// gcc -Wall -O -o dump_udp_ow_4 dump_udp_ow_4.c -lpthread

// maintained by Olaf Wucknitz <wucknitz@mpifr-bonn.mpg.de>



/*
can send test packets like this:
wucknitz@widemap:~$ nc -u4 localhost 4346
(note that tests with netcat may require -4 option)


send TBB data:
observer@lofard3:/media/scratch/observer/wucknitz/lofarD3_local$ socat -b 2140 -u STDIN UDP-DATAGRAM:localhost:4346 < ../TBB_udp/read_udp_gen2.lofard3.2015-03-02T06:17:39




send real beamformed data:

wucknitz@miraculix2:/media/part0/wucknitz/B1508+55_DATA_2017-03-25$ ~/astro_mpi/SOFT/MIRACULIX2/bin/throttle -w 1 -M 10 -s 7824 < dump.bin.2017-03-25-04:00:00.lofarb1.P16011.T15.B1508+55_rcu5 | ~/astro_mpi/SOFT/MIRACULIX2/bin/socat -b 7824 -u STDIN UDP-DATAGRAM:widemap:4346





TBB ports:

#define  PORTSTART  31664 + tbb_no (0-11)

packlen for TBB:

#define  PACKLEN  2140      for transient mode
#define  PACKLEN  2040      for subband mode


BEAMFORMED port:   (depends on station!)
4346 ...

(all on eth1, this may be outdated)
e.g.   tcpdump -i eth1 udp port 4346

 */


// not really sure about these features...
#if 0
#define _XOPEN_SOURCE   600
#define _BSD_SOURCE
#else
#define  _GNU_SOURCE
#endif




// should be the first:  (??)
#include <pthread.h> 

#include  <time.h>

#include  <sys/types.h>
#include  <sys/socket.h>
#include  <assert.h>
#include  <string.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>

#include  <signal.h>

#include <arpa/inet.h>
#include <netinet/in.h>
		   //#include <stdio.h>
		   //#include <sys/types.h>
		   //#include <sys/socket.h>
#include <unistd.h>

#include  <sys/select.h>


#include  <getopt.h>

#include <sys/time.h>
		   //#include <stdio.h>
#include <errno.h>


#define  MAXNSOCK  12

// this is only used for single packet buffers and thus a bit more
// generous
#define  MMAXLEN  10000




pthread_mutex_t region_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t space_available = PTHREAD_COND_INITIALIZER;
pthread_cond_t data_available = PTHREAD_COND_INITIALIZER;
// the last one is also used to stop recording!

/*
  rear is only used in producer, front only in consumer
  size in both
  allbuffs also in both, but there are no conflicts

  region_mutex protects size and the buffer (?)
 */

char *allbuffs;
long size = 0;  /* byte number of full elements */
long front,rear=0;  /* byte queue */
int  maxsock;
int  sock[MAXNSOCK];

int  packlen;

int  do_blocklen= 0;
int  verbose= 0;

int  beamformed_check= 0;


// recording stopped?
// 0: no, 1: for this file, 2: forever
int  stopped= 0;


long bufsize;
long maxsize;

long  maxwrite= 1024*1024;


FILE  *outf;
long  totlen, lasttotlen;
long  packs_seen[MAXNSOCK], packs_dropped[MAXNSOCK], bytes_written[MAXNSOCK],
  beamformed_good_packs[MAXNSOCK];
int  portnos[MAXNSOCK],nsock;
long  beamformed_first_packno[MAXNSOCK], beamformed_last_packno[MAXNSOCK];

long  last_packs_dropped[MAXNSOCK], last_packs_expected[MAXNSOCK], last_packs_seen[MAXNSOCK], last_good_packs[MAXNSOCK];



struct timespec  timeout;


fd_set  allsocks;

char  thisfilename[300], filename[250];

char  *portlist;





// return timestamp for either timestamp as string or yyyy-mm-ddThh:mm:ss
// not yet implemented:  ss.ssss
// <0 for error
double  time_to_timestamp (char  *time)
{
  char  *p;
  double  res;
  struct tm  tm;

  if (strchr (time, 'T')==0)   // then it is a timestamp already
    {
      res= strtod (time, &p);
      if (p[0]!=0)   // still stuff left after conversion
	res= -1;
    }
  else
    {
      memset (&tm, 0, sizeof (tm));
      p= strptime (time, "%Y-%m-%dT%T", &tm);
      if (p==NULL || p[0]!=0)  // not all converted
	res= -1;
      else
	res= timegm (&tm);
    }

  return  res;
}



// currently:  ss.sss
void  timestamp_to_str (double  timestamp, char  *buff, int  len)
{
  int  j;
  time_t  t;
  struct tm  tm;
  
  assert (len>4);
  t= (time_t)timestamp;
  gmtime_r (&t, &tm);
  j= strftime (buff, len-4, "%FT%T", &tm);
  if (j==0)
    {
      fprintf (stderr,
	       "error in strftime() in timestamp_to_str with %e "
	       "for len=%d\n", timestamp, len);
      exit (1);
    }
  sprintf (buff+strlen (buff), ".%03d", (int)((timestamp-t)*1e3));
}


double  realtime ()
{
  struct timeval  tv;
  int  i;

  i= gettimeofday (&tv, NULL);
  if (i!=0)
    {
      perror ("gettimeofday() in realtime()");
      exit (1);
    }

  return tv.tv_sec+1e-6*tv.tv_usec;
}



void  final_statistics ()
{
  long  ntot;
  int  i;

  
  printf ("\ntotal per socket:  (with%s checks for beamformed data)\n",
	  beamformed_check?"":"out");
  for (i= 0; i<nsock; i++)
    {
      if (beamformed_check)
	{
	  ntot= beamformed_last_packno[i]-beamformed_first_packno[i]+1;
	  printf (  "port %5d :  expected packets %9ld\n"
		    "                missed packets %9ld   %10.6f %% of exp\n"
		    "                  seen packets %9ld   %10.6f %% of exp\n"
		    "                  good packets %9ld   %10.6f %% of seen\n"
		    "               dropped packets %9ld   %10.6f %% of seen\n"
		    "               written packets %9ld   %10.6f %% of seen\n"
		    "                                           "
		    "%10.6f %% of exp\n"
		    "                       volume    %7.3f GB\n",
		    portnos[i], ntot,
		    ntot-packs_seen[i], (ntot-packs_seen[i])*100./ntot,
		    packs_seen[i], packs_seen[i]*100./ntot,
		    beamformed_good_packs[i],
		    beamformed_good_packs[i]*100./packs_seen[i],
		    packs_dropped[i], packs_dropped[i]*100./packs_seen[i],
		    packs_seen[i]-packs_dropped[i],
		    (packs_seen[i]-packs_dropped[i])*100./packs_seen[i],
		    (packs_seen[i]-packs_dropped[i])*100./ntot,
		    bytes_written[i]/pow (1024,3));
	}
      else
	{
	  ntot= packs_seen[i];
	  printf (  "port %5d :  seen packets %9ld\n"
		    "           dropped packets %9ld   %10.6f %% of seen\n"
		    "           written packets %9ld   %10.6f %% of seen\n"
		    "                   volume    %7.3f GB\n",
		    portnos[i], ntot,
		    packs_dropped[i], packs_dropped[i]*100./ntot,
		    packs_seen[i]-packs_dropped[i],
		    (packs_seen[i]-packs_dropped[i])*100./ntot,
		    bytes_written[i]/pow (1024,3));
	}
    }
}



// <=0 for no signal
// -1 for timeout
// 0 for regular
void  signal_handler(int signum)
{
  int  i;


  if (signum>0)
    printf ("caught signal %d:  ", signum);

  if (signum<0 && outf==NULL)  // timeout, but no file open
    {
      if (nsock==1 && portnos[0]==0)   // then read stdin
	{
	  printf ("no data on stdin\n");
	  stopped= 2;
	  pthread_cond_signal(&data_available);  // pretend that data is available
	}
      return;
    }

  printf ("\ntotal %7.3f GB  max buff %ld/%ld (%.1f %% full)",
	  totlen/pow(1024,3), maxsize, bufsize,
	  maxsize/(double)bufsize*100.);

  lasttotlen= totlen;

  printf ("\n");

  for (i= 0; i<nsock; i++)
    {
	if (beamformed_check)
	{
	    printf ("port %5d : %8ld exp  %10.6f %% missed  %10.6f %% dropped  "
		    "%7.3f GB\n",
		    portnos[i],
		    beamformed_last_packno[i]-beamformed_first_packno[i]+1,
		    100-packs_seen[i]*100./(
			beamformed_last_packno[i]-beamformed_first_packno[i]+1),
		    packs_dropped[i]*100./packs_seen[i],
		    bytes_written[i]/pow (1024,3));
	    printf ("                           %10.6f %% good\n",
		    beamformed_good_packs[i]*100./packs_seen[i]);

	    printf ("      block: %8ld exp  %10.6f %% missed  "
		    "%10.6f %% dropped\n",
		    beamformed_last_packno[i]-beamformed_first_packno[i]+1 -
		    last_packs_expected[i],
		    100-(packs_seen[i]-last_packs_seen[i])*100. /
		    (beamformed_last_packno[i]-beamformed_first_packno[i]+1
		     -last_packs_expected[i]),
		    (packs_dropped[i]-last_packs_dropped[i])*100./(
			packs_seen[i]-last_packs_seen[i]));
	    printf ("                           %10.6f %% good\n",
		    (beamformed_good_packs[i]-last_good_packs[i])*100./
		    (packs_seen[i]-last_packs_seen[i]));
	    
	    
	    last_packs_expected[i]= beamformed_last_packno[i]-
		beamformed_first_packno[i]+1;
	    last_good_packs[i]= beamformed_good_packs[i];
	}
	else
	{
	    printf ("port %5d : %8ld seen  %10.6f %% dropped  "
		    "%7.3f GB\n",
		    portnos[i], packs_seen[i], // packs_dropped[i],
		    packs_dropped[i]*100./packs_seen[i],
		    bytes_written[i]/pow (1024,3));
	    printf ("      block: %8ld seen  %10.6f %% dropped\n",
		    packs_seen[i]-last_packs_seen[i],
		    //packs_dropped[i]-last_packs_dropped[i],
		    (packs_dropped[i]-last_packs_dropped[i])*100./(
			packs_seen[i]-last_packs_seen[i]));
	}
    last_packs_dropped[i]= packs_dropped[i];
    last_packs_seen[i]= packs_seen[i];
    }

  if (signum==SIGINT || signum==SIGALRM)
    {
      printf ("stopping (%s)\n", signum==SIGINT ? "INT" : "end time");

      stopped= 2;
      pthread_cond_signal(&data_available);  // pretend that data is available

    }
  else if (signum==-1 || signum==SIGHUP)
    {
      if (outf)  
        {
          if (signum<0)
	    {
	      if (nsock==1 && portnos[0]==0)   // then read stdin
		{
		  printf ("no more data on stdin\n");
		  stopped= 2;
		}
	      else
		{
		  printf ("timeout\n");
		  if (stopped==0)  // may already be 2
		      stopped= 1;
		}
	    }
	  else
	      if (stopped==0)  // may already be 2
		  stopped= 1;
	  pthread_cond_signal(&data_available);// pretend that data is available

        }
    }


}


// header_lofar of each packet:
struct __attribute__((__packed__))  header_lofar
{
  uint8_t   version;
  union __attribute__ ((__packed__)) {
    struct __attribute__ ((__packed__)) {
      unsigned int  rsp_id   : 5;
      unsigned int  unused1  : 1;
      unsigned int  error    : 1;
      unsigned int  is200mhz : 1;
      unsigned int  bm       : 2;
      unsigned int  unused2  : 6;
    } source;
    uint8_t  source_int;
  };
  uint8_t   config;
  uint16_t  station;
  uint8_t   num_beamlets, num_slices;
  //  uint32_t  timestamp, sequence;
  int32_t  timestamp, sequence;

  //  int8_t   data[BEAMLETS][SLICES][4];  // 4 : X/Y  R/I
};



long  beamformed_packno (struct header_lofar  *header)
{
  return ((header->timestamp*1000000l*(160+40*header->source.is200mhz)+512)/1024+header->sequence)/16;
}


int  beamformed_checkpack (struct header_lofar  *header)
{
  return header->source.error==0 && header->timestamp!=-1;
}


void *producer ()
{
  char  buff[MMAXLEN+2];
  socklen_t  slen;
  struct sockaddr_in  addr_src;
  fd_set  myallsocks;
  int  i, thissize, thisblock; 



  char  *buff2;




  if (do_blocklen)
      buff2= buff+2;  // we need two bytes for size
  else
      buff2= buff;


  slen= sizeof (addr_src);
  while (1)
    {


      if (totlen-lasttotlen>1e9)
	signal_handler (0);  // regular printouts

      thissize= -999;  // prevent compiler warning on some systems

      if (nsock==1 && portnos[0]==0)   // then read stdin
	{
	  if (stopped)
	    thissize= 0;
	  else
	    {


	      if (size+packlen >= bufsize)   // not enough space
		{
		  // wait till space available in buffer
		  // (for stdin we can wait with reading, don't drop packets)
		  printf ("waiting for buffer space...\n");
		  pthread_mutex_lock(&region_mutex);
		  pthread_cond_wait(&space_available,&region_mutex); 
		  pthread_mutex_unlock(&region_mutex);
		}


	      thissize= fread (buff2, 1, packlen, stdin);
	      if (ferror (stdin))
		perror ("reading from stdin in producer()");

	      if (thissize==0)  // treat as timeout
		{
		  signal_handler (-1);
		  //stopped= 2;
		}
	    }
	  // here we already have read the packet (or thissize==0)
	}
      else
	{


	    if (stopped==2)
	    {
		// close all sockets and return
		assert (nsock!=1 || portnos[0]!=0);  // not reading from stdin
		for (i= 0; i<nsock; i++)
		{
		    int  j;
		    
		    j= close (sock[i]);
		    if (j)
		    {
			perror ("closing socket");
			exit (1);
		    }
		}
		return NULL;
	    }


	  myallsocks= allsocks;  // we have to reset this every call
	  
	  i= pselect (maxsock+1, &myallsocks, NULL, NULL,
		      &timeout,
		      NULL);
	  if (i==-1)
	    {
	      perror ("pselect in producer()");
	      exit (1);
	    }
	  if (i==0)  // timeout
	    {
	      signal_handler (-1);
	    }
	}
      
      for (i= 0; i<nsock; i++)
	{
	  if (nsock!=1 || portnos[0]!=0)  // not reading from stdin
	    {
	      if (FD_ISSET (sock[i], &myallsocks))
		{
		  thissize= recvfrom (sock[i], buff2, MMAXLEN-1, // play safe
				      0,
				      (struct sockaddr *)&addr_src,
				      &slen);
		  if (thissize==-1)
		    {
		      perror ("recvfrom() in producer()");
		      exit (1);
		    }
		  if (thissize>=MMAXLEN)  // this should not happen
		    {
		      fprintf (stderr, "producer(): recvfrom() result %d >= %d"
			       " (should not happen)\n", thissize, MMAXLEN);
		      exit (1);
		    }
		  
		}
	      else
		thissize= 0;
	    }
	  if (thissize)
	  {
	      if (stopped)
		  ;  //printf ("discarding packet\n");
	      else
	      {
		  // now we have a packet either from socket or from stdin
		  if (packlen==0 || thissize==packlen)
		    
		    {

			if (do_blocklen)//add the blocklen if wanted (two bytes)
			    // not tested
			{
			    *(uint16_t*)buff= (uint16_t)thissize;
			    //thissize+= 2;  see below
			}

		      
		      if (do_blocklen)// add the blocklen if wanted (two bytes)
			thissize+= 2;
		      
		      
		      if (beamformed_check)
			{
			  beamformed_last_packno[i]= beamformed_packno (
					(struct header_lofar*)buff2);
			  if (beamformed_first_packno[i]==-1)
			    beamformed_first_packno[i]= 
			      beamformed_last_packno[i];
			  if (beamformed_checkpack ((struct 
						     header_lofar*)buff2))
			    beamformed_good_packs[i]++;
			}
		      
		      packs_seen[i]++;
		      
		      // (no conflicts for rear and allbuffs)
		      
		      if (size+thissize >= bufsize)   // not enough space
			
			{
			  // simply drop this packet
			  packs_dropped[i]++;
			}
		      else // enough space
			
			{
			  if (rear+thissize<=bufsize)   // no wrap necessary
			    memcpy (allbuffs+rear, buff, thissize);
			  else  // we have to wrap
			    {
			      thisblock= bufsize-rear;
			      memcpy (allbuffs+rear, buff, thisblock);
			      memcpy (allbuffs, buff+thisblock,
				      thissize-thisblock);
			    }
			  rear= (rear+thissize)%bufsize;
			  
			  // ... and then lock again here 
			  pthread_mutex_lock(&region_mutex);
			  
			  size+= thissize;

			  if (size>maxsize)
			      maxsize= size;

			  pthread_cond_signal(&data_available);
			  pthread_mutex_unlock(&region_mutex);
			  
			  totlen+= thissize;
			  
			  bytes_written[i]+= thissize;
			}
		      
		    }  // packlen==0 || ...
		  else
		    printf ("received %5d bytes, wrong length in sock %d, "
			    "should be %d\n", thissize, i, packlen);
		}   // not stopped
	    }  // thissize!=0
	}  // for (i= 0; i<nsock ....
    } // while (1)
}



void  init_thisfilestat ()
{
  int  j;

  for (j= 0; j<nsock; j++)
    {
      beamformed_first_packno[j]= -1;
      bytes_written[j]= packs_seen[j]= packs_dropped[j]= 
	beamformed_good_packs[j]= 0;

      last_packs_dropped[j]= last_packs_expected[j]= last_packs_seen[j]= 
	  last_good_packs[j]= 0;
    }

  lasttotlen= totlen= 0;
}


void  start_file (double  timestamp, char  *comment)
{
  char  buff[100];

  printf ("start file\n");
  timestamp_to_str (timestamp, buff, sizeof (buff));
  if (strcmp (filename, "/dev/null")==0)
    {
      strcpy (thisfilename, filename);
      printf ("\nopening %s\n", thisfilename);
    }
  else
    {
      if (comment)
	sprintf (thisfilename, "%s_%s.%s.%s", filename, portlist, comment,
		 buff);
      else
	sprintf (thisfilename, "%s_%s.%s", filename, portlist, buff);
      printf ("\ncreating %s\n", thisfilename);
    }
  
  outf= fopen (thisfilename, "w");
  if (outf==NULL)
    {
      perror ("opening output file in start_file()");
      exit (1);
    }
}


void *consumer ()
{
    long  i, thissize, thisblock;

  while (1)
    {

      pthread_mutex_lock(&region_mutex);
      if (size == 0 && ! stopped)  // no data available
	pthread_cond_wait(&data_available,&region_mutex); 



      if (size==0 && stopped)  // buffer empty and we want to stop
	{
	  int  j;

	  if (outf)
	    {
	      final_statistics ();
	      printf ("closing %s\n", thisfilename);
	      j= fclose (outf);
	      if (j!=0)
		{
		  perror ("closing file in consumer()");
		  exit (1);
		}
	      outf= NULL;
	      init_thisfilestat ();
	    }
	  if (stopped==2)
	    exit (0);

	  stopped= 0;
	  pthread_mutex_unlock(&region_mutex);
	  continue;
	}

      
     
      thissize= size;


      if (outf==NULL)  // no file open, open it now
        // producer can produce in parallel
        {
	  start_file (realtime (), "packets");
        }

      // here I can unlock ...
      // (no conflicts for front and allbuffs)
      pthread_mutex_unlock(&region_mutex);


      // now write to disk

      if (thissize>maxwrite)
	  thissize= maxwrite;



#if 0
      // some delay for tests
      {
	struct timespec  timespec;
	
	timespec.tv_sec= 0;
	timespec.tv_nsec= (long) (1e-3  *1e9);
	nanosleep (&timespec, NULL);
      }
#endif

      if (front+thissize<=bufsize)     // no wrap necessary
        {
          i= fwrite (allbuffs+front, 1, thissize, outf);
	  if (i!=thissize)
	    {
	      perror ("writing file (1) in consumer()");
	      exit (1);
	    }
        }
      else  // we have to wrap
        {
          // first from here to end
          thisblock= bufsize-front;
          i= fwrite (allbuffs+front, 1, thisblock, outf);
	  if (i!=thisblock)
	    {
	      perror ("writing file (2) in consumer()");
	      exit (1);
	    }
          // then the rest
          thisblock= thissize-thisblock;
          i= fwrite (allbuffs, 1, thisblock, outf);
	  if (i!=thisblock)
	    {
	      perror ("writing file (2) in consumer()");
	      exit (1);
	    }
        }
      front= (front+thissize)%bufsize;
      
      
      // ... and then lock again here (for size)
      pthread_mutex_lock(&region_mutex);
      size-= thissize;

      pthread_cond_signal(&space_available);
      pthread_mutex_unlock(&region_mutex);
    }
}

 




int  main (int  argc, char  **argv)
{
  int  i, j, c;
  struct sockaddr_in  addr[MAXNSOCK];

  pthread_t producer_thread; 
  pthread_t consumer_thread; 

  struct option long_options[] =
    {
      {"verbose",  no_argument, &verbose, 1},
      {"len",      required_argument, 0, 'l'},
      {"ports",    required_argument, 0, 'p'},
      {"out",      required_argument, 0, 'o'},
      {"help",     no_argument,       0, 'h'},
      {"Help",     no_argument,       0, 'H'},
      {"sizehead", no_argument,       0, 's'},
      {"bufsize",  required_argument, 0, 'b'},
      {"maxwrite", required_argument, 0, 'm'},
      {"timeout",  required_argument, 0, 't'},
      {"Start",    required_argument, 0, 'S'},
      {"End",      required_argument, 0, 'E'},
      {"duration", required_argument, 0, 'd'},
      {"check",    no_argument,       0, 'c'},
      {0, 0, 0, 0}
    };
  char  *short_options= "hHvl:p:o:sb:m:t:S:E:d:c";

  int  option_index= 0;

  char  stdportlist[]= "4346";
  char  *start_time= NULL, *end_time= NULL;
  double  start_timestamp= 0, end_timestamp= 0;
  double  duration= 0;

  char *cp, *cp2, *cp3, *cp4;


  char  buff2[100];

  double  timeout_sec= 10.0;


  portlist= stdportlist;


  i= gethostname (buff2, sizeof (buff2));
  if (i!=0)
    {
      strcpy (buff2, "unknown");
      fprintf (stderr, "cannot determine hostname, using %s", buff2);
      perror ("gethostname");
    }
  sprintf (filename, "udp.%s", buff2);


  outf= NULL;



  bufsize= 104857600;
  packlen= 0;  // arbitrary

  while (1)
    {
      c = getopt_long (argc, argv, short_options,
		       long_options, &option_index);

      if (c==-1)  // no more options
        {
          // check that there are no other arguments
          if (argc>optind)
            {
              fprintf (stderr, "no other arguments allowed\n");
              c= '?'; // error
            }
          else
            break;
        }
      
      switch (c)
        {
          case 0:  // for the flag
            break;
          case 'v':
            verbose= 1;
            break;
          case 'l':
            if (sscanf (optarg, "%d", &packlen)!=1 ||
                packlen<=0 || packlen>=MMAXLEN)
              {
                fprintf (stderr, "problem with packet length\n");
                c= '?';
              }
	    if (beamformed_check && packlen!=7824)
	      {
		fprintf (stderr, "--check implies --len 7824, cannot use "
			 "other value\n");
		c= '?';
	      }
            break;
          case 'b':
            if (sscanf (optarg, "%ld", &bufsize)!=1 ||
                bufsize<=10000)
              {
                fprintf (stderr,
                         "problem with bufsize\n");
                c= '?';
              }
            break;
          case 'm':
            if (sscanf (optarg, "%ld", &maxwrite)!=1 ||
                maxwrite<=1024)
              {
                fprintf (stderr,
                         "problem with maxwrite\n");
                c= '?';
              }
            break;
          case 't':
            if (sscanf (optarg, "%lf", &timeout_sec)!=1 ||
                timeout_sec<1e-3)
              {
                fprintf (stderr, "problem with timeout\n");
                c= '?';
              }
            break;
          case 'o':
            strncpy (filename, optarg, sizeof (filename)-1);
            filename[sizeof (filename)-1]= 0;
            break;
          case 'p':
            portlist= optarg;
            break;
          case 's':
            do_blocklen= 1;
            break;
          case 'S':
            start_time= optarg;
	    start_timestamp= time_to_timestamp (start_time);
	    if (start_timestamp<0)
              {
                fprintf (stderr, "problem with start time\n");
                c= '?';
              }
            break;
          case 'E':
            end_time= optarg;
	    end_timestamp= time_to_timestamp (end_time);
	    if (end_timestamp<0)
              {
                fprintf (stderr, "problem with end time\n");
                c= '?';
              }
	    if (duration)
	      {
		fprintf (stderr, "cannot use --End and --duration together\n");
		c= '?';
	      }
            break;
	  case 'd':
            if (sscanf (optarg, "%lf", &duration)!=1 ||
                duration<=0)
              {
                fprintf (stderr,
                         "problem with duration\n");
                c= '?';
              }
	    if (end_time)
	      {
		fprintf (stderr, "cannot use --End and --duration together\n");
		c= '?';
	      }
            break;
	  case 'c':
	    if (packlen && packlen!=7824)
	      {
		fprintf (stderr, "--check implies --len 7824, "
			 "cannot use other value\n");
		c= '?';
		break;
	      }
	    packlen= 7824;
	    beamformed_check= 1;
	    break;
          case 'h':
          case 'H':
            break;
          default:  // incl '?'
            c= '?';
        }
      if (c=='?' || c=='h' || c=='H')  // some error or help
        {
          fprintf (stderr, "\n%s  options\n"
                   "    [--ports/-p  portlist]   current: %s\n"
                   "                             e.g.  31664,31665 or 31664x2\n"
		   "                             or 0 for stdin read\n"
                   "    [--out/-o filename]      current: %s\n"
                   "    [--verbose/-v] \n"
                   "    [--len/-l packet_len]    current: %d, 0=arbitrary\n"
                   "    [--sizehead/-s]          write packet lengths as headers\n"
		   "                             (not well tested)\n"
                   "    [--bufsize/-b size]      current: %ld\n"
		   "    [--maxwrite/-m size]     max. write block, current: %ld\n"
                   "    [--timeout/-t sec]       current: %f\n"
                   //		   "    [--out/-o filename]\n"
		   "    [--Start/-S time]        default: now\n"
		   "    [--End/-E time]          default: never\n"
		   "                             time: unix-timestamp or yyyy-mm-ddThh:mm:ss\n"
		   "    [--duration/-d sec]      default: infinity\n"
		   "                             (from start time or first packet)\n"
		   "    [--check/-c]             packet statistics for beamformed data\n"
		   "                             implies --len 7824\n"
// not yet implemented:  ss.ssss
                   "    [--help/-h]\n"
                   "    [--Help/-H]              extended help\n",
		   argv[0], portlist, filename, packlen,
                   bufsize, maxwrite, timeout_sec);
	  if (c=='H')
	    fprintf (stderr,
"\nWe can work in different modes. If --Start is given, start at that time,\n"
"otherwise with first arriving packet. If --End is given, stop at that time.\n"
"If --duration is given, run for that long. This duration either starts at\n"
"--Start or with first packet. --timeout stops recording after that time\n"
"with no packets. If --Start used, timeout can also happen before first\n"
"packet, otherwise only once data have arrived. After timeout the programme\n"
"stops this recording but then waits for next packet and potentially starts\n"
"new file(s). After --duration or at --End, the programme stops.\n"
"We can listen to several ports, but all data will go to one file.\n"
"--ports 0 reads from stdin. It requires --len but cannot use --Start, --End\n"
"or --duration. End of file is treated as timeout.\n"
"Filename is built from --out parameter plus 'start' or 'packet' (depending\n"
"on whether we start at certain time or with first packet) plus UTC timestamp.\n"
"Filename '/dev/null' (this exact spelling) is used directly.\n"
"Packets can be any length, unless --len is given, then only that length is\n"
"accepted (others discarded). For variable packet length we can write the\n"
"lengths as headers (--sizehead). The internal ring buffer size can be set\n"
"with --bufsize. --verbose produces more output.\n"
"Reading and writing have their own threads, data are written in maximum\n"
"blocks given by --maxwrite. (Should be << bufsize, because each block\n"
"is only released after complete write.)\n"
"With --check we compare the number of packets (received and written) with\n"
"the number expected from the packet numbers and determine a completeness.\n"

);

          exit (c=='?');
        }

    }  // while (1)


	
  
  timeout.tv_sec= timeout_sec;
  timeout.tv_nsec= (int)((timeout_sec-timeout.tv_sec)*1e9+0.5);

  if (verbose)
    {
      char  buff[100];

      printf ("packlen %d\n", packlen);
      printf ("filename %s\n", filename);
      printf ("portlist %s\n", portlist);
      printf ("timeout %.6f sec\n", timeout_sec);
      if (start_time)
	{
	  timestamp_to_str (start_timestamp, buff, sizeof (buff));
	  printf ("start time %.3f = %s\n", start_timestamp, buff);
	}
      if (end_time)
	{
	  timestamp_to_str (end_timestamp, buff, sizeof (buff));
	  printf ("end time   %.3f = %s\n", end_timestamp, buff);
	}
      if (duration)
	printf ("duration %.3f sec\n", duration);

      if (beamformed_check)
	printf ("check%s beamformed statistics\n",
		beamformed_check==2 ? "extended" : "");
    }


  allbuffs= malloc (bufsize);
  if (allbuffs==NULL)
    {
      fprintf (stderr, "cannot allocate memory for buffer (%ld bytes)\n",
	       bufsize);
      exit (1);
    }

  maxsize= 0;
  lasttotlen= 0;




  cp2= NULL;  // prevent warnings

  nsock= 0;
  cp= strtok_r (portlist, ",", &cp2);
  while (cp!=NULL)
  {
      cp3= strchr (cp, 'x');
      if (cp3==0)
      {
	if (nsock>=MAXNSOCK)
	  {
	    fprintf (stderr,
		     "number of sockets too large (>%d, allowed max. %d)\n",
		     nsock, MAXNSOCK);
	    exit (1);
	  }
	assert (cp[0]!=0);
	i= strtol (cp, &cp4, 10);
	assert (cp4[0]==0);
	portnos[nsock]= i;
	nsock++;
      }
      else
      {
	  assert (cp[0]!=0);
	  i= strtol (cp, &cp4, 10);
	  assert (cp4[0]=='x');
	  assert (cp3[1]!=0);
	  j= strtol (cp3+1, &cp4, 10);
	  assert (cp4[0]==0);
	  while (j)
	  {
	    if (nsock>=MAXNSOCK)
	      {
		fprintf (stderr,
			 "number of sockets too large (>%d, allowed max. %d)\n",
			 nsock, MAXNSOCK);
		exit (1);
	      }
	    portnos[nsock]= i;
	    i++;
	    nsock++;
	    j--;
	  }
	      
      }
      cp= strtok_r (NULL, ",", &cp2);
  }



  if (verbose)
      for (i= 0; i<nsock; i++)
        printf ("port %d  %d\n", i, portnos[i]);

  if (nsock==1 && portnos[0]==0)   // then read stdin
    {
      if (packlen==0)
	{
	  fprintf (stderr, "Reading from stdin requires --len.\n");
	  exit (1);
	}
      if (start_timestamp || end_timestamp || duration)
	{
	  fprintf (stderr, "Reading from stdin is not compatible with "
		   "--Start, --End, --duration.\n");
	  exit (1);
	}
    }


  if (start_timestamp)
    {
      double  wait_time;
      struct timespec  timespec;


      start_file (start_timestamp, "start");

      wait_time= start_timestamp-realtime ();

      printf ("waiting for %.3f sec...\n", wait_time);

      if (wait_time<0)
	{
	  printf ("negative wait, starting now!\n");
	  if (duration)
	    end_timestamp= realtime ()+duration;
	}
      else
	{
	  if (duration)
	    end_timestamp= start_timestamp+duration;
	  while (wait_time>0)
	    {
	      if (wait_time>=1)
		sleep ((unsigned int)wait_time);
	      else
		{
		  timespec.tv_sec= (time_t)wait_time;
		  timespec.tv_nsec= (long)((wait_time-timespec.tv_sec)*1e9+0.5);
		  nanosleep (&timespec, NULL);
		}
	      wait_time= start_timestamp-realtime ();
	    }
	  if (verbose)
	    printf ("remaining wait_time = %.6f sec\n", wait_time);
	}
      
    }
  else
    if (duration)
      end_timestamp= realtime ()+duration;

  if (end_timestamp)
    {
      double  wait_time;
      struct itimerval  itimer;

      wait_time= end_timestamp-realtime ();
      printf ("running for max %.3f sec...\n", wait_time);
      if (wait_time<0.1)
      {
	  printf ("time is%s negative, do not record at all\n",
		  wait_time>=0 ? " almost" : "");
	  exit (1);
      }
      else
	{
	  itimer.it_value.tv_sec= (long)wait_time;
	  itimer.it_value.tv_usec= (long)((wait_time-
					   itimer.it_value.tv_sec)*1e6+0.5);
	  itimer.it_interval.tv_sec= itimer.it_interval.tv_usec= 0;

	  j= setitimer (ITIMER_REAL, &itimer, NULL);
	  if (j!=0)
	    {
	      perror ("setitimer()");
	      exit (1);
	    }
	  signal (SIGALRM, signal_handler);
	}
    }



  maxsock= -1;
  FD_ZERO (&allsocks);

  if (nsock==1 && portnos[0]==0)   // then read stdin
    printf ("reading from stdin\n");
  else
    {
      printf ("listening to %s\n", portlist);
      for (i= 0; i<nsock; i++)
	{
	  sock[i]= socket (AF_INET, SOCK_DGRAM, IPPROTO_UDP);
	  if (sock[i]==-1)
	    {
	      perror ("socket()");
	      exit (1);
	    }
	  //      printf ("sock %d %d\n", i, sock[i]);
	  if (sock[i]>maxsock)
	    maxsock= sock[i];
	  
	  memset(&addr[i], 0, sizeof(addr[i]));
	  addr[i].sin_family= AF_INET;
	  addr[i].sin_port= htons (portnos[i]);
	  addr[i].sin_addr.s_addr= htonl (INADDR_ANY);
	  
	  j= bind (sock[i], (struct sockaddr *)&addr[i], sizeof (addr[i]));
	  if (j==-1)
	    {
	      perror ("bind()");
	      exit (1);
	    }
	  
	  FD_SET (sock[i], &allsocks);
	}
    }



  init_thisfilestat ();


  signal (SIGINT, signal_handler);
  signal (SIGHUP, signal_handler);

  totlen= 0;

  j= pthread_create(&consumer_thread,NULL,consumer,NULL);
  if (j)
  {
      perror ("pthread_create for consumer");
      exit (1);
  }
  j= pthread_create(&producer_thread,NULL,producer,NULL);
  if (j)
  {
      perror ("pthread_create for producer");
      exit (1);
  }

  j= pthread_join(consumer_thread,NULL);
  if (j)
  {
      perror ("pthread_join");
      exit (1);
  }

  free (allbuffs);

  return 0;
}
