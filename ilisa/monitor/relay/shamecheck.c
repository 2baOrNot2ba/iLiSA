/*****************************************************************************
 *
 * FILE: shamecheck.c
 *
 *   This program can be used to check the reception of shamecasts. It will
 *   connect to the shamecast socket and then print out the names, sizes,
 *   version numbers and time stamps of all shamecast packages arriving.
 *
 * HISTORY
 *
 * who          when           what
 * ---------    -----------    ----------------------------------------------
 * lerner       20 May 2014    Original coding
 *
 *****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>


#define MULTICAST_ADDRESS   0xe0010104
#define MULTICAST_PORT            4242


/*  Define a general shamecast block that we can use for reading arbitrary
    shamecast blocks  */

struct GEN_SHAME_BLOCK {
  char name[16];
  int size;
  int version;
  int timestamp;
  int flag;
  char data[10000];
} gen_block;


const int ttl = 1;



int main(int argc, char *argv[]) {

  struct sockaddr_in from;
  struct ip_mreq multiaddr;
  socklen_t fromlen;
  int sock;
  int itmp;
  char origin[20];
  time_t now;
  char timestring[80];

/*  Set up the socket where we will listen for shamecasts  */

  sock = socket(AF_INET, SOCK_DGRAM, 0);

  itmp = 1;
  if ( setsockopt(sock, SOL_SOCKET, SO_REUSEADDR, &itmp,
                  sizeof(int)) < 0 ) {
    perror("setsockopt SO_REUSEADDR failed");
    exit(1);
  }

  multiaddr.imr_multiaddr.s_addr = htonl(MULTICAST_ADDRESS);
  multiaddr.imr_interface.s_addr = htonl(INADDR_ANY);

  if ( setsockopt(sock, IPPROTO_IP, IP_ADD_MEMBERSHIP, &multiaddr,
                  sizeof(multiaddr)) < 0 ) {
    perror("setsockopt IP_ADD_MEMBERSHIP failed");
    exit(1);
  }

  if ( setsockopt(sock, IPPROTO_IP, IP_MULTICAST_TTL, &ttl,
                  sizeof(ttl)) < 0 ) {
    perror("setsockopt IP_MULTICAST_TTL failed");
    exit(1);
  }

/*  Set up the source address and UDP port  */

  from.sin_family = AF_INET;
  from.sin_addr.s_addr = htonl(MULTICAST_ADDRESS);
  from.sin_port = htons(MULTICAST_PORT);

  if ( bind(sock, (struct sockaddr *) &from,
            sizeof(struct sockaddr_in) ) < 0 ) {
    perror("bind socket failed");
    exit(1);
  }

/*  Start the eternal shamechecking loop --- we read in the shamecast blocks
    into the general block 'gen_block' and then print out information from
    the headers  */

  while ( 1 ) {
    fromlen = sizeof(struct sockaddr_in);
    if ( ( recvfrom(sock, &gen_block, sizeof(gen_block),
                    MSG_WAITALL, (struct sockaddr*) &from, &fromlen)) < 0 )
      continue;
    now = time(NULL);
    strcpy(timestring, ctime(&now));
    strncpy(origin, gen_block.name, 16);
    printf("Received block '%-16s' (%4d) ver. %d  -  %s", origin,
           gen_block.size, gen_block.version, timestring);
  }
}

