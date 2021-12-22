# *iStnMonitor* a software package for monitoring international LOFAR stations remotely

v1.2

2013-09-30

Tobia Carozzi (<tobia@chalmers.se>) and Simon Casey

Onsala Space Observatory, Sweden

Changes in v1.2
---------------

-   The paths to standard Lofar commands have been updated to reflect
    changes in the MAC version 1.16 (currently latest). The paths have
    been placed in global variables, but users shouldn't have to change
    these unless they are using a MAC version prior to 1.16. (Before
    1.16 the path to the environment control scripts were in
    */opt/stationtest/test/envcontroltest/* )

Changes in v1.1
---------------

-   Cacti scripts have been updated to resolve incorrect correspondence
    between graphs labels and variables
-   stnStatMon.py has been updated to reflect changes in *isStatus.py
    *output strings (“cab3” had been removed 21 June 2013)

Introduction
------------

Although international LOFAR stations come equipped with numerous
scripts for monitoring the status of the station, it is often useful to
get the status information out of the LCU network, e.g., on a computer
dedicated to centralized monitoring of observatory resources. The
software package dubbed *iStnMonitor *provides this functionality.

Provides
--------

The software provides a webpage with plots over time of environment
control:

-   Temperature (in container)
-   Humidity (in container)
-   48V power supply ON/OFF status

and in addition

-   swlevel
-   ILT switch ON/OFF status

Prerequisites
-------------

In order to use the software, hardware-wise you need

-   international LOFAR LCU
-   Gateway (GW) PC (usually a computer with 2 ethernet cards: one to
    LCU & one to observatory network)

Software-wise you need python and cacti. It would be convenient if the
GW PC had ubuntu linux.

Software contents
-----------------

The *iStnMonitor *package consists of at least

-   stnStatMon.py

    -   which collates LCU status information
    -   runs on LCU

-   LCUcli.py

    -   which reads status info and converts it for use in cacti
    -   runs on GW PC or other PC on observatory network

-   LOFAR\_cacti\_template\_Fat.xml and LOFAR\_cacti\_templateNoFat.xml

    -   which are xml descriptions of how to plot monitoring data on
        webpage

-   mkRRD.sh

    -   which creates the initial Round Robin Database files

Optionally you could use

-   stnMonitorRelay.py

    -   a script that relays UDP packets from LCU to cacti server via GW
        PC

-   FatterLines.sh

    -   Makes thicker lines in Cacti plots

Installation
------------

### On LCU

-   copy over *stnStatMon.py *and *monitor\_crontabs.txt *to LCU
-   Edit *stnStatMon.py *to reflect the IP addresses and port number of
    your gateway PC (specifically: edit global variable *GW\_PC\_UDP\_IP
    *)
-   Create a crontab entry (crontab -e at prompt) and type\
    *\*/5 \* \* \* \*\$HOME/iStnMonitor/stnStatMon.py &gt;&gt;
    stnStatus.log\
    *for updates every five minutes. Note that \$HOME is where you
    installed the iStnMonitor software. The redirect at the end is an
    example and can edited to suit tastes (and can even be removed), it
    simply provides a log file of the status data that was sent to the
    GW PC, so you should have this data anyways.

### 

### 

### On GW PC (optional)

-   Either copy over *stnMonitorRelay.py *and run this by putting into
    */etc/rc.local*
-   Or configure using *iptables *to reroute packets from LCU to cacti
    server

### 

### On cacti server

This section aims to provide the reader with a guide to install a
web-based platform for viewing monitoring parameters of a LOFAR station.
It is dependant on a package called Cacti which provides an easy to use
interface to a database and graphing package called rrdtool.

First of all you will need to install Cacti, in Ubuntu, *apt-get install
cacti *will install and configure Cacti and its dependencies.

Once Cacti is installed, you can create the empty Round Robin Database
files which will be used for storing the logged data. Execute the
mkRRD.sh script in a directory to which you have write access. This will
create two files, *LOFAR\_status.rrd* and *LOFAR\_environment.rrd*. You
will then need to move these to a suitable permanent location where the
LCUcli.py script has write access to them, and the apache2 process has
read access. On Ubuntu, the default location for RRD files is
*/var/lib/cacti/rra*, so this may be a good location to choose, but
remember to check the permissions, and note down the location you choose
as it will be needed later on.

To enable fatter lines to be displayed in the graphs, some modifications
need making to a few of the Cacti files. On Ubuntu, Cacti is installed
to */usr/share/cacti/site *and from within that directory, execute the
*FatterLines.sh *script. This requires that you have sufficient
permissions to modify files within this directory. If you don't, then
you will need to use the *LOFAR\_cacti\_template\_NoFat.xml* script in
the next steps, otherwise you should use the script
*LOFAR\_cacti\_template\_Fat.xml*.

Assuming you are using a newly installed Cacti, open a web browser and
go to *http://*cactiserver*/cacti *where you will be asked to login. The
default login on Ubuntu is admin/admin, and on the following page you
will be asked to choose a new password.

From the menu on the left, navigate to Import/Export -&gt; Import
Templates. Select the appropriate LOFAR template file, and choose “Use
custom RRA settings from template”.

![Figure 1: Adding a 'LOFAR station'
device](ex_cacti_add_station.png){width="17.59cm"
height="13.085cm"}

With the template imported, you now need to configure a device to use
the template. Click **Management -&gt; Devices** and you should have a
single device listed, ‘localhost’. It’s unlikely you want to monitor the
machine running cacti, so check the box to the right of this, and choose
the delete action then click go to remove this. Now back in **Management
-&gt; Devices**, choose ‘add’. Fill in the details as in Figure 1. The
description can be anything appropriate, probably your station code is a
good idea. Click create.

Next, click ‘Create graphs for this host’. Check the boxes next to ‘LCU
environment’ & ‘LCU status’ and click create, then create again on the
next page.

![Figure 2: Adding graphs to a
tree](ex_cacti_add_graphs.png){width="17.59cm"
height="13.109cm"}

From the left menu, click **Management -&gt; Graph Management**. If no
rows are visible in the table, i.e. you don’t see the two LCU lines as
in Figure 2, click Go in the “Graph Management” search box. Check the
box beside each LCU line, then choose action “Place on a tree (Default
Tree)” and click Go. On the next page, use the default Destination
Branch \[root\] and click continue.

Finally, Cacti needs informing of where to source its data for the
graphs. Click **Management -&gt; Data Sources**. Again if no rows are
visible, click Go in the “Data Sources” search box. For each of the two
data sources, click on the name and on the next page enter the location
to the appropriate RRD file that you created earlier, then click save.

With all this done, you should now be able to click on the Graphs tab,
and be presented with 2 (empty) graphs.

From here on, you need to deploy the LCUcli.py on the machine hosting
cacti. For this, you will need to take a note of the location of the two
.rrd files you created earlier. Again, double check the permissions,
they need to be readable by the httpd process, and writable by the user
id that the LCUcli.py script will be running as.
