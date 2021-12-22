# Setting Up Cacti for a GLOW Station

## Import Template(s)

Import the template for a GLOW station `iStnMonitor/client/LOFAR_cacti_template_GLOW.xml`. 
Make sure to tick the box next to `Use custom RRA settings from the template`, or you Cacti not 
use the high-resolution data but will do increasingly more averaging for old data.

In case you also want to display the data from the barix device in Effelsberg, also import 
the `iStnMonitor/client/LOFAR_cacti_template_DE601-brix.xml` template.

## Set up a new Host for the LOFAR Station

1. Go to `Create` -> `New Graphs` (from the menu on the left), select `Create New Host`.
1. Fill out the form:
    1. Give the new device the name of the station, e.g. `DE601 - Effelsberg`
    1. Set the `Hostname` (IP-address) to `localhost`
    1. Set the `Host Template` to `LOFAR Station GLOW`
    1. Set `Downed Device Detection` to `None`
    1. Set `SNMP Version` to `Not In Use`
    1. Click `Create`
1. Create the Graphs for the new station:
    1. Go to the device page for the station, either by continuing from above, or go to `Management` -> `Devices` from the menu
    1. Select `Create Graphs for this Host`
    1. Select all (or those graphs that you want to use)
    1. Click `Create`
1. Tell Cacti to read the data from the right files
    1. Go to `Management` -> `Data Sources`
    1. Select the correct Host at the top under `Data Sources` -> `Host:`. You should see one entry for each data source (`.rrd` file) that is associated with this host.
    1. For each of the data sources:
        1. Select the data source
        1. Enter the path to the `.rrd` file in the `Data Source Path` field
        1. `Save`
Now all the graphs for the station should be available on the list view of the graphs.

## Set up the Tree View for the Station

All the graphs should be visible on the list view.
1. Got to `Management` -> `Graph Trees`
1. Select the existing "Default" tree od `Add` a new tree
1. Give the tree a useful name / change the name of the "Default" tree
1. Set up the new tree. I do:
    1. Add a graph `LCU environment` (I usually don't change the `Round Robin Archive` value. I don't know what it does here.)
    1. Add a graph `LCU status`
    1. Add a header named `HW-Status` (Keep it a `Manual Ordering`)
    1. In this new header I add (use the `add` link next to the header name):
        1. Add a graph `RCU Modes` 
        1. Add a graph `RSP Temperatures`
        1. Add a graph `TBB Status`
