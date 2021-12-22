sed -i.bak '/\tGRAPH_ITEM_TYPE_LINE3/a\\tGRAPH_ITEM_TYPE_LINE4\t=> "LINE4",\n\tGRAPH_ITEM_TYPE_LINE5\t=> "LINE5",' include/global_arrays.php
sed -i.bak '/GRAPH_ITEM_TYPE_LINE3/a\define("GRAPH_ITEM_TYPE_LINE4",  11);\ndefine("GRAPH_ITEM_TYPE_LINE5",  12);' include/global_constants.php
sed -i.bak 's/LINE\[123\]/LINE\[12345\]/g' graphs.php
sed -i.bak 's/LINE\[123\]/LINE\[12345\]/g' lib/graph_variables.php
sed -i.bak 's/LINE\[123\]/LINE\[12345\]/g' lib/html.php
sed -i.bak 's/LINE\[123\]/LINE\[12345\]/g' lib/rrd.php

