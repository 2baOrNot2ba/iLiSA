===================================
`iLiSA`'s `monitor` `relay` package
===================================

What is it
----------

`iLiSA`'s `relay` package provides an agent that retransmitting the messages
sent out periodically from the `server` to LOFAR monitor clients. It is meant
to run on a gateway node between LOFAR/ASTRON and the local observatory
network.

Setting up
----------

To set up the relay agent, create a `relay_config.py` in the relay package
directory (edit e.g. `relay_config.py.template` and rename appropriately).
Rather than running the relay agent directly, it can make sense to set it up
so it starts during bootup, for instance, setting it up as a rc script.
To do this (on an Ubuntu-like linux) add in the file `/etc/rc.local` the line:

.. code-block:: sh

  sudo nohup /home/tobia/lofar/iStnMonitor/relay/stnMonitorRelay.py &

To manually start it, if it is not already running, execute:

.. code-block:: shell-session

  $ sudo sh /etc/rc.local

Info on the config parameters (found in `relay_config.py`):

* `RELAYSTNSTAT` Should Relay status message be copied as is and sent to
  another address? This option can be used if client cannot be reached from
  LCU.
* `IPto="123.45.67.89` IP address to which gateway should relay to,
  if `RELAYSTNSTAT` is True. This is the client's IP address.
* `SHAMECAST` Should the station status be shamecast?
* `isLogging` Should the last station monitoring status be logged? This can
  be useful for debugging, as it shows if the relay agent is receiving LCU
  status messages.
* `logfilename` filename with the latest station status message.
