Installing
==========
Some notes on installing ``iLiSA``.

Development
-----------
If you intend to develop the ``iLiSA`` code, then after downloading the latest version
using ``git pull`` do the following:

.. code-block:: console

   [iLiSA] python3 setup.py devel --user

Pipeline Only
-------------
Some nodes are intended as a pure data recording unit (DRU) and only
needs the ``pipelines`` module. In other words the node only records
station data in realtime and possibly processes it in quasi-realtime with
a pipeline process. In this case only the ``pipelines`` module needs installation,
so here are the instructions.

After downloading

.. code-block:: console

   [iLiSA] python3 setup.py install --user
