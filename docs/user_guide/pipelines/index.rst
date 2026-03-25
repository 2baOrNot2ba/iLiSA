LOFAR station beamform streams (BFS): pipelines and backend processing
======================================================================

The ``pipelines`` package deals with the so-called *beamformed streams* (BFS), that is,
complex voltages generated when beamforming station-level array and sent out as UDP packets.
This package contain tools to record BFS as files on disk and process these files.

Recording BFS
-------------

BFS UDP packets can be setup to stream to a backend computer called
the *Data Recording Unit (DRU)*. ``iLiSA/pipelines``  provides two S/W recorders:

* ``dump_udp_ow`` (ID: ``ow``) is Olaf Wucknitz's c-code recorder
* ``rec_bf_stream_py`` (ID: ``py``) is ``iLiSA``'s alternative Python-based recorder

These recorders do not need to be run directly, but rather the ``pipelines``
package provides a unified wrapper called ``pl_rec`` that uses these two
recorders. When running ``pl_rec``

.. code-block:: console
    $ pl_rec -h
    usage: pl_rec [-h] [-m] [-t STARTTIME] [-p PORTS] [-b BFDATADIR] [-d DURATION]
                  [-f FILE_DURATION] [-w WHICH] [-r RCUMODE] [-s STNID] [-c] [-v]

    optional arguments:
      -h, --help            show this help message and exit
      -m, --mockrun         Run mock rec
      -t STARTTIME, --starttime STARTTIME
                            Start-time: (iso format) YYYY-mm-ddTHH:MM:SS
      -p PORTS, --ports PORTS
                            List of port number(s)
      -b BFDATADIR, --bfdatadir BFDATADIR
                            Template directory for BF data
      -d DURATION, --duration DURATION
                            Duration of recording in seconds
      -f FILE_DURATION, --file_duration FILE_DURATION
                            Duration of dumped files in seconds
      -w WHICH, --which WHICH
                            Which backend recorder: ow or py
      -r RCUMODE, --rcumode RCUMODE
                            rcumode or spectral window
      -s STNID, --stnid STNID
                            station id
      -c, --compress        Compress recorded data
      -v, --version         Print version of module

one can select which actual recorder is use (default is ``ow``).

BFS file processing
-------------------

The recorders in the previous section dump BFS packets in files on the DRU.
Typically this is produced by running an ``iLiSA`` observing session that
includes  ``rec`` with ``BFS`` option. These BFS files can be processed using
the CLI script ``pl_bfs``, which has several subcommands:

.. code-block:: console
    pl_bfs -h
    usage: pl_bfs [-h] {check,show,bin,bst,num,npy,meta,fix,info} ...

    positional arguments:
      {check,show,bin,bst,num,npy,meta,fix,info}
                            sub-command help
        check               Check BFS file
        show                Show packet from BFS file
        bin                 Convert BFS to binary
        bst                 Convert to BST files
        num                 NOT IMPLEMENTED YET (testing)
        npy                 Convert BFS to np npy
        meta                BFS metadata
        fix                 Strip BFS
        info                Info on BFS filefolder

The BFS contain only a limited set of metadata about the observation associated
with the files. Additional metadata can be contained in either the file headers
(*.h) or filefolder name, both of which are produced by ``iLiSA``.

If ``iLiSA`` was not used to create the BFS files, then either the headers or
the filefolder names will need to be created. The latter is probably easiest
as all is required is symbolic link (or renaming the folder containing the BFS
files) with the format

.. code-block:: console
    $ ln -s <BFS-folder> <iLiSA-filefolder>
    where <iLiSA-filefolder> :=
            <station_id(6c)>[<ANTSET>]_<filenametime>_spw<rcumodes>_sb<subbands>\
            _int<integration>_dur<duration_scan>[_dir<pointing>][_cal|_mod]\
            _<ldat_type>
          where <ANTSET> := <LBA|HBA%<conf1>[%<conf2>]>

It is currently also necessary to uncompress any compressed BSF files and
it is recommend to run ``fix`` on the BFS files before conversion subcommands
to ``pl_bfs`` are use.

bst (subcommand to pl_bfs)
^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``pl_bfs bst`` converts the BFS files into beamform statistics BST-like files.
Normal BSTs are *integer* second integrated BFS of the co-pol channels XX and YY.
The BST-like files created here extend this to sub-second integration and adds the
cross-pol channel XY.

.. code-block:: console
    usage: pl_bfs bst [-h] bfs_ff integration

    positional arguments:
      bfs_ff       BFS filefolder path
      integration  Integration in float seconds

    optional arguments:
      -h, --help   show this help message and exit

