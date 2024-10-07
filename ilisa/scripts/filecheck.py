#!/usr/bin/python3
"""
List size of all lofar datafiles found under the directory passed as argument

This is an example of a iLiSA postprocessing script.
It lives under ~/.config/ilisa/postprocessing/ and can be executed after a
completed ilisa_sched run which will pass it the iLiSA session folder path.

So move this script to <userhome>/.config/ilisa/postprocessing/
and add it to an ilisa_sched yml file to add a postprocessing step.
"""
import sys
import os
sessdir = sys.argv[1]

ldatdirs = {'acc': [], 'bfs': [], 'sst': [], 'xst': []}
ldatfiles = {'acc': {}, 'bfs': {}, 'sst': {}, 'xst': {}}
ls_sessdir = os.listdir(sessdir)
scn_dirs = filter(lambda s: s.startswith('scan_'), os.listdir(sessdir))
for scn_dir in scn_dirs:
    datadirs = os.listdir(os.path.join(sessdir,scn_dir))
    for k in ldatdirs:
        _ = filter(lambda s: s.endswith('_'+k), datadirs)
        ldatdirs[k].extend(list(map(lambda d: os.path.join(scn_dir, d), _)))
for ldatdirtype in ldatdirs:
    for ldatdir in ldatdirs[ldatdirtype]:
        ldatdir_path = os.path.join(sessdir, ldatdir)
        datafiles = os.listdir(ldatdir_path)
        if ldatdirtype == 'bfs':
            datafiles = filter(lambda f: f.endswith('.zst'), datafiles)
        else:
            datafiles = filter(lambda f: f.endswith('.dat'), datafiles)
        ldatfiles[ldatdirtype][ldatdir] = {}
        for f in datafiles:
            datafilepath = os.path.join(ldatdir, f)
            sz = os.path.getsize(os.path.join(sessdir, datafilepath))
            ldatfiles[ldatdirtype][ldatdir][f] = sz
outputf = os.path.join(os.path.dirname(__file__), "out.txt")
with open(outputf, 'w') as f:
    f.write(sessdir+'\n')
    for ldatdirtype in ldatdirs:
        f.write(ldatdirtype+':\n')
        for ldd in ldatdirs[ldatdirtype]:
            f.write('  '+ldd+':\n')
            for df in ldatfiles[ldatdirtype][ldd]:
                sz = ldatfiles[ldatdirtype][ldd][df]
                f.write('    '+df+': '+str(float(sz/1024**2))+'MiB\n')
