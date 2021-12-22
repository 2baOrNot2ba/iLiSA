#!/usr/bin/python

## Script for parsing output of tbbctl
## A. Horneffer 30. Jan 2015

import re

def parse_tbbctl_status(status_data):
    seperator_line = "---  -----  ------  -----  -----  -----  ----  ----  -------  -------  -------  -------  -----"
    line_pattern = re.compile(r" [ \d]\d ")
    tbb_pattern = re.compile(r"( [ \d]\d) +\d +ETH +([\d\.]+)V +([\d\.]+)V +([\d\.]+)V +([\d]+)'C +([\d]+)'C +([\d]+)'C +\w\w +([\d]+)'C +\w\w +([\d]+)'C +\w\w +([\d]+)'C +\w+")
    status_lines = status_data.splitlines()
    tbbstatus = {}
    tbbstatus['badTBBs'] = 0
    tbbstatus['volt12min'] = 999.
    tbbstatus['volt12mean'] = 0.
    tbbstatus['volt12max'] = 0.
    tbbstatus['volt25min'] = 999.
    tbbstatus['volt25mean'] = 0.
    tbbstatus['volt25max'] = 0.
    tbbstatus['volt33min'] = 999.
    tbbstatus['volt33mean'] = 0.
    tbbstatus['volt33max'] = 0.
    tbbstatus['PCBtempmin'] = 999.
    tbbstatus['PCBtempmean'] = 0.
    tbbstatus['PCBtempmax'] = 0.
    tbbstatus['TPtempmin'] = 999.
    tbbstatus['TPtempmean'] = 0.
    tbbstatus['TPtempmax'] = 0.
    tbbstatus['MPtempmin'] = 999.
    tbbstatus['MPtempmean'] = 0.
    tbbstatus['MPtempmax'] = 0.
    tbbstatus['goodlines'] = 0
    seperator_found = False
    for line in status_lines:
        if not seperator_found: 
            if (line != seperator_line):
                seperator_found = True
            continue
        linematch = line_pattern.match(line)
        if (linematch):
            tmatch = tbb_pattern.match(line)        
            if (tmatch):
                volt12 = float(tmatch.group(2))
                volt25 = float(tmatch.group(3))
                volt33 = float(tmatch.group(4))
                PCBtemp = float(tmatch.group(5))
                TPtemp = float(tmatch.group(6))
                MPtemps = [float(v) for v in tmatch.groups()[6:10]]
                tbbstatus['volt12min'] = min(tbbstatus['volt12min'],volt12)
                tbbstatus['volt12max'] = max(tbbstatus['volt12max'],volt12)
                tbbstatus['volt12mean'] += volt12
                tbbstatus['volt25min'] = min(tbbstatus['volt25min'],volt25)
                tbbstatus['volt25max'] = max(tbbstatus['volt25max'],volt25)
                tbbstatus['volt25mean'] += volt25
                tbbstatus['volt33min'] = min(tbbstatus['volt33min'],volt33)
                tbbstatus['volt33max'] = max(tbbstatus['volt33max'],volt33)
                tbbstatus['volt33mean'] += volt33
                tbbstatus['PCBtempmin'] = min(tbbstatus['PCBtempmin'],PCBtemp)
                tbbstatus['PCBtempmax'] = max(tbbstatus['PCBtempmax'],PCBtemp)
                tbbstatus['PCBtempmean'] += PCBtemp
                tbbstatus['TPtempmin'] = min(tbbstatus['TPtempmin'],TPtemp)
                tbbstatus['TPtempmax'] = max(tbbstatus['TPtempmax'],TPtemp)
                tbbstatus['TPtempmean'] += TPtemp
                tbbstatus['MPtempmin']=min(tbbstatus['MPtempmin'],min(MPtemps))
                tbbstatus['MPtempmax']=max(tbbstatus['MPtempmax'],max(MPtemps))
                tbbstatus['MPtempmean'] += sum(MPtemps)
                tbbstatus['goodlines'] += 1
            else:
               tbbstatus['badTBBs'] += 1 
    if tbbstatus['goodlines'] > 0:
        tbbstatus['volt12mean'] /= tbbstatus['goodlines']
        tbbstatus['volt25mean'] /= tbbstatus['goodlines']
        tbbstatus['volt33mean'] /= tbbstatus['goodlines']
        tbbstatus['PCBtempmean'] /= tbbstatus['goodlines']
        tbbstatus['TPtempmean'] /= tbbstatus['goodlines']
        tbbstatus['MPtempmean'] /= (tbbstatus['goodlines']*4.)
    else:
        tbbstatus['volt12mean'] = -1.
        tbbstatus['volt25mean'] = -1.
        tbbstatus['volt33mean'] = -1.
        tbbstatus['PCBtempmean'] = -1.
        tbbstatus['TPtempmean'] = -1.
        tbbstatus['MPtempmean'] = -1.
    return tbbstatus


if __name__ == "__main__":
    import sys
    print("usage: python parse_tbbctl.py <tbbctl-status-out>")
    statusname = sys.argv[1]
    infile = open(statusname,'r')
    statusdata = infile.read()
    infile.close()
    tbbstatus = parse_tbbctl_status(statusdata)
    print("TBB-Status:")
    print(tbbstatus)
    
