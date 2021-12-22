#!/usr/bin/python

## Script for parsing output of rspctl
## A. Horneffer 3. Dec 2013

import re

def parse_rspctl_status(status_data):
    status_lines = status_data.splitlines()
    voltage_pattern = re.compile(r'RSP\[([ \d]+)\] 1.2 V:([\d\. ]+), 2.5 V:([\d\. ]+), 3.3 V:([\d\. ]+)')
    temp_pattern = re.compile(r'RSP\[([ \d]+)\] PCB_temp:([\d\. ]+), BP_temp:([\d\. ]+), Temp AP0:([\d\. ]+), AP1:([\d\. ]+), AP2:([\d\. ]+), AP3:([\d\. ]+)')
    rspstatus = {}
    rspstatus['volt12min'] = 999.
    rspstatus['volt12mean'] = 0.
    rspstatus['volt12max'] = 0.
    rspstatus['volt25min'] = 999.
    rspstatus['volt25mean'] = 0.
    rspstatus['volt25max'] = 0.
    rspstatus['volt33min'] = 999.
    rspstatus['volt33mean'] = 0.
    rspstatus['volt33max'] = 0.
    rspstatus['voltlines'] = 0
    rspstatus['PCBtempmin'] = 999.
    rspstatus['PCBtempmean'] = 0.
    rspstatus['PCBtempmax'] = 0.
    rspstatus['BPtempmin'] = 999.
    rspstatus['BPtempmean'] = 0.
    rspstatus['BPtempmax'] = 0.
    rspstatus['APtempmin'] = 999.
    rspstatus['APtempmean'] = 0.
    rspstatus['APtempmax'] = 0.
    rspstatus['templines'] = 0
    for line in status_lines:
        vmatch = voltage_pattern.match(line)
        if (vmatch):
            volt12 = float(vmatch.group(2))
            volt25 = float(vmatch.group(3))
            volt33 = float(vmatch.group(4))
            rspstatus['volt12min'] = min(rspstatus['volt12min'],volt12)
            rspstatus['volt12max'] = max(rspstatus['volt12max'],volt12)
            rspstatus['volt12mean'] += volt12
            rspstatus['volt25min'] = min(rspstatus['volt25min'],volt25)
            rspstatus['volt25max'] = max(rspstatus['volt25max'],volt25)
            rspstatus['volt25mean'] += volt25
            rspstatus['volt33min'] = min(rspstatus['volt33min'],volt33)
            rspstatus['volt33max'] = max(rspstatus['volt33max'],volt33)
            rspstatus['volt33mean'] += volt33
            rspstatus['voltlines'] += 1
        else:
            tmatch = temp_pattern.match(line)
            if (tmatch):
                PCBtemp = float(tmatch.group(2))
                BPtemp = float(tmatch.group(3))
                APtemps = [float(v) for v in tmatch.groups()[3:7]]
                rspstatus['PCBtempmin'] = min(rspstatus['PCBtempmin'],PCBtemp)
                rspstatus['PCBtempmax'] = max(rspstatus['PCBtempmax'],PCBtemp)
                rspstatus['PCBtempmean'] += PCBtemp
                rspstatus['BPtempmin'] = min(rspstatus['BPtempmin'],BPtemp)
                rspstatus['BPtempmax'] = max(rspstatus['BPtempmax'],BPtemp)
                rspstatus['BPtempmean'] += BPtemp
                rspstatus['APtempmin']=min(rspstatus['APtempmin'],min(APtemps))
                rspstatus['APtempmax']=max(rspstatus['APtempmax'],max(APtemps))
                rspstatus['APtempmean'] += sum(APtemps)
                rspstatus['templines'] += 1
    if rspstatus['voltlines'] >0:
        rspstatus['volt12mean'] /= rspstatus['voltlines']
        rspstatus['volt25mean'] /= rspstatus['voltlines']
        rspstatus['volt33mean'] /= rspstatus['voltlines']
    else:
        rspstatus['volt12mean'] = -1.
        rspstatus['volt25mean'] = -1.
        rspstatus['volt33mean'] = -1.
    if rspstatus['templines'] > 0:
        rspstatus['PCBtempmean'] /= rspstatus['templines']
        rspstatus['BPtempmean'] /= rspstatus['templines']
        rspstatus['APtempmean'] /= (rspstatus['templines']*4.)
    else:
        rspstatus['PCBtempmean'] = -1.
        rspstatus['BPtempmean'] = -1. 
        rspstatus['APtempmean'] = -1.
    return rspstatus

def parse_rspctl_rcu(rcu_data):
    rcu_lines = rcu_data.splitlines()
    rcu_pattern = re.compile(r'RCU\[([ \d]+)\]\.control=0x[\dabcdef]+ => ([ ONF]+), mode:([\d-]+), delay=')
    rcumodes = {'-1':0, '0':0, '1':0, '2':0, '3':0, '4':0, '5':0, '6':0, '7':0}
    for line in rcu_lines:
        rcumatch = rcu_pattern.match(line)
        if (rcumatch):
            rcumodes[rcumatch.group(3)]+=1
    return rcumodes


if __name__ == "__main__":
    import sys
    print("usage: python parse_rspctl.py <rspctl-status-out> <rspctl-rcu-out>")
    statusname = sys.argv[1]
    infile = open(statusname,'r')
    statusdata = infile.read()
    infile.close()
    rspstatus = parse_rspctl_status(statusdata)
    print("RSP-Status:")
    print(rspstatus)
    
    rcuname = sys.argv[2]
    infile = open(rcuname,'r')
    rcudata = infile.read()
    infile.close()
    rcumodes = parse_rspctl_rcu(rcudata)
    print("RCU-Modes:")
    print(rcumodes)
