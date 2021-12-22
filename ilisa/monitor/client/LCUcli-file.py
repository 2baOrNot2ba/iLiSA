#!/usr/bin/python
#LOFAR station monitoring UDP client v1.1
#Receives updates in UDP packets and enters these into RRD files using rrdtool
#
#Simon Casey
#Onsala Space Observatory
#
#
#Revision history
#v1.0	2012-10-12	Initial release
#v1.1	2013-04-03	Added '-t' to rrdupdate to ensure correct field order
#v1.2	2014-??-??	Added handling of RSP data
#v1.3	2015-02-02	Added partial handling of TBB data 

import socket, rrdtool, re, sys

import time, calendar
from optparse import OptionParser

#==================================================
# Adjust settings and file names / pathes below
#==================================================

STN_string="de601ec"

RRD_FILE_env="/home/software/iStnMonitor-data/de601_environment_20130408.rrd"
RRD_FILE_status="/home/software/iStnMonitor-data/de601_status_20130408.rrd"
RRD_FILE_users="/home/software/iStnMonitor-data/de601_users_20130619.rrd"
RRD_FILE_RCUModes="/home/software/iStnMonitor-data/de601_RCUModes_201312003.rrd"
RRD_FILE_MinTemp="/home/software/iStnMonitor-data/de601_RSPmintemp_20131203.rrd"
RRD_FILE_MinVolt="/home/software/iStnMonitor-data/de601_RSPminvolt_20131203.rrd"
RRD_FILE_MeanTemp="/home/software/iStnMonitor-data/de601_RSPmeantemp_20131203.rrd"
RRD_FILE_MeanVolt="/home/software/iStnMonitor-data/de601_RSPmeanvolt_20131203.rrd"
RRD_FILE_MaxTemp="/home/software/iStnMonitor-data/de601_RSPmaxtemp_20131203.rrd"
RRD_FILE_MaxVolt="/home/software/iStnMonitor-data/de601_RSPmaxvlot_20131203.rrd"
RRD_FILE_TBBs="/home/software/iStnMonitor-data/de601_TBBgood_20131203.rrd"
RRD_FILE_TBBMinTemp="/home/software/iStnMonitor-data/de601_TBBmintemp_20131203.rrd"
RRD_FILE_TBBMeanTemp="/home/software/iStnMonitor-data/de601_TBBmeantemp_20131203.rrd"
RRD_FILE_TBBMaxTemp="/home/software/iStnMonitor-data/de601_TBBmaxtemp_20131203.rrd"

#==================================================
# Adjust settings and file names / pathes above
#==================================================


ENV_order = "Temperature:Humidity"
STATUS_order = "Heater:48V:LCU:Lightning:SWLevel:Switch"
USERS_order = "all_users:local_users"
RCU_order = "RCUFail:RCU0:RCU1:RCU2:RCU3:RCU4:RCU5:RCU6:RCU7"
RSPVolt_order = "V12:V25:V33"
RSPTemp_order = "PCBTemp:BPTemp:APTemp"
TBB_order = "BadTBBs:GoodTBBs"
TBBTemp_order = "PCBTemp:TPTemp:MPTemp"

parser = OptionParser();

parser.add_option("-F", "--file", dest="File",
                  help="Logfile from LCUcli.py to read data from.")
parser.add_option("-U", "--update-all", dest="Update",action="store_true",
                  help="Actually update all the Databases")
parser.add_option("--ue", "--update-env", dest="UpEnv",action="store_true",
                  help="Actually update the Environment Database")
parser.add_option("--us", "--update-status",dest="UpStatus",action="store_true",
                  help="Actually update the LCU-Status Database")
parser.add_option("--uu", "--update-users", dest="UpUsers",action="store_true",
                  help="Actually update the Users Database")
parser.add_option("--ur", "--update-rcus", dest="UpRCU",action="store_true",
                  help="Actually update the RCU-Modes Database")
parser.add_option("--update-RSP", dest="UpRSP",action="store_true",
                  help="Actually update the RSP-Hardware-Status Database")
parser.add_option("--update-TBB", dest="UpTBB",action="store_true",
                  help="Actually update the TBB-Hardware-Status Database")

(options, args) = parser.parse_args()

if not options.File:
    print "Please specify a file to read from"
    exit(1)

updateEnv = False
updateStat = False
updateUser = False
updateRCUs = False
updateRSPVolt = False
updateRSPTemp = False
updateTBB = False

if options.Update:
    updateEnv = True
    updateStat = True
    updateUser = True
    updateRCUs = True
    updateRSPVolt = True
    updateRSPTemp = True
    updateTBB = True
if options.UpEnv:
    updateEnv = True
if options.UpStatus:
    updateStat = True
if options.UpUsers:
    updateUser = True
if options.UpRCU:
    updateRCUs = True
if options.UpRSP:
    updateRSPVolt = True
    updateRSPTemp = True
if options.UpTBB:
    updateTBB = True

infile = open(options.File,'r')

if not infile:
    print "Failed to open file:",options.File
    exit(1)

while True:
    message = ""
    line = ""
    line = infile.readline()
    #exit if the last read "line" was empty
    if len(line)<1:
        exit(0)
    while (line and not line.startswith('LOFAR_STN_STATUS (version):')):
        message = message + line
        line = infile.readline()
    if (not message):
        continue
    if re.search('Station: *(.+?)[,\n]', message).group(1) == STN_string or re.search('Station: *(.+?)[,\n]', message).group(1) == 'Unknown':
        timestring = message[0:19]
        timesecs = calendar.timegm(time.strptime(timestring,'%Y-%m-%d-%H:%M:%S'))
        datestamp = str(timesecs)

        L_Temp=re.search('Cab3 Temp: *(.+?)C,', message).group(1)
        L_Hum=re.search('Cab3 Hum: *(.+?)%,', message).group(1)
        if re.search('Heater: *(.+?)[,\n]', message).group(1) == 'ON':
            L_Heater="2"
        else: L_Heater="0"
        if re.search('48V: *(.+?)[,\n]', message).group(1) == 'ON':
            L_48V="1"
        else: L_48V="0"
        if re.search('LCU: *(.+?)[,\n]', message).group(1) == 'ON':
            L_LCU="3"
        else: L_LCU="0"
        if re.search('Lightning: *(.+)[,\n]', message).group(1) == 'N.A.':
            L_Lightning="0"
        else: L_Lightning="4"
        
        if re.search('Switch: *(.+?)[,\n]', message).group(1) == 'ilt':
            L_Switch="1"
        else: L_Switch="0"
        L_SwL = re.search('Software Level: *(.+?)[,\n]', message).group(1)

        if (re.search('All Users: *(.+?)', message)):
            User_LCUall = re.search('All Users: *(.+?)', message).group(1)
            User_LCUlocal = re.search('Local Users: *(.+?)', message).group(1)
        else:
            User_LCUall = "-1"
            User_LCUlocal = "-1"

        if (re.search('RCUmodes', message)):
            RCUmatch = re.search('RCUmodes -1:(\d+) 0:(\d+) 1:(\d+) 2:(\d+) 3:(\d+) 4:(\d+) 5:(\d+) 6:(\d+) 7:(\d+)', message)
        else:
            RCUmatch = False
        if (re.search('RSPvolt ', message)):
            RSPVoltMeanMatch = re.search('RSPvolt V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
            RSPVoltMinMatch = re.search('RSPvoltMin V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
            RSPVoltMaxMatch = re.search('RSPvoltMax V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
        else:
            RSPVoltMeanMatch = False
            RSPVoltMinMatch = False
            RSPVoltMaxMatch = False
        if (re.search('RSPtemps PCBmean', message)):
            RSPTempMeanMatch = re.search('RSPtemps PCBmean: (\d+\.\d\d)C, BPmean: (\d+\.\d\d)C, APmean: (\d+\.\d\d)C', message)
            RSPTempMinMatch = re.search('RSPtemps PCBmin: (\d+\.\d\d)C, BPmin: (\d+\.\d\d)C, APmin: (\d+\.\d\d)C', message)
            RSPTempMaxMatch = re.search('RSPtemps PCBmax: (\d+\.\d\d)C, BPmax: (\d+\.\d\d)C, APmax: (\d+\.\d\d)C', message)
        else:
            RSPTempMeanMatch = False
            RSPTempMinMatch = False
            RSPTempMaxMatch = False

            # Collect TBB Data
        if (re.search('Bad-TBBs', message)):
            TBBgoodmatch = re.search('Bad-TBBs: (\d+), Good-TBBs: (\d+)', message)
            TBBTempMeanMatch = re.search('TBBtemps PCBmean: (\d+\.\d\d)C, TPmean: (\d+\.\d\d)C, MPmean: (\d+\.\d\d)C', message)
            TBBTempMinMatch = re.search('TBBtemps PCBmin: (\d+\.\d\d)C, TPmin: (\d+\.\d\d)C, MPmin: (\d+\.\d\d)C', message)
            TBBTempMaxMatch = re.search('TBBtemps PCBmax: (\d+\.\d\d)C, TPmax: (\d+\.\d\d)C, MPmax: (\d+\.\d\d)C', message)
            TBBVoltMeanMatch = re.search('TBBvolt V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
            TBBVoltMinMatch = re.search('TBBvoltMin V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
            TBBVoltMaxMatch = re.search('TBBvoltMax V1.2: (\d+\.\d\d), V2.5: (\d+\.\d\d), V3.3: (\d+\.\d\d)', message)
        else:
            TBBgoodmatch = False
            TBBTempMeanMatch = False
            TBBTempMinMatch = False
            TBBTempMaxMatch = False
            TBBVoltMeanMatch = False
            TBBVoltMinMatch = False
            TBBVoltMaxMatch = False
    
        sys.stdout.flush()
        
        outstring = "Update for: " + timestring + " -- "
        if updateEnv:
            try:
                rrdtool.update(RRD_FILE_env, '-t', ENV_order, datestamp + ':' + L_Temp + ':' + L_Hum)
            except:
                print "Error updaing " + RRD_FILE_env
            outstring += "Updated Env, "
        else:
            outstring += "Would update Env: " + ENV_order + " --- " + datestamp + ':' + L_Temp + ':' + L_Hum + "\n"

        if updateStat:
            try:
                rrdtool.update(RRD_FILE_status, '-t', STATUS_order, datestamp + ':' + L_Heater + ':' + L_48V + ':' + L_LCU + ':' + L_Lightning + ':' + L_SwL + ':' + L_Switch)
            except:
                print "Error updating " + RRD_FILE_status
            outstring += "Updated Stat, "
        else:
            outstring += "Would update Stat: " + STATUS_order + " --- " + datestamp + ':' + L_Heater + ':' + L_48V + ':' + L_LCU + ':' + L_Lightning + ':' + L_SwL + ':' + L_Switch + "\n"

        if updateUser:
            try:
                rrdtool.update(RRD_FILE_users, '-t', USERS_order, datestamp + ':' + User_LCUall + ':' + User_LCUlocal)
            except:
                print "Error updaing " + RRD_FILE_users
            outstring += "Updated Users "
        else:
            outstring += "Would update Users: " + USERS_order + " --- " + datestamp + ':' + User_LCUall + ':' + User_LCUlocal + '\n' 

        if (RCUmatch):
            if updateRCUs:
                try:
                    rrdtool.update(RRD_FILE_RCUModes, '-t', RCU_order, datestamp + ':' +  RCUmatch.group(1) + ':' +  RCUmatch.group(2) + ':' +  RCUmatch.group(3) + ':' +  RCUmatch.group(4) + ':' +  RCUmatch.group(5) + ':' +  RCUmatch.group(6) + ':' +  RCUmatch.group(7) + ':' +  RCUmatch.group(8) + ':' +  RCUmatch.group(9) )
                except:
                    print "Error updaing " + RRD_FILE_RCUModes
                outstring += "Updated RCUs "
            else:
                outstring +="Would update RCUMode: " + RCU_order + " --- " + datestamp + ':' +  RCUmatch.group(1) + ':' +  RCUmatch.group(2) + ':' +  RCUmatch.group(3) + ':' +  RCUmatch.group(4) + ':' +  RCUmatch.group(5) + ':' +  RCUmatch.group(6) + ':' +  RCUmatch.group(7) + ':' +  RCUmatch.group(8) + ':' +  RCUmatch.group(9) + '\n' 

        if (RSPVoltMeanMatch and RSPVoltMinMatch and RSPVoltMaxMatch):
            if updateRSPVolt:
                try:
                    rrdtool.update(RRD_FILE_MinVolt, '-t', RSPVolt_order, datestamp + ':' + RSPVoltMinMatch.group(1) + ':' + RSPVoltMinMatch.group(2) + ':' + RSPVoltMinMatch.group(3) ) 
                except:
                    print "Error updaing " + RRD_FILE_MinVolt
                try:
                    rrdtool.update(RRD_FILE_MeanVolt, '-t', RSPVolt_order, datestamp + ':' + RSPVoltMeanMatch.group(1) + ':' + RSPVoltMeanMatch.group(2) + ':' + RSPVoltMeanMatch.group(3) ) 
                except:
                    print "Error updaing " + RRD_FILE_MeanVolt
                try:
                    rrdtool.update(RRD_FILE_MaxVolt, '-t', RSPVolt_order, datestamp + ':' + RSPVoltMaxMatch.group(1) + ':' + RSPVoltMaxMatch.group(2) + ':' + RSPVoltMaxMatch.group(3) ) 
                except:
                    print "Error updaing " + RRD_FILE_MaxVolt
                outstring += "Updated RSPVolt "
            else:
                outstring +="Would update RSPVoltMean: " + RSPVolt_order + " --- " + datestamp + ':' + RSPVoltMeanMatch.group(1) + ':' + RSPVoltMeanMatch.group(2) + ':' + RSPVoltMeanMatch.group(3) + '\n' 
                outstring +="Would update RSPVoltMin: " + RSPVolt_order + " --- " + datestamp + ':' + RSPVoltMinMatch.group(1) + ':' + RSPVoltMinMatch.group(2) + ':' + RSPVoltMinMatch.group(3) + '\n' 
                outstring +="Would update RSPVoltMax: " + RSPVolt_order + " --- " + datestamp + ':' + RSPVoltMaxMatch.group(1) + ':' + RSPVoltMaxMatch.group(2) + ':' + RSPVoltMaxMatch.group(3) + '\n' 

        if (RSPTempMaxMatch and RSPTempMeanMatch and RSPTempMinMatch):
            if updateRSPTemp:
                try:
                    rrdtool.update(RRD_FILE_MaxTemp, '-t', RSPTemp_order, datestamp + ':' + RSPTempMaxMatch.group(1) + ':' + RSPTempMaxMatch.group(2) + ':'+ RSPTempMaxMatch.group(3) )
                except:
                    print "Error updaing " + RRD_FILE_MaxTemp
                try:
                    rrdtool.update(RRD_FILE_MeanTemp, '-t', RSPTemp_order, datestamp + ':' + RSPTempMeanMatch.group(1) + ':' + RSPTempMeanMatch.group(2) + ':'+ RSPTempMeanMatch.group(3) )
                except:
                    print "Error updaing " + RRD_FILE_MeanTemp
                try:
                    rrdtool.update(RRD_FILE_MinTemp, '-t', RSPTemp_order, datestamp + ':' + RSPTempMinMatch.group(1) + ':' + RSPTempMinMatch.group(2) + ':'+ RSPTempMinMatch.group(3) )
                except:
                    print "Error updaing " + RRD_FILE_MinTemp
                outstring += "Updated RSPTemp "
            else:
                outstring +="Would update RSPTempMean: " + RSPTemp_order + " --- " + datestamp + ':' + RSPTempMeanMatch.group(1) + ':' + RSPTempMeanMatch.group(2) + ':'+ RSPTempMeanMatch.group(3) + '\n' 
                outstring +="Would update RSPTempMin: " + RSPTemp_order + " --- " + datestamp + ':' + RSPTempMinMatch.group(1) + ':' + RSPTempMinMatch.group(2) + ':'+ RSPTempMinMatch.group(3) + '\n' 
                outstring +="Would update RSPTempMax: " + RSPTemp_order + " --- " + datestamp + ':' + RSPTempMaxMatch.group(1) + ':' + RSPTempMaxMatch.group(2) + ':'+ RSPTempMaxMatch.group(3) + '\n' 
        
        if (TBBgoodmatch and TBBTempMeanMatch and RSPTempMinMatch and TBBTempMaxMatch):
            if updateTBB:
                try:
                    rrdtool.update(RRD_FILE_TBBs, '-t', TBB_order, datestamp + ':' + TBBgoodmatch.group(1) + ':' + TBBgoodmatch.group(2))
                except:
                    print "Error updaing " + RRD_FILE_TBBs
                try:
                    rrdtool.update(RRD_FILE_TBBMeanTemp, '-t', TBBTemp_order, datestamp + ':' + TBBTempMeanMatch.group(1) + ':' + TBBTempMeanMatch.group(2) + ':' + TBBTempMeanMatch.group(3))
                except:
                    print "Error updaing " + RRD_FILE_TBBMeanTemp
                try:
                    rrdtool.update(RRD_FILE_TBBMinTemp, '-t', TBBTemp_order,  datestamp+ ':' + TBBTempMinMatch.group(1) + ':' + TBBTempMinMatch.group(2) + ':' + TBBTempMinMatch.group(3))
                except:
                    print "Error updaing " + RRD_FILE_TBBMinTemp
                try:
                    rrdtool.update(RRD_FILE_TBBMaxTemp, '-t', TBBTemp_order, datestamp + ':' + TBBTempMaxMatch.group(1) + ':' + TBBTempMaxMatch.group(2) + ':' + TBBTempMaxMatch.group(3))
                except:
                    print "Error updaing " + RRD_FILE_TBBMaxTemp
                outstring += "Updated TBB "
            else:
                outstring +="Would update TBBgood: " + TBBTempMeanMatch.group(1) + ':' + TBBTempMeanMatch.group(2) + ':' + TBBTempMeanMatch.group(3)+"\n"
                outstring +="Would update TBBMeanTemp: " + TBBTempMeanMatch.group(1) + ':' + TBBTempMeanMatch.group(2) + ':' + TBBTempMeanMatch.group(3)+"\n"
                outstring +="Would update TBBMinTemp: " + TBBTempMinMatch.group(1) + ':' + TBBTempMinMatch.group(2) + ':' + TBBTempMinMatch.group(3)+"\n"
                outstring +="Would update TBBMaxTemp: " + TBBTempMaxMatch.group(1) + ':' + TBBTempMaxMatch.group(2) + ':' + TBBTempMaxMatch.group(3)+"\n"
        outstring += '\n'
        print outstring
            
        sys.stdout.flush()
    else:
        print "Error parsing message:",message
