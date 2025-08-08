from __future__ import print_function
import os
import sys
import ntpath

#To Run this is script it is necesary first to generate 2 files:
# A) a list with all the folders containing the data that want to be unpacked
# B) a list with a map of folder names and HG/LG. The structure must be [folderName Id IsValidFlag HG LG]. Where IdValidFlag is 1 if the folder wants to be upnacked and 0 otherwise.
#
# 
# Assuming A) is stored in list1 and B) in list2 one can generate .sh files with insturctions to unpack all the folders by running '$ python unpack.py list1 list2' 
# Then one can send jobs using each .sh file to unpack each folder as an individual job.

if len(sys.argv) < 2:
    raise SyntaxError("The script needs an external file containing a target list of directories\n\
                      Example: python unpacking.py listOffolders.list\n The directiories must contain the full path.")

HG_dict = {}
LG_dict = {}

with open(sys.argv[1]) as mapping:
    line = 'line'
    while line:
        line = mapping.readline().rstrip()
        cols = line.split(" ")
        if len(cols) < 5: break
        if int(cols[2]) == 1:
            HG_dict[cols[0]] = int(cols[3])
            LG_dict[cols[0]] = int(cols[4])

job_cnt = 0
with open(sys.argv[2]) as f:
    path = 'path'
    while path:
        path = f.readline().rstrip()
        if path:
            print(path)
            folder = ntpath.basename(path)
            job_name = str(job_cnt) + '.sh'
            job_cnt += 1
            with open(job_name, 'w') as job:
                print('#!/bin/bash\n', file=job)
                print('source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.14.04/x86_64-centos7-gcc48-opt/root/bin/thisroot.sh', file=job)
                print('cd ', path, file=job)
                print('rm -r sfgd_framework', file=job)
                print('cp -r /afs/cern.ch/user/c/cjesus/sfgd_framework .', file=job)
                print('rm *Slot*', file=job)
                print('cd sfgd_framework', file=job)
                print('cd data_preprocessing/build/', file=job)
                print('/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.8.1/Linux-x86_64/bin/cmake ../', file=job) 
                print('make', file=job)
                print('cd ../bin/', file=job)
                print('./TDMunpack -f ../../../*MCR0**.daq*', file=job)
                print('./TDMunpack -f ../../../*MCR1**.daq*', file=job)
                print('./TDMunpack -f ../../../*MCR2**.daq*', file=job)
                print('./TDMunpack -f ../../../*MCR3**.daq*', file=job)
                print('ls ../../../*Slot* > febs_files_list.list', file=job)
                print('./unpack -f d', file=job)
                print('./Calibration  ../../../*_raw.root', HG_dict[folder], LG_dict[folder], file=job)
                print('./EventStructure ../../../*_calib.root', file=job)
                print('rm *Slot*', file=job)
