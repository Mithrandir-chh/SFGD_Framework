## SuperFGD Data-Preprocessing Package
This package is used for the preprocessing of the SuperFGD data. It performs the following tasks:
1. Unpack the binary output of the DAQ to a _raw.root file. 
2. Perform the calibration to obtain the amplitude in units of p.e. and produce a _calib.root file.
3. Prepare the data for analysis by converting it to the Event structure (_events.root).


### Building the package
Assuming that you already have a ROOT installation and ROOTSYS is set.

Build the unpacking package from the /build/ folder by running
   ```
   $ cmake ../
   $ make
   ```

### Using the package
Navigate to the /bin/ folder in order to start running some codes.

```
./TDMunpack -f  /path/to/.daq
```
or 

```
for f in /path/to/daq/folder/*.daq; do ./TDMunpack -f "$f"; done
```

This will create separate .daq files for each FEB on the MCR. You must run this command on all 4 MCRs separately.


```
ls /path/to/data/folder/*Slot* > febs_files_list.list
```
This will create a list of all FEB .daq files for the unpacking code to go through.
   
     
```
./unpack 
```
This will unpack the data from all FEB files into a root file (_raw.root).


```
./Calibration  /path/to/_raw.root
```
This will prompt you to enter the HG and LG factors corresponding to the values set during the beam test for the run. Once you enter them, the code will use the calibration factors stored in /bin/calib_files/ for each setting and for each channel to calibrate the data And produce the _calib.root file.


```
./EventStructure  /path/to/_calib.root
```
This will create the _events.root file with the data in the Event structure. This file is used as the input for the analysis package.

If you are not using FEB 12, use the following variant of event structure.

```
./EventStructure_12free /path/to/_calib.root
```
This will create two TTrees in a single _event.root file: AllEvents and TimeGroupedEvents. AllEvents is simply grouped by spill. TimeGroupedEvents is grouped by event windows in spills. 


Quality of life scripts: 

**Unpack.sh**

'''
./Unpack.sh /path/to/run/folder1 /path/to/run/folder2
'''
Notice this command has a capital "U" instead of the lowercase "u" in ./unpack. 
This bash script will automatically run the entire data_preprocessing pipeline from TDMunpack to EventStructure_12free. 
This script should be operated on run folders, ie. one level above subrun folders and two levels above .daq files. This script can take in multiple folder paths. 

**delet_3D.sh**
'''
./delet_3D.sh
'''
This bash script will search for and delete files or folders that match specified names in the script. 
## This deletion is permanent, be extremely careful when specifying target names.





