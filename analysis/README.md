## SuperFGD Analysis Package
This package is used for the analysis of pre-processed SuperFGD data. It performs the following tasks:
1. Takes in _events.root file from data_preprocess 
2. Perform event filtering and plots 3 view histograms for passed events
3. Reconstruct 3D voxels and fit a trejectory to through passed events

### Building the package
Assuming that you already have a ROOT installation and ROOTSYS is set.

Build the unpacking package from the /build/ folder by running
   ```
   $ cmake ../
   $ make
   ```

### Using the package
Navigate to the /bin/che/ folder in order to start running some codes.

In /bin/che/ folder, here are the important functions that one could use for quick data analysis: 
1. EventDisplayLite
    **purpose**: Filters events by their high-pe counts and display the 3-view histogram for each passed event. 
    **input**: _events.root
    **output**: a folder "EventDisplay_PNG_Filtered" inside of which are 3 view histograms of each high-count event + a _events_display.root with all histograms. 
2. 3D_EventTrajectorySide
    **purpose**: The vanilla reconstruction and fitting function. Reconstructs voxels from hits, apply calibration, and does PCA fitting. Can also do event filtering based on track linearity, planeness, and residual (and outliers, although it can be finicky sometimes). 
    **input**: _events.root + event#
    **output**: a root file with reconstructed voxels and their information in either "SingleTrack" or "Other" folder
3. MuonDecayEvents
    **purpose**: The reconstruction function specialized for Muon Decay detction. Reconstruct decay events by matching the time gap and hysical overlap of normal two events. 
    **input**: _events.root
    **output**: a "_MuonDecay" folder inside of which are reconstructed muon event root files (empty if no detected muon event) + a _MuonDecay_summary.root file that summarize all events analyzed in the operation. 
4. 3D_EventGamma
    **purpose**: The reconstruction function with added spacial veto logic, useful for gamma/beta source analysis.
    **input**: _events.root
    **output**: two folders "_OnlyX23_ZLeft" and "_OnlyX23_ZRight", inside of which are reconstructed event roots that survived the veto. 
5. 3D_AllEventsSide.sh
    **purpose**: a bash script that will automatically perform ./3D_EventTrajectorySide on input _event.root for each .pngs in input folder. 
    **input**: _event.root + EventDiaplay_PNG_Filtered"
    **output**: two folders ""SingleTrack" and "Other" containning reconstructed root files for events that have either passed or not passed PCA thresholds. 
6. All_MuonDecay.sh
    **purpose**: a bash script that automatically perform ./MuonDecayEvents on all _event.root files in a folder. 
    **input**: /path/to/folder/containng/subrun/subfolders/
    **output**: a user specified destination folder for all reconstructed muon decay events. 
7. All_GammaBeta.sh
    **purpose**: a bash script that automatically perform ./3D_EventGamma on all _event.root files in a folder. 
    **input**: /path/to/folder/containng/subrun/subfolders/
    **output**: user specified destination folders for all reconstructed gamma and beta events respectively. 
8. All_process.sh
    **purpose**: a bash script that automatically perform ./EventDisplayLite and ./3D_AllEventsSide.sh on all _event.root in a folder. 
    **input**: /path/to/folder/1-3levels/above/_events/root/files/
    **output**: a user specified destination folder for all reconstructed root files for events that have passed PCA thresholds (So SingleTrack events; The difference between 3D_AllEventsSide.sh and this is that this does not need an input PNG folder and it can take in folders directly instead of _event.roots)
9. VoxelTrackLength
    **purpose**: Calculates the fitted tracklength of the trajectory in each of the voxels in the event.
    **input**: /path/to/folder/containing/3D/roots/ or a _3DEvent.root
    **output**: a folder containing VoxelPathLengths.root for each 3D root file processed + a summary root file, "A_Summary_AllFiles.root", with all voxels, all trajectories, and all track lengths of the 3D root files in the input folder. 
10. GainCalculations
    **purpose**: Calculates channel gain of each channel with channel charges and track lengths of each voxel. 
    **input**: A_Summary_AllFiles.root
    **output**: a folder containing histograms of each channel's gain response for each channel processed + a summary root file, "ChannelGainAnalysis.root", containing gain information for all channels analyzed. 
11. CalibPlots
    **purpose**: Visualized and compiles channel gain information
    **input**: ChannelGainAnalysis.root
    **output**: gain_calibration.txt, GainCalibrationMaps.png, and GainCalibrationMaps.root
12. SingleTrackFilter
    **purpose** A single track filter that operates on 2D histograms from EventDisplayLite by fiting a line through the histograms and checking residual levels. 
    **input**: _events_display.root
    **output**: two folders, "Accepted_events" and "Rejected_events", containing passed and filtered pngs respectively. (Pretty much redundant now with PCA in 3D reconstruction function)

**TO RUN**

1. **EventDisplayLite**
```
./EventDisplayLite -i /path/to/_events.root
```
This function will perform additional filtering on TimeGroupedEvents by setting minimum number of hits over a threshold and then plot out 3-view histograms for each passed event. Output will be one _display.root and a EventDisplay_PNG_Filtered folder with .png for each passed event. 


2. **3D_EventTrajectory** (FOR FLAT DETECTOR CONFIGURATION, CHECK WITH DETECTOR OPERATOR IF DETECTOR IS IN SUCH ORIENTATION. THIS FUNCTION DOES NOT HAVE ADVANCED PCA FITTING <8/4/2025>.)
For a single event of interest use the following
```
./3D_EventTrajectory /path/to/_event.root/ event#
```
    for example 
    if I wish to operate on event#4589, use
    ```
    ./3D_EventTrajectory /path/to/_event.root/ 4589
    ```
This will reconstruct 3D voxels from _event.root for the specified event and output one _GroupRes.root file in a 3D_Dislay_root folder. 

**3D_AllEvents.sh**
Alternatively, to perform reconstruction on all events filtered from EventDisplayLite use the following
```
./3D_AllEvents.sh /path/to/_event.root/ /path/to/EventDisplay_PNG_Filtered/ 
```
This will perform ./3D_EventTrajectory on all events listed in the EventDisplay_PNG_Filtered folder and output one _GroupRes.root file in the 3D_Display_root folder for each event. 

3. **3D_EventTrajectorySide** (FOR SIDE DETECTOR CONFIGURATION, CHECK WITH DETECTOR OPERATOR IF DETECTOR IS IN SUCH ORIENTATION. THIS FUNCTION INCLUDES PCA FITTING.)
For a single event of interest use the following
```
./3D_EventTrajectorySide /path/to/_event.root/ event#
```
    for example 
    if I wish to operate on event#4589, use
    ```
    ./3D_EventTrajectorySide /path/to/_event.root/ 4589
    ```
This will reconstruct 3D voxels from _event.root for the specified event and output one _GroupRes.root file in a 3D_Dislay_root_SingleTrack or _Other folder. 

**3D_AllEventsSide.sh**
Alternatively, to perform reconstruction on all events filtered from EventDisplayLite use the following
```
./3D_AllEventsSide.sh /path/to/_event.root/ /path/to/EventDisplay_PNG_Filtered/ 
```
This will perform ./3D_EventTrajectorySide on all events listed in the EventDisplay_PNG_Filtered folder and output one _GroupRes.root file in the 3D_Display_root folder for each event. 

4. **MuonDecayEvents**
```
./MuonDecayEvents /path/to/_events.root
```
This function will detect muon decay events 



For calculating the track length in each voxel for a given event, use the following: 
```
(for single 3D event)
./VoxelTrackLength /path/to/_3D_Display.root/

(for a folder of 3D events)
./VoxelTrackLength /path/to/_3D_Display_root/
```
This will perform the tracklength calculation for one or many 3D events. When more than one 3D event is analyzed, in addition to a voxelTrackLength.root file for each 3D event processed, an A_Summary_AllFiles.root will also be created, that contains all essential information for each hits recorded in the 3D events processed, all essential track information for each event processed, and voxelTrackLength information for each voxel analyzed. 



For calculating ADC-MeV calibration for each channel, use the following:
```
./GainCalculations /path/to/A_Summary_AllFiles.root/
```