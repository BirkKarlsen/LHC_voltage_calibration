# LHC voltage calibration
Analysis of the beam based voltage calibration performed on the LHC RF system cavity-by-cavity on June 25th 2022. The ```analyse_profiles```-folder contains scripts to perform analysis on beam profiles to find the synchrotron frequencies and the corresponding voltages. The ```simulations```-folder contains scripts to run simulations recreating the conditions in the LHC during the MD, generate beams in the SPS and launch parameter scans in the lxplus. Lastly, ```impedance_models``` contains impedance models of LHC flat-bottom. 

## What you need:
  * Directory called ```data_files/generated_beams/``` and ```data_files/voltage_calibration_sorted/```
  * Python libraries: ```numpy```, ```matplotlib```, ```scipy``` and ```h5py```
