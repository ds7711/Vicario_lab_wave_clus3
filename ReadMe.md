# Vicario Lab Spike Sorting

Use wave_clus3 for fully unsupervised spike sorting. 


## Instructions

1. Add Spike_Sorting folder and _its subfolders_ to Matlab path.
2. Change parameters in vicario_lab_spikesorting.m before on a new computer. 
	- param.CED64_path has to be the directory of CEDS64ML.
	- Criteria for SUA are also defined in param struct. 
3. Put all .smr files in one folder. 
4. Use spike2 scripts to increase the number of supported channels to 300 if original .smr files do not support that. 
	- Create a new folder to store the new .smr files.
5. Run vicario_lab_spikesorting.m and follow the instruction.
6. The detected MUA and SUA will be stored as wavemark channels in original .smr files.
	- MUA has prefix mu. e.g. mu_3 means MUA from channel 3.
	- SUA has prefix su. e.g., sua_3_1 means the 1st classes sorted from channel 3.

## Dependencies

1. Developed with Matlab 2016b (may or may not work on other versions).
2. Spike2 MATLAB SON Interface (included) is required to read and modify .smr files.
3. Wave_clus 3 (as 2018-11-11)  [https://github.com/csn-le/wave_clus](https://github.com/csn-le/wave_clus)
4. Some functions in wave_clus.m has been modified. Including but not limited to:
	- wave_clus_OpeningFcn()
	- load_data_button_Callback()
5. Signal processing toolbox of Matlab may be required.

## Notes

1. One may want to use thresholding to detect spikes before spike sorting.
2. In helper_spike2_scripts, scripts can delete old wavemark channels, add free channels, threshold to detect spikes, and export wavemark channels to matlab format.
3. Alternatively, one could use neo package in Python to directly read .smr files.

## Things to do

1. Support recording with gaps. 


## multiple trials

1. Run wave_clus through vicario_lab_spikesorting. 

    global flag;
    flag.MUA = 0;
    flag.SUA = 0;
    flag.trial = 0;
    
2. For each channel, run wave_clus N times. During each iteration:

    0. Run wave_clus.
        * If no spikes detected, flag.MUA = -1
        * If spikes detected, flga.MUA = 1.
        
        * If SUA classified, flag.SUA = # of SUA
        * Otherwise, flag.SUA = 0
        
        * in all cases, flag.trial += 1
        
    1. If no spikes detected, ()
        * break the loop and continue to next channel.
    2. If spikes detected and SUA detected,
        * break the loop and continue to next channel.
    2. If spikes detected but the number of sorted SUA is 0.
        * loop again.