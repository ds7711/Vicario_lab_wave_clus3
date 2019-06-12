# Vicario Lab Spike Sorting & Data Pre-Processing
 
 Multi-unit, single-unit detection and conversion to `.npz` data for flexible and efficient analysis in Python.
 
 ## Outline

1. Single-unit detection: spike sorting. Use matlab script `vicario_lab_spikesorting.m`, which calls [wave_clus3](https://github.com/csn-le/wave_clus) for fully unsupervised spike sorting.
2. Multi-unit detection: thresholded MUA. Use spike2 script `STD_Threshold2Matlab.s2s`.
3. Export MUA & SUA data stored .smr files into .mat format. Use spike2 script `Spike2Matlab_AllChannel_Batch_v0.1.s2s` or use the export function in `STD_Threshold2Matlab.s2s`.
4. Combine data from master & slave computer based on filenames. Use matlab script 'combine_efe_left_right.m'.
5. For each experiment condition (recording session, corresponds to one pair of .smr files from master and slave computer), combine data from different birds and store the data in .npz format. Use functions in `mdlab.py`.
6. Analyze data in Python based on `mdlab.py`.

## Steps

1. Put all .smr files to be processed in one single folder.
2. Increase the number of "free" channels that each .smr file can store to 300. e.g., use `increase_free_channels.s2s`, the script is in the Spike2Matlab folder.
	- Note: the following steps will not work if the number of "free" channels is not enough.
3. Change parameters & path for `vicario_lab_spikesorting.m` if used on a new computer and run `vicario_lab_spikesorting.m`.
	- The spikesorting often takes a long time and may need to run over night or several days, depending on the amount of data to be processed. It is recommended to first test a few files before running on all files.
	- See "instructions" for details.
4. Run `STD_Threshold2Matlab.s2s` to obtain thresholded multi-unit channels. This script is in the Spike2Matlab folder.
	- Change the `X` at line `var NumStd := X;` to your prefer threshold. The default is 2.
	- Change the `X` at line `var exportflag% := X;`. If `X=1`, .smr files will be exported to .mat files. If `X=0`, the script will only threshold and create multi-unit wavemark channels but NOT export to .mat format.
5. Run spikesorting first and then threshold to make the analysis easier. If one does thresholding and exporting (step 4) first and spikesorting later, single-unit wavemark channels will not be exported to .mat files.
	- If one wants to threshold and spikesorting second, do not export during thresholding. Instead, use `Spike2Matlab_AllChannel_Batch_v0.1.s2s` to export both MUA and SUA data after both procedure is done.
	- One could also use `Spike2Matlab_AllChannel_Batch_v0.1.s2s` to selectively export only SUA or MUA by skipping step 3 or 4.
6. Run `combine_efe_left_right.m` to combine the data from "master" & "slave" computer.
	- Note, the filename from the master and slave computer must follow the following naming convention.
7. Use 

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
	- MUA has prefix mu. e.g. mu_3 means MUA from channel 3. channel 3 corresponds to electrode 2.
	- SUA has prefix su. e.g., sua_3_1 means the 1st classes sorted from channel 3.
	- After data are converted into SpikeData object in Python, representation for MUA can be divided exactly whereas SUA can NOT! It can be used to selectively analyze MUA or SUA data.

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
    flag.MUA = 1;
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
