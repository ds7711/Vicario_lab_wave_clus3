# Vicario Lab Spike Sorting & Data Pre-Processing (Efe's setup)
 
 Multi-unit, single-unit detection and conversion to `.npz` data for flexible and efficient analysis in Python.
 
 ## Outline

1. Single-unit detection: spike sorting. Use matlab script `vicario_lab_spikesorting.m`, which calls [wave_clus3](https://github.com/csn-le/wave_clus) for fully unsupervised spike sorting.
2. Multi-unit detection: thresholded MUA. Use spike2 script `STD_Threshold2Matlab.s2s`.
3. Export MUA & SUA data stored .smr files into .mat format. Use spike2 script `Spike2Matlab_AllChannel_Batch_v0.1.s2s` or use the export function in `STD_Threshold2Matlab.s2s`.
4. Combine data from master & slave computer based on filenames. Use matlab script `combine_efe_left_right.m`.
5. For each experiment condition (recording session, corresponds to one pair of .smr files from master and slave computer), combine data from different birds and store the data in `.npz` format. Use functions in `mdlab.py`.
6. Analyze data in Python based on `mdlab.py`.

## Procedure

### I. Raw data organization

1. Copy all `.smr` files from the Master and Slave computer into one folder.
	- This folder should only contains `.smr` files and no other files/subfolders.
2. Name recording files according to convention `[BirdID]_[M/S][L/R]_[Experiment]_[extra].smr`. 
	- Square brackets are used for notation purprose only. Do NOT use it in your filename. Underscore `_` is used to separate different parts of the filename.
	- `[BirdID]` must consist of EXACTLY 5 characters, with 2 capital letters follwed by 3 arabic numerals. E.g., `RD001`.
	- `[M/S][L/R]`: `M/S` represents Master and Slave computer, respectively. `L/R` represents Left and Right hemisphere. E.g., `ML` represents recording from the Master computer, which is connected to the electrodes in the left hemisphere. `SR` represents recording from the Slave computer, which is connected to the electrodes in the right hemisphere.
	- `[experiment]_[extra]` is optional and can be used to contain extra information about the recording. E.g., experiment name/condition/date.
	- A pair of recording from the Master and Slave computer may look like: `RD001_MR_variant.smr` & `RD001_SL_variant.smr`.
	- The filenames of the data from the Master & Slave computer in the same recording session/condition should be identical except the `[M/S][L/R]` field. This is essential for combining data from Master and Slave computer.
3. Increase the number of "free" channels that each .smr file can store to 300.
	- Both Multi-units & single-units are stored as wavemark channels in `.smr` files. Number of "free" channels in a `.smr` file determines how much wavemark channels can be writtent. Therefore, it is important to ensure each `.smr` file has enough free channels to hold the detected multi-unit and single-unit produced in the following spike-sorting thresholding procedure.
	- Run `increase_free_channels.s2s` in Spike2 to increase the number of free channels to 300 that the `.smr` files in the folder can hold. The script is in the Spike2Matlab folder.
	- Note, if the latest configuration files were used during data acquisition, the number of free channels is around 300 and there is no need to do this step.

### II. Spike-sorting

0. See below for important dependencies for the scripts in this step to work.
1. Run `vicario_lab_spikesorting.m` for conduct fully automated spike-sorting.
2. If the script is used for the 1st time on a new computer, one first has to change the script and matlab path.
	- Add Spike_Sorting folder and _its subfolders_ to Matlab path.
	- Change `param.CED64_path` has to be the directory where CEDS64ML is located.
	- Criteria forvalid SUA can also be changed in param struct. E.g., number of minimum spikes required, threshold for contamination rate (percentage of inter-spike intervals less than 2ms).
3. Detected single-units are stored as wavemark channels in the original `.smr` files.
	- If one runs spike-sorting several times on the same `.smr` files, multiple copies of SUA will be written into the `.smr` files as wavemark channels.
4. In the `.smr` file, detected SUA wavemark channel has prefix `su`. e.g., su_3_1 means the 1st classes sorted from channel 3. 
	- Channel number = elecctrode number + 1.
	- (optional) In the `.smr` file, MUA detected by [wave_clus3](https://github.com/csn-le/wave_clus) has prefix `mu`. `mu_3` means MUA from channel 3. channel 3 corresponds to electrode 2. This is used to distinguish the MUA detected from those in Spike2 (next step). 
5. Note:
	- Spikesorting often takes a long time and may need to run over night or several days, depending on the amount of data to be processed. It is recommended to first test a few files before running on all files. 
	- For long recordings, spike-sorting uses a lot of RAM and insufficient ram will make the whole process significantly slower. A minimum of 16GB RAW is recommended. 32GB will be enough for most cases.

### III. Thresholded Multi-unit detection

1. Run `STD_Threshold2Matlab.s2s` in Spike2 to obtain thresholded multi-units.
	- Change `X` at line `var NumStd := X;` to what threshold is used for MUA detection. The default is 2.
	- The advantage of using `STD_Threshold2Matlab.s2s` over [wave_clus3](https://github.com/csn-le/wave_clus) is that we have explicit control of the threshold used for action potential detection.
2. Detected multi-unit action potentials are also stored as wavemark channels in the original `.smr` files.
	- If thresholding scripts were run several times on the same `.smr` files, multiple copies of MUA will be written into the `.smr` files as wavemark channels.
3. In the `.smr` file, detected MUA wavemark channel has prefix `nw`. `nw-4` means MUA from channel 4, detected by thresholding.
4. 

### IV. Export MUA & SUA (wavemark channels) to `.mat` files.

1. Run `Spike2Matlab_AllChannel_Batch_v0.1.s2s` in Spike2 to export desired channels to `.mat` files.
	- One can specify what channels will be exported by changing parameters in the script. The default is all wavemark channels, trig, IDstim, & sound channels.
	- Because both MUA & SUA are stored as wavemark channels, one must first do both spike-sorting SUA and thresholded MUA for this script to export both to `.mat` files.
2. `STD_Threshold2Matlab.s2s` also has a `var exportflag% := X;` that can be used to export data to `.mat files`. 
	- The default is `X=0`, which does not export.
	- If changed to `X=1`, the script will detect thresholded MUA, write them into wavemark channels, and then export all wavemark channels to `.mat` files. Note, at this step, if spike-sorting has not been completed, the SUA wavemark channels do not exist yet and the exported `.mat` files will only contain MUA data.
3. After this step, all `.mat` files will be in a folder.
4. Run spikesorting first and then thresholding can make the analysis easier.
	- If one wants to threshold and spikesorting second, do not export during thresholding. Instead, use `Spike2Matlab_AllChannel_Batch_v0.1.s2s` to export both MUA and SUA data after both procedure is done.
	- One could also use `Spike2Matlab_AllChannel_Batch_v0.1.s2s` to selectively export only SUA or MUA by skipping spike-sorting or thresholding.

### V. Combine data from Master and Slave and simplify the data.

1. Run `combine_efe_left_right.m` in Matlab to combine the data from "master" & "slave" computer.
	- When prompted, select the folder where exported `.mat` files are stored from the previous step.
	- Note, the filename from the master and slave computer must follow the naming convention specified above (part 1).
	- After combination, a new folder contains keyword `matrix` will be created within the original `.mat` folder.
	- Combined data will be stored in this new `matrix` folder. The filenames of the `.mat` file should not contain `[M/S][L/R]` if data from Master and Slave computer both exist and form a "pair".
2. In case filenames from the Master & Slave computer do not match OR data from Slave computer are missing, the script will only process data from the Master computer and skip those from Slave computer.
	- If only data from the Master computer is processed, the output `.mat` files in the `matrix` folder will contain `M[L/R]`.
3. The newly created `.mat` files inside the `matrix` folder will be used for further analysis.

### VI. Combine data from different bird in the same condition.

1. Organize the `.mat` files from previous step in different folders based on the experiment condition. E.g., each condition should have a folder contains `.mat` data from different birds.
2. Run `mat2npz.py` in Python to combine data from different birds and save data in `.npz` format.
	- After this step, for each experiment condition, you should only have `.npz` file that contains data from all birds.
	- `.npz` file contains the spike trains from each trial, what the stimulus code is, how the stimulus waveform looks like (extracted from the sound channel in the Master computer). For single-unit, it also contains its average spike waveform, which can be used to classify a unit into wide/narrow type.

## Dependencies
### Spike-sorting in Matlab.

1. Developed with Matlab 2016b (may or may not work on other versions).
2. Spike2 MATLAB SON Interface (included) is required to read and modify .smr files.
3. Wave_clus 3 (as 2018-11-11)  [https://github.com/csn-le/wave_clus](https://github.com/csn-le/wave_clus)
4. Some functions in wave_clus.m has been modified. Including but not limited to:
	- wave_clus_OpeningFcn()
	- load_data_button_Callback()
5. Signal processing toolbox of Matlab may be required.

### Data analysis in Python.

1. Developed with Python 2.7 in [Anaconda](https://www.anaconda.com/distribution/).
2. Whenever `mdlab.py` is needed, it has to be included in the same folder as the script that is being executed. e.g., see `mat2npz.py`.

## Notes

2. In helper_spike2_scripts, scripts can delete old wavemark channels, add free channels, threshold to detect spikes, and export wavemark channels to matlab format.
3. Alternatively, one could use neo package in Python to directly read .smr files.

## Things to do

1. Add a tutorial about how to use `mdlab.py` to analyze spike train data.
2. Support recording with gaps. 

## multiple trials (experimental, do not change!)

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
