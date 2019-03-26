%% combine efe's recording
% combine the left and right recordings from two computers 
% does not work under linux because the separator in directory is
% different.

% next to do
% SUA & MUA channel number, line 150 in lr_structure data

%% clear all
clear all; close all; clc;
left_right_kws = {'_MR_', '_ML_', '_SR_', '_SL_'};
left_kw = 'L';
right_kw = 'R';
keyword = '.*mat';  % '.' any character; '*' any number of character

% param for structure data
% param.data_keyword = '*_nw_*';
param.mua_keyword = '*_nw_*';
param.sua_keyword = '*_su_*';
param.sua_num_c = 100.0;
param.right_adjustment = 16;
param.chan2ele = 1;  % channel to electrode adjustment
param.float2int = 100;  % mutiple electrode number by a constant

param.trig_keyword = '_trig';
param.ID_keyword = '_IDstim';
param.sound_keyword = '_sound';
param.bin_size = 0.0005; % bin size in seconds; 
param.master_kw = 'M';
param.slave_kw = 'S';
param.equal_threshold = 5e-4; % if the difference between two values are less than this, they are considered as equal.


%% select the files to be processed
% old files probably cannot be batch processed due to lack of embedded
% duration files
% use File_Detection to get list of files to be processed

starting_path = 'D:\Google_Drive\Lab_Data\';
original_data_folder = uigetdir(starting_path, 'Select the Folder where the data is stored'); 
original_data_folder = [original_data_folder, '\']; % add a "\" for completeness
new_folder = 'matrix_data\';
new_path = [original_data_folder, new_folder];
mkdir(new_path); % create a new folder to store the data

[filenames] = File_Detection(keyword, original_data_folder);

% group the files according to their file name: 
[fn_pairs, unpaired_files] = find_data_pairs(filenames);


%% processed unpaired files
disp('Process unpaired files.......\n')
if ~isempty(unpaired_files{1, 1})
    iii = 1;
    while iii < (length(unpaired_files)+1)
        % go through all unpaired files
        tmp_filename = unpaired_files{iii};
        fn_splited = strsplit(tmp_filename, '_');
        tmp_hem = fn_splited{2};

        % only process files from the master computer
        if tmp_hem(1) == param.master_kw
            % 2nd part: if from the same recording, structure and combine them
            tmp_result = lr_structure_data(original_data_folder, tmp_filename, param); % convert the data into desired format
            
            header = tmp_result.header;
            spiketrains = tmp_result.YYY; 
            TCSE = tmp_result.TCSE;
            stim_codes = tmp_result.stim_codes;
            stim_waveforms = tmp_result.stim_waveforms;
            stim_sr = tmp_result.stim_sr;
            spike_waveforms = tmp_result.spike_waveforms;
            birdids = tmp_result.birdid;
            units = tmp_result.unit;
            
            num_dpts = size(spiketrains);
            num_dpts = num_dpts(2);

            if strcmp(tmp_hem(2), right_kw)  % if from right hemisphere
                header(:, 2) = header(:, 2) + param.right_adjustment; % electrodes in the right hemisphere are represented by bigger numbers
                units = units + param.right_adjustment;
            end

            header(:, 2) = (header(:, 2) - param.chan2ele) * param.float2int; % multiple it by 100 to convert to int
            units = (units - param.chan2ele) * param.float2int;
            % does not remove hemisphere information if unpaired
            % fn_splited{2} = '';
            tmp_combined_fn = strjoin(fn_splited, '_');

            tmp_name = [new_path, 'matrix_', tmp_combined_fn]; % use matrix to indicate matrix representation of the data
            eval(['save ', tmp_name, ' header spiketrains stim_codes stim_waveforms stim_sr spike_waveforms birdids units -v7']); % save the files in the new folder using eval
            disp([num2str(iii), 'th file processed.']);
        else
            disp('Unpaired files from slave PC will not be processed!');
        end

        iii = iii + 1;
    end
end


%% processed paired files 
% two files with the same name except left and right belong to the same group. 
disp('Process paired files.......\n')
iii = 1;
if ~isempty(fn_pairs{1})  % only process if paired files were found
    while iii < (length(fn_pairs)+1)
        tmp_pair = fn_pairs(iii, :);
        tmp_filename_1 = tmp_pair{1};
        fn_splited = strsplit(tmp_filename_1, '_');
        tmp_hem_1 = fn_splited{2};
        tmp_filename_2 = tmp_pair{2};
        fn2_splited = strsplit(tmp_filename_2, '_');
        tmp_hem_2 = fn2_splited{2};

        % 2nd part: if from the same recording, structure and combine them
        result_1 = lr_structure_data(original_data_folder, tmp_filename_1, param); % convert the data into desired format
        header_1 = result_1.header;
        spiketrains_1 = result_1.YYY;
        TCSE_1 = result_1.TCSE;
        
        result_2 = lr_structure_data(original_data_folder, tmp_filename_2, param); % convert the data into desired format
        header_2 = result_2.header;
        spiketrains_2 = result_2.YYY;
        TCSE_2 = result_2.TCSE;
        
        if result_1.master == 1  % if result_1 is from master computer
            % copy the 4th column of header_1 into header_2
            [stim_dict, dur_dict] = header2dict(header_1);
            header_2 = correct_header(header_2, stim_dict, dur_dict);
            stim_codes = result_1.stim_codes;
            stim_waveforms = result_1.stim_waveforms;
            stim_sr = result_1.stim_sr;
        else
            [stim_dict, dur_dict] = header2dict(header_1);
            header_1 = correct_header(header_1, stim_dict, dur_dict);
            stim_codes = result_2.stim_codes;
            stim_waveforms = result_2.stim_waveforms;
            stim_sr = result_2.stim_sr;
        end

        num_dpts = size(spiketrains_1);
        num_dpts = num_dpts(2);
        num_dpts_2 = size(spiketrains_2);
        num_dpts_2 = num_dpts_2(2);
        num_dpts = min(num_dpts, num_dpts_2);

        if strcmp(tmp_hem_1(2), left_kw)
            header_2(:, 2) = header_2(:, 2) + param.right_adjustment; % electrodes in the right hemisphere are represented by bigger numbers
            result_2.unit = result_2.unit + param.right_adjustment; % change unit name associated with spikewaveform
            header = cat(1, header_1, header_2);
            spiketrains = cat(1, spiketrains_1(:, 1:num_dpts), spiketrains_2(:, 1:num_dpts));
        else
            header_1(:, 2) = header_1(:, 2) + param.right_adjustment;
            result_1.unit = result_1.unit + param.right_adjustment;
            header = cat(1, header_2, header_1);
            spiketrains = cat(1, spiketrains_2(:, 1:num_dpts), spiketrains_1(:, 1:num_dpts));
        end
        
        % process spike waveforms
        spike_waveforms = [result_1.spike_waveforms; result_2.spike_waveforms];
        birdids = [result_1.birdid; result_2.birdid];
        units = [result_1.unit; result_2.unit];
        units = (units - param.chan2ele) * param.float2int;

        clearvars header_1 header_2 spiketrains_1 spiketrains_2
        % add sanity check
        % 1. stimulus ending time should be larger than stimulus starting time
        if any(header(:, 6) - header(:, 5) < param.equal_threshold)
            disp('!!!Stimulus starting or ending time is wrong!!!\n')
        end
        % 2. pre_stim + after_stim == isi * 1.25
        % 3. for every electrode, the trial2stim dictionary should be the same
        for zzz = 1 : length(header)
            if stim_dict(header(zzz, 3)) ~= header(zzz, 4)
                disp('!!!Stimulus trial & code does not match!!!\n')
                break;
            end
            if abs(dur_dict(header(zzz, 3)) - (header(zzz, 6) - header(zzz, 5))) > param.equal_threshold
                disp('!!!Stimulus duration & trial does not match!!!\n')
                break;
            end
        end
        % 4. resolution * spikes.shape[1] == pre_stim + after_stim
        st_size = size(spiketrains);
        if abs(st_size(1, 2) * header(1, 9) - (header(1, 7) + header(1, 8))) > param.equal_threshold
            disp([num2str(iii), 'th file, ', tmp_filename_1,', ', tmp_filename_2]);
            disp('!!!Spike train dimension does not match with header!!!\n')
        end

        header(:, 2) = (header(:, 2) - param.chan2ele) * param.float2int; % multiple it by 100 to convert to int
        fn_splited{2} = '';
        tmp_combined_fn = strjoin(fn_splited, '_');

        tmp_name = [new_path, 'matrix_', tmp_combined_fn]; % use matrix to indicate matrix representation of the data
        eval(['save ', tmp_name, ' header spiketrains stim_codes stim_waveforms stim_sr spike_waveforms birdids units -v7']); % save the files in the new folder using eval
        disp([num2str(iii), 'th file processed.']);

        iii = iii + 1;
    end
end
%% clear variables
% clear all;

