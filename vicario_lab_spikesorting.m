
%% implementation details.
% clear everything
close all; clc; clear; 

% CEDLIB path (% change to where CEDS64ML is)
param.CED64_path = 'C:\Users\Mingwen\Dropbox\MD_scripts\Spike_Sorting\CEDMATLAB\CEDS64ML';
wave_clus_param = set_parameters();  % get parameters from wave_clus 

% parameters
param.adc_code = 1;  % code label for the stimulus or recording channel
param.record_1st_idx = 2;  % the first channel that is the recording from neurons, 1st adc channel is stimulus by default
param.int16_factor = 32767;  % change wavefrom from real number to integers
param.sr = 25000;  % default sampling rate assumed by wav clus
param.LEAST_NUM_SPIKE = 100; % the number of spikes in cluster is less than this number, it's not saved
param.min_chans_allowed = 100;  % minimum number of channels allowed in the .smr file
param.cr_refractory = 0.002;  % refractory period used to calculate contamination rate
param.w_pre = wave_clus_param.w_pre;
param.w_post = wave_clus_param.w_post;
param.shift = -wave_clus_param.w_pre;  % spikes need to be aligned by their peak --> determined by number of pre-event points in wav_clus
param.ndpts = wave_clus_param.w_pre + wave_clus_param.w_post;  % number of data points for each spike waveform
param.dRate_wavemark = 1000;  % not fully understood
param.separator = '_';  % separator used to separate channel number from unit number
param.min_fr = 0.1;  % minimum number of spikes per second
param.contamination_rate = 0.02; % for debuging, set to 0.1
param.normalize_wav_mag = 0;  % whether normalize the spikes sorted to avoid clipping
param.max_classes = 7;  % maximum number of classes to check
param.mua_label = 'mu';  % prefix for MUA
param.sua_label = 'su';  % prefix for SUA
param.tmp_fn = '___tmp_ss.mat';  % temporary filename used to interact with wave clus
param.wav_data_keyword = '*.smr';  % '.' any character; '*' any number of character
param.min_trial = 3;  % number of times to try before giving up on a channel

% data path (optional, where data are stored, not necessary for using the scripts)
param.data_directory = 'C:\Users\Mingwen\Desktop\Wave_Clus_ValidationProject';

wav_lists = 1 : param.min_chans_allowed;  % test which channel store waveforms within this range.


%% 2. detect all .smr files in the folder and store them as a cell.
Original_Data_Folder = uigetdir(param.data_directory, 'Select the Folder where the original data is stored');  
FN_struct = dir(strcat(Original_Data_Folder, '\', param.wav_data_keyword)); 
num_files = length(FN_struct);
FileName = cell(num_files, 1);
for i = 1 : num_files
    FileName{i,1}=FN_struct(i,1).name;
end
Data_Store_Folder = Original_Data_Folder;
% Data_Store_Folder = uigetdir(starting_path, 'Select the Folder where you want to store the sorted data');  
cd(Data_Store_Folder); % change the working directory to where the data would be stored for convenience
clearvars wav_data_keyword i FN_struct ans wav_data_keyword % clear unused variables
% create a log file
fid = fopen(fullfile(Data_Store_Folder, 'SpikeSortingLog.txt'), 'wt');


%% 3. Loop through each .smr file.

% load CED64 library
cedpath = param.CED64_path;
if isempty(getenv('CEDS64ML'))
    setenv('CEDS64ML', cedpath);
end
cedpath = getenv('CEDS64ML');
addpath(cedpath);
% load ceds64int.dll
CEDS64LoadLib( cedpath );

% cedpath = param.CED64_path;  % add path to CED code
% addpath(cedpath);
% % load ceds64int.dll
% CEDS64LoadLib( cedpath );


% loop through .smr files
% things to add potentially: remove exisiting wavemark channels
for iii = 1 : num_files
    % load ith .smr file
    cd(Data_Store_Folder); % make sure matlab is working the right folder
    temp_filename = FileName{iii, 1};
    fprintf('%s\n', temp_filename);
    temp_path = ([Original_Data_Folder, '\', temp_filename]); 
    
    % test whether there is a file open, if yes, first close it
    try 
        CEDS64CloseAll();
    catch
        % raise error if previous files are not correct closed
    end
    
    try 
        fhand = CEDS64Open(temp_path, 0);  % write mode

        % if a file fails to open, write it in the log and reload CED64 library
        if (fhand <= 0)  % file cannot be read
            % add the filename to a log file
            fprintf(fid, '\n %s, \t Failed to open, \t %s.\n', temp_filename, datestr(now, 0));
            unloadlibrary ceds64int;

            % load CED64 library
            if isempty(getenv('CEDS64ML'))
                setenv('CEDS64ML', CED64_path);
            end
            cedpath = getenv('CEDS64ML');
            addpath(cedpath);
            % load ceds64int.dll
            CEDS64LoadLib( cedpath ); 
        end
    catch
        fprintf(fid, '%s, \t Failed to open, \t %s.\n', temp_filename, datestr(now, 0));
        continue;
    end
    
    % a) Check if the .smr file could save enough channels.
    max_wav_chans = CEDS64MaxChan(fhand);
    % if channels allowed < 100 --> change the file so that it could
    % contain 100 channels use the spike2 scripts.
    % add later if has time
    if max_wav_chans < param.min_chans_allowed
        fprintf(fid, '%s, \t channels allowed should be increased to 300, \t %s.\n', temp_filename, datestr(now, 0));
        CEDS64Close( fhand );
        continue;
    else  % only do spike sorting if the file has enough free channels
        % test all wavefile channels
        detected_wav_chans = NaN(length(wav_lists), 1);
        kkk = 1;  % index for wavchannels
        for jjj = 1 : length(wav_lists)
            tmp_chan = wav_lists(jjj);
            tmp_type = CEDS64ChanType(fhand, tmp_chan);
            if tmp_type == param.adc_code
                detected_wav_chans(kkk, 1) = tmp_chan;
                kkk = kkk + 1;
            end
        end
        detected_wav_chans = detected_wav_chans(~isnan(detected_wav_chans));
        detected_wav_chans = detected_wav_chans(param.record_1st_idx:end);  % first adc is stimulus channel

        % load the .smr files inside one by one (one or more experiment)
        maxTimeTicks = CEDS64ChanMaxTime( fhand, 1 )+1;
        for mmm = 1 : length(detected_wav_chans)
        % b) List all wave channels. Loop through each wave channel.
            iChan = detected_wav_chans(mmm);
            tic;
            fprintf('\t Channel %d, started at %s\n', iChan, datetime);
            % i. Read the wavechannel.
            try
                [ fRead, data, fTime ] = CEDS64ReadWaveF( fhand, iChan, maxTimeTicks, 0, maxTimeTicks);
                % fVals = readWave(fhand, iChan);
                % fVals = single(fVals);  % store as single to save for memory
            catch
                fprintf(fid, '%s, \t Failed to read channel, \t %s.\n', [temp_filename, '_chan_', num2str(iChan)], datestr(now, 0));
                % CEDS64Close( fhand );
                close all;
                continue;
            end
            
            %% test new wave_clus
            
            % create a structure to pass the filename, fhand, param
            tmp_mat_file = param.tmp_fn;
            tmp_event.param = param;
            tmp_event.fhand = fhand;
            tmp_event.wav_name = temp_filename;
            iChan_str = int2str(iChan);
            tmp_event.iChan = iChan;
            tmp_event.tmp_mat = [Original_Data_Folder, '\', tmp_mat_file];
            
            try 
                %% try clustering, if successful, create wavemark channel
                % ii. Test whether the channel is good.
                    % if a channel is bad, stop current iteration and continue to next wave channel.
                    % otherwise, continue.

                % iii. Call automatic spike-sorting subrouting.
                % iv. Combine all units as the MUA and write into a wavemark channel.
                % v. List all single-units. Loop through each single-unit:
                    % - Check the number of spikes detected.
                    % - Check the inter-spike interval.
                    % - if criterion not met, stop current iteration and continue to next single-unit.
                    % - otherwise, write the single-unit as a wavemark channel.
                    % if multi-units detected but no SUA satisfying the criterion, try the automatic sorting again.

                sr = 1.0 / (CEDS64TimeBase(fhand) * CEDS64ChanDiv( fhand, iChan));
                save(tmp_event.tmp_mat, 'data', 'sr');  % save local file for wave_clus
                
                % set global variable for trying wave_clus multiple times
                global flag;
                flag.MUA = 1;  % assume there is MUA
                flag.SUA = 0;
                flag.trial = 0;
                
                while flag.MUA > 0 && flag.SUA == 0 && flag.trial < param.min_trial
                    % try spike sorting for the first time
                    wave_clus(tmp_event);
                    fprintf('\t Trial %d, *%d* single units were found.\n', flag.trial, flag.SUA);
                    close all;  % close the gui
                end
            
            catch
                fprintf(fid, '%s, \t Failed to cluster or create channel, \t %s.\n', [temp_filename, '_chan_', num2str(iChan)], datestr(now, 0));
                % CEDS64Close( fhand );
                close all;
                fprintf('\t !!! cluster ended unexpectedly...!!! \n');
                
            end
            tmp_ed = toc;
            fprintf('\t Channel %d,   ended at %s\n', iChan, datetime);
            fprintf('\t Elapsed time is %f minutes.\n\n', tmp_ed/60.0);
            clearvars tmp_event;
            continue;

        end
    % c) Once the file is finished, write a log file.
    end
    
    % vi. close the current .smr file and continue to next .smr file.
    CEDS64Close( fhand );
    % delete temporary files
    delete *.s2rx;  % *.dg_01 will be deleted in run_cluster.m
    % warning('off', 'MATLAB:DELETE:FileNotFound');  % suppress this warning
    clearvars data;
    close all;
end

%% 4. Print the message that spike-sorting has been finished. 
CEDS64Close( fhand );
fclose(fid);  % close the log file
% delete temporary files
delete *.dg_01 *.s2rx *.lab;  % *.mat
close all;
CEDS64CloseAll();
