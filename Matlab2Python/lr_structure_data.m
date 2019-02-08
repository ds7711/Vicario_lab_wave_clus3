%% structure_data:
%       takes the exported .mat file from spike2 (thresholded mua) or spike-sorting
%   algorithm (waveclus) and organized the data into a matrix.
%   
%  structure of the matrix:
%       1st column: numerical birdID (color converted to number, see table below for conversion)
%       2nd column: electrode ID (channel number - 1)
%       3rd column: trial number (1 to number of trials)
%       4th column: stimulus code
%       5th column: trial starting time
%       6th column: trial stimulus ending time
%       7th column: pre_stim duration
%       8th column: recording duration starting from stimulus onset
%       9th column: resolution (bin_size); columns to its right are spike data
%       10th to end: 0 indicates no spike and 1 indicates spike
%                   pre-stimulus: 1 / 4 of ISI, maximum 1 seconds;
%                   stimulus: 3/4 of ISI, max stimulus duration + 0.5 sec
%           
%  customized function requirements:         
%       File_Detection: get filenames in the target folder
%       tcse: calculate the Trial, code, starting, ending time
%       Spikes_Between: get the number of spikes between specified range
%

function [result] = lr_structure_data(path, filename, param)
%% load data and define keywords
load([path, filename]);
% find variables that storing the data 
mua_keyword = param.mua_keyword;
sua_keyword = param.sua_keyword;
mua_name_list = who(mua_keyword);
sua_name_list = who(sua_keyword);
data_name_list = [mua_name_list; sua_name_list];

mua_splitter = mua_keyword(2:end-1);
sua_splitter = sua_keyword(2:end-1);

all_fns = who();  % list all filenames & find trig & IDstim
for iii = 1 : length(all_fns)
    str = all_fns{iii, 1};
    if contains(str, param.trig_keyword, 'IgnoreCase', true)
        trig_name = {str};
    end
    if contains(str, param.ID_keyword, 'IgnoreCase', true)
        stimID = {str};
    end
    if contains(str, param.sound_keyword, 'IgnoreCase', true)
        sound = {str};
    end
end
Trig = eval(trig_name{1,1});
stimID = eval(stimID{1, 1});
sound = eval(sound{1, 1});

%% Extract or import stimulus duration matrix 
% [stimdur, stim_name_list] = dur_extract(stimID); 

% test whether there is zero data in the Trig variable, unknown strange
% error from converted matlab files.
if any(Trig.times(3:length(Trig.times)) == 0)
    tmp_logidx = Trig.times(3:length(Trig.times)) == 0;
    tmp_idx = find(tmp_logidx, 1);
    tmp_idx = tmp_idx + 1;
    Trig.times = Trig.times(1 : tmp_idx);
    Trig.codes = Trig.codes(1 : tmp_idx, :);
end

SCode_vec = Trig.codes(:, 1);
if all(SCode_vec==255)
    stimdur = [1, 0];
else
    try 
        stimdur = dur_extract(stimID);
        logical_test = isnan(stimdur); % identify whether the stimulus duration is successfully extracted or not (NaN)
        if all(logical_test(:, 2)) || isempty(logical_test)
    %         stimdur = uiimport('C:\Users\md\Dropbox\Data_Analysis\Field_L_2014_Project\FieldL_Project\stimdur_oddball.xlsx');
    %         stimdur = stimdur.data; % import the stimulus duration data interactively 
            stimdur = csvread(strcat(path, 'stimdur.txt'));
        end
    catch 
        stimdur = csvread(strcat(path, 'stimdur.txt'));
    %     stimdur = uiimport('C:\Users\md\Dropbox\Data_Analysis\Field_L_2014_Project\FieldL_Project\stimdur_oddball.xlsx');
    %     stimdur = stimdur.data; % import the stimulus duration data interactively 
    end
end

%% Check whether the stimulus sequence from TextMark channel and Trig channel are the same
% Because stimuli indicated by stimulus code 2 to 10 in oddball experiment
% are actually the same, the stimulus code needs to be first converted
% back. 
% The reason why multiple stimulus code are used for the same stimulus is
% to manipulate the presenting probability of the stimulus. 
% function stimcode_norm takes struct variable trigger and IDstim as input,
% and outputs the updated Trig for later use. 

% Test whether the stimulus sequence from stimID is the same from 
% StimCode_from_Trig = Trig.codes(2 : 3 : end, 1);
% % stimID.codes(:, 1): stim sequence from TextMark channel 
% if ~isequal(StimCode_from_Trig, stimID.codes(:, 1))
%     warndlg('Error in stimulus code: inconsistency between codes from Trig and from stimID');
%     pause();
% end

%% Obtain the Time Range matrix for Response Period & Control Period; 
TCSE = tcse(Trig, stimdur); 
% num_stim = length(unique(TCSE(:, 2))); 
% TCSE: Trial, Code, StartTiming, EndTiming
% TCSE: 
%       1st column: trial number; 
%       2nd column: Stimulus Code; 
%       3rd column: Stimulus onset time;
%       4th column: Stimulus end time.

% extract sound waveform
% number of unique sounds
unique_sounds = unique(TCSE(:, 2));
num_sounds = length(unique_sounds);
stim_codes = NaN(num_sounds, 1);
sound_waveforms = cell(num_sounds, 1);
max_ndpts = -1;
stim_count = 0;
for i = 1 : length(TCSE)
    % get temporary information
    tmp_item = TCSE(i, :);
    tmp_sound_code = tmp_item(2);
    % if stimulus code not in sound_stims
    if ~ ismember(tmp_sound_code, stim_codes)
        stim_count = stim_count + 1;  % index + 1 
        stim_codes(stim_count, 1) = tmp_sound_code;  % store the stimulus code
        tmp_st = tmp_item(3);
        tmp_ed = tmp_item(4);
        tmp_st_idx = int64(tmp_st / sound.interval);
        tmp_ed_idx = int64(tmp_ed / sound.interval);
        % store the stimulus waveform
        sound_waveforms{stim_count, 1} = sound.values(tmp_st_idx : tmp_ed_idx, 1);
        tmp_ndpts = tmp_ed_idx - tmp_st_idx + 1;
        if max_ndpts < tmp_ndpts
            max_ndpts = tmp_ndpts;
        end
    end
    % if all stimuli found, break the loop
    if stim_count >= num_sounds
       break; 
    end
end
% convert sound_waveforms from a cell array to a matrix
stim_waveforms = NaN(num_sounds, max_ndpts);
for i = 1 : num_sounds
    tmp_waveform = sound_waveforms{i, 1};
    stim_waveforms(i, 1:length(tmp_waveform)) = tmp_waveform;
end


%% !!!!!! ------------------------- Important Parameters, Check before use ------------------------------- !!!!
% After_Effect = 0.1;  % The spikes 0.1 sec after stimulus offset are also considered as responses to the stimulus 
isi = mean(diff(TCSE(:, 3))); 

Pre_Stim = min([1, isi / 4.0]) * 2.0;  % include longer baseline
After_Stim = isi - Pre_Stim / 2.0; % adjust for the longer baseline

% !!!!!! ------------------------- Important Parameters, Check before use ------------------------------- !!!!

name_list_innum = zeros(size(data_name_list));
% DATA: 
%   1st column: electrode name
%   2nd column: trial by trial response firing rate
%   3rd column: average response for each stimulus 

num_chan = length(data_name_list); % number of data channels
Response_Spikes = cell(num_chan, 1);
Control_Spikes = cell(num_chan, 1);
chan_num_vec = zeros(num_chan, 1);

for i = 1 : num_chan
    clear all_fns temp_data_name_splited data spiketiming Response_Time_Range Control_Time_Range Response_FR Control_FR average_control;
    % -------------------------extract the channel number from the data_name-------------------------
    str = data_name_list{i, 1};
    name_list_innum(i, 1) = str2double(str(isstrprop(str, 'digit'))); % data name in numeric form
    data = eval(str);
    
    [temp_data_name_splited, ~] = strsplit(str, {mua_splitter, sua_splitter}, 'CollapseDelimiters', true);

    chanstr = temp_data_name_splited{1, end}; 
    chan_num = chanstr2double(chanstr, param);
%     disp(str);
%     disp(chanstr);
%     disp(chan_num);

    chan_num_vec(i, 1) = chan_num;
    % Calculate the Firing Rate to stimulus & Spontaneous FR during control Period
    % !!! Normlization for the 1st trial needs to be really careful
    % !!! The control duration for the 1st recording is less than 0.26
    % seconds, change the sampling code in the future it in the future
    spiketiming = data.times; 
%     Response_Time_Range = [TCSE(:,3) TCSE(:,4) + repmat(After_Effect, size(TCSE(:, 4)))];
    Response_Time_Range = [TCSE(:,3) - Pre_Stim, TCSE(:,3) + After_Stim];
    [~, Response_Spikes{i, 1}] = Spikes_Between(spiketiming, Response_Time_Range);  % Absolute Spike Timing within the specified response/stimulus window

    Control_Time_Range = [TCSE(:,3) - repmat(Pre_Stim, size(TCSE(:, 3))) TCSE(:,3)];
    [~, Control_Spikes{i, 1}] = Spikes_Between(spiketiming, Control_Time_Range);    % Absolute Spike Timing within the specified Control window
end

%% Transform the absolute spike_timing into relative timing (alligned to the stimulus onset) 
num_trial = size(TCSE, 1); % number of trials 
start_timing_col = 3; % The 3rd column of TCSE stores the absolute onset time of the stimulus
Relative_Response_Spikes = cell(size(Response_Spikes)); % Store Relative Spiketiming cells 
for i = 1 : num_chan
    abs_spike_timing = Response_Spikes{i, 1}; 
    relativ_spike_timing = cell(size(abs_spike_timing));
    for j = 1 : num_trial 
        relativ_spike_timing{j, 1} = abs_spike_timing{j, 1} - TCSE(j, start_timing_col); % relative spike timing = absolute spiketiming - starting timing of corresponding trials
    end
    Relative_Response_Spikes{i, 1} = relativ_spike_timing; 
end

%% Count spikes within defined time bin 
bin_size = param.bin_size; % bin size in seconds; 

edges = -1 * Pre_Stim : bin_size : After_Stim; 

% x = zeros(num_chan * num_trial, 5);
xxx = [];
YYY = false(num_chan * num_trial, length(edges));

k = 1;
for i = 1 : num_chan
    temp_chan_data = Relative_Response_Spikes{i, 1};
    for j = 1 : num_trial
        temp_chan_trial_data = temp_chan_data{j, 1};
        if length(temp_chan_trial_data) < 3
            temp_chan_trial_data = [temp_chan_trial_data; NaN; NaN]; %#ok<AGROW>
        end
%         temp_Y = histc(temp_chan_trial_data, edges) ./ bin_size;
        temp_Y = histc(temp_chan_trial_data, edges);
%         disp([num2str(i), ', ', num2str(j), ', ', num2str(size(temp_Y))])
        YYY(k, :) = logical(temp_Y');
        k = k + 1;
    end
    xxx = [xxx; [repmat(chan_num_vec(i, 1), num_trial, 1), TCSE]]; %#ok<AGROW>
end
YYY = YYY(:, 1 : end - 1); % YYY is spike trains
%clearvars -except xxx YYY edges stim_dur chan_num_vec data_name_list TCSE bin_size filename Pre_Stim After_Stim; 

%% extract numerical birdID from filename
expression = '\w+\d+';
birdid = regexp(filename, expression);
birdid = filename(birdid:end);
if or(birdid(1:2) == 'MD', birdid(1:2) == 'md')
    birdid = birdid(3 : 7);
end
tmp1 = num2str(double(birdid(1))); % letter is converted to number; e.g., A --> 65
tmp2 = num2str(double(birdid(2)));
num_birdid = str2num([tmp1, tmp2, birdid(3:5)]); %#ok<ST2NM>

rep_ones = ones(size(xxx, 1), 1);
header = [rep_ones * num_birdid,... 
                xxx, Pre_Stim * rep_ones, After_Stim * rep_ones,...
                rep_ones * bin_size];
result.header = header;
result.YYY = YYY;
result.TCSE = TCSE;
result.stim_codes = stim_codes;
result.stim_waveforms = stim_waveforms;
result.stim_sr = 1.0 / sound.interval;
if all(Trig.codes(:, 1)==255)
    result.master = 0;  % whether it is master computer or slave computer
else
    result.master = 1;
end
% clearvars -except header YYY filename path;
end