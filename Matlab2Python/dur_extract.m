% Duration Extract 
% Input: 
%       ID_Stim from exported Spike2 File
% Output: 
%       stim_duration matrix: 
%               1st column: stimulus code; 2nd column: stimulus duration
%       stim_name_list: 
%               The corresponding stimulus name for the stimulus code in
%               stim_duration matrix

function [stim_dur, stim_name_list] = dur_extract(stimID)
stimcode = stimID.codes(:,1);
stimname = stimID.text;

stimcode_set = sort(unique(stimcode));
num_stim = size(stimcode_set, 1);
dur_keyword = 'Duration: ';
dur_keyword_length = length(dur_keyword);
stim_dur = zeros(num_stim, 2);
stim_name_list = cell(num_stim, 1);
for i = 1 : num_stim
    temp_stim = double(stimcode_set(i, 1)); 
    temp_stim = double(temp_stim);  % convert to double
    temp_stim_idx = find(temp_stim == stimcode, 1, 'first');
    stim_name_list{i, 1} = stimname(temp_stim_idx, :);
    temp_stim_name = stim_name_list{i, 1}; 
    temp_start_idx = strfind(temp_stim_name, dur_keyword) + dur_keyword_length;
    temp_end_idx = strfind(temp_stim_name, ';'); % find the index of comma ";"
    temp_end_idx = temp_end_idx(1) - 1; % the first comman denotes the ending of duration
    temp_stim_dur = str2double(temp_stim_name(temp_start_idx : temp_end_idx));
    stim_dur(i, :) = [temp_stim, temp_stim_dur];
end


end