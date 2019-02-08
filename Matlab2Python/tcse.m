% Preprocess and organize the timing, stimcode, & trial information into
% the format of TCSE; 
% TCSE: Trial, Code, StartTiming, EndTiming
% TCSE 1st column: trial number (1st to last trial); 2nd column: Stimulus Code; 3rd column:
% Stimulus onset time; 4th column: Stimulus end time. 
% input: Trig (trig construct of the Matlab file exported from spike2, )
%        stimdur (stimulus duration matrix) 
% output: TCSE

function TCSE = tcse(Trig, stimdur)

% new code
threshold = 0.001; % if the difference between two time are less than 1ms, treat them as equal.

% initialize an empty list to store time records
time_records = NaN(length(Trig.times), 1);
code_records = NaN(length(Trig.codes), 1);

% calculate interval between times
time_diffs = diff(Trig.times);

% put in the first time & code
time_records(1, 1) = Trig.times(1, 1);
code_records(1, 1) = Trig.codes(1, 1);

record_idx = 2;  % index of the valid records in Trig
for i = 1 : length(time_records) - 1
    
    if time_diffs(i, 1) < threshold % if current time = previous, to next loop
        continue;
        
    else  % otherwise, add time & code to the records
        time_records(record_idx, 1) = Trig.times(i+1, 1); 
        code_records(record_idx, 1) = Trig.codes(i, 1);
        record_idx = record_idx + 1;
    end
    
end

STime_vec = time_records(2 : 2 : end, 1);
SCode_vec = code_records(2:2:end, 1); % obtain the stimulus code for each trial
% remove NaN elements
STime_vec = STime_vec(~isnan(STime_vec));
SCode_vec = SCode_vec(~isnan(SCode_vec));
if length(STime_vec) ~= length(SCode_vec)
    fprintf('Potential errors in data analysis!!!!!!!!!!!\n'); 
    num_trials = min(length(STime_vec), length(SCode_vec));
    STime_vec = STime_vec(1:num_trials, 1);
    SCode_vec = SCode_vec(1:num_trials, 1);
end
% shrink to the short ISI

ETime_vec=zeros(size(STime_vec));

% if Trig.times(2) - Trig.times(1) > 0.1  % if the first timestamp is somehow missing, occassionally in the slave computer
%     STime_vec = Trig.times(2:3:end, 1); % obtain the timing of stimulus onset
% else  % of normal, use normal stuff
%     STime_vec = Trig.times(3:3:end, 1); % obtain the timing of stimulus onset
% end
% ETime_vec=zeros(size(STime_vec));
% SCode_vec = Trig.codes(2:3:end, 1); % obtain the stimulus code for each trial



%%%%%%%%%%%%%%% only uses in Efe's data specific %%%%%%%%%%%%%%%%%%%%%%%
% test if the codes are all the same: 255
if all(SCode_vec==255)
    SCode_vec = ones(length(SCode_vec), 1);
end
%%%%%%%%%%%%%%% only uses in Efe's data specific %%%%%%%%%%%%%%%%%%%%%%%
keyset = uint8(stimdur(:, 1));
valueset = stimdur(:, 2);
stimdur_map = containers.Map(keyset, valueset);
Num_Trials = length(STime_vec); % number of trials 

% Convert the stimulus code back 
clear temp;

for i=1:Num_Trials
    % temp=stimdur(SCode_vec(i,1), 2);  % Get the stimcode for trial i and then obtain its duration
    temp = stimdur_map(SCode_vec(i,1));
    ETime_vec(i,1)=STime_vec(i,1)+temp;  % Stimulus End time for each trial
end

clear i temp; 

Trial_vec=1:1:Num_Trials;
Trial_vec=Trial_vec';
TCSE=[Trial_vec double(SCode_vec(1:Num_Trials)) STime_vec ETime_vec];  
% TCSE: Trial, Code, StartTiming, EndTiming
% TCSE 1st column: trial number; 2nd column: Stimulus Code; 3rd column:
% Stimulus onset time; 4th column: Stimulus end time  
% clear SCode_vec STime_vec ETime_vec  % They're stored in TCSE all together
end