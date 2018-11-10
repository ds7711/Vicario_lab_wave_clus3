% add_gap_timing
% return the spike timing with gap duration included

function abs_cluster_class = add_gap_timing(gap_logic_idx, cluster_class, SAMPLING_RATE)
% get the change point of gap_logic_index, which reflectss the recording and gap start timing stamps; 
% add 1 before gap_logic_idx so that the real start timing points beomces a
% recording start point (recorded as -1)
change_point = diff([1, gap_logic_idx]); 
recording_sts = find(change_point == -1); % -1 mark the recording start index
gap_sts = find(change_point == 1); % 1 mark the gap starting index
% if the raw_wave_data ends with recording and is not followed with a gap,
% modify gap_st to get the same length
if length(recording_sts) > length(gap_sts) 
    gap_sts = [gap_sts, gap_sts(end) + max(diff(gap_sts))];  % if the data stop with recording, create and artificial gap starting index at the end
end

% recording duration in number of time stamps
recording_dur_insteps = gap_sts - recording_sts;
% recording_dur_insteps = recording_dur_insteps - 1; % real recording duration is 1 time stamp shorter
% gap duration in number of time stamps
gap_dur_insteps = recording_sts(2 : end) - gap_sts(1 : end - 1);
gap_dur_insteps = gap_dur_insteps + 1; % real gap duration is 1 time stamp longer

cum_recording_duration = cumsum(recording_dur_insteps) ./ SAMPLING_RATE;
cum_gap_duration = cumsum(gap_dur_insteps) ./ SAMPLING_RATE;
cum_gap_duration = [0, cum_gap_duration]; % no gap before the 1st recording block
num_spikes = size(cluster_class, 1); 
abs_cluster_class = cluster_class; % correct spike timing by adding gap duration back
for iii = 1 : num_spikes
    block_number_of_spike = find(cum_recording_duration > cluster_class(iii, 2), 1, 'first'); % the least required number of recording durations defines which recording block the spike belong to
    abs_cluster_class(iii, 2) = cluster_class(iii, 2) + cum_gap_duration(block_number_of_spike);
end
end