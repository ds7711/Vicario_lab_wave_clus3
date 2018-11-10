% gap_detection
% input: continous waveform vector (dimension: 1 x N)
% output: 
%   gap_label: if gap exists, 1; if no, 0
%   logical index of continous gaps

function [gap_label, gap_logic_idx] = gap_detection(data_withzero)
gap_detection_threshold = 15; % if number of continous zeros is more than this threshold, say there's gap in the data
idx_all_zeros = ~data_withzero;  % logical index of all zeros (real indices)
max_continous_zero = 15; % The duration of the gap has to > max_continous_zero / samplingrate
conv_vector = ones(1, max_continous_zero); 
conv_idxzeros = conv(single(idx_all_zeros), conv_vector, 'same'); % only continous 1s could consistently generate big running sum
gap_logic_idx = conv_idxzeros >= max_continous_zero / 2; % if half of the values inside window is 1, then that zero belongs to a gap
gap_label = gap_logic_idx * ones(size(gap_logic_idx')) > gap_detection_threshold; 