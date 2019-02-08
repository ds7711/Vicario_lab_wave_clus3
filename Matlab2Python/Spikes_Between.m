% Spikes_Between take the absolute spiketiming vector and user specified time
% range as input; 
% The input time_range should be a matrix with its 1st column storing the starting
% time poing and 2nd column storing end time point. 
% The input spiketiming is just absolute spike timing during the whole
% recording session (should be a column vector). 
% The output:
%       1: Mean firing rate for each trial
%       2: is a cell which has the same number of elements as the rows of
%          time range
function [FR, Spikes] = Spikes_Between(spiketiming,time_range)
% num_spikes=length(spiketiming);
num_ranges=size(time_range,1);
Spikes=cell(num_ranges,1);
FR = zeros(num_ranges, 1);
clear ST ET; 
for i=1:num_ranges
    ST=time_range(i,1);
    ET=time_range(i,2);
%     Index=find(spiketiming(:,1)>ST&spiketiming(:,1)<ET);
%     Spikes{i,1}=spiketiming(Index);
    Spikes{i,1}=spiketiming(spiketiming(:,1)>ST & spiketiming(:,1)<ET); 
    FR(i, 1) = length(Spikes{i,1}) / (ET - ST);
end
end