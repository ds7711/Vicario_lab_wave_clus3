% unit_stats
% Get the comman statistics about a channel:
%       1. Contamination rate: upper band for a specific time range
%       2. Total number of spikes in the units: 
%       3. Firing rate: 
% input: 
%       spiketiming: a column vector of spike times
%       Optional input: 
%           cr_refractory: default is 2ms
% output: 
%       contamination_rate: 
%       tot_spikes: total number of spikes
%       firing_rate: 

function [contamination_rate, tot_spikes, firing_rate] = unit_stats(spiketiming, varargin)
if nargin > 1
    cr_refractory = varargin{1}; % if not specified, default refractory period is set to be 2ms
else 
    cr_refractory = 0.002;
end
tot_spikes = size(spiketiming, 1);
% firing rate in seconds (not precise because the last spike may come way earlier before the recording ends)
firing_rate = tot_spikes / max(spiketiming); 
isi = diff(spiketiming);
contamination_rate = sum(isi < cr_refractory) / (tot_spikes - 1);
end