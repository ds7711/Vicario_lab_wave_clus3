% ISI distribution 
% input: 
%       YYY: binary representation of spike trains 
%       resolution: bin size used to enumerate spikes (most data is binned at 0.0001 second)
% output: 
%       isi: a vector storing all the inter-spike-interval values
function isi = spike_isi(YYY, resolution)
num_trial = size(YYY, 1);
isi = [];
for iii = 1 : num_trial
    tmp_idx = find(YYY(iii, :) == 1); % 1 represents spikes, index of spikes in iTH trial
    isi = [isi, diff(tmp_idx)];
end
isi = isi * resolution; 
end