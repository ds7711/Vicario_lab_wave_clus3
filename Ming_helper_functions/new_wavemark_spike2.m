% create new wavemark output variables 
function [num_classes, varargout] = new_wavemark_spike2(results, eventdata)

% unpack files
param = eventdata.param;
fhand = eventdata.fhand;
iChan = eventdata.iChan;

LEAST_NUM_SPIKE = param.LEAST_NUM_SPIKE;
cr_refractory = param.cr_refractory;

timebase = CEDS64TimeBase(fhand);
max_time = CEDS64MaxTime( fhand ) * timebase;

spikes = results.spikes;
times = results.times; % in seconds
classes = results.classes;
% combine classes and times information into cluster_classes
cluster_class = [classes; times]';
clearvars results times classes;
% start calculation

num_spikes = length(cluster_class);
if num_spikes / max_time < param.min_fr  % if total number of spikes detected is small, jump to next channel
    return;  % stop the function and return
else
    % write multi-unit spikes using keyword mu
    chan_name = [param.mua_label, param.separator, num2str(iChan)];
    create_wavemark(fhand, cluster_class(:, 2), spikes, iChan, chan_name, param);
end

num_classes = min(param.max_classes, length(unique(cluster_class(:, 1)))); % only check the first few classes


% write single-unit spikes using keyword su
if num_classes

    num_sorted_sua = 0;  % number of SUA detected in the current clustering
    for kkk = 1 : (num_classes-1)  % 0 indicates units not assigned to any class, will not be considered

         % get the logic index of one class of spikes
        temp_logic_idx_class = cluster_class(:, 1) == kkk;
        temp_spiketiming = cluster_class(temp_logic_idx_class, 2);
        temp_spikewave = spikes(temp_logic_idx_class, :);
        % test whether it is a good single-unit
        [contamination_rate, tot_spikes, ~] = unit_stats(temp_spiketiming, cr_refractory);

        % only write wavemark if statisfying all criterions
        if contamination_rate < param.contamination_rate && tot_spikes > LEAST_NUM_SPIKE && tot_spikes / max_time > param.min_fr
            % create the wavemark data
            chan_name = [param.sua_label, param.separator, num2str(iChan), param.separator, num2str(kkk)];
            fillret = create_wavemark(fhand, temp_spiketiming, temp_spikewave, iChan, chan_name, param);
            
            if fillret > 0
                num_sorted_sua = num_sorted_sua + 1;
            end
        end

        clearvars temp_logic_idx_class temp_spiketiming temp_spikewave tot_spikes firing_rate

    end
    if num_sorted_sua > 0  % if more than 1 SUA was detected, continue to next channel
        return;
    end
end

end
