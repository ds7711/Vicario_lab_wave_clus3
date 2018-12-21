% create new wavemark output variables 
function [num_sorted_sua, varargout] = new_wavemark_spike2(results, eventdata)

global flag; % global parameter for controling wave_clus

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

num_sorted_sua = -1;  % flag that the file cannot be clustered or no spikes detected

num_spikes = length(cluster_class);  % if this is the first time that current channel is sorted, write MUA

if flag.trial == 0  % only create the 
    if num_spikes / max_time < param.min_fr  % if total number of spikes detected is small, jump to next channel
        return;  % stop the function and return
    else
        % write multi-unit spikes using keyword mu
        chan_name = [param.mua_label, param.separator, num2str(iChan)];
        fillret = create_wavemark(fhand, cluster_class(:, 2), spikes, iChan, chan_name, param);
        if fillret == 0 % fillret == 0 means wavemark channel was successfully created
            flag.MUA = 1;
            flag.trial = flag.trial + 1;
        end
    end
    
else
    flag.trial = flag.trial + 1;
end


unq_classes = unique(cluster_class(:, 1));
num_classes = min(param.max_classes, length(unq_classes)); % only check the first few classes


% write single-unit spikes using keyword su
if num_classes
    num_sorted_sua = 0;  % number of SUA detected in the current clustering
    for kkk = 1 : num_classes  % 0 indicates units not assigned to any class, will not be considered
        
        tmp_class = unq_classes(kkk);
        if tmp_class ~= 0  % only processs non-zero class
             % get the logic index of one class of spikes
            temp_logic_idx_class = cluster_class(:, 1) == tmp_class;
            temp_spiketiming = cluster_class(temp_logic_idx_class, 2);
            temp_spikewave = spikes(temp_logic_idx_class, :);
            % test whether it is a good single-unit
            [contamination_rate, tot_spikes, ~] = unit_stats(temp_spiketiming, cr_refractory);

            % only write wavemark if statisfying all criterions
            if contamination_rate < param.contamination_rate && tot_spikes > LEAST_NUM_SPIKE && tot_spikes / max_time > param.min_fr
                % create the wavemark data
                chan_name = [param.sua_label, param.separator, num2str(iChan), param.separator, num2str(tmp_class)];
                fillret = create_wavemark(fhand, temp_spiketiming, temp_spikewave, iChan, chan_name, param);

                if fillret == 0
                    num_sorted_sua = num_sorted_sua + 1;
                    flag.SUA = flag.SUA + 1;
                end
            end
        end

        clearvars temp_logic_idx_class temp_spiketiming temp_spikewave tot_spikes firing_rate

    end
    if num_sorted_sua > 0  % if more than 1 SUA was detected, continue to next channel
        return;
    end
end

end
