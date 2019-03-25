function chan_num = chanstr2double(chanstr, param)
% convert a channel string to a number.
sua_num_c = param.sua_num_c;

if ~contains(chanstr, '_')
% if MUA, no '_' in between,
% directly convert to an integer
chan_num = str2double(chanstr(isstrprop(chanstr, 'digit'))); % channel #
else
% if SUA, with '_' in between,
% convert to an integer
chan_units = strsplit(chanstr, '_');
chan_num = str2double(chan_units{1,1}) + str2double(chan_units{1,2}) / sua_num_c;
chan_num = chan_num + 1.0 / sua_num_c;  % add a small constant to avoid 0.
end
end