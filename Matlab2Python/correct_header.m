function header = correct_header(header, stim_dict, dur_dict)
% correct the stim in header and return the new header
for i = 1 : length(header)
    header(i, 4) = stim_dict(header(i, 3));
    header(i, 6) = header(i, 6) + dur_dict(header(i, 3));
end
end