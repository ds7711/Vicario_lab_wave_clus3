function [stim_dict, dur_dict] = header2dict(bird_header)
% find the relationship between trial and stim id from the header.
% header should come from one single bird.
logidx = bird_header(:, 2) == bird_header(1, 2);
bird_header = bird_header(logidx, :);
dur = diff(bird_header(:, 5:6), 1, 2);
stim_dict = containers.Map(bird_header(:, 3),bird_header(:, 4));
dur_dict = containers.Map(bird_header(:, 3), dur);
end