% calculate logic index based on the multi-dimension label and conditions
% input: 
%       data_label: M x N matrix, each column corresponds to a different
%           feature, each row corresponds to one sample
%       varargin: 
%               conditions stored in varargin one-by-one, from left to right. 
%           single number: the value in data label has to equal to it
%           NaN: skip this feature
%           string: in the format of '< X', '> X', '== X', '<= X', '>= X'. 
% output: 
%       logidx: 
%           the logic index that of labels that satify all specified
%       conditions (multi-conditions are evaluated by using "&" operation)
%       
function logidx = logic_idx(data_label, varargin)
% fprintf('Number of arguments: %d\n',nargin) % for debugging
logidx = ones(size(data_label, 1), 1); % initialize the logidx
num_conditions = nargin - 1; % number of specified conditions, data_label is not a condition
if num_conditions > size(data_label, 2) 
    logidx = logidx * 0;
    errordlg('Number of conditions > number of features!!!') % Error: too many conditions
end
for iii = 1 : num_conditions
    para = varargin{iii}; % obtain the condition
    if isnan(para) % if NaN, skip this feature
        continue;
    elseif ischar(para) % if string, use 'eval' to evaluate
        logidx = logidx .* (eval(['data_label(:, iii)', para]));
    else % if a number, use "equal" operation
        logidx = logidx .* (data_label(:, iii) == para);
    end
end
logidx = 1 == logidx;
end
