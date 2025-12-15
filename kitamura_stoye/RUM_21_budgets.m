function [budgets,shares,income,Z,obs_list] = RUM_21_budgets(obs_subset)
%% Halevy loader wrapper
% Use the Halevy et al. (2016) adjusted data instead of the UK survey
% files. This is a thin wrapper around halevy/RUM_21_budgets_halevy.m so
% the rest of the pipeline can keep calling RUM_21_budgets.

%% Globals (only used for a consistency check)
global classes

if nargin < 1
    obs_subset = [];
end

% Add the Halevy helper to the path and call it
this_dir  = fileparts(mfilename('fullpath'));
halevy_dir = fullfile(this_dir,'..','halevy');
addpath(halevy_dir);

[budgets,shares,income,Z,obs_list] = RUM_21_budgets_halevy(obs_subset);

% Sanity check: enforce two goods unless the caller explicitly overrides
if isempty(classes)
    classes = size(budgets,2);
elseif classes ~= size(budgets,2)
    error('classes=%d but Halevy data has %d goods. Set number_classes=2.',classes,size(budgets,2));
end

end
