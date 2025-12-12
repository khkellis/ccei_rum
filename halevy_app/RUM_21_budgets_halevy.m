function [budgets,shares,income,Z,obs_list] = RUM_21_budgets_halevy(obs_subset)
%% Code Description: Halevy et al. (2016) loader (2 goods, unit income)
% Builds budgets, expenditure shares, income, and instrument placeholders
% from the experimental data in ../data/Halevy et al (2016) - Data.csv.
% - Prices are recovered from intercepts: p_x = 1/X_intercept,
%   p_y = 1/Y_intercept (income normalized to 1).
% - Shares are expenditure shares on the two goods for each subject.
% - Income is total expenditure per subject (should be 1 up to rounding).
% Inputs:
%   obs_subset (optional): vector of Observation IDs to keep (e.g., [1 2]).
%       If empty/omitted, uses all observations (22 budgets).
% Outputs:
%   budgets: [num_budgets x 2] price matrix, median-expenditure normalized
%   shares:  {num_budgets} cell, each [N_t x 2] expenditure shares
%   income:  {num_budgets} cell, each [N_t x 1] total expenditure
%   Z:       {num_budgets} cell, instrument placeholder (empty; exogenous)
%   obs_list: vector of observation IDs used (order matches budgets)
%
% Note: assumes two goods (classes = 2) and exogenous prices; intended for
% use with the series (estimator=1) or kernel (estimator=0) estimators.

%% Globals (kept for consistency with the original loader)
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores

if nargin < 1
    obs_subset = [];
end

% Read data (path is relative to halevy_app/)
data = readtable('../data/Halevy_Data_adjusted.csv');

% Trim to requested observation IDs (budgets)
if isempty(obs_subset)
    obs_list = unique(data.Observation);
else
    obs_list = obs_subset(:).';
    data = data(ismember(data.Observation, obs_list), :);
end

num_budgets = numel(obs_list);
budgets = zeros(num_budgets, 2);
shares = cell(num_budgets, 1);
income = cell(num_budgets, 1);
Z = cell(num_budgets, 1);

for ii = 1:num_budgets
    obs_id = obs_list(ii);
    d = data(data.Observation == obs_id, :);

    % Recover prices from intercepts (income = 1)
    px = median(1 ./ d.("X_intercept"));
    py = median(1 ./ d.("Y_intercept"));

    % Expenditure on each good and total expenditure
    exp_x = px .* d.X;
    exp_y = py .* d.Y;
    tot = exp_x + exp_y;

    % Shares and income for this budget
    shares{ii} = [exp_x ./ tot, exp_y ./ tot];
    income{ii} = tot;
    Z{ii} = []; % Exogenous prices: instrument not used

    % Normalize prices by median expenditure (tot ~ 1 here)
    budgets(ii, :) = [px, py] ./ median(tot);
end

end
