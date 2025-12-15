%% Column generation main (Halevy data)
% Self-contained entrypoint that uses the Halevy adjusted dataset with the
% column-generation routines. Defaults to exogenous series estimator with
% polynomial_degree = 0 (income effectively exogenous).

clear
tic

warning('off','MATLAB:nargchk:deprecated');
warning('off','MATLAB:nearlySingularMatrix');

%% Paths
this_dir  = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
addpath(this_dir); % column-generation helpers + solver wrappers
addpath(fullfile(repo_root,'kitamura_stoye'));
addpath(fullfile(repo_root,'halevy'));
addpath(repo_root);

if exist('solve_lp','file') ~= 2
    error('Optimization wrappers (solve_lp/solve_qp/solve_bilp) not on path. Current path: %s', path);
end

%% User parameters
budget_length     = [];     % use all budgets by default (no rolling windows)
number_classes    = 2;      % Halevy has two goods
polynomial_degree = 0;      % 0 because income is effectively exogenous
estimator         = 1;      % 1 = Basis (exogenous), 2 = Endogenous
genAX             = 1;      % regenerate X and Tau_Set
tau_set_size      = 500;    % patterns used for tightening
indJ              = 1;      % compute bootstrap
cores             = 10;      % parallel workers
bs_reps           = 1000;    % bootstrap reps (keep modest for Halevy)
tau_ind           = 0;      % tau=0 for exogenous case

%% Globals
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores tau_val

%% Flags
flag_estimator = estimator;
tau_val        = tau_ind;

flag_genAX = genAX;
nameAX = strcat('Input/RUM_',num2str(budget_length),'LengthBudget_',num2str(number_classes),'.mat');
if exist(nameAX) ~= 2 && genAX == 0
    flag_genAX = 1;
end

%% Parameters
budget_l  = budget_length;
classes   = number_classes;
poly_degree = polynomial_degree;
bootstrap_reps = bs_reps;
if length(budget_l) == 0
    seed = 10000*(bootstrap_reps/1000 + 1) + poly_degree*1000  + classes*10 + flag_estimator;
else
    seed = 10000*(bootstrap_reps/1000 + 1) + poly_degree*1000  + budget_l*100 + classes*10 + flag_estimator;
end

if feature('numCores') < cores
    num_cores = feature('numCores');
else
    num_cores = cores;
end

%% Parallel pool
v_=ver;
[installedToolboxes{1:length(v_)}] = deal(v_.Name);
Parellel = all(ismember('Parallel Computing Toolbox',installedToolboxes));
if Parellel == 1
    delete(gcp('nocreate'));
    parpool(num_cores);
end

%% Load Halevy budgets
[budgets_all,shares_all,income_all,instrument_all] = RUM_21_budgets_halevy;

periods    = size(budgets_all,1);
if isempty(budget_length)
    budget_length = periods; % no rolling windows for Halevy
end
start_year = 1;             % placeholders; Halevy data are not time-indexed
end_year   = budget_length;
budget_n   = 1;
budget_l   = budget_length;
budgets    = {budgets_all}; % single window over all budgets

%% Compute X and pattern sets
fprintf('Generating X and patterns using all %d budgets (no rolling windows)\n',budget_l);
X = cell(1,1);
Tau_Set = cell(1,1);
X{1} = RUMCG_31_genX(budgets{1});
Tau_Set{1} = RUMCG_PatternSet(tau_set_size,X{1});

%% Test statistic and bootstrap
if indJ
    share_slice = shares_all(1:budget_l);
    income_slice = income_all(1:budget_l);
    inst_slice = instrument_all(1:budget_l);
    if estimator == 1
        [pi_hat,Jstat,Jstat_bs,CV,pval,tau,eta_hat] = ...
            RUMCG_51_BasisStatistics(budgets{1},share_slice,income_slice,X{1},Tau_Set{1});
    elseif estimator == 2
        [pi_hat,Jstat,Jstat_bs,CV,pval,tau,eta_hat] = ...
            RUMCG_71_EndogenousStatistic(budgets{1},share_slice,income_slice,inst_slice,X{1},Tau_Set{1});
    else
        error('Estimator %d not supported in column-generation code.', estimator);
    end
    fprintf('Completed bootstrap for full Halevy sample.\n');
end

%% Save
Time_taken = toc;
c = clock;
datetime_str = sprintf('%d%02d%02d_%02d%02d',c(1),c(2),c(3),c(4),c(5));
name = sprintf('../../halevy/Halevy_RUMCG_len%d_classes%d_poly%d_bs%d_est%d_tau%d_%s', ...
    budget_l, classes, poly_degree, bootstrap_reps, estimator, tau_ind, datetime_str);
save(name);
fprintf('Saved results to %s (elapsed %.2f seconds)\n', name, Time_taken);
