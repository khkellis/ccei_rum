%% Code Description: Main file for Halevy et al. (2016) application
% Runs the Kitamura-Stoye RUM test on the experimental data with two goods.
% Allows selecting a subset of budgets (Observation IDs) for lighter runs.

clear
tic

warning('off','MATLAB:nargchk:deprecated')

% Add path to core Kitamura-Stoye code
addpath(genpath('../kitamura_stoye_code'))

%% User-specified parameters
obs_subset        = [];  % e.g., [1 2] for a two-budget smoke test; [] = all 22
polynomial_degree = 0;   % 0 gives constant series estimator (income is 1)
estimator         = 1;   % 0 = Kernel; 1 = Exogenous Series; 2 = Endogenous Series
genAX             = 1;   % Regenerate X/A (set to 0 to load cached matrices)
indJ              = 1;   % Run bootstrap
cores             = 10;   % Number of cores for parallel computing
bs_reps           = 200; % Bootstrap repetitions (keep small for laptop runs)
tau_ind           = 0;   % Set equal to 0/1 to set tau equal to 0/not 0.

%% Global parameters
global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores tau_val

%% Flags and sizes
classes        = 2;              % Two goods: X and Y
flag_estimator = estimator;
tau_val        = tau_ind;
flag_genAX     = genAX;
poly_degree    = polynomial_degree;
bootstrap_reps = max(100,bs_reps); % keep positive

% Load budgets/shares/income from Halevy data
[budgets_all,shares_all,income_all,instrument_all,obs_list] = RUM_21_budgets_halevy(obs_subset);
budget_l    = numel(obs_list);  % Number of budgets used
budget_n    = 1;                % Single window (all selected budgets together)
start_year  = 1;
end_year    = budget_l;
periods     = budget_l;

% Seed for reproducibility
seed = 10000*(bootstrap_reps/1000) + poly_degree*1000  + budget_l*100 + classes*10 + flag_estimator;

% Pre-specify cores
if feature('numCores') < cores
    num_cores = feature('numCores');
else
    num_cores = cores;
end

% Parallel computing setup if available
v_=ver;
[installedToolboxes{1:length(v_)}] = deal(v_.Name);
Parellel = all(ismember('Parallel Computing Toolbox',installedToolboxes));
if Parellel == 1
    delete(gcp);
    parpool(num_cores);
end

%% Solver check (CVX only)
RUM_11_cvx;

%% Map budgets to cell for downstream code
budgets = cell(1, budget_n);
budgets{1} = budgets_all;

%% Compute X and A
nameAX= strcat('Input/RUM_',num2str(budget_l),'LengthBudget_',num2str(classes),'.mat');
if exist(nameAX) ~= 2 && genAX == 0
    flag_genAX = 1;
end

if flag_genAX == 0
    load(nameAX)
else
    X{budget_n} = [];
    A{budget_n} = [];
    fprintf('Starting to generate X and A using %d budgets \\n',budget_l)
    for ii = 1:budget_n
        X{ii} = RUM_31_genX(budgets{ii});
        A{ii} = RUM_41_genA(X{ii});
        fprintf('Completed budget %d out of %d. \\n',ii,budget_n);
    end
    save(nameAX,'X','A')
end
Time_taken_AX = toc;

%% Pihat, Jstat, Critical Value, Pr(Jstat = 0)
if indJ
pi_hat{budget_n} = [];
eta_hat{budget_n} = [];
Jstat_bs{budget_n} = [];
Jstat   = zeros(budget_n,1);
CV      = zeros(budget_n,2);
pval    = zeros(budget_n,1);
tau     = zeros(budget_n,1);

if flag_estimator == 0
    for ii = 1:budget_n
        [pi_hat{ii},Jstat(ii,1),Jstat_bs{ii},CV(ii,:), pval(ii,1)] ...
            =  RUM_61_KernelStatistics(budgets{ii},cell(shares_all),cell(income_all),X{ii},A{ii});
        fprintf('Completed bootstrap %d out of %d. \\n',ii,budget_n);
    end
elseif flag_estimator == 1
    for ii = 1:budget_n
        [pi_hat{ii},Jstat(ii,1),Jstat_bs{ii},CV(ii,:), pval(ii,1),tau(ii,1),eta_hat{ii}] ...
            =  RUM_51_BasisStatistics(budgets{ii},cell(shares_all),cell(income_all),X{ii},A{ii});
        fprintf('Completed bootstrap %d out of %d. \\n',ii,budget_n);
    end
elseif flag_estimator == 2
    for ii = 1:budget_n
        [pi_hat{ii},Jstat(ii,1),Jstat_bs{ii},CV(ii,:), pval(ii,1), tau(ii,1),eta_hat{ii}] ...
            =  RUM_71_EndogenousStatistic(budgets{ii},cell(shares_all),cell(income_all),cell(instrument_all),X{ii},A{ii});
        fprintf('Completed bootstrap %d out of %d. \\n',ii,budget_n);
    end
end
end

%% Save
Time_taken = toc;
c = clock;
datetime = strcat(num2str(c(1)),num2str(c(2)),num2str(c(3)),'_',num2str(c(4)),num2str(c(5)));
prefix = 'Halevy_';
if flag_estimator == 1
     name = strcat('Output/',prefix,'RUM_',num2str(budget_l),'Budgets_',num2str(classes),'Classes_',num2str(poly_degree),'PolyDegree_',num2str(bootstrap_reps),'BSreps_Series_',num2str(tau_ind),'taunot0_',datetime);
elseif flag_estimator == 0
     name = strcat('Output/',prefix,'RUM_',num2str(budget_l),'Budgets_',num2str(classes),'Classes_',num2str(bootstrap_reps),'BSreps_Kernel_',num2str(tau_ind),'taunot0_',datetime);
else
     name = strcat('Output/',prefix,'RUM_',num2str(budget_l),'Budgets_',num2str(classes),'Classes_',num2str(poly_degree),'PolyDegree_',num2str(bootstrap_reps),'BSreps_Series_Endogenous',num2str(tau_ind),'taunot0_',datetime);
end
save(name);
