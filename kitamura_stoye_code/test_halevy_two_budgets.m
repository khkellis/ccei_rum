%% Smoke test: Halevy data with two budgets (Observations 1 and 2)
% Runs the series estimator with polynomial_degree = 0 on two budgets to
% keep runtime small. Requires CVX for the X-generation feasibility checks.

clear

global classes budget_l budget_n start_year end_year periods poly_degree
global flag_genAX flag_estimator bootstrap_reps seed num_cores tau_val

% Minimal globals
classes        = 2;
budget_l       = 2;
budget_n       = 1;
start_year     = 1;
end_year       = budget_l;
periods        = budget_l;
poly_degree    = 0;
flag_genAX     = 1;
flag_estimator = 1;   % Series estimator (exogenous)
bootstrap_reps = 20;  % Small for a quick check
seed           = 12345;
num_cores      = 1;
tau_val        = 0;

% Load budgets/shares/income for Observation IDs 1 and 2
[budgets, shares, income, Z, obs_list] = RUM_21_budgets_halevy([1 2]);

% Generate X and A
X = RUM_31_genX(budgets);
A = RUM_41_genA(X);

% Compute test statistic
[pi_hat, Jstat, Jstat_bs, CV, pval, tau, eta_hat] = ...
    RUM_51_BasisStatistics(budgets, shares, income, X, A);

disp('Obs IDs used:'); disp(obs_list);
disp('Pi_hat (first 10 rows):'); disp(pi_hat(1:min(10,end)));
disp('J-stat:'); disp(Jstat);
disp('p-value (approx):'); disp(pval);
