# Kitamura–Stoye RUM Testing Code

MATLAB implementation of the bootstrap tests for random utility models (RUM) and GARP from Kitamura & Stoye (Econometrica 2018). The code builds budget sets from UK Family Expenditure Survey data (1975–1999), constructs GARP-consistent pseudo-agents, and bootstraps test statistics under kernel, series, or endogenous series estimators.

## Layout
- `kitamura_stoye_code/` MATLAB source plus `Input/` (data, precomputed matrices) and `Output/` (saved runs).
- `theory/` Paper PDFs and notes that explain the econometric procedure.
- `data/` Additional external data (not used directly by the MATLAB scripts).
- `test.m` Scratch script opened in the IDE.

## Requirements
- MATLAB with Statistics and Optimization Toolbox (`linprog`/`quadprog`); Parallel Computing Toolbox is optional.
- Input CSVs named `data75.csv` … `data99.csv` in `kitamura_stoye_code/Input/` plus the price file `pall.csv`. Prebuilt `RUM_*LengthBudget_*.mat` files speed up reruns of the X/A construction.

## How to run
1) Open MATLAB with the working directory set to `kitamura_stoye_code/`.  
2) Edit the user parameters at the top of `kitamura_stoye_code/RUM_00_Main.m`:
   - `budget_length`: rolling window size (e.g., 7).
   - `number_classes`: number of commodity bundles (3–5 supported; adjust `RUM_21_budgets.m` to change aggregation).
   - `polynomial_degree`: basis degree for series estimators.
   - `estimator`: `0` kernel, `1` exogenous series, `2` endogenous series.
   - `genAX`: set `1` to (re)generate patches `X` and pseudo-agents `A`; leave `0` to reuse saved matrices in `Input/`.
   - `indJ`: set `1` to compute the bootstrap test statistic.
   - `cores`, `bs_reps`, `tau_ind`: compute resources and bootstrap tuning.
3) Run `RUM_00_Main`. Outputs (bootstrap draws, p-values, settings) are saved to a timestamped file in `kitamura_stoye_code/Output/`. Existing `Output/Tau*` folders contain prior runs.

Computation is intensive: large `budget_length`, higher `number_classes`, and high `bs_reps` can take many hours. Parallelization is used when available.

### Halevy et al. (2016) data (2 goods, unit income)
- Separate folder: see `halevy_app/` (keeps the application distinct from the original code).
- Solvers: Halevy scripts now rely solely on CVX. Run `cvx_setup` once to install CVX, then `RUM_11_cvx` only calls `cvx_startup` (no reinstall).
- Loader: `halevy_app/RUM_21_budgets_halevy.m` ingests `data/Halevy et al (2016) - Data.csv`. `obs_subset` lets you keep a subset of budgets (e.g., `[1 2]`).
- Main: `halevy_app/RUM_00_Main_halevy.m` sets `classes = 2`, `polynomial_degree = 0`, and defaults to the series estimator. Set `obs_subset` inside the script for the budgets you want; run the script to generate X/A, bootstrap, and save results under `Output/` with a `Halevy_` prefix.
- Quick smoke test: `halevy_app/test_halevy_two_budgets.m` runs the series estimator on budgets 1 and 2 with a tiny bootstrap to validate the pipeline.

## What each file does (MATLAB)
- `RUM_00_Main.m` Orchestrates the entire workflow: checks solver availability, builds budgets, constructs patches `X` and pseudo-agents `A`, then computes/bootstraps the test statistic.
- `RUM_11_cvx.m` Legacy name; now simply checks for Optimization Toolbox (preferred) and only initializes CVX if it is already installed.
- `RUM_21_budgets.m` Loads yearly FES microdata, filters to couples with one car and children, trims extremes, builds expenditure shares and normalized budgets, and outputs the instrument (`Z`) for endogenous estimation.
- `RUM_31_genX.m` Generates all feasible patches for a budget set; `RUM_32_genXdroprow.m` drops infeasible patches via `linprog` feasibility checks.
- `RUM_41_genA.m` Builds the matrix of GARP-consistent pseudo-agents; `RUM_42_genAup.m` and `RUM_43_genAfw.m` implement the tree search and cycle detection.
- Basis (series) estimator: `RUM_51_BasisStatistics.m` (main wrapper), `RUM_52_BasisRandomize.m` (bootstrap resampling), `RUM_53_BasisPihat.m` (estimate patch choice probabilities), `RUM_54_BasisTeststat.m` (compute tightened J-statistic).
- Kernel estimator: `RUM_61_*` analogues for kernel smoothing (`Statistics`, `Randomize`, `Pihat`, `Teststat`).
- Endogenous estimator: `RUM_71_*` analogues that use the instrument `Z`.
- `correct_pvals.m` Rounds bootstrap draws and recomputes p-values; helpful when J-statistic values are effectively zero.
- `Input/` Data and cached X/A matrices; `Output/` saved bootstrap results; `cvx/` legacy CVX binaries (no longer needed for normal runs).

## Notes
- Changing the price series or budget construction requires regenerating `X` and `A` (`genAX = 1`).
- The code assumes data years 1975–1999; adjust `start_year`/`end_year` in `RUM_00_Main.m` if you add data.
- The `theory/` PDFs describe the econometric test; keep them handy to interpret the outputs.
