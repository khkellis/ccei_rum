function [ drop_row ] = RUM_32_genXdroprow(budgets,X,dimx,ii,J)
%% Code Description: RUM_32_genXdroprow (CVX version for Halevy)
% Determines whether to drop a particular row of X using a CVX feasibility
% check: find x >= 0 such that A*x <= b. No Optimization Toolbox calls.

warning('off','MATLAB:nargchk:deprecated')

% Tolerances to mitigate floating-point issues in binding/strict checks.
tol_eq = 1e-9;      % For equality: |p*x - 1| <= tol_eq
tol_strict = 1e-6;  % For strict above/below: buffer away from 1

A = [];
b = [];
for tt=1:J %%translate (0,-1,1) into linear inequalities
    if X(ii,tt)==0
        % |p*x - 1| <= tol_eq
        A=[A;budgets(tt,:);-budgets(tt,:)];
        b=[b;1+tol_eq;-(1-tol_eq)];
    elseif X(ii,tt)==1
        % p*x >= 1 + tol_strict
        A=[A;-budgets(tt,:)];
        b=[b;-(1+tol_strict)];
    else
        % p*x <= 1 - tol_strict
        A=[A;budgets(tt,:)];
        b=[b;1-tol_strict];
    end
end
A=[A;-eye(dimx)];
b=[b;zeros(dimx,1)];

% Feasibility check via CVX.
Aineq = A;
bineq = b;

cvx_begin quiet
    variable x(dimx)
    minimize(0)
    Aineq * x <= bineq;
    x >= 0;
cvx_end

is_feasible = strcmpi(cvx_status,'Solved') || strcmpi(cvx_status,'Inaccurate/Solved');
drop_row = ~is_feasible;
end
