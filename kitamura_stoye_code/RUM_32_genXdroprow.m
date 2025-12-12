function [ drop_row ] = RUM_32_genXdroprow(budgets,X,dimx,ii,J)
%% Code Description: RUM_33_genXdroprow
% This code determines whether or not we drop a particular row of X.  Due
% to problems with parfor, genX module needs to call this one rather than
% embed it inside genX module.

warning('off','MATLAB:nargchk:deprecated')

% Tolerance to handle floating-point issues in tight/loose budget checks.
tol_eq = 1e-9;      % For equality (binding) constraints: |p*x - 1| <= tol_eq
tol_strict = 1e-6;  % For strict inequalities: enforce a buffer away from 1

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
% Feasibility check: Ax <= b with x >= 0 using linprog (CVX removed).
f = zeros(dimx,1);
Aineq = A;
bineq = b;
lb = zeros(dimx,1);
try
    opts = optimoptions('linprog','Display','none');
catch
    opts = [];
end
[~,~,exitflag] = linprog(f,Aineq,bineq,[],[],lb,[],opts);
drop_row = ~(exitflag == 1 || exitflag == 2);
end
