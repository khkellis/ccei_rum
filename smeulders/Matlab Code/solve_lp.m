function [x,fval,exitflag] = solve_lp(f,Aineq,bineq,Aeq,beq,lb,ub)
%SOLVE_LP Thin wrapper over linprog to mimic cplexlp signature.
if nargin < 8, ub = []; end
if nargin < 7, lb = []; end
if nargin < 6, beq = []; end
if nargin < 5, Aeq = []; end
if nargin < 4, bineq = []; end
if nargin < 3, Aineq = []; end

if exist('linprog','file') ~= 2
    error('Optimization Toolbox (linprog) not available.');
end

opts = optimoptions('linprog','Display','none','Algorithm','dual-simplex');
[x,fval,exitflag] = linprog(f,Aineq,bineq,Aeq,beq,lb,ub,opts);
end
