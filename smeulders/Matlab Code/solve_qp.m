function [x,fval,exitflag] = solve_qp(Q,f,Aineq,bineq,Aeq,beq,lb,ub)
%SOLVE_QP Thin wrapper over quadprog to mimic cplexqp signature.
if nargin < 8, ub = []; end
if nargin < 7, lb = []; end
if nargin < 6, beq = []; end
if nargin < 5, Aeq = []; end
if nargin < 4, bineq = []; end
if nargin < 3, Aineq = []; end

if exist('quadprog','file') ~= 2
    error('Optimization Toolbox (quadprog) not available.');
end

opts = optimoptions('quadprog','Display','none','Algorithm','interior-point-convex');
[x,fval,exitflag] = quadprog(Q,f,Aineq,bineq,Aeq,beq,lb,ub,[],opts);
end
