function [x,fval,exitflag] = solve_bilp(f,Aineq,bineq,Aeq,beq)
%SOLVE_BILP Binary ILP wrapper over intlinprog (all binaries).
if nargin < 5, beq = []; end
if nargin < 4, Aeq = []; end
if nargin < 3, bineq = []; end
if nargin < 2, Aineq = []; end

if exist('intlinprog','file') ~= 2
    error('Optimization Toolbox (intlinprog) not available.');
end

n = numel(f);
intcon = 1:n;
lb = zeros(n,1);
ub = ones(n,1);
opts = optimoptions('intlinprog','Display','none');
[x,fval,exitflag] = intlinprog(f,intcon,Aineq,bineq,Aeq,beq,lb,ub,opts);
end
