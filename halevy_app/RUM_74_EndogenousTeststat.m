function [nuhat,fval] =  RUM_74_EndogenousTeststat(A,pi_hat,N_j,tau,Omega) 
%% Code Description (CVX version for Halevy)
% nuhat solves min (A*nu - pihat)'Omega(A*nu - pihat) subject to nu >= 0
% using CVX.

warning('off','MATLAB:nargchk:deprecated')
%% Setup
[I,H]=size(A);
if nargin<5
    Omega=eye(I);
end
if nargin<4
    tau=0;
end

% Check to see if Omega is PD and symmetric
[~,p] = chol(Omega);
if p > 0
    error('Omega is not positive definite')
end
if Omega ~= Omega.'
    error('Omega is not symmetric')
end
% N = sum N_j
N = sum(N_j);

% Lower-bound for v
low=tau.*ones(H,1)/H;

cvx_begin quiet
    variable v(H)
    minimize( quad_form(A*v - pi_hat, Omega) )
    subject to
        v >= low;
cvx_end

if ~(strcmpi(cvx_status,'Solved') || strcmpi(cvx_status,'Inaccurate/Solved'))
    error('CVX failed to solve for nuhat (status: %s)',cvx_status);
end

%% Output
nuhat =A*v;
fval  =N*(A*v-pi_hat).'*Omega*(A*v-pi_hat);
fval  = round(fval,6,'decimals');

end
