function [nuhat,fval] =  RUM_54_BasisTeststat(A,pi_hat,N,poly_degree,tau,Omega) 
%% Code Description (CVX version for Halevy)
% nuhat solves min (A*nu - pihat)'Omega(A*nu - pihat) subject to nu >= 0
% using CVX. Avoids quadprog/Optimization Toolbox entirely.

warning('off','MATLAB:nargchk:deprecated')

%% Setup
[I,H]=size(A);
if nargin<6
    Omega=eye(I);
end
if nargin<5
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
% N = sum N_j/poly_degree
% We rescale N by poly_degree as per series estimation section in KS
% That is, we replace the vactor N with N/max_jK(j), where K(j) is the
% polynomial degree in year j.  This needs to be adjusted if one uses a
% different polynomial degree for each year, in which case we use the
% largest polynomial degree.
denom = max(poly_degree,1); % avoid divide-by-zero when poly_degree = 0 in Halevy runs
N = sum(N)/denom;

% Lower-bound for v
low=tau.*ones(H,1)/H;

% CVX quadratic program
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
