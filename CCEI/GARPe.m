function pass = GARPe(q,p,e,fosd)
% GARP(e)    : Check GARP(e) condition
% Written by : Taisuke Imai
% 
% - Input
%   [q]    Quantities
%   [p]    Prices
%   [e]    Slack
%   [fosd] =1 if testing F-GARP (Nishimura et al. 2017; Polisson et al. 2015)
% - Output
%   [pass] 1: pass, 0: fail

% Size of the input matrix
T = size(q,1);  % Number of observations

% Mirror-image choices
qM = q(:,[2 1]);

% Directly revealed preference relations
R0 = eye(T);      % Element (i,j) = 1 if i is directly revealed preferred to j
P0 = zeros(T,T);  % Element (i,j) = 1 if i is directly strictly revealed preferred to j
for i = 1:T
    if fosd == 0
        R0(i,:) = (e*p(i,:)*q(i,:)'>=p(i,:)*q')';
        P0(i,:) = (e*p(i,:)*q(i,:)'>p(i,:)*q')';
    elseif fosd == 1
        R0(i,:) = ((e*p(i,:)*q(i,:)'>=p(i,:)*q')' | (e*p(i,:)*q(i,:)'>=p(i,:)*qM')');
        P0(i,:) = ((e*p(i,:)*q(i,:)'>p(i,:)*q')' | (e*p(i,:)*q(i,:)'>p(i,:)*qM')');
    end
end

% Revealed preference relations
% Construct transitive closure of R0 with Warshall's algorithm
R = R0;
for i = 1:T
    for j=1:T
        % Make changes only when R(i,j) = 1
        R(i,:) = max([R(i,:); R(i,j)*R(j,:)]);
    end
end

% Check GARP condition
pass = 1;
for i = 1:T
    for j = 1:T
        if R(i,j) == 1 && P0(j,i) == 1
            pass = 0;
            break;
        end
    end
    if pass == 0
        break;
    end
end

