function e = CCEI(q,p,fosd)
% CCEI       : Calculate Afriat's Critical Cost Efficiency Index using a binary search procedure
% Written by : Taisuke Imai
% 
% - Input
%   [q]    Quantities
%   [p]    Prices
%   [fosd] =1 if testing F-GARP (Nishimura et al. 2017; Polisson et al. 2015)
% - Output
%   [e] CCEI value (=1 if GARP is satisfied)

if GARP(q,p,fosd) == 1
    % CCEI = 1 if GARP/F-GARP is satisfied
    e = 1;
else
    eH  = 1;      % Upper bound
    eL  = 0;      % Lower bound
    tol = 1e-4;  % Tolerance for stopping criterion
    
    while eH-eL >= tol
        
        current_e = (eH+eL)/2;                  % Evaluate this value
        pass      = GARPe(q,p,current_e,fosd);  % Check if Garp(v) is satisfied
        
        if pass == 1
            % Update lower bound
            eL = current_e;
        elseif pass == 0
            % Update upper bound
            eH = current_e;
        end
    end
    e = current_e;
end

