%% Solver setup (CVX optional)
% Legacy hook kept for compatibility. We now prefer MATLAB's Optimization
% Toolbox (`linprog`/`quadprog`) and only touch CVX if it is already on the
% path. No license files or bundled binaries are needed anymore.

% If Optimization Toolbox is present, just announce and return.
has_opt_toolbox = exist('linprog','file') == 2 && exist('quadprog','file') == 2;
if has_opt_toolbox
    fprintf('Optimization Toolbox detected (linprog/quadprog); skipping CVX setup.\\n');
    return
end

% If CVX is already installed on the MATLAB path, initialize it quietly.
if exist('cvx_setup','file') == 2
    cvx_setup;
    if exist('cvx_startup','file') == 2
        cvx_startup;
    end
    fprintf('CVX initialized because Optimization Toolbox was not detected.\\n');
    return
end

warning(['Neither Optimization Toolbox nor CVX was detected on the path. ', ...
         'Please install the Optimization Toolbox (preferred) or add an existing CVX installation if you want to run the code.']);
