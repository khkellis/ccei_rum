function RUM_11_cvx()
%% Solver setup for the Halevy application (CVX-only, no reinstall)
% Assumes CVX has already been installed via cvx_setup. We only call
% cvx_startup to add CVX to the path and never trigger installation logic.

warning('off','MATLAB:nargchk:deprecated')

% Try to locate CVX without invoking cvx_setup.
has_cvx = exist('cvx_startup','file') == 2 || exist('cvx_begin','file') == 2;

% If CVX is not on the path, add the bundled copy (../cvx) if present.
if ~has_cvx
    here = fileparts(mfilename('fullpath'));
    cvx_root = fullfile(here,'..','cvx');
    if exist(fullfile(cvx_root,'cvx_startup.m'),'file') == 2
        addpath(cvx_root);
        addpath(genpath(cvx_root));
        has_cvx = true;
    end
end

if ~has_cvx
    error(['CVX is required for the Halevy code. Please install CVX and ', ...
           'run cvx_setup once manually, then re-run this script.']);
end

% Initialize CVX without re-running cvx_setup (avoids reinstall prompts).
if exist('cvx_startup','file') == 2
    cvx_startup;
end

fprintf('Using existing CVX installation (cvx_startup only; no reinstall).\\n');
end
