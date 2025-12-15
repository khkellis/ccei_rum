function setup_cplex(cplex_root)
%SETUP_CPLEX Add the CPLEX MATLAB API to the path.
%   Call without arguments to auto-detect a CPLEX_Studio installation using
%   common environment variables (CPLEX_STUDIO_DIR, CPLEXROOT) and default
%   install locations. Optionally pass the root install directory
%   explicitly, e.g. setup_cplex('/Applications/CPLEX_Studio2211').
%
%   This is a convenience helper; it leaves the path untouched if CPLEX
%   cannot be found.

if nargin < 1 || isempty(cplex_root)
    cplex_root = first_existing({ ...
        getenv('CPLEX_STUDIO_DIR'), ...
        getenv('CPLEX_STUDIO_DIR128'), ...
        getenv('CPLEXROOT'), ...
        '/Applications/CPLEX_Studio2211', ...
        '/Applications/CPLEX_Studio201', ...
        '/Applications/CPLEX_Studio2010', ...
        '/Applications/CPLEX_Studio1210', ...
        '/opt/ibm/ILOG/CPLEX_Studio2211', ...
        '/opt/ibm/ILOG/CPLEX_Studio201', ...
        '/Users/Keaton/Applications/CPLEX_Studio_Community2212'
    });
end

if isempty(cplex_root) || ~isfolder(cplex_root)
    warning('setup_cplex:missing','CPLEX root not found. Set CPLEX_STUDIO_DIR or pass the path to setup_cplex.');
    return;
end

if ismac
    api_path = fullfile(cplex_root,'cplex','matlab','x86-64_osx');
elseif ispc
    api_path = fullfile(cplex_root,'cplex','matlab','x64_win64');
else
    api_path = fullfile(cplex_root,'cplex','matlab','x86-64_linux');
end

if ~isfolder(api_path)
    warning('setup_cplex:apipath','CPLEX MATLAB API folder not found at %s',api_path);
    return;
end

addpath(api_path);

if exist('cplexlp','file') == 2
    disp(['CPLEX MATLAB API added to path: ', api_path]);
else
    warning('setup_cplex:nocplex','API path added but cplex functions are not visible. Check your installation.');
end

end

function path_out = first_existing(candidates)
path_out = '';
for ii = 1:numel(candidates)
    c = candidates{ii};
    if isempty(c)
        continue;
    end
    if isfolder(c)
        path_out = c;
        return;
    end
end
end
