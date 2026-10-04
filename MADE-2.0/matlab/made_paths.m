function made_paths(family)
%MADE_PATHS Set the MATLAB path for one MADE 2.0 coil family.
%
%   MADE_PATHS(family) with family 'TF', 'CS' or 'PFC' removes every folder
%   of MADE-2.0/matlab left on the path by an earlier run (another family,
%   test shims, legacy) and adds:
%     - shared/   functions common to all families (axisymmetric field
%                 xbr/xbz/xlm, REBCO fits Ic_sst33/jc_ybco, hot-spot ODE,
%                 fatigue crack growth fcgr/ParametriY, axisymmetric FE
%                 solver and field cases);
%     - the folders of the family: TF = MADE-2.0/matlab and its subfolders
%       except CS, PFC, shared, legacy; CS = CS/ except CS/test; PFC = PFC/
%       except PFC/test.
%   The three families keep their own pipeline functions with the same
%   names (read_machine_input, scan_wp_designs, ...): only one family is on
%   the path at a time, so the wrong one can never shadow the right one.

root = fileparts(mfilename('fullpath'));
family = upper(family);
if ~any(strcmp(family, {'TF', 'CS', 'PFC'}))
    error('made_paths:family', 'family must be ''TF'', ''CS'' or ''PFC'' (got %s).', family);
end

% 1. remove every MADE folder already on the path (except root itself)
old = strsplit(path, pathsep);
old = old(strncmp(old, root, numel(root)) & ~strcmp(old, root));
if ~isempty(old), rmpath(old{:}); end

% 2. shared functions, then the family
add = [list(fullfile(root, 'shared'), {}), family_dirs(root, family)];
addpath(strjoin(add, pathsep));
addpath(root);
end

function d = family_dirs(root, family)
switch family
    case 'TF'
        d = list(root, {fullfile(root, 'CS'), fullfile(root, 'PFC'), fullfile(root, 'shared'), ...
            fullfile(root, 'legacy')});
        d = d(~strcmp(d, root));
    case 'CS'
        d = list(fullfile(root, 'CS'), {fullfile(root, 'CS', 'test')});
    case 'PFC'
        d = list(fullfile(root, 'PFC'), {fullfile(root, 'PFC', 'test')});
end
end

function d = list(top, excl)
% top and all its subfolders, minus the excluded subtrees
d = strsplit(genpath(top), pathsep);
d = d(~cellfun(@isempty, d));
for k = 1:numel(excl)
    d = d(~(strcmp(d, excl{k}) | strncmp(d, [excl{k} filesep], numel(excl{k}) + 1)));
end
end
