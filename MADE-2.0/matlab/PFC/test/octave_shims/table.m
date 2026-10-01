function t = table(varargin)
% Octave test shim: a one-row "table" is a struct; table() is [].
if nargin == 0
    t = []; return
end
names = {};
vals = varargin;
for i = 1:numel(varargin)-1
    if ischar(varargin{i}) && strcmp(varargin{i}, 'VariableNames')
        names = varargin{i+1}; vals = varargin(1:i-1); break
    end
end
if numel(vals) ~= numel(names)
    error('table shim: %d values for %d variable names (MATLAB would reject this)', numel(vals), numel(names));
end
t = struct();
for i = 1:numel(names)
    t.(names{i}) = vals{i};
end
end
