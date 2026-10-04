here = fileparts(mfilename('fullpath')); cs = fileparts(here);
addpath(fileparts(cs)); made_paths('CS'); addpath(fullfile(here, 'octave_shims'));
fid = fopen(fullfile(here,'params_template.txt')); p = struct();
while true, l = fgetl(fid); if ~ischar(l), break; end, t = strsplit(strtrim(l),' '); p.(t{1}) = str2double(t{2}); end, fclose(fid);
g = compute_operating_params(p);
fem_cs_selftest(p, g);
