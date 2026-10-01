% Run the FEM verification on a saved design point (Octave dev helper).
here = fileparts(mfilename('fullpath')); cs = fileparts(here);
addpath(genpath(cs)); addpath(fileparts(cs)); addpath(fullfile(here, 'octave_shims'));
S = getenv('CS_ROWS'); load(S);                 % D, p, g
d = csvread(fullfile(here, 'baseline_VNS_07_2026_V4_LG.csv'));
geom.R = d(:,1); geom.Z = d(:,2); geom.dr = d(:,3); geom.dz = d(:,4); geom.MAt_signed = d(:,5:end);
tic; res = fem_cs_verify(D(1), p, g, geom, struct('scenarios', 1:2)); toc
