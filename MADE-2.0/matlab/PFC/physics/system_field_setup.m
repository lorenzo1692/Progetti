function sys = system_field_setup(geom, row, Fz_scenarios, nfr, nfz)
%SYSTEM_FIELD_SETUP Filament model of every OTHER conductor of the machine, for system-based sizing.
%
%   sys = SYSTEM_FIELD_SETUP(geom, row, Fz_scenarios, nfr, nfz) prepares what
%   SYSTEM_FIELD_EVAL needs to give the background field that the rest of
%   the machine (CS, the other PF coils, plasma) produces on the PF coil in
%   geom row `row`, in every plasma scenario: each other conductor block
%   (geom.dr x geom.dz) becomes nfr x nfz filaments (default 4 x 8) carrying
%   its signed scenario Ampere-turns (geom.MAt_signed), the same model as
%   COMPUTE_SCENARIO_FORCES. Fz_scenarios (optional, [MN], rows as in geom,
%   one column per scenario) is kept so that the vertical stress of each
%   scenario can be evaluated with its own net axial force.
%
%   Used only when the input option system_sizing = 1: the default (0) sizes
%   the coil on its own field, which needs no scenario definition.

if nargin < 3, Fz_scenarios = []; end
if nargin < 4, nfr = 4; end
if nargin < 5, nfz = 8; end
nrow = numel(geom.R);
nper = nfr*nfz;
others = setdiff(1:nrow, row);
Rf = zeros(numel(others)*nper, 1); Zf = Rf;
Ifil = zeros(numel(others)*nper, size(geom.MAt_signed, 2));
k = 0;
for r = others
    rr = geom.R(r) + ((1:nfr)-0.5)/nfr*geom.dr(r) - geom.dr(r)/2;
    zz = geom.Z(r) + ((1:nfz)-0.5)/nfz*geom.dz(r) - geom.dz(r)/2;
    [RR, ZZ] = meshgrid(rr, zz);
    idx = k + (1:nper);
    Rf(idx) = RR(:); Zf(idx) = ZZ(:);
    Ifil(idx, :) = repmat(geom.MAt_signed(r, :)/nper, nper, 1);
    k = k + nper;
end
sys = struct('Rf', Rf, 'Zf', Zf, 'Ifil', Ifil, 'Zrow', geom.Z(row), ...
    'MAt_row', geom.MAt_signed(row, :), 'Fz_row', []);
if ~isempty(Fz_scenarios)
    sys.Fz_row = Fz_scenarios(row, :);
end
end
