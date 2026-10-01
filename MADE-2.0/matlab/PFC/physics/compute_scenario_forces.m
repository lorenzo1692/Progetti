function [Fz, Fr] = compute_scenario_forces(geom, nfr, nfz)
%COMPUTE_SCENARIO_FORCES Net axial and radial Lorentz force on every conductor, for every plasma scenario.
%
%   [Fz, Fr] = COMPUTE_SCENARIO_FORCES(geom, nfr, nfz) returns, in MN, the
%   net axial force Fz and the net radial (outward) force Fr on each row of
%   the combined CS+PFC+plasma table geom (READ_COIL_GEOMETRY), one column
%   per scenario, from the field of all the OTHER conductors acting on the
%   current of each one (a conductor's own field exerts no net axial force
%   on itself). Every conductor block (geom.dr x geom.dz) is represented by
%   nfr x nfz uniform filaments (default 4 x 8) carrying the scenario
%   Ampere-turns (geom.MAt_signed); fields via XBR/XBZ.
%
%   Purpose: an in-pipeline alternative to the external axial-force file
%   FZ_PF.xlsx read by READ_AXIAL_FORCE (the design models only need
%   max(abs(Fz),[],2)*1e6 [N]). Also used as a consistency check of that
%   file: the two should agree within the filament discretization.

if nargin < 2, nfr = 4; end
if nargin < 3, nfz = 8; end
nrow = numel(geom.R);
nscen = size(geom.MAt_signed, 2);
nper = nfr*nfz;

Rf = zeros(nrow*nper, 1); Zf = Rf; owner = Rf;
Ifil = zeros(nrow*nper, nscen);
k = 0;
for r = 1:nrow
    rr = geom.R(r) + ((1:nfr)-0.5)/nfr*geom.dr(r) - geom.dr(r)/2;
    zz = geom.Z(r) + ((1:nfz)-0.5)/nfz*geom.dz(r) - geom.dz(r)/2;
    [RR, ZZ] = meshgrid(rr, zz);
    idx = k + (1:nper);
    Rf(idx) = RR(:); Zf(idx) = ZZ(:); owner(idx) = r;
    Ifil(idx, :) = repmat(geom.MAt_signed(r, :)/nper, nper, 1);
    k = k + nper;
end

Fz = zeros(nrow, nscen); Fr = zeros(nrow, nscen);
for r = 1:nrow
    mine = owner == r;
    others = ~mine;
    Rp = Rf(mine); Zp = Zf(mine);
    Sc = 1e-32*ones(sum(others), 1);
    BR = xbr(Rp, Zp, Rf(others), Zf(others), Sc)*Ifil(others, :);
    BZ = xbz(Rp, Zp, Rf(others), Zf(others), Sc)*Ifil(others, :);
    Im = Ifil(mine, :);
    Fz(r, :) = sum(-BR.*Im.*(2*pi*Rp), 1)*1e-6;
    Fr(r, :) = sum( BZ.*Im.*(2*pi*Rp), 1)*1e-6;
end
end
