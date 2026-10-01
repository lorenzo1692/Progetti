function [BR, BZ] = fem_field_cases(Rp, Zp, Rf, Zf, Ifil)
%FEM_FIELD_CASES Radial/vertical field at points from many circular filaments, several current cases at once.
%
%   [BR, BZ] = FEM_FIELD_CASES(Rp, Zp, Rf, Zf, Ifil) returns the field [T]
%   at the points (Rp, Zp) (column vectors) produced by circular
%   filaments at (Rf, Zf) (column vectors) carrying, in load case k, the
%   currents Ifil(:,k) [A] (nfil x ncase). The influence matrices come from
%   XBR/XBZ (elliptic integrals, shared with TF/CS) and are evaluated in
%   chunks of points, so memory stays bounded for fine meshes.

npts = numel(Rp);
nfil = numel(Rf);
ncase = size(Ifil, 2);
BR = zeros(npts, ncase); BZ = zeros(npts, ncase);
Sc = 1e-32*ones(nfil, 1);
chunk = max(1, floor(2e6/nfil));
for i0 = 1:chunk:npts
    ii = i0:min(npts, i0+chunk-1);
    KBR = xbr(Rp(ii), Zp(ii), Rf, Zf, Sc);
    KBZ = xbz(Rp(ii), Zp(ii), Rf, Zf, Sc);
    BR(ii, :) = KBR*Ifil;
    BZ(ii, :) = KBZ*Ifil;
end
end
