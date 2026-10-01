function [Bcheck] = estimate_peak_field(Rm, h, n_t, n_l, Iop, shape_cable)
%ESTIMATE_PEAK_FIELD Crude self-field estimate at the coil centre, for pre-filtering candidates.
%
%   Bcheck = ESTIMATE_PEAK_FIELD(Rm, h, n_t, n_l, Iop, shape_cable) lays
%   out an n_t (axial) x n_l (radial) grid of turns centred on (Rm, 0),
%   each carrying Iop, and returns the field magnitude at the winding's
%   own inner-most centre point (Rp=Ri, Zp=0) via XBR/XBZ. Used by
%   SEARCH/GENERATE_COMBINATIONS to discard turns/layers/Iop combinations
%   whose field falls outside [min_B, max_B] before running the full
%   per-candidate physics chain in SCAN_WP_DESIGNS - the same role as the
%   infinite-solenoid B_check pre-filter in the CS pipeline, but computed
%   with the exact filament field model instead of a smeared estimate.
%
%   Relocated and renamed from the PF legacy archive (bmax_check.m, part
%   of MADE_PF.7z), same algorithm. Actively called by PF_opt_VNS.m
%   (unlike the TF/CS ports, where the equivalent function was found
%   unused and left as legacy/orphan - see manuale PFC).

dy = h/n_t;
dx = (shape_cable == 200)*dy + (shape_cable ~= 200)*(dy/2);
Ri = Rm - dx*n_l/2;

n = 0;
data = zeros(n_t*n_l, 2);
for i = 1:n_t
    for j = 1:n_l
        n = n + 1;
        Zc = dy/2 + dy*(i-1);
        Rc = Ri + dx/2 + dx*(j-1);
        data(n, :) = [Zc, Rc];
    end
end

Zp = 0;
Rp = Ri;

Zc = data(:, 1);
Rc = data(:, 2);

% Move any source turn that coincides with the observation point, to
% avoid a singular field evaluation.
mask = (Rc == Rp & Zc == Zp);
Rc(mask) = Rc(mask) - 0.01;

Sc = 1e-30;

Br = xbr(Rp, Zp, Rc, Zc, Sc);
Bz = xbz(Rp, Zp, Rc, Zc, Sc);

Br_ = sum(Br*Iop);
Bz_ = sum(Bz*Iop);

Bcheck = sqrt(Bz_^2 + Br_^2);
end
