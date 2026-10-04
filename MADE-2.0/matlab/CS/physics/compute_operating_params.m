function g = compute_operating_params(p)
%COMPUTE_OPERATING_PARAMS Derive CS module geometry and scaling shared by every candidate.
%
%   g = COMPUTE_OPERATING_PARAMS(p) takes the raw machine parameters (as
%   returned by READ_MACHINE_INPUT) and returns a struct g with the
%   quantities every candidate design shares: physical constants, the
%   scaled module geometry (WP_h, Re, insulation/spacer thicknesses),
%   the axial stack height, and the guessed inner/outer radius ratio used
%   before the winding pack is actually sized.
%
%   Relocated from the CS legacy archive (CS_opt_VNS.m, the block computed
%   once per WP_h0/Re_0 pair before the turns/layers/Iop scan). The
%   original swept WP_h0 and Re_0 in an outer for-loop that, in
%   CS_opt_VNS.m, only ever iterated one value each - here WP_h0 and Re_0
%   are plain scalar inputs, matching how the TF pipeline treats machine
%   geometry as fixed per run (see manuale CS).

g = struct();
g.Mu_0 = 4*pi*1e-7;

g.WP_h = p.WP_h0 * p.Increm;
g.Re = p.Re_0 * p.Increm;
g.grins_w = p.grins_w * p.Increm;
g.grins_h = p.grins_h * p.Increm;
g.tins = p.tins * p.Increm;
g.spacer = p.spacer * p.Increm;
g.r_SC = p.r_SC * p.Increm;

g.h_stack = g.WP_h*p.n_moduli + (p.n_moduli-1)*g.spacer;

if p.shape_cable == 200
    g.guess_shape = 1;    % RIS
else
    g.guess_shape = 1.5;  % rectangular CICC
end

g.Re_WP_outer = g.Re - g.grins_w; % fixed outer radius reference for the stress check (Re_grades(1) in the legacy driver)
end
