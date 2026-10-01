function g = compute_operating_params(p)
%COMPUTE_OPERATING_PARAMS Derive PFC constants and geometry-independent scaling shared by every candidate.
%
%   g = COMPUTE_OPERATING_PARAMS(p) takes the raw machine/search
%   parameters (as returned by READ_MACHINE_INPUT) and returns a struct g
%   with the physical constants and derived quantities every PF coil and
%   every WP_h/turns/layers/Iop candidate shares.
%
%   Relocated from the PF legacy archive (PF_opt_VNS.m, the block computed
%   once at the top of the script, before the per-coil loop).

g = struct();
g.Mu_0 = 4*pi*1e-7;

g.S_hoop_amm = p.S_hoop_allow/p.SF_hoop; % [Pa]
g.S_T_amm = g.S_hoop_amm*p.Tresca_factor; % [Pa]
end
