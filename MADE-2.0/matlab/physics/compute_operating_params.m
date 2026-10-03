function g = compute_operating_params(p)
%COMPUTE_OPERATING_PARAMS Derive TF coil operating parameters.
%
%   g = COMPUTE_OPERATING_PARAMS(p) takes the raw machine parameters (as
%   returned by READ_MACHINE_INPUT) and returns a struct g with the
%   derived operating quantities shared by every candidate design: coil
%   geometry envelope, total TF current, peak field, and the
%   plasma-side-limited discharge time estimate.

g = struct();
g.Mu_0 = 4e-7*pi;
g.theta_TF = 2*pi/p.n_TF;

g.RTFi = p.R0 - p.R0/p.A - p.wb;                        % [m] TF inner-leg outer radius
g.RTFo = (p.R0 + p.R0/p.A)*(1/p.ripple)^(1/p.n_TF);      % [m] TF outer-leg inner radius
g.R_VV = g.RTFi*p.R_VV_margin_factor;                    % [m]

% Bending-free shell-model factors: the original script recomputed
% 0.5*log(RTFo/RTFi) four times (inductance, hoop tension x2, DTF sizing)
% with the exact same inputs; computed once here and reused everywhere.
g.k_bf = 0.5*log(g.RTFo/g.RTFi);
g.r_bf = sqrt(g.RTFo*g.RTFi);

amp_corr_R0 = 1/(p.R0/g.RTFi) * ...
    (1 + 1/((p.R0/g.RTFi)^p.n_TF - 1) + 1/((g.RTFo/p.R0)^p.n_TF - 1));

NI_MA = (2*pi*p.R0*p.B0/g.Mu_0)/p.n_TF*1e-6;             % [MA] total TF current

B_PHI_TF = g.Mu_0*p.n_TF*NI_MA*1e6/(2*pi*g.RTFi - p.dr_plasma_side); % Max field on TF
g.B_PHI_0 = B_PHI_TF*amp_corr_R0;
g.NI = NI_MA*1e6;                                        % [A] total TF current

g.R_TF_Outerleg = g.RTFo;                                % Outer-leg inner radius
g.R_TF_Innerleg = g.RTFi;                                % Inner-leg outer radius

g.Tau_discharge1 = p.B0*g.NI*p.n_TF*(p.R0/p.A)^2/(g.R_VV*p.S_VV);

g.B_PHI_TF = B_PHI_TF*p.corr_B_WP;                       % Max field on TF, corrected
end
