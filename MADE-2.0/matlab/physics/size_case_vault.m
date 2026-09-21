function out = size_case_vault(in)
%SIZE_CASE_VAULT Size the TF case nose thickness (DTF) against Tresca limits.
%
%   out = SIZE_CASE_VAULT(in) grows the case nose thickness DTF until both
%   the case vault and the innermost jacket layer satisfy their Tresca
%   stress allowables, replacing the corresponding while-loop of the
%   original monolithic script.
%
%   Carries original fix #1: S_c_VT and S_rm_JT are explicitly initialized
%   before the loop (instead of implicitly reusing whatever value a
%   previous, unrelated turns/layers combination happened to leave behind)
%   and fix #4: the loop is capped at in.max_iter passes.
%
%   in fields (all scalars unless noted):
%     Rj_, Ri_, CASE_w, theta_TF, A_WP, A_SC_tot, A_JT_tot, WP_h,
%     lateral_w, E_case, E_cbl, E_jckt, p_rs, S_rm, n_TF, RTFo, RTFi,
%     n_spire1, Iop, Mu_0, dr_plasma_side, S_amm_VT, S_amm_JT,
%     safety_membrane, DTF0, DTF_step, max_iter
%     Ke_cavo_rad, n_turns : vectors over 1:n_layers
%
%   out fields:
%     Rk_, DTF, S_T_VT, S_T_JT, A_CASE, iterations

DTF = in.DTF0;
S_T_JT = Inf;
S_T_VT = Inf;
S_c_VT = Inf;
S_rm_JT = Inf;
it = 0;

while S_T_VT > in.S_amm_VT || S_c_VT > in.S_amm_VT/in.safety_membrane || ...
        S_T_JT > in.S_amm_JT || S_rm_JT > in.S_amm_JT/in.safety_membrane
    it = it + 1;
    if it > in.max_iter
        error('size_case_vault:not_converged', ...
            'DTF sizing did not converge after %d iterations: check the input parameters.', it);
    end

    DTF = DTF + in.DTF_step;                        % If Tresca is exceeded, increase the TF nose thickness
    Rk_ = in.Rj_ - DTF;                              % Innermost Case radius
    CASE_w_l = 2*Rk_*tan(in.theta_TF/2);             % Case low part width
    A_tot = (in.CASE_w + CASE_w_l)*(in.Ri_ - Rk_)/2;
    A_CASE = A_tot - in.A_WP;

    Ke_WP_rad = sum(1./(in.Ke_cavo_rad .* in.n_turns))^-1;
    K_ps_rad = in.E_case/in.dr_plasma_side*(2*in.Ri_*tan(in.theta_TF/2));
    K_lat_rad = in.E_case*(in.lateral_w/2)/in.WP_h;
    K_vault_rad = in.E_case*CASE_w_l/DTF;
    Ke_case_rad = (1/(Ke_WP_rad + 2*K_lat_rad) + 1/K_vault_rad + 1/K_ps_rad)^-1;
    dcr_WP_rad = Ke_WP_rad/Ke_case_rad;
    S_rm_JT = in.S_rm*dcr_WP_rad;                    % Radial membrane stress, innermost jacket layer correction

    h_unit = 1;
    k_steel_tor = in.E_case*(h_unit*(in.Ri_-Rk_)*2*pi*(in.Ri_+Rk_)/2);   % full steel casing
    k_SC_tor = in.E_cbl*(h_unit*(in.Ri_-in.Rj_)*2*pi*(in.Ri_+in.Rj_)/2)*(1-(in.A_WP-in.A_SC_tot)/in.A_WP);
    k_JT_tor = in.E_jckt*(h_unit*(in.Ri_-in.Rj_)*2*pi*(in.Ri_+in.Rj_)/2)*(1-(in.A_WP-in.A_JT_tot)/in.A_WP);
    k_vault_tor = in.E_case*(h_unit*(in.Rj_-Rk_)*2*pi*(in.Rj_+Rk_)/2);
    dcr_vault_tor = k_steel_tor/(k_SC_tor + k_JT_tor + k_vault_tor);

    beta = Rk_/in.Ri_;
    S_c_VT = 2/(1-beta^2)*in.p_rs*dcr_vault_tor;

    k_bf = 0.5*log(in.RTFo/in.RTFi);                                     % k bending free
    T_bf = 0.5*(k_bf*in.n_TF*(in.n_spire1*in.Iop)^2*in.Mu_0/(2*pi));      % Hoop tension along TF longitudinal axis
    S_z = T_bf/(in.A_JT_tot + A_CASE);

    S_T_VT = S_z + S_c_VT;   % Vault Tresca stress
    S_T_JT = S_z + S_rm_JT;  % Jacket Tresca stress
end

out.Rk_ = Rk_;
out.DTF = DTF;
out.S_T_VT = S_T_VT;
out.S_T_JT = S_T_JT;
out.A_CASE = A_CASE;
out.iterations = it;
end
