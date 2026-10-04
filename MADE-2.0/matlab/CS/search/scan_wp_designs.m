function DATA = scan_wp_designs(p, g, env)
%SCAN_WP_DESIGNS Explore every turns/layers/Iop candidate CS design.
%
%   DATA = SCAN_WP_DESIGNS(p, g, env) pre-filters every combination of
%   turns-per-grade (env.combT), layers-per-grade (env.combL) and
%   operating current (env.Iop_range) against a crude infinite-solenoid
%   field window [p.min_B, p.max_B], then for every surviving candidate
%   sizes the conductor (SIZE_CONDUCTOR_CICC), the jacket geometry
%   (SIZE_CICC_CABLE), the field/forces (EMAG_FIELD_FORCES), the final
%   Tresca stress (EQV_STRESS_COIL_CICC) and the fatigue life (FCGR), and
%   keeps every candidate that passes the geometric consistency and
%   stress-limit checks. Returns them as a table, one row per feasible
%   design point.
%
%   This function is a direct, line-by-line port of the scan loop in
%   CS_opt_VNS.m (see manuale CS for the full derivation), restructured
%   the same way the TF port restructured its own monolithic script:
%   duplicated sizing logic delegated to SIZE_CONDUCTOR_CICC/
%   SIZE_CICC_CABLE, everything else kept as close to the original control
%   flow as possible. Known fidelity note (see manuale CS): k_Bmax is
%   algebraically always 1 given the currently-active (non bmax_check)
%   code path.
%
%   Deliberate deviation from the legacy driver: the axial force in the
%   final stress check is each candidate's own stack compression FZmax
%   from EMAG_FIELD_FORCES with the full stack evaluated (saved as column
%   Fz_MN), not the fixed 35 MN the legacy CS_opt_VNS.m overwrote it with.

% fatigue (shared/fatigue/fcgr.m): JK2LB jacket, ITER CS criteria
% (C0, m: ITER DDD CS p. 6-81; surface flaw 1e-6 m^2 x SFa 2, SFK 1.5, SFN 2)
FCGR_CS = struct('C0', 1.75e-13, 'm', 3.7, 'mw', 0.5, 'residual_stress', 200, 'KIC', 200, ...
    'flaw_aspect', 3, 'flaw_type', 1, 'flaw_area', 1e-6, 'SFa', 2, 'SFK', 1.5, 'SFN', 2);

Mu_0 = g.Mu_0;
n_moduli = p.n_moduli;
n_grades = p.n_grades;
h_stack = g.h_stack;

%% 1. Pre-filter turns/layers/Iop combinations against the crude B window
good = 0;
max_possible = size(env.combT,1) * size(env.combL,1) * numel(env.Iop_range);
comb_nli = zeros(max_possible, 2*n_grades + 1);

for t = 1:size(env.combT,1)
    for l = 1:size(env.combL,1)
        turns = env.combT(t,:);
        layers = env.combL(l,:);
        N_total = sum(turns .* layers) * n_moduli;

        for i = 1:numel(env.Iop_range)
            Iop = env.Iop_range(i);
            B_check = Mu_0 * Iop * N_total / h_stack; % infinite-solenoid estimate
            if B_check < p.min_B || B_check > p.max_B
                continue
            end
            good = good + 1;
            comb_nli(good, :) = [layers, turns, Iop];
        end
    end
end
comb_nli = sortrows(comb_nli(1:good, :), 2*n_grades+1);

%% 2. Full physics chain per candidate
DATA = table();
counter = 0;

for comb = 1:size(comb_nli,1)
    coil = struct();
    coil.n_layers = comb_nli(comb, 1:n_grades);
    coil.n_turns  = comb_nli(comb, n_grades+1 : 2*n_grades);
    coil.Iop      = comb_nli(comb, end);

    coil.N_spire = n_moduli * sum(coil.n_layers .* coil.n_turns);

    %% Flux/inductance (per module)
    k_Bmax = 1; % algebraically always 1 on this code path (see fidelity note above)

    Re_ = g.Re_WP_outer; % = Re - grins_w
    Ri_ = g.Re / g.guess_shape;

    Rp_L = (Ri_+Re_)/2; Zp_L = g.WP_h/2;
    coil.L = xlm(Rp_L, Zp_L, Rp_L, Zp_L, Re_-Ri_, g.WP_h, 1) * (coil.N_spire/n_moduli)^2;

    coil.Tau_discharge = ceil(coil.L * coil.Iop / p.V_MAX * p.SF_V);

    coil.n_spire = coil.n_layers .* coil.n_turns;
    coil.conc_spire = cumsum(coil.n_spire);
    coil.B_grades  = n_moduli * coil.conc_spire * coil.Iop * Mu_0 / h_stack * k_Bmax;
    coil.B_partial = n_moduli * coil.n_spire    * coil.Iop * Mu_0 / h_stack * k_Bmax;

    %% Per-grade conductor sizing
    coil.type_cable = cell(1, n_grades);
    coil.SC_mat = cell(1, n_grades);
    coil.N_Cu = zeros(1, n_grades); coil.N_Sc = zeros(1, n_grades); coil.N_tot = zeros(1, n_grades);
    coil.S_Cable = zeros(1, n_grades); coil.S_REBCO = zeros(1, n_grades); coil.S_Cu_HTS = zeros(1, n_grades);
    coil.THS = zeros(1, n_grades);
    coil.shape_cable = zeros(1, n_grades);

    for var = 1:n_grades
        [coil.type_cable(var), coil.N_Cu(var), coil.N_Sc(var), coil.N_tot(var), ...
         coil.S_Cable(var), coil.S_REBCO(var), coil.S_Cu_HTS(var), coil.THS(var), mat_code] = ...
            size_conductor_cicc(coil.B_grades(var), coil.Iop, coil.Tau_discharge, p.WP_SC_type);

        switch mat_code
            case 0, coil.SC_mat{var} = 'Nb3Sn'; coil.shape_cable(var) = 201;
            case 1, coil.SC_mat{var} = 'Nb_Ti'; coil.shape_cable(var) = 201;
            case 2, coil.SC_mat{var} = 'REBCO'; coil.shape_cable(var) = 201;
        end
    end

    if min(coil.N_Sc) <= 1
        continue % invalid design
    end

    %% Per-grade jacket geometry sizing
    coil.PB = zeros(1, n_grades);
    coil.Cond_h = zeros(1, n_grades); coil.Cond_w = zeros(1, n_grades);
    coil.SC_w = zeros(1, n_grades); coil.SC_h = zeros(1, n_grades);
    coil.JT = zeros(1, n_grades); coil.S_CICC = zeros(1, n_grades); coil.S_JT = zeros(1, n_grades);
    coil.Ri_grades = zeros(1, n_grades+1);
    coil.Phi = zeros(1, n_grades);

    Re_grades = zeros(1, n_grades+1);
    Re_grades(1) = g.Re_WP_outer;

    for var = 1:n_grades
        coil.PB(var) = coil.B_grades(var)^2 / (2*Mu_0);
        coil.Cond_h(var) = (g.WP_h - 2*g.grins_h) / coil.n_turns(var);
        r_Sc_var = g.r_SC;

        sizing_in = struct();
        sizing_in.Cond_h = coil.Cond_h(var);
        sizing_in.S_Cable_var = coil.S_Cable(var);
        sizing_in.r_SC = r_Sc_var;
        sizing_in.tins = g.tins;
        sizing_in.shape_cable = coil.shape_cable(var);
        sizing_in.type_cable = coil.type_cable(var);
        sizing_in.PB = coil.PB(var);
        sizing_in.Re_this_grade = Re_grades(var);
        sizing_in.Re_WP_outer = g.Re_WP_outer;
        sizing_in.n_layers = coil.n_layers(var);
        sizing_in.S_hoop_allow = p.S_hoop_allow;
        sizing_in.SF_hoop = p.SF_hoop;
        sizing_in.JT_min = p.JT_min;
        sizing_in.JT_step = p.JT_step;
        sizing_in.max_iter = p.max_sizing_iterations;

        sizing_out = size_cicc_cable(sizing_in);

        coil.SC_w(var) = sizing_out.SC_w;
        coil.SC_h(var) = sizing_out.SC_h;
        coil.JT(var) = sizing_out.JT;
        coil.Cond_w(var) = sizing_out.Cond_w;
        coil.S_CICC(var) = sizing_out.S_CICC;
        coil.S_JT(var) = sizing_out.S_JT;
        coil.Ri_grades(var) = sizing_out.Ri_grade;
        Re_grades(var+1) = sizing_out.Re_grade_next;

        coil.Phi(var) = (pi/3) * (Re_grades(var)^2 + coil.Ri_grades(var)^2 + coil.Ri_grades(var)*Re_grades(var)) * coil.B_partial(var);
    end

    %% Winding pack geometry and consistency checks
    coil.Ri = coil.Ri_grades(1) - g.grins_w;
    coil.WP_w = g.Re - coil.Ri;
    coil.Rm = (g.Re + coil.Ri)/2;
    coil.total_length = coil.n_spire * coil.Rm * 2*pi;

    coil.A_Cable_tot = sum(coil.S_CICC .* coil.n_spire);
    coil.A_JT_tot = sum(coil.S_JT .* coil.n_spire);
    coil.r_cable = coil.Cond_h ./ coil.Cond_w;

    if min(coil.SC_w) <= p.SC_w_min || min(coil.r_cable) < p.r_cable_min || ...
       min(coil.JT) < p.JT_min_check || coil.Ri < 0 || coil.Ri < p.Ri_min || ...
       max(coil.r_cable) > p.r_cable_max
        continue
    end

    coil.Phi_TOT = sum(coil.Phi);
    coil.L = n_moduli * coil.L;
    coil.E = 0.5 * coil.L * coil.Iop^2 * 1e-6; % [MJ]

    %% Field/forces (innermost grade) and final stress check
    var = n_grades;
    [Bsum, ~, ~, ~, FZmax, ~, ~, ~, ~, ~, ~] = emag_field_forces( ...
        coil.Cond_h, coil.Cond_w, coil.Ri_grades, coil.n_turns, coil.n_layers, ...
        n_grades, n_moduli, coil.Iop, g.spacer, 1);

    if Bsum < p.min_B || Bsum > p.max_B
        continue
    end

    coil.PB(var) = Bsum^2 / (2*Mu_0);

    % Fz per candidate, not hardcoded: FZmax from EMAG_FIELD_FORCES with
    % the whole stack evaluated (full_stack=1) is the axial compression
    % through the stack mid-plane (sum of the Fz of the lower half of the
    % stack), in MN - see EQV_STRESS_COIL_CICC's unit note. With
    % full_stack=0 only one module is evaluated and FZmax would be the
    % force on half a module (about 26x too small for the VNS baseline, as
    % the FEM stack verification showed). The legacy driver overwrote this
    % value with a fixed 35 MN; per user decision (01/10/2026) each
    % candidate now uses its own computed compression.
    coil.F_z_MN = FZmax;
    [S_hoop, ~, S_ver, S_T] = eqv_stress_coil_cicc( ...
        FZmax, coil.Ri, g.Re, coil.Cond_h(var), coil.Cond_w(var), coil.JT(var), ...
        coil.SC_h(var), coil.SC_w(var), g.tins, coil.type_cable(var), coil.S_CICC(var), coil.S_JT(var), coil.PB(var)); %#ok<ASGLU>

    coil.S_hoop_max = max(S_hoop);
    coil.S_T_max = max(S_T);
    coil.eps_Sc = coil.SC_w(var) / (2*coil.Ri) * 100; % [%]

    coil.plasma_cycles = zeros(1, n_grades);
    coil.plasma_cycles(var) = fcgr(coil.JT(var), coil.Cond_w(var), coil.S_hoop_max, FCGR_CS);

    %% Acceptance
    if coil.S_T_max < p.S_T_allow/p.SF_T && coil.S_hoop_max < p.S_hoop_allow/p.SF_hoop && ...
       coil.S_hoop_max > 0 && coil.S_T_max > 0

        counter = counter + 1;

        coil.B_max = Bsum;
        coil.B_dim = max(coil.B_grades);
        coil.Iop_kA = coil.Iop * 1e-3;
        coil.MAt = coil.Iop_kA * sum(coil.n_layers .* coil.n_turns) * 1e-3;
        coil.Jeng = coil.Iop / (coil.Cond_h(1) * coil.Cond_w(1)) * 1e-6; % [A/mm^2], outermost grade

        row = table( ...
            coil.plasma_cycles, coil.S_hoop_max, coil.S_T_max, ...
            coil.Phi, coil.Phi_TOT, coil.L, coil.B_max, ...
            coil.Iop_kA, coil.MAt, coil.Ri, g.Re, ...
            coil.n_layers, coil.n_turns, coil.WP_w, g.WP_h, g.spacer, ...
            coil.type_cable, coil.SC_mat, ...
            coil.Cond_w, coil.Cond_h, coil.JT, coil.r_cable, ...
            coil.N_Cu, coil.N_Sc, coil.B_dim, coil.S_REBCO, coil.S_Cu_HTS, ...
            coil.S_Cable, coil.THS, coil.E, ...
            coil.Jeng, coil.eps_Sc, coil.total_length, coil.F_z_MN, ...
            'VariableNames', { ...
                'plasma_cycles', 'S_hoop_max', 'S_T_max', ...
                'Phi', 'Phi_TOT', 'L', 'B_grades', ...
                'Iop_kA', 'Itot_MAt', 'Ri', 'Re', ...
                'n_layers', 'n_turns', 'WP_w', 'WP_h', 'spacer', ...
                'type_cable', 'SC_mat', ...
                'Cond_w', 'Cond_h', 'JT', 'r_cable', ...
                'N_Cu', 'N_Sc', 'B_dim', 'S_REBCO', 'S_Cu_HTS', ...
                'S_Cable', 'THS', 'E', ...
                'Jeng', 'eps_Sc', 'length_module', 'Fz_MN' ...
            });

        DATA = [DATA; row]; %#ok<AGROW>
    end
end
end
