function DATA = scan_wp_designs(p, g, env, WP_h, R_center, MAt_target, FZ, L_factor, n_PF_label)
%SCAN_WP_DESIGNS Explore every turns/layers/Iop candidate design for one PFC coil at one WP_h.
%
%   DATA = SCAN_WP_DESIGNS(p, g, env, WP_h, R_center, MAt_target, FZ,
%   L_factor, n_PF_label) does, in two stages, what PF_opt_VNS.m's own
%   nested loops did inline:
%
%   1. Pre-filter every combination of turns-per-grade (env.combT),
%      layers-per-grade (env.combL) and operating current (env.Iop_range)
%      first against the Ampere-turns actually required by the plasma
%      scenario (MAt_target, +MAt_tolerance window - not just a feasibility
%      bound: PFC coils must supply close to a specific target, unlike
%      TF/CS's pure feasibility-window filters), then against a crude
%      self-field estimate window [min_B, max_B] (PHYSICS/
%      ESTIMATE_PEAK_FIELD). Surviving candidates are sorted by total
%      Ampere-turns, as in the legacy driver.
%
%   2. For every surviving candidate, sizes the conductor (SIZE_CONDUCTOR_
%      CICC), the jacket geometry against the ring stress model
%      (SIZE_CICC_CABLE), re-derives the field once the geometry has
%      converged (PHYSICS/EMAG_FIELD_FORCES), and keeps every candidate
%      that passes the geometric consistency and stress-limit checks.
%      FCGR fatigue life is also computed for every accepted candidate
%      (wired in for consistency with the CS pipeline - PF_opt_VNS.m
%      itself never called it, see manuale PFC).
%
%   Returns DATA as a table, one row per feasible design point, for this
%   PF coil and this WP_h alone - MAIN_WP_PFC_DESIGN.m concatenates the
%   tables from every WP_h in the sweep into one coil-level result table,
%   matching what the legacy driver did across its own separate per-WP_h
%   output files.
%
%   This function is a direct, line-by-line port of the scan loop in
%   PF_opt_VNS.m (see manuale PFC for the full derivation). Two behaviors
%   are preserved unchanged and flagged here (see manuale PFC, "decisioni
%   aperte"): FZ is an externally supplied scenario value (READ_AXIAL_
%   FORCE), not the FZmax this same physics chain computes internally from
%   the coil's own geometry; and the check_B outer convergence loop's
%   n_grades>1 conc_spire construction reproduces the legacy switch
%   verbatim, including its case-2/case-3 asymmetry (never exercised by
%   the n_grades=1 VNS baseline).

n_grades = p.n_grades;
Mu_0 = g.Mu_0;

%% 1. Pre-filter turns/layers/Iop combinations against the MAt target and the crude B window
good = 0;
max_possible = size(env.combT,1)*size(env.combL,1)*numel(env.Iop_range);
valid_comb = zeros(max_possible, 2*n_grades + 2); % [layers, turns, Iop, MAt_comb]

for t = 1:size(env.combT,1)
    for l = 1:size(env.combL,1)
        turns = env.combT(t,:);
        layers = env.combL(l,:);

        for i = 1:numel(env.Iop_range)
            Iop = env.Iop_range(i);
            MAt_comb = sum(layers.*turns)*Iop;
            if MAt_comb < MAt_target || MAt_comb > p.MAt_tolerance*MAt_target
                continue
            end

            Bcheck = estimate_peak_field(R_center, WP_h - 2*p.grins_h, sum(turns), sum(layers), Iop, p.shape_cable);
            if Bcheck <= p.min_B || Bcheck >= p.max_B
                continue
            end

            good = good + 1;
            valid_comb(good, :) = [layers, turns, Iop, MAt_comb];
        end
    end
end
valid_comb = sortrows(valid_comb(1:good, :), 2*n_grades+2);

%% 2. Full physics chain per candidate
DATA = table();
counter = 0;

for comb = 1:size(valid_comb,1)
    coil = struct();
    coil.n_layers = valid_comb(comb, 1:n_grades);
    coil.n_turns  = valid_comb(comb, n_grades+1 : 2*n_grades);
    coil.Iop      = valid_comb(comb, 2*n_grades+1);

    coil.n_spire = coil.n_layers .* coil.n_turns;
    coil.N_spire = p.n_moduli*sum(coil.n_spire);

    coil.n_t = sum(coil.n_turns);
    coil.n_l = sum(coil.n_layers);
    coil.n_g = n_grades;

    r_coil = 2;
    coil.dy = (WP_h - 2*p.grins_h)/coil.n_t;
    coil.dx = (p.shape_cable == 200)*coil.dy + (p.shape_cable ~= 200)*(coil.dy/r_coil);

    coil.Ri = R_center - coil.dx*coil.n_l/2;
    coil.Re = R_center + coil.dx*coil.n_l/2;

    [coil.Bsum, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = emag_field_forces( ...
        coil.dy, coil.dx, coil.Ri, coil.n_t, coil.n_l, coil.n_g, coil.Iop);

    coil.k_Bmax = coil.Bsum/(Mu_0*coil.N_spire*coil.Iop/(WP_h - 2*p.grins_h));

    %% Inductance / energy (monoturn coupling scaled by L_factor for the discharge voltage only)
    Rp = (coil.Ri + coil.Re)/2; Zp = 0;
    DRc = coil.Re - coil.Ri - 2*p.grins_w; DZc = WP_h - 2*p.grins_h;
    coil.L = xlm(Rp, Zp, Rp, Zp, DRc, DZc, 1)*coil.N_spire^2;
    coil.E = 0.5*coil.L*coil.Iop^2*1e-6; % [MJ]

    coil.Tau_discharge = p.Tau_discharge_fixed;
    coil.V_ = ceil(coil.L*L_factor*coil.Iop/coil.Tau_discharge);

    %% Field per grade (concatenated spire count, from outermost to innermost)
    % Reproduces the legacy switch verbatim, including its case-2/case-3
    % asymmetry - not exercised by n_grades=1 (see function header).
    switch n_grades
        case 1
            coil.conc_spire = sum(coil.n_spire);
        case 2
            coil.conc_spire = [sum(coil.n_spire), coil.n_spire(2)];
        case 3
            coil.conc_spire = [sum(coil.n_spire), sum(coil.n_spire)-coil.n_spire(1), coil.n_spire(3)];
        otherwise
            error('scan_wp_designs:unsupported_n_grades', 'n_grades must be 1, 2 or 3 (got %g).', n_grades);
    end

    coil.B_grades = zeros(1, n_grades);
    for var = 1:n_grades
        coil.B_grades(var) = p.n_moduli*coil.conc_spire(var)*coil.Iop*Mu_0/(WP_h - 2*p.grins_w)*coil.k_Bmax;
    end

    %% Per-grade conductor + jacket sizing, iterated until the assumed field converges
    coil.t_ins = p.tins*ones(1, n_grades);
    coil.r_SC = 0.002*p.Increm*ones(1, n_grades);
    coil.Cond_h = zeros(1, n_grades); coil.Cond_w = zeros(1, n_grades);
    coil.type_cable = cell(1, n_grades); coil.SC_mat = cell(1, n_grades);
    coil.N_Cu = zeros(1, n_grades); coil.N_sc = zeros(1, n_grades); coil.N_tot = zeros(1, n_grades);
    coil.S_Cable = zeros(1, n_grades); coil.S_REBCO = zeros(1, n_grades); coil.S_Cu_HTS = zeros(1, n_grades);
    coil.THS = zeros(1, n_grades);
    coil.JT = zeros(1, n_grades); coil.SC_w = zeros(1, n_grades); coil.SC_h = zeros(1, n_grades);
    coil.Ri_grades = zeros(1, n_grades); coil.Re_grades = zeros(1, n_grades);
    coil.S_hoop = []; coil.S_T = [];
    invalid = false;

    check_B = false;
    while ~check_B
        for var = 1:n_grades
            B_local = coil.B_grades(var) + p.B_background;
            [type, N_Cu, N_Sc, N_tot, S_Cable, S_REBCO, S_Cu_HTS, THS, mat_code] = ...
                size_conductor_cicc(B_local, coil.Iop, coil.Tau_discharge, p.WP_SC_type);

            coil.type_cable(var) = type;
            coil.N_Cu(var) = N_Cu; coil.N_sc(var) = N_Sc; coil.N_tot(var) = N_tot;
            coil.S_Cable(var) = S_Cable; coil.S_REBCO(var) = S_REBCO; coil.S_Cu_HTS(var) = S_Cu_HTS;
            coil.THS(var) = THS;
            switch mat_code
                case 0, coil.SC_mat{var} = 'Nb3Sn';
                case 1, coil.SC_mat{var} = 'Nb_Ti';
                case 2, coil.SC_mat{var} = 'REBCO';
            end

            coil.Cond_h(var) = (WP_h - 2*p.grins_h)/coil.n_turns(var);

            sizing_in = struct();
            sizing_in.Cond_h = coil.Cond_h(var);
            sizing_in.S_Cable_var = coil.S_Cable(var);
            sizing_in.r_SC = coil.r_SC(var);
            sizing_in.tins = coil.t_ins(var);
            sizing_in.type_cable = coil.type_cable(var);
            sizing_in.Cond_w1 = coil.Cond_h(1); % Cond_w(1) in the legacy driver - equal to this grade's Cond_h when n_grades=1
            sizing_in.n_layers_var = coil.n_layers(var);
            sizing_in.R_center = R_center;
            sizing_in.dy = coil.Cond_h(var);
            sizing_in.n_t = coil.n_t; sizing_in.n_l = coil.n_l; sizing_in.n_g = coil.n_g;
            sizing_in.Iop = coil.Iop;
            sizing_in.FZ = FZ;
            sizing_in.WP_h = WP_h;
            sizing_in.var = var;
            sizing_in.JT_min = p.JT_min;
            sizing_in.JT_step = p.JT_step;
            sizing_in.S_hoop_allow = g.S_hoop_amm;
            sizing_in.Tresca_factor = p.Tresca_factor;
            sizing_in.max_iter = p.max_sizing_iterations;

            sized = size_cicc_cable(sizing_in);

            coil.SC_w(var) = sized.SC_w; coil.SC_h(var) = sized.SC_h;
            coil.JT(var) = sized.JT; coil.Cond_w(var) = sized.Cond_w;
            coil.S_CICC(var) = sized.S_CICC; coil.S_JT(var) = sized.S_JT;
            coil.Ri_grades(var) = sized.Ri_grade; coil.Re_grades(var) = sized.Re_grade;
            coil.Bsum = sized.Bsum; coil.Bmin(var) = sized.Bmin; coil.Bmax(var) = sized.Bmax;
            coil.S_hoop(var) = sized.S_hoop; coil.S_T(var) = sized.S_T;
        end

        if coil.B_grades(var) > coil.Bsum
            coil.B_grades(var) = coil.Bsum;
        else
            check_B = true;
        end
    end

    %% Global winding-pack geometry
    WP_w0 = coil.Re_grades(end) - coil.Ri_grades(end);
    coil.Ri_grades(1) = R_center - WP_w0/2;
    for var = 1:n_grades
        coil.Re_grades(var) = coil.Ri_grades(var) + coil.Cond_w(var)*coil.n_layers(var);
        if var < n_grades
            coil.Ri_grades(var+1) = coil.Re_grades(var) + p.ins_grades;
        end
    end

    coil.Ri = coil.Ri_grades(1) - p.grins_w;
    coil.Re = coil.Re_grades(n_grades) + p.grins_w;

    coil.A_Cable_tot = sum(coil.S_CICC.*coil.n_spire);
    coil.A_JT_tot = sum(coil.S_JT.*coil.n_spire);
    coil.r_cable = coil.Cond_h./coil.Cond_w;

    if coil.Ri < 0 || min(coil.SC_w) < p.SC_w_min || ...
       min(coil.r_cable) < p.r_cable_min || max(coil.r_cable) > p.r_cable_max || any(coil.JT < p.JT_min)
        continue % Skip invalid design
    end

    coil.SH = max(abs(coil.S_hoop));
    coil.ST = max(coil.S_T);

    if ~(coil.ST < g.S_T_amm && coil.SH < g.S_hoop_amm && coil.SH > 0 && coil.ST > 0 && (coil.Re - coil.Ri) < p.WP_radial_build_max)
        continue % Skip design failing the mechanical/geometric acceptance check
    end

    counter = counter + 1;

    Jeng = (coil.Iop*1e-3)./(coil.Cond_w.*coil.Cond_h)*1e-3; % [A/mm^2]
    coil.plasma_cycles = zeros(1, n_grades);
    for var = 1:n_grades
        coil.plasma_cycles(var) = fcgr(coil.JT(var), coil.Cond_w(var), coil.S_hoop(var)*1e-6);
    end

    S_Cable_mm2 = round(1e6*coil.S_Cable);
    S_REBCO_mm2 = round(1e6*coil.S_REBCO);
    S_Cu_HTS_mm2 = round(1e6*coil.S_Cu_HTS);

    row = table(n_PF_label, WP_h, coil.L, coil.Bsum, Jeng, coil.Iop*1e-3, coil.N_spire*coil.Iop, ...
        coil.Ri, coil.Re, Rp, coil.n_layers, coil.n_turns, coil.N_spire, coil.Re-coil.Ri, WP_h, ...
        coil.SC_mat, coil.Cond_w, coil.Cond_h, coil.JT, coil.r_cable, coil.N_Cu, coil.N_sc, ...
        coil.B_grades, S_REBCO_mm2, S_Cu_HTS_mm2, S_Cable_mm2, coil.THS, coil.E, ...
        coil.plasma_cycles, coil.SH*1e-6, coil.ST*1e-6, L_factor, coil.Tau_discharge, coil.V_, ...
        'VariableNames', { ...
            'n_PF', 'WP_h0', 'L', 'Btot', 'Jeng', 'Iop_kA', 'Itot', ...
            'Ri', 'Re', 'Rp', 'n_layers', 'n_turns', 'N_spire', 'WP_w', 'WP_h', ...
            'SC_mat', 'Cond_w', 'Cond_h', 'JT', 'r_cable', 'N_Cu', 'N_SC', ...
            'B_local', 'S_REBCO_mm2', 'S_Cu_HTS_mm2', 'S_Cable_mm2', 'THS', 'E', ...
            'plasma_cycles', 'SH_MPa', 'ST_MPa', 'L_tot_factor', 'Tau_discharge', 'V_max'});

    DATA = [DATA; row]; %#ok<AGROW>
end
end
