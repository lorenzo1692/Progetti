function [DATA, cal] = scan_wp_designs(p, g, env, combT)
%SCAN_WP_DESIGNS Explore every case-wedge x turns/layers candidate design.
%
%   DATA = SCAN_WP_DESIGNS(p, g, env, combT) plays out each candidate
%   winding-pack layout (case wedge thickness from env.lateral_w_min:
%   p.lateral_w_step:env.lateral_w_max, crossed with every turns/layers
%   combination in combT from GENERATE_COMBINATIONS), sizes its CICC/
%   jacket cross-sections (SIZE_CICC_CABLE) and its case nose
%   (SIZE_CASE_VAULT), and keeps every candidate that satisfies the
%   geometric and Tresca stress checks. Returns them as a table, one row
%   per feasible design point, with the same columns as the original
%   monolithic script (minus the vestigial "dp" design-point index, since
%   this solver now handles one machine input per run).
%
%   Prints a self-overwriting progress line to the console (total
%   candidates, analyzed, passed, remaining, elapsed/ETA), updated ~500
%   times over the whole scan so it stays cheap even for very large scans.
%
%   This function is a direct, line-by-line port of the scan loop from
%   legacy/scripts/WP_TF_VNS_Design_Point_2026.m: it calls CICC(...) (conductor/)
%   function exactly as before, and keeps the same array layout and
%   control flow, only delegating the three duplicated jacket-sizing
%   while-loops to SIZE_CICC_CABLE and the case/vault sizing while-loop to
%   SIZE_CASE_VAULT.

maxdim = p.maxdim;
theta_TF = g.theta_TF;
Mu_0 = g.Mu_0;

tins_const = p.turn_insulation_nominal*p.Increm;   % Turn insulation
if ~isfield(p, 'shape_cable') || ~isscalar(p.shape_cable) || ~any(p.shape_cable == [200 201])
    error('scan_wp_designs:shape_cable', 'p.shape_cable must be 200 (RIS) or 201 (Rect).');
end

% field model used to size the grades (see the sizing / field loop below):
%   0 'smeared'    Ampere x corr_B_WP, linear across the WP (old model)
%   1 'calibrated' k(W, Iop) and per-layer profile calibrated at the start
%                  with the discrete model (WP_FIELD_CALIBRATION); with
%                  p.field_verify = 1 (default) the accepted candidates are
%                  then checked with the discrete model and re-sized if the
%                  calibrated field is off by more than field_tol
%   2 'discrete'   smeared first guess, then discrete model and re-sizing
%                  until converged
field_model = 1; if isfield(p, 'field_model') && ~isempty(p.field_model), field_model = p.field_model; end
if ischar(field_model) || isstring(field_model)
    field_model = find(strcmpi(field_model, {'smeared', 'calibrated', 'discrete'})) - 1;
end
if isempty(field_model) || ~isscalar(field_model) || ~any(field_model == [0 1 2])
    error('scan_wp_designs:field_model', 'p.field_model must be 0 (smeared), 1 (calibrated) or 2 (discrete).');
end
% jacket stress model: 0 = analytic formula with the SCF table and
% SCF_transition_provisional (default); 1 = fast surrogate calibrated on
% the 2D FE (JACKET_STRESS_SURROGATE): the jacket thickness of every grade
% is SIZED with it (SIZE_JACKET_SURROGATE: smallest JT with
% jacket_margin x Pm <= Sm and jacket_margin x (Pm+Pb) <= 1.5 Sm), iterated
% with the case nose (docs/DIMENSIONAMENTO_JACKET.md)
scf_model = 0; if isfield(p, 'scf_model') && ~isempty(p.scf_model), scf_model = p.scf_model; end
if ~any(scf_model == [0 1])
    error('scan_wp_designs:scf_model', 'p.scf_model must be 0 (SCF table) or 1 (surrogate).');
end
if scf_model == 1 && p.shape_cable == 200
    error('scan_wp_designs:scf_model_ris', ['scf_model = 1: the jacket stress surrogate is calibrated on ' ...
        'rectangular cables only; use scf_model = 0 for RIS (shape_cable = 200).']);
end
Sm_jacket = p.S_amm_JT; if isfield(p, 'Sm_jacket') && ~isempty(p.Sm_jacket), Sm_jacket = p.Sm_jacket; end
jacket_margin = 1.0;   % default set from the recalibration (docs/DIMENSIONAMENTO_JACKET.md), in progress if isfield(p, 'jacket_margin') && ~isempty(p.jacket_margin), jacket_margin = p.jacket_margin; end
JT_max = 0.010; if isfield(p, 'JT_max') && ~isempty(p.JT_max), JT_max = p.JT_max; end
jacket_max_passes = 4; if isfield(p, 'jacket_max_passes') && ~isempty(p.jacket_max_passes), jacket_max_passes = p.jacket_max_passes; end
n_jacket_rejected = 0;
% maximum cell aspect ratio Cond_w/Cond_h (input files without it: 2)
max_cable_aspect_ratio = 2;
if isfield(p, 'max_cable_aspect_ratio') && ~isempty(p.max_cable_aspect_ratio), max_cable_aspect_ratio = p.max_cable_aspect_ratio; end
field_verify = 1; if isfield(p, 'field_verify') && ~isempty(p.field_verify), field_verify = p.field_verify; end
use_discrete_field = field_model == 2 || (field_model == 1 && field_verify);
cal = [];                 % returned: the start-of-scan field calibration (field_model = 1)
if field_model == 1
    cal = wp_field_calibration(p, g, env);
end
field_tol = 0.05;     if isfield(p, 'field_tol') && ~isempty(p.field_tol), field_tol = p.field_tol; end           % [T]
field_max_iter = 6;   if isfield(p, 'field_max_iter') && ~isempty(p.field_max_iter), field_max_iter = p.field_max_iter; end
n_field_rejected = 0;

counter = 0;
% DATA is intentionally left undefined here: like the original script, it
% is created on the first successful candidate via the indexed assignment
% below (DATA(1,:) = row), which lets MATLAB infer its column structure
% from that first row. If no candidate ever succeeds, it is set to an
% empty table just before returning (see bottom of this function).

% Progress reporting: total candidates = every (lateral_w, turns/layers
% combination) pair, regardless of how far each one gets before a
% "continue" - so the count and the ETA are stable from the first candidate.
n_lateral_w = numel(env.lateral_w_min:p.lateral_w_step:env.lateral_w_max);
n_combT_tot = 0;
for i_ = 1:numel(combT)
    n_combT_tot = n_combT_tot + size(combT{i_}, 1);
end
total_candidates = n_lateral_w * n_combT_tot;
examined = 0;
progress_every = max(1, round(total_candidates/500)); % ~500 console updates over the whole scan
progress_msg_len = 0;
progress_t0 = tic;
fprintf('Scanning %d candidate(s)...\n', total_candidates);

for lateral_w = env.lateral_w_min:p.lateral_w_step:env.lateral_w_max
    for i = 1:numel(combT)
        for j = 1:size(combT{i}, 1)

            examined = examined + 1;
            if mod(examined, progress_every) == 0 || examined == total_candidates
                elapsed = toc(progress_t0);
                remaining = total_candidates - examined;
                rate = examined / max(elapsed, eps);
                eta_s = remaining / max(rate, eps);
                msg = sprintf('  %d/%d analyzed (%.1f%%) | %d passed | %d remaining | elapsed %s | ETA %s', ...
                    examined, total_candidates, 100*examined/total_candidates, counter, remaining, ...
                    format_hms(elapsed), format_hms(eta_s));
                fprintf('%s%s', repmat(char(8), 1, progress_msg_len), msg);
                progress_msg_len = length(msg);
                if examined == total_candidates
                    fprintf('\n');
                end
            end

            n_layers = env.layers_comb(i);
            n_layers0 = n_layers;     % turns/layers count before any grade is added back

            n_turns = zeros(1, maxdim);
            n_turns(1:n_layers0) = combT{i}(j, :);

            % Turns in each grade
            n_spire_ = zeros(1, n_layers);
            n_spire_(1) = sum(n_turns(1:n_layers0));
            for var = 2:n_layers
                n_spire_(var) = n_spire_(var-1) - n_turns(var-1);
            end

            Iop = ceil(g.NI/n_spire_(1));
            if Iop < p.Iop_min || Iop > p.Iop_max
                continue
            end

            % Sizing / field loop: size the candidate with the field of
            % each grade. First pass: the calibrated field k(W, Iop) x
            % Ampere x profile (field_model = 1, default) or the smeared
            % Ampere x corr_B_WP linear field (0 and 2). Then, with
            % field_model = 2 or with 1 and field_verify, compute the real
            % per-layer peak on the sized WP with the discrete 2D model
            % (WP_PEAK_FIELD_FAST, validated vs ANSYS) and re-size with it
            % until the field each grade was sized for matches its peak
            % within p.field_tol. The smeared model's error on the peak
            % grows with the toroidal narrowness of the WP (+1 to +2 T on
            % the plasma side, more on the low-field grades), so no
            % constant corr_B_WP can fix it. Candidates failing any check
            % are rejected on the pass where they fail (a higher field only
            % makes cables bigger).
            n_layers_start = n_layers; n_turns_start = n_turns;
            n_spire_start = n_spire_; Iop_start = Iop;
            B_field_layer = []; reject = false; field_ok = false;
            JT_req = zeros(1, maxdim);   % jacket thickness required by the surrogate (scf_model = 1), per layer
            jacket_passes = 0;
            for field_it = 1:field_max_iter + (scf_model == 1)*jacket_max_passes
                n_layers = n_layers_start; n_turns = n_turns_start;
                n_spire_ = n_spire_start; Iop = Iop_start; WP_w0 = [];
                % Preallocations (maxdim-sized arrays are large enough for any
                % candidate; n_layers-sized arrays auto-grow, exactly as in
                % the original script, if extra grades get appended below)
                type_cable = repmat({'---'}, 1, maxdim);
                Cond_h = zeros(1, maxdim);   Cond_w = zeros(1, maxdim);
                S_Cable = zeros(1, maxdim);  JT = zeros(1, maxdim);
                r_cable = zeros(1, maxdim);  tins = zeros(1, maxdim);
                N_Sc = zeros(1, maxdim);     N_Cu = zeros(1, maxdim);
                THS = zeros(1, maxdim);      S_Cu_HTS = zeros(1, maxdim);
                S_REBCO = zeros(1, maxdim);  Ke_cavo_rad = zeros(1, maxdim);
                B_grade = zeros(1, maxdim);  % field each grade's cable was sized at
                Ke_cavo_tor = zeros(1, maxdim);

                N_tot = zeros(1, n_layers);
                SC_w = zeros(1, n_layers);   SC_h = zeros(1, n_layers);
                R_J = zeros(1, n_layers);
                S_CICC = zeros(1, n_layers); S_JT = zeros(1, n_layers);
                Ri = zeros(1, n_layers);     Re = zeros(1, n_layers);

                % Definition of the magnetic field peak in each grade (linear
                % behavior from B(Re)=Bmax to B(Ri)=0)
                Re(1) = g.R_TF_Innerleg - p.dr_plasma_side - p.GoundIns; % WP inner-leg outer radius (non-insulated)
                B_TF = g.B_PHI_TF;
                B_layers = B_TF .* (n_spire_ ./ n_spire_(1));

                % TF inductance, shell model
                Ntot_turns = p.n_TF*n_spire_(1);
                L = Mu_0*(Ntot_turns*g.k_bf)^2*g.r_bf/2 * ...
                    (besseli(0, g.k_bf) + 2*besseli(1, g.k_bf) + besseli(2, g.k_bf)) / p.n_TF;

                Tau_discharge2 = L*Iop/p.V_MAX; % [s] - discharge in groups of n coils
                Tau_discharge = max([g.Tau_discharge1, Tau_discharge2, 4]);

                % Field each grade is sized for: first pass = smeared model;
                % then the per-layer peak of the discrete model on the
                % previous pass's sized WP (layers added by the turn
                % reduction take the field of the last layer)
                B_cal_layer = nan(1, n_layers);
                if field_model == 1
                    W_cand = 2*Re(1)*tan(theta_TF/2) - lateral_w*2 - p.GoundIns*2;
                    B_cal_layer = wp_field_calibrated(cal, W_cand, Iop, n_spire_(1:n_layers));
                end
                if isempty(B_field_layer) && field_model == 1
                    B_size_layer = B_cal_layer;
                elseif isempty(B_field_layer)
                    B_size_layer = B_layers;
                else
                    B_size_layer = B_field_layer([1:min(n_layers, numel(B_field_layer)), ...
                        numel(B_field_layer)*ones(1, n_layers - numel(B_field_layer))]);
                end
                jump_grade = pick_jump_grades(p.n_grades, B_size_layer, p.grade_B_target_2, p.grade_B_target_3);
                grade_end = max(jump_grade, [jump_grade(2:end)-1, n_layers]);   % repeated jump indices -> one layer

                for ig = 1:numel(jump_grade)
                    var = jump_grade(ig);
                    B = max(B_size_layer(var:grade_end(ig)));
                    [type_cable(var), N_Cu(var), N_Sc(var), N_tot(var), S_Cable(var), ...
                        S_REBCO(var), S_Cu_HTS(var), THS(var)] = cicc(B, Iop, Tau_discharge, p.WP_SC_type, ...
                        p.THS_max_LTS, p.THS_max_HTS);
                    B_grade(var:n_layers) = B;

                    type_cable(var:n_layers) = type_cable(var);
                    N_Cu(var:n_layers) = N_Cu(var);
                    N_Sc(var:n_layers) = N_Sc(var);
                    N_tot(var:n_layers) = N_tot(var);
                    S_Cable(var:n_layers) = S_Cable(var);
                    S_REBCO(var:n_layers) = S_REBCO(var);
                    S_Cu_HTS(var:n_layers) = S_Cu_HTS(var);
                    THS(var:n_layers) = THS(var);
                end

                % Define cable cross-sections in each layer
                T_bf = 0.5*(g.k_bf*p.n_TF*(n_spire_(1)*Iop)^2*Mu_0/(2*pi)); % Hoop tension along TF longitudinal axis

                WP_w0(1) = 2*Re(1)*tan(theta_TF/2) - lateral_w*2 - p.GoundIns*2; % maximum toroidal WP envelope
                S_z_JT = T_bf/(WP_w0(1)^2)/2;
                n_turns_add = 0;
                % Magnetic pressure of the thin-WP formulas (jacket sizing,
                % jacket radial stress, vault): from the field grade 1 is
                % sized for, i.e. the plasma-side J x B field of the chosen
                % field model - the calibrated or discrete peak when
                % available (field_model 1/2), Ampere x corr_B_WP otherwise
                % (field_model 0: B_grade(1) = B_TF, as before)
                p_rs = B_grade(1)^2/(2*Mu_0);

                for var = 1:n_layers
                    Cond_w(var) = WP_w0(1)/n_turns(1);
                    E_cbl = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);

                    sized = size_grade_cable(Cond_w(var), S_Cable(var), p.r_SC_min, p.r_SC_max, tins_const, ...
                        p.E_jckt, E_cbl, p.E_ins, p.shape_cable, p_rs, S_z_JT, ...
                        p.S_amm_JT, p.safety_membrane, max(p.min_JT, JT_req(var) - p.JT_step), p.JT_step, p.max_sizing_iterations);
                    tins(var) = tins_const;
                    Cond_h(var) = sized.Cond_h; JT(var) = sized.JT;
                    SC_w(var) = sized.SC_w;     SC_h(var) = sized.SC_h; R_J(var) = sized.R_J;
                    Ke_cavo_rad(var) = sized.Ke_rad; Ke_cavo_tor(var) = sized.Ke_tor;
                    S_CICC(var) = sized.S_CICC; S_JT(var) = sized.S_JT;

                    Ri(var) = Re(var) - Cond_h(var);          % Inner radius of the i-th grade of the WP
                    Re(var+1) = Ri(var) - p.INS_grades;       % Outer radius of the (i+1)-th grade of the WP
                    WP_w0(var) = Cond_w(var)*n_turns(var);
                    check_w = 2*Ri(var)*tan(theta_TF/2);

                    shrink_iter = 0;
                    while (check_w - (WP_w0(var) + p.GoundIns*2))/2 <= p.toroidal_gap
                        shrink_iter = shrink_iter + 1;
                        if shrink_iter > p.max_sizing_iterations
                            error('scan_wp_designs:turn_shrink_not_converged', ...
                                'Turn-count reduction for grade %d did not converge: check the input parameters.', var);
                        end
                        n_turns(var) = n_turns(var) - 2;
                        Cond_w(var) = WP_w0(1)/n_turns(1);
                        E_cbl = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);

                        sized = size_grade_cable(Cond_w(var), S_Cable(var), p.r_SC_min, p.r_SC_max, tins_const, ...
                            p.E_jckt, E_cbl, p.E_ins, p.shape_cable, p_rs, S_z_JT, ...
                            p.S_amm_JT, p.safety_membrane, max(p.min_JT, JT_req(var) - p.JT_step), p.JT_step, p.max_sizing_iterations);
                        Cond_h(var) = sized.Cond_h; JT(var) = sized.JT;
                        SC_w(var) = sized.SC_w;     SC_h(var) = sized.SC_h; R_J(var) = sized.R_J;
                        Ke_cavo_rad(var) = sized.Ke_rad; Ke_cavo_tor(var) = sized.Ke_tor;
                        S_CICC(var) = sized.S_CICC; S_JT(var) = sized.S_JT;

                        Ri(var) = Re(var) - Cond_h(var);
                        Re(var+1) = Ri(var) - p.INS_grades;
                        WP_w0(var) = Cond_w(var)*n_turns(var);
                        check_w = 2*Ri(var)*tan(theta_TF/2);
                        n_turns_add = n_spire_(1) - sum(n_turns(1:n_layers0));
                    end
                end

                if n_turns(n_layers0) <= 0
                    reject = true; break
                end

                if var == n_layers && n_turns_add >= 1
                    n_layers_add = ceil(n_turns_add/n_turns(n_layers0));
                    n_layers = n_layers + n_layers_add;

                    same_parity = (mod(n_turns(n_layers0), 2) == 0) == (mod(n_turns_add/n_layers_add, 2) == 0);
                    if same_parity
                        pluss = round(n_turns_add/n_layers_add);
                    else
                        pluss = round(n_turns_add/n_layers_add) + 1;
                    end

                    n_turns(n_layers0+1 : n_layers0+n_layers_add) = pluss;

                    var_ = var;
                    for var = var_+1:n_layers
                        Cond_w(var) = WP_w0(1)/n_turns(1);
                        S_Cable(var) = S_Cable(var-1);
                        type_cable(var) = type_cable(var-1);
                        THS(var) = THS(var-1); B_grade(var) = B_grade(var-1);
                        N_Sc(var) = N_Sc(var-1); N_Cu(var) = N_Cu(var-1);
                        S_REBCO(var) = S_REBCO(var-1); S_Cu_HTS(var) = S_Cu_HTS(var-1);
                        E_cbl = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);

                        sized = size_grade_cable(Cond_w(var), S_Cable(var), p.r_SC_min, p.r_SC_max, tins_const, ...
                            p.E_jckt, E_cbl, p.E_ins, p.shape_cable, p_rs, S_z_JT, ...
                            p.S_amm_JT, p.safety_membrane, max(p.min_JT, JT_req(var) - p.JT_step), p.JT_step, p.max_sizing_iterations);
                        tins(var) = tins_const;
                        Cond_h(var) = sized.Cond_h; JT(var) = sized.JT;
                        SC_w(var) = sized.SC_w;     SC_h(var) = sized.SC_h; R_J(var) = sized.R_J;
                        Ke_cavo_rad(var) = sized.Ke_rad; Ke_cavo_tor(var) = sized.Ke_tor;
                        S_CICC(var) = sized.S_CICC; S_JT(var) = sized.S_JT;

                        Ri(var) = Re(var) - Cond_h(var);
                        Re(var+1) = Ri(var) - p.INS_grades;
                        WP_w0(var) = Cond_w(var)*n_turns(var);
                    end
                end

                if n_layers > maxdim || n_layers < 0
                    reject = true; break
                end

                % Recompute B per grade
                n_spire_ = zeros(1, n_layers);
                n_spire_(1) = sum(n_turns(1:n_layers));
                for var = 2:n_layers
                    n_spire_(var) = n_spire_(var-1) - n_turns(var);
                end
                B_layers = B_TF .* (n_spire_/n_spire_(1));
                Iop = ceil(g.NI/n_spire_(1));

                % WP data
                WP_w = WP_w0(1);
                WP_h = sum(Cond_h(1:n_layers));
                A_WP = sum(Cond_h(1:n_layers).*Cond_w(1:n_layers).*n_turns(1:n_layers));
                A_SC_tot = sum(S_Cable(1:n_layers).*n_turns(1:n_layers));
                A_JT_tot = sum(S_JT(1:n_layers).*n_turns(1:n_layers));
                Ri_ = g.R_TF_Innerleg;
                % WP inner radius: cells + inter-layer insulation + ground insulation
                % (same radial stack as the Re/Ri recursion above, the section plot
                % and export_ansys_input's WPH; before, the (n_layers-1)*INS_grades
                % gaps were missing, so the nose DTF was over-estimated - by 6 mm
                % on the 13-layer design 7 checked against FEM).
                Rj_ = Ri_ - WP_h - (n_layers-1)*p.INS_grades - p.dr_plasma_side - p.GoundIns*2;

                check_w_arr = 2*Ri(1:n_layers)*tan(theta_TF/2); % maximum toroidal Case envelope
                if min((check_w_arr - (WP_w0(1:n_layers)+p.GoundIns*2))/2) < p.toroidal_gap
                    reject = true; break
                end

                % Geometric check on the obtained cable dimensions
                r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);
                if min(SC_w) <= p.min_SC_w || min(r_cable(1:n_layers)) < p.min_cable_aspect_ratio || ...
                        max(r_cable(1:n_layers)) > max_cable_aspect_ratio || min(JT(1:n_layers)) < p.min_JT
                    reject = true; break
                end

                % Hot-spot temperature check (CICC): the allowable depends on
                % the cable type (LTS/HTS), so it cannot be folded into a
                % single scalar threshold.
                ths_exceeded = false;
                for var = 1:n_layers
                    if strcmp(type_cable{var}, 'HTS')
                        ths_limit = p.THS_max_HTS;
                    else
                        ths_limit = p.THS_max_LTS;
                    end
                    if THS(var) > ths_limit
                        ths_exceeded = true;
                        break
                    end
                end
                if ths_exceeded
                    reject = true; break
                end

                % Primary radial stress (Pm+Pb) - evaluated at every layer,
                % keeping the worst case, instead of only the last layer.
                %
                % NOTE: even checking every layer, a lumped stiffness-network
                % model like this one cannot reproduce the local bending stress
                % a true 2D FEM shows at a grade transition (see
                % validation/TF_FEM_benchmark_2026_findings.md: the FEM Jacket
                % peak sits at the grade1/grade2 row boundary, roughly 2x higher
                % than this formula predicts there). p.SCF_transition_provisional
                % is an explicit, clearly-flagged empirical multiplier applied
                % only to layers adjacent to a grade change, calibrated against
                % that single FEM benchmark point - a placeholder for the
                % physics-based local-bending correction still to be developed,
                % not a validated general law. Revisit once more FEM points are
                % available.
                p_rs = B_grade(1)^2/(2*Mu_0);     % as above
                is_transition = false(1, n_layers);
                for k = 2:numel(jump_grade)
                    is_transition(max(jump_grade(k)-1, 1)) = true;
                    is_transition(jump_grade(k)) = true;
                end

                S_rm_per_layer = zeros(1, n_layers);
                E_cbl_per_layer = zeros(1, n_layers);
                sigma_nom = zeros(1, n_layers);
                xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778];
                yyy_LTS = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85];
                yyy_HTS = [1.01,1.02,1.02,1.04,1.06,1.07,1.11,1.13,1.16,1.18,1.22,1.27,1.29,1.29,1.31];
                for var = 1:n_layers
                    if strcmp(type_cable{var}, 'HTS')
                        E_cbl_var = p.E_cbl_HTS;
                        yyy = yyy_HTS;
                    else
                        E_cbl_var = p.E_cbl_LTS;
                        yyy = yyy_LTS;
                    end
                    param = Cond_w(var)/SC_w(var);
                    if param > 1 && param < 2.8
                        pf = polyfit(xxx, yyy, 5);
                        scf = polyval(pf, param);
                    else
                        scf = 1.5;
                    end

                    K_jckt = 2*p.E_jckt*JT(var)/Cond_h(var);
                    dcr_jckt = K_jckt/Ke_cavo_rad(var);
                    r_steel = (Cond_w(var)-2*tins(var))/(2*JT(var));
                    S_rm_var = p_rs*r_steel*scf*dcr_jckt; % Radial membrane stress, this Jacket layer
                    if is_transition(var)
                        S_rm_var = S_rm_var*p.SCF_transition_provisional;
                    end
                    S_rm_per_layer(var) = S_rm_var;
                    E_cbl_per_layer(var) = E_cbl_var;
                    sigma_nom(var) = p_rs*r_steel*dcr_jckt;   % nominal radial stress, no SCF (surrogate input)
                end
                [S_rm, worst_var] = max(S_rm_per_layer);
                E_cbl = E_cbl_per_layer(worst_var);

                % Case/vault sizing (fix #1 + fix #4 live in size_case_vault.m)
                ctx = struct();
                ctx.Rj_ = Rj_; ctx.Ri_ = Ri_; ctx.CASE_w = env.CASE_w; ctx.theta_TF = theta_TF;
                ctx.A_WP = A_WP; ctx.A_SC_tot = A_SC_tot; ctx.A_JT_tot = A_JT_tot;
                ctx.WP_h = WP_h; ctx.lateral_w = lateral_w;
                ctx.E_case = p.E_case; ctx.E_cbl = E_cbl; ctx.E_jckt = p.E_jckt;
                ctx.p_rs = p_rs; ctx.S_rm = S_rm;
                ctx.n_TF = p.n_TF; ctx.RTFo = g.RTFo; ctx.RTFi = g.RTFi;
                ctx.n_spire1 = n_spire_(1); ctx.Iop = Iop; ctx.Mu_0 = Mu_0;
                ctx.dr_plasma_side = p.dr_plasma_side;
                ctx.S_amm_VT = p.S_amm_VT; ctx.S_amm_JT = p.S_amm_JT; ctx.safety_membrane = p.safety_membrane;
                ctx.Ke_cavo_rad = Ke_cavo_rad(1:n_layers); ctx.n_turns = n_turns(1:n_layers);
                ctx.DTF0 = p.DTF_initial; ctx.DTF_step = p.DTF_step; ctx.max_iter = p.max_sizing_iterations;
                % with the jacket stress surrogate the jacket is checked
                % below, layer by layer, and does not drive the nose
                ctx.jacket_in_loop = scf_model == 0;

                cv = size_case_vault(ctx);
                Rk_ = cv.Rk_; S_T_VT = cv.S_T_VT; S_T_JT = cv.S_T_JT;

                JT_Pm = NaN; JT_PmPb = NaN;
                JT_crit_layer = worst_var;      % critical jacket layer (analytic formula)
                if scf_model == 1
                    % jacket thickness of every grade sized with the fast
                    % surrogate calibrated on the 2D FE (SIZE_JACKET_SURROGATE,
                    % JACKET_STRESS_SURROGATE): smallest JT per grade with
                    % jacket_margin x Pm <= Sm and jacket_margin x (Pm+Pb) <= 1.5 Sm.
                    % A thicker jacket changes the WP height, its stiffness and
                    % the axial stress S_z, so the candidate is sized again with
                    % the new JT (next pass of this loop) until the jacket
                    % no longer grows; then the field check below runs as usual.
                    B_lay = B_size_layer([1:min(n_layers, numel(B_size_layer)), ...
                        numel(B_size_layer)*ones(1, n_layers - numel(B_size_layer))]);
                    gidx = numel(jump_grade)*ones(1, n_layers);
                    for ig = numel(jump_grade):-1:1
                        gidx(jump_grade(ig):min(grade_end(ig), n_layers)) = ig;
                    end
                    E_cbl_l = zeros(1, n_layers);
                    for var = 1:n_layers
                        E_cbl_l(var) = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);
                    end
                    js = size_jacket_surrogate(struct('A_cable', S_Cable(1:n_layers), 'Cond_w', Cond_w(1:n_layers), ...
                        'tins', tins_const, 'E_cbl', E_cbl_l, 'n_turns', n_turns(1:n_layers), 'B_layer', B_lay, ...
                        'grade', gidx, 'JT0', JT(1:n_layers), 'E_jckt', p.E_jckt, 'E_ins', p.E_ins, ...
                        'r_SC_min', p.r_SC_min, 'r_SC_max', p.r_SC_max, 'p_rs', p_rs, 'Iop', Iop, 'S_z', cv.S_z, ...
                        'Sm', Sm_jacket, 'margin', jacket_margin, 'JT_step', p.JT_step, 'JT_max', JT_max));
                    if ~js.ok
                        n_jacket_rejected = n_jacket_rejected + 1;
                        reject = true; break
                    end
                    if any(js.JT > JT(1:n_layers) + 1e-9)
                        % thicker jacket needed: size the candidate again with it
                        JT_req(1:n_layers) = max(JT_req(1:n_layers), js.JT);
                        JT_req(n_layers+1:end) = JT_req(n_layers);
                        jacket_passes = jacket_passes + 1;
                        continue
                    end
                    JT_Pm = max(js.Pm); JT_PmPb = max(js.PmPb); JT_crit_layer = js.crit_layer;
                    S_T_JT = JT_Pm;
                end

                if ~(S_T_VT < p.S_amm_VT && S_T_JT < p.S_amm_JT && S_T_JT > 0 && S_T_VT > 0)   % vault vs its own allowable (review C07)
                    reject = true; break
                end

                % Real peak field of the sized candidate vs the field each
                % layer's cable was sized for (B_grade); smeared model only:
                % accept as before
                if ~use_discrete_field
                    B_peak_layer = nan(1, n_layers);  % not computed
                    field_ok = true; break
                end
                fr = struct('Iop', Iop, 'n_layers', n_layers, 'n_turns', n_turns(1:n_layers), ...
                    'Cond_w', Cond_w(1:n_layers), 'Cond_h', Cond_h(1:n_layers), 'JT', JT(1:n_layers), ...
                    'Ri_', Ri_, 'shape_cable', p.shape_cable);
                B_peak_layer = wp_peak_field_fast(fr, p);
                % converged: every grade sized for its own peak (within
                % field_tol, neither under- nor over-sized) and no layer
                % above the field its cable was sized for
                ge = max(jump_grade, [jump_grade(2:end)-1, n_layers]);
                dB = zeros(1, numel(jump_grade));
                for ig = 1:numel(jump_grade)
                    dB(ig) = max(B_peak_layer(jump_grade(ig):ge(ig))) - B_grade(jump_grade(ig));
                end
                if all(abs(dB) <= field_tol) && all(B_peak_layer <= B_grade(1:n_layers) + field_tol)
                    field_ok = true; break
                end
                B_field_layer = B_peak_layer;
            end
            if reject || ~field_ok
                n_field_rejected = n_field_rejected + (~reject);
                continue
            end
            counter = counter + 1;
            n_cond = n_spire_(1);
            WP_w = WP_w0(1);
            JENG = Iop/(min(Cond_w(1:n_layers))*min(Cond_h(1:n_layers)))*1e-6;
            E = 0.5*L*Iop^2*1e-6;
            S_T_VT = S_T_VT*1e-6;
            S_T_JT = S_T_JT*1e-6;
            JT_Pm = JT_Pm*1e-6; JT_PmPb = JT_PmPb*1e-6;   % [MPa], surrogate jacket stresses (NaN with scf_model = 0)
            Nose = Rj_-Rk_;
            R_0 = p.R0;
            radial_build = Ri_-Rk_;

            % shape_cable is saved with the solution (200 RIS / 201 Rect): every
            % downstream step (section plots, FEM surrogate, ANSYS export)
            % reads the conductor shape from the row, so a reloaded RIS
            % design cannot silently become Rect under a different input file.
            shape_cable = p.shape_cable;
            % peak field on the conductor (discrete model; NaN with the
            % smeared model), per layer padded to maxdim like n_turns, and
            % the number of sizing/field passes
            B_peak = max(B_peak_layer);
            B_peak_layers = nan(1, maxdim); B_peak_layers(1:n_layers) = B_peak_layer;
            % calibrated estimate of the first pass (NaN unless field_model = 1)
            B_cal = max(B_cal_layer);
            field_iter = field_it;
            row = table(S_T_VT,S_T_JT,R_0,g.B_PHI_0,B_TF,Iop,JENG,L,E,Ri_,Rj_,Rk_,radial_build,Nose,WP_h,WP_w,...
                lateral_w,n_cond,n_layers,n_turns,type_cable,Cond_w,Cond_h,JT,r_cable,N_Sc,N_Cu,S_Cable,S_REBCO,S_Cu_HTS,THS,B_grade,Tau_discharge, ...
                shape_cable, B_peak, B_peak_layers, field_iter, B_cal, JT_Pm, JT_PmPb, JT_crit_layer, ...
                'VariableNames', {'S_T_VT','S_T_JT','R_0','B_PHI_0','B_TF','Iop','JENG','L','E','Ri_','Rj_','Rk_', ...
                'radial_build','Nose','WP_h','WP_w','lateral_w','n_cond','n_layers','n_turns','type_cable','Cond_w', ...
                'Cond_h','JT','r_cable','N_Sc','N_Cu','S_Cable','S_REBCO','S_Cu_HTS','THS','B_grade','Tau_discharge', ...
                'shape_cable','B_peak','B_peak_layers','field_iter','B_cal','JT_Pm','JT_PmPb','JT_crit_layer'});
            DATA(counter,:) = row;
        end
    end
end

if progress_msg_len > 0 && examined < total_candidates
    fprintf('\n'); % make sure the cursor isn't left mid-progress-line on an early return path
end

if use_discrete_field && n_field_rejected > 0
    fprintf(['%d candidate(s) passed every check but their grade fields did not converge to the ' ...
        'discrete peak within %.3g T in %d passes: rejected.\n'], n_field_rejected, field_tol, field_max_iter);
end
if scf_model == 1 && n_jacket_rejected > 0
    fprintf('%d candidate(s) rejected: no jacket thickness up to JT_max = %.1f mm satisfies the surrogate criteria.\n', ...
        n_jacket_rejected, 1e3*JT_max);
end
if counter == 0
    DATA = table();
else
    % ranking: smallest radial build of the inner leg (Ri_ - Rk_) first, so
    % design #1 is the best starting point
    [~, order] = sort(DATA.radial_build);
    DATA = DATA(order, :);
end
end

function s = format_hms(seconds)
%FORMAT_HMS Render a duration in seconds as HH:MM:SS for the progress line.
if ~isfinite(seconds), seconds = 0; end
seconds = max(0, round(seconds));
h = floor(seconds/3600); m = floor(mod(seconds,3600)/60); sec = mod(seconds,60);
s = sprintf('%02d:%02d:%02d', h, m, sec);
end

function jump_grade = pick_jump_grades(n_grades, B_layers, target2, target3)
%PICK_JUMP_GRADES Grade indices at which the cable type/current design changes.
if n_grades == 3
    jump_grade = zeros(1,3);
    [~, jump_grade(1)] = min(abs(B_layers - B_layers(1)));
    [~, jump_grade(2)] = min(abs(B_layers - target2));
    [~, jump_grade(3)] = min(abs(B_layers - target3));
elseif n_grades == 2
    jump_grade = zeros(1,2);
    [~, jump_grade(1)] = min(abs(B_layers - B_layers(1)));
    [~, jump_grade(2)] = min(abs(B_layers - target3));
else
    error('pick_jump_grades:unsupported_n_grades', 'n_grades must be 2 or 3 (got %g).', n_grades);
end
end

function E_cbl = pick_E_cbl(cable_type, E_cbl_HTS, E_cbl_LTS)
%PICK_E_CBL Equivalent cable Young's modulus for the given cable type.
if strcmp(cable_type, 'HTS')
    E_cbl = E_cbl_HTS;
else
    E_cbl = E_cbl_LTS;
end
end

function sized = size_grade_cable(Cond_w, S_Cable_var, r_SC_min, r_SC_max, tins, E_jckt, E_cbl, E_ins, ...
    shape_cable, p_rs, S_z_JT, S_amm_JT, safety_membrane, min_JT, JT_step, max_iter)
%SIZE_GRADE_CABLE Thin convenience wrapper around SIZE_CICC_CABLE.
in.Cond_w = Cond_w; in.S_Cable_var = S_Cable_var;
in.r_SC_min = r_SC_min; in.r_SC_max = r_SC_max; in.tins = tins;
in.E_jckt = E_jckt; in.E_cbl = E_cbl; in.E_ins = E_ins; in.shape_cable = shape_cable;
in.p_rs = p_rs; in.S_z_JT = S_z_JT; in.S_amm_JT = S_amm_JT; in.safety_membrane = safety_membrane;
in.min_JT = min_JT; in.JT_step = JT_step; in.max_iter = max_iter;
sized = size_cicc_cable(in);
end
