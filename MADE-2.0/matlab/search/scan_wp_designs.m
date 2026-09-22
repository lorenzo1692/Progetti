function DATA = scan_wp_designs(p, g, env, combT)
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
%   This function is a direct, line-by-line port of the scan loop from
%   WP_TF_VNS_Design_Point_2026.m: it calls the external CICC(...) sizing
%   function exactly as before, and keeps the same array layout and
%   control flow, only delegating the three duplicated jacket-sizing
%   while-loops to SIZE_CICC_CABLE and the case/vault sizing while-loop to
%   SIZE_CASE_VAULT.

maxdim = p.maxdim;
theta_TF = g.theta_TF;
Mu_0 = g.Mu_0;

tins_const = p.turn_insulation_nominal*p.Increm;   % Turn insulation

counter = 0;
% DATA is intentionally left undefined here: like the original script, it
% is created on the first successful candidate via the indexed assignment
% below (DATA(1,:) = row), which lets MATLAB infer its column structure
% from that first row. If no candidate ever succeeds, it is set to an
% empty table just before returning (see bottom of this function).

for lateral_w = env.lateral_w_min:p.lateral_w_step:env.lateral_w_max
    for i = 1:numel(combT)
        for j = 1:size(combT{i}, 1)

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

            jump_grade = pick_jump_grades(p.n_grades, B_layers, p.grade_B_target_2, p.grade_B_target_3);

            for var = jump_grade
                B = B_layers(var);
                [type_cable(var), N_Cu(var), N_Sc(var), N_tot(var), S_Cable(var), ...
                    S_REBCO(var), S_Cu_HTS(var), THS(var)] = cicc(B, Iop, Tau_discharge, p.WP_SC_type);

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
            p_rs = B_TF^2/(2*Mu_0); % Magnetic pressure, thin WP

            for var = 1:n_layers
                Cond_w(var) = WP_w0(1)/n_turns(1);
                E_cbl = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);

                sized = size_grade_cable(Cond_w(var), S_Cable(var), p.r_SC_min, p.r_SC_max, tins_const, ...
                    p.E_jckt, E_cbl, p.E_ins, p.shape_cable, p_rs, S_z_JT, ...
                    p.S_amm_JT, p.safety_membrane, p.min_JT, p.JT_step, p.max_sizing_iterations);
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
                        p.S_amm_JT, p.safety_membrane, p.min_JT, p.JT_step, p.max_sizing_iterations);
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
                continue
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
                    E_cbl = pick_E_cbl(type_cable{var}, p.E_cbl_HTS, p.E_cbl_LTS);

                    sized = size_grade_cable(Cond_w(var), S_Cable(var), p.r_SC_min, p.r_SC_max, tins_const, ...
                        p.E_jckt, E_cbl, p.E_ins, p.shape_cable, p_rs, S_z_JT, ...
                        p.S_amm_JT, p.safety_membrane, p.min_JT, p.JT_step, p.max_sizing_iterations);
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
                continue
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
            Rj_ = Ri_ - WP_h - p.dr_plasma_side - p.GoundIns*2;

            check_w_arr = 2*Ri(1:n_layers)*tan(theta_TF/2); % maximum toroidal Case envelope
            if min((check_w_arr - (WP_w0(1:n_layers)+p.GoundIns*2))/2) < p.toroidal_gap
                continue
            end

            % Geometric check on the obtained cable dimensions
            r_cable(1:n_layers) = Cond_w(1:n_layers)./Cond_h(1:n_layers);
            if min(SC_w) <= p.min_SC_w || min(r_cable(1:n_layers)) < p.min_cable_aspect_ratio || min(JT(1:n_layers)) < p.min_JT
                continue
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
            p_rs = B_TF^2/(2*Mu_0);
            is_transition = false(1, n_layers);
            for k = 2:numel(jump_grade)
                is_transition(max(jump_grade(k)-1, 1)) = true;
                is_transition(jump_grade(k)) = true;
            end

            S_rm_per_layer = zeros(1, n_layers);
            E_cbl_per_layer = zeros(1, n_layers);
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

            cv = size_case_vault(ctx);
            Rk_ = cv.Rk_; S_T_VT = cv.S_T_VT; S_T_JT = cv.S_T_JT;

            if S_T_VT < p.S_amm_JT && S_T_JT < p.S_amm_JT && S_T_JT > 0 && S_T_VT > 0
                counter = counter + 1;
                n_cond = n_spire_(1);
                WP_w = WP_w0(1);
                JENG = Iop/(min(Cond_w(1:n_layers))*min(Cond_h(1:n_layers)))*1e-6;
                E = 0.5*L*Iop^2*1e-6;
                S_T_VT = S_T_VT*1e-6;
                S_T_JT = S_T_JT*1e-6;
                Nose = Rj_-Rk_;
                R_0 = p.R0;
                radial_build = Ri_-Rk_;

                row = table(S_T_VT,S_T_JT,R_0,g.B_PHI_0,B_TF,Iop,JENG,L,E,Ri_,Rj_,Rk_,radial_build,Nose,WP_h,WP_w,...
                    lateral_w,n_cond,n_layers,n_turns,type_cable,Cond_w,Cond_h,JT,r_cable,N_Sc,N_Cu,S_Cable,S_REBCO,S_Cu_HTS,THS,Tau_discharge, ...
                    'VariableNames', {'S_T_VT','S_T_JT','R_0','B_PHI_0','B_TF','Iop','JENG','L','E','Ri_','Rj_','Rk_', ...
                    'radial_build','Nose','WP_h','WP_w','lateral_w','n_cond','n_layers','n_turns','type_cable','Cond_w', ...
                    'Cond_h','JT','r_cable','N_Sc','N_Cu','S_Cable','S_REBCO','S_Cu_HTS','THS','Tau_discharge'});
                DATA(counter,:) = row;
            end
        end
    end
end

if counter == 0
    DATA = table();
end
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
