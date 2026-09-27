function R = verify_ris_rect_fem(p, out_dir, which)
%VERIFY_RIS_RECT_FEM Reproducible RIS / Rect check of geometry, mesh and FEM.
%
%   R = VERIFY_RIS_RECT_FEM(p, out_dir) builds two design points and, for
%   each, checks the conductor-shape identity end to end:
%     'Rect' - design point 7 (13 layers x 8 turns, EM_2D007, validated
%              against ANSYS);
%     'RIS'  - a round-in-square test point sized with SIZE_CICC_CABLE
%              (shape_cable = 200): 10 layers x 7 turns, 46 mm square
%              cells, three cable areas (not a scan result: a controlled
%              fixture for the geometry/FEM path).
%   For each case it reports: the shape actually used (read from the row,
%   with p.shape_cable deliberately set to the OTHER code to prove a
%   current default cannot change it), cable/jacket areas vs sizing,
%   meshed cable area, material counts, Jacobians, contact convergence,
%   residuals (total and primary), force balance, axial balance, the
%   Lorentz resultant, and writes figures: section plot, GPS section, and
%   (RIS) a zoom on one conductor with the mesh and the contact interface.
%
%   p       - machine parameters: in MATLAB
%             p = read_machine_input('input/WP_TF_input_template.xlsx');
%   out_dir - folder for figures and the text report (created)
%   which   - optional {'Rect','RIS'} subset
%
%   Returns R.(case) with the checks and the FEM output.

if nargin < 2 || isempty(out_dir), out_dir = fullfile(pwd, 'verify_ris_rect'); end
if nargin < 3 || isempty(which), which = {'Rect', 'RIS'}; end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
fid = fopen(fullfile(out_dir, 'verify_ris_rect_report.txt'), 'w');
cleanup = onCleanup(@() fclose(fid));
say = @(varargin) fprintf_both(fid, varargin{:});

R = struct();
for w = 1:numel(which)
    name = which{w};
    [row, opts] = fixture(name, p);
    % prove the solution's own shape wins over the current default
    pp = p; pp.shape_cable = 401 - row.shape_cable;       % 200 <-> 201
    say('\n===== %s (row.shape_cable = %d, p.shape_cable deliberately = %d) =====\n', ...
        name, row.shape_cable, pp.shape_cable);
    tg = wp_turn_geometry(row, pp);
    say('shape used: %s\n', tg.shape_name);
    for k = unique([1, round(row.n_layers/2), row.n_layers])
        say(['  layer %2d: cell %.2fx%.2f mm, cable %.2fx%.2f mm r=%.2f mm, jacket outer r=%.2f mm, ' ...
             'JT_min %.2f mm, A_cable %.1f mm2 (sized %.1f), A_jacket %.1f mm2\n'], k, ...
            1e3*tg.cell_w(k), 1e3*tg.cell_h(k), 1e3*tg.cab_w(k), 1e3*tg.cab_h(k), 1e3*tg.cab_r(k), ...
            1e3*tg.jk_r(k), 1e3*tg.jt_min(k), 1e6*tg.A_cable(k), 1e6*tg.A_cable_row(k), 1e6*tg.A_jacket(k));
    end

    out = wp_mech_surrogate(row, pp, opts);
    c = out.checks;
    m = out.mesh;
    say('mesh: %d nodes, %d Q8, %d T6; materials (Q8+T6 count): jacket %d, turn/layer ins %d, filler %d, cable %d, case %d\n', ...
        size(m.xy,1), size(m.q8,1), size(m.t6,1), sum(m.q8_mat==2), sum(m.q8_mat==3)+sum(m.t6_mat==3), ...
        sum(m.q8_mat==7), sum(m.t6_mat==1)+sum(m.t6_mat==6), sum(m.t6_mat==4));
    say('checks: %s\n', c.summary);
    fn = fieldnames(c);
    for i = 1:numel(fn)
        v = c.(fn{i});
        if isnumeric(v) || islogical(v), say('  %-22s %s\n', fn{i}, mat2str(double(v), 5)); end
    end
    say('contact (total):   %s\n', contact_line(out.sol.contact));
    if ~isempty(out.primary)
        say('contact (primary): %s\n', contact_line(out.primary.contact));
    end
    say('Lorentz resultant: Fx = %.4g N/m, Fy = %.4g N/m (from the assembled load vector)\n', ...
        out.load_resultant.Fx_lorentz, out.load_resultant.Fy_lorentz);
    say('jacket: max Pm %.0f MPa, max Pm+Pb %.0f MPa, max peak %.0f MPa; eps_z %.4e\n', ...
        max(out.layer.Pm)/1e6, max(out.layer.PmPb)/1e6, max(out.layer.peak)/1e6, out.eps_z);
    say('FoM: %s (valid = %d)\n', out.fom.status, out.valid);

    % figures
    title_str = sprintf('%s verification point', name);
    plot_wp_section(row, pp, title_str);
    save_fig(gcf, fullfile(out_dir, sprintf('%s_section.png', name)));
    fg = plot_wp_gps_section(out, [title_str ' - GPS']);
    save_fig(fg, fullfile(out_dir, sprintf('%s_gps_section.png', name)));
    fz = plot_conductor_zoom(out, 1, sprintf('%s: turn 1 mesh and interfaces', name));
    save_fig(fz, fullfile(out_dir, sprintf('%s_conductor_zoom.png', name)));
    R.(name) = struct('row', row, 'tg', tg, 'checks', c, 'out', out);
end
end

function [row, opts] = fixture(name, p)
opts = struct('verbose', true);
if isfield(p, 'verify_load_sequence'), opts.load_sequence = p.verify_load_sequence; end
switch name
    case 'Rect'
        % design point 7 (EM_2D007), as validation/validate_mech_surrogate_2026.m
        row = struct('Iop', 64628, 'n_layers', 13, 'n_turns', 8*ones(1,13), ...
            'Cond_w', repmat(0.04300426499905924,1,13), ...
            'Cond_h', [repmat(0.034167368089179834,1,4) repmat(0.026104457386404455,1,3) repmat(0.024226622311504957,1,6)], ...
            'JT', repmat(0.003,1,13), 'Ri_', 1.259425287356322, 'Rk_', 0.7300827089713595);
        row.type_cable = repmat({'LTS'}, 1, 13);
        row.shape_cable = 201;
        opts.T_bf = 36.266e6;
    case 'RIS'
        nl = 10; nt = 7; Cw = 0.046;
        A = [repmat(900e-6,1,4) repmat(700e-6,1,3) repmat(600e-6,1,3)];
        tins = p.turn_insulation_nominal*p.Increm;
        Ch = zeros(1,nl); JT = zeros(1,nl);
        for k = 1:nl
            in = struct('Cond_w', Cw, 'S_Cable_var', A(k), 'r_SC_min', p.r_SC_min, 'r_SC_max', p.r_SC_max, ...
                'tins', tins, 'E_jckt', p.E_jckt*1e9, 'E_cbl', p.E_cbl_LTS*1e9, 'E_ins', p.E_ins*1e9, ...
                'shape_cable', 200, 'p_rs', 0, 'S_z_JT', 0, 'S_amm_JT', p.S_amm_JT, ...
                'safety_membrane', p.safety_membrane, 'min_JT', p.min_JT, 'JT_step', p.JT_step, 'max_iter', 1e5);
            sz = size_cicc_cable(in);
            Ch(k) = sz.Cond_h; JT(k) = sz.JT;
        end
        Ri_ = 1.259425287356322;
        WP_h = sum(Ch);
        Rj_ = Ri_ - WP_h - (nl-1)*p.INS_grades - p.dr_plasma_side - 2*p.GoundIns;
        row = struct('Iop', 96000, 'n_layers', nl, 'n_turns', nt*ones(1,nl), 'Cond_w', Cw*ones(1,nl), ...
            'Cond_h', Ch, 'JT', JT, 'Ri_', Ri_, 'Rk_', Rj_ - 0.105, 'S_Cable', A);
        row.type_cable = repmat({'LTS'}, 1, nl);
        row.shape_cable = 200;
    otherwise
        error('verify_ris_rect_fem:case', 'Unknown case %s', name);
end
g = compute_operating_params(p);
row.B_TF = g.B_PHI_TF;
end

function s = contact_line(ci)
st = ci.step(end);
s = sprintf(['steps %d, iterations %d, converged %d, last changes %d of %d pairs, dN %.2g, residual %.3g ' ...
    '(first solve %.3g, %d refinement steps), balance %.3g, contact violation force/total normal %.3g ' ...
    '(normal %.3g, friction %.3g), bonded cable pairs %d (tension fraction %.3g), max penetration %.2g mm'], ...
    numel(ci.step), ci.iter, ci.converged, st.n_changes, ci.n_pairs, st.dN, ci.residual, ...
    ci.residual_first_solve, ci.refinement_steps, ci.force_balance, ci.final_violation_force, ...
    ci.final_violation_normal, ci.final_violation_friction, ci.n_cable_bonded, ci.bond_tension_fraction, 1e3*ci.max_penetration);
end

function fprintf_both(fid, varargin)
fprintf(varargin{:});
fprintf(fid, varargin{:});
end

function save_fig(f, file)
try
    print(f, '-dpng', '-r150', file);
catch err
    warning('verify_ris_rect_fem:print', 'Could not save %s: %s', file, err.message);
end
end
