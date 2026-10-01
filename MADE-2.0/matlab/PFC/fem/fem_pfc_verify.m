function res = fem_pfc_verify(row, p, geom, opts)
%FEM_PFC_VERIFY Axisymmetric FE verification of a PFC design point against the analytic design models.
%
%   res = FEM_PFC_VERIFY(row, p, geom, opts) takes one design point from
%   SEARCH/SCAN_WP_DESIGNS (a one-row results table), the machine
%   parameters p, the combined CS+PFC+plasma geometry/scenario struct geom
%   (READ_COIL_GEOMETRY) and checks the solution with an independent model:
%
%   1. Field: every turn of the coil is a uniform-current block (n x n
%      filaments, p.fem_subfil), and the CS, the other PF coils and the
%      plasma ring (geom rows, block dr x dz, p.fem_bg_fil_r x
%      p.fem_bg_fil_z filaments each) are added with the currents of every
%      plasma scenario (columns of geom.MAt_scenario).
%   2. Loads: Lorentz body force density J x B at every Gauss point,
%      f_r = J*Bz, f_z = -J*Br (J = turn current / conductor cell area).
%   3. Mechanics: axisymmetric FE (FEM_AXISYM_SOLVE) of the winding pack,
%      homogenized orthotropic material (FEM_WP_MATERIAL). Vertical support
%      from p.fem_bc (1: inertia relief, 2: bottom face).
%   4. Comparison with the design models: hoop/vertical/Tresca stress of
%      EQV_STRESS_COIL_RING_CICC (legacy and with the vertical-area fix) and
%      peak field/forces of EMAG_FIELD_FORCES.
%
%   Load cases: case 0 reproduces the design assumptions (this coil alone,
%   at the design current Iop, no background); cases 1..ns are the plasma
%   scenarios with the coil current scaled to the scenario Ampere-turns.
%
%   opts (optional struct): .coil_row  row of the coil in geom (default
%   p.cs_modules + row.n_PF); .plot (default false); .scenarios  subset of
%   scenario columns (default all).
%
%   res fields: wp, mat, mesh, cases (struct array with field/force/stress
%   summaries per load case), analytic (design-model results), table
%   (text summary printed to the command window).

if nargin < 4, opts = struct(); end
p = fem_defaults(p);
wp = fem_wp_geometry(row, p);
mat = fem_wp_material(wp, p);

coil_row = p.cs_modules + wp.n_PF;
if isfield(opts, 'coil_row'), coil_row = opts.coil_row; end
do_plot = isfield(opts, 'plot') && opts.plot;
nscen_all = size(geom.MAt_scenario, 2);
scen = 1:nscen_all;
if isfield(opts, 'scenarios'), scen = opts.scenarios; end

%% Mesh of the turns region (local z from 0 to H, coil mid-plane at Z(coil_row))
mr = p.fem_mesh_r; mz = p.fem_mesh_z;
nr = wp.n_l*mr; nz = wp.n_t*mz;
mesh = fem_mesh_rect(wp.Ri, wp.Re, 0, wp.H, nr, nz);
ne = size(mesh.elems, 1);
Zc_coil = geom.Z(coil_row);
zshift = Zc_coil - wp.H/2;

% Gauss point coordinates (global z)
gp = [-1 -1; 1 -1; 1 1; -1 1]/sqrt(3);
rg = zeros(ne, 4); zg = zeros(ne, 4);
for g = 1:4
    xi = gp(g,1); eta = gp(g,2);
    N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    for e = 1:ne
        xe = mesh.nodes(mesh.elems(e,:), :);
        rg(e,g) = N*xe(:,1);
        zg(e,g) = N*xe(:,2) + zshift;
    end
end

%% Filaments: self (sub-filaments of every turn) + background coils
ns = p.fem_subfil;
[Rf_s, Zf_s] = turn_filaments(wp, ns, zshift);
n_self = numel(Rf_s);
n_per_turn = ns*ns;

bg_rows = setdiff(1:numel(geom.R), coil_row);
nbr = p.fem_bg_fil_r; nbz = p.fem_bg_fil_z;
Rf_b = []; Zf_b = []; bg_id = [];
for k = bg_rows
    rr = geom.R(k) + ((1:nbr)-0.5)/nbr*geom.dr(k) - geom.dr(k)/2;
    zz = geom.Z(k) + ((1:nbz)-0.5)/nbz*geom.dz(k) - geom.dz(k)/2;
    [RR, ZZ] = meshgrid(rr, zz);
    Rf_b = [Rf_b; RR(:)]; Zf_b = [Zf_b; ZZ(:)]; %#ok<AGROW>
    bg_id = [bg_id; k*ones(numel(RR), 1)]; %#ok<AGROW>
end

ncase = 1 + numel(scen);
Ifil = zeros(n_self + numel(Rf_b), ncase);
case_name = cell(1, ncase);
I_turn = zeros(1, ncase);
Ifil(1:n_self, 1) = wp.Iop/n_per_turn;
I_turn(1) = wp.Iop;
case_name{1} = 'design (self, Iop)';
for c = 1:numel(scen)
    s = scen(c);
    I_turn(c+1) = geom.MAt_signed(coil_row, s)/wp.N_spire;
    Ifil(1:n_self, c+1) = I_turn(c+1)/n_per_turn;
    for k = bg_rows
        idx = n_self + find(bg_id == k);
        Ifil(idx, c+1) = geom.MAt_signed(k, s)/numel(idx);
    end
    case_name{c+1} = sprintf('scenario %d', s);
end

Rp = rg(:); Zp = zg(:);
[BRg, BZg] = fem_field_cases(Rp, Zp, [Rf_s; Rf_b], [Zf_s; Zf_b], Ifil);

%% Solve every load case
bc_mode = p.fem_bc;
if bc_mode == 1
    bc = struct('type', 'inertia_relief');
else
    bc = struct('type', 'dofs', 'dofs', 2*mesh.bottom');
end

tc_layer = ceil(mesh.ir/mr);
cases = struct([]);
for c = 1:ncase
    J = I_turn(c)/wp.A_cell;                       % [A/m^2]
    fgp = zeros(ne, 4, 2);
    BRc = reshape(BRg(:, c), ne, 4); BZc = reshape(BZg(:, c), ne, 4);
    fgp(:,:,1) = J*BZc;
    fgp(:,:,2) = -J*BRc;
    sol = fem_axisym_solve(mesh.nodes, mesh.elems, mat.D, fgp, bc);

    sth_avg = sol.sel(:,3); sz_avg = sol.sel(:,2); sr_avg = sol.sel(:,1);
    s_th = sth_avg*mat.f_hoop;                       % steel (jacket) hoop stress
    s_z = sz_avg*mat.f_z;                            % jacket-wall vertical stress
    s_r = sr_avg;                                    % smeared radial stress
    S_T_sum = abs(s_th + s_z);                       % design-model definition
    S_T_tresca = zeros(ne, 1); S_vm = zeros(ne, 1);
    for e = 1:ne
        sg = [s_th(e), s_z(e), s_r(e)];
        S_T_tresca(e) = max(sg) - min(sg);
        S_vm(e) = sqrt(0.5*((sg(1)-sg(2))^2 + (sg(2)-sg(3))^2 + (sg(3)-sg(1))^2));
    end

    layer_hoop = zeros(1, wp.n_l);
    for l = 1:wp.n_l
        sl = s_th(tc_layer == l);
        [~, im] = max(abs(sl));
        layer_hoop(l) = sl(im);   % signed value of largest magnitude (tension +, compression -)
    end

    % global equilibrium check: hoop tension of the FE section vs the radial load
    T_fe = sum(sth_avg.*sol.vol./(2*pi*sol.rcen));     % integral of sigma_theta over the section [N]
    T_load = sol.Fnet(1)/(2*pi);                        % hoop tension required by the radial load [N]

    cases(c).name = case_name{c};
    cases(c).I_turn = I_turn(c);
    cases(c).Bz_max = max(mean(BZc, 2)); cases(c).Bz_min = min(mean(BZc, 2));
    cases(c).B_peak = max(sqrt(mean(BRc, 2).^2 + mean(BZc, 2).^2)); % peak of the cell-averaged field (point values near a filament are noisy)
    cases(c).Fr_MN = sol.Fnet(1)*1e-6;
    cases(c).Fz_net_MN = sol.Fnet(2)*1e-6;
    cases(c).react_MN = sol.react*1e-6;
    [~, ih] = max(abs(s_th));
    cases(c).hoop_max_MPa = s_th(ih)*1e-6;      % signed hoop stress of largest magnitude (compression < 0)
    cases(c).hoop_layer_MPa = layer_hoop*1e-6;
    cases(c).vert_max_MPa = max(abs(s_z))*1e-6;
    cases(c).radial_avg_max_MPa = max(abs(s_r))*1e-6;
    cases(c).S_T_sum_MPa = max(S_T_sum)*1e-6;
    cases(c).S_T_tresca_MPa = max(S_T_tresca)*1e-6;
    cases(c).S_vm_MPa = max(S_vm)*1e-6;
    cases(c).Fz_relief = sol.Fz_relief;
    cases(c).u_max_mm = max(abs(sol.u))*1e3;
    cases(c).T_hoop_MN = T_fe*1e-6;
    cases(c).T_check = T_fe/T_load;
    cases(c).sol = sol;
    cases(c).s_th = s_th; cases(c).s_z = s_z; cases(c).s_r = s_r;
    cases(c).BR = BRc; cases(c).BZ = BZc;
end

%% Design-model results for the same geometry
Ri_sz = wp.Ri; 
[Bsum, Bmin, Bmax, FRmax, FZmax, ~, ~, Fr_leg, Fz_leg] = emag_field_forces( ...
    wp.Cond_h, wp.Cond_w, Ri_sz, wp.n_t, wp.n_l, 1, wp.Iop); %#ok<ASGLU>
tc = {wp.type_cable};
ring_leg = struct('E_jckt', p.E_jckt, 'E_cbl_LTS', p.E_cbl_LTS, 'E_cbl_HTS', p.E_cbl_HTS, ...
    'E_ins', p.E_ins, 'ni', p.ring_nu, 'Fz_area_fix', 0);
ring_fix = ring_leg; ring_fix.Fz_area_fix = 1;
WP_h_coil = wp.H + 2*p.grins_h;
[Sh_a, ~, Sv_leg, ST_leg] = eqv_stress_coil_ring_cicc(wp.FZ_ext, wp.Re, wp.Ri, wp.Cond_h, wp.Cond_w, wp.JT, ...
    wp.SC_h, wp.SC_w, wp.tins, tc, wp.S_CICC, wp.S_JT, 1, wp.Iop, Bmin, Bmax, WP_h_coil, wp.n_l, ring_leg);
[~, ~, Sv_fix, ST_fix] = eqv_stress_coil_ring_cicc(wp.FZ_ext, wp.Re, wp.Ri, wp.Cond_h, wp.Cond_w, wp.JT, ...
    wp.SC_h, wp.SC_w, wp.tins, tc, wp.S_CICC, wp.S_JT, 1, wp.Iop, Bmin, Bmax, WP_h_coil, wp.n_l, ring_fix);

% vertical load = the self-field squeeze computed by EMAG_FIELD_FORCES (fz_source = 1), annulus-area form
ring_sq = ring_fix;
[~, ~, Sv_sq, ST_sq] = eqv_stress_coil_ring_cicc(abs(FZmax)*1e6, wp.Re, wp.Ri, wp.Cond_h, wp.Cond_w, wp.JT, ...
    wp.SC_h, wp.SC_w, wp.tins, tc, wp.S_CICC, wp.S_JT, 1, wp.Iop, Bmin, Bmax, WP_h_coil, wp.n_l, ring_sq);
analytic.S_ver_squeeze_MPa = max(Sv_sq)*1e-6;
analytic.S_T_squeeze_MPa = max(ST_sq)*1e-6;

analytic.Bsum = Bsum; analytic.Bmin = Bmin; analytic.Bmax = Bmax;
analytic.FRmax_MN = FRmax; analytic.FZmax_half_MN = FZmax;
analytic.FR_total_MN = sum(Fr_leg)*1e-6; analytic.FZ_total_MN = sum(Fz_leg)*1e-6;
analytic.hoop_layer_MPa = Sh_a*1e-6;
analytic.hoop_max_MPa = max(abs(Sh_a))*1e-6;
analytic.S_ver_legacy_MPa = max(Sv_leg)*1e-6;
analytic.S_ver_fixed_MPa = max(Sv_fix)*1e-6;
analytic.S_T_legacy_MPa = max(ST_leg)*1e-6;
analytic.S_T_fixed_MPa = max(ST_fix)*1e-6;
analytic.FZ_ext_MN = wp.FZ_ext*1e-6;

res.wp = wp; res.mat = mat; res.mesh = mesh; res.cases = cases; res.analytic = analytic;
res.opts = struct('coil_row', coil_row, 'bc', bc_mode, 'scenarios', scen);

print_summary(res, p);
if do_plot
    fem_plot_results(res, 1);
end
end

function [Rf, Zf] = turn_filaments(wp, ns, zshift)
% ns x ns uniform sub-filaments in every conductor cell of the winding pack
Rf = zeros(wp.n_l*wp.n_t*ns*ns, 1); Zf = Rf;
k = 0;
off = ((1:ns)-0.5)/ns;
for l = 1:wp.n_l
    for t = 1:wp.n_t
        r0 = wp.Ri + (l-1)*wp.Cond_w;
        z0 = (t-1)*wp.Cond_h + zshift;
        for a = 1:ns
            for b = 1:ns
                k = k + 1;
                Rf(k) = r0 + off(a)*wp.Cond_w;
                Zf(k) = z0 + off(b)*wp.Cond_h;
            end
        end
    end
end
end

function print_summary(res, p)
wp = res.wp; an = res.analytic; cs = res.cases;
fprintf('\nFEM verification - PF%d, %d layers x %d turns, Iop = %.1f kA, JT = %.2f mm, WP %.3f x %.3f m\n', ...
    wp.n_PF, wp.n_l, wp.n_t, wp.Iop*1e-3, wp.JT*1e3, wp.Re-wp.Ri, wp.H);
fprintf('  Mesh %d x %d elements, homogenized WP: Eth=%.1f GPa  Er=%.1f GPa  Ez=%.1f GPa, jacket factors f_hoop=%.2f f_z=%.2f\n', ...
    res.mesh.nr, res.mesh.nz, res.mat.Eth*1e-9, res.mat.Er*1e-9, res.mat.Ez*1e-9, res.mat.f_hoop, res.mat.f_z);
fprintf('  Design model: Bpeak %.2f T (Bz %.2f..%.2f), hoop max %.0f MPa, FZ ext %.2f MN, emag half-stack FZ %.1f MN\n', ...
    an.Bsum, an.Bmin, an.Bmax, an.hoop_max_MPa, an.FZ_ext_MN, an.FZmax_half_MN);
fprintf('  Design-model vertical stress / Tresca [MPa]: legacy (fz_source 0) %.0f / %.0f, area fix %.0f / %.0f, emag squeeze (fz_source 1) %.0f / %.0f\n', ...
    an.S_ver_legacy_MPa, an.S_T_legacy_MPa, an.S_ver_fixed_MPa, an.S_T_fixed_MPa, an.S_ver_squeeze_MPa, an.S_T_squeeze_MPa);
fprintf('  %-22s %8s %9s %9s %9s %9s %9s %9s %9s\n', 'load case', 'I [kA]', 'Bpk [T]', 'Fr [MN]', 'Fz [MN]', 'hoop MPa', 'vert MPa', 'Tresca', 'vonMises');
    fprintf('  (hoop = signed value of largest magnitude: negative = compression, i.e. net inward radial force)\n');
for c = 1:numel(cs)
    fprintf('  %-22s %8.1f %9.2f %9.2f %9.2f %9.0f %9.0f %9.0f %9.0f\n', cs(c).name, cs(c).I_turn*1e-3, cs(c).B_peak, ...
        cs(c).Fr_MN, cs(c).Fz_net_MN, cs(c).hoop_max_MPa, cs(c).vert_max_MPa, cs(c).S_T_tresca_MPa, cs(c).S_vm_MPa);
end
fprintf('  (hoop/vertical: jacket stress recovered from the smeared FE stress; Tresca/von Mises from [hoop, vertical, radial])\n');
end
