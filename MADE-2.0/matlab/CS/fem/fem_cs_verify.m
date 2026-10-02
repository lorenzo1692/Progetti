function res = fem_cs_verify(row, p, g, geom, opts)
%FEM_CS_VERIFY Axisymmetric FE verification of a CS design point: the whole stack of modules.
%
%   res = FEM_CS_VERIFY(row, p, g, geom, opts) takes one design point from
%   SEARCH/SCAN_WP_DESIGNS (a one-row results table), the machine
%   parameters p, the derived parameters g (COMPUTE_OPERATING_PARAMS) and,
%   optionally, the machine geometry / scenario struct geom, and checks the
%   solution with an independent model of the complete stack:
%
%   1. Geometry: p.n_moduli identical modules, each with ground
%      insulation, separated by steel spacer plates (FEM_CS_STACK). All
%      modules are homogenized orthotropic (FEM_WP_MATERIAL), so the
%      stack carries the vertical load path through the plates; a single
%      module can later be replaced by a turn-by-turn model.
%   2. Field: every turn of every module is a block of n x n filaments
%      (p.fem_subfil) carrying its module current; with geom the PF coils
%      and the plasma ring (geom rows after the CS modules) are added as
%      blocks of p.fem_bg_fil_r x p.fem_bg_fil_z filaments.
%   3. Loads: Lorentz body force density J x B at every Gauss point of the
%      winding-pack elements (f_r = J*Bz, f_z = -J*Br).
%   4. Mechanics: axisymmetric FE (FEM_AXISYM_SOLVE). Support p.fem_bc: 1
%      inertia relief, 2 stack bottom face fixed axially, optionally with a
%      stack preload p.fem_preload_MN pushing down on the top face.
%   5. Comparison with the design models: hoop / vertical / Tresca stress
%      of EQV_STRESS_COIL_CICC (with the Fz stored in the design point) and
%      the field / forces of EMAG_FIELD_FORCES for the full stack.
%
%   Load cases: case 1 reproduces the design assumptions (this stack alone,
%   all modules at the design current Iop, no background); with geom, cases
%   2..ns+1 are the plasma scenarios, module m carrying the scenario
%   Ampere-turns of geom row m (rows 1..n_moduli are the CS modules, the
%   same convention as the scenario file) divided by its turns.
%
%   Module m = p.fem_detail_module > 0 is modeled turn by turn (insulation,
%   jacket wall, cable of every conductor cell); the others stay homogenized.
%   Its hoop / vertical / Tresca are then read from the jacket elements.
%
%   opts (optional struct): .plot (default false); .scenarios subset of
%   scenario columns (default all); .zc axial center of the stack (default
%   the mean Z of the CS rows of geom, or 0).
%
%   res fields: wp, mat, stack, cases (struct array: per-case, per-module
%   loads and stresses, force through every plate), analytic (design-model
%   results). A summary is printed to the command window.

if nargin < 4, geom = []; end
if nargin < 5, opts = struct(); end
p = fem_defaults(p);
wp = fem_cs_geometry(row, p, g);
nm = wp.n_mod;

have_geom = ~isempty(geom) && numel(geom.R) >= nm;
zc = 0;
if have_geom, zc = mean(geom.Z(1:nm)); end
if isfield(opts, 'zc'), zc = opts.zc; end
do_plot = isfield(opts, 'plot') && opts.plot;
scen = [];
if have_geom
    scen = 1:size(geom.MAt_signed, 2);
    if isfield(opts, 'scenarios'), scen = opts.scenarios; end
end

st = fem_cs_stack(wp, p, zc);
mesh = st.mesh; mat = st.mat;
ne = size(mesh.elems, 1);
ew = find(st.matid == 1 | st.matid == 4);   % loaded elements: homogenized winding pack, cable of the detailed module
nw = numel(ew);
is_cable = (st.matid(ew) == 4);

%% Field points: centroid of every loaded element (the load is uniform in the element)
Rn = mesh.nodes(:, 1); Zn = mesh.nodes(:, 2);
rg = mean(reshape(Rn(mesh.elems(ew, :)), nw, 4), 2);
zg = mean(reshape(Zn(mesh.elems(ew, :)), nw, 4), 2);

%% Filaments: self (sub-filaments of every turn of every module) + background
ns = p.fem_subfil; n_per_turn = ns*ns;
[Rf_s, Zf_s, mod_s] = stack_filaments(wp, st, ns);
n_self = numel(Rf_s);

Rf_b = []; Zf_b = []; bg_id = [];
bg_rows = [];
if have_geom
    bg_rows = nm+1:numel(geom.R);
    nbr = p.fem_bg_fil_r; nbz = p.fem_bg_fil_z;
    for k = bg_rows
        rr = geom.R(k) + ((1:nbr)-0.5)/nbr*geom.dr(k) - geom.dr(k)/2;
        zz = geom.Z(k) + ((1:nbz)-0.5)/nbz*geom.dz(k) - geom.dz(k)/2;
        [RR, ZZ] = meshgrid(rr, zz);
        Rf_b = [Rf_b; RR(:)]; Zf_b = [Zf_b; ZZ(:)]; %#ok<AGROW>
        bg_id = [bg_id; k*ones(numel(RR), 1)]; %#ok<AGROW>
    end
end

ncase = 1 + numel(scen);
Ifil = zeros(n_self + numel(Rf_b), ncase);
I_mod = zeros(nm, ncase);           % turn current of every module and case
case_name = cell(1, ncase);
I_mod(:, 1) = wp.Iop;
Ifil(1:n_self, 1) = wp.Iop/n_per_turn;
case_name{1} = 'design (self, Iop)';
for c = 1:numel(scen)
    s = scen(c);
    I_mod(:, c+1) = geom.MAt_signed(1:nm, s)/wp.N_mod;
    Ifil(1:n_self, c+1) = I_mod(mod_s, c+1)/n_per_turn;
    for k = bg_rows
        idx = n_self + find(bg_id == k);
        Ifil(idx, c+1) = geom.MAt_signed(k, s)/numel(idx);
    end
    case_name{c+1} = sprintf('scenario %d', s);
end

[BRg, BZg] = fem_field_cases(rg, zg, [Rf_s; Rf_b], [Zf_s; Zf_b], Ifil);

%% Boundary conditions
if p.fem_bc == 1
    bc = struct('type', 'inertia_relief');
    if p.fem_preload_MN ~= 0
        warning('fem_cs_verify:preload_ignored', 'fem_preload_MN needs fem_bc = 2 (bottom support); ignored.');
    end
else
    bc = struct('type', 'dofs', 'dofs', 2*mesh.bottom');
    if p.fem_preload_MN ~= 0
        rt = mesh.rv(:);
        rm = [rt(1); (rt(1:end-1)+rt(2:end))/2; rt(end)];
        w = pi*(rm(2:end).^2 - rm(1:end-1).^2);
        fext = zeros(2*size(mesh.nodes, 1), 1);
        fext(2*mesh.top) = -p.fem_preload_MN*1e6*w/sum(w);
        bc.fext = fext;
    end
end

%% Solve every load case
cases = struct([]);
for c = 1:ncase
    A_cond = wp.A_cell*ones(nw, 1);                         % current area: whole cell (homogenized) or the cable rectangle
    A_cond(is_cable) = wp.SC_w*wp.SC_h;
    Jw = I_mod(st.modid(ew), c)./A_cond;                    % [A/m^2], per loaded element
    fgp = zeros(ne, 4, 2);
    BRc = repmat(BRg(:, c), 1, 4); BZc = repmat(BZg(:, c), 1, 4);
    fgp(ew, :, 1) = repmat(Jw, 1, 4).*BZc;
    fgp(ew, :, 2) = -repmat(Jw, 1, 4).*BRc;
    sol = fem_axisym_solve(mesh.nodes, mesh.elems, st.D, fgp, bc);

    sel = sol.sel;
    % stress in the steel: homogenized elements are recovered from the smeared stress,
    % the jacket elements of a detailed module are read directly
    ej = find(st.matid == 1 | st.matid == 5);
    hom = (st.matid(ej) == 1);
    s_th = sel(ej, 3); s_z = sel(ej, 2); s_r = sel(ej, 1);
    s_th(hom) = s_th(hom)*mat.f_hoop;
    s_z(hom) = s_z(hom)*mat.f_z;
    S_T_sum = abs(s_th + s_z);
    S_T_tresca = max([s_th, s_z, s_r], [], 2) - min([s_th, s_z, s_r], [], 2);
    mid = st.modid(ej);
    s_th_all = nan(ne, 1); s_th_all(ej) = s_th;

    mods = struct([]);
    for m = 1:nm
        im = (mid == m);
        [~, ih] = max(abs(s_th(im))); sh = s_th(im); 
        fm = (st.modid == m);
        vm = sol.vol(fm);
        mods(m).I_turn = I_mod(m, c);
        mods(m).Fr_MN = sum(vm.*mean(fgp(fm, :, 1), 2))*1e-6;
        mods(m).Fz_MN = sum(vm.*mean(fgp(fm, :, 2), 2))*1e-6;
        ml = (st.modid(ew) == m);
        mods(m).B_peak = max(sqrt(BRc(ml, 1).^2 + BZc(ml, 1).^2));
        mods(m).hoop_max_MPa = sh(ih)*1e-6;                % signed value of largest magnitude
        mods(m).vert_max_MPa = max(abs(s_z(im)))*1e-6;
        mods(m).S_T_sum_MPa = max(S_T_sum(im))*1e-6;
        mods(m).S_T_tresca_MPa = max(S_T_tresca(im))*1e-6;
    end

    % axial force through every spacer plate (compression positive)
    F_plate = zeros(1, nm-1);
    for k = 1:nm-1
        er = find(mesh.iz == st.iz_plate(k));
        F_plate(k) = -sum(sel(er, 2).*st.area(er))*1e-6;
    end

    T_fe = sum(sel(:, 3).*sol.vol./(2*pi*sol.rcen));       % hoop tension of the FE section [N]
    T_load = sol.Fnet(1)/(2*pi);

    cases(c).name = case_name{c};
    cases(c).modules = mods;
    cases(c).F_plate_MN = F_plate;
    cases(c).Fr_MN = sol.Fnet(1)*1e-6;
    cases(c).Fz_net_MN = sol.Fnet(2)*1e-6;
    cases(c).react_MN = sol.react*1e-6;
    cases(c).Fz_relief = sol.Fz_relief;
    [~, imx] = max(abs([mods.hoop_max_MPa])); cases(c).hoop_max_MPa = mods(imx).hoop_max_MPa;
    cases(c).vert_max_MPa = max([mods.vert_max_MPa]);
    cases(c).S_T_tresca_MPa = max([mods.S_T_tresca_MPa]);
    cases(c).u_max_mm = max(abs(sol.u))*1e3;
    cases(c).T_check = T_fe/T_load;
    % ANSYS model (STR/emag_loads/Load_import.lgw) preload variables: cumulative module Fz from the top / from the bottom
    Fzm = [mods.Fz_MN];
    cases(c).check_F_pos = max(cumsum(Fzm(end:-1:1)));
    cases(c).check_F_neg = min(cumsum(Fzm));
    cases(c).sol = sol; cases(c).s_th = s_th_all; cases(c).ej = ej;
end

%% Design-model results for the same geometry
mu0 = 4e-7*pi;
[Bsum, ~, ~, FRmax, FZmax, ~, ~, Fr_t, Fz_t] = emag_field_forces( ...
    wp.Cond_h, wp.Cond_w, wp.Ri, wp.n_t, wp.n_l, 1, nm, wp.Iop, wp.spacer, 1);
n_per = wp.N_mod;
for m = 1:nm
    rng = (m-1)*n_per+1:m*n_per;
    analytic.Fr_mod_MN(m) = sum(Fr_t(rng))*1e-6;
    analytic.Fz_mod_MN(m) = sum(Fz_t(rng))*1e-6;
end
Fz_row = rowget(row, 'Fz_MN');
PB = rowget(row, 'B_grades')^2/(2*mu0);                    % B_grades column = peak field of the scan
[Sh, ~, Sv, ST] = eqv_stress_coil_cicc(Fz_row, rowget(row, 'Ri'), g.Re, wp.Cond_h, wp.Cond_w, wp.JT, ...
    wp.SC_h, wp.SC_w, wp.tins, {wp.type_cable}, wp.S_CICC, wp.S_JT, PB);
analytic.Bsum_stack = Bsum; analytic.FRmax_MN = FRmax; analytic.FZmax_half_MN = FZmax;
analytic.Fz_design_MN = Fz_row;
analytic.hoop_max_MPa = max(Sh); analytic.S_ver_MPa = max(Sv); analytic.S_T_MPa = max(ST);
analytic.B_scan = rowget(row, 'B_grades');

res.wp = wp; res.mat = mat; res.stack = st; res.cases = cases; res.analytic = analytic;
res.opts = struct('zc', zc, 'bc', p.fem_bc, 'scenarios', scen, 'preload_MN', p.fem_preload_MN);

print_summary(res);
if do_plot
    fem_cs_plot(res, 1);
end
end

function [Rf, Zf, mod_id] = stack_filaments(wp, st, ns)
% ns x ns uniform sub-filaments in every conductor cell of every module
n = wp.n_mod*wp.n_l*wp.n_t*ns*ns;
Rf = zeros(n, 1); Zf = Rf; mod_id = Rf;
off = ((1:ns)-0.5)/ns;
k = 0;
for m = 1:wp.n_mod
    for l = 1:wp.n_l
        for t = 1:wp.n_t
            r0 = wp.Ri + (l-1)*wp.Cond_w;
            z0 = st.z_turns0(m) + (t-1)*wp.Cond_h;
            for a = 1:ns
                for b = 1:ns
                    k = k + 1;
                    Rf(k) = r0 + off(a)*wp.Cond_w;
                    Zf(k) = z0 + off(b)*wp.Cond_h;
                    mod_id(k) = m;
                end
            end
        end
    end
end
end

function v = rowget(row, name)
v = row.(name);
if iscell(v), v = v{1}; end
end

function print_summary(res)
wp = res.wp; an = res.analytic; cs = res.cases; nm = wp.n_mod;
fprintf('\nFEM verification - CS stack of %d modules, %d layers x %d turns, Iop = %.1f kA, JT = %.2f mm, turns %.3f x %.3f m per module\n', ...
    nm, wp.n_l, wp.n_t, wp.Iop*1e-3, wp.JT*1e3, wp.Re-wp.Ri, wp.H);
if res.stack.detail > 0
    fprintf('  Module %d modeled turn by turn (jacket stress read directly); the other modules are homogenized.\n', res.stack.detail);
end
fprintf('  Mesh %d x %d elements, homogenized WP: Eth=%.1f GPa  Er=%.1f GPa  Ez=%.1f GPa, jacket factors f_hoop=%.2f f_z=%.2f\n', ...
    res.stack.mesh.nr, res.stack.mesh.nz, res.stack.mat.Eth*1e-9, res.stack.mat.Er*1e-9, res.stack.mat.Ez*1e-9, res.stack.mat.f_hoop, res.stack.mat.f_z);
fprintf('  Design model: hoop %.0f MPa, vertical %.0f MPa (Fz = %.2f MN), Tresca %.0f MPa, B peak %.2f T; emag full stack: B %.2f T\n', ...
    an.hoop_max_MPa, an.S_ver_MPa, an.Fz_design_MN, an.S_T_MPa, an.B_scan, an.Bsum_stack);
fprintf('  Design model force per module [MN]  Fr: %s   Fz: %s\n', sprintf('%.2f ', an.Fr_mod_MN), sprintf('%.2f ', an.Fz_mod_MN));
for c = 1:numel(cs)
    fprintf('\n  Case: %s   (FE hoop-tension check %.3f, max displacement %.2f mm)\n', cs(c).name, cs(c).T_check, cs(c).u_max_mm);
    fprintf('  %-8s %8s %9s %9s %9s %9s %9s %9s\n', 'module', 'I [kA]', 'Bpk [T]', 'Fr [MN]', 'Fz [MN]', 'hoop MPa', 'vert MPa', 'Tresca');
    for m = 1:nm
        md = cs(c).modules(m);
        fprintf('  %-8d %8.1f %9.2f %9.2f %9.2f %9.0f %9.0f %9.0f\n', m, md.I_turn*1e-3, md.B_peak, md.Fr_MN, md.Fz_MN, ...
            md.hoop_max_MPa, md.vert_max_MPa, md.S_T_tresca_MPa);
    end
    fprintf('  Compression through the plates (bottom to top) [MN]: %s\n', sprintf('%.2f ', cs(c).F_plate_MN));
    fprintf('  ANSYS CHECK_F definitions (cumulative module Fz from top, max / from bottom, min) [MN]: %.2f / %.2f\n', cs(c).check_F_pos, cs(c).check_F_neg);
end
fprintf('\n  (module Fr/Fz: Lorentz load on the module; hoop: signed value of largest magnitude, negative = compression;\n');
fprintf('   hoop/vertical: jacket stress recovered from the smeared FE stress; plates: axial force through each spacer, positive = compression)\n');
end
