function out = wp_mech_surrogate(row, p, opts)
%WP_MECH_SURROGATE 2D finite-element equilibrium of the TF inner-leg section.
%
%   out = WP_MECH_SURROGATE(row, p) solves, from the design point alone, the
%   same problem as the ANSYS 2D generalized-plane-strain model of the
%   inner leg (one coil sector) with a light, self-meshing linear FE model:
%     - every turn with its real rounded geometry (curved 8-node quads):
%       jacket (cable fillet r_SC, outer radius r_SC+JT), turn insulation,
%       corner fillers and inter-layer insulation; cable interiors,
%       ground insulation, case (flat plasma side, radial flanks at
%       +-pi/n_TF, nose = arc of radius Rk_/cos(pi/n_TF), WP cavity = convex
%       hull of the layers + ground insulation) and wedge insulation as
%       6-node triangles of a conforming Delaunay mesh;
%     - wedging: zero normal displacement on the flank lines, sliding free
%       (the wedge insulation is only normal-constrained, as in the FEM);
%     - frictional unilateral contact (Coulomb, mu = 0.2 as the FEM)
%       between WP and case and between every cable and its jacket;
%     - loads: Lorentz force J x B in every cable (WP_FIELD_AT_POINTS),
%       cool-down T_ref -> T_op (orthotropic insulation, filler, cable and
%       steel CTEs), and the vertical force T_bf per leg as generalized
%       plane strain.
%
%   Validated against two ANSYS runs (validation/validate_mech_surrogate_2026.m,
%   TF_FEM_benchmark_2026_findings.md): axial strain within 0.4%; linearized
%   jacket-wall stresses (Pm, Pm+Pb) on ~1400 sections each, mean ratio
%   0.99-1.01 with 4-7% scatter; case stress-classification lines within
%   -4..+2% (nose, vault, side walls); jacket peak in the cable fillet per
%   turn within -11..+8% (global maximum -6..-9%).
%
%   row - design point (DATA row or struct): Iop, n_layers, n_turns,
%         Cond_w, Cond_h, JT, type_cable, Ri_, Rk_
%   p   - machine parameters (READ_MACHINE_INPUT); the surrogate-specific
%         ones (category "Mechanical surrogate" of the input template)
%         fall back to the benchmark values if absent
%   opts (optional) - overrides: T_bf [N] (default from
%         COMPUTE_OPERATING_PARAMS), r_SC [m], n_cable, n_arc, n_thk,
%         h_fine, mu_case, mu_cable, mu_flank, interface_case
%         ('sliding'|'bonded'), interface_cable ('contact'|'bonded'), verbose
%
%   out - per-turn jacket results (out.turn: Pm, PmPb over all wall and
%         fillet sections, peak = Tresca in the fillet), per-layer maxima
%         (out.layer), case stress-classification lines (out.case_scl),
%         figure of merit vs the allowables (out.fom), nodal-averaged
%         stresses per material (out.nodal), element Tresca for plotting
%         (out.elem_SINT), mesh and solution.

if nargin < 3, opts = struct(); end
opts = set_default(opts, 'n_cable', get_or(p, 'surrogate_n_cable', 4));
opts = set_default(opts, 'n_arc', get_or(p, 'surrogate_n_arc', 8));
opts = set_default(opts, 'n_thk', get_or(p, 'surrogate_n_thk', 2));
opts = set_default(opts, 'h_fine', get_or(p, 'surrogate_h_fine', 3e-3));
opts = set_default(opts, 'h_mid', 6e-3);
opts = set_default(opts, 'h_coarse', 12e-3);
opts = set_default(opts, 'verbose', true);
opts = set_default(opts, 'interface_case', 'sliding');   % 'sliding' | 'bonded'
opts = set_default(opts, 'mu_case', get_or(p, 'mu_case', 0.2));    % WP / case friction
opts = set_default(opts, 'mu_cable', get_or(p, 'mu_cable', 0.2));  % cable / jacket friction
opts = set_default(opts, 'interface_cable', 'contact');              % 'contact' | 'bonded'
opts = set_default(opts, 'mu_flank', get_or(p, 'mu_flank', 0));    % flank: the wedge insulation is
% only normal-constrained on its outer face, so it slides with the case (as in the FEM)
t_start = tic;

geo = surr_geometry(row, p, opts);
mat = surr_materials(p, opts);
mesh = surr_mesh(geo, opts);
if strcmpi(opts.interface_case, 'sliding')
    mesh = split_cavity_interface(mesh, geo);
else
    mesh.cpair = zeros(0,2); mesh.cnorm = zeros(0,2);
end
if strcmpi(opts.interface_cable, 'contact')
    mesh = split_cable_interface(mesh, geo);
else
    mesh.kpair = zeros(0,2); mesh.knorm = zeros(0,2);
end
t_mesh = toc(t_start);

if ~isfield(opts, 'T_bf') || isempty(opts.T_bf)
    g = compute_operating_params(p);
    NI = sum(geo.n_turns)*row.Iop;
    opts.T_bf = 0.5*g.k_bf*p.n_TF*NI^2*(4e-7*pi)/(2*pi);
end

sol = surr_solve(mesh, mat, geo, row.Iop, opts.T_bf, opts);
t_solve = toc(t_start) - t_mesh;

out = surr_postprocess(mesh, sol, geo);
out.fom = surr_fom(out, p);
out.geo = geo; out.mesh = mesh; out.sol = sol; out.T_bf = opts.T_bf;
out.time = struct('mesh', t_mesh, 'solve', t_solve, 'total', toc(t_start));
if opts.verbose
    fprintf(['wp_mech_surrogate: %d nodes, %d quads, %d triangles, eps_z = %.3e | ' ...
        'mesh %.1f s, solve+post %.1f s\n'], size(mesh.xy,1), size(mesh.q8,1), ...
        size(mesh.t6,1), sol.eps_z, t_mesh, toc(t_start)-t_mesh);
end
end

%% ======================================================================
function s = set_default(s, name, val)
if ~isfield(s, name) || isempty(s.(name)), s.(name) = val; end
end

%% ======================================================================
function geo = surr_geometry(row, p, opts)
nl = row.n_layers;
geo.n_layers = nl;
geo.n_turns = row.n_turns(1:nl);
geo.Cond_w = row.Cond_w(1:nl);
geo.Cond_h = row.Cond_h(1:nl);
geo.JT = row.JT(1:nl);
geo.is_hts = strcmp(row.type_cable(1:nl), 'HTS');
geo.tins = p.turn_insulation_nominal*p.Increm;
geo.GIT = p.GoundIns;
geo.INS = p.INS_grades;
geo.Ri_ = row.Ri_;
geo.half = pi/p.n_TF;
geo.n_TF = p.n_TF;
geo.R_in = row.Rk_/cos(geo.half);            % nose arc radius (FEM WEDGE_INTR_INS)
geo.t_w = get_or(p, 'wedge_insulation', 1.5e-3);

Y2 = row.Ri_ - p.dr_plasma_side;              % top of cavity
Re = zeros(1, nl+1); Ri = zeros(1, nl);
Re(1) = Y2 - geo.GIT;
for k = 1:nl
    Ri(k) = Re(k) - geo.Cond_h(k);
    Re(k+1) = Ri(k) - geo.INS;
end
geo.Re = Re(1:nl); geo.Ri = Ri;
geo.half_w = geo.n_turns.*geo.Cond_w/2;

% Cavity = convex hull of the layers, widened by the ground insulation
hx = []; hy = [];
for k = 1:nl
    ytop = Re(k) + (k == 1)*geo.GIT;
    ybot = Ri(k) - (k == nl)*geo.GIT;
    w = geo.half_w(k) + geo.GIT;
    hx = [hx, -w, w, w, -w]; %#ok<AGROW>
    hy = [hy, ybot, ybot, ytop, ytop]; %#ok<AGROW>
end
kh = convhull(hx, hy);
geo.cav = [hx(kh(:)); hy(kh(:))]';            % closed polygon (first = last)

% Turn cells and cable rectangles
geo.cells = zeros(0, 6);                      % [x0 x1 y0 y1 layer col]
for k = 1:nl
    for j = 1:geo.n_turns(k)
        x0 = -geo.half_w(k) + (j-1)*geo.Cond_w(k);
        geo.cells(end+1,:) = [x0, x0+geo.Cond_w(k), Ri(k), Re(k), k, j];
    end
end
c = geo.cells; L = c(:,5);
off = geo.tins + geo.JT(L)';
geo.cable = [ (c(:,1)+c(:,2))/2, (c(:,3)+c(:,4))/2, ...
              c(:,2)-c(:,1)-2*off, c(:,4)-c(:,3)-2*off ];   % [xc yc w h]
if any(geo.cable(:,3) <= 0) || any(geo.cable(:,4) <= 0)
    error('wp_mech_surrogate:geometry', 'Cable size <= 0: check Cond_w/Cond_h vs JT and turn insulation.');
end
% cable corner fillet: r_SC = JT clamped to [r_SC_min, r_SC_max] (as
% size_cicc_cable / export_ansys_input); opts.r_SC overrides
if isfield(opts, 'r_SC') && ~isempty(opts.r_SC)
    rl = opts.r_SC.*ones(1, nl);
else
    rl = min(max(geo.JT, get_or(p, 'r_SC_min', 2e-3)), get_or(p, 'r_SC_max', 6e-3));
end
geo.r_SC = rl;
r = min(rl(L)', 0.49*min(geo.cable(:,3), geo.cable(:,4)));
geo.cable_rr = [geo.cable(:,1)-geo.cable(:,3)/2, geo.cable(:,1)+geo.cable(:,3)/2, ...
                geo.cable(:,2)-geo.cable(:,4)/2, geo.cable(:,2)+geo.cable(:,4)/2, r];
geo.cable_area = geo.cable(:,3).*geo.cable(:,4) - (4-pi)*r.^2;
end

function v = get_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end

%% ======================================================================
function mat = surr_materials(p, opts)
% Material table (SI). Defaults = the ANSYS benchmark material set
% (TFBM_model.txt): 1 LTS cable, 2 jacket, 3 insulation (glass-epoxy,
% orthotropic, local n = through-thickness), 4 case, 5 wedge insulation,
% 6 HTS cable. Alphas are secant 293 K -> T_op.
T_ref = get_or(opts, 'T_ref', get_or(p, 'T_ref', 293));
T_op  = get_or(opts, 'T_op',  get_or(p, 'T_op', 4.2));
dT = T_op - T_ref;
E_j = get_or(p, 'E_jckt', 205)*1e9; E_c = get_or(p, 'E_case', 205)*1e9;
a_steel = get_or(opts, 'alpha_steel', get_or(p, 'alpha_steel', 1.038e-5));
a_cbl   = get_or(opts, 'alpha_cable', get_or(p, 'alpha_cable', 5.54e-6));
ins.En = get_or(p, 'E_ins', 12)*1e9;
ins.Et = get_or(opts, 'E_ins_t', get_or(p, 'E_ins_t', 20))*1e9;
ins.G  = get_or(opts, 'G_ins', get_or(p, 'G_ins', 6))*1e9;
ins.nu_tn = get_or(p, 'nu_ins_tn', 0.33); ins.nu_tz = get_or(p, 'nu_ins_tz', 0.17);
nu_s = get_or(p, 'nu_steel', 0.29); nu_c = get_or(p, 'nu_cable', 0.3);
ins.an = get_or(opts, 'alpha_ins_n', get_or(p, 'alpha_ins_n', 2.422e-5));
ins.at = get_or(opts, 'alpha_ins_t', get_or(p, 'alpha_ins_t', 8.65e-6));

mat = cell(1, 7);
mat{1} = iso(get_or(p, 'E_cbl_LTS', 10)*1e9, nu_c, a_cbl, dT);
mat{2} = iso(E_j, nu_s, a_steel, dT);
mat{3} = ortho(ins, dT);
mat{4} = iso(E_c, nu_s, a_steel, dT);
mat{5} = ortho(ins, dT);
mat{6} = iso(get_or(p, 'E_cbl_HTS', 120)*1e9, nu_c, a_cbl, dT);
mat{7} = iso(get_or(p, 'E_filler', 7)*1e9, get_or(p, 'nu_filler', 0.3), ...
    get_or(opts, 'alpha_filler', get_or(p, 'alpha_filler', 1.73e-5)), dT);   % corner filler
end

function m = iso(E, nu, a, dT)
S = [1 -nu -nu 0; -nu 1 -nu 0; -nu -nu 1 0; 0 0 0 2*(1+nu)]/E;
m.S = S; m.eth = a*dT*[1;1;1;0]; m.ortho = false;
end

function m = ortho(ins, dT)
S = [1/ins.En, -ins.nu_tn/ins.Et, -ins.nu_tn/ins.Et, 0;
     -ins.nu_tn/ins.Et, 1/ins.Et, -ins.nu_tz/ins.Et, 0;
     -ins.nu_tn/ins.Et, -ins.nu_tz/ins.Et, 1/ins.Et, 0;
     0 0 0 1/ins.G];
m.S = S; m.eth = dT*[ins.an; ins.at; ins.at; 0]; m.ortho = true;
end

function [D, eth] = mat_D(m, phi)
% 4x4 stiffness for [exx eyy ezz gxy] in global axes; phi = angle of the
% local n axis from global x (only matters for orthotropic materials).
D0 = inv(m.S);
if ~m.ortho || phi == 0
    D = D0; eth = m.eth; return
end
c = cos(phi); s = sin(phi);
T = [c^2, s^2, 0, c*s; s^2, c^2, 0, -c*s; 0 0 1 0; -2*c*s, 2*c*s, 0, c^2-s^2];
D = T'*D0*T;
eth = T\m.eth;
end

%% ======================================================================
function mesh = surr_mesh(geo, opts)
% Every turn cell is meshed with curved 8-node quads that follow the real
% rounded geometry (as the FEM): jacket ring (cable fillet r_SC, outer
% radius r_SC+JT), turn-insulation ring (outer radius r_SC+JT+tins), corner
% fillers (fans between the ring and the square cell corner) and, below
% the cell, the inter-layer insulation strip. Adjacent layers are merged
% where their nodes coincide and tied (quadratic interpolation) where they
% do not. The cable interiors, the ground insulation, the case and the
% wedge insulation are 6-node triangles of one Delaunay mesh whose
% midside nodes on curved quad edges are snapped to the quads' ones.
tol = 1e-9;
nt = size(geo.cells,1);
lmax = get_or(opts, 'wall_element_max', 6e-3);   % max straight-wall element length
na = 2*ceil(get_or(opts, 'n_arc', 4)/2);           % elements per 90 deg fillet (even)
nth = get_or(opts, 'n_thk', 1);                    % elements through the jacket thickness

QX = cell(nt,1); QI = cell(nt,1);
for t = 1:nt
    c = geo.cable_rr(t,:);                         % [xa xb ya yb r]
    k = geo.cells(t,5); J = geo.JT(k); ti = geo.tins; r = c(5);
    Ro = r + J + ti;
    nsx = max(opts.n_cable, ceil((c(2)-c(1)-2*r)/lmax));
    nsy = max(opts.n_cable, ceil((c(4)-c(3)-2*r)/lmax));
    [P, Nn, wall, kk] = ring_stations(c(1), c(2), c(3), c(4), r, nsx, nsy, na);
    ne = numel(wall);
    X = zeros(0,16); info = zeros(0,6);
    dl = [(0:nth)*J/nth, J + ti];                  % jacket sub-layers, then insulation
    for e = 1:ne
        s0 = 2*e-1; sm = 2*e; s1 = mod(2*e, 2*ne) + 1;
        for L = 1:nth+1
            d0 = dl(L); d1 = dl(L+1); ring = 1 + (L > nth);
            Q = [P(s0,:)+d1*Nn(s0,:); P(s1,:)+d1*Nn(s1,:); P(s1,:)+d0*Nn(s1,:); P(s0,:)+d0*Nn(s0,:); ...
                 P(sm,:)+d1*Nn(sm,:); P(s1,:)+(d0+d1)/2*Nn(s1,:); P(sm,:)+d0*Nn(sm,:); P(s0,:)+(d0+d1)/2*Nn(s0,:)];
            X(end+1,:) = reshape(Q', 1, []); %#ok<AGROW>
            info(end+1,:) = [t, ring, wall(e), kk(e), atan2(Nn(sm,2), Nn(sm,1)), L]; %#ok<AGROW>
        end
        if wall(e) >= 5                           % corner filler fan outside the ring
            C = P(s0,:) - r*Nn(s0,:);              % fillet centre
            A0 = C + Ro*Nn(s0,:); A1 = C + Ro*Nn(s1,:); Am = C + Ro*Nn(sm,:);
            Q0 = C + Ro*Nn(s0,:)/max(abs(Nn(s0,:))); Q1 = C + Ro*Nn(s1,:)/max(abs(Nn(s1,:)));
            Q = [A1; A0; Q0; Q1; Am; (A0+Q0)/2; (Q0+Q1)/2; (Q1+A1)/2];
            X(end+1,:) = reshape(Q', 1, []); %#ok<AGROW>
            info(end+1,:) = [t, 3, wall(e), kk(e), 0, 0]; %#ok<AGROW>
        end
    end
    % inter-layer insulation strip under the cell (not under the last layer)
    if k < geo.n_layers
        yb = geo.cells(t,3); ys = yb - geo.INS;
        xs = unique(round([c(1)-J-ti, c(1)+r+(0:2*nsx)*(c(2)-c(1)-2*r)/(2*nsx), c(2)+J+ti, ...
            C_fan_x(c, Ro, na)]/tol))*tol;
        for i = 1:2:numel(xs)-2
            xa = xs(i); xm = xs(i+1); xb = xs(i+2);
            Q = [xa ys; xb ys; xb yb; xa yb; xm ys; xb (ys+yb)/2; xm yb; xa (ys+yb)/2];
            X(end+1,:) = reshape(Q', 1, []); %#ok<AGROW>
            info(end+1,:) = [t, 4, 0, 0, pi/2, 0]; %#ok<AGROW>
        end
    end
    QX{t} = X; QI{t} = info;
end
QX = cell2mat(QX); QI = cell2mat(QI);
% orientation (counter-clockwise corners)
for e = 1:size(QX,1)
    Q = reshape(QX(e,:), 2, [])';
    a = polyarea_signed(Q(1:4,:));
    if a < 0, Q = Q([1 4 3 2 8 7 6 5],:); QX(e,:) = reshape(Q', 1, []); end
end
qxy = reshape(QX', 2, [])';
[qn, ~, iq] = unique(round(qxy/tol), 'rows');
qn = qn*tol;
q8 = reshape(iq, 8, [])';
mesh.q8_turn = QI(:,1); mesh.q8_ring = QI(:,2); mesh.q8_wall = QI(:,3); mesh.q8_k = QI(:,4);
mesh.q8_lay = QI(:,6);                               % jacket sub-layer (1 = cable face)
matmap = [2 3 7 3];                                  % jacket, turn ins, filler, layer ins
mesh.q8_mat = matmap(QI(:,2))';
mesh.q8_phi = QI(:,5).*(QI(:,2) == 2 | QI(:,2) == 4);

% ---- ties between adjacent layers where nodes do not coincide -------
[tie, tied_edge] = layer_ties(geo, qn, q8, mesh.q8_ring, QI(:,1));
mesh.tie = tie;

% boundary edges of the quad region (used once), excluding tied interfaces
ed = [q8(:,[1 2 5]); q8(:,[2 3 6]); q8(:,[3 4 7]); q8(:,[4 1 8])];
ed = ed(ed(:,1) ~= ed(:,2), :);                      % collapsed fan sides
key = sort(ed(:,1:2), 2);
[~, ia, ic] = unique(key, 'rows');
cnt = accumarray(ic, 1);
bnd = ia(cnt == 1);
bedge = ed(bnd,:);
if ~isempty(tied_edge)
    bedge = bedge(~ismember(sort(bedge(:,1:2),2), sort(tied_edge,2), 'rows'), :);
end

% ---- Delaunay points ------------------------------------------------
hf = opts.h_fine;
Pq = qn(unique(bedge(:,1:2)),:);
P = [Pq; sample_polyline(geo.cav, hf); sample_outer(geo, hf); interior_points(geo, opts)];
[~, iu] = unique(round(P/tol), 'rows', 'stable');
P = P(iu,:);
[~, locA] = ismember(round(qn(bedge(:,1),:)/tol), round(P/tol), 'rows');
[~, locB] = ismember(round(qn(bedge(:,2),:)/tol), round(P/tol), 'rows');
if exist('delaunayTriangulation', 'file') || exist('delaunayTriangulation', 'class')
    DT = delaunayTriangulation(P(:,1), P(:,2), [locA, locB]);
    tri = DT.ConnectivityList; P = DT.Points;
else
    tri = delaunay(P(:,1), P(:,2));
end
cxy = (P(tri(:,1),:) + P(tri(:,2),:) + P(tri(:,3),:))/3;
area = 0.5*((P(tri(:,2),1)-P(tri(:,1),1)).*(P(tri(:,3),2)-P(tri(:,1),2)) - ...
            (P(tri(:,3),1)-P(tri(:,1),1)).*(P(tri(:,2),2)-P(tri(:,1),2)));
[cable_of, in_quads] = tri_region(geo, cxy);
keep = in_domain(geo, cxy) & ~in_quads & abs(area) > 1e-14;
tri = tri(keep,:); cxy = cxy(keep,:); area = area(keep); cable_of = cable_of(keep);
flip = area < 0; tri(flip,[2 3]) = tri(flip,[3 2]);

TE = sort([tri(:,[1 2]); tri(:,[2 3]); tri(:,[3 1])], 2);
missing = ~ismember(sort([locA locB], 2), TE, 'rows');
if any(missing)
    error('wp_mech_surrogate:mesh', '%d quad-boundary edges not conforming in the triangle mesh.', sum(missing));
end

% triangle materials
nt6 = size(tri,1);
tmat = zeros(nt6,1); tphi = zeros(nt6,1); tturn = zeros(nt6,1);
nR = [cos(geo.half), -sin(geo.half)]; nL = [-cos(geo.half), -sin(geo.half)];
dR = cxy*nR'; dL = cxy*nL';
in_wedge = dR > -geo.t_w | dL > -geo.t_w;
in_cav = inpolygon(cxy(:,1), cxy(:,2), geo.cav(:,1), geo.cav(:,2));
tmat(in_wedge) = 5; tphi(in_wedge & dR > -geo.t_w) = -geo.half; tphi(in_wedge & dL > -geo.t_w) = geo.half;
tmat(~in_wedge & ~in_cav) = 4;
ins = ~in_wedge & in_cav;
tmat(ins) = 3; tphi(ins) = ins_normal(geo, cxy(ins,:));
isc = cable_of > 0;
tturn(isc) = cable_of(isc);
tmat(isc) = 1 + 5*geo.is_hts(geo.cells(cable_of(isc),5))';

% ---- nodes: quads + triangle corners + triangle midsides ---------------
[~, pq] = ismember(round(P/tol), round(qn/tol), 'rows');
nid = pq; newp = find(pq == 0);
nid(newp) = size(qn,1) + (1:numel(newp))';
xy = [qn; P(newp,:)];
T3 = nid(tri);
bk = sort(bedge(:,1:2), 2);
E6 = sort([T3(:,[1 2]); T3(:,[2 3]); T3(:,[3 1])], 2);
[ue, ~, ie] = unique(E6, 'rows');
[isq, loc] = ismember(ue, bk, 'rows');
mid = zeros(size(ue,1),1);
mid(isq) = bedge(loc(isq), 3);
nnew = sum(~isq);
mid(~isq) = size(xy,1) + (1:nnew)';
xy = [xy; (xy(ue(~isq,1),:) + xy(ue(~isq,2),:))/2];
M = reshape(mid(ie), [], 3);
mesh.t6 = [T3, M];
mesh.xy = xy; mesh.q8 = q8;
mesh.t6_mat = tmat; mesh.t6_phi = tphi; mesh.t6_turn = tturn;

onR = abs(xy*nR') < 1e-7; onL = abs(xy*nL') < 1e-7;
mesh.flank_nodes = [find(onR); find(onL)];
mesh.flank_normal = [repmat(nR, sum(onR), 1); repmat(nL, sum(onL), 1)];

Aq = 0;
g = [-sqrt(3/5) 0 sqrt(3/5)]; w = [5 8 5]/9;
for e = 1:size(q8,1)
    Xe = xy(q8(e,:),:);
    for i = 1:3
        for j = 1:3
            [~, dJ] = elem_B(Xe, @q8_dN, [g(i) g(j)]); Aq = Aq + dJ*w(i)*w(j);
        end
    end
end
mesh.area_check = [sum(abs(area)) + Aq, domain_area(geo)];
end

function x = C_fan_x(c, Ro, na)
% x of the fan points on the bottom edge of the cell (both bottom corners)
r = c(5);
th = linspace(-pi/2, 0, na+1); th = th(th <= -pi/4 + 1e-12);
xr = (c(2)-r) + Ro*cos(th)./max(abs([cos(th); sin(th)]));
th = linspace(pi, 3*pi/2, na+1); th = th(th >= 5*pi/4 - 1e-12);
xl = (c(1)+r) + Ro*cos(th)./max(abs([cos(th); sin(th)]));
% fan midpoints on the edge
x = [xr, xl];
x = [x, (xr(1:end-1)+xr(2:end))/2, (xl(1:end-1)+xl(2:end))/2];
end

function a = polyarea_signed(Q)
a = 0.5*sum(Q(:,1).*Q([2:end 1],2) - Q([2:end 1],1).*Q(:,2));
end

function [tie, tied_edge] = layer_ties(geo, qn, q8, ring, qturn)
% Tie constraints on the interface between the insulation strip under
% layer k and the top of layer k+1 where the nodes do not coincide:
% each unmatched node is tied to the quadratic edge of the other side.
tie = struct('slave', {}, 'master', {}, 'w', {});
tied_edge = zeros(0,2);
tolx = 1e-9;
for k = 1:geo.n_layers-1
    y0 = geo.Re(k+1);
    upper = ismember(qturn, find(geo.cells(:,5) == k)) & ring == 4;     % strips of layer k
    lower = ismember(qturn, find(geo.cells(:,5) == k+1)) & ring ~= 4;   % cells of layer k+1
    EU = edges_on_line(qn, q8(upper,:), y0);
    EL = edges_on_line(qn, q8(lower,:), y0);
    if isempty(EU) || isempty(EL), continue, end
    xo = [max(min(qn(EU(:,1),1)), min(qn(EL(:,1),1))), min(max(qn(EU(:,2),1)), max(qn(EL(:,2),1)))];
    nU = unique(EU(:)); nL = unique(EL(:));
    if isempty(setxor(nU(in_range(qn(nU,1), xo)), nL(in_range(qn(nL,1), xo)))), continue, end
    for side = 1:2
        if side == 1, S = nL; Em = EU; else, S = nU; Em = EL; end
        S = S(in_range(qn(S,1), xo));
        S = setdiff(S, Em(:));
        for s = S'
            x = qn(s,1);
            e = find(qn(Em(:,1),1) <= x + tolx & qn(Em(:,2),1) >= x - tolx, 1);
            if isempty(e), continue, end
            xa = qn(Em(e,1),1); xb = qn(Em(e,2),1);
            xi = 2*(x - xa)/(xb - xa) - 1;
            N = [xi*(xi-1)/2, xi*(xi+1)/2, 1-xi^2];
            tie(end+1) = struct('slave', s, 'master', Em(e,:), 'w', N); %#ok<AGROW>
        end
    end
    inU = in_range(qn(EU(:,1),1), xo) & in_range(qn(EU(:,2),1), xo);
    inL = in_range(qn(EL(:,1),1), xo) & in_range(qn(EL(:,2),1), xo);
    tied_edge = [tied_edge; EU(inU,1:2); EL(inL,1:2)]; %#ok<AGROW>
end
end

function E = edges_on_line(qn, Q, y0)
% quad edges lying on y = y0, as [left right mid] node ids
cand = [Q(:,[1 2 5]); Q(:,[2 3 6]); Q(:,[3 4 7]); Q(:,[4 1 8])];
on = abs(qn(cand(:,1),2) - y0) < 1e-9 & abs(qn(cand(:,2),2) - y0) < 1e-9 & cand(:,1) ~= cand(:,2);
E = cand(on,:);
sw = qn(E(:,1),1) > qn(E(:,2),1);
E(sw,[1 2]) = E(sw,[2 1]);
E = unique(E, 'rows');
end

function tf = in_range(x, xo)
tf = x >= xo(1) - 1e-9 & x <= xo(2) + 1e-9;
end

function [P, Nn, wall, kk] = ring_stations(xa, xb, ya, yb, r, nsx, nsy, na)
% Stations (element ends and midpoints) around the rounded cable boundary,
% counter-clockwise from the start of the bottom straight wall.
% wall: 1 bottom, 2 right, 3 top, 4 left (straight), 5..8 fillets BR,TR,TL,BL
segs = { 'L', [xa+r ya], [xb-r ya], [0 -1], nsx, 1;
         'A', [xb-r ya+r], [-pi/2 0], [], na, 5;
         'L', [xb ya+r], [xb yb-r], [1 0], nsy, 2;
         'A', [xb-r yb-r], [0 pi/2], [], na, 6;
         'L', [xb-r yb], [xa+r yb], [0 1], nsx, 3;
         'A', [xa+r yb-r], [pi/2 pi], [], na, 7;
         'L', [xa yb-r], [xa ya+r], [-1 0], nsy, 4;
         'A', [xa+r ya+r], [pi 3*pi/2], [], na, 8 };
P = zeros(0,2); Nn = zeros(0,2); wall = zeros(0,1); kk = zeros(0,1);
for s = 1:size(segs,1)
    n = segs{s,5};
    tt = (0:2*n-1)'/(2*n);
    if segs{s,1} == 'L'
        A = segs{s,2}; B = segs{s,3};
        P = [P; A + tt*(B-A)]; Nn = [Nn; repmat(segs{s,4}, 2*n, 1)]; %#ok<AGROW>
    else
        C = segs{s,2}; th = segs{s,3}(1) + tt*(segs{s,3}(2)-segs{s,3}(1));
        P = [P; C + r*[cos(th) sin(th)]]; Nn = [Nn; [cos(th) sin(th)]]; %#ok<AGROW>
    end
    wall = [wall; segs{s,6}*ones(n,1)]; kk = [kk; (1:n)']; %#ok<AGROW>
end
end

function [cable_of, in_quads] = tri_region(geo, C)
% cable_of: turn whose cable contains the point (0 if none); in_quads:
% inside a cell (outside its cable) or inside an inter-layer strip
n = size(C,1);
cable_of = zeros(n,1); in_quads = false(n,1);
for t = 1:size(geo.cells,1)
    c = geo.cable_rr(t,:); cl = geo.cells(t,:);
    inc = C(:,1) > cl(1) & C(:,1) < cl(2) & C(:,2) > cl(3) & C(:,2) < cl(4);
    if cl(5) < geo.n_layers
        ins = C(:,1) > cl(1) & C(:,1) < cl(2) & C(:,2) > cl(3) - geo.INS & C(:,2) <= cl(3);
        in_quads = in_quads | ins;
    end
    if ~any(inc), continue, end
    idx = find(inc);
    in_cab = in_rrect(C(idx,:), c(1), c(2), c(3), c(4), c(5));
    cable_of(idx(in_cab)) = t;
    in_quads(idx(~in_cab)) = true;
end
end

function tf = in_rrect(C, xa, xb, ya, yb, r)
xc = (xa+xb)/2; yc = (ya+yb)/2; hw = (xb-xa)/2; hh = (yb-ya)/2;
dx = abs(C(:,1)-xc); dy = abs(C(:,2)-yc);
tf = dx < hw & dy < hh;
cr = dx > hw - r & dy > hh - r;
tf(cr) = (dx(cr)-(hw-r)).^2 + (dy(cr)-(hh-r)).^2 < r^2;
end

function mesh = split_cable_interface(mesh, geo)
% Duplicate the nodes on the cable face of every jacket for the cable
% triangles: cable and jacket then interact only through unilateral
% frictional contact (as the FEM), solved with the other contacts.
xy = mesh.xy;
ce = find(mesh.t6_turn > 0);
pairs = zeros(0,2); nrm = zeros(0,2);
for t = unique(mesh.t6_turn(ce))'
    el = ce(mesh.t6_turn(ce) == t);
    nodes = unique(mesh.t6(el,:));
    c = geo.cable_rr(t,:);
    [d, n] = rr_boundary(xy(nodes,:), c(1), c(2), c(3), c(4), c(5));
    on = nodes(abs(d) < 1e-8);
    non = n(abs(d) < 1e-8, :);
    newid = size(xy,1) + (1:numel(on))';
    xy = [xy; xy(on,:)]; %#ok<AGROW>
    T = mesh.t6(el,:);
    [tf, loc] = ismember(T, on);
    T(tf) = newid(loc(tf));
    mesh.t6(el,:) = T;
    pairs = [pairs; on, newid]; %#ok<AGROW>
    nrm = [nrm; -non]; %#ok<AGROW>                   % from the jacket into the cable
end
mesh.xy = xy; mesh.kpair = pairs; mesh.knorm = nrm;
end

function [d, n] = rr_boundary(C, xa, xb, ya, yb, r)
% signed distance to the rounded-rectangle boundary (>0 outside) and outward normal
xc = (xa+xb)/2; yc = (ya+yb)/2; hw = (xb-xa)/2 - r; hh = (yb-ya)/2 - r;
qx = abs(C(:,1)-xc) - hw; qy = abs(C(:,2)-yc) - hh;
sx = sign(C(:,1)-xc); sy = sign(C(:,2)-yc);
d = zeros(size(C,1),1); n = zeros(size(C,1),2);
cor = qx > 0 & qy > 0;
L = hypot(qx(cor), qy(cor));
d(cor) = L - r; n(cor,:) = [sx(cor).*qx(cor)./L, sy(cor).*qy(cor)./L];
side = ~cor & qx >= qy;
d(side) = qx(side) - r; n(side,:) = [sx(side), zeros(sum(side),1)];
tb = ~cor & ~side;
d(tb) = qy(tb) - r; n(tb,:) = [zeros(sum(tb),1), sy(tb)];
end

function mesh = split_cavity_interface(mesh, geo)
% Duplicate the nodes on the cavity boundary on the case side, so that the
% ground insulation and the case are connected only in the normal
% direction (frictionless unilateral contact, solved by active set).
xy = mesh.xy; V = geo.cav;
nseg = size(V,1)-1;
on = false(size(xy,1),1); segs = cell(size(xy,1),1);
for s = 1:nseg
    A = V(s,:); B = V(s+1,:); AB = B - A;
    t = ((xy(:,1)-A(1))*AB(1) + (xy(:,2)-A(2))*AB(2))/(AB*AB');
    d = hypot(xy(:,1)-A(1)-t*AB(1), xy(:,2)-A(2)-t*AB(2));
    hit = find(d < 1e-8 & t > -1e-9 & t < 1+1e-9);
    on(hit) = true;
    for h = hit', segs{h} = [segs{h}, s]; end
end
ids = find(on);
case_el = mesh.t6_mat == 4;
used_by_case = false(size(xy,1),1);
tc = mesh.t6(case_el,:); used_by_case(tc(:)) = true;
ids = ids(used_by_case(ids));
newid = zeros(size(xy,1),1);
newid(ids) = size(xy,1) + (1:numel(ids))';
mesh.xy = [xy; xy(ids,:)];
tc(ismember(tc, ids)) = newid(tc(ismember(tc, ids)));
mesh.t6(case_el,:) = tc;
% constrained pairs (insulation node, case node) with outward normal
pairs = zeros(0,2); nrm = zeros(0,2);
cw = polyarea(V(:,1), V(:,2)) > 0;          % orientation sign
for a = ids'
    for s = segs{a}
        AB = V(s+1,:) - V(s,:); n = [AB(2), -AB(1)]/norm(AB);
        if ~cw, n = -n; end
        % make n point from the cavity into the case
        mid = (V(s,:)+V(s+1,:))/2;
        if inpolygon(mid(1)+1e-6*n(1), mid(2)+1e-6*n(2), V(:,1), V(:,2)), n = -n; end
        pairs(end+1,:) = [a, newid(a)]; %#ok<AGROW>
        nrm(end+1,:) = n; %#ok<AGROW>
    end
end
mesh.cpair = pairs; mesh.cnorm = nrm;
end

function Q = sample_polyline(V, h)
Q = zeros(0,2);
for s = 1:size(V,1)-1
    L = norm(V(s+1,:)-V(s,:));
    n = max(1, ceil(L/h));
    tt = (0:n-1)'/n;
    Q = [Q; V(s,:) + tt*(V(s+1,:)-V(s,:))]; %#ok<AGROW>
end
end

function Q = sample_outer(geo, h)
% Outer boundary (flank lines, top, nose arc) + the offset lines between
% the case and the wedge insulation, sampled at matching positions.
hf = geo.half; R = geo.R_in; tw = geo.t_w;
top = geo.Ri_;
r_top = top/cos(hf);                         % flank reaches y = Ri_
n = max(2, ceil((r_top - R)/h));
r = linspace(R, r_top, n+1)';
QF = zeros(0,2); side = zeros(0,1); Q = zeros(0,2);
for sgn = [1 -1]
    d = [sgn*sin(hf), cos(hf)];              % flank direction
    nrm = [sgn*cos(hf), -sin(hf)];           % outward normal
    F = r*d;                                 % flank line (constrained)
    % case/wedge-insulation interface: parallel line at distance tw inside,
    % from its intersection with the nose arc to its intersection with y = Ri_
    s_arc = sqrt(R^2 - tw^2); s_top = (top - tw*sin(hf))/cos(hf);
    O = linspace(s_arc, s_top, n+1)'*d - tw*nrm;
    QF = [QF; F]; side = [side; sgn*ones(size(F,1),1)]; %#ok<AGROW>
    Q = [Q; F; O]; %#ok<AGROW>
end
% top segment
xt = top*tan(hf);
nt = max(2, ceil(2*xt/h));
Q = [Q; [linspace(-xt, xt, nt+1)', top*ones(nt+1,1)]];
% nose arc
na = max(4, ceil(2*hf*R/h));
a = linspace(-hf, hf, na+1)';
Q = [Q; [R*sin(a), R*cos(a)]];
end

function P = interior_points(geo, opts)
% Graded interior points: fine near the cavity and the outer boundary.
hs = [opts.h_fine, opts.h_mid, opts.h_coarse];
dlim = [3*opts.h_fine, 4*opts.h_mid, Inf];
xmax = geo.Ri_*tan(geo.half);
P = zeros(0,2);
for lev = 1:3
    h = hs(lev);
    [X, Y] = meshgrid(-xmax:h:xmax, geo.R_in*cos(geo.half):h:geo.Ri_);
    C = [X(:), Y(:)];
    C = C(in_domain(geo, C), :);
    d_cav = dist_polyline(C, geo.cav);
    d_out = dist_outer(geo, C);
    in_cav = inpolygon(C(:,1), C(:,2), geo.cav(:,1), geo.cav(:,2));
    d = min(d_cav, d_out);
    lo = [0, dlim(1:end-1)];
    sel = d >= lo(lev) & d < dlim(lev) & d > 0.45*h & ~in_cav;
    P = [P; C(sel,:)]; %#ok<AGROW>
    if lev == 1
        % inside the cavity, outside the turns (ground-insulation fills)
        dc = dist_cells(geo, C);           % fills of stepped WPs only
        sel2 = in_cav & dc > 1.5*h & d_cav > 1.5*h;   % only in wide fills, not in thin bands
        P = [P; C(sel2,:)]; %#ok<AGROW>
    end
end
end

function tf = in_domain(geo, C)
nR = [cos(geo.half), -sin(geo.half)]; nL = [-cos(geo.half), -sin(geo.half)];
tf = C(:,2) <= geo.Ri_ + 1e-12 & sqrt(sum(C.^2,2)) >= geo.R_in - 1e-12 & ...
     C*nR' <= 1e-12 & C*nL' <= 1e-12;
end

function d = dist_cells(geo, C)
d = inf(size(C,1),1);
for t = 1:size(geo.cells,1)
    c = geo.cells(t,:);
    dx = max([c(1)-C(:,1), zeros(size(C,1),1), C(:,1)-c(2)], [], 2);
    dy = max([c(3)-C(:,2), zeros(size(C,1),1), C(:,2)-c(4)], [], 2);
    d = min(d, hypot(dx, dy));
end
end

function d = dist_polyline(C, V)
d = inf(size(C,1),1);
for s = 1:size(V,1)-1
    A = V(s,:); B = V(s+1,:); AB = B - A;
    t = ((C(:,1)-A(1))*AB(1) + (C(:,2)-A(2))*AB(2))/(AB*AB');
    t = min(max(t,0),1);
    d = min(d, hypot(C(:,1)-A(1)-t*AB(1), C(:,2)-A(2)-t*AB(2)));
end
end

function d = dist_outer(geo, C)
nR = [cos(geo.half), -sin(geo.half)]; nL = [-cos(geo.half), -sin(geo.half)];
d = min([abs(C*nR' + geo.t_w), abs(C*nL' + geo.t_w), abs(C*nR'), abs(C*nL'), geo.Ri_ - C(:,2), ...
         sqrt(sum(C.^2,2)) - geo.R_in], [], 2);
end

function A = domain_area(geo)
% sector between the nose arc and y = Ri_ (numerical, fine polygon)
a = linspace(-geo.half, geo.half, 2001)';
arc = geo.R_in*[sin(a), cos(a)];
xt = geo.Ri_*tan(geo.half);
V = [arc; xt, geo.Ri_; -xt, geo.Ri_];
A = polyarea(V(:,1), V(:,2));
end

function phi = ins_normal(geo, C)
% Through-thickness direction of the ground / inter-layer insulation:
% vertical (pi/2) in the strips between layers and above/below the WP,
% horizontal (0) beside the layers.
phi = pi/2*ones(size(C,1),1);
for k = 1:geo.n_layers
    ytop = geo.Re(k) + (k == 1)*geo.GIT + (k > 1)*geo.INS/2;
    ybot = geo.Ri(k) - (k == geo.n_layers)*geo.GIT - (k < geo.n_layers)*geo.INS/2;
    beside = C(:,2) < ytop & C(:,2) > ybot & abs(C(:,1)) > geo.half_w(k);
    phi(beside) = 0;
end
end

%% ======================================================================
function sol = surr_solve(mesh, mat, geo, Iop, T_bf, opts)
xy = mesh.xy; N = size(xy,1);
ndof = 2*N + 1; iz = ndof;

% Gauss rules
g = sqrt(3/5); w1 = 5/9; w2 = 8/9;
[GX, GY] = meshgrid([-g 0 g], [-g 0 g]); [WX, WY] = meshgrid([w1 w2 w1], [w1 w2 w1]);
gq = [GX(:), GY(:)]; wq = WX(:).*WY(:);
gt = [1/6 1/6; 2/3 1/6; 1/6 2/3]; wt = [1 1 1]'/6;

nq = size(mesh.q8,1); nt = size(mesh.t6,1);
nnz_est = nq*17^2 + nt*13^2;
I = zeros(nnz_est,1); Jc = I; V = I; ptr = 0;
F = zeros(ndof,1);

% Lorentz body force at the cable Gauss points (triangles; one field call)
ec = find(mesh.t6_turn > 0);
Pg = zeros(numel(ec)*3, 2);
Nt = t6_N(gt);
for a = 1:numel(ec)
    Pg(3*a-2:3*a,:) = Nt*xy(mesh.t6(ec(a),:),:);
end
[Bx, By] = wp_field_at_points(Pg(:,1), Pg(:,2), geo.cable(:,1), geo.cable(:,2), ...
    geo.cable(:,3), geo.cable(:,4), Iop, geo.n_TF);
Jd = Iop./geo.cable_area;                  % current density per turn (rounded cable)
fbx = zeros(nt, 3); fby = zeros(nt, 3);
for a = 1:numel(ec)
    Jt = Jd(mesh.t6_turn(ec(a)));
    fbx(ec(a),:) = -Jt*By(3*a-2:3*a)';      % f = J x B, J along +z
    fby(ec(a),:) =  Jt*Bx(3*a-2:3*a)';
end

% quads (jacket and insulation rings): identical turns share element
% matrices, cached by the element shape relative to its first node
cache = containers.Map();
sol.qD = cell(nq,1); sol.qeth = cell(nq,1);
for e = 1:nq
    nodes = mesh.q8(e,:); X = xy(nodes,:);
    [D, eth] = mat_D(mat{mesh.q8_mat(e)}, mesh.q8_phi(e));
    sol.qD{e} = D; sol.qeth{e} = eth;
    rel = round((X - X(1,:))'/1e-9);
    key = [sprintf('%d_', rel(:)), sprintf('%d_%.6f', mesh.q8_mat(e), mesh.q8_phi(e))];
    if isKey(cache, key)
        kk = cache(key); Ke = kk{1}; fe = kk{2};
    else
        [Ke, fe] = elem_K(X, D, eth, @q8_N, @q8_dN, gq, wq);
        cache(key) = {Ke, fe};
    end
    dofs = [reshape([2*nodes-1; 2*nodes], 1, []), iz];
    [ii, jj] = ndgrid(dofs, dofs);
    m = numel(Ke);
    I(ptr+1:ptr+m) = ii(:); Jc(ptr+1:ptr+m) = jj(:); V(ptr+1:ptr+m) = Ke(:); ptr = ptr + m;
    F(dofs) = F(dofs) + fe;
end
sol.tD = cell(nt,1); sol.teth = cell(nt,1);
for e = 1:nt
    nodes = mesh.t6(e,:); X = xy(nodes,:);
    [D, eth] = mat_D(mat{mesh.t6_mat(e)}, mesh.t6_phi(e));
    sol.tD{e} = D; sol.teth{e} = eth;
    [Ke, fe] = elem_K(X, D, eth, @t6_N, @t6_dN, gt, wt);
    if mesh.t6_turn(e) > 0
        fe = fe + elem_body(X, @t6_N, @t6_dN, gt, wt, fbx(e,:)', fby(e,:)');
    end
    dofs = [reshape([2*nodes-1; 2*nodes], 1, []), iz];
    [ii, jj] = ndgrid(dofs, dofs);
    m = numel(Ke);
    I(ptr+1:ptr+m) = ii(:); Jc(ptr+1:ptr+m) = jj(:); V(ptr+1:ptr+m) = Ke(:); ptr = ptr + m;
    F(dofs) = F(dofs) + fe;
end
K = sparse(I(1:ptr), Jc(1:ptr), V(1:ptr), ndof, ndof);
K = (K + K')/2;                 % exact symmetry -> Cholesky in backslash

% Axial force (generalized plane strain)
F(iz) = F(iz) + T_bf;

% Layer-to-layer ties where the nodes do not coincide (penalty)
kt = 1e3*max(diag(K));
if isfield(mesh, 'tie') && ~isempty(mesh.tie)
    Ii = []; Jj = []; Vv = [];
    for q = 1:numel(mesh.tie)
        nd = [mesh.tie(q).slave, mesh.tie(q).master];
        cw = [1, -mesh.tie(q).w];
        for d = 0:1                                  % x then y
            dofs = 2*nd - 1 + d;
            [ii, jj] = ndgrid(dofs, dofs);
            Ii = [Ii; ii(:)]; Jj = [Jj; jj(:)]; Vv = [Vv; kt*reshape(cw'*cw, [], 1)]; %#ok<AGROW>
        end
    end
    K = K + sparse(Ii, Jj, Vv, ndof, ndof);
end

% Contacts (penalty + active set, Coulomb friction):
%  - flanks: bilateral zero normal displacement (as the FEM), friction mu_flank
%  - cavity: WP ground insulation / case, unilateral, friction mu_case
kp = 1e3*max(diag(K));
[U, cinfo] = solve_contacts(K, F, mesh, kp, opts);
sol.contact = cinfo;
sol.U = U(1:end-1); sol.eps_z = U(end);
sol.Fz_check = T_bf;
end

function [U, info] = solve_contacts(K, F, mesh, kp, opts)
% Penalty contact with Coulomb friction (stick/slip return mapping on the
% total slip, as a one-step implicit FEM contact solution):
%  - normal: closed pairs get a penalty spring kp on the gap; a pair opens
%    when its gap becomes positive (tension) and closes when it penetrates;
%    flank pairs are bilateral (always closed);
%  - tangential: stick = penalty spring kp; a pair slips when the stick
%    force exceeds mu*N and then carries the sliding force mu*N; it sticks
%    again if the slip reverses. Stops when fewer than 0.5% of the pairs
%    change state and the contact forces are stable within 1%.
ndof = size(K,1);
fa = mesh.flank_nodes; fnv = mesh.flank_normal;
ca = mesh.cpair(:,1); cb = mesh.cpair(:,2); cnv = mesh.cnorm;
ka = mesh.kpair(:,1); kb = mesh.kpair(:,2); knv = mesh.knorm;
a = [fa; ca; ka]; b = [zeros(size(fa)); cb; kb]; n = [fnv; cnv; knv];
bil = [true(size(fa)); false(size(ca)); false(size(ka))];
mu = [opts.mu_flank*ones(size(fa)); opts.mu_case*ones(size(ca)); opts.mu_cable*ones(size(ka))];
info.group = [ones(size(fa)); 2*ones(size(ca)); 3*ones(size(ka))];
t = [-n(:,2), n(:,1)];
np_ = numel(a);
closed = true(np_,1);
slip = zeros(np_,1);             % 0 stick, +-1 slipping along +-t
N = zeros(np_,1);
info.hist = [];
maxit = get_or(opts, 'contact_maxit', 40);
done = false;
for it = 1:maxit
    [Kn_i, Kn_j, Kn_v] = pair_matrix(a(closed), b(closed), n(closed,:), kp);
    st = closed & slip == 0 & mu > 0;
    [Kt_i, Kt_j, Kt_v] = pair_matrix(a(st), b(st), t(st,:), kp);
    Kc = K + sparse([Kn_i; Kt_i], [Kn_j; Kt_j], [Kn_v; Kt_v], ndof, ndof);
    Fc = F;
    sl = closed & slip ~= 0;
    fs = -mu(sl).*max(N(sl),0).*slip(sl);          % sliding force on a along t
    Fc = Fc + accumarray([2*a(sl)-1; 2*a(sl)], [fs.*t(sl,1); fs.*t(sl,2)], [ndof 1]);
    g = sl & b > 0;
    fg = -mu(g).*max(N(g),0).*slip(g);
    Fc = Fc - accumarray([2*b(g)-1; 2*b(g)], [fg.*t(g,1); fg.*t(g,2)], [ndof 1]);
    U = spd_solve(Kc, Fc);
    ua = [U(2*a-1), U(2*a)];
    ub = zeros(np_,2); gb = b > 0; ub(gb,:) = [U(2*b(gb)-1), U(2*b(gb))];
    gn = sum((ub - ua).*n, 2);                     % >0: open
    gt = sum((ua - ub).*t, 2);
    N_new = -kp*gn; N_new(~closed) = 0;
    closed_new = bil | (closed & gn <= 0) | (~closed & gn < 0);
    Tst = -kp*gt;
    slip_new = slip;
    q = closed_new & mu > 0 & slip == 0 & abs(Tst) > mu.*max(N_new,0);
    slip_new(q) = sign(gt(q));
    q = closed_new & slip ~= 0 & gt.*slip < 0;
    slip_new(q) = 0;
    slip_new(~closed_new) = 0;
    dN = max(abs(N_new - N))/max([abs(N_new); 1]);
    nchg = sum(closed_new ~= closed) + sum(slip_new ~= slip);
    info.hist(end+1,:) = [it, sum(closed_new), sum(slip_new ~= 0), nchg, dN];
    done = nchg <= max(2, 5e-3*np_) && dN < 1e-2;
    closed = closed_new; slip = slip_new; N = N_new;
    if done, break, end
end
info.iter = it; info.closed = closed; info.slip = slip; info.N = N;
info.n_cavity_closed = sum(closed(info.group == 2)); info.n_cavity = sum(info.group == 2);
info.n_cable_closed = sum(closed(info.group == 3)); info.n_cable = sum(info.group == 3);
info.converged = done;
if ~done
    warning('wp_mech_surrogate:contact', ...
        'contact iterations did not settle in %d iterations (%d state changes, dN = %.1e)', maxit, nchg, dN);
end
end

function [I, J, V] = pair_matrix_w(a, b, v, k)
% as pair_matrix, with a stiffness per pair
m = numel(a);
if m == 0, I = zeros(0,1); J = I; V = I; return, end
g = b > 0;
dofs = [2*a-1, 2*a, 2*max(b,1)-1, 2*max(b,1)];
coef = [v(:,1), v(:,2), -v(:,1).*g, -v(:,2).*g];
I = zeros(16*m,1); J = I; V = I; q = 0;
for r = 1:4
    for c = 1:4
        I(q+1:q+m) = dofs(:,r); J(q+1:q+m) = dofs(:,c);
        V(q+1:q+m) = k.*coef(:,r).*coef(:,c); q = q + m;
    end
end
end

function U = spd_solve(K, F)
% Sparse Cholesky with fill-reducing ordering (backslash may miss the
% symmetry after penalty assembly and fall back to a much slower LU)
K = (K + K')/2;
[R, flag, q] = chol(K, 'vector');
if flag == 0
    U = zeros(size(F));
    U(q) = R\(R'\F(q));
else
    U = K\F;
end
end

function [I, J, V] = pair_matrix(a, b, v, k)
% penalty on (u_a - u_b).v for every pair; b = 0 -> ground
m = numel(a);
if m == 0, I = zeros(0,1); J = I; V = I; return, end
g = b > 0;
dofs = [2*a-1, 2*a, 2*max(b,1)-1, 2*max(b,1)];
coef = [v(:,1), v(:,2), -v(:,1).*g, -v(:,2).*g];
I = zeros(16*m,1); J = I; V = I; q = 0;
for r = 1:4
    for c = 1:4
        I(q+1:q+m) = dofs(:,r); J(q+1:q+m) = dofs(:,c);
        V(q+1:q+m) = k*coef(:,r).*coef(:,c); q = q + m;
    end
end
end

function [Ke, fth] = elem_K(X, D, eth, Nf, dNf, gp, w)
n = size(X,1);
Ke = zeros(2*n+1); fth = zeros(2*n+1,1);
for q = 1:size(gp,1)
    [B, detJ] = elem_B(X, dNf, gp(q,:));
    Ke = Ke + (B'*D*B)*detJ*w(q);
    fth = fth + (B'*D*eth)*detJ*w(q);
end
end

function fe = elem_body(X, Nf, dNf, gp, w, fx, fy)
n = size(X,1);
fe = zeros(2*n+1,1);
for q = 1:size(gp,1)
    [~, detJ] = elem_B(X, dNf, gp(q,:));
    Nq = Nf(gp(q,:));
    fe(1:2:2*n) = fe(1:2:2*n) + Nq'*fx(q)*detJ*w(q);
    fe(2:2:2*n) = fe(2:2:2*n) + Nq'*fy(q)*detJ*w(q);
end
end

function [B, detJ] = elem_B(X, dNf, xi)
dN = dNf(xi);                   % 2 x n (d/dxi; d/deta)
Jm = dN*X;                      % 2x2
detJ = det(Jm);
dNx = Jm\dN;                    % 2 x n (d/dx; d/dy)
n = size(X,1);
B = zeros(4, 2*n+1);
B(1,1:2:2*n) = dNx(1,:);
B(2,2:2:2*n) = dNx(2,:);
B(3,end) = 1;                   % eps_z (generalized plane strain)
B(4,1:2:2*n) = dNx(2,:);
B(4,2:2:2*n) = dNx(1,:);
end

function N = q8_N(xi)
x = xi(:,1); e = xi(:,2);
N = [-(1-x).*(1-e).*(1+x+e)/4, -(1+x).*(1-e).*(1-x+e)/4, -(1+x).*(1+e).*(1-x-e)/4, ...
     -(1-x).*(1+e).*(1+x-e)/4, (1-x.^2).*(1-e)/2, (1+x).*(1-e.^2)/2, ...
     (1-x.^2).*(1+e)/2, (1-x).*(1-e.^2)/2];
end

function dN = q8_dN(xi)
x = xi(1); e = xi(2);
dN = [ (1-e)*(2*x+e)/4, (1-e)*(2*x-e)/4, (1+e)*(2*x+e)/4, (1+e)*(2*x-e)/4, ...
       -x*(1-e), (1-e^2)/2, -x*(1+e), -(1-e^2)/2;
       (1-x)*(x+2*e)/4, (1+x)*(2*e-x)/4, (1+x)*(x+2*e)/4, (1-x)*(2*e-x)/4, ...
       -(1-x^2)/2, -(1+x)*e, (1-x^2)/2, -(1-x)*e ];
end

function N = t6_N(xi)
L2 = xi(:,1); L3 = xi(:,2); L1 = 1 - L2 - L3;
N = [L1.*(2*L1-1), L2.*(2*L2-1), L3.*(2*L3-1), 4*L1.*L2, 4*L2.*L3, 4*L3.*L1];
end

function dN = t6_dN(xi)
L2 = xi(1); L3 = xi(2); L1 = 1 - L2 - L3;
dN = [ -(4*L1-1), 4*L2-1, 0, 4*(L1-L2), 4*L3, -4*L3;
       -(4*L1-1), 0, 4*L3-1, -4*L2, 4*L2, 4*(L1-L3) ];
end

%% ======================================================================
function out = surr_postprocess(mesh, sol, geo)
xy = mesh.xy; U = sol.U; ez = sol.eps_z;

% Nodal-averaged stresses per material (like the FEM's material-restricted
% nodal averaging), evaluated directly at the element nodes.
q8nat = 0.999*[-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0];   % inset: collapsed fan quads
t6nat = [0 0; 1 0; 0 1; .5 0; .5 .5; 0 .5];
nmat = 7; N = size(xy,1);
acc = zeros(N, 4, nmat); cnt = zeros(N, nmat);
for e = 1:size(mesh.q8,1)
    nodes = mesh.q8(e,:); m = mesh.q8_mat(e);
    S = elem_stress_at(xy(nodes,:), U, nodes, ez, sol.qD{e}, sol.qeth{e}, @q8_dN, q8nat);
    acc(nodes,:,m) = acc(nodes,:,m) + S; cnt(nodes,m) = cnt(nodes,m) + 1;
end
for e = 1:size(mesh.t6,1)
    nodes = mesh.t6(e,:); m = mesh.t6_mat(e);
    S = elem_stress_at(xy(nodes,:), U, nodes, ez, sol.tD{e}, sol.teth{e}, @t6_dN, t6nat);
    acc(nodes,:,m) = acc(nodes,:,m) + S; cnt(nodes,m) = cnt(nodes,m) + 1;
end
names = {'cable_LTS','jacket','insulation','case','wedge_ins','cable_HTS','filler'};
for m = 1:nmat
    has = cnt(:,m) > 0;
    Sm = acc(has,:,m)./cnt(has,m);
    tr = tresca(Sm);
    out.nodal.(names{m}) = struct('node', find(has), 'S', Sm, 'SINT', tr);
    if any(has)
        [mx, i] = max(tr); nid = find(has);
        out.([names{m} '_SINT_max']) = mx;
        out.([names{m} '_SINT_max_xy']) = xy(nid(i),:);
    else
        out.([names{m} '_SINT_max']) = NaN;
        out.([names{m} '_SINT_max_xy']) = [NaN NaN];
    end
end

% Linearized jacket-wall stresses (Pm, Pm+Pb, Tresca) per turn, on
% through-thickness sections of the jacket ring: every element end along
% the four straight walls, and the middle of every fillet (45 deg).
% Peak: nodal-averaged Tresca on the cable face of the jacket (fillets).
nturn = size(geo.cells,1);
wallnames = {'bottom','right','top','left'};
T = struct('layer', num2cell(geo.cells(:,5)), 'col', num2cell(geo.cells(:,6)));
jq = find(mesh.q8_ring == 1);
% inner-face nodal averaging for the peak
accp = zeros(size(xy,1), 4); cntp = zeros(size(xy,1), 1);
for e = jq(mesh.q8_lay(jq) == 1)'
    nodes = mesh.q8(e,:);
    S = elem_stress_at(xy(nodes,:), U, nodes, ez, sol.qD{e}, sol.qeth{e}, @q8_dN, [1 1; -1 1; 0 1]);
    id = nodes([3 4 7]);
    accp(id,:) = accp(id,:) + S; cntp(id) = cntp(id) + 1;
end
has = cntp > 0;
peakS = zeros(size(xy,1),1); peakS(has) = tresca(accp(has,:)./cntp(has));
sec = zeros(0, 10);    % [turn wall station kind(0 end,1 interior,2 fillet) outer_xy inner_xy Pm PmPb]
for t = 1:nturn
    et = jq(mesh.q8_turn(jq) == t);
    res = struct();
    for w = 1:4
        ew = et(mesh.q8_wall(et) == w);
        ks = unique(mesh.q8_k(ew))';
        ns = numel(ks) + 1;
        Pm = zeros(1, ns); PmPb = zeros(1, ns); Sm = zeros(ns, 4);
        for s = 1:ns
            if s <= numel(ks), kq = ks(s); a = -1; else, kq = ks(end); a = 1; end
            e = stack_of(mesh, ew, kq);
            [Pm(s), PmPb(s), Sm(s,:), en] = lin_section(xy, U, ez, mesh, sol, e, a);
            sec(end+1,:) = [t, w, s, (s > 1 && s < ns), en(1,:), en(2,:), Pm(s), PmPb(s)]; %#ok<AGROW>
        end
        res.(wallnames{w}) = struct('Pm', Pm, 'PmPb', PmPb, 'S_membrane', Sm);
    end
    fil = zeros(4, 3);                           % [Pm PmPb 0] at the fillet middles
    for w = 5:8
        ef = et(mesh.q8_wall(et) == w);
        ks = unique(mesh.q8_k(ef))'; n = numel(ks);
        if mod(n, 2) == 0, kq = ks(n/2); a = 1; else, kq = ks((n+1)/2); a = 0; end
        e = stack_of(mesh, ef, kq);
        [fil(w-4,1), fil(w-4,2), ~, en] = lin_section(xy, U, ez, mesh, sol, e, a);
        sec(end+1,:) = [t, w, 0, 2, en(1,:), en(2,:), fil(w-4,1), fil(w-4,2)]; %#ok<AGROW>
    end
    str_in_Pm = []; str_in_Pb = []; str_all_Pm = []; str_all_Pb = [];
    for w = 1:4
        r = res.(wallnames{w}); in = 2:numel(r.Pm)-1;
        str_in_Pm = [str_in_Pm, r.Pm(in)]; str_in_Pb = [str_in_Pb, r.PmPb(in)]; %#ok<AGROW>
        str_all_Pm = [str_all_Pm, r.Pm]; str_all_Pb = [str_all_Pb, r.PmPb]; %#ok<AGROW>
    end
    res.fillet = struct('Pm', fil(:,1)', 'PmPb', fil(:,2)');
    T(t).walls = res;
    T(t).Pm_straight = max(str_in_Pm);           % straight walls, away from the fillets
    T(t).PmPb_straight = max(str_in_Pb);
    T(t).Pm = max([str_all_Pm, fil(:,1)']);      % all sections incl. fillet start/middle
    T(t).PmPb = max([str_all_Pb, fil(:,2)']);
    inner_nodes = unique(mesh.q8(et(mesh.q8_lay(et) == 1), [3 4 7]));
    [T(t).peak, ip] = max(peakS(inner_nodes));
    T(t).peak_xy = xy(inner_nodes(ip),:);
end
out.turn = T;
out.sections = sec;
L = [T.layer];
for k = 1:geo.n_layers
    out.layer.Pm(k) = max([T(L==k).Pm]);
    out.layer.PmPb(k) = max([T(L==k).PmPb]);
    out.layer.Pm_straight(k) = max([T(L==k).Pm_straight]);
    out.layer.PmPb_straight(k) = max([T(L==k).PmPb_straight]);
    out.layer.peak(k) = max([T(L==k).peak]);
end
out.case_scl = case_scl(mesh, sol, geo);

% element Tresca (at the element centre) for plotting
out.elem_SINT.q8 = zeros(size(mesh.q8,1),1);
for e = 1:size(mesh.q8,1)
    nodes = mesh.q8(e,:);
    out.elem_SINT.q8(e) = tresca(elem_stress_at(xy(nodes,:), U, nodes, ez, sol.qD{e}, sol.qeth{e}, @q8_dN, [0 0]));
end
out.elem_SINT.t6 = zeros(size(mesh.t6,1),1);
for e = 1:size(mesh.t6,1)
    nodes = mesh.t6(e,:);
    out.elem_SINT.t6(e) = tresca(elem_stress_at(xy(nodes,:), U, nodes, ez, sol.tD{e}, sol.teth{e}, @t6_dN, [1/3 1/3]));
end
out.eps_z = ez;
end

function e = stack_of(mesh, list, kq)
% jacket elements of one station, ordered from the outer face to the cable face
e = list(mesh.q8_k(list) == kq);
[~, o] = sort(mesh.q8_lay(e), 'descend'); e = e(o);
end

function [Pm, PmPb, m, ends] = lin_section(xy, U, ez, mesh, sol, elist, a)
% Linearization across the jacket thickness (all sub-layers, eta = -1
% outer .. +1 inner of each) at xi = a: membrane = mean, bending = surface
% value of the linear part (Tresca of both).
Sall = zeros(0,4); Pall = zeros(0,2);
nat = [a*ones(5,1), linspace(-1,1,5)'];
for e = elist(:)'
    nodes = mesh.q8(e,:);
    S = elem_stress_at(xy(nodes,:), U, nodes, ez, sol.qD{e}, sol.qeth{e}, @q8_dN, nat);
    Sall = [Sall; S]; Pall = [Pall; q8_N(nat)*xy(nodes,:)]; %#ok<AGROW>
end
d = [0; cumsum(sqrt(sum(diff(Pall).^2, 2)))];
t = d(end); z = d - t/2;
m = trapz(d, Sall)/t;
b = 6/t^2*trapz(d, Sall.*z);
Pm = tresca(m);
PmPb = max(tresca(m + b), tresca(m - b));
ends = Pall([1 end],:);
end

function scl = case_scl(mesh, sol, geo)
% Stress classification lines through the case (right half, symmetric):
% nose centreline and mid-nose, vault diagonal, lateral wall at 1/4, 1/2,
% 3/4 of the WP height, plasma-side plate. Linearized Pm and Pm+Pb (Tresca).
xy = mesh.xy; ce = find(mesh.t6_mat == 4);
Tc = mesh.t6(ce, 1:3);
V = geo.cav(1:end-1,:);
Y1 = min(V(:,2)); Y2 = max(V(:,2));
xb = max(V(abs(V(:,2)-Y1) < 1e-9, 1));                 % cavity bottom corner x
h = geo.half; tw = geo.t_w; R = geo.R_in;
d = [sin(h), cos(h)]; nr = [cos(h), -sin(h)];
xoff = @(y) (y*sin(h) - tw)/cos(h);                  % case flank (offset line) at height y
corner = sqrt(R^2 - tw^2)*d - tw*nr;                 % arc / case-flank corner
lines = {'nose centreline', [0 Y1], [0 R];
         'nose mid',        [xb/2 Y1], [xb/2 sqrt(R^2-(xb/2)^2)];
         'vault diagonal',  [xb Y1], corner};
for f = [0.25 0.5 0.75]
    y = Y1 + f*(Y2 - Y1);
    lines(end+1,:) = {sprintf('side wall %.0f%%', 100*f), [cav_x_at(V, y) y], [xoff(y) y]}; %#ok<AGROW>
end
lines(end+1,:) = {'plasma-side plate', [0 Y2], [0 geo.Ri_]};
scl = struct('name', lines(:,1), 'P0', lines(:,2), 'P1', lines(:,3));
ns = 21;
for q = 1:numel(scl)
    P0 = scl(q).P0; P1 = scl(q).P1;
    s = linspace(0.002, 0.998, ns)';
    P = P0 + s*(P1 - P0);
    [idx, bc] = tsearchn(xy, Tc, P);
    S = nan(ns, 4);
    for k = 1:ns
        if isnan(idx(k)), continue, end
        e = ce(idx(k)); nodes = mesh.t6(e,:);
        S(k,:) = elem_stress_at(xy(nodes,:), sol.U, nodes, sol.eps_z, sol.tD{e}, sol.teth{e}, @t6_dN, bc(k,[2 3]));
    end
    ok = ~any(isnan(S), 2);
    L = norm(P1 - P0); z = (s - 0.5)*L;
    m = trapz(s(ok)*L, S(ok,:))/(L*(s(find(ok,1,'last')) - s(find(ok,1))));
    b = 6/L^2*trapz(s(ok)*L, S(ok,:).*z(ok));
    scl(q).Pm = tresca(m);
    scl(q).PmPb = max(tresca(m + b), tresca(m - b));
    scl(q).S_membrane = m;
    scl(q).valid_fraction = mean(ok);
end
end

function x = cav_x_at(V, y)
% right-most intersection of the horizontal line y with the cavity polygon
x = -Inf;
W = [V; V(1,:)];
for s = 1:size(W,1)-1
    A = W(s,:); B = W(s+1,:);
    if (A(2)-y)*(B(2)-y) <= 0 && A(2) ~= B(2)
        x = max(x, A(1) + (y - A(2))*(B(1)-A(1))/(B(2)-A(2)));
    end
end
end

function fom = surr_fom(out, p)
% Primary-stress figure of merit (ITER / ASME III style): Pm <= Sm and
% Pm+Pb <= 1.5 Sm on the jacket walls and on the case SCLs, with
% Sm = S_amm_JT (jacket) and S_amm_VT (case).
[fom.jacket_Pm, it] = max([out.turn.Pm]);
[fom.jacket_PmPb, ib] = max([out.turn.PmPb]);
fom.jacket_Pm_at = [out.turn(it).layer, out.turn(it).col];
fom.jacket_PmPb_at = [out.turn(ib).layer, out.turn(ib).col];
[fom.jacket_peak, ip] = max([out.turn.peak]);
fom.jacket_peak_at = [out.turn(ip).layer, out.turn(ip).col];
[fom.case_Pm, ic] = max([out.case_scl.Pm]);
[fom.case_PmPb, id] = max([out.case_scl.PmPb]);
fom.case_Pm_at = out.case_scl(ic).name; fom.case_PmPb_at = out.case_scl(id).name;
Smj = get_or(p, 'S_amm_JT', 667e6); Smc = get_or(p, 'S_amm_VT', 667e6);
fom.util = [fom.jacket_Pm/Smj, fom.jacket_PmPb/(1.5*Smj), fom.case_Pm/Smc, fom.case_PmPb/(1.5*Smc)];
fom.util_names = {'jacket Pm/Sm', 'jacket (Pm+Pb)/1.5Sm', 'case Pm/Sm', 'case (Pm+Pb)/1.5Sm'};
fom.max_util = max(fom.util);
fom.ok = fom.max_util <= 1;
end

function S = elem_stress_at(X, U, nodes, ez, D, eth, dNf, nat)
ue = reshape([U(2*nodes-1)'; U(2*nodes)'], [], 1);
S = zeros(size(nat,1), 4);
for q = 1:size(nat,1)
    B = elem_B(X, dNf, nat(q,:));
    S(q,:) = (D*(B*[ue; ez] - eth))';
end
end

function t = tresca(S)
% S: n x 4 [sx sy sz txy] -> Tresca stress intensity
c = (S(:,1)+S(:,2))/2; r = sqrt(((S(:,1)-S(:,2))/2).^2 + S(:,4).^2);
P = [c+r, c-r, S(:,3)];
t = max(P,[],2) - min(P,[],2);
end
