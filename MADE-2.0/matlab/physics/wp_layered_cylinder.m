function out = wp_layered_cylinder(row, p, opts)
%WP_LAYERED_CYLINDER Layered thick-cylinder model of the TF inner leg (jacket and case together).
%
%   out = WP_LAYERED_CYLINDER(row, p, opts) solves the inner-leg section of
%   a wedged TF system as a 360-degree cylinder made of concentric rings,
%   in generalized plane strain, under the Lorentz body force of the
%   winding pack and the vertical force T_bf. Equations and assumptions:
%   docs/MODELLO_A_STRATI.md (and the web manual, section 5c).
%
%   Rings, from the bore outwards: nose (case steel, Rk..Rj), ground
%   insulation, the WP layers separated by the inter-layer insulation,
%   ground insulation, plasma-side plate (case steel). Every ring is
%   homogeneous and cylindrically orthotropic; a WP ring contains, along
%   the circumference of one sector, the two lateral case walls, the two
%   lateral ground insulations and the turns of the layer:
%     radial modulus   E_r = sum f_k E_r,k          (components in parallel)
%     hoop modulus     1/E_t = sum f_k / E_t,k       (components in series)
%     axial modulus    E_z = sum f_k E_z,k           (parallel)
%   with f_k the fraction of the sector chord 2 r tan(pi/n_TF). Steel rings
%   are isotropic (Poisson coupling kept), mixed rings have no Poisson
%   coupling.
%
%   In each ring the radial displacement is
%     u(r) = A r^k + B r^-k + u_p(r),   k = sqrt(C22/C11)
%   with u_p the particular solution of the Lorentz body force
%   f(r) = -j B(r) = a1 r + a_1/r (and of the axial strain eps0 for steel
%   rings, none for isotropic rings). Continuity of u and sigma_r at every
%   interface, sigma_r = 0 at the bore and at the plasma-side surface and
%   the axial equilibrium int sigma_z dA = n_TF T_bf give 2N+1 linear
%   equations in the 2N constants A, B and eps0.
%
%   The jacket of every layer is then recovered from the strains of its
%   ring with a cell model (cable, jacket and insulation in series or in
%   parallel): side walls take the radial strain of the ring, top and
%   bottom walls the hoop strain of the turn, all walls the axial strain.
%   Their Tresca intensity is the jacket membrane stress Pm (no local
%   bending: see opts.bending).
%
%   row: design point (Iop, n_layers, n_turns, Cond_w, Cond_h, JT, Ri_,
%   Rk_, type_cable, shape_cable 201). p: input parameters (n_TF, E_jckt,
%   E_case, E_ins, E_cbl_LTS/HTS, nu_steel, dr_plasma_side, GoundIns,
%   INS_grades, turn_insulation_nominal). opts: T_bf (default: tool
%   formula), nu (0.3), bending ('none' or a struct, see below), n_sample
%   (points per ring, 9).
%
%   out: rings (r_in, r_out, type, C, constants), sample (r, u, sigma_r,
%   sigma_t, sigma_z per ring), eps0, layer (per WP layer: ring stresses at
%   mid-height, jacket wall stresses, Pm), nose (linearized Pm, PmPb and
%   the hoop stress at the bore), plate, walls (lateral case walls per
%   layer), check (equilibrium residuals).

if nargin < 3, opts = struct(); end
nu = getf(opts, 'nu', getf(p, 'nu_steel', 0.3));
% share of the lateral case walls in the RADIAL stiffness of a WP ring:
% 1 = WP bonded to the walls (they carry radial load in parallel), 0 = WP
% free to slide on the walls (its radial load stays in the WP down to the
% nose; the walls still close the hoop path)
wall_r = getf(opts, 'wall_radial', 1);
% radial stress of the WP: 'ring' = from the ring solution (WP bonded to
% the case, load shared with the lateral walls); 'column' = the WP is a
% column sliding on the lateral walls: the radial pressure on layer k is
% the Lorentz force of the layers above per unit width of each layer
% (where a layer is narrower than the one above, the overhanging turns
% rest on the case shoulder), p_k = sum_{j<k} F_j/W_j + F_k/(2 W_k)
radial = getf(opts, 'radial', 'ring');
ns = getf(opts, 'n_sample', 9);
Mu0 = 4e-7*pi;
nTF = p.n_TF; th = 2*pi/nTF;
nl = row.n_layers;
nt = row.n_turns(1:nl); Cw = row.Cond_w(1:nl); Ch = row.Cond_h(1:nl); JT = row.JT(1:nl);
tg = wp_turn_geometry(row, p);
tins = tg.tins;
Ej = p.E_jckt; Ec = p.E_case; Ei = p.E_ins;
Ecb = p.E_cbl_LTS*ones(1, nl);
if isfield(row, 'type_cable')
    tc = row.type_cable;
    for k = 1:nl, if strcmp(tc{k}, 'HTS'), Ecb(k) = p.E_cbl_HTS; end, end
end
ins = p.INS_grades; git = p.GoundIns; dps = p.dr_plasma_side;
Ri_ = row.Ri_; Rk_ = row.Rk_;

% ---- cell (turn) moduli: cable, jacket, insulation ----------------------
SCw = Cw - 2*JT - 2*tins; SCh = Ch - 2*JT - 2*tins;
% radial loading: side walls (jacket, insulation) in parallel with the
% central column, where insulation, jacket top/bottom walls and cable are
% in series
Emid_r = 1./(2*tins./(Ch*Ei) + 2*JT./(Ch*Ej) + SCh./(Ch.*Ecb));      % central column, series
Er_t = (2*JT*Ej + 2*tins*Ei)./Cw + (SCw./Cw).*Emid_r;
% hoop (toroidal) loading: top/bottom walls in parallel with the central
% band, where insulation, jacket side walls and cable are in series
Emid_t = 1./(2*tins./(Cw*Ei) + 2*JT./(Cw*Ej) + SCw./(Cw.*Ecb));
Et_t = (2*JT*Ej + 2*tins*Ei)./Ch + (SCh./Ch).*Emid_t;
Acell = Cw.*Ch;
Ez_t = (tg.A_jacket*Ej + tg.A_cable.*Ecb + (Acell - tg.A_jacket - tg.A_cable)*Ei)./Acell;

% ---- rings (inner to outer) ---------------------------------------------
Re = zeros(1, nl); Rii = zeros(1, nl);
Re(1) = Ri_ - dps - git;
for k = 1:nl
    Rii(k) = Re(k) - Ch(k);
    if k < nl, Re(k+1) = Rii(k) - ins; end
end
Rj = Rii(nl) - git;
R = struct('r1', {}, 'r2', {}, 'type', {}, 'layer', {}, 'C', {}, 'jz', {});
R(end+1) = ring(Rk_, Rj, 'nose', 0, steel(Ec, nu), 0);
R(end+1) = ring(Rj, Rii(nl), 'ground', nl, mixed(Rj, Rii(nl), nt(nl)*Cw(nl), git, 0, Ei, Ei, Ei, Ec, th), 0);
for k = nl:-1:1
    W = nt(k)*Cw(k);
    I_ring = nTF*nt(k)*row.Iop;
    R(end+1) = ring(Rii(k), Re(k), 'layer', k, ...
        mixed(Rii(k), Re(k), W, git, 1, Er_t(k), Et_t(k), Ez_t(k), Ec, th, Ei, wall_r), ...
        I_ring/(pi*(Re(k)^2 - Rii(k)^2))); %#ok<AGROW>
    if k > 1
        R(end+1) = ring(Re(k), Rii(k-1), 'interlayer', k, ...
            mixed(Re(k), Rii(k-1), max(nt(k), nt(k-1))*Cw(k), git, 0, Ei, Ei, Ei, Ec, th), 0); %#ok<AGROW>
    end
end
R(end+1) = ring(Re(1), Ri_ - dps, 'ground', 1, mixed(Re(1), Ri_ - dps, nt(1)*Cw(1), git, 0, Ei, Ei, Ei, Ec, th), 0);
R(end+1) = ring(Ri_ - dps, Ri_, 'plate', 0, steel(Ec, nu), 0);
N = numel(R);

% ---- Lorentz body force f(r) = -j B(r) = a1 r + am1 / r -------------------
Ienc = 0;                                     % current inside the inner radius of the ring
for i = 1:N
    j = R(i).jz;
    % B(r) = mu0 (Ienc + j pi (r^2 - r1^2)) / (2 pi r)
    R(i).a1 = -j*Mu0*j/2;
    R(i).am1 = -j*Mu0*(Ienc - j*pi*R(i).r1^2)/(2*pi);
    Ienc = Ienc + j*pi*(R(i).r2^2 - R(i).r1^2);
end
g = compute_operating_params(p);
T_bf = getf(opts, 'T_bf', 0.5*(g.k_bf*nTF*(sum(nt)*row.Iop)^2*Mu0/(2*pi)));
Nz = nTF*T_bf;

% ---- linear system in [A_1 B_1 ... A_N B_N eps0] -------------------------
nu_ = 2*N + 1;
M = zeros(nu_); rhs = zeros(nu_, 1); e = 0;
% sigma_r = 0 at the bore
[a, b] = sr_row(R(1), R(1).r1, N, 1); e = e + 1; M(e,:) = a; rhs(e) = -b;
for i = 1:N-1
    r = R(i).r2;
    [a1u, b1u] = u_row(R(i), r, N, i); [a2u, b2u] = u_row(R(i+1), r, N, i+1);
    e = e + 1; M(e,:) = a1u - a2u; rhs(e) = b2u - b1u;
    [a1s, b1s] = sr_row(R(i), r, N, i); [a2s, b2s] = sr_row(R(i+1), r, N, i+1);
    e = e + 1; M(e,:) = a1s - a2s; rhs(e) = b2s - b1s;
end
[a, b] = sr_row(R(N), R(N).r2, N, N); e = e + 1; M(e,:) = a; rhs(e) = -b;
% axial equilibrium: sum int sigma_z 2 pi r dr = Nz
[xg, wg] = gauss_pts(6);
az = zeros(1, nu_); bz = 0;
for i = 1:N
    rr = (R(i).r1 + R(i).r2)/2 + (R(i).r2 - R(i).r1)/2*xg;
    ww = (R(i).r2 - R(i).r1)/2*wg;
    for q = 1:numel(rr)
        [aa, bb] = sz_row(R(i), rr(q), N, i);
        az = az + aa*2*pi*rr(q)*ww(q); bz = bz + bb*2*pi*rr(q)*ww(q);
    end
end
e = e + 1; M(e,:) = az; rhs(e) = Nz - bz;
x = M\rhs;
eps0 = x(end);

% ---- results --------------------------------------------------------------
out = struct('rings', R, 'eps0', eps0, 'T_bf', T_bf);
smp = struct('r', {}, 'u', {}, 'sr', {}, 'st', {}, 'sz', {}, 'er', {}, 'et', {});
for i = 1:N
    rr = linspace(R(i).r1, R(i).r2, ns);
    s = struct('r', rr, 'u', 0*rr, 'sr', 0*rr, 'st', 0*rr, 'sz', 0*rr, 'er', 0*rr, 'et', 0*rr);
    for q = 1:ns
        [s.u(q), s.er(q), s.et(q), s.sr(q), s.st(q), s.sz(q)] = state(R(i), rr(q), x, i, eps0);
    end
    smp(i) = s;
end
out.sample = smp;
% column pressure per layer (radial Lorentz force per unit length of each
% layer, per coil, from the smeared field)
F_lay = zeros(1, nl); W_lay = nt.*Cw;
for i = 1:N
    if ~strcmp(R(i).type, 'layer'), continue, end
    rr = (R(i).r1 + R(i).r2)/2 + (R(i).r2 - R(i).r1)/2*xg; ww = (R(i).r2 - R(i).r1)/2*wg;
    F_lay(R(i).layer) = -sum((R(i).a1*rr + R(i).am1./rr).*2*pi.*rr.*ww)/nTF;
end
% friction on the lateral walls (Janssen): the hoop (wedging) compression
% sigma_t of the WP ring presses the WP on the two lateral walls; with the
% WP / case friction coefficient mu, each layer of height h_k hands over to
% the case up to 2 mu |sigma_t| h_k per unit length, i.e. 2 mu |sigma_t|
% h_k / W_k of column pressure (never more than the pressure itself)
mu = getf(opts, 'mu_case', getf(p, 'mu_case', 0.2));
st_lay = zeros(1, nl);
for i = 1:N
    if ~strcmp(R(i).type, 'layer'), continue, end
    [~, ~, ~, ~, st_lay(R(i).layer)] = state(R(i), (R(i).r1 + R(i).r2)/2, x, i, eps0);
end
p_col = zeros(1, nl); p_in = 0;
for k = 1:nl
    dp_load = F_lay(k)/W_lay(k);
    dp_fric = 2*mu*abs(min(st_lay(k), 0))*Ch(k)/W_lay(k);
    p_col(k) = max(p_in + dp_load/2 - dp_fric/2, 0);
    p_in = max(p_in + dp_load - dp_fric, 0);
end
out.F_layer = F_lay; out.p_column = p_col;
% WP layers: ring stresses at mid-height and jacket walls
lay = struct('k', {}, 'r', {}, 'sr', {}, 'st', {}, 'sz', {}, 'er', {}, 'et', {}, ...
    'J_side', {}, 'J_top', {}, 'Pm', {}, 'wall', {});
for i = 1:N
    if ~strcmp(R(i).type, 'layer'), continue, end
    k = R(i).layer; rm = (R(i).r1 + R(i).r2)/2;
    [~, er, et, sr, st, sz] = state(R(i), rm, x, i, eps0);
    % hoop: the ring components are in series -> same hoop stress; strain
    % of the turn = sigma_t / Et_turn
    et_turn = st/Et_t(k);
    if strcmp(radial, 'column')
        er = -p_col(k)/Er_t(k);             % radial strain of the turn under the column pressure
    end
    sm_r = Emid_r(k)*er;                    % stress of the central column (radial)
    sm_t = Emid_t(k)*et_turn;               % stress of the central band (hoop)
    % each wall: strain imposed in its 'parallel' direction, stress
    % imposed in its 'series' direction, axial strain eps0 (generalized
    % plane strain) -> 3D isotropic Hooke for the other two stresses
    [sp, sz_s] = wall_hooke(er, sm_t, eps0, Ej, nu);
    J_side = [sp, sm_t, sz_s];              % side walls: radial (parallel), hoop (series), axial
    [sp, sz_t] = wall_hooke(et_turn, sm_r, eps0, Ej, nu);
    J_top  = [sm_r, sp, sz_t];              % top/bottom walls: radial (series), hoop (parallel), axial
    Pm = max(tresca3(J_side), tresca3(J_top));
    [sp, sz_w] = wall_hooke(er, st, eps0, Ec, nu);
    wall = [sp, st, sz_w];                  % lateral case wall: radial (parallel), hoop (series)
    lay(k) = struct('k', k, 'r', rm, 'sr', sr, 'st', st, 'sz', sz, 'er', er, 'et', et, ...
        'J_side', J_side, 'J_top', J_top, 'Pm', Pm, 'wall', wall);
end
out.layer = lay;
out.Pm_jacket = [lay.Pm];
out.cell = struct('Er', Er_t, 'Et', Et_t, 'Ez', Ez_t, 'Emid_r', Emid_r, 'Emid_t', Emid_t);
% nose and plate: linearization through the thickness (like the FE SCLs)
out.nose = linearize(smp(1));
out.plate = linearize(smp(N));
% lateral wall Pm per layer
out.walls_Pm = arrayfun(@(L) tresca3(L.wall), lay);
% checks: radial force balance of the WP and axial force
Fr = 0;
for i = 1:N
    if R(i).jz == 0, continue, end
    rr = (R(i).r1 + R(i).r2)/2 + (R(i).r2 - R(i).r1)/2*xg; ww = (R(i).r2 - R(i).r1)/2*wg;
    Fr = Fr + sum((R(i).a1*rr + R(i).am1./rr).*2*pi.*rr.*ww);
end
% hoop balance of the whole section (free surfaces): int sigma_t dr = int f r dr
Ht = 0;
for i = 1:N
    Ht = Ht + trapz(smp(i).r, smp(i).st);
end
Hf = 0;
for i = 1:N
    rr = (R(i).r1 + R(i).r2)/2 + (R(i).r2 - R(i).r1)/2*xg; ww = (R(i).r2 - R(i).r1)/2*wg;
    Hf = Hf + sum((R(i).a1*rr + R(i).am1./rr).*rr.*ww);
end
out.check = struct('residual', norm(M*x - rhs)/max(norm(rhs), eps), 'F_centering_per_coil', Fr/nTF, ...
    'N_z_per_coil', Nz/nTF, 'hoop_balance', abs(Ht - Hf)/max(abs(Hf), eps));
end

% ===========================================================================
function s = ring(r1, r2, type, layer, C, jz)
s = struct('r1', r1, 'r2', r2, 'type', type, 'layer', layer, 'C', C, 'jz', jz);
end

function C = steel(E, nu)
c = E/((1 + nu)*(1 - 2*nu));
C = c*[1-nu nu nu; nu 1-nu nu; nu nu 1-nu];
end

function C = mixed(r1, r2, W, git, is_wp, Er_w, Et_w, Ez_w, Ec, th, Ei, wall_r)
% ring of one sector: lateral case walls + lateral ground insulation +
% the WP (or the insulation across the cavity width W)
if nargin < 11, Ei = Er_w; end
if nargin < 12, wall_r = 1; end
rm = (r1 + r2)/2;
chord = 2*rm*tan(th/2);
w_ins = 2*git; w_wall = max(chord - W - w_ins, 0);
f = [w_wall, w_ins, W]/chord;
if is_wp
    Ek = [Ec Ei Er_w; Ec Ei Et_w; Ec Ei Ez_w];
else
    Ek = [Ec Ei Er_w; Ec Ei Et_w; Ec Ei Ez_w];     % W filled with insulation (Er_w = Et_w = Ez_w = E_ins)
end
Er = f*(Ek(1,:).*[wall_r 1 1])';
Et = 1/(f*(1./Ek(2,:))');
Ez = f*Ek(3,:)';
C = diag([Er Et Ez]);
end

function [kk, P, Q, logq, Qe, logqe] = coeffs(Rg)
% homogeneous exponent and particular-solution coefficients of one ring:
% u_p = P r^3 + Q (r or r ln r) + eps0 * Qe (r or r ln r)
C = Rg.C; C11 = C(1,1); C22 = C(2,2); dC = C(1,3) - C(2,3);
kk = sqrt(C22/C11);
P = -Rg.a1/(9*C11 - C22);
logq = abs(C11 - C22) < 1e-6*C11;
if logq, Q = -Rg.am1/(2*C11); else, Q = -Rg.am1/(C11 - C22); end
logqe = logq;
if logq, Qe = -dC/(2*C11); else, Qe = -dC/(C11 - C22); end
end

function [ub, dub, up, dup, ue, due] = basis(Rg, r)
[kk, P, Q, logq, Qe] = coeffs(Rg);
ub = [r^kk, r^(-kk)]; dub = [kk*r^(kk-1), -kk*r^(-kk-1)];
if logq
    up = P*r^3 + Q*r*log(r); dup = 3*P*r^2 + Q*(log(r) + 1);
    ue = Qe*r*log(r); due = Qe*(log(r) + 1);
else
    up = P*r^3 + Q*r; dup = 3*P*r^2 + Q;
    ue = Qe*r; due = Qe;
end
end

function [a, b] = u_row(Rg, r, N, i)
[ub, ~, up, ~, ue] = basis(Rg, r);
a = zeros(1, 2*N + 1); a(2*i-1:2*i) = ub; a(end) = ue; b = up;
end

function [a, b] = sr_row(Rg, r, N, i)
[ub, dub, up, dup, ue, due] = basis(Rg, r);
C = Rg.C;
a = zeros(1, 2*N + 1);
a(2*i-1:2*i) = C(1,1)*dub + C(1,2)*ub/r;
a(end) = C(1,1)*due + C(1,2)*ue/r + C(1,3);
b = C(1,1)*dup + C(1,2)*up/r;
end

function [a, b] = sz_row(Rg, r, N, i)
[ub, dub, up, dup, ue, due] = basis(Rg, r);
C = Rg.C;
a = zeros(1, 2*N + 1);
a(2*i-1:2*i) = C(3,1)*dub + C(3,2)*ub/r;
a(end) = C(3,1)*due + C(3,2)*ue/r + C(3,3);
b = C(3,1)*dup + C(3,2)*up/r;
end

function [u, er, et, sr, st, sz] = state(Rg, r, x, i, eps0)
[ub, dub, up, dup, ue, due] = basis(Rg, r);
AB = x(2*i-1:2*i);
u = ub*AB + up + eps0*ue;
er = dub*AB + dup + eps0*due;
et = u/r;
C = Rg.C;
s = C*[er; et; eps0];
sr = s(1); st = s(2); sz = s(3);
end

function L = linearize(s)
% membrane (mean) and bending (linear part) of sigma_r, sigma_t, sigma_z
% through the ring thickness; Tresca of membrane and of membrane +/- bending
t = s.r(end) - s.r(1); z = s.r - (s.r(1) + s.r(end))/2;
S = [s.sr; s.st; s.sz];
m = trapz(s.r, S, 2)/t;
bb = 6/t^2*trapz(s.r, S.*z, 2);
L = struct('Pm', tresca3(m'), 'PmPb', max(tresca3((m + bb)'), tresca3((m - bb)')), ...
    'membrane', m', 'bending', bb', 'st_inner', s.st(1), 'st_outer', s.st(end));
end

function [s_par, s_z] = wall_hooke(e_par, s_ser, e_z, E, nu)
% isotropic Hooke with the strain e_par in one in-plane direction, the
% stress s_ser in the other in-plane direction and the axial strain e_z:
%   e_par = (s_par - nu (s_ser + s_z))/E,  e_z = (s_z - nu (s_par + s_ser))/E
A = [1 -nu; -nu 1];
x = A\[E*e_par + nu*s_ser; E*e_z + nu*s_ser];
s_par = x(1); s_z = x(2);
end

function t = tresca3(s)
t = max(s) - min(s);
end

function [x, w] = gauss_pts(n)
b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D)); w = 2*V(1,i)'.^2;
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
