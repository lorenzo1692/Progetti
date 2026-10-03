function s = tf3d_section_from_design(row, p, opts)
%TF3D_SECTION_FROM_DESIGN Beam section of the TF coil for the 3D global model, from the 2D design.
%
%   s = TF3D_SECTION_FROM_DESIGN(row, p) integrates the inner-leg cross
%   section of the design point (x toroidal, y radial, as the 2D sizing):
%
%   - case: the sector |x| <= y tan(pi/n_TF) between the bore arc R_bore =
%     Rk_/cos(pi/n_TF) and the plasma-side face y = Ri_, minus the WP
%     cavity: the convex hull of the layers (width n_turns*Cond_w, height
%     Cond_h, inter-grade insulation INS_grades) widened by the ground
%     insulation GoundIns, as the 2D FE (WP_MECH_SURROGATE) and ANSYS;
%   - jackets: the jacket steel of every layer (WP_TURN_GEOMETRY), smeared
%     over the layer band, weighted by E_jckt/E_case.
%
%   Cable, insulation and filler carry no load in the beam (as the steel
%   area of the 2D sizing, S_z = T_bf/(A_jacket + A_case)).
%
%   s: E, nu (case steel), A [m^2] (steel, case-equivalent), A_case,
%   A_jacket, y_c (radius of the centroid), Iz (bending in the plane of the
%   coil, about the toroidal axis: integral of (y - y_c)^2), Iy (bending
%   out of the plane, about the radial axis: integral of x^2), J (torsion,
%   Bredt formula for the closed case box with its mean wall thickness),
%   Asy, Asz (shear areas, A_case/2 each), c_r, c_t (largest radial and
%   toroidal distance of the steel from the centroid, for the fibre
%   stresses), t_vault (radial steel thickness crossing the coil mid
%   plane: nose + plasma-side plate + 2 JT per layer; thickness of the
%   vault shells of the global model), R_bore, Nose, radial_build.
%   opts: n_grid (points per direction, default 500).
%
%   This is the section of the inner leg. The global model uses it along
%   the whole coil, as STR_360 uses one HREC section (TF3D_GLOBAL_MODEL).

if nargin < 3, opts = struct(); end
ng = 500; if isfield(opts, 'n_grid'), ng = opts.n_grid; end
if ~isstruct(row), row = table2struct(row); end
nl = row.n_layers;
nt = row.n_turns(1:nl); Cw = row.Cond_w(1:nl); Ch = row.Cond_h(1:nl); JT = row.JT(1:nl);
th2 = pi/p.n_TF;
R_bore = row.Rk_/cos(th2);
Ri = row.Ri_;
git = p.GoundIns; ins = p.INS_grades; dps = p.dr_plasma_side;
tg = wp_turn_geometry(row, p);
E_case = modulus(p.E_case); E_jk = modulus(p.E_jckt);

% layer bands (top = plasma side) and cavity half widths
Re = Ri - dps - git;                       % top of layer 1
y_top = zeros(1, nl); y_bot = y_top;
for k = 1:nl
    y_top(k) = Re; y_bot(k) = Re - Ch(k);
    Re = y_bot(k) - ins;
end
W = nt.*Cw;                                 % layer width
hw = W/2 + git;                             % cavity half width at layer k
Rj = y_bot(nl) - git;                       % bottom of the cavity

% grid over the sector
xm = Ri*tan(th2);
y0 = R_bore*cos(th2);
x = linspace(-xm, xm, ng); y = linspace(y0, Ri, ng);
dx = x(2) - x(1); dy = y(2) - y(1);
[X, Y] = meshgrid(x, y);
in_sector = abs(X) <= Y*tan(th2) & hypot(X, Y) >= R_bore & Y <= Ri;
% cavity: convex hull of the layers widened by the ground insulation
hx = []; hy = [];
for k = 1:nl
    ytop = y_top(k) + (k == 1)*git; ybot = y_bot(k) - (k == nl)*git;
    hx = [hx, -hw(k), hw(k), hw(k), -hw(k)]; %#ok<AGROW>
    hy = [hy, ybot, ybot, ytop, ytop]; %#ok<AGROW>
end
kh = convhull(hx, hy); hull = [hx(kh(:)); hy(kh(:))]';
cav = inpolygon(X, Y, hull(:,1), hull(:,2));
w = zeros(size(X));
for k = 1:nl
    band = Y >= y_bot(k) & Y <= y_top(k);
    fk = nt(k)*tg.A_jacket(k)/(W(k)*Ch(k))*E_jk/E_case;     % smeared jacket, case-equivalent
    w(band & abs(X) <= W(k)/2) = fk;
end
case_m = in_sector & ~cav;
w(case_m) = 1;
w(~in_sector) = 0;
dA = dx*dy;
A = sum(w(:))*dA;
A_case = sum(case_m(:))*dA;
yc = sum(w(:).*Y(:))*dA/A;
Iz = sum(w(:).*(Y(:) - yc).^2)*dA;
Iy = sum(w(:).*X(:).^2)*dA;
st = w > 0;
c_r = max(abs(Y(st) - yc)); c_t = max(abs(X(st)));

% torsion: Bredt, closed box of the case with its mean wall thickness
A_cav = sum(cav(:) & in_sector(:))*dA;
side = Ri/cos(th2) - R_bore;
P_out = 2*xm + 2*side + 2*th2*R_bore;
P_cav = sum(sqrt(sum(diff(hull).^2, 2)));
P_m = (P_out + P_cav)/2;
A_m = (A_case + 2*A_cav)/2;                 % area enclosed by the mid line of the walls
t_eff = A_case/P_m;
J = 4*A_m^2*t_eff/P_m;

Nose = Rj - R_bore;
t_vault = Nose + dps + sum(2*JT)*E_jk/E_case;
nu = 0.3; if isfield(p, 'nu_steel') && ~isempty(p.nu_steel), nu = p.nu_steel; end
s = struct('E', E_case, 'nu', nu, 'A', A, 'A_case', A_case, 'A_jacket', sum(nt.*tg.A_jacket), ...
    'y_c', yc, 'Iz', Iz, 'Iy', Iy, 'J', J, 'Asy', A_case/2, 'Asz', A_case/2, 'c_r', c_r, 'c_t', c_t, ...
    't_vault', t_vault, 't_case_mean', t_eff, 'R_bore', R_bore, 'Nose', Nose, 'radial_build', Ri - R_bore, ...
    'Rj', Rj, 'A_m', A_m);
end

function E = modulus(E)
% input files give the moduli in GPa
if E < 1e6, E = E*1e9; end
end
