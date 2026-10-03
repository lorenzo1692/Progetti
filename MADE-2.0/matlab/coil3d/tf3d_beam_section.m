function [sec, info] = tf3d_beam_section(row, p, opts)
%TF3D_BEAM_SECTION Beam section of the global model from the MADE TF design point.
%
%   [sec, info] = TF3D_BEAM_SECTION(row, p) returns the section of the 3D
%   beams of TF3D_GLOBAL_MODEL for the design point row, instead of the
%   generic hollow square of STR_360. The section is the inner-leg
%   equatorial cross-section sized by MADE:
%
%     - case: wedge (trapezoid) between Rk_ and Ri_, half-width r*tan(pi/n_TF);
%     - winding pack: one rectangle per layer (n_turns*Cond_w x Cond_h,
%       with ground insulation around the pack), as SCAN_WP_DESIGNS
%       stacks them from Ri_ - dr_plasma_side - GoundIns inwards, with the
%       LONGITUDINAL modulus of the turn (parallel of jacket, cable and
%       insulation areas of WP_TURN_GEOMETRY) - the beam axis follows the
%       conductor;
%     - ground insulation with E_ins.
%   Properties are modulus-weighted with E_ref = E_case (sec.E):
%     A, Iz (bending in the coil plane: radial lever arm), Iy (out-of-plane
%     bending: toroidal lever arm) about the elastic centroid;
%     J from Bredt's formula on the closed case (plasma-side plate, two
%     lateral walls, nose; the WP is neglected in torsion);
%     Asy (in-plane shear: lateral walls), Asz (out-of-plane shear:
%     plasma-side plate + nose).
%   The same section is used along the whole coil (as STR_360): the outer
%   leg is assumed at least as stiff as the inner leg; opts.sec_scale_outer
%   is left for a later version.
%
%   info: areas, E_WP_long, centroid radius, eccentricity w.r.t. the current
%   centre of the WP (the beam axis), polygons for plotting.

if nargin < 3, opts = struct(); end
if ~isstruct(row), row = table2struct(row); end
Eu = 1; if p.E_case < 1e6, Eu = 1e9; end          % input moduli in GPa
E_case = p.E_case*Eu; E_j = p.E_jckt*Eu; E_ins = p.E_ins*Eu;
nu = 0.3; if isfield(opts, 'nu'), nu = opts.nu; end
th2 = pi/p.n_TF;
nl = row.n_layers; nt = row.n_turns(1:nl); cw = row.Cond_w(1:nl); ch = row.Cond_h(1:nl);
Ri = row.Ri_; Rk = row.Rk_; GI = p.GoundIns; dps = p.dr_plasma_side;

% turn longitudinal modulus per layer
tg = wp_turn_geometry(row, p);
E_cbl = zeros(1, nl);
tc = row.type_cable; if ischar(tc), tc = {tc}; end
for k = 1:nl
    if strcmp(tc{min(k, numel(tc))}, 'HTS'), E_cbl(k) = p.E_cbl_HTS*Eu; else, E_cbl(k) = p.E_cbl_LTS*Eu; end
end
A_cell = cw.*ch;
A_ins_turn = A_cell - tg.A_jacket - tg.A_cable;
E_turn = (E_j*tg.A_jacket + E_cbl.*tg.A_cable + E_ins*A_ins_turn)./A_cell;

% shapes: [polygon in (r, x) with x toroidal], modulus
shapes = {};
shapes{end+1} = {[Rk -Rk*tan(th2); Ri -Ri*tan(th2); Ri Ri*tan(th2); Rk Rk*tan(th2)], E_case};
r_wp_top = Ri - dps - GI;                             % outer radius of layer 1
Re = r_wp_top;
W = nt.*cw;
for k = 1:nl
    r2 = Re; r1 = Re - ch(k);
    w = W(k)/2;
    % ground insulation band around this layer (lateral), WP layer
    shapes{end+1} = {rect(r1, r2, w + GI), E_ins - E_case}; %#ok<AGROW>
    shapes{end+1} = {rect(r1, r2, w), E_turn(k) - E_ins}; %#ok<AGROW>
    if k < nl
        gap = p.INS_grades;
        shapes{end+1} = {rect(r1 - gap, r1, max(W(k), W(k+1))/2 + GI), E_ins - E_case}; %#ok<AGROW>
        Re = r1 - gap;
    else
        Re = r1;
    end
end
r_bot = Re - GI;
shapes{end+1} = {rect(r_wp_top, r_wp_top + GI, W(1)/2 + GI), E_ins - E_case};      % top GI band
shapes{end+1} = {rect(r_bot, Re, W(nl)/2 + GI), E_ins - E_case};                                           % bottom GI band

% modulus-weighted moments (E_ref = E_case)
A = 0; Sr = 0; Sx = 0; Irr = 0; Ixx = 0;
for i = 1:numel(shapes)
    [a, sr, sx, irr, ixx] = poly_moments(shapes{i}{1});
    n = shapes{i}{2}/E_case;
    A = A + n*a; Sr = Sr + n*sr; Sx = Sx + n*sx; Irr = Irr + n*irr; Ixx = Ixx + n*ixx;
end
rc = Sr/A; xc = Sx/A;
Iz = Irr - A*rc^2;                                    % radial lever arm: in-plane bending
Iy = Ixx - A*xc^2;                                    % toroidal lever arm: out-of-plane bending

% torsion (Bredt) and shear areas from the case walls
Rj = row.Rj_; nose = Rj - Rk;
lat = row.lateral_w;
WPh = row.WP_h; WPw = row.WP_w;
t_ps = dps + GI;
b_top = 2*(Ri - t_ps/2)*tan(th2); b_bot = 2*(Rk + nose/2)*tan(th2);
h_m = (Ri - t_ps/2) - (Rk + nose/2);
A_m = (b_top + b_bot)/2*h_m;
L_side = hypot(h_m, (b_top - b_bot)/2);
J = 4*A_m^2/(b_top/t_ps + b_bot/nose + 2*L_side/lat);
Asy = 2*lat*(Ri - Rk);
Asz = t_ps*b_top + nose*b_bot;

sec = struct('E', E_case, 'nu', nu, 'A', A, 'Iy', Iy, 'Iz', Iz, 'J', J, 'Asy', Asy, 'Asz', Asz);
r_current = (r_wp_top + Re)/2;                        % mid-depth of the WP
info = struct('A_transformed', A, 'r_centroid', rc, 'r_wp_centre', r_current, ...
    'eccentricity', rc - r_current, 'E_turn_long', E_turn, 'nose', nose, ...
    'lateral_w', lat, 't_plasma_side', t_ps, 'WP_h', WPh, 'WP_w', WPw, 'shapes', {shapes});
end

function P = rect(r1, r2, hw)
P = [r1 -hw; r2 -hw; r2 hw; r1 hw];
end

function [a, sr, sx, irr, ixx] = poly_moments(P)
% area, first and second moments of a simple polygon (r, x), about r = 0, x = 0
r = P(:,1); x = P(:,2); r2 = r([2:end 1]); x2 = x([2:end 1]);
c = r.*x2 - r2.*x;
a = sum(c)/2;
sr = sum((r + r2).*c)/6;
sx = sum((x + x2).*c)/6;
irr = sum((r.^2 + r.*r2 + r2.^2).*c)/12;
ixx = sum((x.^2 + x.*x2 + x2.^2).*c)/12;
if a < 0, a = -a; sr = -sr; sx = -sx; irr = -irr; ixx = -ixx; end
end
