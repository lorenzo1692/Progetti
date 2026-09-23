function out = compute_discrete_field_profile(row, p)
%COMPUTE_DISCRETE_FIELD_PROFILE Exact 2D discrete-coil field at every turn.
%
%   out = COMPUTE_DISCRETE_FIELD_PROFILE(row, p) replaces the smooth/
%   smeared field model used during the scan (B_layers, a linear
%   interpolation between B_TF and 0 assuming n_TF continuous current
%   sheets via Ampere's law) with an exact 2D magnetostatic calculation in
%   the same idealization as the generalized-plane-strain FEM model
%   (infinite straight conductors along Z):
%
%   - Coil 1 (the one analysed): every turn is a uniform-current
%     RECTANGLE with the cable's own cross-section (SC_w x SC_h, i.e. the
%     cell minus jacket and turn insulation), using the closed-form
%     Biot-Savart field of a rectangular current block. This includes the
%     turn's own SELF-FIELD (~mu0*I/(2*pi*a), about 1 T for a 60 kA cable)
%     and the field of its neighbours at their true size, which a line
%     filament model misses.
%   - The other p.n_TF-1 coils are far enough away to be line filaments.
%
%   The field is evaluated on a 3x3 grid on each cable (corners, edge
%   midpoints, centre) and the per-turn PEAK is reported: that is the
%   value that drives Ic degradation and that the FEM reports as the peak
%   BSUM on the conductor. On the TF_2D FEM benchmark (FEM peak 13.489 T)
%   this model gives 13.57 T (+0.6%); the previous centroid/filament model
%   without self-field gave 12.92 T (-4.2%) and under-predicted the low
%   field layers by up to ~50%.
%
%   Only meant for ONE chosen design point (postprocessing/verification),
%   not inside the combinatorial scan.
%
%   row - one row of a DATA table or an equivalent struct (see
%         PLOT_WP_SECTION). Used: Iop, n_layers, n_turns, Cond_w, Cond_h,
%         JT, Ri_, B_TF.
%   p   - machine parameter struct. Required: n_TF, dr_plasma_side,
%         GoundIns, INS_grades, turn_insulation_nominal, Increm.
%
%   out fields (coil 1's turns, one entry per turn, ordered by layer then
%   by position within the layer; plot frame: x = toroidal, y = radial):
%     x, y       - turn centroid [m]
%     layer      - layer (1..n_layers) of each turn
%     B_smooth   - smeared/linear model's field for that layer [T]
%     B_center   - field at the cable centre [T]
%     B_discrete - peak field over the cable cross-section [T]
%     ripple     - B_discrete./B_smooth (>1: smeared model under-predicts)

n_layers = row.n_layers;
n_turns  = row.n_turns(1:n_layers);
Cond_w   = row.Cond_w(1:n_layers);
Cond_h   = row.Cond_h(1:n_layers);
JT       = row.JT(1:n_layers);

theta_TF = 2*pi/p.n_TF;
Mu_0 = 4e-7*pi;
tins = p.turn_insulation_nominal*p.Increm;
I_turn = row.Iop; % same current in every turn, every coil

% Per-layer radial positions, same recursion as search/scan_wp_designs.m
Re = zeros(1, n_layers+1);
Re(1) = row.Ri_ - p.dr_plasma_side - p.GoundIns;
Ri = zeros(1, n_layers);
for k = 1:n_layers
    Ri(k) = Re(k) - Cond_h(k);
    Re(k+1) = Ri(k) - p.INS_grades;
end

% Turn centroids (toroidal xt, radial rt) and cable sizes, winding order
n_tot = sum(n_turns);
xt = zeros(1, n_tot); rt = zeros(1, n_tot); layer_of_turn = zeros(1, n_tot);
sc_w = zeros(1, n_tot); sc_h = zeros(1, n_tot);
i0 = 0;
for k = 1:n_layers
    idx = i0 + (1:n_turns(k));
    xt(idx) = -Cond_w(k)*n_turns(k)/2 + Cond_w(k)/2 + (0:n_turns(k)-1)*Cond_w(k);
    rt(idx) = Ri(k) + Cond_h(k)/2;
    layer_of_turn(idx) = k;
    sc_w(idx) = Cond_w(k) - 2*JT(k) - 2*tins;
    sc_h(idx) = Cond_h(k) - 2*JT(k) - 2*tins;
    i0 = i0 + n_turns(k);
end
if any(sc_w <= 0) || any(sc_h <= 0)
    error('compute_discrete_field_profile:bad_geometry', ...
        'Cable size <= 0 (Cond_w/Cond_h too small for JT + turn insulation).');
end

% Field points: 3x3 grid on each cable. Global frame: X = radial, Y = toroidal
% (coil 1 centred on the X axis, layers flat as in the FEM).
[a, b] = meshgrid([-0.5 0 0.5], [-0.5 0 0.5]);
a = a(:); b = b(:);                          % 9x1 (a: toroidal, b: radial)
PX = rt + b.*sc_h;  PY = xt + a.*sc_w;       % 9 x n_tot
PX = PX(:); PY = PY(:);
owner = repmat(1:n_tot, 9, 1); owner = owner(:);

% Coil 1: analytic rectangular conductors (self-field included)
BX = zeros(size(PX)); BY = zeros(size(PX));
for j = 1:n_tot
    [bx, by] = rect_field(PX, PY, rt(j), xt(j), sc_h(j), sc_w(j), I_turn, Mu_0);
    BX = BX + bx; BY = BY + by;
end

% Coils 2..n_TF: line filaments at the rotated centroids
for c = 2:p.n_TF
    th = (c-1)*theta_TF;
    Xs = rt*cos(th) - xt*sin(th);
    Ys = rt*sin(th) + xt*cos(th);
    dx = PX - Xs; dy = PY - Ys;               % implicit expansion
    d2 = dx.^2 + dy.^2;
    BX = BX + Mu_0*I_turn/(2*pi) * sum(-dy./d2, 2);
    BY = BY + Mu_0*I_turn/(2*pi) * sum( dx./d2, 2);
end
Bmag = sqrt(BX.^2 + BY.^2);

B_peak = accumarray(owner, Bmag, [n_tot 1], @max)';
B_center = Bmag(5:9:end)';                   % a=0, b=0 is the 5th grid point

% Smeared/linear model at each layer, as in scan_wp_designs.m
n_spire_ = zeros(1, n_layers);
n_spire_(1) = sum(n_turns);
for k = 2:n_layers
    n_spire_(k) = n_spire_(k-1) - n_turns(k);
end
B_layers = row.B_TF .* (n_spire_/n_spire_(1));
B_smooth = B_layers(layer_of_turn);

out.x = xt;
out.y = rt;
out.layer = layer_of_turn;
out.B_smooth = B_smooth;
out.B_center = B_center;
out.B_discrete = B_peak;
out.ripple = B_peak ./ B_smooth;
end

function [Bx, By] = rect_field(x, y, xc, yc, w, h, I, Mu_0)
% Closed-form 2D field of a uniform current I (along Z) in the rectangle
% centred at (xc,yc), width w along x, height h along y.
J = I/(w*h);
ua = x - (xc - w/2); ub = x - (xc + w/2);
va = y - (yc - h/2); vb = y - (yc + h/2);
k = Mu_0*J/(2*pi);
By =  k*(G(ua,va) - G(ua,vb) - G(ub,va) + G(ub,vb));
Bx = -k*(G(va,ua) - G(vb,ua) - G(va,ub) + G(vb,ub));
end

function g = G(u, v)
% Primitive of the line-current kernel: 1/2*v*ln(u^2+v^2) + u*atan(v/u)
t1 = 0.5*v.*log(u.^2 + v.^2); t1(v == 0) = 0;
t2 = u.*atan(v./u);           t2(u == 0) = 0;
g = t1 + t2;
end
