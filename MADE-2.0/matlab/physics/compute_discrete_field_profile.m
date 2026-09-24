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
%   - Coil 1 (the one analysed): every cable carries a uniform current
%     density over its real ROUNDED cross-section (SC_w x SC_h, i.e. the
%     cell minus jacket and turn insulation, corner fillet r_SC), decomposed
%     exactly into rectangles with the closed-form Biot-Savart field of a
%     rectangular current block (see WP_FIELD_AT_POINTS). This includes the
%     turn's own SELF-FIELD (~mu0*I/(2*pi*a), about 1 T for a 60 kA cable)
%     and the field of its neighbours at their true size, which a line
%     filament model misses.
%   - The other p.n_TF-1 coils are far enough away to be line filaments.
%
%   The field is evaluated on 48 points along each cable's rounded outline
%   (where the self-field peak lies) plus a 3x3 interior grid, and the
%   per-turn PEAK is reported: that is the value that drives Ic degradation
%   and that the FEM reports as the peak BSUM on the conductor. Against
%   ANSYS: 14.476 T vs 14.482 T (design 7, EM_2D007) and 13.488 T vs
%   13.489 T (TF_2D benchmark); the old centroid/filament model without
%   self-field gave 12.92 T (-4.2%) and under-predicted the low-field
%   layers by up to ~50%.
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

% Cable corner fillet, as the FEM and size_cicc_cable: r_SC = JT clamped
% to [r_SC_min, r_SC_max]
rmin = 2e-3; rmax = 6e-3;
if isfield(p, 'r_SC_min'), rmin = p.r_SC_min; end
if isfield(p, 'r_SC_max'), rmax = p.r_SC_max; end
r_l = min(max(JT, rmin), rmax);
rc = min(r_l(layer_of_turn), 0.49*min(sc_w, sc_h));

% Field points on every cable: its whole boundary (48 points on the rounded
% outline, where the peak of the self-field lies) plus a 3x3 interior grid
% (review C08: 9 points alone do not locate the peak)
np_b = 48;
[PXl, PYl] = deal(zeros(np_b + 9, n_tot));
[a, b] = meshgrid([-0.5 0 0.5], [-0.5 0 0.5]);
for t = 1:n_tot
    [bx, by] = rr_outline(sc_w(t)*(1-1e-6), sc_h(t)*(1-1e-6), rc(t), np_b);
    PXl(:,t) = [xt(t) + bx; xt(t) + a(:)*sc_w(t)];
    PYl(:,t) = [rt(t) + by; rt(t) + b(:)*sc_h(t)];
end
owner = repmat(1:n_tot, np_b + 9, 1);
[Bx, By] = wp_field_at_points(PXl(:), PYl(:), xt, rt, sc_w, sc_h, I_turn, p.n_TF, rc);
Bmag = sqrt(Bx.^2 + By.^2);
B_peak = accumarray(owner(:), Bmag, [n_tot 1], @max)';
B_center = Bmag(sub2ind(size(PXl), (np_b + 5)*ones(1, n_tot), 1:n_tot))';
B_center = B_center(:)';

% Smeared/linear model at each layer, as in scan_wp_designs.m
n_spire_ = zeros(1, n_layers);
n_spire_(1) = sum(n_turns);
for k = 2:n_layers
    n_spire_(k) = n_spire_(k-1) - n_turns(k);
end
B_layers = row.B_TF .* (n_spire_/n_spire_(1));
B_smooth = B_layers(layer_of_turn);

out.r_SC = rc;
out.x = xt;
out.y = rt;
out.layer = layer_of_turn;
out.B_smooth = B_smooth;
out.B_center = B_center;
out.B_discrete = B_peak;
out.ripple = B_peak ./ B_smooth;
end

function [x, y] = rr_outline(w, h, r, n)
% n points evenly spaced along the outline of a w x h rectangle with corner radius r
L = 2*(w - 2*r) + 2*(h - 2*r) + 2*pi*r;
s = (0:n-1)'/n*L;
x = zeros(n,1); y = zeros(n,1);
seg = [w-2*r, pi*r/2, h-2*r, pi*r/2, w-2*r, pi*r/2, h-2*r, pi*r/2];
c = [0 cumsum(seg)];
cx = [w/2-r, -w/2+r, -w/2+r, w/2-r]; cy = [h/2-r, h/2-r, -h/2+r, -h/2+r];
for i = 1:n
    k = find(s(i) >= c(1:end-1) & s(i) < c(2:end), 1); u = s(i) - c(k);
    switch k
        case 1, x(i) = w/2-r - u; y(i) = h/2;                          % top, right to left
        case 3, x(i) = -w/2; y(i) = h/2-r - u;                          % left, top to bottom
        case 5, x(i) = -w/2+r + u; y(i) = -h/2;                         % bottom
        case 7, x(i) = w/2; y(i) = -h/2+r + u;                          % right
        otherwise                                                        % fillets
            q = k/2; th = q*pi/2 + u/max(r, eps); iq = mod(q, 4) + 1;
            x(i) = cx(iq) + r*cos(th); y(i) = cy(iq) + r*sin(th);
    end
end
end
