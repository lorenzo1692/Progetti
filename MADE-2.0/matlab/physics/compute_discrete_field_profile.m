function out = compute_discrete_field_profile(row, p)
%COMPUTE_DISCRETE_FIELD_PROFILE Exact 2D discrete-coil field at every turn.
%
%   out = COMPUTE_DISCRETE_FIELD_PROFILE(row, p) replaces the smooth/
%   smeared field model used during the scan (B_layers, a linear
%   interpolation between B_TF and 0 assuming n_TF continuous current
%   sheets via Ampere's law) with an exact 2D magnetostatic calculation:
%   every turn, in every one of the p.n_TF coils, is an infinite line
%   current perpendicular to this cross-section - the same idealization
%   the generalized-plane-strain FEM model uses (PLANE183, KEYOPT(3)=5).
%   The field at each turn is the vector sum of the Biot-Savart
%   contribution of every OTHER turn (its own filament is excluded, as
%   usual for the "external" field driving Ic degradation); this captures
%   both the toroidal ripple from having only p.n_TF discrete coils and
%   the local ripple from the WP's own discrete turns, neither of which
%   the smeared model sees.
%
%   Only meant to be run on ONE chosen design point (postprocessing/
%   verification), not inside the combinatorial scan: cost scales as
%   (n_TF * turns_per_coil)^2, still well under a second for a realistic
%   WP, but far too slow to repeat for every scanned candidate - the scan
%   keeps using the smeared model on purpose.
%
%   row - one row of a DATA table or an equivalent struct (see
%         PLOT_WP_SECTION for the field list); only Iop, n_layers,
%         n_turns, Cond_w, Cond_h, Ri_ are used here.
%   p   - machine parameter struct. Required: n_TF, dr_plasma_side,
%         GoundIns, INS_grades.
%
%   out fields (all for coil 1's turns, one entry per turn, ordered by
%   layer then by position within the layer):
%     x, y      - global Cartesian position of each turn [m]
%     layer     - which layer (1..n_layers) each turn belongs to
%     B_smooth  - the smeared/linear model's field at that turn's radius
%     B_discrete- the exact discrete Biot-Savart field magnitude [T]
%     ripple    - B_discrete./B_smooth (>1 means the smeared model under-predicts)

n_layers = row.n_layers;
n_turns  = row.n_turns(1:n_layers);
Cond_w   = row.Cond_w(1:n_layers);
Cond_h   = row.Cond_h(1:n_layers);

theta_TF = 2*pi/p.n_TF;
Mu_0 = 4e-7*pi;

% Per-layer radial positions, same recursion as search/scan_wp_designs.m
Re = zeros(1, n_layers+1);
Re(1) = row.Ri_ - p.dr_plasma_side - p.GoundIns;
Ri = zeros(1, n_layers);
for k = 1:n_layers
    Ri(k) = Re(k) - Cond_h(k);
    Re(k+1) = Ri(k) - p.INS_grades;
end

% Every turn's local toroidal offset (x_local, from coil 1's centerline)
% and radius (mid-height of its cell), in winding order.
x_local = [];
r_of_turn = [];
layer_of_turn = [];
for k = 1:n_layers
    y0 = Ri(k) + Cond_h(k)/2;
    w_total = Cond_w(k)*n_turns(k);
    x0 = -w_total/2 + Cond_w(k)/2;
    xs = x0 + (0:n_turns(k)-1)*Cond_w(k);
    x_local = [x_local, xs]; %#ok<AGROW>
    r_of_turn = [r_of_turn, repmat(y0, 1, n_turns(k))]; %#ok<AGROW>
    layer_of_turn = [layer_of_turn, repmat(k, 1, n_turns(k))]; %#ok<AGROW>
end
n_turns_tot = numel(x_local);

% Global position of every turn, replicated at each of the n_TF coils.
% x_local is an arc length at radius r_of_turn, so the angular offset
% from the coil's own centerline is x_local/r_of_turn (consistent with
% treating the WP cross-section as flat/local, as the rest of this solver
% already does for CASE_w etc.).
theta_local = x_local ./ r_of_turn;
X_all = zeros(p.n_TF, n_turns_tot);
Y_all = zeros(p.n_TF, n_turns_tot);
for c = 1:p.n_TF
    theta_c = (c-1)*theta_TF + theta_local;
    X_all(c,:) = r_of_turn .* cos(theta_c);
    Y_all(c,:) = r_of_turn .* sin(theta_c);
end

I_turn = row.Iop; % same current in every turn, every coil (all wind the same way)
Xs = X_all(:); Ys = Y_all(:);

Xf = X_all(1,:); Yf = Y_all(1,:);
Bx = zeros(1, n_turns_tot);
By = zeros(1, n_turns_tot);
for i = 1:n_turns_tot
    dx = Xf(i) - Xs;
    dy = Yf(i) - Ys;
    d2 = dx.^2 + dy.^2;
    d2(d2 < 1e-12) = Inf; % exclude this turn's own singular self-term
    Bx(i) = sum(Mu_0*I_turn/(2*pi) .* (-dy)./d2);
    By(i) = sum(Mu_0*I_turn/(2*pi) .* ( dx)./d2);
end
B_discrete = sqrt(Bx.^2 + By.^2);

% Smeared/linear model's field at each turn's radius, for comparison -
% same formula as the "Recompute B per grade" section of scan_wp_designs.m
n_spire_ = zeros(1, n_layers);
n_spire_(1) = sum(n_turns);
for k = 2:n_layers
    n_spire_(k) = n_spire_(k-1) - n_turns(k);
end
B_layers = row.B_TF .* (n_spire_/n_spire_(1));
B_smooth = B_layers(layer_of_turn);

out.x = Xf;
out.y = Yf;
out.layer = layer_of_turn;
out.B_smooth = B_smooth;
out.B_discrete = B_discrete;
out.ripple = B_discrete ./ B_smooth;
end
