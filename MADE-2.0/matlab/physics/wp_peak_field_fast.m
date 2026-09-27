function [B_layer, info] = wp_peak_field_fast(row, p)
%WP_PEAK_FIELD_FAST Per-layer peak field on the conductor, fast enough for the scan.
%
%   [B_layer, info] = WP_PEAK_FIELD_FAST(row, p) returns, for every layer,
%   the peak field magnitude [T] over the cables of that layer, with the
%   same 2D discrete-coil model as COMPUTE_DISCRETE_FIELD_PROFILE (coil 1:
%   uniform current density on the real cable section, rounded or round,
%   closed-form Biot-Savart incl. self-field; other coils: line filaments;
%   validated against ANSYS within 0.1% on the peak BSUM), but evaluated
%   only where the peak of a layer can be:
%     - the edge turn of each layer and its neighbour, and the central
%       turn(s): the coil section is symmetric about its own axis, so one
%       edge suffices; the peak of a layer is at its edge turn (low field
%       layers, return field of the case side) or at its centre (plasma
%       side layer);
%     - 24 points on the outline of each of these cables (the peak of the
%       self-field is on the cable surface).
%   Checked against the full profile on the cases of
%   validation/check_peak_field_model.m.
%
%   It replaces, inside SCAN_WP_DESIGNS, the smeared model (Ampere field
%   of n_TF current sheets times the constant p.corr_B_WP, linear across
%   the WP), whose error on the peak grows with the toroidal narrowness of
%   the WP (1-2 T on the plasma side, more on the low-field grades).
%
%   row - struct (or table row) with Iop, n_layers, n_turns, Cond_w, Cond_h,
%         JT, Ri_ and shape_cable (as WP_TURN_GEOMETRY)
%   p   - machine parameters (n_TF, dr_plasma_side, GoundIns, INS_grades,
%         turn_insulation_nominal, Increm, r_SC_min, r_SC_max)
%
%   B_layer - 1 x n_layers peak field per layer [T]
%   info    - struct: B_peak (overall peak), layer_peak (its layer)

nl = row.n_layers;
n_turns = row.n_turns(1:nl);
Cond_w = row.Cond_w(1:nl); Cond_h = row.Cond_h(1:nl);
tg = wp_turn_geometry(row, p);

% turn centroids, same radial recursion as the scan and the section plot
Re = row.Ri_ - p.dr_plasma_side - p.GoundIns;
n_tot = sum(n_turns);
xt = zeros(1, n_tot); rt = zeros(1, n_tot); lay = zeros(1, n_tot); pos = zeros(1, n_tot);
i0 = 0;
for k = 1:nl
    idx = i0 + (1:n_turns(k));
    xt(idx) = -Cond_w(k)*n_turns(k)/2 + Cond_w(k)/2 + (0:n_turns(k)-1)*Cond_w(k);
    rt(idx) = Re - Cond_h(k)/2;
    lay(idx) = k; pos(idx) = 1:n_turns(k);
    Re = Re - Cond_h(k) - p.INS_grades;
    i0 = i0 + n_turns(k);
end
cw = tg.cab_w(lay); ch = tg.cab_h(lay); cr = tg.cab_r(lay);

% candidate turns: edge, next to edge, centre (one or two)
nt_l = n_turns(lay);
cand = pos == 1 | pos == 2 | pos == floor((nt_l+1)/2) | pos == floor(nt_l/2) + 1;
ic = find(cand);
np_b = 48;
PX = zeros(np_b, numel(ic)); PY = PX;
for q = 1:numel(ic)
    t = ic(q);
    [bx, by] = outline(cw(t)*(1-1e-6), ch(t)*(1-1e-6), cr(t)*(1-1e-6), np_b);
    PX(:,q) = xt(t) + bx; PY(:,q) = rt(t) + by;
end
[Bx, By] = field_near_far(PX(:), PY(:), xt, rt, cw, ch, cr, row.Iop, p.n_TF);
Bq = max(reshape(sqrt(Bx.^2 + By.^2), np_b, []), [], 1);
B_layer = accumarray(lay(ic)', Bq(:), [nl 1], @max)';
[info.B_peak, info.layer_peak] = max(B_layer);
end

function [Bx, By] = field_near_far(px, py, xc, yc, w, h, r, I, n_TF)
% Same field as WP_FIELD_AT_POINTS, vectorized: every cable of coil 1 acts
% on every point as ONE rectangle of equal area and aspect (exact for a
% sharp rectangle, a multipole error that decays fast with distance for a
% rounded or round one); for the points closer than NEAR cell sizes to a
% cable that rectangle is replaced by the exact decomposition of the
% rounded/round section (WP_FIELD_AT_POINTS's rectangles). Coils 2..n_TF:
% line filaments.
Mu_0 = 4e-7*pi; near = 2.5;
px = px(:); py = py(:); xc = xc(:)'; yc = yc(:)'; w = w(:)'; h = h(:)'; r = r(:)';
A = w.*h - (4 - pi)*r.^2;
we = w.*sqrt(A./(w.*h)); he = h.*sqrt(A./(w.*h));        % equivalent sharp rectangle
[Bx, By] = rect_field(px, py, xc, yc, we, he, I*ones(size(xc)));
% exact decomposition near each cable (correction = exact - equivalent)
[Rx, Ry, Rw, Rh, RI, Ro] = decompose(xc, yc, w, h, r, I);
dist = max(abs(px - xc)./w, abs(py - yc)./h);            % points x cables, in cell sizes
[ip, jc] = find(dist < near);
if ~isempty(ip)
    [cBx, cBy] = rect_field(px(ip), py(ip), xc(jc)', yc(jc)', we(jc)', he(jc)', -I*ones(numel(jc),1), true);
    % pairs (point, sub-rectangle of the same cable)
    [~, ~, g] = unique(Ro); nsub = accumarray(g(:), 1);
    first = cumsum([1; nsub(1:end-1)]);
    k = nsub(jc); rep = repelem((1:numel(ip))', k);
    off = (1:sum(k))' - repelem(cumsum([0; k(1:end-1)]), k);
    q = first(jc(rep)) + off - 1;
    [eBx, eBy] = rect_field(px(ip(rep)), py(ip(rep)), Rx(q), Ry(q), Rw(q), Rh(q), RI(q), true);
    cBx = cBx + accumarray(rep, eBx, [numel(ip) 1]); cBy = cBy + accumarray(rep, eBy, [numel(ip) 1]);
    Bx = Bx + accumarray(ip, cBx, [numel(px) 1]); By = By + accumarray(ip, cBy, [numel(px) 1]);
end
theta = 2*pi/n_TF;
for c = 2:n_TF
    a = (c-1)*theta;
    dx = px - (xc*cos(a) - yc*sin(a)); dy = py - (xc*sin(a) + yc*cos(a));
    d2 = dx.^2 + dy.^2;
    Bx = Bx + Mu_0*I/(2*pi)*sum(-dy./d2, 2);
    By = By + Mu_0*I/(2*pi)*sum( dx./d2, 2);
end
end

function [Bx, By] = rect_field(px, py, xc, yc, w, h, I, paired)
% field of uniform-current rectangles: all points x all rectangles (summed
% over rectangles), or element-wise pairs when paired is true
if nargin < 8 || ~paired
    xc = xc(:)'; yc = yc(:)'; w = w(:)'; h = h(:)'; I = I(:)';
else
    xc = xc(:); yc = yc(:); w = w(:); h = h(:); I = I(:);
end
k = 4e-7*pi*I./(w.*h)/(2*pi);
ua = px - (xc - w/2); ub = px - (xc + w/2);
va = py - (yc - h/2); vb = py - (yc + h/2);
By = k.*(G(ua,va) - G(ua,vb) - G(ub,va) + G(ub,vb));
Bx = -k.*(G(va,ua) - G(vb,ua) - G(va,ub) + G(vb,ub));
if nargin < 8 || ~paired
    Bx = sum(Bx, 2); By = sum(By, 2);
end
end

function [Rx, Ry, Rw, Rh, RI, Ro] = decompose(xc, yc, w, h, r, I)
% rounded / round cable as rectangles, as WP_FIELD_AT_POINTS (8 strips per
% quarter disk, area preserving); Ro = owner cable
ns = 8; R = zeros(0, 6);
for j = 1:numel(xc)
    rj = min(r(j), 0.5*min(w(j), h(j)));
    J = I/(w(j)*h(j) - (4 - pi)*rj^2);
    if rj <= 0
        R(end+1,:) = [xc(j), yc(j), w(j), h(j), I, j]; %#ok<AGROW>
        continue
    end
    R(end+1,:) = [xc(j), yc(j), w(j)-2*rj, h(j), J*(w(j)-2*rj)*h(j), j]; %#ok<AGROW>
    for sx = [-1 1]
        R(end+1,:) = [xc(j)+sx*(w(j)/2-rj/2), yc(j), rj, h(j)-2*rj, J*rj*(h(j)-2*rj), j]; %#ok<AGROW>
    end
    yk = linspace(0, rj, ns+1);
    Fk = @(y) 0.5*(y.*sqrt(rj^2 - y.^2) + rj^2*asin(min(y/rj, 1)));
    for k = 1:ns
        wk = (Fk(yk(k+1)) - Fk(yk(k)))/(yk(k+1) - yk(k));
        yb = (yk(k) + yk(k+1))/2; hb = yk(k+1) - yk(k);
        for sx = [-1 1]
            for sy = [-1 1]
                R(end+1,:) = [xc(j)+sx*(w(j)/2-rj+wk/2), yc(j)+sy*(h(j)/2-rj+yb), wk, hb, J*wk*hb, j]; %#ok<AGROW>
            end
        end
    end
end
R = R(R(:,3) > 0 & R(:,4) > 0, :);
Rx = R(:,1); Ry = R(:,2); Rw = R(:,3); Rh = R(:,4); RI = R(:,5); Ro = R(:,6);
end

function g = G(u, v)
% primitive of the line-current kernel: 1/2*v*ln(u^2+v^2) + u*atan(v/u)
t1 = 0.5*v.*log(u.^2 + v.^2); t1(v == 0) = 0;
t2 = u.*atan(v./u);           t2(u == 0) = 0;
g = t1 + t2;
end

function [x, y] = outline(w, h, r, n)
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
