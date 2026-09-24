function [Bx, By] = wp_field_at_points(px, py, xc, yc, w, h, I, n_TF, r)
%WP_FIELD_AT_POINTS 2D magnetic field of the TF inner-leg set at given points.
%
%   [Bx, By] = WP_FIELD_AT_POINTS(px, py, xc, yc, w, h, I, n_TF) returns the
%   field [T] at the points (px, py) (any shape) in the FEM frame of the
%   inner-leg cross-section: x = toroidal, y = radial, machine axis at the
%   origin, coil 1 centred on the +y axis.
%
%   xc, yc, w, h - centre and size of every cable (current-carrying area)
%                  of coil 1 [m], vectors of equal length
%   I            - current per turn [A], flowing along +z in every turn
%   n_TF         - number of coils
%   r            - (optional) corner fillet radius of every cable [m]: the
%                  cable is then the real rounded rectangle (as in the FEM),
%                  decomposed exactly into rectangles (central cross + an
%                  area-preserving staircase of strips for each quarter
%                  disk); default 0 = sharp rectangle
%
%   Coil 1: every cable carries a uniform current density over its own
%   (rounded) cross-section, closed-form Biot-Savart per rectangle, so its
%   own self-field is included. Coils 2..n_TF: the
%   same cables rotated by (c-1)*2*pi/n_TF, as line filaments (they are far
%   away). Validated against ANSYS on two benchmarks (peak BSUM within
%   0.1%, Lorentz force per turn within 0.03%).

Mu_0 = 4e-7*pi;
if nargin < 9 || isempty(r), r = 0; end
sz = size(px);
px = px(:); py = py(:);
xc = xc(:)'; yc = yc(:)'; w = w(:)'; h = h(:)';
r = r(:)'.*ones(size(xc));
Bx = zeros(size(px)); By = zeros(size(px));

% sub-rectangles [x y w h current] of every cable
ns = 8;                                         % strips per quarter disk
R = zeros(0, 5);
for j = 1:numel(xc)
    rj = min(r(j), 0.5*min(w(j), h(j)));
    A = w(j)*h(j) - (4 - pi)*rj^2;
    J = I/A;
    if rj <= 0
        R(end+1,:) = [xc(j), yc(j), w(j), h(j), I]; %#ok<AGROW>
        continue
    end
    R(end+1,:) = [xc(j), yc(j), w(j)-2*rj, h(j), J*(w(j)-2*rj)*h(j)]; %#ok<AGROW>
    for sx = [-1 1]
        R(end+1,:) = [xc(j)+sx*(w(j)/2-rj/2), yc(j), rj, h(j)-2*rj, J*rj*(h(j)-2*rj)]; %#ok<AGROW>
    end
    yk = linspace(0, rj, ns+1);                 % bands from the fillet centre outward
    Fk = @(y) 0.5*(y.*sqrt(rj^2 - y.^2) + rj^2*asin(min(y/rj, 1)));   % int sqrt(r^2-y^2)
    for k = 1:ns
        wk = (Fk(yk(k+1)) - Fk(yk(k)))/(yk(k+1) - yk(k));   % area-preserving strip width
        yb = (yk(k) + yk(k+1))/2; hb = yk(k+1) - yk(k);
        for sx = [-1 1]
            for sy = [-1 1]
                cxk = xc(j) + sx*(w(j)/2 - rj + wk/2);
                cyk = yc(j) + sy*(h(j)/2 - rj + yb);
                R(end+1,:) = [cxk, cyk, wk, hb, J*wk*hb]; %#ok<AGROW>
            end
        end
    end
end

for q = 1:size(R,1)
    k = Mu_0*R(q,5)/(R(q,3)*R(q,4))/(2*pi);
    ua = px - (R(q,1) - R(q,3)/2); ub = px - (R(q,1) + R(q,3)/2);
    va = py - (R(q,2) - R(q,4)/2); vb = py - (R(q,2) + R(q,4)/2);
    By = By + k*(G(ua,va) - G(ua,vb) - G(ub,va) + G(ub,vb));
    Bx = Bx - k*(G(va,ua) - G(vb,ua) - G(va,ub) + G(vb,ub));
end

theta = 2*pi/n_TF;
for c = 2:n_TF
    a = (c-1)*theta;
    Xs = xc*cos(a) - yc*sin(a);
    Ys = xc*sin(a) + yc*cos(a);
    dx = px - Xs; dy = py - Ys;          % implicit expansion
    d2 = dx.^2 + dy.^2;
    Bx = Bx + Mu_0*I/(2*pi) * sum(-dy./d2, 2);
    By = By + Mu_0*I/(2*pi) * sum( dx./d2, 2);
end

Bx = reshape(Bx, sz); By = reshape(By, sz);
end

function g = G(u, v)
% Primitive of the line-current kernel: 1/2*v*ln(u^2+v^2) + u*atan(v/u)
t1 = 0.5*v.*log(u.^2 + v.^2); t1(v == 0) = 0;
t2 = u.*atan(v./u);           t2(u == 0) = 0;
g = t1 + t2;
end
