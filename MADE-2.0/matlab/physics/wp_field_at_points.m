function [Bx, By] = wp_field_at_points(px, py, xc, yc, w, h, I, n_TF)
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
%
%   Coil 1: every cable is a uniform-current rectangle (closed-form
%   Biot-Savart, so its own self-field is included). Coils 2..n_TF: the
%   same cables rotated by (c-1)*2*pi/n_TF, as line filaments (they are far
%   away). Validated against ANSYS on two benchmarks (peak BSUM within
%   0.1%, Lorentz force per turn within 0.03%).

Mu_0 = 4e-7*pi;
sz = size(px);
px = px(:); py = py(:);
xc = xc(:)'; yc = yc(:)'; w = w(:)'; h = h(:)';
Bx = zeros(size(px)); By = zeros(size(px));

for j = 1:numel(xc)
    J = I/(w(j)*h(j));
    k = Mu_0*J/(2*pi);
    ua = px - (xc(j) - w(j)/2); ub = px - (xc(j) + w(j)/2);
    va = py - (yc(j) - h(j)/2); vb = py - (yc(j) + h(j)/2);
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
