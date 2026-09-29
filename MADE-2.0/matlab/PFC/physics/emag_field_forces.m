function [Bsum, Bmin, Bmax, FRmax, FZmax, Br, Bz, Fr, Fz, Rc, Zc] = emag_field_forces(dy, dx, Ri_, n_t, n_l, n_g, Iop)
%EMAG_FIELD_FORCES Field map and radial/vertical Lorentz forces of one PFC grade.
%
%   [Bsum, Bmin, Bmax, FRmax, FZmax, Br, Bz, Fr, Fz, Rc, Zc] =
%   EMAG_FIELD_FORCES(dy, dx, Ri_, n_t, n_l, n_g, Iop) places every turn of
%   the coil on an n_t (axial) x n_l (radial) grid, all carrying Iop, and
%   evaluates the radial/vertical field (via XBR/XBZ, elliptic-integral
%   filament formulas) and the resulting Lorentz forces at one radial face
%   of each turn: the inner face for the inner half of the layers, the
%   outer face for the outer half (same convention as the field/hot-spot
%   check point of a single-layer solenoid). FZmax and FRmax are in MN.
%
%   Relocated from the PF legacy archive (emag.m, part of MADE_PF.7z),
%   same algorithm. Single-module, single-grade version (no n_moduli/
%   spacer stacking): matches how PF_opt_VNS.m always calls it, one grade
%   and one module at a time. Requires XBR.M and XBZ.M (elliptic-integral
%   filament field), shared with the TF/CS pipelines and already on the
%   path from MADE-2.0/matlab/.

n = 0;
data = zeros(n_t*n_l, 5);

for g = 1:n_g
    for i = 1:n_t
        for j = 1:n_l
            n = n + 1;
            data(n, 1) = dy(g)*(i-1);                       % Zc
            data(n, 2) = Ri_(g) + dx(g)/2 + dx(g)*(j-1);    % Rc
            data(n, 3) = Iop;
            if j <= n_l/2 + 1
                delta_r = -dx(g)/2;
            else
                delta_r = dx(g)/2;
            end
            data(n, 4) = Ri_(g) + dx(g)/2 + dx(g)*(j-1) + delta_r; % Rp
            data(n, 5) = dy(g)*(i-1);                                % Zp
        end
    end
end

Rp = data(:, 4); Zp = data(:, 5);
Ic = data(:, 3);

Rc = data(:, 2); Zc = data(:, 1);
Sc = 1e-32*ones(numel(data(:, 2)), 1);
Br = sum(xbr(Rp, Zp, Rc, Zc, Sc), 2).*Ic;
Bz = sum(xbz(Rp, Zp, Rc, Zc, Sc), 2).*Ic;
Fr = Bz.*(2*pi).*Rp.*Ic;
Fz = -Br.*(2*pi).*Rp.*Ic;

FZmax = sum(Fz(1:ceil(end/2)))*1e-6;
FRmax = sum(Fr)*1e-6;
Bsum = max(sqrt(Bz.^2 + Br.^2));
Bmin = min(Bz);
Bmax = max(Bz);
end
