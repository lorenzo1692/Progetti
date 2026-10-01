function [Bsum, Bmin, Bmax, FRmax, FZmax, Br, Bz, Fr, Fz, Rc, Zc] = emag_field_forces(Cond_h, Cond_w, Ri_, n_turns_, n_layers_, n_grades, n_moduli, Iop, spacer, full_stack)
%EMAG_FIELD_FORCES Field map and radial/axial forces of the CS module stack.
%
%   [Bsum, Bmin, Bmax, FRmax, FZmax, Br, Bz, Fr, Fz, Rc, Zc] =
%   EMAG_FIELD_FORCES(Cond_h, Cond_w, Ri_, n_turns_, n_layers_, n_grades,
%   n_moduli, Iop, spacer, full_stack) places every turn of the coil (all
%   modules carrying the same Iop) and evaluates the radial/axial field
%   (via XBR/XBZ, elliptic-integral filament formulas) and the resulting
%   radial/axial Lorentz forces. FZmax and FRmax are in MN (sum of Fz/Fr
%   in N, scaled by 1e-6 - see EQV_STRESS_COIL_CICC for why this matters).
%   full_stack=0 evaluates only the middle module (module 3) for speed;
%   full_stack=1 evaluates every module. FZmax is the sum of the Fz of the
%   first half of the evaluated turns: with full_stack=1 it is the axial
%   compression through the stack mid-plane (use this for the vertical
%   stress), with full_stack=0 it is only the force on half of one module.
%
%   Relocated from the CS legacy archive (emag.m), same algorithm. The
%   original accepted an extra `input` (scenario currents) argument that
%   the active code path never used (I(1:n_moduli)=Iop unconditionally;
%   the per-scenario line was already commented out) - dropped here, see
%   manuale CS.
%
%   Requires XBR.M and XBZ.M (elliptic-integral filament field), shared
%   with the TF pipeline and already on the path from MADE-2.0/matlab/.

I(1:n_moduli) = Iop;

n = 0;
data = zeros(sum(n_turns_ .* n_layers_)*n_moduli, 4);
h_m = Cond_h(1)*n_turns_(1)+spacer;

for m = 1:n_moduli
    for g = 1:n_grades
        for i = 1:n_turns_(g)
            for j = 1:n_layers_(g)
                n = n + 1;
                data(n, 1) = Cond_h(g) * (i - 1)-h_m*n_moduli+h_m*(m-1);    % Zc
                data(n, 2) = Ri_(g) + Cond_w(g) / 2 + Cond_w(g) * (j - 1);  % Rc
                data(n, 3) = I(m);
                data(n, 4) = Ri_(g) + Cond_w(g) / 2 + Cond_w(g) * (j - 1) - Cond_w(g) / 2;  % Rp
            end
        end
    end
end

mod_plot = 3;
switch full_stack
    case 0
        range = 1+n/n_moduli*(mod_plot-1):n/n_moduli*(mod_plot);
    case 1
        range = 1:n;
end
Rp = data(range, 4);    Zp = data(range, 1);
Ic = data(range, 3);

Rc = data(:, 2);    Zc = data(:, 1);
Sc = 1e-32*ones(numel(data(:, 2)), 1);
Br = sum(xbr(Rp, Zp, Rc, Zc, Sc),2).*Ic;  Bz = sum(xbz(Rp, Zp, Rc, Zc, Sc),2).*Ic;
Fr = Bz .* (2*pi) .* Rp .* Ic;     Fz = -Br .* (2*pi) .* Rp .* Ic;
FZmax = sum(Fz(1:ceil(end/2))) * 1e-6;
FRmax = sum(Fr) * 1e-6;
Bsum  = max(sqrt(Bz.^2 + Br.^2));
Bmin  = min(Bz);
Bmax  = max(Bz);
end
