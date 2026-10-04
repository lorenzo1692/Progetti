function [S_hoop, S_rad, S_ver, S_T] = eqv_stress_coil_ring_cicc(FZ, Re_grades, Ri_grades, Cond_h, Cond_w, JT, SC_h, SC_w, tins, type_cable, S_CICC, S_JT, var, Iop, Bmin, Bmax, WP_h, n_l, opts)
%EQV_STRESS_COIL_RING_CICC Thick-ring hoop/radial/vertical stress on a PFC jacket.
%
%   [S_hoop, S_rad, S_ver, S_T] = EQV_STRESS_COIL_RING_CICC(FZ, Re_grades,
%   Ri_grades, Cond_h, Cond_w, JT, SC_h, SC_w, tins, type_cable, S_CICC,
%   S_JT, var, Iop, Bmin, Bmax, WP_h, n_l) returns the equivalent hoop,
%   radial, vertical and Tresca (S_T=|S_hoop+S_ver|) stresses [Pa] on the
%   jacket of a Poloidal Field coil, from a thick-cylinder (ring) model.
%
%   Unlike MADE-2.0/matlab/CS/physics/eqv_stress_coil_cicc.m (a "coil"
%   formula assuming a single uniform field PB through the winding pack),
%   this ring formula distinguishes the field at the coil's inner and
%   outer radius (Bmin, Bmax) and switches model depending on the coil's
%   aspect ratio:
%     - WP_h <= WP_w (a "flat" ring, wider than it is tall) or a nonzero
%       Bmin: linear field gradient (Bmax at Ri, Bmin at Re) via the Kk/Mm
%       terms below. This matches a coil dominated by its own field (peak
%       Bz on the inner radius, return flux on the outer one); the formula
%       assumes it, it does not check where Bmax/Bmin actually are (swapping
%       them raises the hoop stress of a typical PF design by ~25%, see
%       manuale PFC, revisione modelli);
%     - otherwise (a "tall" ring, closer to a solenoid): the same
%       thin-walled, uniform-PB formula as EQV_STRESS_COIL_CICC.
%   FZ (the coil's axial/vertical force, [N]) is an external input here -
%   see PFC/search/scan_wp_designs.m for where it comes from (a
%   per-scenario external file, not the FZmax this module's own
%   EMAG_FIELD_FORCES would compute - see manuale PFC, "fidelity note:
%   axial force").
%
%   opts (optional struct) overrides the legacy hardcoded constants:
%     E_jckt, E_cbl_LTS, E_cbl_HTS, E_ins [GPa], ni (Poisson ratio) and
%     Fz_area_fix (1 = corrected, default; 0 = legacy, see below). Other
%     defaults are the legacy values (205, 0.1, 120, 20 GPa, 1/3).
%
%   VERTICAL STRESS: S_ver = FZ/(pi*(Re^2-Ri^2)) * K_rv, the axial force over
%   the annulus area. The legacy driver computed FZ/(Re^2-Ri^2)*pi (the
%   factor pi multiplying instead of dividing), which overestimates S_ver by
%   pi^2 = 9.87; it is still available with Fz_area_fix = 0 to reproduce
%   legacy results.
%
%   Relocated unchanged from the PF legacy archive
%   (eqv_stress_coil_ring_cicc.m, part of MADE_PF.7z), the formula
%   actually used by PF_opt_VNS.m (as opposed to
%   eqv_stress_coil_ring_lasso.m, used only by the separate, unported
%   "LASSO" winding variant - see manuale PFC).

%% Materials [GPa]
if nargin < 19 || isempty(opts)
    opts = struct();
end
E_jckt = get_opt(opts, 'E_jckt', 205);
E_cbl_HTS = get_opt(opts, 'E_cbl_HTS', 120);
E_cbl_LTS = get_opt(opts, 'E_cbl_LTS', 0.1);
E_ins = get_opt(opts, 'E_ins', 20);
ni = get_opt(opts, 'ni', 1/3);
Fz_area_fix = get_opt(opts, 'Fz_area_fix', 1);
J = Iop/(Cond_h*Cond_w);

%% Hoop stress
WP_w = Re_grades - Ri_grades;
a_ = Ri_grades;
b_ = Re_grades;
alpha = b_/a_;
e = linspace(a_, b_, n_l)./a_;
Kk = (alpha*Bmax - Bmin)*J*a_/(alpha-1);
Mm = (Bmax - Bmin)*J*a_/(alpha-1);
PB = Bmax^2/(2*4*pi*10^-7);

k_tot = Cond_h*Cond_w*E_jckt;
S_Cable = S_CICC - S_JT;
S_ins = (Cond_h*Cond_w) - S_CICC;

k_jkt_LTS = (S_JT*E_jckt + S_Cable*E_cbl_LTS + S_ins*E_ins);
k_jkt_HTS = (S_JT*E_jckt + S_Cable*E_cbl_HTS + S_ins*E_ins);

if WP_h <= WP_w || abs(Bmin) > 0
    if strcmp(type_cable{1, var}, 'LTS')
        S_hoop__ = Kk .* (2+ni)./(3*(alpha+1)) .* (alpha^2+alpha+1+alpha^2./e.^2 - e.*((1+2*ni)*(alpha+1)/(2+ni))) ...
                   - Mm .* (3+ni)/8 .* (alpha^2+1+alpha^2./e.^2 - (1+3*ni)/(3+ni)*e.^2);
        S_hoop = S_hoop__ * 1.1 * (k_tot/k_jkt_LTS);
    else
        S_hoop__ = Kk .* (2+ni)./(3*(alpha+1)) .* (alpha^2+alpha+1+alpha^2./e.^2 - e.*((1+2*ni)*(alpha+1)/(2+ni))) ...
                   - Mm .* (3+ni)/8 .* (alpha^2+1+alpha^2./e.^2 - (1+3*ni)/(3+ni)*e.^2);
        S_hoop = S_hoop__ * 1.1 * (k_tot/k_jkt_HTS);
    end
else
    c1 = 2/((alpha-1)^2);
    c2 = 7*alpha/(9*(alpha+1));
    c3 = alpha^2+alpha+1+(alpha^2./e.^2) - 5/7*(alpha+1)*e;
    c4 = (5/12)*(alpha^2+1+(alpha^2./e.^2) - 3/5*e.^2);
    if strcmp(type_cable{1, var}, 'LTS')
        S_hoop = c1*(c2*c3-c4)*max(PB)*1.0*(k_tot/k_jkt_LTS);
    else
        S_hoop = c1*(c2*c3-c4)*max(PB)*1.0*(k_tot/k_jkt_HTS);
    end
end

%% Radial stress
e = 1;
S_rad__ = Kk*(2+ni)/(3*(alpha+1))*(alpha^2+alpha+1+alpha^2/e^2 - e*(alpha+1)) ...
          - Mm*(3+ni)/8*(alpha^2+1+alpha^2/e^2 - e^2);
if strcmp(type_cable{1, var}, 'LTS')
    S_rad = S_rad__*(k_tot/k_jkt_LTS);
else
    S_rad = S_rad__*(k_tot/k_jkt_HTS);
end

%% Vertical stress
E_cbl = strcmp(type_cable{1, var}, 'LTS')*E_cbl_LTS + strcmp(type_cable{1, var}, 'HTS')*E_cbl_HTS;

Ke_cavo = (2*E_jckt*JT/(Cond_h-2*tins)) + (2*tins*E_ins/Cond_h) + ...
          ((1/(E_cbl*SC_w/SC_h)) + (2/(E_jckt*Cond_w/JT)) + (2/(E_ins*Cond_w/tins)))^-1;

K_jckt = 2*JT/(Cond_h-2*tins)*E_jckt;
dcr_jckt = K_jckt/Ke_cavo;
r_steel = (Cond_w-2*tins)/(2*JT);
K_rv = r_steel*dcr_jckt;

if Fz_area_fix
    S_ver = ones(size(S_hoop))*(FZ/((Re_grades^2 - Ri_grades^2)*pi))*max(K_rv);
else
    S_ver = ones(size(S_hoop))*(FZ/((Re_grades^2 - Ri_grades^2))*pi)*max(K_rv);
end

%% Tresca stress
S_T = abs(S_hoop + S_ver);
end

function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
