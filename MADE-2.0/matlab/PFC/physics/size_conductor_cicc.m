function [type_cable, N_Cu, N_Sc, N_tot, S_Cable, S_REBCO, S_Cu_HTS, THS, mat] = size_conductor_cicc(B_local, Iop, Tau_discharge, WP_SC_type)
%SIZE_CONDUCTOR_CICC Size one PFC grade's conductor (SC + Cu strands) against Ic and hot-spot limits.
%
%   [type_cable, N_Cu, N_Sc, N_tot, S_Cable, S_REBCO, S_Cu_HTS, THS, mat] =
%   SIZE_CONDUCTOR_CICC(B_local, Iop, Tau_discharge, WP_SC_type) picks the
%   superconductor (Nb3Sn / NbTi / REBCO, depending on WP_SC_type and
%   B_local) via the material Ic_* laws, sets N_Sc from Ic (rounded up to
%   the nearest of a fixed catalogue of strand counts, LTS only), then
%   bisects N_Cu (number of copper strands) against a hot-spot temperature
%   limit (HEAT_BALANCE_CICC_ODE) and returns the equivalent cable area
%   S_Cable.
%
%   WP_SC_type: 100 = full LTS, 101 = full HTS, 102 = Hybrid (HTS above
%   15 T, LTS below; LTS further splits Nb3Sn/NbTi at 6.2 T).
%
%   Relocated and renamed from the PF legacy archive (cicc.m, part of
%   MADE_PF.7z). Kept PF's own constants throughout (d_fili=0.82 mm,
%   T_dim=4.2 K, theta=0, the 15 T HTS/Hybrid threshold, the N_fili
%   catalogue) - these differ from every other domain's cicc-equivalent
%   (TF/CS use T_dim=20 K and no strand-count catalogue), so this is NOT a
%   copy of size_conductor_cicc.m from CS: it is PFC's own conductor law
%   (see manuale PFC).
%
%   FIDELITY NOTE (preserved as-is, see manuale PFC): for an LTS design
%   (mat 0 or 1), N_Sc is rounded UP to the nearest entry of N_fili, a
%   fixed catalogue of 25 standard strand counts (1500 down to 144). If
%   the required N_Sc exceeds every catalogue entry, N_fili>N_Sc is empty
%   and N_Sc becomes empty too, which SIZE_CICC_CABLE below will simply
%   propagate as an unusable candidate rather than raising a dedicated
%   error - exactly the legacy behavior.
%
%   Fix (same as CS, 01/10/2026): the legacy cicc.m called
%   Ic_sst33(T_dim, B_local, ...), i.e. with B and T swapped relative to
%   Ic_sst33's (B, T, theta, opt) signature. Here it is called as
%   Ic_sst33(B_local, T_dim, ...). Changes results only for
%   WP_SC_type=101 or the HTS branch of 102 (never for the LTS default).

CunonCu = 1;
S_tapes = 4*1e-7;    % [m^2] REBCO tape cross-section
d_fili = 0.00082;    % [m] strand diameter
d_cc = 0.005;        % [m] cooling-channel diameter
cos_theta = 0.97;    % twist correction on the strand bundle cross-section
VF = 0.8;            % void fraction in the conductor
N_fili = [1500 1440 1350 1296 1200 1152 1080 972 960 900 864 810 768 720 675 648 540 486 360 324 300 216 180 162 144]; % [-] catalogue of standard strand counts

theta = 0;
T_dim = 4.2; % [K] LTS design temperature (also used as the REBCO temperature in Ic_sst33)

N_Cu = 0;
type_cable = {'X'};
S_REBCO = 0;
S_Cu_HTS = 0;
S_Cable_0 = 0;

if WP_SC_type == 101
    type_cable = {'HTS'};
    Ic_sc = Ic_sst33(B_local, T_dim, theta, [3,4]);
    mat = 2; % 0 = Nb3Sn, 1 = NbTi, 2 = REBCO
    Tlim = 150;
elseif WP_SC_type == 100
    type_cable = {'LTS'};
    Ic_sc = Ic_Nb3Sn_PF(B_local, d_fili*1e3);
    mat = 0;
    Tlim = 250;
    if B_local < 6.2
        Ic_sc = Ic_NbTi_PF(B_local, d_fili*1e3);
        mat = 1;
        Tlim = 250;
    end
elseif WP_SC_type == 102
    if B_local > 15
        Ic_sc = Ic_sst33(B_local, T_dim, theta, [3,4]); %#ok<NASGU>
        mat = 2;
        Tlim = 150;
    else
        type_cable = {'LTS'};
        Ic_sc = Ic_Nb3Sn_PF(B_local, d_fili*1e3);
        mat = 0;
        Tlim = 250;
        if B_local < 6.2
            Ic_sc = Ic_NbTi_PF(B_local, d_fili*1e3);
            mat = 1;
            Tlim = 250;
        end
    end
else
    error('size_conductor_cicc:unknown_WP_SC_type', 'Unknown WP_SC_type code: %g', WP_SC_type);
end

N_Sc = ceil(Iop/Ic_sc); % Number of SC wires/tapes needed

if mat ~= 2
    above = N_fili > N_Sc;
    catalogue_hits = N_fili(above);
    N_Sc = min(catalogue_hits); % rounds up to the nearest standard strand count - see fidelity note above
end

if N_Sc > 1500
    N_tot = N_Sc; S_Cable = S_Cable_0; THS = NaN;
    return
end

%% Bisect N_Cu against the hot-spot temperature limit
div = 10;
N_Cu0 = linspace(1, 1000, div);
while true
    THS = zeros(size(N_Cu0,2), 1);
    for j = 1:size(N_Cu0,2)
        THS(j) = heat_balance_cicc_ode(N_Sc, N_Cu0(j), d_fili, CunonCu, Iop, B_local, Tau_discharge, mat, d_cc, VF, cos_theta, S_tapes);
    end
    [~, indx] = min(abs(THS-Tlim));
    if THS(indx) > Tlim
        b = round(N_Cu0(min(max(1,indx+1), div))); a = round(N_Cu0(indx));
    else
        a = round(N_Cu0(max(1,indx-1))); b = round(N_Cu0(indx));
    end
    if abs(Tlim-THS(indx)) <= 5 || b-a <= 1 || a > 980
        N_Cu = ceil(N_Cu0(indx));
        THS = THS(indx);
        break
    end
    N_Cu0 = linspace(a, b, div);
end

N_tot = N_Sc + N_Cu;

if mat == 2
    S_REBCO = N_Sc*S_tapes;
    S_Cu_HTS = N_Cu*pi*d_fili^2/4;
    S_Cable_0 = S_REBCO + S_Cu_HTS + (pi*d_cc^2/4);
elseif N_Sc < max(N_fili)
    S_Cable_0 = (N_tot*pi*d_fili^2/4/cos_theta/VF + (pi*d_cc^2/4));
end

D_eqv = (S_Cable_0*4/pi)^0.5;          % Equivalent cable diameter from area
A_w = ((D_eqv+0.0004)^2 - D_eqv^2)*pi/4; % Steel wrapping area (0.4 mm thick)
S_Cable = pi/4*D_eqv^2 + A_w;            % Corrected LTS cable equivalent area
end
