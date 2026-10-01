function [plasma_cycles] = fcgr(JT, Cond_w, S_hoop)
%FCGR Fatigue crack growth (Paris/Walker) life of the PFC jacket, in plasma cycles.
%
%   plasma_cycles = FCGR(JT, Cond_w, S_hoop) grows a semi-elliptical
%   surface flaw (ITER aspect ratio r=3) under the given jacket thickness
%   JT [m], conductor width Cond_w [m] and hoop stress S_hoop [MPa] using
%   the Paris law with Walker mean-stress correction, and returns the
%   number of plasma cycles to grow the flaw to 80% of the jacket
%   thickness/width (whichever governs).
%
%   Relocated and renamed from the PF legacy archive (fgcr.m -> FCGR,
%   part of MADE_PF.7z; the legacy filename was a typo of "fcgr", renamed
%   here for consistency with MADE-2.0/matlab/CS/physics/fcgr.m). Uses
%   316LN stainless steel Paris-law constants (C0=3.86e-11, m=2.394,
%   Walker exponent mw=0.5) - PFC's own jacket material, distinct from the
%   JK2LB constants used by CS/physics/fcgr.m (see manuale PFC).
%
%   NOT called by PF_opt_VNS.m (present in the legacy archive but
%   unexercised there); wired into PFC/search/scan_wp_designs.m anyway for
%   architectural consistency with the CS pipeline, per the 29/09/2026
%   decision recorded in manuale PFC.

%% Constants and initial parameters
type_flaw = 1;  % 0 == "embedded" - 1 == "surface"

r = 3;          % Geometry factor (ITER standard)
SFa = 1;        % Safety factor for area
SFK = 1;        % Safety factor for stress intensity factor
SFN = 1;        % Safety factor for cycles
residual_stress = 200; % Residual stress [MPa]

% Material properties for 316LN stainless steel
C0 = 3.86e-11;  % Paris law constant [m/cycle*(MPa*sqrt(m))^m]
m = 2.394;      % Paris law exponent
mw = 0.5;       % Walker exponent

%% Convert inputs
Sa = max(S_hoop); % Maximum hoop stress [MPa]
w = Cond_w(1);     % Conductor width [m]
t = JT(1);         % Jacket thickness [m]
a_f = t*0.8;       % Final crack depth [m]
c_f = w*0.8;       % Final crack length [m]

%% Initial crack dimensions
if type_flaw == 0
    surf = 5e-6; % [m^2]
    a_0 = sqrt(surf*SFa/pi/r);
elseif type_flaw == 1
    surf = 2e-6; % [m^2]
    a_0 = sqrt(2*surf*SFa/pi/r);
end
c_0 = a_0*r;

%% Stress intensity factors and ratios
Smin = residual_stress;
Smax = residual_stress + Sa;
R = Smin/Smax;
KIC = 200/SFK; % [MPa*sqrt(m)]

n_coeff = -m*(mw-1);
C = C0*(1-R)^(-n_coeff);

%% Crack growth loop
da = 1e-6;
a = a_0;
c = c_0;
Ncicli = 0;
Ke = 0;

while a < a_f && Ke < KIC && c < c_f
    [Ka, Kc] = ParametriY(a, c, w, t, Sa);
    Ke = min(Ka, Kc);

    if Ke >= KIC
        break
    end

    da_dN = C*Ke^m;
    dN = da/da_dN;
    dc = (c/a)*da;

    a = a + da;
    c = c + dc;
    Ncicli = Ncicli + dN;
end

cycles = Ncicli/SFN*1e-3; % [thousands of cycles]
plasma_cycles = cycles/2;
end
