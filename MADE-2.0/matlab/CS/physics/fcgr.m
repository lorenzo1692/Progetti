function [plasma_cycles] = fcgr(JT, Cond_w, S_hoop)
%FCGR Fatigue crack growth rate: number of plasma cycles to failure.
%
%   plasma_cycles = FCGR(JT, Cond_w, S_hoop) grows a surface flaw in the
%   jacket (JK2LB steel, Paris/Walker law) under the given hoop stress
%   S_hoop [MPa] until it reaches its critical size or the fracture
%   toughness KIC, and returns the resulting number of plasma cycles
%   (thousands, already divided by the safety factor SFN and by 2 for
%   ramp-up/ramp-down).
%
%   Inputs:
%     JT      - jacket thickness [m]
%     Cond_w  - conductor width [m]
%     S_hoop  - hoop stress array [MPa]; the max is used
%
%   Relocated unchanged from the CS legacy archive (fcgr.m). Calls
%   PARAMETRIY. Note: CS_DEMO_DesignExplorer_A45.m duplicated this logic
%   inline with 316LN constants instead of JK2LB - that duplicate was not
%   ported (see manuale CS, duplicazione #7).

%% Constants
type_flaw = 1;  % 0 = "embedded", 1 = "surface"
r = 3;          % Geometry factor (ITER standard)
SFa = 2;        % Safety factor for area
SFK = 1.5;      % Safety factor for stress intensity factor
SFN = 2;        % Safety factor for cycles
residual_stress = 200; % [MPa]
mw = 0.5;       % Walker exponent

C0 = 1.75e-13;  % JK2LB (ITER DDD CS p. 6-81), Paris law constant [m/cycle*(MPa*sqrt(m))^m]
m = 3.7;        % JK2LB, Paris law exponent

%% Inputs
Sa = max(S_hoop);  % [MPa]
w = Cond_w;        % [m]
t = JT;             % [m]
a_f = t * 0.8;      % Final crack depth [m]
c_f = w * 0.8;      % Final crack length [m]

%% Initial crack dimensions
if type_flaw == 0
    surf = 5e-6;    % [m^2]
    a_0 = sqrt(surf * SFa / pi / r);
elseif type_flaw == 1
    surf = 1e-6;    % [m^2]
    a_0 = sqrt(2 * surf * SFa / pi / r);
end
c_0 = a_0 * r;

%% Stress intensity factors and ratios
Smin = residual_stress;
Smax = residual_stress + Sa;
R = Smin / Smax;
KIC = 200 / SFK;   % [MPa*sqrt(m)]

n_coeff = -m * (mw - 1);
C = C0 * (1 - R) ^ (-n_coeff);

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
        break;
    end

    da_dN = C * Ke ^ m;
    dN = da / da_dN;
    dc = (c / a) * da;

    a = a + da;
    c = c + dc;
    Ncicli = Ncicli + dN;
end

cycles = Ncicli / SFN * 1e-3;
plasma_cycles = cycles / 2;
end
