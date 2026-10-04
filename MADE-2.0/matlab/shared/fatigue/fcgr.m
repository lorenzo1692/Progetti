function [plasma_cycles] = fcgr(JT, Cond_w, S_hoop, fp)
%FCGR Fatigue crack growth (Paris/Walker) life of a CICC jacket, in thousands of plasma cycles.
%
%   plasma_cycles = FCGR(JT, Cond_w, S_hoop, fp) grows a semi-elliptical
%   flaw under the given jacket thickness JT [m], conductor width Cond_w
%   [m] and hoop stress S_hoop [MPa] using the Paris law with Walker
%   mean-stress correction, and returns the number of plasma cycles [k]
%   to grow the flaw to 80% of the jacket thickness/width (whichever
%   governs).
%
%   fp (optional struct, see FCGR_OPTIONS) holds the material/flaw
%   parameters read from the input workbook; fields: C0, m, mw,
%   residual_stress, KIC, flaw_aspect, flaw_type, flaw_area and, optional,
%   the safety factors SFa (flaw area), SFK (KIC), SFN (cycles), default 1.
%   Shared by CS (JK2LB, ITER safety factors, see CS/search/scan_wp_designs)
%   and PFC (FCGR_OPTIONS, from the input workbook). Without fp the
%   legacy 316LN values are used (C0=3.86e-11, m=2.394, mw=0.5, residual
%   stress 200 MPa, KIC 200 MPa*sqrt(m), r=3, surface flaw of 2e-6 m^2).
%
%   Relocated and renamed from the PF legacy archive (fgcr.m -> FCGR,
%   part of MADE_PF.7z; the legacy filename was a typo of "fcgr"). The
%   hardcoded constants of the legacy function are now the defaults of fp.
%   PF_opt_VNS.m never called it; it is wired into PFC/search/
%   scan_wp_designs.m and PHYSICS/SIZE_CICC_CABLE.m as a sizing option
%   (input parameter fcgr_mode: 0 = off, 1 = compute only, 2 = jacket
%   sized until plasma_cycles >= plasma_cycles_min).

if nargin < 4 || isempty(fp)
    fp = struct('C0', 3.86e-11, 'm', 2.394, 'mw', 0.5, 'residual_stress', 200, ...
        'KIC', 200, 'flaw_aspect', 3, 'flaw_type', 1, 'flaw_area', 2e-6);
end

SFa = getf(fp, 'SFa', 1);  % Safety factor for area
SFK = getf(fp, 'SFK', 1);  % Safety factor for stress intensity factor
SFN = getf(fp, 'SFN', 1);  % Safety factor for cycles

r = fp.flaw_aspect;
C0 = fp.C0;
m = fp.m;
mw = fp.mw;
residual_stress = fp.residual_stress;

%% Convert inputs
Sa = max(S_hoop); % Maximum hoop stress [MPa]
w = Cond_w(1);     % Conductor width [m]
t = JT(1);         % Jacket thickness [m]
a_f = t*0.8;       % Final crack depth [m]
c_f = w*0.8;       % Final crack length [m]

%% Initial crack dimensions
if fp.flaw_type == 0
    a_0 = sqrt(fp.flaw_area*SFa/pi/r);     % embedded flaw
else
    a_0 = sqrt(2*fp.flaw_area*SFa/pi/r);   % surface flaw
end
c_0 = a_0*r;

%% Stress intensity factors and ratios
Smin = residual_stress;
Smax = residual_stress + Sa;
R = Smin/Smax;
KIC = fp.KIC/SFK;

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

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
