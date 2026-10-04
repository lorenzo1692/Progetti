function fp = fcgr_options(p)
%FCGR_OPTIONS Build the FCGR parameter struct from the machine input parameters.
%
%   fp = FCGR_OPTIONS(p) collects the fatigue crack growth parameters read
%   from the input workbook (category "Fatigue (FCGR)", see
%   PFC/input/WP_PFC_input_template.xlsx) into the struct FCGR expects.
%   Preset values for the two jacket steels used in the MADE domains:
%     316LN : C0 = 3.86e-11, m = 2.394, mw = 0.5   (PFC legacy fgcr.m)
%     JK2LB : C0 = 1.75e-13, m = 3.7,   mw = 0.5   (CS, ITER DDD CS p.6-81)

fp = struct();
fp.C0 = p.fcgr_C0;                          % Paris law constant [m/cycle*(MPa*sqrt(m))^m]
fp.m = p.fcgr_m;                            % Paris law exponent
fp.mw = p.fcgr_mw;                          % Walker exponent
fp.residual_stress = p.fcgr_residual_stress; % [MPa]
fp.KIC = p.fcgr_KIC;                        % [MPa*sqrt(m)]
fp.flaw_aspect = p.fcgr_flaw_aspect;        % c0/a0 (ITER standard: 3)
fp.flaw_type = p.fcgr_flaw_type;            % 0 = embedded, 1 = surface
fp.flaw_area = p.fcgr_flaw_area;            % [m^2] initial flaw area
end
