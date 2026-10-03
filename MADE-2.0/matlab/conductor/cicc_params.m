function cp = cicc_params(p)
%CICC_PARAMS Strand/cable constants shared by CICC and its postprocessing.
%   cp = CICC_PARAMS(p) also applies the overrides of the machine input p
%   (Tau_delay_mode, Tau_delay, v_quench_LTS).
%
%   cp = CICC_PARAMS() is the single place for the constants CICC uses to
%   size the cable and to run HEAT_BALANCE_CICC_ODE, so that
%   postprocess/plot_hotspot_transient.m re-runs the hot-spot transient
%   with exactly the same assumptions as the sizing did.

cp.CunonCu   = 1;          % Cu/non-Cu ratio of the SC strands
cp.S_tapes   = 4e-7;       % [m^2] REBCO tape cross-section
cp.d_fili    = 0.001;      % [m] strand diameter
cp.theta     = 30;         % [deg] REBCO field angle for Ic_sst33
cp.T_dim     = 20;         % [K] HTS design temperature
cp.d_cc      = 0.005;      % [m] cooling channel diameter
cp.cos_theta = 0.97;       % twist-pitch factor on the strand cross-section
cp.VF        = 0.7;        % cable void fraction
cp.N_fili_max = 1500;      % max number of SC strands/tapes per cable
cp.B_NbTi_max = 6;         % [T] LTS: NbTi below this field, Nb3Sn above
cp.B_hybrid_HTS = 15;      % [T] hybrid WP (WP_SC_type 102): HTS above this field
cp.Tau_delay = 1;          % [s] fixed quench delay before the dump (Tau_delay_mode 'fixed')
% Tau_delay_mode (input: Tau_delay_mode 0/1):
%   'fixed' - Tau_delay for every grade;
%   'iter'  - QUENCH_DELAY_TIME with the ITER TF protection assumptions:
%             dump 2 s after the coil voltage reaches 0.1 V (R. Zanino et
%             al., "Quench analysis of an ITER TF coil"); HTS grades keep
%             Tau_delay (no ITER data for HTS).
cp.Tau_delay_mode = 'fixed';
cp.quench.V_th = 0.1;      % [V] detection voltage (ITER TF)
cp.quench.t_hold = 2.0;    % [s] delay after detection: holding + dump actuation (ITER TF)
cp.quench.t_act = 0;       % [s] included in t_hold for the ITER assumptions
cp.quench.v_q_LTS = 5;     % [m/s] normal-zone propagation speed, LTS CICC (input v_quench_LTS)
cp.quench.v_q_HTS = 0.05;  % [m/s] HTS (study only: 'iter' mode uses Tau_delay for HTS)
cp.quench.t_det_HTS = [];  % [s] fixed HTS detection time; [] = from v_q_HTS
cp.quench.RRR_seg = 300;   % as HEAT_BALANCE_CICC_ODE
cp.quench.RRR_nonseg = 100;
cp.quench.T_ref_LTS = 6.8; % [K] resistivity temperature of the normal zone (= ODE start)
cp.quench.T_ref_HTS = 20;
cp.THS_max_LTS = 250;      % [K] default hot-spot limit, LTS
cp.THS_max_HTS = 150;      % [K] default hot-spot limit, HTS

% overrides from the machine input file (optional argument p)
if nargin >= 1 && ~isempty(p)
    if isfield(p, 'Tau_delay_mode') && ~isempty(p.Tau_delay_mode)
        modes = {'fixed', 'iter'};
        if ~any(p.Tau_delay_mode == [0 1])
            error('cicc_params:Tau_delay_mode', 'Tau_delay_mode must be 0 (fixed) or 1 (ITER), got %g.', p.Tau_delay_mode);
        end
        cp.Tau_delay_mode = modes{p.Tau_delay_mode + 1};
    end
    if isfield(p, 'Tau_delay') && ~isempty(p.Tau_delay), cp.Tau_delay = p.Tau_delay; end
    if isfield(p, 'v_quench_LTS') && ~isempty(p.v_quench_LTS), cp.quench.v_q_LTS = p.v_quench_LTS; end
end
end
