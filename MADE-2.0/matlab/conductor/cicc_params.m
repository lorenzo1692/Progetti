function cp = cicc_params()
%CICC_PARAMS Strand/cable constants shared by CICC and its postprocessing.
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
cp.Tau_delay = 1;          % [s] quench detection + discharge delay (Tau_delay_mode 'fixed')
% Tau_delay_mode: 'fixed' uses cp.Tau_delay; 'detection' computes it for
% every copper amount with QUENCH_DELAY_TIME (t_det + t_hold + t_act).
% PROVISIONAL values below: to be agreed and referenced before use.
cp.Tau_delay_mode = 'fixed';
cp.quench.V_th = 0.1;      % [V] detection voltage threshold
cp.quench.t_hold = 1.0;    % [s] holding / validation time
cp.quench.t_act = 0.5;     % [s] dump-circuit actuation (breaker opening)
cp.quench.v_q_LTS = 5;     % [m/s] normal-zone propagation speed, LTS CICC
cp.quench.v_q_HTS = 0.05;  % [m/s] normal-zone propagation speed, HTS cable
cp.quench.t_det_HTS = [];  % [s] fixed HTS detection time (non-voltage detection); [] = from v_q_HTS
cp.quench.RRR_seg = 300;   % as HEAT_BALANCE_CICC_ODE
cp.quench.RRR_nonseg = 100;
cp.quench.T_ref_LTS = 6.8; % [K] resistivity temperature of the normal zone (= ODE start)
cp.quench.T_ref_HTS = 20;
cp.THS_max_LTS = 250;      % [K] default hot-spot limit, LTS
cp.THS_max_HTS = 150;      % [K] default hot-spot limit, HTS
end
