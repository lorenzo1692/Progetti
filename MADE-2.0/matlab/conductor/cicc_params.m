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
cp.Tau_delay = 1;          % [s] quench detection + discharge delay
cp.THS_max_LTS = 250;      % [K] default hot-spot limit, LTS
cp.THS_max_HTS = 150;      % [K] default hot-spot limit, HTS
end
