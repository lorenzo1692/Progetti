function [t_delay, d] = quench_delay_time(Iop, A_Cu_seg, A_Cu_nonseg, B, mat, qp)
%QUENCH_DELAY_TIME Delay between quench onset and start of the dump (proposal).
%
%   [t_delay, d] = QUENCH_DELAY_TIME(Iop, A_Cu_seg, A_Cu_nonseg, B, mat, qp)
%   returns the time during which the hot spot is heated at full current,
%   i.e. the Tau_delay of HEAT_BALANCE_CICC_ODE, as the sum of
%
%     t_det    time for the normal zone to build up the detection voltage,
%     t_hold   holding (validation) time of the detection system,
%     t_act    actuation time of the dump circuit (breaker opening).
%
%   Detection: the normal zone grows from the quench point with two fronts
%   at speed v_q, so its length is 2*v_q*t and its voltage
%       V(t) = 2*v_q*t * Iop * r',   r' = rho_Cu(T_ref, B) / A_Cu
%   (current fully in the copper, resistivity at the low temperature T_ref:
%   the normal zone is hotter, so this is conservative, i.e. a longer t_det).
%   V = V_th gives
%       t_det = V_th * A_Cu / (2 * v_q * Iop * rho_Cu(T_ref, B)).
%   t_det grows with the copper: a cable with more stabilizer is detected
%   later. Called inside the copper search of CICC, the delay is therefore
%   consistent with the copper it protects.
%
%   mat: 0 Nb3Sn, 1 NbTi (LTS, qp.v_q_LTS), 2 REBCO (HTS, qp.v_q_HTS).
%   qp fields (see cicc_params.m): V_th [V], t_hold [s], t_act [s],
%   v_q_LTS, v_q_HTS [m/s], RRR_seg, RRR_nonseg, T_ref_LTS, T_ref_HTS [K],
%   t_det_HTS [s] (optional: if not empty, fixed detection time for HTS,
%   e.g. non-voltage detection such as optical fibres).
%
%   d: struct with t_det, t_hold, t_act, v_q, r_prime [ohm/m], V_per_m [V/m].

if mat == 2
    v_q = qp.v_q_HTS; T_ref = qp.T_ref_HTS;
else
    v_q = qp.v_q_LTS; T_ref = qp.T_ref_LTS;
end

rho_seg = rho_cu(T_ref, B, qp.RRR_seg);
rho_non = rho_cu(T_ref, B, qp.RRR_nonseg);
% parallel of segregated and non-segregated copper, per unit length
G = A_Cu_seg/rho_seg + A_Cu_nonseg/rho_non;          % [m/ohm]
r_prime = 1/G;                                        % [ohm/m]

if mat == 2 && isfield(qp, 't_det_HTS') && ~isempty(qp.t_det_HTS)
    t_det = qp.t_det_HTS;
else
    t_det = qp.V_th/(2*v_q*Iop*r_prime);
end
t_delay = t_det + qp.t_hold + qp.t_act;

d = struct('t_det', t_det, 't_hold', qp.t_hold, 't_act', qp.t_act, 'v_q', v_q, ...
    'r_prime', r_prime, 'V_per_m', Iop*r_prime);
end

function rho = rho_cu(T, B, RRR)
% Copper resistivity with magnetoresistance: same fit as HEAT_BALANCE_CICC_ODE
rho1 = (1.171e-17*T^4.49)/(1 + 4.5e-7*T^3.35*exp(-(50/T)^6.428));
rho2 = 1.69e-8/RRR + rho1 + 0.4531*(1.69e-8*rho1)/(RRR*rho1 + 1.69e-8);
A = log10(1.553e-8*B/rho2);
a = -2.662 + 0.3168*A + 0.6229*A^2 - 0.1839*A^3 + 0.01827*A^4;
rho = rho2*(1 + 10^a);
end
