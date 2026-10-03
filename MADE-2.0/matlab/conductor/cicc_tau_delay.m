function t_delay = cicc_tau_delay(cp, Iop, B, mat, N_Sc, N_Cu)
%CICC_TAU_DELAY Tau_delay for the hot-spot transient (cicc_params.m: Tau_delay_mode).
t_delay = cp.Tau_delay;
if any(strcmpi(cp.Tau_delay_mode, {'iter', 'detection'})) && mat ~= 2
    A_seg = N_Cu*pi*cp.d_fili^2/4;
    A_non = N_Sc*pi*cp.d_fili^2/(4*(1 + cp.CunonCu));
    t_delay = quench_delay_time(Iop, A_seg, A_non, B, mat, cp.quench);
elseif strcmpi(cp.Tau_delay_mode, 'detection')
    t_delay = quench_delay_time(Iop, N_Cu*pi*cp.d_fili^2/4, 0, B, mat, cp.quench);
end
end
