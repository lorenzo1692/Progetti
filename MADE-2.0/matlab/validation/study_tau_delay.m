%% study_tau_delay.m
% Effect of the quench delay criterion (cicc_params.m: Tau_delay_mode) on
% the copper sized by CICC, for one LTS and one HTS grade. Compares fixed
% delays with the detection-based delay of QUENCH_DELAY_TIME for a range of
% normal-zone propagation speeds. Prints a table; no file is written.
%
% Run from MADE-2.0/matlab:  study_tau_delay

this_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(this_dir); made_paths('TF');

Iop = 60e3; tau_dis = 20;   % [A], [s] illustrative grade
cases = {
    % label, B [T], WP_SC_type, mode, Tau_delay or v_q
    'LTS fixed 1 s',        12, 100, 'fixed', 1
    'LTS fixed 3 s',        12, 100, 'fixed', 3
    'LTS det. v_q 1 m/s',   12, 100, 'detection', 1
    'LTS det. v_q 5 m/s',   12, 100, 'detection', 5
    'LTS det. v_q 20 m/s',  12, 100, 'detection', 20
    'HTS fixed 1 s',        18, 101, 'fixed', 1
    'HTS det. v_q 0.2 m/s', 18, 101, 'detection', 0.2
    'HTS det. v_q 0.05 m/s',18, 101, 'detection', 0.05
    'HTS det. v_q 0.01 m/s',18, 101, 'detection', 0.01
    };
fprintf('%-24s %6s %6s %6s %8s %8s %8s\n', 'case', 'N_Sc', 'N_Cu', 'THS', 't_det', 't_delay', 'Cu mm2');
for i = 1:size(cases, 1)
    [lab, B, typ, mode, val] = cases{i, :};
    cp = cicc_params(); cp.Tau_delay_mode = mode;
    if strcmp(mode, 'fixed'), cp.Tau_delay = val;
    elseif typ == 101, cp.quench.v_q_HTS = val;
    else, cp.quench.v_q_LTS = val; end
    [~, N_Cu, N_Sc, ~, ~, ~, ~, THS, mat] = cicc(B, Iop, tau_dis, typ, [], [], cp);
    A_seg = N_Cu*pi*cp.d_fili^2/4;
    A_non = (mat ~= 2)*N_Sc*pi*cp.d_fili^2/(4*(1 + cp.CunonCu));
    if strcmp(mode, 'fixed')
        t_det = NaN; t_del = cp.Tau_delay;
    else
        [t_del, d] = quench_delay_time(Iop, A_seg, A_non, B, mat, cp.quench); t_det = d.t_det;
    end
    fprintf('%-24s %6d %6d %6g %8.2f %8.2f %8.0f\n', lab, N_Sc, N_Cu, THS, t_det, t_del, (A_seg + A_non)*1e6);
end
