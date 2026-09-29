function [Ic, FA] = Ic_sst33(B, T, theta, opt)
%IC_SST33 REBCO (SuperOx SST) tape critical current fit with angular dependence.
%
%   [Ic, FA] = IC_SST33(B, T, theta, opt) with theta=0 perpendicular field,
%   theta=90 parallel field, opt=[B-T dataset, angular-law] (default [3,4]
%   = SuperOx'23 4mm-wide tape data, both B-T and angular).
%
%   Relocated unchanged from the CS legacy archive (Ic_sst33.m).
%
%   CALLER NOTE (open question, see manuale CS - decisione aperta): in the
%   legacy cicc.m the call was Ic_sst33(T_dim, B_local, theta, [3,4]) with
%   T_dim=20 and B_local the actual field - i.e. B and T are swapped
%   relative to this signature. size_conductor_cicc.m (the CS/physics port
%   of cicc.m) reproduces that same call exactly, unchanged, pending
%   confirmation - it is NOT exercised by the WP_SC_type=100 (full LTS)
%   default path used by CS_opt_VNS.m.

if nargin < 4
    opt = [3, 4];
end

switch opt(1) % B-T dependence
    case 1 % first approach
        FF = jc_ybco(B, T, 'sst2')*3.3e-9;
    case 2 % more accurate
        FF = jc_ybco(B, T, 'sst')/jc_ybco(0, 77, 'sst')*43*3.3;
    case 3 % SO'23 4 mm wide
        FF = jc_ybco(B, T, 'sst')/jc_ybco(0, 77, 'sst')*43*4*1.92;
    otherwise
        error('Ic_sst33:unknown_BT_option', 'Unknown B-T dependence option: %g', opt(1));
end

switch opt(2) % angular dependence
    case 1 % interp 4.2 K curve
        data = load('SP_angular_4.2_DU.mat');
        FA = interp1(data.angle, data.alpha, theta/180*pi);
    case 2 % analytical
        FA = 1 + 3.222*exp((theta-90)/8.234) + 1.278*exp((theta-90)/0.5901);
    case 3 % interp var T curves (B influence neglected)
        data = load('SST_angular_low_T_Wimbush+DU.mat');
        FA = griddata(data.angle, data.temp, data.Icnorm, theta, T);
    case 4 % SO'23 analytical
        FA = 1 + 4.482*exp((theta-90)/10.03);
    otherwise
        error('Ic_sst33:unknown_angular_option', 'Unknown angular dependence option: %g', opt(2));
end

Ic = FF.*FA;
end
