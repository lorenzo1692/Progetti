function [Ic, FA] = Ic_sst33(B, T, theta, opt)
%IC_SST33 REBCO (SuperOx SST) tape critical current fit with angular dependence.
%
%   [Ic, FA] = IC_SST33(B, T, theta, opt) with theta=0 perpendicular field,
%   theta=90 parallel field, opt=[B-T dataset, angular-law] (default [3,4]
%   = SuperOx'23 4mm-wide tape data, both B-T and angular).
%
%   Relocated unchanged from the PF legacy archive (Ic_sst33.m, part of
%   MADE_PF.7z). Byte-for-byte the same fit as MADE-2.0/matlab/CS/physics/
%   materials/Ic_sst33.m; kept as PFC's own copy per the 29/09/2026
%   decision to keep the CS and PFC physics trees independent for now (see
%   manuale PFC) - consolidate into a shared MADE-2.0/matlab/materials/
%   later if/when the two are merged.
%
%   CALLER NOTE: the legacy cicc.m called this as Ic_sst33(T_dim, B_local, ...),
%   swapping B and T. Fixed in size_conductor_cicc.m (as in CS): always call
%   as Ic_sst33(B, T, theta, opt).

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
