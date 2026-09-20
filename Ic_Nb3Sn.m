function [Ic_A,Jc_sc] = Ic_Nb3Sn(B,d_fili,type)
% Ic_Nb3Sn
% Calcola Ic e Jc_sc per diversi strand Nb3Sn.
%
% INPUT:
%   B       = campo magnetico [T]
%   d_fili  = diametro strand [m]
%   type    = tipo di strand:
%             'DTT_TF_KAT'
%             'ITER'
%             'WST'
%
% OUTPUT:
%   Ic_A    = corrente critica [A]
%   Jc_sc   = densità di corrente critica sulla sola area non-Cu [A/m2]

    if nargin < 3
        type = 'DTT_TF_KAT';
    end

    %% Common constants
    strand_d = d_fili;
    CunonCu = 1;
    strand_A = pi*strand_d^2/(4*(1+CunonCu)); % superconducting area in the strand

    c_ = 1.0;

    %% Select strand parameters
    switch upper(type)

        case 'DTT_TF_KAT'

            Ca1 = 36.04;        % Deviatoric strain
            Ca2 = 0.00;         % Deviatoric strain
            eps_0a = 0.0022;    % Hydrostatic strain
            eps_m = 0;          % Thermal pre-strain

            Bc20max = 27.6;     % Maximum upper critical field [T]
            Tc0max = 16.0;      % Maximum critical temperature [K]
            C = 25304;          % Pre-constant [AT]

            p = 0.50;
            q = 1.75;

            int_eps = -0.55/100; % intrinsic strain

        case 'ITER'

            Ca1 = 36.35;        % Deviatoric strain
            Ca2 = 0.00;         % Deviatoric strain
            eps_0a = 0.21/100;  % Hydrostatic strain
            eps_m = -0.33/100;  % Thermal pre-strain

            Bc20max = 29.70;    % Maximum upper critical field [T]
            Tc0max = 15.35;     % Maximum critical temperature [K]
            C = 28300;          % Pre-constant [AT]

            p = 0.5;
            q = 2.0;

            int_eps = -0.55/100; % intrinsic strain

        case 'WST'

            Ca1 = 50.06;        % Deviatoric strain
            Ca2 = 0.00;         % Deviatoric strain
            eps_0a = 0.00312;   % Hydrostatic strain
            eps_m = -0.00059;   % Thermal pre-strain

            Bc20max = 33.24;    % Maximum upper critical field [T]
            Tc0max = 16.34;     % Maximum critical temperature [K]
            C = 83075*strand_A; % Pre-constant [AT]

            p = 0.593;
            q = 2.156;

            int_eps = -0.36/100; % intrinsic strain

        otherwise
            error('Unknown Nb3Sn strand type: %s', type);

    end

    %% Input operating temperature
    T = 4.2 + 1.5; % [K]

    %% Scaling-law fit
    eps_sh = Ca2*eps_0a/(sqrt(Ca1^2 - Ca2^2));

    s_eps = 1 + ...
        (Ca1*(sqrt(eps_sh^2 + eps_0a^2) - ...
        sqrt((int_eps - eps_sh)^2 + eps_0a^2)) - Ca2*int_eps) ...
        /(1 - Ca1*eps_0a);

    Bc0_eps = Bc20max*s_eps;
    Tc0_eps = Tc0max*(s_eps)^(1/3);

    t = T/Tc0_eps;

    BcT_eps = Bc0_eps*(1 - t^(1.52));

    b = B./BcT_eps;

    hT = (1 - t^(1.52))*(1 - t^2);

    fPb = b.^p .* (1 - b).^q;

    %% Critical values
    Ic_A = c_ .* (C./B) .* s_eps .* fPb .* hT;

    Jc_sc = Ic_A ./ strand_A;

end