function[Ic_A,Jc_sc] = Ic_Nb3Sn_WST(B,d_fili)
    %% constants
    strand_d = d_fili;
    CunonCu = 1;
    strand_A = pi*strand_d^2/(4*(1+CunonCu)); % area di superconduttore nello strand 
    c_ = 1.0;
    Ca1 = 50.06;        % Deviatoric strain
    Ca2 = 0.00;         % Deviatoric strain
    eps_0a = 0.00312;   % Hydrostatic strain
    eps_m = -0.00059;   % Thermal pre-strain
    Bc20max = 33.24;    % Maximum upper critical feld [T]
    Tc0max = 16.34;     % Maximum critical temperature[K]
    C = 83075*strand_A; % Pre-constant [AT]
    p = 0.593;
    q = 2.156;
    %% input
    T = 4.2+1.5; % temp
    % B = 14; % field    
    int_eps = -0.36/100; % intrinsic strain __ R&W -0.36 __ W&R -0.55
    % app_eps = int_eps-eps_m; % applied strain 
    %% fit functions
    eps_sh = Ca2*eps_0a/(sqrt(Ca1^2-Ca2^2));
    s_eps = 1+(Ca1*(sqrt(eps_sh^2+eps_0a^2)-sqrt((int_eps-eps_sh)^2+eps_0a^2))-Ca2*int_eps)/(1-Ca1*eps_0a);
    Bc0_eps = Bc20max*s_eps;
    Tc0_eps = Tc0max*(s_eps)^(1/3);
    t = T/Tc0_eps;
    BcT_eps = Bc0_eps*(1-t^(1.52));
    TcB_eps = Tc0max*(s_eps)^(1/3)*(1-B./Bc0_eps).^(1/1.52);
    b = B./BcT_eps;
    hT = (1-t^(1.52))*(1-t^2);
    fPb = b^p*(1-b)^q;
    %% critical values 
    Ic_A = c_*(C./B).*s_eps*fPb*hT;
    Jc_sc = Ic_A./strand_A;
    Je_strand = Ic_A./(pi*strand_d^2/4);
end