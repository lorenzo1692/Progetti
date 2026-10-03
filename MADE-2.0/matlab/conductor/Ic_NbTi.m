function[Ic_A] = Ic_NbTi(B,d_fili)
%     %% fit Bottura
%     c = 1.05;
%     Bc20_T = 14.61;
%     Tc0_K = 9.03;
%     C0 = 1;
%     alpha = 1;
%     beta = 1.54;
%     gamma = 2.1;
%     n = 1.7;
%     T = 6.2; % operating temperature    
%     %% fit functions
%     t = T/Tc0_K;
%     b = B/Bc20_T;
%     Ic_A = c*C0/B*b^alpha*(1-b)^beta*(1-t^n)^gamma;
% end
    %% Fit 2 componenti
    d = d_fili;
    CunonCu = 1.9;
    % B = 7;
    T = 6.2;
    n = 1.7;
    Bc20_T = 15.19;
    Tc0_K = 8.907;
    t = T/Tc0_K;
    tt = 1-t^n;
    b = B/Bc20_T;
    C0	= 3.00E+04;
    C1	= 0.45;
    a1	= 3.2;
    b1	= 2.43;
    C2	= 0.55;
    a2	= 0.65;
    b2	= 2;
    g1	= 1.8;
    g2	= 1.8;
    G = (a1/(a1+b1))^a1;
    GG = (b1/(a1+b1))^b1;
    GGG = G*GG;
    F = (a2/(a2+b2))^a2;
    FF = (b2/(a2+b2))^b2;
    FFF = F*FF;    
    Jc = C0*C1/(B*GGG)*(b/tt)^a1*(1-b/tt)^b1*tt^g1+C0*C2/(B*FFF)*(b/tt)^a2*(1-b/tt)^b2*tt^g2;    
    Ic_A = Jc*pi*d^2/(4*(1+CunonCu));
end