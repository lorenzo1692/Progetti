function [Ic_A] = Ic_NbTi_PF(B, strand_d)
%IC_NBTI_PF Critical current of the NbTi strand (2-component fit).
%
%   Ic_A = IC_NBTI_PF(B, strand_d) returns the strand critical current [A]
%   at fixed operating temperature (5.72 K) for a strand of diameter
%   strand_d [mm] in a field B [T].
%
%   Relocated unchanged from the PF legacy archive (Ic_NbTi_PF.m, part of
%   MADE_PF.7z). Same fit family, constants and structure as
%   MADE-2.0/matlab/CS/physics/materials/Ic_NbTi_CS.m, but at a different
%   operating temperature (PF: 4.22+1.5 = 5.72 K, CS: 5.2 K) - kept as a
%   separate file for that reason (see manuale PFC). The legacy file also
%   carried a commented-out single-component Bottura fit, dropped here as
%   dead code (never active).

CunonCu = 1.9;
T = 4.22 + 1.5; % [K]
n = 1.7;
Bc20_T = 15.19;
Tc0_K = 8.907;
t = T/Tc0_K;
tt = 1 - t^n;
b = B/Bc20_T;

C0 = 3.00e4;
C1 = 0.45; a1 = 3.2; b1 = 2.43;
C2 = 0.55; a2 = 0.65; b2 = 2;
g1 = 1.8; g2 = 1.8;

G = (a1/(a1+b1))^a1;
GG = (b1/(a1+b1))^b1;
GGG = G*GG;
F = (a2/(a2+b2))^a2;
FF = (b2/(a2+b2))^b2;
FFF = F*FF;

Jc = C0*C1/(B*GGG)*(b/tt)^a1*(1-b/tt)^b1*tt^g1 + C0*C2/(B*FFF)*(b/tt)^a2*(1-b/tt)^b2*tt^g2;
Ic_A = Jc*pi*strand_d^2/(4*(1+CunonCu));
end
