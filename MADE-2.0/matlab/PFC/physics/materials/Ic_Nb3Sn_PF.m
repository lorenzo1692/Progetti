function [Ic_A] = Ic_Nb3Sn_PF(B, strand_d)
%IC_NB3SN_PF Critical current of the Nb3Sn strand (Bottura strain-scaling fit).
%
%   Ic_A = IC_NB3SN_PF(B, strand_d) returns the strand critical current [A]
%   at fixed operating temperature (6.2 K) and applied strain (-0.55%) for
%   a strand of diameter strand_d [mm] in a field B [T].
%
%   Relocated unchanged from the PF legacy archive (Ic_Nb3Sn_PF.m, part of
%   MADE_PF.7z). Same fit family as MADE-2.0/matlab/CS/physics/materials/
%   Ic_Nb3Sn_CS_DTT.m, with PF's own constants (Bc20max=32.97 T,
%   Tc0max=16.06 K, C=20918.1 AT) - kept as a separate file since the
%   constants differ from every other domain's Nb3Sn fit (see manuale
%   PFC). Dropped two locals (Jc_sc, Je_strand) that the legacy function
%   computed but never returned.

CunonCu = 1;
strand_A = pi*strand_d^2/(4*(1+CunonCu)); %#ok<NASGU> % strand SC area, kept for parity with the legacy source
c_ = 1.08;
Ca1 = 44.48;
Ca2 = 0.00;
eps_0a = 0.00256;
eps_m = -0.00049;
Bc20max = 32.97; % [T]
Tc0max = 16.06;  % [K]
C = 20918.1;     % [AT]
p = 0.63;
q = 2.1;

T = 6.2;      % [K] operating temperature
eps = -0.55;  % [%] applied strain
int_eps = eps/100 + eps_m;

eps_sh = Ca2*eps_0a/sqrt(Ca1^2 - Ca2^2);
s_eps = 1 + (Ca1*(sqrt(eps_sh^2+eps_0a^2) - sqrt((int_eps-eps_sh)^2+eps_0a^2)) - Ca2*int_eps)/(1 - Ca1*eps_0a);
Bc0_eps = Bc20max*s_eps;
Tc0_eps = Tc0max*s_eps^(1/3);

t = T/Tc0_eps;
BcT_eps = Bc0_eps*(1 - t^1.52);

b = B./BcT_eps;
hT = (1 - t^1.52)*(1 - t^2);
fPb = b.^p.*(1-b).^q;

Ic_A = c_*(C./B).*s_eps.*fPb.*hT;
end
