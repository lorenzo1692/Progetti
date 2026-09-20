function[Jc,Ic_REBCO] = Jc_REBCO(B,S_tapes)
% B = 13;
% S_tapes = 0.004*100e-6; % mm2 
T = 20;
Tc = 92.83; % K
Birr0 = 120; % T B ortogonale
C = 12510; % A T
p = 0.5;
q = 1.7;
a = 1.52;
b = 2.33;
Birr = Birr0*(1-T/Tc)^a;
Ic_REBCO = C/B*(Birr/Birr0)^b*(B/Birr)^p*(1-B/Birr)^q;
Jc = Ic_REBCO/S_tapes;
end

% JB_REBCO = (0.9000*(B)^2-54.100*(B)+1.1455e3); 


%%
