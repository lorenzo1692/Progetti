function [S_hoop, S_rad, S_ver, S_T] = eqv_stress_coil_cicc(Fz, Ri_grades, Re_grades, Cond_h, Cond_w, JT, SC_h, SC_w, tins, type_cable, S_CICC, S_JT, PB)
%EQV_STRESS_COIL_CICC Thick-cylinder (Lame) hoop/radial/vertical stress on the CICC jacket.
%
%   [S_hoop, S_rad, S_ver, S_T] = EQV_STRESS_COIL_CICC(...) returns the
%   equivalent hoop, radial, vertical and Tresca (S_T=S_hoop+S_ver)
%   stresses [MPa] on the innermost jacket layer of a CICC-wound solenoid
%   module, from a thick-cylinder model plus stress-concentration factors.
%
%   UNITS: Fz must be passed in MN (not N) - since areas are in m^2,
%   Fz[MN]/Area[m^2] already comes out in MPa (1 MN/m^2 = 1 MPa), matching
%   the legacy driver's convention (e.g. emag_field_forces.m's FZmax is in
%   MN). Passing Fz in N here silently gives stresses 1e6 too high.
%
%   Relocated unchanged from the CS legacy archive (eqv_stress_coil_cicc.m).

%% Materials [GPa]
E_jckt = 205;
E_cbl_HTS = 120;
E_cbl_LTS = 0.1;
E_ins = 20;

%% Hoop
a_ = Ri_grades;
b_ = Re_grades;
a = b_/a_;
e = linspace(Ri_grades, Re_grades, 10)./a_;
c1 = 2/((a-1)^2);
c2 = 7*a/(9*(a+1));
c3 = a^2+a+1+(a^2./e.^2)-5/7.*(a+1).*e;
c4 = (5/12)*(a^2+1+(a^2./e.^2)-3/5.*e.^2);

k_tot = (Cond_h*Cond_w)*(E_jckt); % as if fully steel, consistent with S_hoop's reference stiffness
S_Cable = S_CICC-S_JT;
S_ins = (Cond_h*Cond_w)-S_CICC;
k_jkt_LTS = (S_JT*(E_jckt)+S_Cable*(E_cbl_LTS)+S_ins*(E_ins));
k_jkt_HTS = (S_JT*(E_jckt)+S_Cable*(E_cbl_HTS)+S_ins*(E_ins));

if strcmp(type_cable, 'LTS')
    S_hoop = c1*(c2*c3-c4)*(PB)*(k_tot/k_jkt_LTS)*1.1e-6;
else
    S_hoop = c1*(c2*c3-c4)*(PB)*(k_tot/k_jkt_HTS)*1.1e-6;
end

%% Radial
c1 = 2/((a-1)^2);
c2 = 7*a/(9*(a+1));
c3 = a^2+a+1-(a^2./e.^2)-(a+1).*e;
c4 = (5/12)*(a^2+1-(a^2./e.^2)-e.^2);

if strcmp(type_cable, 'LTS')
    S_rad = c1*(c2*c3-c4)*max(PB)*(k_tot/k_jkt_LTS)*1e-6;
else
    S_rad = c1*(c2*c3-c4)*max(PB)*1e-6*(k_tot/k_jkt_HTS);
end

%% Vertical
E_cbl = E_cbl_LTS;
Ke_cavo = 2*E_jckt*JT/Cond_h+2*tins*E_ins/Cond_h+(1/(E_cbl*SC_w/SC_h)+2/(E_jckt*Cond_w/JT)+2/(E_ins*Cond_w/tins))^-1;
K_jckt = 2*(JT/Cond_h*E_jckt);
dcr_jckt = K_jckt./Ke_cavo;
r_steel = (Cond_w-2.*tins)./(2.*JT);
K_rvLTS = r_steel.*dcr_jckt;

E_cbl = E_cbl_HTS;
Ke_cavo = 2*E_jckt*JT/Cond_h+2*tins*E_ins/Cond_h+(1/(E_cbl*SC_w/SC_h)+2/(E_jckt*Cond_w/JT)+2/(E_ins*Cond_w/tins))^-1;
K_jckt = 2*JT/Cond_h*E_jckt;
dcr_jckt = K_jckt/Ke_cavo;
r_steel = (Cond_w-2*tins)/(2*JT);
K_rvHTS = r_steel*dcr_jckt;

if strcmp(type_cable, 'LTS')
    S_ver = ones(size(S_hoop)).*Fz./((Re_grades^2-Ri_grades^2)*pi)*max(K_rvLTS);
else
    S_ver = ones(size(S_hoop)).*Fz./((Re_grades^2-Ri_grades^2)*pi)*max(K_rvHTS);
end

%% Tresca
S_T = S_hoop+S_ver;
end
