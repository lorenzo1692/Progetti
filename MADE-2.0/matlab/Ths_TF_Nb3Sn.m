% HOT SPOT TEMPERATURE exercise (Wilson - Superconducing Magnets) (modified) for DEMO
% USE unities from SI (A,T,m2)
%% PRIMARY FUNCION:
function[Ths] = Ths_TF_Nb3Sn(Iop,Bp,Cus,Cuin,NB3SN,Tau_discharge)
%     V_MAX = 10000; % [V]
    RRR_seg = 300;
    RRR_inSC = 100;
    Tc = 6.8; % [K] 
    Tau_delay = 3;% [s] 
    
%     Tau_discharge_ = L*Iop/V_MAX*3;% [s]
%     if Tau_discharge_ > Tau_discharge
%         Tau_discharge = Tau_discharge_;
%     end
   
    t = (0:Tau_discharge/100:Tau_discharge); % [s] time 
    dt=t(2)-t(1);
    Cu_cross_section=Cus+Cuin; % total Cu cross-section [mm^2];
    %% Preallocations:
    Icu=zeros(1,max(size(t)));
    B=zeros(1,max(size(t)));
    energy=zeros(1,max(size(t)));
    Jcu=zeros(1,max(size(t)));
    f=zeros(1,max(size(t)));
    deltaT=zeros(1,max(size(t))+1);
    T_calc=zeros(1,max(size(t))+1);
    for i=1:max(size(t))
           if (t(i) <= Tau_delay)
                 Icu(i)=(Iop);
                 B(i)=Bp;
           else
                 Icu(i)=(Iop)*exp(-((t(i)-Tau_delay)/Tau_discharge));
                 B(i)=(Bp)*exp(-((t(i)-Tau_delay)/Tau_discharge)); % B(i)=(B_grades);%
                %Icu(i)=exp(-(t(i)/Tau_discharge(j)));
           end
           Jcu(i)=Icu(i)/Cu_cross_section; %  (/Atot);
    end
    %% temperature calculation algorithm
    deltaT(1)=0;
    T_calc(1)=Tc;
    for t_ind=1:max(size(t))
       energy(t_ind)=(Jcu(t_ind)^2)*rho_calculationB_parallel_RRR(t(t_ind),T_calc(t_ind),B(t_ind),RRR_inSC,RRR_seg,Cuin,Cus,Cu_cross_section)...
						*Cu_cross_section*dt; % energy generated during an interval dt

       f(t_ind)=((Cu_spec_heat_calculationB(T_calc(t_ind),Cu_cross_section))+...
           (Nb3Sn_spec_heat_calculationB(T_calc(t_ind),NB3SN)));

       % Sarebbe una dipendenza dalla Temperatura che trasformo in una dipendenza dal tempo (visto che T=T(t));
       deltaT(t_ind+1)=energy(t_ind)/f(t_ind);  % Temperature increment for the balance;
       T_calc(t_ind+1)=T_calc(t_ind)+deltaT(t_ind+1);
    end
    Ths=round(max(T_calc)); 
    % plot
	%     figure
	%     plot(t,T_calc(1:max(size(t))));
	%     ylabel('T_c_a_l_c [K]','fontSize',11)
	%     xlabel('t [s]','fontsize',12)
	%     title('HOT spot Temperature DEMO TF','fontSize',11)
	%     grid on
	%     keyt(1)={['Hot spot Temperature = ',num2str(round(max(T_calc))),' [K]']};
	%     keyt(2)={['Cu Area = ',num2str(round(Cu_cross_section*1e6)),' [mm^2]']};
	%     keyt(3)={['Nb_3Sn Area = ',num2str(round(NB3SN*1e6)),' [mm^2]']};
	%     %keyt(4)={['316LN steel Area = ',num2str(LNsteel_cross_section*1e6),' [mm^2]']};
	%     %keyt(5)={['He Area = ',num2str(He_cross_section*1e6),' [mm^2] Pressure =',num2str(pressure),' [bar]']};
	%     keyt(6)={['Tau_d_e_l = ',num2str(round(Tau_delay)),' [s] Tau_d_i_s = ',num2str(round(Tau_discharge)),' [s]']};
	%     %keyt(7)={['J_c = ',num2str(Jop/1000),' [kA] B_p = ',num2str(B_grades),' [T]']};
	%     text(30,120,keyt,'HorizontalAlignment','left','FontSize',12);
    %% SUBFUNCTIONS:
    function [rhoCu_parallel]=rho_calculationB_parallel_RRR(t,T,B,RRR_inSC,RRR_seg,Cu_cross_sect_in_SC,Cu_cross_sect_seg,Cu_cross_section)
        rhoCu_parallel=Cu_cross_section*(rho_calculationB(t,T,B,RRR_seg)*rho_calculationB(t,T,B,RRR_inSC))/...
        (rho_calculationB(t,T,B,RRR_inSC)*Cu_cross_sect_seg+rho_calculationB(t,T,B,RRR_seg)*Cu_cross_sect_in_SC);
    end
    %%
    function [rhoCu]=rho_calculationB(t,T,B,RRR)
        %preallocations:
        A=zeros(max(size(t)),max(size(T)));
        a=zeros(max(size(t)),max(size(T)));
        rhoCu=zeros(max(size(t)),max(size(T)));
        rho1=(1.171.*(10^-17).*(T.^4.49))./...
           (1+(4.5.*(10^-7).*(T.^3.35).*(exp(-(50./T).^6.428)))); % [ohm m]
        rho2=((1.69.*(10^-8)./RRR)+rho1+0.4531.*...
            ((1.69.*(10^-8).*rho1)./(RRR.*rho1+1.69.*(10^-8))));  % [ohm m]
        for j=1:max(size(T))
            for k=1:max(size(t)) %t
              A(k,j)=log10(1.553.*(10^-8).*B(k)./rho2(j));
              a(k,j)=-2.662+(0.3168.*A(k,j))+(0.6229.*(A(k,j).^2))-(0.1839.*(A(k,j).^3))+(0.01827.*(A(k,j).^4));
              rhoCu(k,j)=(rho2(j).*(1+(10.^a(k,j))));  % [ohm m]
            end
        end
    end
    %%
    function [Cp_Cu]=Cu_spec_heat_calculationB(T,Cu_cross_section)
        % specific heat Cu parameters from Dresner ("Stability of superconductors", 1995):
        density=8960; % [Kg/m^3]
        cp300=3.454e6; % [J/K/m^3] known data point at 300K
        gamma=0.011;  % [J/K^2/Kg]
        beta=0.0011;  % [J/K^4/Kg]
        c_plow=(beta.*(T.^3))+(gamma.*T); % [J/K-Kg] low temperature range
        Cp_Cu=(1./((1/cp300)+(1./(c_plow.*density))))*Cu_cross_section; % [J/K/m^3] volumetric specific heat for the whole temperature range
    end
    %%
    function [Cp_Nb3Sn]=Nb3Sn_spec_heat_calculationB(T,NB3SN)
        % Nb3Sn specific heat material parameters from Dataset2 ITER DRG1
        %preallocation:
        Cp_Nb3Sn=zeros(1,max(size(T)));
        gamma_Nb=0.1; % [J/K^2/Kg]
        beta_Nb=0.001; % [J/K^4/Kg]
        density_Nb=8040; % [Kg/m^3]
        Cp300_Nb=210; % [J/K/Kg]
        % Phenomenological fit based on the Debye model:
        Cp_low_NC=(beta_Nb.*(T.^3))+(gamma_Nb.*T); % [J/K/Kg] NORMAL
        % Cp_low_SC=((beta_Nb+(15*gamma_Nb/(Tc0_Nb^2))).*(T.^3))+(gamma_Nb*(B./Bc20_Nb).*T); % [J/K/Kg] SUPERCONDUCTING applies for T<Tc
        % I use Cp_low_NC
        for kk=1:max(size(T))
             %if T(i)<Tc
                 %Cp_Nb3Sn(i)=1./((1/Cp300_Nb)+(1./Cp_low_SC(i)));   
                 %else
                 Cp_Nb3Sn(kk)=1./((1/Cp300_Nb)+(1./Cp_low_NC(kk)));
                 %end
        end
        Cp_Nb3Sn=Cp_Nb3Sn*density_Nb*NB3SN;%;
    end
end