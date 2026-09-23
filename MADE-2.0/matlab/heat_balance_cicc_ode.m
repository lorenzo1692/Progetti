function [THS, t, TF] = heat_balance_cicc_ode(N_Sc,N_Cu,d_fili,CunonCu,Iop0,B0,Tau_discharge,mat,d_cc,VF,costheta,S_tapes,Tau_delay)
%HEAT_BALANCE_CICC_ODE Adiabatic hot-spot transient of a CICC after a quench.
%
%   THS = HEAT_BALANCE_CICC_ODE(...) returns the peak hot-spot temperature [K].
%   [THS, t, TF] = HEAT_BALANCE_CICC_ODE(...) also returns the full
%   temperature history TF(t), used by postprocess/plot_hotspot_transient.m.
%   Tau_delay [s] is optional (default 1 s, see cicc_params.m).

    if nargin < 13 || isempty(Tau_delay)
        Tau_delay = 1;
    end
    tspan = [0, (-Tau_discharge*log(1/Iop0)+Tau_delay)];

    % Initial temperature -> Tc
    if mat == 2 
        T0 = 20.0;
    else
        T0 = 6.8;
    end

    % ODE solver
    options = odeset('RelTol', 1e-6);
    [t, TF] = ode45(@temperatureODE, tspan, T0, options);

    % Results
    THS = round(max(TF));

    %  if THS <=250
    %     figure
    %     plot(t,TF(1:max(size(t))));
    %     ylabel('T_c_a_l_c [K]','fontSize',10)
    %     xlabel('t [s]','fontsize',10)
    %     title('HOT spot Temperature','fontSize',10)
    %     grid on
    %     keyt(1)={['Hot spot Temperature = ',num2str(round(max(TF))),' [K]']};
    %     keyt(2)={['Cu Area non seg = ',num2str((N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu)))*1e6),' [mm^2]']};
    %     keyt(3)={['Cu Area seg = ',num2str((N_Cu * pi * d_fili^2 / 4)*1e6),' [mm^2]']};
    %     keyt(4)={['Nb_3Sn Area = ',num2str(N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu))*1e6),' [mm^2]']};
    %     % keyt(4)={['316LN steel Area = ',num2str(LNsteel_cross_section*1e6),' [mm^2]']};
    %     % keyt(5)={['He Area = ',num2str(He_cross_section*1e6),' [mm^2] Pressure =',num2str(pressure),' [bar]']};
    %     keyt(5)={['Nsc = ',num2str(N_Sc),' NCu = ',num2str(N_Cu)]};
    %     keyt(6)={['Tau_d_e_l = ',num2str(Tau_delay),' [s] Tau_d_i_s = ',num2str(Tau_discharge),' [s]']};
    %     keyt(7)={['Iop = ',num2str(Iop0/1000),' [kA] Bop = ',num2str(B0),' [T]']};
    % 
    %     text(t(10),(max(TF)/1.55),keyt,'HorizontalAlignment','left','FontSize',10);
    % end

    function dTdt = temperatureODE(t, T)
        % ODE for temperature evolution

        % Calculate parameters

        if mat == 2
            SC_cross_sect = N_Sc * S_tapes;
            Cu_cross_sect_non_seg = 0;
        else
            SC_cross_sect = N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu));
            Cu_cross_sect_non_seg = N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu));
        end
        
        Cu_cross_sect_seg = N_Cu * pi * d_fili^2 / 4;
        Cu_cross_section = Cu_cross_sect_seg + Cu_cross_sect_non_seg;
        S_cc = d_cc^2 / 4 * pi;

        % Time-dependent current and magnetic field
        if t <= Tau_delay
            Iop = Iop0;
            B = B0;
        else
            Iop = Iop0 * exp(-((t - Tau_delay) / Tau_discharge));
            B = B0 * exp(-((t - Tau_delay) / Tau_discharge));
        end

        % Calculate resistances
        [rho_Nb3Sn, rho_NbTi, rho_Steel, rhoCu_non_seg, rhoCu_seg] = calculateResistances(T, B);

        % Equivalent resistance
        rho_Cu_parallel = (Cu_cross_section * rhoCu_non_seg * rhoCu_seg) / (rhoCu_seg * Cu_cross_sect_non_seg + rhoCu_non_seg * Cu_cross_sect_seg);
        Cu_res = rho_Cu_parallel / Cu_cross_section;
        Nb3Sn_res = rho_Nb3Sn / SC_cross_sect;
        NbTi_res = rho_NbTi / SC_cross_sect;
        Steel_res = rho_Steel / SC_cross_sect;
        cable_res_Nb3Sn = (1 / Cu_res + 1 / Nb3Sn_res)^(-1);
        cable_res_NbTi = (1 / Cu_res + 1 / NbTi_res)^(-1);
        cable_res_Steel = (1 / Cu_res + 1 / Steel_res)^(-1);
        % rho_cable_Nb3Sn = cable_res_Nb3Sn * ((Cu_cross_section + SC_cross_sect) / VF / costheta + S_cc);
        % rho_cable_NbTi = cable_res_NbTi * ((Cu_cross_section + SC_cross_sect) / VF / costheta + S_cc);
        % rho_cable_Steel = cable_res_Steel * ((Cu_cross_section + SC_cross_sect) / VF / costheta + S_cc);
        rho_cable_Nb3Sn = cable_res_Nb3Sn * ((Cu_cross_section + SC_cross_sect));
        rho_cable_NbTi = cable_res_NbTi * ((Cu_cross_section + SC_cross_sect));
        rho_cable_Steel = cable_res_Steel * ((Cu_cross_section + SC_cross_sect));

        % Specific heat
        [Cp_Cu, Cp_v_Nb3Sn, Cp_v_NbTi, Cp_v_Steel] = calculateSpecificHeat(T, SC_cross_sect, Cu_cross_section);

        % Select material parameters
        if mat == 0 
            Cp_v = Cp_v_Nb3Sn;
            rho_cable = rho_cable_Nb3Sn;
        elseif mat == 1 
            Cp_v = Cp_v_NbTi; 
            rho_cable = rho_cable_NbTi;
        elseif mat == 2 
            Cp_v = Cp_v_Steel;
            rho_cable = rho_cable_Steel;   
        end

        % Heat generation and absorption
        % Q_gen = (Iop / ((Cu_cross_section + SC_cross_sect) / VF / costheta + S_cc))^2 * rho_cable; 
        Q_gen = (Iop / ((Cu_cross_section + SC_cross_sect)))^2 * rho_cable; 
        Q_abs = Cp_v;
        deltaT = Q_gen / Q_abs;
        
        % Output
        dTdt = deltaT;
    end

    function [rho_Nb3Sn, rho_NbTi, rho_Steel, rhoCu_non_seg, rhoCu_seg] = calculateResistances(T, B)
        % Calculate resistances for Nb3Sn, NbTi, Cu non-segregated, and Cu segregated
        rho_Nb3Sn = (-1e-4 * T^2 + 0.0938 * T + 22.601) * 1e-8; % [ohm m]
        rho_NbTi = (0.0558 * T + 55.668) * 1e-8; % [ohm m]
        rho_Steel = (76.2063 + 0.071375*(T-273) - 2.3109e-5*(T-273)^2) * 1e-8;  % [ohm m]

        RRR = 300;
        rho1 = (1.171 * (10^-17) * (T^4.49)) / (1 + (4.5 * (10^-7) * (T^3.35) * (exp(-(50 / T)^6.428))));
        rho2 = ((1.69 * (10^-8) / RRR) + rho1 + 0.4531 * ((1.69 * (10^-8) * rho1) / (RRR * rho1 + 1.69 * (10^-8))));
        A = log10(1.553 * (10^-8) * B / rho2);
        a = -2.662 + (0.3168 * A) + (0.6229 * (A^2)) - (0.1839 * (A^3)) + (0.01827 * (A^4));
        rhoCu_seg = (rho2 * (1 + (10^a)));  % [ohm m]

        RRR = 100;
        rho1 = (1.171 * (10^-17) * (T^4.49)) / (1 + (4.5 * (10^-7) * (T^3.35) * (exp(-(50 / T)^6.428))));
        rho2 = ((1.69 * (10^-8) / RRR) + rho1 + 0.4531 * ((1.69 * (10^-8) * rho1) / (RRR * rho1 + 1.69 * (10^-8))));
        A = log10(1.553 * (10^-8) * B / rho2);
        a = -2.662 + (0.3168 * A) + (0.6229 * (A^2)) - (0.1839 * (A^3)) + (0.01827 * (A^4));
        rhoCu_non_seg = (rho2 * (1 + (10^a)));  % [ohm m]
    end

    function [Cp_Cu, Cp_v_Nb3Sn, Cp_v_NbTi, Cp_v_Steel] = calculateSpecificHeat(T, SC_cross_sect, Cu_cross_section)
        % Calculate specific heat for Cu, Nb3Sn, and NbTi
        density = 8960;  % [Kg/m^3]
        cp300 = 3.454e6; % [J/K/m^3] known data point at 300K
        gamma = 0.011;   % [J/K^2/Kg]
        beta = 0.0011;   % [J/K^4/Kg]
        c_plow = (beta * (T^3)) + (gamma * T); % [J/K-Kg] low temperature range
        Cp_Cu = (1 / ((1 / cp300) + (1 / (c_plow * density)))); % [J/K/m^3] volumetric specific heat for the whole temperature range        

        gamma_Nb = 0.1; % [J/K^2/Kg]
        beta_Nb = 0.001; % [J/K^4/Kg]
        density_Nb = 8040; % [Kg/m^3]
        Cp300_Nb = 210; % [J/K/Kg]
        Cp_low_NC = (beta_Nb * (T^3)) + (gamma_Nb * T); % [J/K/Kg] NORMAL
        Cp_Nb3Sn = 1 / ((1 / Cp300_Nb) + (1 / Cp_low_NC));
        Cp_Nb3Sn = Cp_Nb3Sn * density_Nb;

        gamma_Nb = 0.145; % [J/K^2/Kg]
        beta_Nb = 0.0023; % [J/K^4/Kg]
        density_Nb = 6000; % [Kg/m^3]
        Cp300_Nb = 400; % [J/K/Kg]
        Cp_low_NC = (beta_Nb * (T^3)) + (gamma_Nb * T); % [J/K/Kg] NORMAL
        Cp_NbTi = 1 / ((1 / Cp300_Nb) + (1 / Cp_low_NC));
        Cp_NbTi = Cp_NbTi * density_Nb;

        gamma_Steel = 0.48; % [J/K^2/Kg]
        beta_Steel = 0.00075; % [J/K^4/Kg]
        density_Steel = 7900; % [Kg/m^3]
        Cp300_Steel = 500; % [J/K/Kg]  
        Cp_low_NC = (beta_Steel * (T^3)) + (gamma_Steel * T); % [J/K/Kg] NORMAL
        Cp_Steel = 1 / ((1 / Cp300_Steel) + (1 / Cp_low_NC));
        Cp_Steel = Cp_Steel * density_Steel;

        % Equivalent Specif Heat (volumetric)
        Cp_v_Nb3Sn = (Cp_Nb3Sn * SC_cross_sect + Cp_Cu * Cu_cross_section) / (Cu_cross_section + SC_cross_sect);
        Cp_v_NbTi  = (Cp_NbTi * SC_cross_sect + Cp_Cu * Cu_cross_section) / (Cu_cross_section + SC_cross_sect);
        Cp_v_Steel  = (Cp_Steel * SC_cross_sect + Cp_Cu * Cu_cross_section) / (Cu_cross_section + SC_cross_sect);      
    end
end
