function THS = heat_balance_cicc_ode(N_Sc, N_Cu, d_fili, CunonCu, Iop0, B0, Tau_discharge, mat, d_cc, VF, costheta, S_tapes)
%HEAT_BALANCE_CICC_ODE Adiabatic hot-spot temperature of the CICC during a fast discharge.
%
%   THS = HEAT_BALANCE_CICC_ODE(...) integrates the adiabatic hot-spot
%   energy balance (Joule heating vs. specific heat, exponential current/
%   field decay after Tau_delay) with ode45 and returns the peak
%   temperature THS [K] reached during the discharge.
%
%   Relocated unchanged from the CS legacy archive (heat_balance_cicc_ode.m).
%   Used by SIZE_CONDUCTOR_CICC to bisect the number of copper strands
%   N_Cu against a hot-spot temperature limit.

Tau_delay = 0.5; % [s] detection + breaker opening delay
tspan = [0, (-Tau_discharge*log(1/Iop0) + Tau_delay)];

if mat == 2
    T0 = 15.0;
else
    T0 = 6.8;
end

options = odeset('RelTol', 1e-6);
[~, TF] = ode45(@temperatureODE, tspan, T0, options);

THS = round(max(TF));

    function dTdt = temperatureODE(t, T)
        if mat == 2
            SC_cross_sect = N_Sc * S_tapes;
            Cu_cross_sect_non_seg = 0;
        else
            SC_cross_sect = N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu));
            Cu_cross_sect_non_seg = N_Sc * pi * d_fili^2 / (4 * (1 + CunonCu));
        end

        Cu_cross_sect_seg = N_Cu * pi * d_fili^2 / 4;
        Cu_cross_section = Cu_cross_sect_seg + Cu_cross_sect_non_seg;

        if t <= Tau_delay
            Iop = Iop0;
            B = B0;
        else
            Iop = Iop0 * exp(-((t - Tau_delay) / Tau_discharge));
            B = B0 * exp(-((t - Tau_delay) / Tau_discharge));
        end

        [rho_Nb3Sn, rho_NbTi, rho_Steel, rhoCu_non_seg, rhoCu_seg] = calculateResistances(T, B);

        rho_Cu_parallel = (Cu_cross_section * rhoCu_non_seg * rhoCu_seg) / (rhoCu_seg * Cu_cross_sect_non_seg + rhoCu_non_seg * Cu_cross_sect_seg);
        Cu_res = rho_Cu_parallel / Cu_cross_section;
        Nb3Sn_res = rho_Nb3Sn / SC_cross_sect;
        NbTi_res = rho_NbTi / SC_cross_sect;
        Steel_res = rho_Steel / SC_cross_sect;
        cable_res_Nb3Sn = (1 / Cu_res + 1 / Nb3Sn_res)^(-1);
        cable_res_NbTi = (1 / Cu_res + 1 / NbTi_res)^(-1);
        cable_res_Steel = (1 / Cu_res + 1 / Steel_res)^(-1);
        rho_cable_Nb3Sn = cable_res_Nb3Sn * (Cu_cross_section + SC_cross_sect);
        rho_cable_NbTi = cable_res_NbTi * (Cu_cross_section + SC_cross_sect);
        rho_cable_Steel = cable_res_Steel * (Cu_cross_section + SC_cross_sect);

        [Cp_Cu, Cp_v_Nb3Sn, Cp_v_NbTi, Cp_v_Steel] = calculateSpecificHeat(T, SC_cross_sect, Cu_cross_section); %#ok<ASGLU>

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

        Q_gen = (Iop / (Cu_cross_section + SC_cross_sect))^2 * rho_cable;
        Q_abs = Cp_v;
        dTdt = Q_gen / Q_abs;
    end

    function [rho_Nb3Sn, rho_NbTi, rho_Steel, rhoCu_non_seg, rhoCu_seg] = calculateResistances(T, B)
        rho_Nb3Sn = (-1e-4 * T^2 + 0.0938 * T + 22.601) * 1e-8; % [ohm m]
        rho_NbTi = (0.0558 * T + 55.668) * 1e-8;                % [ohm m]
        rho_Steel = (76.2063 + 0.071375*(T-273) - 2.3109e-5*(T-273)^2) * 1e-8; % [ohm m]

        RRR = 300;
        rho1 = (1.171e-17 * T^4.49) / (1 + 4.5e-7 * T^3.35 * exp(-(50/T)^6.428));
        rho2 = (1.69e-8/RRR) + rho1 + 0.4531*((1.69e-8*rho1)/(RRR*rho1 + 1.69e-8));
        A = log10(1.553e-8 * B / rho2);
        a = -2.662 + 0.3168*A + 0.6229*A^2 - 0.1839*A^3 + 0.01827*A^4;
        rhoCu_seg = rho2 * (1 + 10^a); % [ohm m]

        RRR = 100;
        rho1 = (1.171e-17 * T^4.49) / (1 + 4.5e-7 * T^3.35 * exp(-(50/T)^6.428));
        rho2 = (1.69e-8/RRR) + rho1 + 0.4531*((1.69e-8*rho1)/(RRR*rho1 + 1.69e-8));
        A = log10(1.553e-8 * B / rho2);
        a = -2.662 + 0.3168*A + 0.6229*A^2 - 0.1839*A^3 + 0.01827*A^4;
        rhoCu_non_seg = rho2 * (1 + 10^a); % [ohm m]
    end

    function [Cp_Cu, Cp_v_Nb3Sn, Cp_v_NbTi, Cp_v_Steel] = calculateSpecificHeat(T, SC_cross_sect, Cu_cross_section)
        density = 8960; cp300 = 3.454e6; gamma = 0.011; beta = 0.0011;
        c_plow = beta*T^3 + gamma*T;
        Cp_Cu = 1 / ((1/cp300) + (1/(c_plow*density)));

        gamma_Nb = 0.1; beta_Nb = 0.001; density_Nb = 8040; Cp300_Nb = 210;
        Cp_low_NC = beta_Nb*T^3 + gamma_Nb*T;
        Cp_Nb3Sn = 1 / ((1/Cp300_Nb) + (1/Cp_low_NC)) * density_Nb;

        gamma_Nb = 0.145; beta_Nb = 0.0023; density_Nb = 6000; Cp300_Nb = 400;
        Cp_low_NC = beta_Nb*T^3 + gamma_Nb*T;
        Cp_NbTi = 1 / ((1/Cp300_Nb) + (1/Cp_low_NC)) * density_Nb;

        gamma_Steel = 0.48; beta_Steel = 0.00075; density_Steel = 7900; Cp300_Steel = 500;
        Cp_low_NC = beta_Steel*T^3 + gamma_Steel*T;
        Cp_Steel = 1 / ((1/Cp300_Steel) + (1/Cp_low_NC)) * density_Steel;

        Cp_v_Nb3Sn = (Cp_Nb3Sn*SC_cross_sect + Cp_Cu*Cu_cross_section) / (Cu_cross_section + SC_cross_sect);
        Cp_v_NbTi  = (Cp_NbTi*SC_cross_sect  + Cp_Cu*Cu_cross_section) / (Cu_cross_section + SC_cross_sect);
        Cp_v_Steel = (Cp_Steel*SC_cross_sect + Cp_Cu*Cu_cross_section) / (Cu_cross_section + SC_cross_sect);
    end
end
