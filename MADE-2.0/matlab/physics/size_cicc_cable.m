function out = size_cicc_cable(in)
%SIZE_CICC_CABLE Size the CICC/jacket cross-section for one WP grade.
%
%   out = SIZE_CICC_CABLE(in) replaces the three near-identical
%   turns/jacket sizing while-loops of the original monolithic script with
%   a single reusable function. For a rectangular cable (shape_cable=201)
%   it grows the jacket thickness JT until the radial membrane stress
%   drops within the allowable, exactly as before; a round/RIS cable
%   (shape_cable=200) uses the original closed-form sizing.
%
%   Also carries the original fix #4: the iterative sizing is capped at
%   in.max_iter passes and raises a clear error instead of looping forever
%   if the input parameters never let it converge.
%
%   in fields (all scalars):
%     Cond_w, S_Cable_var, r_SC, tins, E_jckt, E_cbl, E_ins, shape_cable,
%     p_rs, S_z_JT, S_amm_JT, safety_membrane, min_JT, max_iter
%
%   out fields:
%     SC_w, SC_h, R_J, Cond_h, JT, Ke_rad, Ke_tor, S_CICC, S_JT, iterations

switch in.shape_cable
    case 200
        SC_w = 2*sqrt(in.S_Cable_var/pi);
        SC_h = SC_w;                       % SC cable height
        R_J  = in.r_SC;                    % Jacket corner curvature radius
        Cond_h = in.Cond_w;                % RIS_ cable height
        JT = (in.Cond_w - 2*in.tins - SC_w)/2;
        [Ke_rad, Ke_tor] = cavo_stiffness(in.E_jckt, in.E_cbl, in.E_ins, ...
            JT, in.tins, Cond_h, in.Cond_w, SC_w, SC_h);
        iterations = 0;

    case 201
        JT = in.min_JT;
        S_rm = Inf;
        iterations = 0;
        while S_rm > in.S_amm_JT/in.safety_membrane
            iterations = iterations + 1;
            if iterations > in.max_iter
                error('size_cicc_cable:not_converged', ...
                    'JT sizing did not converge after %d iterations: check the input parameters.', iterations);
            end
            JT = JT + 1e-4;
            SC_w = in.Cond_w - JT*2 - in.tins*2;                          % SC cable width
            SC_h = (in.S_Cable_var + (4-pi)*in.r_SC^2)/SC_w;              % SC cable height
            R_J  = in.r_SC + JT;                                          % Jacket corner curvature radius
            Cond_h = SC_h + JT*2 + in.tins*2;                             % Rectangular cable height
            [Ke_rad, Ke_tor] = cavo_stiffness(in.E_jckt, in.E_cbl, in.E_ins, ...
                JT, in.tins, Cond_h, in.Cond_w, SC_w, SC_h);
            K_jckt = 2*JT/Cond_h*in.E_jckt;
            dcr_jckt = K_jckt/Ke_rad;
            r_steel = (in.Cond_w - 2*in.tins)/(2*JT);
            S_rm = in.p_rs*r_steel*dcr_jckt + in.S_z_JT;                  % Radial membrane stress, innermost jacket layer
        end

    otherwise
        error('size_cicc_cable:unknown_shape', 'Unknown shape_cable code: %g', in.shape_cable);
end

S_CICC = ((in.Cond_w - 2*in.tins)*(Cond_h - 2*in.tins)) - ((4-pi)*R_J^2); % Non-insulated cable cross-section
S_JT   = S_CICC - in.S_Cable_var;                                        % Steel cross-section

out.SC_w = SC_w;
out.SC_h = SC_h;
out.R_J = R_J;
out.Cond_h = Cond_h;
out.JT = JT;
out.Ke_rad = Ke_rad;
out.Ke_tor = Ke_tor;
out.S_CICC = S_CICC;
out.S_JT = S_JT;
out.iterations = iterations;
end

function [Ke_rad, Ke_tor] = cavo_stiffness(E_jckt, E_cbl, E_ins, JT, tins, Cond_h, Cond_w, SC_w, SC_h)
Ke_rad = 2*E_jckt*JT/Cond_h + 2*tins*E_ins/Cond_h + ...
    (1/(E_cbl*SC_w/SC_h) + 2/(E_jckt*Cond_w/JT) + 2/(E_ins*Cond_w/tins))^-1;
Ke_tor = 2*E_jckt*JT/Cond_w + 2*tins*E_ins/Cond_w + ...
    (1/(E_cbl*SC_h/SC_w) + 2/(E_jckt*JT/Cond_w) + 2/(E_ins*tins/Cond_w))^-1;
end
