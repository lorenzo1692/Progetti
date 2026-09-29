function out = size_cicc_cable(in)
%SIZE_CICC_CABLE Size one CICC grade's jacket geometry against the hoop stress limit.
%
%   out = SIZE_CICC_CABLE(in) grows the jacket thickness JT (rectangular
%   CICC, shape_cable=201) until the radial hoop stress from
%   EQV_STRESS_COIL_CICC (evaluated at Fz=0, i.e. the sizing-time check
%   used by the legacy driver) drops within in.S_hoop_allow/in.SF_hoop,
%   and derives the resulting cable/jacket geometry (SC_w, SC_h, Cond_w,
%   S_CICC, S_JT, Ri_grade). A round/RIS cable (shape_cable=200) uses the
%   original closed-form sizing instead.
%
%   in fields (all scalars unless noted):
%     Cond_h, S_Cable_var, r_SC, tins, shape_cable, type_cable (cell str),
%     PB, Re_this_grade, Re_WP_outer, n_layers, S_hoop_allow, SF_hoop,
%     JT_min, JT_step, max_iter
%
%   Re_this_grade is this grade's own current outer radius (Re_grades(var)
%   in the legacy driver); Re_WP_outer is the winding pack's fixed
%   outermost radius (Re_grades(1), constant across grades) - the legacy
%   stress check uses Re_WP_outer, not Re_this_grade, as the thick-cylinder
%   outer radius reference, while Ri_grade is derived from Re_this_grade.
%   Keep the two distinct: mixing them silently changes the stress result.
%
%   out fields:
%     SC_w, SC_h, JT, Cond_w, R_J, S_CICC, S_JT, Ri_grade, Re_grade_next,
%     iterations
%
%   Relocated and renamed from the CS legacy archive (CS_opt_VNS.m, the
%   per-grade shape_cable switch inside the main scan loop), consolidating
%   part of duplicazione #1/#6 of the manuale CS. In the ported
%   CS_opt_VNS.m path shape_cable is always 201 (every material branch set
%   coil.shape_cable(var)=201) - the shape_cable=200 branch is carried
%   over unchanged for fidelity/future use, but not currently exercised.

switch in.shape_cable
    case 200 % round / RIS
        SC_w = 2*sqrt(in.S_Cable_var/pi);
        SC_h = SC_w;
        R_J = in.r_SC;
        JT = (in.Cond_h - 2*in.tins - SC_w)/2;
        Cond_w = in.Cond_h;
        S_CICC = (in.Cond_h - 2*in.tins)*(Cond_w - 2*in.tins) - (4-pi)*R_J^2;
        S_JT = S_CICC - in.S_Cable_var;
        Ri_grade = in.Re_this_grade - Cond_w*in.n_layers;
        iterations = 0;

    case 201 % rectangular CICC
        JT = in.JT_min;
        S_hoop = Inf;
        iterations = 0;
        while max(S_hoop) > in.S_hoop_allow/in.SF_hoop
            iterations = iterations + 1;
            if iterations > in.max_iter
                error('size_cicc_cable:not_converged', ...
                    'JT sizing did not converge after %d iterations: check the input parameters.', iterations);
            end
            JT = JT + in.JT_step;
            SC_h = in.Cond_h - 2*JT - 2*in.tins;
            SC_w = (in.S_Cable_var + (4-pi)*in.r_SC^2)/SC_h;
            R_J = in.r_SC + JT;
            Cond_w = SC_w + 2*JT + 2*in.tins;
            S_CICC = (in.Cond_h - 2*in.tins)*(Cond_w - 2*in.tins) - (4-pi)*R_J^2;
            S_JT = S_CICC - in.S_Cable_var;
            Ri_grade = in.Re_this_grade - Cond_w*in.n_layers;

            [S_hoop, ~, ~, ~] = eqv_stress_coil_cicc(0, Ri_grade, in.Re_WP_outer, ...
                in.Cond_h, Cond_w, JT, SC_h, SC_w, in.tins, in.type_cable, S_CICC, S_JT, in.PB);
        end

    otherwise
        error('size_cicc_cable:unknown_shape', 'Unknown shape_cable code: %g', in.shape_cable);
end

out.SC_w = SC_w;
out.SC_h = SC_h;
out.JT = JT;
out.Cond_w = Cond_w;
out.R_J = R_J;
out.S_CICC = S_CICC;
out.S_JT = S_JT;
out.Ri_grade = Ri_grade;
out.Re_grade_next = Ri_grade - in.tins;
out.iterations = iterations;
end
