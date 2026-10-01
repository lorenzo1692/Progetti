function out = size_cicc_cable(in)
%SIZE_CICC_CABLE Size one PFC grade's jacket geometry against the ring Tresca stress (and optionally fatigue) limit.
%
%   out = SIZE_CICC_CABLE(in) grows the jacket thickness JT (rectangular
%   CICC, the only shape_cable PF_opt_VNS.m exercises) until BOTH the hoop
%   stress and the Tresca stress from EQV_STRESS_COIL_RING_CICC drop
%   within [in.S_hoop_allow, in.S_hoop_allow*in.Tresca_factor], and derives
%   the resulting cable/jacket geometry (SC_w, SC_h, Cond_w, S_CICC, S_JT,
%   Ri_grade, Re_grade).
%
%   Fatigue sizing option (in.fcgr_mode, input parameter fcgr_mode):
%     0 / 1 : jacket sized on stress only (1 only affects reporting in
%             SCAN_WP_DESIGNS, which then computes plasma_cycles);
%     2     : once the stress limits are met, JT keeps growing until the
%             FCGR life reaches in.plasma_cycles_min [k cycles].
%
%   Unlike MADE-2.0/matlab/CS/physics/size_cicc_cable.m, this sizing loop
%   recomputes the field (EMAG_FIELD_FORCES) and the full ring stress
%   (hoop AND Tresca, not just hoop at Fz=0) on every JT increment: the
%   ring formula needs the field at both the inner and outer radius
%   (Bmin, Bmax), which move every time Cond_w changes, and the axial
%   force FZ is an externally supplied scenario value applied from the
%   very first iteration (see manuale PFC, "fidelity note: axial force").
%   This mirrors PF_opt_VNS.m's own control flow exactly; the CS-style
%   split ("size hoop-only, check Tresca separately after") does not apply
%   here.
%
%   in fields (all scalars unless noted):
%     Cond_h, S_Cable_var, r_SC, tins, type_cable (1x1 cell string),
%     Cond_w1 (Cond_w of grade 1, used for Ri_grade when var>1; for var=1
%     the width sized in the current iteration is used, as in PF_opt_VNS.m),
%     n_layers_var, R_center, dy, n_t, n_l, n_g, Iop, FZ, WP_h, var,
%     JT_min, JT_step, S_hoop_allow, Tresca_factor, max_iter
%   optional fields: ring_opts (struct for EQV_STRESS_COIL_RING_CICC),
%     fcgr_mode, plasma_cycles_min, fp (struct for FCGR)
%
%   out fields:
%     feasible, JT, SC_h, SC_w, R_J, Cond_w, S_CICC, S_JT, Ri_grade,
%     Re_grade, Bsum, Bmin, Bmax, FRmax, FZmax, BR, BZ, FR, FZ, S_hoop,
%     S_rad, S_ver, S_T, plasma_cycles (NaN unless fatigue was evaluated),
%     iterations
%
%   feasible=false when the jacket grows until the cable no longer fits
%   (SC_h or SC_w <= 0): the legacy loop had no such guard (it relied on
%   the geometric checks that follow to reject the candidate).
%
%   Relocated and renamed from the PF legacy archive (PF_opt_VNS.m, the
%   shape_cable==201 branch inside the main scan loop), consolidating the
%   per-candidate JT-growing while-loop into a single reusable function.

ring_opts = [];
if isfield(in, 'ring_opts'), ring_opts = in.ring_opts; end
fcgr_mode = 0;
if isfield(in, 'fcgr_mode'), fcgr_mode = in.fcgr_mode; end
fp = [];
if isfield(in, 'fp'), fp = in.fp; end

JT = in.JT_min;
S_hoop = Inf;
S_T = Inf;
plasma_cycles = NaN;
iterations = 0;
feasible = true;

SC_h = NaN; SC_w = NaN; R_J = NaN; Cond_w = NaN; S_CICC = NaN; S_JT = NaN;
Ri_grade = NaN; Re_grade = NaN;
Bsum = NaN; Bmin = NaN; Bmax = NaN; FRmax = NaN; FZmax = NaN;
BR = NaN; BZ = NaN; FR = NaN; FZ = NaN; S_rad = NaN; S_ver = NaN;

keep_growing = true;
while keep_growing
    iterations = iterations + 1;
    if iterations > in.max_iter
        error('size_cicc_cable:not_converged', ...
            'JT sizing did not converge after %d iterations: check the input parameters.', iterations);
    end

    JT = JT + in.JT_step;
    SC_h = in.Cond_h - 2*JT - 2*in.tins;
    if SC_h <= 0
        feasible = false;
        break
    end
    SC_w = (in.S_Cable_var + (4-pi)*in.r_SC^2)/SC_h;
    R_J = in.r_SC + JT;
    Cond_w = SC_w + 2*JT + 2*in.tins;
    S_CICC = (in.Cond_h - 2*in.tins)*(Cond_w - 2*in.tins) - (4-pi)*R_J^2;
    S_JT = S_CICC - in.S_Cable_var;

    if in.var == 1
        Cond_w_first = Cond_w; % grade 1 sizes its own width in this very iteration
    else
        Cond_w_first = in.Cond_w1;
    end
    Ri_grade = in.R_center - Cond_w_first*in.n_layers_var/2;
    Re_grade = Ri_grade + Cond_w*in.n_layers_var;

    [Bsum, Bmin, Bmax, FRmax, FZmax, BR, BZ, FR, FZ, ~, ~] = emag_field_forces( ...
        in.dy, Cond_w, Ri_grade, in.n_t, in.n_l, in.n_g, in.Iop);

    [S_hoop, S_rad, S_ver, S_T] = eqv_stress_coil_ring_cicc(in.FZ, Re_grade, Ri_grade, ...
        in.Cond_h, Cond_w, JT, SC_h, SC_w, in.tins, in.type_cable, S_CICC, S_JT, in.var, ...
        in.Iop, Bmin, Bmax, in.WP_h, in.n_layers_var, ring_opts);

    stress_ok = max(abs(S_hoop)) <= in.S_hoop_allow && max(abs(S_T)) <= in.S_hoop_allow*in.Tresca_factor;
    if ~stress_ok
        continue
    end

    if fcgr_mode == 2
        plasma_cycles = fcgr(JT, Cond_w, max(abs(S_hoop))*1e-6, fp);
        keep_growing = plasma_cycles < in.plasma_cycles_min;
    else
        keep_growing = false;
    end
end

out.feasible = feasible;
out.JT = JT;
out.SC_h = SC_h;
out.SC_w = SC_w;
out.R_J = R_J;
out.Cond_w = Cond_w;
out.S_CICC = S_CICC;
out.S_JT = S_JT;
out.Ri_grade = Ri_grade;
out.Re_grade = Re_grade;
out.Bsum = Bsum; out.Bmin = Bmin; out.Bmax = Bmax;
out.FRmax = FRmax; out.FZmax = FZmax;
out.BR = BR; out.BZ = BZ; out.FR = FR; out.FZ = FZ;
out.S_hoop = S_hoop; out.S_rad = S_rad; out.S_ver = S_ver; out.S_T = S_T;
out.plasma_cycles = plasma_cycles;
out.iterations = iterations;
end
