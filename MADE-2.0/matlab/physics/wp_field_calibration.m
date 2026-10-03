function cal = wp_field_calibration(p, g, env)
%WP_FIELD_CALIBRATION Peak-field factor and layer profile vs WP width, once per machine.
%
%   cal = WP_FIELD_CALIBRATION(p, g, env) replaces the constant
%   p.corr_B_WP of the smeared field model with a calibration made at the
%   start of the scan, with the discrete 2D field model validated against
%   ANSYS (WP_PEAK_FIELD_FAST):
%
%     k(W, Iop)   = peak field on the conductor / Ampere field, as a
%                   function of the toroidal WP width W (the main effect:
%                   the narrower the WP, the more concentrated the current)
%                   and of the operating current Iop (cable self-field);
%     f(s; W, Iop) = per-layer peak / peak, as a function of the fraction
%                   s of the turns from that layer inward (s = 1 on the
%                   plasma-side layer). The smeared model assumes f = s,
%                   which under-predicts the inner layers by up to 2-3 T.
%
%   The scan then sizes layer k at  B_k = k(W, Iop) * B_Ampere * f(s_k).
%
%   Grid: W = every toroidal WP width the scan will use (one per case
%   wedge thickness lateral_w), Iop = p.field_cal_n_iop (default 5)
%   values between the feasible current limits. For each point, the
%   reference WPs are built as the scan builds them: every admissible
%   turns-per-layer count nt of env.turns_comb in the first layer, deeper
%   layers reduced by 2 turns at a time until they fit the sector with the
%   toroidal gap and the remaining turns added as further layers; cable
%   sized with CICC at the Ampere field x 1.1 and SIZE_CICC_CABLE for the
%   cell height, one grade; nt giving a cell outside the scan's geometric limits (cable
%   width and aspect ratio, cell size, every layer inside the sector with
%   the toroidal gap) are skipped. k and f are the MEAN over the admissible nt (the spread is
%   reported in cal.k_spread): the residual dependence on the layout is
%   left to the optional discrete verification of the scan.
%
%   B_Ampere is g.B_PHI_TF / p.corr_B_WP, the same expression the smeared
%   model starts from, so k replaces corr_B_WP exactly.
%
%   cal fields: W [1 x nW], Iop [1 x nI], k [nW x nI], k_spread [nW x nI],
%   s (profile abscissa, 1 x ns), f [nW x nI x ns], B_amp, n_ref (number
%   of reference WPs per point), elapsed.

t0 = tic;
Mu_0 = g.Mu_0;
theta = g.theta_TF;
B_amp = g.B_PHI_TF/p.corr_B_WP;
tins = p.turn_insulation_nominal*p.Increm;
Re1 = g.R_TF_Innerleg - p.dr_plasma_side - p.GoundIns;

lw = env.lateral_w_min:p.lateral_w_step:env.lateral_w_max;
W = 2*Re1*tan(theta/2) - 2*lw - 2*p.GoundIns;
n_iop = 5; if isfield(p, 'field_cal_n_iop') && ~isempty(p.field_cal_n_iop), n_iop = p.field_cal_n_iop; end
Imin = max(p.Iop_min, g.NI/max(env.n_spire_fsbl));
Imax = min(p.Iop_max, g.NI/min(env.n_spire_fsbl));
Iop = linspace(Imin, Imax, n_iop);
s = linspace(0, 1, 21);

nW = numel(W); nI = numel(Iop);
k = nan(nW, nI); ks = nan(nW, nI); f = nan(nW, nI, numel(s)); nref = zeros(nW, nI);
ref = zeros(0, 8);   % W Iop nt nl Cond_w Cond_h depth k, one row per reference WP

% one reference cable per current (grade 1 at ~ the real peak)
S_cab = zeros(1, nI);
for j = 1:nI
    N = ceil(g.NI/Iop(j));
    L = Mu_0*(p.n_TF*N*g.k_bf)^2*g.r_bf/2*(besseli(0, g.k_bf) + 2*besseli(1, g.k_bf) + besseli(2, g.k_bf))/p.n_TF;
    tau = max([g.Tau_discharge1, L*Iop(j)/p.V_MAX, 4]);
    [~, ~, ~, ~, S_cab(j)] = cicc(1.1*B_amp, Iop(j), tau, p.WP_SC_type, p.THS_max_LTS, p.THS_max_HTS);
end

p_rs = (1.1*B_amp)^2/(2*Mu_0);           % as the scan: pressure from the (calibrated) plasma-side field
for i = 1:nW
    for j = 1:nI
        N = ceil(g.NI/Iop(j));
        T_bf = axial_load_factor(p)*0.5*(g.k_bf*p.n_TF*(N*Iop(j))^2*Mu_0/(2*pi));
        S_z_JT = T_bf/(W(i)^2)/2;
        kk = []; ff = [];
        for nt = env.turns_comb
            cw = W(i)/nt;
            in = struct('Cond_w', cw, 'S_Cable_var', S_cab(j), 'r_SC_min', p.r_SC_min, 'r_SC_max', p.r_SC_max, ...
                'tins', tins, 'E_jckt', p.E_jckt, 'E_cbl', p.E_cbl_LTS, 'E_ins', p.E_ins, ...
                'shape_cable', p.shape_cable, 'p_rs', p_rs, 'S_z_JT', S_z_JT, 'S_amm_JT', p.S_amm_JT, ...
                'safety_membrane', p.safety_membrane, 'min_JT', p.min_JT, 'JT_step', p.JT_step, ...
                'max_iter', p.max_sizing_iterations);
            try
                sz = size_cicc_cable(in);
            catch
                continue
            end
            ch = sz.Cond_h;
            ar_max = 2; if isfield(p, 'max_cable_aspect_ratio') && ~isempty(p.max_cable_aspect_ratio), ar_max = p.max_cable_aspect_ratio; end
            if sz.SC_w <= p.min_SC_w || cw/ch < p.min_cable_aspect_ratio || cw/ch > ar_max || ...
                    cw < p.min_size_CICC || cw > p.max_size_CICC
                continue
            end
            % layers as the scan builds them: nt turns in the first layer;
            % a deeper layer that does not fit the sector with the toroidal
            % gap loses 2 turns at a time, and the turns left over go into
            % further layers at the back
            nt_l = zeros(1, 0); rem = N; Re_k = Re1;
            while rem > 0 && numel(nt_l) < p.maxdim
                Ri_k = Re_k - ch;
                n = min(nt, rem);
                while n > 0 && (2*Ri_k*tan(theta/2) - (n*cw + 2*p.GoundIns))/2 < p.toroidal_gap
                    n = n - 2;
                end
                if n <= 0 || Ri_k <= 0, break, end
                nt_l(end+1) = n; rem = rem - n; Re_k = Ri_k - p.INS_grades; %#ok<AGROW>
            end
            nl = numel(nt_l);
            if rem > 0 || nl < 2, continue, end
            row = struct('Iop', Iop(j), 'n_layers', nl, 'n_turns', nt_l, 'Cond_w', cw*ones(1, nl), ...
                'Cond_h', ch*ones(1, nl), 'JT', sz.JT*ones(1, nl), 'Ri_', g.R_TF_Innerleg, ...
                'shape_cable', p.shape_cable);
            B = wp_peak_field_fast(row, p);
            sk = fliplr(cumsum(fliplr(nt_l)))/N;                   % turns from layer k inward
            [su, iu] = unique(sk);
            fk = interp1([0 su], [B(end)/max(B), B(iu)/max(B)], s, 'linear', 'extrap');
            kk(end+1) = max(B)/B_amp; %#ok<AGROW>
            ref(end+1,:) = [W(i), Iop(j), nt, nl, cw, ch, nl*ch + (nl-1)*p.INS_grades, kk(end)]; %#ok<AGROW>
            ff(end+1,:) = fk; %#ok<AGROW>
        end
        if isempty(kk), continue, end
        k(i,j) = mean(kk); ks(i,j) = max(kk) - min(kk); nref(i,j) = numel(kk);
        f(i,j,:) = mean(ff, 1);
    end
end
% fill points without an admissible reference layout from the nearest
% valid one along Iop, then along W (the scan may still visit them)
for i = 1:nW
    ok = ~isnan(k(i,:));
    if any(ok) && ~all(ok)
        idx = find(ok);
        for j = find(~ok)
            [~, m] = min(abs(idx - j)); jj = idx(m);
            k(i,j) = k(i,jj); ks(i,j) = ks(i,jj); f(i,j,:) = f(i,jj,:);
        end
    end
end
for j = 1:nI
    ok = ~isnan(k(:,j));
    if any(ok) && ~all(ok)
        idx = find(ok);
        for i = find(~ok)'
            [~, m] = min(abs(idx - i)); ii = idx(m);
            k(i,j) = k(ii,j); ks(i,j) = ks(ii,j); f(i,j,:) = f(ii,j,:);
        end
    end
end
if any(isnan(k(:)))
    error('wp_field_calibration:empty', ...
        'No admissible reference WP for the field calibration: check Iop_min/max, min/max_size_CICC.');
end
cal = struct('W', W, 'lateral_w', lw, 'Iop', Iop, 'k', k, 'k_spread', ks, 's', s, 'f', f, ...
    'B_amp', B_amp, 'n_ref', nref, 'ref', ref, 'elapsed', toc(t0));
fprintf(['Field calibration: %d widths x %d currents, k = peak/Ampere %.3f-%.3f ' ...
    '(corr_B_WP = %.3f), layout spread <= %.3f, %.0f s\n'], nW, nI, min(k(:)), max(k(:)), ...
    p.corr_B_WP, max(ks(:)), cal.elapsed);
end
