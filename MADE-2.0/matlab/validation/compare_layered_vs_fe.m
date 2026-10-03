function C = compare_layered_vs_fe(row, p, name, opts_lc, fe)
%COMPARE_LAYERED_VS_FE Side-by-side review of the layered model against the 2D FE.
%
%   C = COMPARE_LAYERED_VS_FE(row, p, name) runs the 2D FE (WP_MECH_SURROGATE,
%   primary load case only: Lorentz + vertical force, no cool-down) and the
%   layered model (WP_LAYERED_CYLINDER) on the same design and compares,
%   quantity by quantity, where the loads go:
%     1. hoop (wedging) force across the centre plane x = 0, split into
%        nose / WP region (lateral walls + WP) / plasma-side plate;
%     2. axial force, same split (FE: also lateral walls and WP apart);
%     3. membrane stress components on the case SCLs (nose centreline,
%        plasma-side plate, lateral wall at 50%);
%     4. per WP layer: hoop and radial stress averaged over the layer at the
%        centre plane (FE) against the ring hoop stress and the column
%        pressure (model), and the jacket Pm.
%   C = COMPARE_LAYERED_VS_FE(row, p, name, opts_lc, fe) reuses a FE result.
%   Units in C: forces per unit length [N/m] (hoop) and [N] (axial), stresses [Pa].
%   docs/MODELLO_A_STRATI.md, section 10.

if nargin < 4 || isempty(opts_lc), opts_lc = struct(); end
if nargin < 5 || isempty(fe)
    T_ref = 293; if isfield(p, 'T_ref'), T_ref = p.T_ref; end
    fe = wp_mech_surrogate(row, p, struct('classify', 0, 'T_op', T_ref));   % no cool-down: total = primary
end
lc = wp_layered_cylinder(row, p, opts_lc);

% ---- radii -------------------------------------------------------------------
R = lc.rings; Rb = R(1).r1; Rj = R(1).r2; Ri = row.Ri_; Rp = R(end).r1;   % plate from Rp to Ri
E = fe.elem; x = E.xy(:,1); y = E.xy(:,2); A = E.area; S = E.S; m = E.mat;
is_case = m == 4;
reg = zeros(size(m));                                  % 1 nose, 2 walls, 3 plate, 4 WP, 5 wedge insulation
reg(is_case & y < Rj) = 1; reg(is_case & y >= Rj & y <= Rp) = 2; reg(is_case & y > Rp) = 3;
reg(~is_case & m ~= 5) = 4; reg(m == 5) = 5;

% ---- 1. hoop force across x = 0 -------------------------------------------------
band = 4e-3; nb = max(ceil((Ri - Rb)/1e-3), 50);
edges = linspace(Rb, Ri, nb + 1); yc = (edges(1:end-1) + edges(2:end))/2;
sel = abs(x) < band;
sx_bin = nan(1, nb); sy_bin = nan(1, nb);
for i = 1:nb
    q = sel & y >= edges(i) & y < edges(i+1);
    if any(q), sx_bin(i) = sum(S(q,1).*A(q))/sum(A(q)); sy_bin(i) = sum(S(q,2).*A(q))/sum(A(q)); end
end
ok = ~isnan(sx_bin);
sx_bin = interp1(yc(ok), sx_bin(ok), yc, 'linear', 'extrap');
sy_bin = interp1(yc(ok), sy_bin(ok), yc, 'linear', 'extrap');
dy = edges(2) - edges(1);
H_fe = [sum(sx_bin(yc < Rj)), sum(sx_bin(yc >= Rj & yc <= Rp)), sum(sx_bin(yc > Rp))]*dy;
H_lc = zeros(1, 3);
for i = 1:numel(R)
    s = lc.sample(i); h = trapz(s.r, s.st);
    if strcmp(R(i).type, 'nose'), H_lc(1) = H_lc(1) + h;
    elseif strcmp(R(i).type, 'plate'), H_lc(3) = H_lc(3) + h;
    else, H_lc(2) = H_lc(2) + h; end
end

% ---- 2. axial force ------------------------------------------------------------------
Fz_fe = zeros(1, 5);
for k = 1:5, Fz_fe(k) = sum(S(reg == k, 3).*A(reg == k)); end
Fz_lc = zeros(1, 3);
for i = 1:numel(R)
    s = lc.sample(i); f = trapz(s.r, s.sz.*s.r)*2*pi/p.n_TF;
    if strcmp(R(i).type, 'nose'), Fz_lc(1) = Fz_lc(1) + f;
    elseif strcmp(R(i).type, 'plate'), Fz_lc(3) = Fz_lc(3) + f;
    else, Fz_lc(2) = Fz_lc(2) + f; end
end

% ---- 3. case SCL membrane components (FE: [sx sy sz txy]; model: [sr st sz]) ---
scl = fe.case_scl; nm = {scl.name};
g = @(t) scl(strcmp(nm, t)).S_membrane;
nl = row.n_layers; mid = round(nl/2);
mem = struct('name', {'nose centreline', 'plasma-side plate', 'side wall 50%'}, ...
    'fe_hoop_radial_axial', {g('nose centreline')([1 2 3]), g('plasma-side plate')([1 2 3]), g('side wall 50%')([1 2 3])}, ...
    'lc_hoop_radial_axial', {lc.nose.membrane([2 1 3]), lc.plate.membrane([2 1 3]), lc.layer(mid).wall([2 1 3])});

% ---- 4. per layer --------------------------------------------------------------------
lay = struct('k', {}, 'st_fe', {}, 'st_lc', {}, 'sr_fe', {}, 'sr_lc', {}, 'Pm_fe', {}, 'Pm_lc', {});
pc = lc.p_column;
for i = 1:numel(R)
    if ~strcmp(R(i).type, 'layer'), continue, end
    k = R(i).layer; q = yc >= R(i).r1 & yc <= R(i).r2;
    sr_lc = lc.layer(k).sr; if ~isempty(pc), sr_lc = -pc(k); end
    lay(k) = struct('k', k, 'st_fe', mean(sx_bin(q)), 'st_lc', lc.layer(k).st, 'sr_fe', mean(sy_bin(q)), ...
        'sr_lc', sr_lc, 'Pm_fe', fe.layer.Pm(k), 'Pm_lc', lc.Pm_jacket(k));
end

C = struct('name', name, 'Rb', Rb, 'Rj', Rj, 'Rp', Rp, 'Ri', Ri, 'nose_mm', 1e3*(Rj - Rb), ...
    'H_fe', H_fe, 'H_lc', H_lc, 'Fz_fe', Fz_fe, 'Fz_lc', Fz_lc, 'Fz_total_fe', fe.Fz_integral, ...
    'mem', mem, 'layer', lay, 'profile', struct('y', yc, 'sx', sx_bin, 'sy', sy_bin), ...
    'lc_sample', lc.sample, 'lc_types', {{R.type}}, 'valid', fe.valid);

% ---- report --------------------------------------------------------------------------
fprintf('\n=== %s: nose %.0f mm (FE geometry), FE valid %d ===\n', name, C.nose_mm, fe.valid);
fprintf('1. Hoop force across x = 0 [MN/m]   nose / walls+WP / plate / total\n');
fprintf('   FE     %7.2f %7.2f %7.2f | %7.2f   shares %.0f / %.0f / %.0f %%\n', H_fe/1e6, sum(H_fe)/1e6, 100*H_fe/sum(H_fe));
fprintf('   model  %7.2f %7.2f %7.2f | %7.2f   shares %.0f / %.0f / %.0f %%\n', H_lc/1e6, sum(H_lc)/1e6, 100*H_lc/sum(H_lc));
fprintf('2. Axial force [MN]   nose / walls+WP / plate / total (FE: walls %.2f, WP %.2f, wedge ins %.2f)\n', ...
    Fz_fe(2)/1e6, Fz_fe(4)/1e6, Fz_fe(5)/1e6);
Fz3 = [Fz_fe(1), Fz_fe(2) + Fz_fe(4) + Fz_fe(5), Fz_fe(3)];
fprintf('   FE     %7.2f %7.2f %7.2f | %7.2f   shares %.0f / %.0f / %.0f %%\n', Fz3/1e6, sum(Fz3)/1e6, 100*Fz3/sum(Fz3));
fprintf('   model  %7.2f %7.2f %7.2f | %7.2f   shares %.0f / %.0f / %.0f %%\n', Fz_lc/1e6, sum(Fz_lc)/1e6, 100*Fz_lc/sum(Fz_lc));
fprintf('3. Membrane stresses [MPa]          hoop    radial   axial  (FE / model)\n');
for i = 1:numel(mem)
    a = mem(i).fe_hoop_radial_axial/1e6; b = mem(i).lc_hoop_radial_axial/1e6;
    fprintf('   %-18s %6.0f/%-6.0f %6.0f/%-6.0f %6.0f/%-6.0f\n', mem(i).name, a(1), b(1), a(2), b(2), a(3), b(3));
end
fprintf('4. Layers: hoop (centre plane, layer mean) | radial | jacket Pm  [MPa]  FE / model\n');
for k = 1:numel(lay)
    L = lay(k);
    fprintf('   L%-2d %6.0f/%-6.0f %6.0f/%-6.0f %6.0f/%-6.0f\n', k, L.st_fe/1e6, L.st_lc/1e6, L.sr_fe/1e6, L.sr_lc/1e6, ...
        L.Pm_fe/1e6, L.Pm_lc/1e6);
end
end
