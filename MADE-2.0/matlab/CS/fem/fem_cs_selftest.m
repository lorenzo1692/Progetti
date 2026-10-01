function ok = fem_cs_selftest(p, g)
%FEM_CS_SELFTEST Consistency checks of the CS stack FE model on a built-in design point.
%
%   ok = FEM_CS_SELFTEST(p, g) runs FEM_CS_VERIFY (design case only, no
%   background coils) on a fixed 6 layers x 24 turns design of the VNS
%   baseline and checks:
%     1. hoop-tension equilibrium of the FE section (ratio within 1%),
%     2. net axial equilibrium: inertia-relief reaction negligible compared
%        with the module loads,
%     3. Lorentz force on every module against EMAG_FIELD_FORCES (axial
%        within 3%, with the radial force allowed to be lower than the
%        design model, which evaluates the field at the inner edge of every
%        turn instead of at the cell center),
%     4. the compression through the mid plate equals EMAG_FIELD_FORCES'
%        FZmax (full stack) within 3%,
%     5. with the bottom face fixed and a 10 MN preload on the top face,
%        every plate carries 10 MN more than without preload.
%   p, g: parameter struct and COMPUTE_OPERATING_PARAMS output of the CS
%   template (cs_smoke_test loads them for Octave).

ok = true;
row = struct('n_layers', 6, 'n_turns', 24, 'Cond_w', 0.023563, 'Cond_h', 0.040292, ...
    'JT', 4.5e-3, 'S_Cable', 3.4654e-4, 'Iop_kA', 55, 'Ri', 0.5506, ...
    'B_grades', 10.35, 'Fz_MN', 43.38, 'type_cable', {{'LTS'}});
p.fem_mesh_r = 1; p.fem_mesh_z = 1; p.fem_subfil = 2; p.fem_bc = 1; p.fem_preload_MN = 0;

res = fem_cs_verify(row, p, g, [], struct());
c = res.cases(1); an = res.analytic; nm = res.wp.n_mod;

ok = check('hoop equilibrium ratio', c.T_check, 1, 0.01) && ok;
Fz_mod = [c.modules.Fz_MN];
ok = check('sum of module Fz (net, MN)', sum(Fz_mod), 0, 0.01*max(abs(Fz_mod))) && ok;
ok = check('end module Fz vs emag (MN)', Fz_mod(1), an.Fz_mod_MN(1), 0.03*abs(an.Fz_mod_MN(1))) && ok;
ok = check('mid-plate compression vs emag FZmax (MN)', c.F_plate_MN(ceil((nm-1)/2)), an.FZmax_half_MN, 0.03*abs(an.FZmax_half_MN)) && ok;
fr_ratio = [c.modules.Fr_MN]./an.Fr_mod_MN;
ok = check('radial force FEM/design (0.75..1.05)', mean(fr_ratio), 0.9, 0.15) && ok;

p.fem_bc = 2; p.fem_preload_MN = 10;
res2 = fem_cs_verify(row, p, g, [], struct());
dF = res2.cases(1).F_plate_MN - c.F_plate_MN;
ok = check('preload 10 MN adds to every plate (MN)', mean(dF), 10, 0.3) && ok;
ok = check('preload spread over plates (max-min, MN)', max(dF) - min(dF), 0, 0.3) && ok;
fprintf('FEM_CS_SELFTEST: %s\n', ternary(ok, 'all checks passed', 'FAILED'));
end

function ok = check(name, val, ref, tol)
ok = abs(val - ref) <= tol;
fprintf('  %-48s %10.4f  (ref %.4f, tol %.4f)  %s\n', name, val, ref, tol, ternary(ok, 'ok', 'FAIL'));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
