function ok = fem_selftest(verbose)
%FEM_SELFTEST Verification of the axisymmetric FE solver against closed-form solutions.
%
%   ok = FEM_SELFTEST() runs the benchmarks below and returns true when all
%   pass (also printing the errors). Run it after any change to the FEM
%   modules; it needs no input file.
%
%   1. Thick ring with a uniform radial body force and free faces,
%      plane strain (axial displacement fixed at both ends): the closed-form
%      solution u(r) = -f r^2/(3(l+2m)) + C1 r/2 + C2/r gives sigma_r and
%      sigma_theta; the FE result is compared at every element centroid.
%   2. Force equilibrium: with inertia relief the reaction at the single
%      constrained dof and the net applied axial force must vanish; with a
%      rigid bottom support the reaction must equal the applied axial force.
%   3. Mesh convergence of test 1 (error decreases with refinement).
%   4. Uniform axial body force on a column (support at the bottom face):
%      the resultant of sigma_z across a section equals the load above it.
%   5. COMPUTE_SCENARIO_FORCES: the axial forces of an isolated system of
%      conductors sum to zero and parallel coaxial currents attract.

if nargin < 1, verbose = true; end
ok = true;

E = 205e9; nu = 0.3;
lam = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu));
Diso = [lam+2*mu, lam, lam, 0; lam, lam+2*mu, lam, 0; lam, lam, lam+2*mu, 0; 0 0 0 mu];

%% Test 1 and 3: Lame-type ring with radial body force
a = 0.5; b = 0.8; h = 0.2; f0 = 5e7; % [m], [m], [m], [N/m^3]
errs = zeros(1, 3);
meshes = [4 2; 8 2; 16 2];
for k = 1:3
    nr = meshes(k, 1); nz = meshes(k, 2);
    m = fem_mesh_rect(a, b, 0, h, nr, nz);
    ne = size(m.elems, 1);
    fgp = zeros(ne, 4, 2); fgp(:, :, 1) = f0;
    bc.type = 'dofs';
    bc.dofs = [2*m.bottom', 2*m.top']; % axial dof of the end faces
    sol = fem_axisym_solve(m.nodes, m.elems, Diso, fgp, bc);

    % closed form: (l+2m) g' = -f0, g = u'+u/r, sigma_r(a)=sigma_r(b)=0
    L2 = lam + 2*mu;
    % u = -f0 r^2/(3 L2) + C1 r/2 + C2/r
    sr = @(r, C1, C2) L2*(-2*f0*r/(3*L2) + C1/2 - C2./r.^2) + lam*(-f0*r/(3*L2) + C1/2 + C2./r.^2);
    A = [L2*(1/2) + lam*(1/2), L2*(-1/a^2) + lam*(1/a^2); L2*(1/2) + lam*(1/2), L2*(-1/b^2) + lam*(1/b^2)];
    rhs = -[sr(a, 0, 0); sr(b, 0, 0)];
    C = A\rhs;
    r = sol.rcen;
    u = -f0*r.^2/(3*L2) + C(1)*r/2 + C(2)./r;
    up = -2*f0*r/(3*L2) + C(1)/2 - C(2)./r.^2;
    s_r_an = L2*up + lam*u./r;
    s_t_an = lam*up + L2*u./r;
    err_r = max(abs(sol.sel(:, 1) - s_r_an))/max(abs(s_t_an));
    err_t = max(abs(sol.sel(:, 3) - s_t_an))/max(abs(s_t_an));
    errs(k) = max(err_r, err_t);
    if verbose
        fprintf('  [Lame ring] mesh %2dx%d : max error sigma_r %.3f%%, sigma_theta %.3f%% (of peak hoop %.1f MPa)\n', ...
            nr, nz, 100*err_r, 100*err_t, max(s_t_an)*1e-6);
    end
end
if errs(3) > 0.02
    ok = false; fprintf('  FAIL: Lame benchmark error %.2f%% > 2%% on the finest mesh\n', 100*errs(3));
end
if ~(errs(3) < errs(1))
    ok = false; fprintf('  FAIL: no mesh convergence (%.3f%% -> %.3f%%)\n', 100*errs(1), 100*errs(3));
end

%% Test 2: equilibrium with inertia relief and with rigid support
m = fem_mesh_rect(a, b, 0, h, 6, 6);
ne = size(m.elems, 1);
fgp = zeros(ne, 4, 2);
fgp(:, :, 1) = 3e7;
fgp(:, :, 2) = repmat(linspace(-4e7, 9e7, 4), ne, 1); % non-uniform axial load, net != 0
bc = struct('type', 'inertia_relief');
sol = fem_axisym_solve(m.nodes, m.elems, Diso, fgp, bc);
if abs(sol.react) > 1e-6*abs(sol.Fnet(2))
    ok = false; fprintf('  FAIL: inertia relief reaction %.3g N (net axial force %.3g N)\n', sol.react, sol.Fnet(2));
elseif verbose
    fprintf('  [equilibrium] inertia relief: reaction %.2e N for a net axial force of %.3e N\n', sol.react, sol.Fnet(2));
end
bc2.type = 'dofs'; bc2.dofs = 2*m.bottom';
sol2 = fem_axisym_solve(m.nodes, m.elems, Diso, fgp, bc2);
if abs(sol2.react + sol2.Fnet(2)) > 1e-8*abs(sol2.Fnet(2))
    ok = false; fprintf('  FAIL: rigid support reaction %.6g N vs applied %.6g N\n', sol2.react, sol2.Fnet(2));
elseif verbose
    fprintf('  [equilibrium] rigid support: reaction %.6e N, applied %.6e N\n', sol2.react, sol2.Fnet(2));
end

%% Test 4: column loaded by a uniform axial body force, supported at the bottom
nzc = 10;
m = fem_mesh_rect(a, b, 0, h, 3, nzc);
ne = size(m.elems, 1);
g0 = -2e7; % downward body force density [N/m^3]
fgp = zeros(ne, 4, 2); fgp(:, :, 2) = g0;
bc3.type = 'dofs'; bc3.dofs = 2*m.bottom';
sol3 = fem_axisym_solve(m.nodes, m.elems, Diso, fgp, bc3);
% The resultant of sigma_z across any section equals the load above it
% (sigma_z in tension positive): check the row of elements at mid-height.
jz = nzc/2;
row = find(m.iz == jz);
zc = (jz - 0.5)*h/nzc;
F_sec = sum(sol3.sel(row, 2).*sol3.vol(row))/(h/nzc);
F_above = g0*pi*(b^2 - a^2)*(h - zc);
rel = abs(F_sec - F_above)/abs(F_above);
if verbose
    fprintf('  [column] section force at z=%.3f m: FE %.4e N, load above %.4e N (diff %.2f%%)\n', zc, F_sec, F_above, 100*rel);
end
if rel > 0.05
    ok = false; fprintf('  FAIL: column section force differs by %.2f%%\n', 100*rel);
end

%% Test 5: scenario forces (Newton's third law and sign of the interaction)
gm.R = [1.0; 2.5; 1.5]; gm.Z = [1.0; -1.0; 0.2]; gm.dr = [0.2; 0.3; 0.2]; gm.dz = [0.3; 0.3; 0.3];
gm.MAt_signed = [2e6 -1e6; 2e6 3e6; -1.5e6 0.5e6];
Fzs = compute_scenario_forces(gm, 3, 3);
if max(abs(sum(Fzs, 1))) > 1e-9*max(abs(Fzs(:)))
    ok = false; fprintf('  FAIL: net axial force of the system is not zero (%.3g MN)\n', max(abs(sum(Fzs, 1))));
elseif verbose
    fprintf('  [forces] sum of the axial forces over all conductors: %.2e MN (largest single force %.2f MN)\n', max(abs(sum(Fzs, 1))), max(abs(Fzs(:))));
end
gm2.R = [1.0; 1.0]; gm2.Z = [1.0; -1.0]; gm2.dr = [0.1; 0.1]; gm2.dz = [0.1; 0.1]; gm2.MAt_signed = [1e6; 1e6];
Fz2 = compute_scenario_forces(gm2, 2, 2);
if ~(Fz2(1) < 0 && Fz2(2) > 0)
    ok = false; fprintf('  FAIL: two coaxial coils with parallel currents must attract each other\n');
end

if verbose
    if ok, fprintf('FEM selftest: PASS\n'); else, fprintf('FEM selftest: FAIL\n'); end
end
end
