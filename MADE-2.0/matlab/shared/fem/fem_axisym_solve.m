function sol = fem_axisym_solve(nodes, elems, D, fgp, bc)
%FEM_AXISYM_SOLVE Axisymmetric (r,z) linear-elastic finite element solver, bilinear quads.
%
%   sol = FEM_AXISYM_SOLVE(nodes, elems, D, fgp, bc) solves the stress
%   problem of a ring (body of revolution, no torsion) discretized on its
%   (r,z) cross-section with 4-node isoparametric elements and 2x2 Gauss
%   integration. The hoop direction is the third axis, so the strain
%   vector is [eps_r; eps_z; eps_theta; gamma_rz] with eps_theta = u/r.
%
%   Inputs
%     nodes  nn x 2   node coordinates [r z] in m (r > 0)
%     elems  ne x 4   node numbers, counter-clockwise (r-z- , r+z- , r+z+ , r-z+)
%     D      4x4 or 4x4xne   constitutive matrix [Pa], stress order
%            [sigma_r; sigma_z; sigma_theta; tau_rz] (see FEM_WP_MATERIAL)
%     fgp    ne x 4 x 2   body force density at the 4 Gauss points of every
%            element, [f_r f_z] in N/m^3 (e.g. the Lorentz force J x B);
%            may be empty (no load, e.g. for a stiffness-only call)
%     bc     struct with
%              .type  'inertia_relief' : the net axial force is removed by a
%                     uniform opposite body force (equivalent to carrying
%                     it as a stress spread over the whole volume), and one
%                     node is fixed axially to remove the rigid-body mode;
%                     'dofs' : constrain the dofs in bc.dofs (zero value).
%              .dofs  vector of constrained dofs, dof = 2*(node-1)+{1 (r), 2 (z)}
%              .fext  (optional) 2nn x 1 extra nodal forces [N], added after the
%                     inertia relief (CS copy only: used for the stack preload)
%
%   Output sol:
%     u        2nn x 1  displacements [u1 w1 u2 w2 ...] [m]
%     sgp      ne x 4 x 4  stress at the Gauss points [sr sz sth trz] [Pa]
%     sel      ne x 4      element-average stress [sr sz sth trz] [Pa]
%     vol      ne x 1      element volume (full revolution) [m^3]
%     rcen     ne x 1      element centroid radius [m]
%     Fnet     1 x 2       applied net force [Fr Fz] over the full revolution
%                          (radial component is the integrated radial load;
%                          for a ring only the axial one is a resultant) [N]
%     Fz_relief            uniform axial body force added by inertia relief [N/m^3]
%     react    sum of the axial reactions at the constrained dofs [N]
%
%   All forces are for the complete 360 degree ring (2*pi*r in every
%   integral), so axial resultants are directly in N.

nn = size(nodes, 1);
ne = size(elems, 1);
ndof = 2*nn;

if size(D, 3) == 1
    D = repmat(D, [1 1 ne]);
end
if isempty(fgp)
    fgp = zeros(ne, 4, 2);
end

gp = [-1 -1; 1 -1; 1 1; -1 1]/sqrt(3);
Nw = 1; % Gauss weights (2x2 rule: all 1)

ii = zeros(64*ne, 1); jj = zeros(64*ne, 1); vv = zeros(64*ne, 1);
f = zeros(ndof, 1);
vol = zeros(ne, 1); rcen = zeros(ne, 1);
Bstore = zeros(4, 8, 4, ne);
cnt = 0;

for e = 1:ne
    en = elems(e, :);
    xe = nodes(en, :); % 4 x 2 [r z]
    edofs = reshape([2*en-1; 2*en], 1, 8);
    ke = zeros(8, 8);
    fe = zeros(8, 1);
    for g = 1:4
        xi = gp(g, 1); eta = gp(g, 2);
        N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
        dNxi = 0.25*[-(1-eta), (1-eta), (1+eta), -(1+eta)];
        dNeta = 0.25*[-(1-xi), -(1+xi), (1+xi), (1-xi)];
        J = [dNxi; dNeta]*xe;   % [dr/dxi dz/dxi; dr/deta dz/deta]
        detJ = det(J);
        dN = J\[dNxi; dNeta];   % row 1: dN/dr, row 2: dN/dz
        r = N*xe(:, 1);
        B = zeros(4, 8);
        for i = 1:4
            B(1, 2*i-1) = dN(1, i);
            B(2, 2*i)   = dN(2, i);
            B(3, 2*i-1) = N(i)/r;
            B(4, 2*i-1) = dN(2, i);
            B(4, 2*i)   = dN(1, i);
        end
        Bstore(:, :, g, e) = B;
        wgt = 2*pi*r*detJ*Nw;
        ke = ke + B'*D(:, :, e)*B*wgt;
        fg = squeeze(fgp(e, g, :));
        for i = 1:4
            fe(2*i-1) = fe(2*i-1) + N(i)*fg(1)*wgt;
            fe(2*i)   = fe(2*i)   + N(i)*fg(2)*wgt;
        end
        vol(e) = vol(e) + wgt;
        rcen(e) = rcen(e) + r*wgt;
    end
    rcen(e) = rcen(e)/vol(e);
    [I, Jc] = meshgrid(edofs, edofs);
    ii(cnt+1:cnt+64) = I(:); jj(cnt+1:cnt+64) = Jc(:); vv(cnt+1:cnt+64) = ke(:);
    cnt = cnt + 64;
    f(edofs) = f(edofs) + fe;
end
K = sparse(ii, jj, vv, ndof, ndof);

%% Net applied force and boundary conditions
Fnet = [sum(f(1:2:end)), sum(f(2:2:end))];
Fz_relief = 0;
bctype = 'dofs';
if isfield(bc, 'type'), bctype = bc.type; end

switch bctype
    case 'inertia_relief'
        Fz_relief = -Fnet(2)/sum(vol);       % uniform axial body force [N/m^3]
        fr = zeros(ndof, 1);
        for e = 1:ne
            en = elems(e, :);
            xe = nodes(en, :);
            edofs = reshape([2*en-1; 2*en], 1, 8);
            fe = zeros(8, 1);
            for g = 1:4
                xi = gp(g, 1); eta = gp(g, 2);
                N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
                dNxi = 0.25*[-(1-eta), (1-eta), (1+eta), -(1+eta)];
                dNeta = 0.25*[-(1-xi), -(1+xi), (1+xi), (1-xi)];
                J = [dNxi; dNeta]*xe;
                wgt = 2*pi*(N*xe(:, 1))*det(J)*Nw;
                for i = 1:4
                    fe(2*i) = fe(2*i) + N(i)*Fz_relief*wgt;
                end
            end
            fr(edofs) = fr(edofs) + fe;
        end
        f = f + fr;
        fixed = 2; % axial dof of node 1: rigid-body translation only (reaction ~ 0)
    case 'dofs'
        fixed = bc.dofs(:)';
    otherwise
        error('fem_axisym_solve:unknown_bc', 'Unknown bc.type "%s".', bctype);
end

if isfield(bc, 'fext') && ~isempty(bc.fext)
    f = f + bc.fext(:);                   % extra nodal forces (e.g. stack preload)
end

free = setdiff(1:ndof, fixed);
u = zeros(ndof, 1);
u(free) = K(free, free)\f(free);

R = K*u - f;
react = sum(R(fixed(mod(fixed, 2) == 0)));

%% Stress recovery
sgp = zeros(ne, 4, 4);
sel = zeros(ne, 4);
for e = 1:ne
    en = elems(e, :);
    edofs = reshape([2*en-1; 2*en], 1, 8);
    ue = u(edofs);
    for g = 1:4
        s = D(:, :, e)*(Bstore(:, :, g, e)*ue);
        sgp(e, g, :) = s;
    end
    sel(e, :) = squeeze(mean(sgp(e, :, :), 2))';
end

sol.u = u;
sol.sgp = sgp;
sol.sel = sel;
sol.vol = vol;
sol.rcen = rcen;
sol.Fnet = Fnet;
sol.Fz_relief = Fz_relief;
sol.react = react;
sol.f = f;
end
