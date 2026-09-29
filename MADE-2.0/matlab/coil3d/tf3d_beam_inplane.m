function res = tf3d_beam_inplane(X, conn, sec, F, strips, supports, n_TF)
%TF3D_BEAM_INPLANE In-plane beam model of one TF coil with cyclic-symmetric intercoil strips.
%
%   res = TF3D_BEAM_INPLANE(X, conn, sec, F, strips, supports, n_TF) solves
%   the global structural model of the TF system for loads in the plane of
%   the coils and identical on all coils (TF Lorentz forces, the in-plane
%   PF forces), as the ANSYS model STR_360 (BEAM188 along the centreline of
%   each coil, SHELL181 strips between adjacent coils at the inner leg -
%   the wedged vault - and at the outer intercoil structure). With cyclic
%   symmetry one coil is enough:
%
%   - beam: 2D Timoshenko frame, 3 DOF per node (u_R radial, u_z vertical,
%     theta rotation about the toroidal axis); out-of-plane DOF are not
%     loaded and decouple for a doubly symmetric section;
%   - strip between coil nodes a, b and their images on the next coil:
%     plane-stress membrane of thickness t, length h = |X_b - X_a| and width
%     the chord c = 2 R sin(pi/n_TF). Cyclic symmetry gives the hoop strain
%     eps_h = u_R/R (both edges move radially by the same amount) and the
%     strain along the coil eps_s from the relative displacement of a and b;
%     no in-plane shear. Plus plate bending from the rotation gradient along
%     the strip. Each coil carries the energy of one strip per segment.
%
%   Inputs:
%     X        N x 2 node coordinates [R z] [m]
%     conn     E x 2 beam elements (node indices)
%     sec      struct: E, nu, A, I (bending about the toroidal axis), As
%              (shear area)
%     F        N x 2 nodal forces [F_R F_z] [N]
%     strips   S x 3 [node_a node_b thickness]
%     supports K x 2 [node dof] (dof 1 = u_R, 2 = u_z, 3 = theta)
%     n_TF     number of coils
%
%   res: U (N x 3), elem (E x 6 end forces in local axes at I and J:
%   [N_I V_I M_I N_J V_J M_J], N > 0 tension, V and M in the element
%   local frame, x from I to J, y = rotated +90 deg in the (R, z) plane),
%   strip (S x 4 [N_hoop N_s M_s eps_h] per unit width), reactions.

N = size(X, 1); ndof = 3*N;
Kt = zeros(0,1); Ki = Kt; Kj = Kt;
E = sec.E; G = E/(2*(1 + sec.nu));
kloc_all = cell(size(conn,1), 1); T_all = kloc_all;
for e = 1:size(conn, 1)
    a = conn(e,1); b = conn(e,2);
    d = X(b,:) - X(a,:); L = norm(d); c = d(1)/L; s = d(2)/L;
    phi = 12*E*sec.I/(G*sec.As*L^2);
    k1 = E*sec.A/L; k2 = 12*E*sec.I/(L^3*(1+phi)); k3 = 6*E*sec.I/(L^2*(1+phi));
    k4 = (4+phi)*E*sec.I/(L*(1+phi)); k5 = (2-phi)*E*sec.I/(L*(1+phi));
    kl = [ k1 0 0 -k1 0 0; 0 k2 k3 0 -k2 k3; 0 k3 k4 0 -k3 k5; ...
          -k1 0 0 k1 0 0; 0 -k2 -k3 0 k2 -k3; 0 k3 k5 0 -k3 k4];
    R3 = [c s 0; -s c 0; 0 0 1]; T = blkdiag(R3, R3);
    ke = T'*kl*T;
    dofs = [3*a-2:3*a, 3*b-2:3*b];
    [ii, jj] = ndgrid(dofs, dofs);
    Ki = [Ki; ii(:)]; Kj = [Kj; jj(:)]; Kt = [Kt; ke(:)]; %#ok<AGROW>
    kloc_all{e} = kl; T_all{e} = T;
end
% intercoil strips
Dm = E/(1 - sec.nu^2);
for q = 1:size(strips, 1)
    a = strips(q,1); b = strips(q,2); t = strips(q,3);
    d = X(b,:) - X(a,:); h = norm(d); es = d/h;
    Rm = (X(a,1) + X(b,1))/2; cw = 2*Rm*sin(pi/n_TF);
    % strains in terms of [uRa uza tha uRb uzb thb]
    Bh = [1/(2*Rm) 0 0 1/(2*Rm) 0 0];
    Bs = [-es(1) -es(2) 0 es(1) es(2) 0]/h;
    Dmem = Dm*t*[1 sec.nu; sec.nu 1];
    Bm = [Bh; Bs];
    ke = cw*h*(Bm'*Dmem*Bm);
    Bk = [0 0 -1 0 0 1]/h;                               % rotation gradient along the strip
    ke = ke + cw*h*(Dm*t^3/12)*(Bk'*Bk);
    dofs = [3*a-2:3*a, 3*b-2:3*b];
    [ii, jj] = ndgrid(dofs, dofs);
    Ki = [Ki; ii(:)]; Kj = [Kj; jj(:)]; Kt = [Kt; ke(:)]; %#ok<AGROW>
end
K = sparse(Ki, Kj, Kt, ndof, ndof);
f = zeros(ndof, 1); f(1:3:end) = F(:,1); f(2:3:end) = F(:,2);
fixed = 3*(supports(:,1) - 1) + supports(:,2);
free = setdiff(1:ndof, fixed);
u = zeros(ndof, 1);
u(free) = K(free, free)\f(free);
res.U = reshape(u, 3, N)';
r = K*u - f; res.reactions = [supports, r(fixed)];
res.elem = zeros(size(conn,1), 6);
for e = 1:size(conn, 1)
    a = conn(e,1); b = conn(e,2);
    dl = T_all{e}*u([3*a-2:3*a, 3*b-2:3*b]);
    fe = kloc_all{e}*dl;                                  % end forces on the element
    res.elem(e,:) = [-fe(1) fe(2) fe(3) fe(4) fe(5) fe(6)];   % N_I > 0 tension
end
res.strip = zeros(size(strips,1), 4);
for q = 1:size(strips, 1)
    a = strips(q,1); b = strips(q,2); t = strips(q,3);
    d = X(b,:) - X(a,:); h = norm(d); es = d/h; Rm = (X(a,1) + X(b,1))/2;
    ua = res.U(a,1:2); ub = res.U(b,1:2);
    eh = (ua(1) + ub(1))/(2*Rm); esv = dot(ub - ua, es)/h;
    res.strip(q,:) = [Dm*t*(eh + sec.nu*esv), Dm*t*(esv + sec.nu*eh), ...
        Dm*t^3/12*(res.U(b,3) - res.U(a,3))/h, eh];
end
end
