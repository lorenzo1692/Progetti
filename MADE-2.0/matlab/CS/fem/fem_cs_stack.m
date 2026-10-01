function st = fem_cs_stack(wp, p, zc)
%FEM_CS_STACK Mesh and materials of the CS stack: n_mod modules separated by spacer plates.
%
%   st = FEM_CS_STACK(wp, p, zc) builds one structured (r,z) mesh of the
%   whole stack, centered axially on z = zc. From bottom to top every module
%   is: bottom ground insulation (wp.grins_h), turns region (n_t*fem_mesh_z
%   elements), top ground insulation; modules are separated by a steel
%   spacer plate of height wp.spacer. Radially each module is: inner ground
%   insulation (wp.grins_w), turns region (n_l*fem_mesh_r elements), outer
%   ground insulation. The spacer plates span the full radial width.
%
%   Materials: turns region = homogenized winding pack (FEM_WP_MATERIAL),
%   ground insulation = isotropic E_ins, plates = isotropic E_jckt.
%
%   st fields: mesh, D (4x4xne), matid (1 winding pack, 2 insulation,
%   3 plate), modid (module number of every winding-pack element, else 0),
%   z_turns0 (z of the bottom of the turns of each module), iz_plate (axial
%   element row of each plate), area (annulus area of each element [m^2]),
%   mat (homogenized winding-pack material).

mat = fem_wp_material(wp, p);
mr = p.fem_mesh_r; mz = p.fem_mesh_z;
gw = wp.grins_w; gh = wp.grins_h;

% radial grid lines and band of every radial element (1 ins, 2 turns, 3 ins)
nrt = wp.n_l*mr;
rv = [wp.Ri - gw, linspace(wp.Ri, wp.Re, nrt+1), wp.Re + gw];
rband = [1, 2*ones(1, nrt), 3];

% axial grid lines and band of every axial element (1 ins, 2 turns, 3 plate)
H_stack = wp.n_mod*wp.WP_h + (wp.n_mod-1)*wp.spacer;
z = zc - H_stack/2;
zv = z; aband = []; amod = [];
z_turns0 = zeros(1, wp.n_mod); iz_plate = zeros(1, wp.n_mod-1);
for m = 1:wp.n_mod
    [zv, aband, amod] = add_band(zv, aband, amod, gh, 1, 1, m);
    z_turns0(m) = zv(end);
    [zv, aband, amod] = add_band(zv, aband, amod, wp.H, wp.n_t*mz, 2, m);
    [zv, aband, amod] = add_band(zv, aband, amod, gh, 1, 1, m);
    if m < wp.n_mod
        [zv, aband, amod] = add_band(zv, aband, amod, wp.spacer, 1, 3, 0);
        iz_plate(m) = numel(aband);
    end
end

mesh = fem_mesh_grid(rv, zv);
ne = size(mesh.elems, 1);

Dwp = mat.D;
Dins = iso_D(p.E_ins*1e9, p.fem_nu_ins);
Dpl = iso_D(p.E_jckt*1e9, p.fem_nu_steel);
D = zeros(4, 4, ne);
matid = zeros(ne, 1); modid = zeros(ne, 1);
for e = 1:ne
    rb = rband(mesh.ir(e)); ab = aband(mesh.iz(e));
    if ab == 3
        matid(e) = 3; D(:, :, e) = Dpl;
    elseif ab == 2 && rb == 2
        matid(e) = 1; D(:, :, e) = Dwp; modid(e) = amod(mesh.iz(e));
    else
        matid(e) = 2; D(:, :, e) = Dins;
    end
end

st.mesh = mesh; st.D = D; st.matid = matid; st.modid = modid;
st.z_turns0 = z_turns0; st.iz_plate = iz_plate; st.mat = mat;
rvc = rv(:);
st.area = pi*(rvc(mesh.ir+1).^2 - rvc(mesh.ir).^2);
st.rv = rv; st.zv = zv; st.aband = aband; st.amod = amod;
end

function [zv, aband, amod] = add_band(zv, aband, amod, h, n, type, m)
z0 = zv(end);
zv = [zv, z0 + (1:n)*h/n];
aband = [aband, type*ones(1, n)];
amod = [amod, m*ones(1, n)];
end

function D = iso_D(E, nu)
lam = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu));
D = [lam+2*mu, lam, lam, 0; lam, lam+2*mu, lam, 0; lam, lam, lam+2*mu, 0; 0 0 0 mu];
end
