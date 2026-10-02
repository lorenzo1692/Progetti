function st = fem_cs_stack(wp, p, zc)
%FEM_CS_STACK Mesh and materials of the CS stack: n_mod modules separated by spacer plates.
%
%   st = FEM_CS_STACK(wp, p, zc) builds one structured (r,z) mesh of the
%   whole stack, centered axially on z = zc. From bottom to top every module
%   is: bottom ground insulation (wp.grins_h), turns region, top ground
%   insulation; modules are separated by a steel spacer plate of height
%   wp.spacer. Radially each module is: inner ground insulation
%   (wp.grins_w), turns region, outer ground insulation. The plates span
%   the full radial width.
%
%   Turns region: homogenized winding pack (FEM_WP_MATERIAL, n_l*fem_mesh_r
%   by n_t*fem_mesh_z elements per module). If p.fem_detail_module = m > 0,
%   module m is instead modeled turn by turn: every conductor cell is
%   turn insulation (tins) - jacket wall (JT) - cable (SC_w x SC_h, with
%   p.fem_detail_cable_r x p.fem_detail_cable_z elements) - jacket wall -
%   turn insulation, as a rectangular cell (the corner radius of the
%   jacket is neglected). Because the mesh is a tensor grid, the radial
%   grid lines of the detailed cells run through all modules (the other
%   modules stay homogenized, with elements of the same radial size);
%   axial refinement stays inside module m.
%
%   Materials: ground insulation / turn insulation = isotropic E_ins,
%   plates and jacket = isotropic E_jckt, cable = isotropic E_cbl (LTS or
%   HTS) in the detailed module.
%
%   st fields: mesh, D (4x4xne), matid (1 homogenized winding pack,
%   2 insulation, 3 plate, 4 cable, 5 jacket), modid (module number of the
%   elements of the turns regions, else 0), z_turns0 (z of the bottom of
%   the turns of each module), iz_plate (axial element row of each plate),
%   area (annulus area of each element [m^2]), mat (homogenized material),
%   detail (detailed module, 0 if none), lay, tur (layer / turn index of
%   every turns-region element, else 0).

mat = fem_wp_material(wp, p);
mr = p.fem_mesh_r; mz = p.fem_mesh_z;
gw = wp.grins_w; gh = wp.grins_h;
md = p.fem_detail_module;
if md > wp.n_mod || md < 0
    error('fem_cs_stack:bad_detail_module', 'fem_detail_module must be between 0 and %d.', wp.n_mod);
end

%% radial grid: rcode 0 ground insulation; in the turns zone 1 ins, 2 jacket, 3 cable, 4 homogenized
if md > 0
    ncr = p.fem_detail_cable_r;
    one = [wp.tins, wp.JT, wp.SC_w/ncr*ones(1, ncr), wp.JT, wp.tins];
    code1 = [1, 2, 3*ones(1, ncr), 2, 1];
    dr = repmat(one, 1, wp.n_l); rcode_t = repmat(code1, 1, wp.n_l);
    rlay = repelem(1:wp.n_l, numel(one));
else
    dr = wp.Cond_w/mr*ones(1, wp.n_l*mr); rcode_t = 4*ones(1, wp.n_l*mr);
    rlay = repelem(1:wp.n_l, mr);
end
rv = [wp.Ri - gw, wp.Ri, wp.Ri + cumsum(dr), wp.Re + gw];
rcode = [0, rcode_t, 0];
rlayer = [0, rlay, 0];

%% axial grid: acode 1 ground ins, 2 turns (homogenized), 3 plate, 4 turns (detailed cell sub-band codes in asub)
H_stack = wp.n_mod*wp.WP_h + (wp.n_mod-1)*wp.spacer;
zv = zc - H_stack/2;
acode = []; amod = []; asub = []; aturn = [];
z_turns0 = zeros(1, wp.n_mod); iz_plate = zeros(1, wp.n_mod-1);
for m = 1:wp.n_mod
    [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, gh, 1, 1, m, 0, 0);
    z_turns0(m) = zv(end);
    if m == md
        nz_c = p.fem_detail_cable_z;
        one = [wp.tins, wp.JT, wp.SC_h/nz_c*ones(1, nz_c), wp.JT, wp.tins];
        code1 = [1, 2, 3*ones(1, nz_c), 2, 1];
        for t = 1:wp.n_t
            for k = 1:numel(one)
                [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, one(k), 1, 4, m, code1(k), t);
            end
        end
    else
        for t = 1:wp.n_t
            [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, wp.Cond_h, mz, 2, m, 0, t);
        end
    end
    [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, gh, 1, 1, m, 0, 0);
    if m < wp.n_mod
        [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, wp.spacer, 1, 3, 0, 0, 0);
        iz_plate(m) = numel(acode);
    end
end

mesh = fem_mesh_grid(rv, zv);
ne = size(mesh.elems, 1);

Ec = p.E_cbl_LTS*1e9;
if strcmp(wp.type_cable, 'HTS'), Ec = p.E_cbl_HTS*1e9; end
Dwp = mat.D;
Dins = iso_D(p.E_ins*1e9, p.fem_nu_ins);
Dpl = iso_D(p.E_jckt*1e9, p.fem_nu_steel);
Dcab = iso_D(Ec, p.fem_nu_cable);
D = zeros(4, 4, ne);
matid = zeros(ne, 1); modid = zeros(ne, 1); lay = zeros(ne, 1); tur = zeros(ne, 1);
for e = 1:ne
    rc = rcode(mesh.ir(e)); iz = mesh.iz(e); ac = acode(iz);
    if ac == 3
        matid(e) = 3;
    elseif ac == 1 || rc == 0
        matid(e) = 2;
    else
        modid(e) = amod(iz); lay(e) = rlayer(mesh.ir(e)); tur(e) = aturn(iz);
        if ac == 2
            matid(e) = 1;                                   % homogenized winding pack
        else                                                % detailed cell
            sa = asub(iz);
            if sa == 1 || rc == 1,    matid(e) = 2;
            elseif sa == 3 && rc == 3, matid(e) = 4;
            else,                      matid(e) = 5;
            end
        end
    end
    switch matid(e)
        case 1, D(:, :, e) = Dwp;
        case 2, D(:, :, e) = Dins;
        case 3, D(:, :, e) = Dpl;
        case 4, D(:, :, e) = Dcab;
        case 5, D(:, :, e) = Dpl;
    end
end

st.mesh = mesh; st.D = D; st.matid = matid; st.modid = modid; st.lay = lay; st.tur = tur;
st.z_turns0 = z_turns0; st.iz_plate = iz_plate; st.mat = mat; st.detail = md;
rvc = rv(:);
st.area = pi*(rvc(mesh.ir+1).^2 - rvc(mesh.ir).^2);
st.rv = rv; st.zv = zv; st.aband = acode; st.amod = amod;
end

function [zv, acode, amod, asub, aturn] = add_band(zv, acode, amod, asub, aturn, h, n, type, m, sub, turn)
z0 = zv(end);
zv = [zv, z0 + (1:n)*h/n];
acode = [acode, type*ones(1, n)];
amod = [amod, m*ones(1, n)];
asub = [asub, sub*ones(1, n)];
aturn = [aturn, turn*ones(1, n)];
end

function D = iso_D(E, nu)
lam = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu));
D = [lam+2*mu, lam, lam, 0; lam, lam+2*mu, lam, 0; lam, lam, lam+2*mu, 0; 0 0 0 mu];
end
