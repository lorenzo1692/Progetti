function mesh = fem_mesh_rect(r0, r1, z0, z1, nr, nz)
%FEM_MESH_RECT Structured quadrilateral mesh of a rectangular (r,z) cross-section.
%
%   mesh = FEM_MESH_RECT(r0, r1, z0, z1, nr, nz) meshes [r0,r1] x [z0,z1]
%   with nr x nz bilinear quadrilaterals. Node (ir,iz), ir=1..nr+1,
%   iz=1..nz+1, has number (iz-1)*(nr+1)+ir; element (jr,jz) has number
%   (jz-1)*nr+jr and nodes ordered counter-clockwise starting at its
%   (r-, z-) corner, as FEM_AXISYM_SOLVE expects.
%
%   mesh fields: nodes, elems, ir (ne x 1 radial element index), iz (ne x 1
%   axial element index), nr, nz, bottom (node numbers at z=z0), top
%   (node numbers at z=z1).

r = linspace(r0, r1, nr+1);
z = linspace(z0, z1, nz+1);
[R, Z] = meshgrid(r, z);
R = R'; Z = Z'; % (nr+1) x (nz+1): R(ir,iz), so that R(:) runs over ir first
nodes = [R(:), Z(:)];

ne = nr*nz;
elems = zeros(ne, 4);
jr_all = zeros(ne, 1); jz_all = zeros(ne, 1);
e = 0;
for jz = 1:nz
    for jr = 1:nr
        e = e + 1;
        n1 = (jz-1)*(nr+1) + jr;
        n2 = n1 + 1;
        n4 = n1 + (nr+1);
        n3 = n4 + 1;
        elems(e, :) = [n1 n2 n3 n4];
        jr_all(e) = jr; jz_all(e) = jz;
    end
end

mesh.nodes = nodes;
mesh.elems = elems;
mesh.ir = jr_all;
mesh.iz = jz_all;
mesh.nr = nr;
mesh.nz = nz;
mesh.bottom = (1:nr+1)';
mesh.top = (nz*(nr+1) + (1:nr+1))';
end
