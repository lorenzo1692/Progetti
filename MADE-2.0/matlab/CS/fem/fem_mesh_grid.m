function mesh = fem_mesh_grid(rv, zv)
%FEM_MESH_GRID Structured quadrilateral mesh on arbitrary radial / axial grid lines.
%
%   mesh = FEM_MESH_GRID(rv, zv) meshes the rectangle spanned by the
%   increasing vectors rv (radial grid lines) and zv (axial grid lines) with
%   (numel(rv)-1) x (numel(zv)-1) bilinear quadrilaterals. Numbering and
%   element node order are the same as FEM_MESH_RECT (node (ir,iz) has
%   number (iz-1)*(nr+1)+ir, element (jr,jz) has number (jz-1)*nr+jr,
%   counter-clockwise from its (r-, z-) corner).
%
%   mesh fields: nodes, elems, ir, iz (element radial/axial indices), nr, nz,
%   rv, zv, bottom, top (node numbers at the first / last axial grid line).

rv = rv(:)'; zv = zv(:)';
nr = numel(rv) - 1; nz = numel(zv) - 1;
[R, Z] = meshgrid(rv, zv);
R = R'; Z = Z';
nodes = [R(:), Z(:)];

ne = nr*nz;
elems = zeros(ne, 4);
jr_all = zeros(ne, 1); jz_all = zeros(ne, 1);
e = 0;
for jz = 1:nz
    for jr = 1:nr
        e = e + 1;
        n1 = (jz-1)*(nr+1) + jr;
        elems(e, :) = [n1, n1+1, n1+1+(nr+1), n1+(nr+1)];
        jr_all(e) = jr; jz_all(e) = jz;
    end
end

mesh.nodes = nodes; mesh.elems = elems;
mesh.ir = jr_all; mesh.iz = jz_all;
mesh.nr = nr; mesh.nz = nz; mesh.rv = rv; mesh.zv = zv;
mesh.bottom = (1:nr+1)';
mesh.top = (nz*(nr+1) + (1:nr+1))';
end
