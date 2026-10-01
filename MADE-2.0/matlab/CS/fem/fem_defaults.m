function p = fem_defaults(p)
%FEM_DEFAULTS Complete the machine parameter struct with the CS FEM defaults.
%
%   p = FEM_DEFAULTS(p) adds every FEM parameter that the CS input workbook
%   does not define, so the FEM module works with the current template.
%   Values already in p are never overwritten.
%
%   fem_mesh_r, fem_mesh_z  elements per conductor cell, radial / axial
%   fem_subfil              n x n current sub-filaments per turn (self field)
%   fem_bg_fil_r/z          filaments per background coil block (PF, plasma)
%   fem_bc                  1 inertia relief, 2 stack bottom face fixed axially
%   fem_nu_wp/ins/steel     Poisson ratios: winding pack, insulation, spacer plates
%   fem_preload_MN          stack compression applied on the top face (needs fem_bc=2)

d = struct( ...
    'fem_mesh_r', 2, 'fem_mesh_z', 2, 'fem_subfil', 3, ...
    'fem_bg_fil_r', 4, 'fem_bg_fil_z', 8, 'fem_bc', 1, ...
    'fem_nu_wp', 0.3, 'fem_nu_ins', 0.3, 'fem_nu_steel', 0.3, ...
    'fem_preload_MN', 0, ...
    'E_jckt', 205, 'E_cbl_LTS', 0.1, 'E_cbl_HTS', 120, 'E_ins', 20);
names = fieldnames(d);
for k = 1:numel(names)
    if ~isfield(p, names{k})
        p.(names{k}) = d.(names{k});
    end
end
end
