function p = fem_defaults(p)
%FEM_DEFAULTS Complete the machine parameter struct with the FEM module defaults.
%
%   p = FEM_DEFAULTS(p) adds every parameter of the "FEM verification" and
%   material categories that the input workbook does not define (older
%   workbooks have none of the fem_* rows), so the FEM module can be used
%   with any PFC input file. Values already in p are never overwritten.

d = struct( ...
    'fem_mesh_r', 2, 'fem_mesh_z', 2, 'fem_subfil', 4, ...
    'fem_bg_fil_r', 4, 'fem_bg_fil_z', 8, 'fem_bc', 1, 'fem_nu_wp', 0.3, ...
    'E_jckt', 205, 'E_cbl_LTS', 0.1, 'E_cbl_HTS', 120, 'E_ins', 20, ...
    'ring_nu', 1/3, 'ring_Fz_area_fix', 0, 'fz_source', 0, 'tins', 0.001);
names = fieldnames(d);
for k = 1:numel(names)
    if ~isfield(p, names{k})
        p.(names{k}) = d.(names{k});
    end
end
end
