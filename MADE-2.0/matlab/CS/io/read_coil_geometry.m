function geom = read_coil_geometry(xlsx_path)
%READ_COIL_GEOMETRY Load the combined CS+PFC+plasma geometry and per-scenario Ampere-turns table.
%
%   geom = READ_COIL_GEOMETRY(xlsx_path) reads a machine geometry file
%   (columns: R, Z, dr, dz, then one column per plasma scenario) such as
%   Baseline_VNS_07_2026_V3_CREATE.xlsx, and returns a struct geom with:
%     geom.R, geom.Z, geom.dr, geom.dz  - one value per conductor row
%     geom.MAt_signed                  - scenario currents with their sign
%                                         (needed for the field of the other
%                                         coils), same layout as MAt_scenario
%     geom.MAt_scenario                - abs(scenario currents), one row
%                                         per conductor, one column per
%                                         scenario
%     geom.MAt                         - max(MAt_scenario, [], 2): the
%                                         governing (worst-case) Ampere-
%                                         turns requirement per conductor
%
%   Row convention: rows 1:n_moduli are the CS modules; any further row (PF
%   coils, plasma filament) is a background conductor of the FEM field.
%   Copy of MADE-2.0/matlab/PFC/io/read_coil_geometry.m, used by
%   FEM_CS_VERIFY (the CS main reads only the Ampere-turns bound from the
%   same file through READ_SCENARIO_CURRENTS).

if ~isfile(xlsx_path)
    error('read_coil_geometry:file_not_found', ...
        'Geometry file not found: %s', xlsx_path);
end

raw = table2array(readtable(xlsx_path));
if size(raw, 2) < 5
    error('read_coil_geometry:unexpected_shape', ...
        ['Geometry file %s has %d column(s); expected at least 5 ' ...
         '(R, Z, dr, dz, then one column per plasma scenario).'], xlsx_path, size(raw,2));
end

geom = struct();
geom.R = raw(:, 1);
geom.Z = raw(:, 2);
geom.dr = raw(:, 3);
geom.dz = raw(:, 4);
geom.MAt_signed = raw(:, 5:end);
geom.MAt_scenario = abs(geom.MAt_signed);
geom.MAt = max(geom.MAt_scenario, [], 2);
end
