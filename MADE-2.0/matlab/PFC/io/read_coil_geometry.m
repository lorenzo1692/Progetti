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
%   Row convention (unchanged from the legacy driver): rows
%   1:p.cs_modules are the CS modules, rows p.cs_modules+1:p.cs_modules+6
%   are PF1..PF6, and any further row is the plasma filament - the same
%   file COMPUTE_COUPLING_MATRIX reads to build the mutual inductance
%   matrix (it needs R, Z, dr, dz for every conductor, this function's own
%   geom.R/Z/dr/dz).
%
%   This is a separate file from the machine parameter workbook (see
%   READ_MACHINE_INPUT / WP_PFC_input_template.xlsx): the PF legacy driver
%   (PF_opt_VNS.m) reads a plain numeric table here, not a Name/Value
%   sheet, so it is kept as its own simple loader, mirroring
%   MADE-2.0/matlab/CS/io/read_scenario_currents.m.

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
