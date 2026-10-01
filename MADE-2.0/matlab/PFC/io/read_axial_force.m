function FZ_max = read_axial_force(xlsx_path)
%READ_AXIAL_FORCE Load the per-scenario axial (vertical) force table and reduce it to a bound.
%
%   FZ_max = READ_AXIAL_FORCE(xlsx_path) reads a per-coil, per-scenario
%   axial force file (such as FZ_PF.xlsx: one row per conductor, same row
%   convention as READ_COIL_GEOMETRY, one column per plasma scenario, in
%   MN) and returns FZ_max [N]: the maximum absolute axial force required
%   across all scenarios, one value per row.
%
%   FZ_max feeds directly into PHYSICS/EQV_STRESS_COIL_RING_CICC as the
%   coil's vertical force - it is NOT the FZmax that
%   PHYSICS/EMAG_FIELD_FORCES computes internally from the coil's own
%   geometry (see manuale PFC, "fidelity note: axial force"): the legacy
%   driver always used this externally supplied, scenario-derived value
%   for the final stress check, and that behavior is preserved unchanged.
%
%   This is a separate file from the machine parameter workbook and from
%   the geometry/Ampere-turns file (see READ_MACHINE_INPUT /
%   READ_COIL_GEOMETRY), mirroring how PF_opt_VNS.m reads it as its own
%   plain numeric table.

if ~isfile(xlsx_path)
    error('read_axial_force:file_not_found', ...
        'Axial force file not found: %s', xlsx_path);
end

raw = table2array(readtable(xlsx_path));
FZ_scenario = abs(raw);
FZ_max = max(FZ_scenario, [], 2)*1e6; % [MN] -> [N]
end
