function MAt = read_scenario_currents(xlsx_path, units)
%READ_SCENARIO_CURRENTS Load the per-scenario Ampere-turns table and reduce it to a bound.
%
%   MAt = READ_SCENARIO_CURRENTS(xlsx_path) reads a machine scenario file
%   (columns: R, Z, dr, dz, then one column per plasma scenario) such as
%   Baseline_VNS_22_07_2025.xlsx or Baseline_CEFTR-10-12-2025.xlsx, and
%   returns MAt: the maximum absolute Ampere-turns required across all
%   scenarios, one value per row of the table, in A. units (optional struct
%   from CS_SELECT_INPUTS) gives the unit of the file's Ampere-turns
%   (units.cur_scale to A); default A.
%
%   This is a separate file from the machine parameter workbook (see
%   READ_MACHINE_INPUT / WP_CS_input_template.xlsx): the CS legacy driver
%   (CS_opt_VNS.m) reads a plain numeric table here, not a Name/Value
%   sheet, so it is kept as its own simple loader rather than forced into
%   the Excel-parameter convention.
%
%   GENERATE_COMBINATIONS uses max(MAt(1:6)) to bound the feasible number
%   of layers - the "first 6 scenarios" convention is carried over
%   unchanged from the legacy driver (see manuale CS).

if ~isfile(xlsx_path)
    error('read_scenario_currents:file_not_found', ...
        'Scenario file not found: %s', xlsx_path);
end

raw = table2array(readtable(xlsx_path));
if size(raw, 2) < 6
    error('read_scenario_currents:unexpected_shape', ...
        ['Scenario file %s has %d column(s); expected at least 6 ' ...
         '(R, Z, dr, dz, then per-scenario currents).'], xlsx_path, size(raw,2));
end

cur_scale = 1;
if nargin > 1 && isfield(units, 'cur_scale'), cur_scale = units.cur_scale; end
MAt_scenario = abs(raw(:, 5:end))*cur_scale;
MAt = max(MAt_scenario, [], 2);
end
