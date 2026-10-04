function p = read_machine_input(xlsx_path)
%READ_MACHINE_INPUT Load the machine/search parameters from the input Excel file.
%
%   p = READ_MACHINE_INPUT(xlsx_path) reads the "Input" sheet of xlsx_path
%   (columns: Category, Name, Value, Unit, Description) and returns a
%   struct p with one field per parameter Name, holding its Value. This is
%   the single entry point through which every independent parameter of
%   the WP_PFC solver is made available to the rest of the code.
%
%   The workbook is meant to be reviewed and edited by hand before every
%   run (see MADE-2.0/matlab/PFC/input/WP_PFC_input_template.xlsx). Do not
%   rename or remove rows: this function matches parameters by the exact
%   text in the "Name" column.
%
%   Mirrors MADE-2.0/matlab/CS/io/read_machine_input.m (the CS version).
%   The per-coil geometry/scenario tables (CS+PFC+plasma R/Z/dr/dz and
%   Ampere-turns, and the axial-force scenario) are separate files, read
%   by READ_COIL_GEOMETRY / READ_AXIAL_FORCE instead.
%
%   The "FEM verification" parameters (fem_*) are optional: a workbook
%   without them is completed with the defaults of FEM/FEM_DEFAULTS when
%   the FEM module is used.

if ~isfile(xlsx_path)
    error('read_machine_input:file_not_found', ...
        'Input file not found: %s', xlsx_path);
end

T = readtable(xlsx_path, 'Sheet', 'Input', 'TextType', 'string');

required_cols = ["Name", "Value"];
missing_cols = setdiff(required_cols, string(T.Properties.VariableNames));
if ~isempty(missing_cols)
    error('read_machine_input:missing_columns', ...
        'Input sheet is missing required column(s): %s', strjoin(missing_cols, ', '));
end

p = struct();
for k = 1:height(T)
    name = strtrim(T.Name(k));
    if name == "" || ismissing(name)
        continue
    end
    field = matlab.lang.makeValidName(name);
    value = T.Value(k);
    if ~isnumeric(value) || isnan(value)
        error('read_machine_input:invalid_value', ...
            'Parameter "%s" has a non-numeric or missing Value in %s.', name, xlsx_path);
    end
    if isfield(p, field)
        error('read_machine_input:duplicate_parameter', ...
            'Parameter "%s" appears more than once in %s.', name, xlsx_path);
    end
    p.(field) = value;
end

check_required_parameters(p, xlsx_path);
end

function check_required_parameters(p, xlsx_path)
required = ["cs_modules","n_moduli","WP_SC_type","shape_cable","Increm", ...
    "grins_w","grins_h","tins","ins_grades","spacer", ...
    "n_grades","WP_h_min","WP_h_max","WP_h_step", ...
    "min_size_CICC","max_size_CICC","Iop_min","Iop_max","Iop_step", ...
    "min_B","max_B","MAt_tolerance", ...
    "S_hoop_allow","SF_hoop","Tresca_factor", ...
    "E_jckt","E_cbl_LTS","E_cbl_HTS","E_ins","ring_nu","ring_Fz_area_fix","fz_source", ...
    "fcgr_mode","plasma_cycles_min","fcgr_C0","fcgr_m","fcgr_mw","fcgr_residual_stress", ...
    "fcgr_KIC","fcgr_flaw_aspect","fcgr_flaw_type","fcgr_flaw_area", ...
    "JT_min","JT_step","max_sizing_iterations","Tau_discharge_fixed", ...
    "SC_w_min","r_cable_min","r_cable_max","WP_radial_build_max","B_background","system_sizing"];

missing = required(~isfield(p, required));
if ~isempty(missing)
    error('read_machine_input:missing_parameters', ...
        ['The following required parameters are missing from %s:\n  %s\n', ...
         'Check them against MADE-2.0/matlab/PFC/input/WP_PFC_input_template.xlsx.'], ...
        xlsx_path, strjoin(missing, ', '));
end
end
