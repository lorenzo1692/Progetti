function p = read_machine_input(xlsx_path)
%READ_MACHINE_INPUT Load the machine/search parameters from the input Excel file.
%
%   p = READ_MACHINE_INPUT(xlsx_path) reads the "Input" sheet of xlsx_path
%   (columns: Category, Name, Value, Unit, Description) and returns a
%   struct p with one field per parameter Name, holding its Value. This is
%   the single entry point through which every independent parameter of
%   the WP_CS solver is made available to the rest of the code.
%
%   The workbook is meant to be reviewed and edited by hand before every
%   run (see MADE-2.0/matlab/CS/input/WP_CS_input_template.xlsx). Do not
%   rename or remove rows: this function matches parameters by the exact
%   text in the "Name" column.
%
%   Mirrors MADE-2.0/matlab/io/read_machine_input.m (the TF version); the
%   machine scenario table (Ampere-turns per plasma scenario) is a
%   separate file, read by READ_SCENARIO_CURRENTS instead.

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
required = ["n_moduli","WP_h0","Re_0","grins_w","grins_h","tins","spacer","r_SC", ...
    "WP_SC_type","shape_cable","V_MAX","SF_V","Increm", ...
    "n_grades","min_size_CICC","max_size_CICC","Iop_min","Iop_max","Iop_steps", ...
    "min_B","max_B","use_pancake","JT_step","JT_min","max_sizing_iterations", ...
    "E_jckt","E_cbl_HTS","E_cbl_LTS","E_ins", ...
    "S_hoop_allow","SF_hoop","S_T_allow","SF_T","F_z_check_MN", ...
    "SC_w_min","r_cable_min","r_cable_max","JT_min_check","Ri_min"];

missing = required(~isfield(p, required));
if ~isempty(missing)
    error('read_machine_input:missing_parameters', ...
        ['The following required parameters are missing from %s:\n  %s\n', ...
         'Check them against MADE-2.0/matlab/CS/input/WP_CS_input_template.xlsx.'], ...
        xlsx_path, strjoin(missing, ', '));
end
end
