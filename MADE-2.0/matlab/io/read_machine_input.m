function p = read_machine_input(xlsx_path)
%READ_MACHINE_INPUT Load the machine/search parameters from the input Excel file.
%
%   p = READ_MACHINE_INPUT(xlsx_path) reads the "Input" sheet of xlsx_path
%   (columns: Category, Name, Value, Unit, Description) and returns a
%   struct p with one field per parameter Name, holding its Value. This is
%   the single entry point through which every independent parameter of
%   the WP_TF solver is made available to the rest of the code.
%
%   The workbook is meant to be reviewed and edited by hand before every
%   run (see MADE-2.0/matlab/input/WP_TF_input_template.xlsx). Do not
%   rename or remove rows: this function matches parameters by the exact
%   text in the "Name" column.

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
required = ["n_TF","B0","R0","A","wb","ripple","dr_plasma_side", ...
    "E_jckt","E_cbl_HTS","E_cbl_LTS","E_case","E_ins", ...
    "WP_SC_type","shape_cable","corr_B_WP", ...
    "S_amm_JT","S_amm_VT","S_amm_Cu","safety_membrane","S_VV", ...
    "V_MAX","toroidal_gap","GoundIns","INS_grades", ...
    "maxdim","Increm","n_grades","grade_B_target_2","grade_B_target_3", ...
    "Iop_min","Iop_max","min_size_CICC","max_size_CICC","min_JT", ...
    "WP_radial_build_max_est","case_wedge_frac_min","case_wedge_frac_max", ...
    "lateral_w_step","max_sizing_iterations","SCF_transition_provisional"];

missing = required(~isfield(p, required));
if ~isempty(missing)
    error('read_machine_input:missing_parameters', ...
        ['The following required parameters are missing from %s:\n  %s\n', ...
         'Check them against MADE-2.0/matlab/input/WP_TF_input_template.xlsx.'], ...
        xlsx_path, strjoin(missing, ', '));
end
end
