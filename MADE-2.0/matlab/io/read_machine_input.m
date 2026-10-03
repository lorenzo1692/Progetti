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

% Raw cells, not READTABLE: on some MATLAB versions / workbooks the table
% import returned a "Name" variable shorter than the table height ("Index
% exceeds the number of array elements" in tabular/dotParenReference).
C = read_raw_cells(xlsx_path);
if isempty(C)
    error('read_machine_input:empty_sheet', 'Sheet "Input" of %s is empty.', xlsx_path);
end
hdr = cellfun(@cell_text, C(1,:), 'UniformOutput', false);
iN = find(strcmp(hdr, 'Name'), 1); iV = find(strcmp(hdr, 'Value'), 1);
missing_cols = setdiff({'Name', 'Value'}, hdr(~cellfun(@isempty, hdr)));
if isempty(iN) || isempty(iV)
    error('read_machine_input:missing_columns', ...
        'Input sheet is missing required column(s): %s', strjoin(missing_cols, ', '));
end

p = struct();
for k = 2:size(C, 1)
    name = cell_text(C{k, iN});
    if isempty(name)
        continue
    end
    field = matlab.lang.makeValidName(name);
    value = cell_number(C{k, iV});
    if isnan(value)
        error('read_machine_input:invalid_value', ...
            'Parameter "%s" (row %d) has a non-numeric or missing Value in %s.', name, k, xlsx_path);
    end
    if isfield(p, field)
        error('read_machine_input:duplicate_parameter', ...
            'Parameter "%s" (row %d) appears more than once in %s.', name, k, xlsx_path);
    end
    p.(field) = value;
end

check_required_parameters(p, xlsx_path);
end

function check_required_parameters(p, xlsx_path)
required = {'n_TF','B0','R0','A','wb','ripple','dr_plasma_side','R_VV_margin_factor', ...
    'E_jckt','E_cbl_HTS','E_cbl_LTS','E_case','E_ins', ...
    'WP_SC_type','shape_cable','corr_B_WP', ...
    'S_amm_JT','S_amm_VT','S_amm_Cu','safety_membrane','S_VV','THS_max_LTS','THS_max_HTS', ...
    'V_MAX','toroidal_gap','GoundIns','INS_grades','turn_insulation_nominal', ...
    'r_SC_min','r_SC_max','min_SC_w','min_cable_aspect_ratio', ...
    'maxdim','Increm','n_grades','grade_B_target_2','grade_B_target_3', ...
    'Iop_min','Iop_max','min_size_CICC','max_size_CICC','min_JT','JT_step', ...
    'DTF_initial','DTF_step', ...
    'WP_radial_build_max_est','case_wedge_frac_min','case_wedge_frac_max', ...
    'lateral_w_step','max_sizing_iterations','SCF_transition_provisional'};

missing = required(~isfield(p, required));
if ~isempty(missing)
    error('read_machine_input:missing_parameters', ...
        ['The following required parameters are missing from %s:\n  %s\n', ...
         'Check them against MADE-2.0/matlab/input/WP_TF_input_template.xlsx.'], ...
        xlsx_path, strjoin(missing, ', '));
end
end

function C = read_raw_cells(xlsx_path)
% all cells of the "Input" sheet as a cell array (text, numbers, empties)
if exist('readcell', 'file') || exist('readcell', 'builtin')
    C = readcell(xlsx_path, 'Sheet', 'Input');
else
    [~, ~, C] = xlsread(xlsx_path, 'Input');   %#ok<XLSRD> MATLAB before R2019a
end
end

function t = cell_text(x)
% trimmed text of a cell ('' for empty / missing / numeric cells)
if ischar(x) || isstring(x)
    t = strtrim(char(x));
    if strcmpi(t, '<missing>'), t = ''; end
else
    t = '';
end
end

function v = cell_number(x)
% numeric value of a cell (NaN for empty / missing / non-numeric text)
if isnumeric(x) && isscalar(x)
    v = double(x);
elseif islogical(x) && isscalar(x)
    v = double(x);
elseif ischar(x) || isstring(x)
    v = str2double(x);
else
    v = NaN;
end
end
