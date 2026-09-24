function [row, p] = load_design_point(results_file, idx, input_file)
%LOAD_DESIGN_POINT Reload one design point saved by main_WP_TF_design.m.
%
%   [row, p] = LOAD_DESIGN_POINT(results_file, idx, input_file) points
%   PLOT_WP_SECTION / EXPORT_ANSYS_INPUT at a solution saved earlier,
%   without re-running the scan:
%
%   - If results_file has a sibling .mat with the same base name (saved
%     automatically by main_WP_TF_design.m alongside the .xlsx), it is
%     loaded directly: DATA and p come back exactly as MATLAB produced
%     them, arrays included. input_file is then optional (only needed to
%     use a DIFFERENT machine input than the one that produced this run).
%   - Otherwise, results_file (the *_results_*.xlsx itself) is read with
%     READTABLE and its array-valued columns (n_turns, Cond_w, Cond_h, JT,
%     r_cable, N_Sc, N_Cu, S_Cable, S_REBCO, S_Cu_HTS, THS, type_cable),
%     which MATLAB's WRITETABLE spread across "<name>_1".."<name>_k"
%     spreadsheet columns, are reassembled into arrays. input_file is then
%     REQUIRED, since the xlsx alone does not carry the machine parameters.
%
%   idx is the row/design-point index (as shown by BROWSE_SOLUTIONS).

[base_dir, base_name] = fileparts(results_file);
mat_file = fullfile(base_dir, [base_name '.mat']);

if isfile(mat_file)
    S = load(mat_file, 'DATA', 'p');
    row = S.DATA(idx, :);
    if nargin >= 3 && ~isempty(input_file)
        p = read_machine_input(input_file); % explicit override
    else
        p = S.p;
    end
    return
end

if nargin < 3 || isempty(input_file)
    error('load_design_point:input_file_required', ...
        ['No %s found next to %s, so the machine parameters cannot be recovered ' ...
         'from the results file alone: pass the machine input xlsx as the third argument.'], ...
        [base_name '.mat'], results_file);
end
p = read_machine_input(input_file);

T = readtable(results_file, 'VariableNamingRule', 'preserve');
array_cols = {'n_turns','Cond_w','Cond_h','JT','r_cable','N_Sc','N_Cu','S_Cable','S_REBCO','S_Cu_HTS','THS','B_grade','type_cable'};
row = struct();
for c = 1:width(T)
    name = T.Properties.VariableNames{c};
    is_array_piece = false;
    for a = 1:numel(array_cols)
        prefix = [array_cols{a} '_'];
        if startsWith(name, prefix) && ~isnan(str2double(extractAfter(name, prefix)))
            is_array_piece = true;
            break
        end
    end
    if ~is_array_piece
        row.(name) = T.(name)(idx);
    end
end
for a = 1:numel(array_cols)
    base = array_cols{a};
    pieces = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, [base '_']));
    if isempty(pieces)
        continue
    end
    piece_idx = cellfun(@(s) str2double(extractAfter(s, [base '_'])), pieces);
    [~, order] = sort(piece_idx);
    pieces = pieces(order);
    if strcmp(base, 'type_cable')
        row.(base) = cellfun(@(v) T.(v){idx}, pieces, 'UniformOutput', false);
    else
        row.(base) = cellfun(@(v) T.(v)(idx), pieces);
    end
end
end
