function [input_file, scenario_file, work_dir, units] = cs_select_inputs(default_template)
%CS_SELECT_INPUTS Choose the machine parameter workbook and find its scenario file next to it.
%
%   [input_file, scenario_file, work_dir] = CS_SELECT_INPUTS(default_template)
%   asks for the machine parameter workbook (default: a WP_CS_input*.xlsx in
%   the current folder, else default_template) and then looks for the plasma
%   scenario Ampere-turns file (e.g. Baseline_VNS_22_07_2025.xlsx) in the
%   SAME folder as the workbook:
%     - exactly one candidate: proposed, Enter confirms it;
%     - several: listed, choose by number;
%     - none: asks for its path, until a valid file is given.
%   Candidates are the .xlsx files of that folder named Baseline*, or, if
%   there are none, any other .xlsx that is not the workbook itself, an
%   Excel lock file (~$...) or a previous result (*_results_*, *_design_*).
%   work_dir (the folder of the workbook) is where the run should start:
%   results are saved there.
%
%   units (struct) holds the units of the scenario file, asked explicitly
%   (a wrong guess silently scales every field and force): lengths R, Z, dr,
%   dz in m or mm, Ampere-turns in A, kA or MA. A suggestion deduced from the
%   magnitude of the data is shown as the default. Fields: len_name,
%   len_scale (to m), cur_name, cur_scale (to A). Pass it to
%   READ_SCENARIO_CURRENTS and READ_COIL_GEOMETRY.

start_dir = pwd;
local = dir(fullfile(start_dir, 'WP_CS_input*.xlsx'));
local = local(~strncmp({local.name}, '~$', 2));
if ~isempty(local)
    default_input = fullfile(local(1).folder, local(1).name);
else
    default_input = default_template;
end

fprintf('Machine parameter file [%s]:\n', default_input);
answer = clean_path(input('Press Enter to use it, or type another path: ', 's'));
if isempty(answer)
    input_file = default_input;
else
    input_file = answer;
end
if ~isfile(input_file)
    error('cs_select_inputs:input_not_found', 'Machine parameter file not found: %s', input_file);
end
d = dir(input_file);
input_file = fullfile(d(1).folder, d(1).name);
work_dir = d(1).folder;

%% scenario file in the same folder
all_x = dir(fullfile(work_dir, '*.xlsx'));
names = {all_x.name};
keep = ~strncmp(names, '~$', 2) & ~strcmpi(names, d(1).name) ...
    & cellfun(@isempty, strfind(names, '_results_')) & cellfun(@isempty, strfind(names, '_design_'));
cand = names(keep);
is_base = ~cellfun(@isempty, regexpi(cand, '^baseline', 'once'));
if any(is_base), cand = cand(is_base); end

fprintf('\nScenario Ampere-turns file, looked for in %s\n', work_dir);
scenario_file = '';
if numel(cand) == 1
    fprintf('Found: %s\n', cand{1});
    answer = clean_path(input('Press Enter to use it, or type another path: ', 's'));
    if isempty(answer), scenario_file = fullfile(work_dir, cand{1}); else, scenario_file = answer; end
elseif numel(cand) > 1
    for k = 1:numel(cand), fprintf('  %d) %s\n', k, cand{k}); end
    answer = clean_path(input(sprintf('Number [1..%d] or another path: ', numel(cand)), 's'));
    n = str2double(answer);
    if isempty(answer)
        scenario_file = fullfile(work_dir, cand{1});
    elseif ~isnan(n) && n >= 1 && n <= numel(cand) && n == round(n)
        scenario_file = fullfile(work_dir, cand{n});
    else
        scenario_file = answer;
    end
else
    fprintf('No scenario file (Baseline*.xlsx) in that folder.\n');
end
while ~isfile(scenario_file)
    if ~isempty(scenario_file)
        fprintf('File not found: %s\n', scenario_file);
    end
    scenario_file = clean_path(input('Path of the scenario file: ', 's'));
end

units = ask_units(scenario_file);
end

function units = ask_units(scenario_file)
raw = table2array(readtable(scenario_file));
if size(raw, 2) < 5
    error('cs_select_inputs:unexpected_shape', ...
        'Scenario file %s has %d column(s); expected R, Z, dr, dz and one column per scenario.', scenario_file, size(raw, 2));
end
geo = raw(:, 1:4); cur = raw(:, 5:end);
len_max = max(abs(geo(:)));
cur_max = max(abs(cur(:)));
fprintf('\nScenario file contents: R, Z, dr, dz up to %.6g; Ampere-turns up to %.6g (%d rows, %d scenarios).\n', ...
    len_max, cur_max, size(raw, 1), size(raw, 2) - 4);

if len_max > 50, len_guess = 'mm'; else, len_guess = 'm'; end
if cur_max >= 1e5, cur_guess = 'A'; elseif cur_max >= 100, cur_guess = 'kA'; else, cur_guess = 'MA'; end

units.len_name = ask_choice(sprintf('Lengths (R, Z, dr, dz) are in [m / mm] (suggested %s): ', len_guess), {'m', 'mm'}, len_guess);
units.cur_name = ask_choice(sprintf('Ampere-turns are in [A / kA / MA] (suggested %s): ', cur_guess), {'A', 'kA', 'MA'}, cur_guess);
if strcmp(units.len_name, 'mm'), units.len_scale = 1e-3; else, units.len_scale = 1; end
switch units.cur_name
    case 'A',  units.cur_scale = 1;
    case 'kA', units.cur_scale = 1e3;
    case 'MA', units.cur_scale = 1e6;
end
fprintf('-> lengths in %s, Ampere-turns in %s: R, Z up to %.4g m, currents up to %.4g MA\n', units.len_name, units.cur_name, ...
    len_max*units.len_scale, cur_max*units.cur_scale*1e-6);
if len_max*units.len_scale > 30 || len_max*units.len_scale < 0.05 || cur_max*units.cur_scale > 1e9 || cur_max*units.cur_scale < 1e4
    warning('cs_select_inputs:implausible_units', ...
        'These units give lengths of %.3g m and currents of %.3g A: check them before running.', ...
        len_max*units.len_scale, cur_max*units.cur_scale);
end
end

function c = ask_choice(prompt, options, default)
c = '';
while isempty(c)
    a = strtrim(input(prompt, 's'));
    if isempty(a), a = default; end
    k = find(strcmpi(a, options), 1);
    if isempty(k)
        fprintf('  Please answer one of: %s\n', strjoin(options, ', '));
    else
        c = options{k};
    end
end
end

function s = clean_path(s)
s = strtrim(s);
s = regexprep(s, '^["'']|["'']$', '');
end
