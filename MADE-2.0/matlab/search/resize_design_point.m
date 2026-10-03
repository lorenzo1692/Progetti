function [row2, info] = resize_design_point(row, p, opts)
%RESIZE_DESIGN_POINT Size again the candidate of a design point with the current input p.
%
%   [row2, info] = RESIZE_DESIGN_POINT(row, p) runs the sizing chain of
%   SCAN_WP_DESIGNS (cables, jackets, case nose, field passes, every check)
%   on the one candidate the design point row came from: the same lateral
%   case width row.lateral_w and the same starting turns/layers combination
%   row.n_turns0. Use it when an input changes for a chosen design, for
%   example the axial load factor from the 3D global model
%   (COUPLE_GLOBAL_SIZING): the result is exactly what a full scan with the
%   new p would give for that candidate, in a few seconds.
%
%   Results saved before the column n_turns0 existed: the candidate is
%   found among the combinations of the same lateral width whose turns
%   count is one of the layer counts of the row, matching the layout
%   (n_layers and n_turns) of the row (info.matched_by = 'layout').
%
%   opts: cal (field calibration of the scan, saves its recomputation),
%   g, env (as computed by the main, optional).
%   row2: the re-sized design point (same type as a DATA row), [] if the
%   candidate is now rejected. info: ok, matched_by, n_found, elapsed.

if nargin < 3, opts = struct(); end
t0 = tic;
r = row; if ~isstruct(r), r = table2struct(r); end
if isfield(opts, 'g') && ~isempty(opts.g), g = opts.g; else, g = compute_operating_params(p); end
if isfield(opts, 'env') && ~isempty(opts.env), env = opts.env; else, env = wp_envelope(p, g); end
env.lateral_w_min = r.lateral_w; env.lateral_w_max = r.lateral_w;
combT_all = generate_combinations(env);
combT = cell(size(combT_all));
nl = r.n_layers; nt_row = r.n_turns(1:nl);
if isfield(r, 'n_turns0') && any(r.n_turns0 > 0)
    nt0 = r.n_turns0(r.n_turns0 > 0);
    matched_by = 'n_turns0';
    for i = 1:numel(combT_all)
        c = combT_all{i};
        if size(c, 2) == numel(nt0)
            combT{i} = c(all(c == nt0(:)', 2), :);
        else
            combT{i} = zeros(0, size(c, 2));
        end
    end
else
    matched_by = 'layout';
    for i = 1:numel(combT_all)
        c = combT_all{i};
        combT{i} = c(ismember(c(:, 1), unique(nt_row)), :);
    end
end
sopts = struct('quiet', true);
if isfield(opts, 'cal'), sopts.cal = opts.cal; end
D = scan_wp_designs(p, g, env, combT, sopts);
n = height_of(D);
row2 = []; k = [];
if n > 0
    if strcmp(matched_by, 'n_turns0')
        k = 1;
    else
        for i = 1:n
            ri = pick(D, i);
            if ri.n_layers == nl && isequal(ri.n_turns(1:nl), nt_row), k = i; break, end
        end
    end
    if ~isempty(k), row2 = D(k, :); end
end
info = struct('ok', ~isempty(row2), 'matched_by', matched_by, 'n_found', n, 'elapsed', toc(t0));
end

function n = height_of(D)
if isstruct(D), n = numel(D); elseif isempty(D), n = 0; else, n = height(D); end
end

function r = pick(D, i)
if isstruct(D), r = D(i); else, r = table2struct(D(i, :)); end
end
