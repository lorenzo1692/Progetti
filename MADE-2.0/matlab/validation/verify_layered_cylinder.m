function R = verify_layered_cylinder(p, out_file, opts)
%VERIFY_LAYERED_CYLINDER Layered-cylinder model vs the 2D FE reference results.
%
%   R = VERIFY_LAYERED_CYLINDER(p, out_file) evaluates WP_LAYERED_CYLINDER
%   on every design of the 2D FE reference data and compares:
%     - jacket Pm (primary, max over the turns of each layer) with
%       validation/results/jacket_surrogate_calibration.csv (23 runs, the
%       jacket-thickness variants rebuilt with JACKET_JT_VARIANT);
%     - case SCLs (nose centreline, plasma-side plate, lateral walls) and
%       jacket Pm with validation/results/mech_reference_fe.csv
%       (RUN_MECH_REFERENCE_FE), when present.
%   Prints per design and overall ratios model/FE and writes them to
%   out_file (text). Designs: validation/results/mech_reference_designs.mat.
%   opts: passed to WP_LAYERED_CYLINDER.

if nargin < 2, out_file = ''; end
if nargin < 3, opts = struct(); end
this_dir = fileparts(mfilename('fullpath'));
B = load(fullfile(this_dir, 'results', 'mech_reference_designs.mat'));
names = {B.base.name};
fid = 1; if ~isempty(out_file), fid = fopen(out_file, 'w'); end
say = @(varargin) both(fid, varargin{:});

% ---- jacket: calibration data set ---------------------------------------
cal = read_csv(fullfile(this_dir, 'results', 'jacket_surrogate_calibration.csv'));
runs = unique(cal.name, 'stable');
R = struct('name', {}, 'k', {}, 'Pm_FE', {}, 'Pm_LC', {});
say('Jacket Pm (primary), layered cylinder / 2D FE, per run: mean, min, max ratio; design max ratio\n');
allr = []; dmax = [];
for i = 1:numel(runs)
    nm = runs{i};
    row = design_row(nm, B, names, p);
    if isempty(row), continue, end
    sel = strcmp(cal.name, nm);
    k = cal.k(sel); PmFE = cal.Pm_P(sel);
    o = wp_layered_cylinder(row, p, opts);
    PmLC = o.Pm_jacket(k)';
    r = PmLC./PmFE;
    allr = [allr; r]; %#ok<AGROW>
    dmax(end+1) = max(PmLC)/max(PmFE); %#ok<AGROW>
    say('  %-18s %2d layers  mean %.3f  min %.3f  max %.3f | design max %.3f (FE %.0f, model %.0f MPa)\n', nm, numel(k), ...
        mean(r), min(r), max(r), dmax(end), max(PmFE)/1e6, max(PmLC)/1e6);
    R(end+1) = struct('name', nm, 'k', k', 'Pm_FE', PmFE', 'Pm_LC', PmLC'); %#ok<AGROW>
end
say('ALL: %d layers, ratio mean %.3f, rms error %.1f %%, min %.3f, max %.3f | design max ratio: min %.3f, max %.3f\n', ...
    numel(allr), mean(allr), 100*sqrt(mean((allr - 1).^2)), min(allr), max(allr), min(dmax), max(dmax));

% ---- case: reference FE runs ----------------------------------------------
ref_file = fullfile(this_dir, 'results', 'mech_reference_fe.csv');
if exist(ref_file, 'file')
    ref = read_csv(ref_file);
    rr = unique(ref.name, 'stable');
    say('\nCase SCLs (primary Pm), layered cylinder vs 2D FE [MPa]\n');
    say('  %-18s %14s %14s %14s %14s\n', 'design', 'nose centre', 'nose mid', 'plate', 'side wall 50%');
    cn = []; cp = []; cw = [];
    for i = 1:numel(rr)
        nm = rr{i}; row = design_row(nm, B, names, p);
        if isempty(row), continue, end
        o = wp_layered_cylinder(row, p, opts);
        s = strcmp(ref.name, nm) & strcmp(ref.kind, 'scl');
        lab = ref.label(s); Pm = ref.Pm_P(s);
        g = @(t) Pm(strcmp(lab, t));
        mid = round(numel(o.layer)/2);
        say('  %-18s %6.0f / %-6.0f %6.0f / %-6.0f %6.0f / %-6.0f %6.0f / %-6.0f\n', nm, g('nose centreline')/1e6, ...
            o.nose.Pm/1e6, g('nose mid')/1e6, o.nose.Pm/1e6, g('plasma-side plate')/1e6, o.plate.Pm/1e6, ...
            g('side wall 50%')/1e6, o.walls_Pm(mid)/1e6);
        cn(end+1) = o.nose.Pm/g('nose centreline'); cp(end+1) = o.plate.Pm/g('plasma-side plate'); %#ok<AGROW>
        cw(end+1) = o.walls_Pm(mid)/g('side wall 50%'); %#ok<AGROW>
    end
    say('  ratio model/FE: nose %.3f .. %.3f (mean %.3f), plate %.3f .. %.3f, side wall %.3f .. %.3f\n', ...
        min(cn), max(cn), mean(cn), min(cp), max(cp), min(cw), max(cw));
end
if fid ~= 1, fclose(fid); end
end

function row = design_row(nm, B, names, p)
% base design, or a jacket-thickness variant '<base>@jt+<dJT>' / the graded
% d10 variant '<base>@jt+0/0.5/1.0'
parts = strsplit(nm, '@');
i = find(strcmp(names, parts{1}), 1);
row = [];
if isempty(i), return, end
row = B.base(i).row;
if numel(parts) > 1
    v = sscanf(strrep(strrep(parts{2}, 'jt+', ''), '/', ' '), '%f')'*1e-3;
    if numel(v) == 1
        dJT = v;
    else                                   % one value per grade (cable change)
        tg = wp_turn_geometry(row, p);
        A = tg.A_cable; gr = cumsum([1, abs(diff(A)) > 1e-12*max(A)]);
        dJT = v(gr);
    end
    row = jacket_jt_variant(row, p, dJT);
end
end

function T = read_csv(f)
fid = fopen(f); hdr = strsplit(fgetl(fid), ',');
fmt = '';
first = fgetl(fid); frewind(fid); fgetl(fid);
vals = strsplit(first, ',');
for i = 1:numel(vals)
    if isnan(str2double(vals{i})), fmt = [fmt '%s']; else, fmt = [fmt '%f']; end %#ok<AGROW>
end
C = textscan(fid, fmt, 'Delimiter', ',');
fclose(fid);
T = struct();
for i = 1:numel(hdr), T.(hdr{i}) = C{i}; end
end

function both(fid, varargin)
fprintf(varargin{:});
if fid ~= 1, fprintf(fid, varargin{:}); end
end
