function R = verify_tf3d_ansys_mag(ansys_dir, out_dir, cfg)
%VERIFY_TF3D_ANSYS_MAG MATLAB centreline field and forces vs the ANSYS 3D MAG model.
%
%   R = VERIFY_TF3D_ANSYS_MAG(ansys_dir, out_dir, cfg) reads the results of
%   the ANSYS 3D magnetic model MAG_360 (SOURC36 WP sources on 12 TF coils
%   + PF/CS/plasma rings, BIOT) from ansys_dir:
%     B_TF_1.csv  - x y z Bx By Bz |B| at the centreline points of coil 1
%     FL_TF_1.csv - x y z Fx Fy Fz |F| nodal forces [N] on them
%     poloidal_1_2026.csv (optional) - rc zc dr dz I(t1) I(t2) PF rings
%   rebuilds the same sources with TF3D_CENTRELINE_LOADS on the same
%   centreline points and compares field and forces point by point.
%
%   cfg (defaults = VNS 2026 input ParametriTF_VNS_LTS_09072026):
%     WPH 0.356 m, WPW 0.4368 m (mean of the high- and low-field widths, as
%     the SOURC36 real constants), NI 100 x 62.3 kA, n_TF 12, pf_time 1,
%     with_pf true, n_gauss 4.
%   Writes out_dir/verify_tf3d_ansys_mag.txt and two figures.

if nargin < 2 || isempty(out_dir), out_dir = fullfile(pwd, 'verify_tf3d_ansys'); end
if nargin < 3, cfg = struct(); end
d = struct('WPH', 0.356, 'WPW', (0.546 + 0.3276)/2, 'NI', 100*62300, 'n_TF', 12, ...
    'pf_time', 1, 'with_pf', true, 'n_gauss', 4);
fn = fieldnames(cfg); for i = 1:numel(fn), d.(fn{i}) = cfg.(fn{i}); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

Ba = read_num(fullfile(ansys_dir, 'B_TF_1.csv'));
Fa = read_num(fullfile(ansys_dir, 'FL_TF_1.csv'));
P = Ba(:,1:3);
pf = zeros(0, 5);
pf_file = fullfile(ansys_dir, 'poloidal_1_2026.csv');
if d.with_pf && exist(pf_file, 'file')
    T = read_num(pf_file);
    pf = [T(:,1:4), T(:,4+d.pf_time)];
end
opts = struct('n_gauss', d.n_gauss);
out = tf3d_centreline_loads(P, d.WPH, d.WPW, d.NI, d.n_TF, pf, opts);
% current direction: the one that reproduces the sign of the toroidal
% field on the inner leg (the ANSYS source real constants carry a sign
% convention of their own)
[~, i0] = min(P(:,2));
if sign(out.B(i0,1)) ~= sign(Ba(i0,4))
    out = tf3d_centreline_loads(P, d.WPH, d.WPW, -d.NI, d.n_TF, pf, opts);
    d.NI = -d.NI;
end
Bm = out.B; Fm = out.F;
dB = sqrt(sum((Bm - Ba(:,4:6)).^2, 2));
fid = fopen(fullfile(out_dir, 'verify_tf3d_ansys_mag.txt'), 'w');
say = @(varargin) both(fid, varargin{:});
say('MATLAB centreline loads vs ANSYS MAG_360 (%d points, PF %s, NI %.4g A)\n', size(P,1), ...
    mat2str(~isempty(pf)), d.NI);
say('|B|: ANSYS max %.4f T, MATLAB max %.4f T; max |B_M - B_A| %.4f T (%.2f %% of max |B|), rms %.4f T\n', ...
    max(Ba(:,7)), max(out.Bmag), max(dB), 100*max(dB)/max(Ba(:,7)), sqrt(mean(dB.^2)));
names = {'x (toroidal)', 'y (radial)', 'z (vertical)'};
for c = 1:3
    say('B%s: max |diff| %.4f T | F%s: ANSYS sum %.4g N, MATLAB sum %.4g N (%+.2f %%), max |diff| %.4g N\n', ...
        names{c}(1), max(abs(Bm(:,c) - Ba(:,3+c))), names{c}(1), sum(Fa(:,3+c)), sum(Fm(:,c)), ...
        100*(sum(Fm(:,c))/sum(Fa(:,3+c)) - 1), max(abs(Fm(:,c) - Fa(:,3+c))));
end
fclose(fid);
fig1 = figure('Name', 'B along the centreline'); hold on; grid on
cols = lines(3);
for c = 1:3
    plot(out.s, Ba(:,3+c), '-', 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', ['ANSYS B' names{c}(1)]);
    plot(out.s, Bm(:,c), '--', 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', ['MATLAB B' names{c}(1)]);
end
xlabel('curvilinear abscissa [m]'); ylabel('B [T]'); legend('Location', 'best'); title('Centreline field, coil 1');
fig2 = figure('Name', 'Nodal forces'); hold on; grid on
for c = 1:3
    plot(out.s, Fa(:,3+c)/1e6, '-', 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', ['ANSYS F' names{c}(1)]);
    plot(out.s, Fm(:,c)/1e6, '--', 'Color', cols(c,:), 'LineWidth', 1.5, 'DisplayName', ['MATLAB F' names{c}(1)]);
end
xlabel('curvilinear abscissa [m]'); ylabel('nodal force [MN]'); legend('Location', 'best'); title('Nodal EM forces, coil 1');
try
    print(fig1, '-dpng', '-r120', fullfile(out_dir, 'centreline_B.png'));
    print(fig2, '-dpng', '-r120', fullfile(out_dir, 'centreline_F.png'));
catch
end
R = struct('P', P, 'B_ansys', Ba(:,4:6), 'F_ansys', Fa(:,4:6), 'B_matlab', Bm, 'F_matlab', Fm, 's', out.s, 'cfg', d);
end

function A = read_num(file)
% numeric table with comma, tab or blank separators
txt = fileread(file);
txt = regexprep(txt, '[,;\t]', ' ');
A = sscanf(txt, '%f');
lines = regexp(strtrim(fileread(file)), '\r?\n', 'split');
nc = numel(sscanf(regexprep(lines{1}, '[,;\t]', ' '), '%f'));
A = reshape(A, nc, [])';
end

function both(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end
