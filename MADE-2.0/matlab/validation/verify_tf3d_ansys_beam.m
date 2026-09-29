function R = verify_tf3d_ansys_beam(ansys_dir, out_dir, cfg)
%VERIFY_TF3D_ANSYS_BEAM MATLAB in-plane beam model vs the ANSYS global beam model STR_360.
%
%   R = VERIFY_TF3D_ANSYS_BEAM(ansys_dir, out_dir, cfg) reads from ansys_dir
%   the export of validation/ansys/POST_TF_BEAM_EXPORT.mac (BEAM_TF1_EL.csv,
%   NODE_TF1_U.csv, SHELL_TF1_EL.csv) and the loads FL_TF_1.csv, rebuilds
%   the same structural model with TF3D_BEAM_INPLANE (coil 1, cyclic
%   symmetry) and compares axial force, in-plane bending moment, nodal
%   displacements and strip membrane forces.
%
%   Model as STR_360.dat (run of 29/09/2026): BEAM188 HREC section, width
%   cfg.W = 0.675 m, wall cfg.t_wall = 0.05 m, E 205 GPa, nu 0.3; SHELL181
%   strips cfg.t_shell = 0.14 m on every segment of the straight inner leg
%   (vault) and on the outer-leg nodes with R in [3.5, 5] m, upper and lower
%   groups separately, each closed by a strip between its first and last
%   node (as in STR_360.dat; STR_360_improved removes that closure,
%   cfg.ois_wrap = false); loads FL_TF_1 row j applied to node j (radial
%   and vertical components, the toroidal one is not applied in STR_360);
%   vertical support at the gravity-support node (cfg.gs_node, default the
%   node of coil 1 listed first in REACT_TF.txt).

if nargin < 2 || isempty(out_dir), out_dir = fullfile(pwd, 'verify_tf3d_beam'); end
if nargin < 3, cfg = struct(); end
d = struct('W', 0.675, 't_wall', 0.05, 'E', 205e9, 'nu', 0.3, 't_shell', 0.14, 'n_TF', 12, ...
    'ois_R', [3.5 5.0], 'ois_wrap', true, 'gs_node', []);
fn = fieldnames(cfg); for i = 1:numel(fn), d.(fn{i}) = cfg.(fn{i}); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

Bm = read_num(fullfile(ansys_dir, 'BEAM_TF1_EL.csv'));
Nd = read_num(fullfile(ansys_dir, 'NODE_TF1_U.csv'));
Sh = read_num(fullfile(ansys_dir, 'SHELL_TF1_EL.csv'));
FL = read_num(fullfile(ansys_dir, 'FL_TF_1.csv'));
id = Nd(:,1); X = [hypot(Nd(:,2), Nd(:,3)), Nd(:,4)];
nn = numel(id);
% beam elements from the end-node coordinates
conn = zeros(size(Bm,1), 2);
for e = 1:size(Bm,1)
    [~, conn(e,1)] = min(sum((Nd(:,2:4) - Bm(e,2:4)).^2, 2));
    [~, conn(e,2)] = min(sum((Nd(:,2:4) - Bm(e,5:7)).^2, 2));
end
% section: hollow square (BEAM188 HREC)
W = d.W; ti = W - 2*d.t_wall;
sec = struct('E', d.E, 'nu', d.nu, 'A', W^2 - ti^2, 'I', (W^4 - ti^4)/12, 'As', 2*(W - d.t_wall)*d.t_wall);
% loads: row j of FL_TF_1 on node j (node order), radial = global Y for coil 1
[~, ord] = sort(id);
F = zeros(nn, 2); m = min(nn, size(FL,1));
F(ord(1:m), :) = [FL(1:m,5), FL(1:m,6)];
% strips
tol = 1e-3;
Rmin = min(X(:,1));
il = find(abs(X(:,1) - Rmin) < tol); il = sort_by_id(il, id);
strips = zeros(0, 3);
half = floor(numel(il)/2);
for j = [1:half-1, half+1:numel(il)-1]
    strips(end+1,:) = [il(j) il(j+1) d.t_shell]; %#ok<AGROW>
end
strips(end+1,:) = [il(1) il(end) d.t_shell];
for side = [1 -1]
    g = find(X(:,1) >= d.ois_R(1) & X(:,1) <= d.ois_R(2) & side*X(:,2) >= 0);
    g = sort_by_id(g, id);
    for j = 1:numel(g)-1
        strips(end+1,:) = [g(j) g(j+1) d.t_shell]; %#ok<AGROW>
    end
    if d.ois_wrap && numel(g) > 2
        strips(end+1,:) = [g(1) g(end) d.t_shell]; %#ok<AGROW>
    end
end
if isempty(d.gs_node)
    txt = fileread(fullfile(ansys_dir, 'REACT_TF.txt'));
    tok = regexp(txt, '\n\s+(\d+)\s+[-0-9.E+]+\s+[-0-9.E+]+\s+[-0-9.E+]+', 'tokens', 'once');
    d.gs_node = str2double(tok{1});
end
gs = find(id == d.gs_node);
res = tf3d_beam_inplane(X, conn, sec, F, strips, [gs 2], d.n_TF);

% comparison
fid = fopen(fullfile(out_dir, 'verify_tf3d_ansys_beam.txt'), 'w');
say = @(varargin) both(fid, varargin{:});
say('MATLAB in-plane beam vs ANSYS STR_360: %d beams, %d strips per coil (vault %d), GS node %d\n', ...
    size(conn,1), size(strips,1), sum(abs(X(strips(:,1),1) - Rmin) < tol), d.gs_node);
say('applied: sum F_R %.4g N, sum F_z %.4g N; MATLAB GS reaction %.4g N (ANSYS REACT_TF per coil)\n', ...
    sum(F(:,1)), sum(F(:,2)), res.reactions(1,3));
Na = Bm(:,11); Nm = res.elem(:,1);
Ma = Bm(:,13); Mm = res.elem(:,3);
sgn = sign(Mm'*Ma); if sgn == 0, sgn = 1; end
s_el = [0; cumsum(sqrt(sum(diff(X(conn(:,1),:)).^2, 2)))];
say('axial force: ANSYS %.2f..%.2f MN, MATLAB %.2f..%.2f MN, max |diff| %.2f MN (%.1f %% of max)\n', ...
    min(Na)/1e6, max(Na)/1e6, min(Nm)/1e6, max(Nm)/1e6, max(abs(Nm - Na))/1e6, 100*max(abs(Nm - Na))/max(abs(Na)));
say('in-plane moment MZ: ANSYS %.3f..%.3f MN m, MATLAB (sign %+d) %.3f..%.3f MN m, max |diff| %.3f MN m\n', ...
    min(Ma)/1e6, max(Ma)/1e6, sgn, min(sgn*Mm)/1e6, max(sgn*Mm)/1e6, max(abs(sgn*Mm - Ma))/1e6);
uRa = Nd(:,6); uza = Nd(:,7);            % coil 1 in the plane X = 0: global Y radial
say('radial displacement: ANSYS %.2f..%.2f mm, MATLAB %.2f..%.2f mm, max |diff| %.2f mm\n', ...
    1e3*min(uRa), 1e3*max(uRa), 1e3*min(res.U(:,1)), 1e3*max(res.U(:,1)), 1e3*max(abs(res.U(:,1) - uRa)));
say('vertical displacement: ANSYS %.2f..%.2f mm, MATLAB %.2f..%.2f mm, max |diff| %.2f mm\n', ...
    1e3*min(uza), 1e3*max(uza), 1e3*min(res.U(:,2)), 1e3*max(res.U(:,2)), 1e3*max(abs(res.U(:,2) - uza)));
Rs = hypot(Sh(:,2), Sh(:,3)); v = Rs < 2;
iv = abs(X(strips(:,1),1) - Rmin) < tol;
say('vault strips: hoop N ANSYS N11 %.1f..%.1f MN/m, MATLAB %.1f..%.1f MN/m; along-leg N22 ANSYS %.1f..%.1f, MATLAB %.1f..%.1f MN/m\n', ...
    min(Sh(v,5))/1e6, max(Sh(v,5))/1e6, min(res.strip(iv,1))/1e6, max(res.strip(iv,1))/1e6, ...
    min(Sh(v,6))/1e6, max(Sh(v,6))/1e6, min(res.strip(iv,2))/1e6, max(res.strip(iv,2))/1e6);
say('OIS strips: N11 ANSYS %.1f..%.1f, N22 ANSYS %.1f..%.1f MN/m; MATLAB hoop %.1f..%.1f, along-leg %.1f..%.1f MN/m\n', ...
    min(Sh(~v,5))/1e6, max(Sh(~v,5))/1e6, min(Sh(~v,6))/1e6, max(Sh(~v,6))/1e6, ...
    min(res.strip(~iv,1))/1e6, max(res.strip(~iv,1))/1e6, min(res.strip(~iv,2))/1e6, max(res.strip(~iv,2))/1e6);
fclose(fid);
f1 = figure('Name', 'Beam axial force and moment');
subplot(2,1,1); plot(s_el, Na/1e6, '-', s_el, Nm/1e6, '--', 'LineWidth', 1.5); grid on
ylabel('N [MN]'); legend('ANSYS', 'MATLAB'); title('Axial force along coil 1');
subplot(2,1,2); plot(s_el, Ma/1e6, '-', s_el, sgn*Mm/1e6, '--', 'LineWidth', 1.5); grid on
ylabel('M_{in-plane} [MN m]'); xlabel('curvilinear abscissa [m]'); legend('ANSYS', 'MATLAB');
f2 = figure('Name', 'Displacements');
[~, o] = sort(id); sn = [0; cumsum(sqrt(sum(diff(X(o,:)).^2, 2)))];
plot(sn, 1e3*uRa(o), '-', sn, 1e3*res.U(o,1), '--', sn, 1e3*uza(o), '-', sn, 1e3*res.U(o,2), '--', 'LineWidth', 1.5);
grid on; ylabel('[mm]'); xlabel('curvilinear abscissa [m]');
legend('u_R ANSYS', 'u_R MATLAB', 'u_z ANSYS', 'u_z MATLAB');
try
    print(f1, '-dpng', '-r120', fullfile(out_dir, 'beam_N_M.png'));
    print(f2, '-dpng', '-r120', fullfile(out_dir, 'beam_U.png'));
catch
end
R = struct('X', X, 'conn', conn, 'strips', strips, 'res', res, 'N_ansys', Na, 'M_ansys', Ma, ...
    'U_ansys', [uRa uza], 'cfg', d, 'moment_sign', sgn);
end

function g = sort_by_id(g, id)
[~, o] = sort(id(g)); g = g(o);
end

function A = read_num(file)
txt = fileread(file);
txt = regexprep(txt, '(\d)([-+]\d{3})', '$1E$2');     % Fortran exponents without E
lines = regexp(strtrim(txt), '\r?\n', 'split');
first = regexprep(lines{1}, '[,;\t]', ' ');
nc = numel(sscanf(first, '%f'));
A = sscanf(regexprep(txt, '[,;\t]', ' '), '%f');
A = reshape(A(1:nc*floor(numel(A)/nc)), nc, [])';
end

function both(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end
