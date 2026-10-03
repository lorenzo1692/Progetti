function R = verify_tf3d_ansys_global(ansys_dir, out_dir, cfg)
%VERIFY_TF3D_ANSYS_GLOBAL MATLAB 3D beam + shell global model vs ANSYS STR_360.
%
%   R = VERIFY_TF3D_ANSYS_GLOBAL(ansys_dir, out_dir, cfg) rebuilds the
%   ANSYS global model STR_360 with TF3D_GLOBAL_MODEL / TF3D_GLOBAL_SOLVE
%   (all n_TF coils, 6 DOF beams, 4-node shells) on the same centreline
%   nodes and loads, and compares with the export of
%   validation/ansys/POST_TF_BEAM_EXPORT.mac found in ansys_dir:
%     NODE_TF1_U.csv, BEAM_TF1_EL.csv, SHELL_TF1_EL.csv, REACT_TF.txt and
%     the loads FL_TF_1.csv (row j applied to node j, as STR_360).
%
%   cfg: W 0.675, t_wall 0.05, t_shell 0.14, n_TF 12, ois_wrap true (the
%   extra OIS element of STR_360), with_toroidal false (STR_360 of
%   29/09/2026 applies only the radial and vertical components; true adds
%   the toroidal one for an out-of-plane run).
%   Writes out_dir/verify_tf3d_ansys_global.txt and figures.

if nargin < 2 || isempty(out_dir), out_dir = fullfile(pwd, 'verify_tf3d_global'); end
if nargin < 3, cfg = struct(); end
d = struct('W', 0.675, 't_wall', 0.05, 't_shell', 0.14, 'n_TF', 12, 'ois_wrap', true, 'with_toroidal', false);
fn = fieldnames(cfg); for i = 1:numel(fn), d.(fn{i}) = cfg.(fn{i}); end
if ~exist(out_dir, 'dir'), [~] = mkdir(out_dir); end

Nd = read_num(fullfile(ansys_dir, 'NODE_TF1_U.csv'));
Bm = read_num(fullfile(ansys_dir, 'BEAM_TF1_EL.csv'));
Sh = read_num(fullfile(ansys_dir, 'SHELL_TF1_EL.csv'));
FL = read_num(fullfile(ansys_dir, 'FL_TF_1.csv'));
[~, o] = sort(Nd(:,1)); Nd = Nd(o,:);                  % node number = order along the coil
P = Nd(:, 2:4); M = size(P, 1);
F1 = zeros(M, 3); m_ = min(M, size(FL, 1));
F1(1:m_, 2:3) = FL(1:m_, 5:6);
if d.with_toroidal, F1(1:m_, 1) = FL(1:m_, 4); end
mo = struct('n_TF', d.n_TF, 'W', d.W, 't_wall', d.t_wall, 't_shell', d.t_shell, 'ois_wrap', d.ois_wrap);
txt = fileread(fullfile(ansys_dir, 'REACT_TF.txt'));
tok = regexp(txt, '\n\s+(\d+)\s+[-0-9.E+]+', 'tokens', 'once');
gsA = str2double(tok{1});
mdl = tf3d_global_model(P, F1, mo);
if mdl.gs_node ~= find(Nd(:,1) == gsA)
    fprintf('note: gravity-support rule picks node %d, ANSYS %d: using the ANSYS node\n', Nd(mdl.gs_node,1), gsA);
    mo.gs_node = find(Nd(:,1) == gsA);
    mdl = tf3d_global_model(P, F1, mo);
end
t0 = tic; res = tf3d_global_solve(mdl); ts = toc(t0);

fid = fopen(fullfile(out_dir, 'verify_tf3d_ansys_global.txt'), 'w');
say = @(varargin) both(fid, varargin{:});
c1 = 1:M;
say('MATLAB 3D global model vs ANSYS STR_360: %d coils x %d nodes, %d beams, %d shells, %d DOF, solved in %.2f s (residual %.1e)\n', ...
    d.n_TF, M, size(mdl.beam,1), size(mdl.shell,1), res.info.ndof, ts, res.info.residual);
rz = arrayfun(@(q) q.frame_force(3), res.reactions);
say('GS reaction (vertical) per coil: MATLAB %.4g N, ANSYS -1.8378e6 N (REACT_TF); sum of applied Fz %.4g N\n', ...
    mean(rz), sum(mdl.F(:,3)));
% beams of coil 1 (element e = nodes e, e+1)
names = {'FX (axial)', 'MY', 'MZ (in-plane)', 'TQ (torsion)', 'SFZ', 'SFY'};
idxA = [11 12 13 14 15 16];
idxM = [1 5 6 4 3 2];                  % MATLAB local [N Vy Vz T My Mz]
for q = 1:6
    a = Bm(:, idxA(q)); b = res.beam_f(c1, idxM(q));
    s = sign(a'*b); if s == 0, s = 1; end
    say('beam %-14s ANSYS %10.4g .. %10.4g | MATLAB (x%+d) %10.4g .. %10.4g | max |diff| %.3g (%.1f %% of max |ANSYS|)\n', ...
        names{q}, min(a), max(a), s, min(s*b), max(s*b), max(abs(s*b - a)), 100*max(abs(s*b - a))/max(max(abs(a)), eps));
end
for q = [3 6]
    a = Bm(:, idxA(q)); b = res.beam_f(c1, idxM(q)); s = sign(a'*b); if s == 0, s = 1; end
    [~, w] = max(abs(s*b - a));
    say('   largest %s difference at element %d (R %.2f m, z %.2f m): ANSYS %.4g, MATLAB %.4g\n', names{q}, w, ...
        hypot(P(w,1), P(w,2)), P(w,3), a(w), s*b(w));
end
% beam section stresses (BEAM188 SMISC 31-35 at node I: SDIR = N/A,
% SBYT/SBYB = bending about local z (in-plane moment) at the +/- y fibres,
% SBZT/SBZB = bending about local y at the +/- z fibres)
sec = mdl.sec; cy = d.W/2;
sd = res.beam_f(c1,1)/sec.A;
sby = -res.beam_f(c1,6)*cy/sec.Iz;
sbz = res.beam_f(c1,5)*cy/sec.Iy;
cmpS = {'SDIR', 23, sd; 'SBYT', 24, sby; 'SBZT', 26, sbz};
for q = 1:3
    a = Bm(:, cmpS{q,2}); b = cmpS{q,3}; s = sign(a'*b); if s == 0, s = 1; end
    say('beam stress %-5s ANSYS %9.4g .. %9.4g MPa | MATLAB (x%+d) %9.4g .. %9.4g MPa | max |diff| %.3g MPa\n', cmpS{q,1}, ...
        min(a)/1e6, max(a)/1e6, s, min(s*b)/1e6, max(s*b)/1e6, max(abs(s*b - a))/1e6);
end
smax_A = Bm(:,23) + max(abs(Bm(:,24:25)), [], 2) + max(abs(Bm(:,26:27)), [], 2);
smax_M = sd + abs(sby) + abs(sbz);
say('beam max fibre stress |SDIR| + |SBY| + |SBZ|: ANSYS %.1f MPa, MATLAB %.1f MPa\n', max(smax_A)/1e6, max(smax_M)/1e6);
U = res.U(c1, :);
lab = {'UX (toroidal)', 'UY (radial)', 'UZ (vertical)', 'ROTX', 'ROTY', 'ROTZ'};
for q = 1:6
    a = Nd(:, 4+q); b = U(:, q);
    say('node %-14s ANSYS %10.4g .. %10.4g | MATLAB %10.4g .. %10.4g | max |diff| %.3g\n', lab{q}, min(a), max(a), min(b), max(b), max(abs(b - a)));
end
% shells attached to coil 1. The ANSYS export gives, as 'centroid', the
% mid-point of the element edge on coil 1 (the same for the shell on each
% side of the coil); the side is told by the element number: STR_360
% creates the shells of the pair (1, 2) before those of the pair
% (n_TF, 1), region by region (vault, OIS upper, OIS lower).
sc = mdl.shell_coil; nTF = d.n_TF;
e1 = zeros(size(mdl.shell, 1), 2);
e1(sc == 1, :) = mdl.shell(sc == 1, [1 4]);
e1(sc == nTF, :) = mdl.shell(sc == nTF, [2 3]);
keep = sc == 1 | sc == nTF;
Xm = nan(size(mdl.shell, 1), 3);
Xm(keep, :) = (mdl.X(e1(keep,1),:) + mdl.X(e1(keep,2),:))/2;
regA = zeros(size(Sh,1), 1);
RA = hypot(Sh(:,2), Sh(:,3));
regA(abs(RA - min(hypot(P(:,1), P(:,2)))) < 1e-2) = 1;
regA(regA == 0 & Sh(:,4) >= 0) = 2; regA(regA == 0) = 3;
sideA = zeros(size(Sh,1), 1);
for g = 1:3
    k = find(regA == g); [~, o2] = sort(Sh(k,1)); k = k(o2);
    h = floor(numel(k)/2);
    sideA(k(1:h)) = 1; sideA(k(h+1:end)) = nTF;
end
mi = zeros(size(Sh,1), 1); dist = mi;
for k = 1:size(Sh, 1)
    cand = find(sc == sideA(k) & mdl.shell_region == regA(k));
    [dist(k), j] = min(sum((Xm(cand,:) - Sh(k, 2:4)).^2, 2));
    mi(k) = cand(j);
end
say('shells: %d ANSYS elements matched by region, side and coil-1 edge (max distance %.2g m, %d distinct)\n', ...
    numel(mi), sqrt(max(dist)), numel(unique(mi)));
reg = mdl.shell_region(mi);
% the closure element of the vault (coil-1 edge between the last and the
% first node) has, in this model, the edge oriented from node M to node 1;
% STR_360 builds it from node 1 to node M: normal flipped, so bending
% moments M11, M22 and the shear Q13 change sign
flip = ones(numel(mi), 1);
for k = 1:numel(mi)
    nd1 = mod(mdl.shell(mi(k), :) - 1, mdl.M) + 1;
    if any(nd1 == 1) && any(nd1 == mdl.M) && reg(k) == 1, flip(k) = -1; end
end
res.shell_r(mi, [4 5 7]) = res.shell_r(mi, [4 5 7]).*flip;   % M11, M22, Q13 (M12 and Q23 keep their sign: both local y and z flip)
rn = {'vault', 'OIS upper', 'OIS lower'};
sn = {'N11 (hoop)', 'N22 (along coil)', 'N12', 'M11', 'M22', 'M12', 'Q13', 'Q23'};
% STR_360 builds the shells between coil n_TF and coil 1 with the node
% order reversed (I on coil 1): local x and normal are flipped there, so
% the two sides of coil 1 are compared separately
side = mdl.shell_coil(mi);
sd = {sprintf('%d-1', d.n_TF), '1-2'};
for g = 1:3
    for si = 1:2
        k = reg == g & side == (si == 2) + (si == 1)*d.n_TF;
        if ~any(k), continue, end
        for q = [1 2 4 5 7 8]
            a = Sh(k, 4+q); b = res.shell_r(mi(k), q);
            s = sign(a'*b); if s == 0, s = 1; end
            say('shell %-9s %-5s %-16s ANSYS %9.4g .. %9.4g | MATLAB (x%+d) %9.4g .. %9.4g | max |diff| %.3g (%.1f %%)\n', ...
                rn{g}, sd{si}, sn{q}, min(a), max(a), s, min(s*b), max(s*b), max(abs(s*b - a)), ...
                100*max(abs(s*b - a))/max(max(abs(a)), eps));
        end
    end
end
% shell stress intensity at the top and bottom surface (plane stress:
% sigma = N/t +- 6 M/t^2, Tresca including sigma_3 = 0)
t = d.t_shell; r = res.shell_r(mi, :);
SI = zeros(numel(mi), 2);
for sgn = [1 -1]
    sx = r(:,1)/t + sgn*6*r(:,4)/t^2; sy = r(:,2)/t + sgn*6*r(:,5)/t^2; txy = r(:,3)/t + sgn*6*r(:,6)/t^2;
    c = (sx + sy)/2; rr = sqrt(((sx - sy)/2).^2 + txy.^2);
    p = [c + rr, c - rr, zeros(size(c))];
    SI(:, (sgn < 0) + 1) = max(p, [], 2) - min(p, [], 2);
end
for g = 1:3
    k = reg == g;
    say('shell %-9s SINT max (top/bottom): ANSYS %.1f MPa, MATLAB %.1f MPa\n', rn{g}, max(max(Sh(k,13:14)))/1e6, max(max(SI(k,:)))/1e6);
end
fclose(fid);
s_el = [0; cumsum(sqrt(sum(diff(P).^2, 2)))];
f1 = figure('Name', 'Global model: beam forces');
subplot(3,1,1); plot(s_el, Bm(:,11)/1e6, '-', s_el, res.beam_f(c1,1)/1e6, '--', 'LineWidth', 1.4); grid on; ylabel('N [MN]'); legend('ANSYS', 'MATLAB');
title('Coil 1, beam end I');
sM = sign(Bm(:,13)'*res.beam_f(c1,6)); if sM == 0, sM = 1; end
subplot(3,1,2); plot(s_el, Bm(:,13)/1e6, '-', s_el, sM*res.beam_f(c1,6)/1e6, '--', 'LineWidth', 1.4); grid on; ylabel('M in-plane [MN m]');
subplot(3,1,3); plot(s_el, Nd(:,6)*1e3, '-', s_el, U(:,2)*1e3, '--', s_el, Nd(:,7)*1e3, '-', s_el, U(:,3)*1e3, '--', 'LineWidth', 1.4); grid on
ylabel('u [mm]'); xlabel('curvilinear abscissa [m]'); legend('U_R ANSYS', 'U_R MATLAB', 'U_z ANSYS', 'U_z MATLAB');
try, print(f1, '-dpng', '-r120', fullfile(out_dir, 'global_beam.png')); catch, end
R = struct('model', mdl, 'res', res, 'ansys_beam', Bm, 'ansys_node', Nd, 'ansys_shell', Sh, 'shell_match', mi, 'cfg', d);
end

function A = read_num(file)
txt = fileread(file);
txt = regexprep(txt, '(\d)([-+]\d{3})', '$1E$2');
lines = regexp(strtrim(txt), '\r?\n', 'split');
nc = numel(sscanf(regexprep(lines{1}, '[,;\t]', ' '), '%f'));
A = sscanf(regexprep(txt, '[,;\t]', ' '), '%f');
A = reshape(A(1:nc*floor(numel(A)/nc)), nc, [])';
end

function both(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end
