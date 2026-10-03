function m = tf3d_global_model(P, F1, opts)
%TF3D_GLOBAL_MODEL 3D beam + shell model of the whole TF system (all coils).
%
%   m = TF3D_GLOBAL_MODEL(P, F1, opts) builds, for TF3D_GLOBAL_SOLVE, the
%   MATLAB counterpart of the ANSYS global model STR_360:
%
%   - every TF coil is a closed chain of 3D beams along its centreline (the
%     ANSYS BEAM188 with a hollow rectangular HREC section W x W, wall t,
%     standing for the case); coil 1 is given by P (M x 3, in loop order,
%     coil 1 in the plane x = 0 with y the major radius, z vertical, as the
%     MATLAB 3D module and ANSYS CSYS 0); coil c is coil 1 rotated about z
%     by (c-1)*360/n_TF degrees (EGEN in CSYS 1). Beam orientation vector:
%     the toroidal direction of the coil (ANSYS node K = node + 1 m in x for
%     coil 1), so local z is toroidal and the bending about z is in the
%     plane of the coil;
%   - shells (SHELL181, thickness t_shell) between every pair of adjacent
%     coils, one element across the gap, on consecutive centreline nodes of
%     two regions: the straight inner leg (the wedged vault / case nose:
%     nodes with R within vault_tol of the smallest R) and the outer
%     intercoil structure (OIS: nodes with R in ois_R, upper and lower parts
%     separately). Element nodes [coil c, node j; coil c+1, node j; coil
%     c+1, node j+1; coil c, node j+1] (local x along the chord, hoop);
%   - gravity support of every coil: toroidal and vertical displacement
%     fixed at the node nearest to (max R, 0.8 * min z) among the nodes
%     with R <= gs_R (the ANSYS rule NODE(XXMAX, YYMAX, ZZMIN*0.8) on R in
%     [0, 3] m), radial free;
%   - loads: F1 (M x 3, global, coil 1) nodal forces, rotated to every coil
%     (same load on all coils), or F_all (n_TF*M x 3) in opts.F_all for
%     coil-by-coil loads.
%
%   opts (defaults = STR_360 of 29/09/2026): n_TF 12, W 0.675, t_wall
%   0.05, E 205e9, nu 0.3, t_shell 0.14, vault_tol 1e-3, ois_R [3.5 5],
%   ois_wrap false (true reproduces the extra element STR_360 builds
%   between the first and last node of each OIS group, which are not
%   adjacent along the coil), gs_R 3, gs_node [] (index in coil 1 to
%   override the rule), with_vault/with_ois true. Design-driven options
%   (TF3D_GLOBAL_FROM_DESIGN): sec (beam section struct E, nu, A, Iy, Iz,
%   J, Asy, Asz, replacing the HREC W x t_wall), t_vault and t_ois (shell
%   thickness of the vault and of the OIS, default t_shell).
%
%   m: model struct for TF3D_GLOBAL_SOLVE plus M, n_TF, node_of(c, k),
%   shell_region (1 vault, 2 OIS upper, 3 OIS lower), shell_coil, gs_node.

if nargin < 3, opts = struct(); end
d = struct('n_TF', 12, 'W', 0.675, 't_wall', 0.05, 'E', 205e9, 'nu', 0.3, 't_shell', 0.14, ...
    'vault_tol', 1e-3, 'ois_R', [3.5 5], 'ois_wrap', false, 'gs_R', 3, 'gs_node', [], ...
    'with_vault', true, 'with_ois', true, 'F_all', [], 'sec', [], 't_vault', [], 't_ois', []);
fn = fieldnames(opts); for i = 1:numel(fn), d.(fn{i}) = opts.(fn{i}); end
M = size(P, 1); n = d.n_TF;
node = @(c, k) (c-1)*M + k;
X = zeros(n*M, 3); F = zeros(n*M, 6);
beam = zeros(n*M, 2); bk = zeros(n*M, 3);
for c = 1:n
    th = (c-1)*2*pi/n; Rz = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
    X(node(c, 1:M), :) = P*Rz';
    if isempty(d.F_all)
        F(node(c, 1:M), 1:3) = F1*Rz';
    end
    beam(node(c, 1:M), :) = [node(c, 1:M)' node(c, [2:M 1])'];
    bk(node(c, 1:M), :) = repmat([1 0 0]*Rz', M, 1);
end
if ~isempty(d.F_all), F(:, 1:3) = d.F_all; end
W = d.W; t = d.t_wall; Wi = W - 2*t;
sec = struct('E', d.E, 'nu', d.nu, 'A', W^2 - Wi^2, 'Iy', (W^4 - Wi^4)/12, 'Iz', (W^4 - Wi^4)/12, ...
    'J', (W - t)^3*t, 'Asy', 2*t*(W - t), 'Asz', 2*t*(W - t));
if ~isempty(d.sec)
    for f = {'E', 'nu', 'A', 'Iy', 'Iz', 'J', 'Asy', 'Asz'}
        sec.(f{1}) = d.sec.(f{1});
    end
end
if isempty(d.t_vault), d.t_vault = d.t_shell; end
if isempty(d.t_ois), d.t_ois = d.t_shell; end
% shell regions on coil 1 (indices along the loop)
R1 = hypot(P(:,1), P(:,2)); z1 = P(:,3);
pairs = zeros(0, 3);                                  % [j j+1 region]
nxt = [2:M 1];
if d.with_vault
    inV = abs(R1 - min(R1)) < d.vault_tol;
    j = find(inV & inV(nxt));
    pairs = [pairs; j nxt(j)' ones(numel(j), 1)];
end
if d.with_ois
    for reg = [2 3]
        if reg == 2, inO = R1 >= d.ois_R(1) & R1 <= d.ois_R(2) & z1 >= 0;
        else,        inO = R1 >= d.ois_R(1) & R1 <= d.ois_R(2) & z1 <= 0; end
        j = find(inO & inO(nxt));
        pairs = [pairs; j nxt(j)' reg*ones(numel(j), 1)]; %#ok<AGROW>
        if d.ois_wrap
            g = find(inO);
            pairs = [pairs; g(1) g(end) reg]; %#ok<AGROW>
        end
    end
end
np = size(pairs, 1);
shell = zeros(n*np, 4); reg = zeros(n*np, 1); sc = zeros(n*np, 1);
for c = 1:n
    c2 = mod(c, n) + 1;
    r = (c-1)*np + (1:np);
    shell(r, :) = [node(c, pairs(:,1)) node(c2, pairs(:,1)) node(c2, pairs(:,2)) node(c, pairs(:,2))];
    reg(r) = pairs(:,3); sc(r) = c;
end
% gravity supports
if isempty(d.gs_node)
    cand = find(R1 <= d.gs_R);
    tgt = [max(R1(cand)) 0.8*min(z1(cand))];
    [~, i] = min((R1(cand) - tgt(1)).^2 + (z1(cand) - tgt(2)).^2);
    gs = cand(i);
else
    gs = d.gs_node;
end
sup = struct('node', {}, 'frame', {}, 'mask', {});
for c = 1:n
    k = node(c, gs); x = X(k,:);
    er = [x(1:2) 0]/norm(x(1:2)); ez = [0 0 1]; et = cross(ez, er);
    sup(c) = struct('node', k, 'frame', [er; et; ez], 'mask', [0 1 1 0 0 0]);
end
t_reg = [d.t_vault d.t_ois d.t_ois];
m = struct('X', X, 'beam', beam, 'beam_k', bk, 'sec', sec, 'shell', shell, ...
    'shell_t', reshape(t_reg(max(reg, 1)), [], 1), 'shell_E', d.E, 'shell_nu', d.nu, 'sup', sup, 'F', F, ...
    'M', M, 'n_TF', n, 'shell_region', reg, 'shell_coil', sc, 'gs_node', gs, 'opts', d);
m.node_of = node;
end
