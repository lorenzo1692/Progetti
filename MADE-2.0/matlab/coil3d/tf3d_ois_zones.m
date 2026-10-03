function Z = tf3d_ois_zones(P, spec, opts)
%TF3D_OIS_ZONES Define and check the outer intercoil structures before a global run.
%
%   Z = TF3D_OIS_ZONES(P, spec, opts) turns a user specification of the
%   outer intercoil structures (OIS) into node ranges of the coil-1
%   centreline P (M x 3, plane x = 0, y = R, z vertical, loop starting at
%   the inner mid-plane, as RUN_TF3D_GLOBAL builds it), checks them, and
%   optionally shows them and asks for confirmation.
%
%   spec: K x 3 [s_start s_end t] - each row one OIS panel between adjacent
%   coils: position as fractions of the centreline length (0 = inner
%   mid-plane, ~0.5 = outer mid-plane, increasing upwards on the inner leg),
%   t = panel thickness [m]. Use TF3D_OIS_ZONES(P, [], opts) to only plot the
%   centreline with the s ticks, to choose the zones.
%
%   opts: n_TF (required), W (case width [m], for the toroidal gap check),
%   vault_idx (nodes of the wedged vault, excluded), pf (K x 5 [rc zc dr dz I]
%   rings, clash check), clearance (min distance from PF coils [m], default
%   0.1), plot (default true if a display is available), confirm (default
%   true if plot: ask y/n in the command window).
%
%   Checks (errors stop the run, warnings are printed):
%     error   s outside [0 1], s_start >= s_end, t <= 0, zone with < 2 nodes,
%             zones overlapping each other or the vault;
%     warning panel width between adjacent cases (2 R sin(pi/n_TF) - W) <= 0
%             somewhere in the zone (cases touch: wedge, not an OIS),
%             PF coil closer than clearance to the zone,
%             thickness larger than the free width between the cases.
%
%   Z: struct array idx (nodes of coil 1, in loop order), t, s, R range,
%   z range, length [m], gap [min max] free width between cases [m].

if nargin < 3, opts = struct(); end
M = size(P, 1);
s = [0; cumsum(sqrt(sum(diff(P).^2, 2)))]; Ltot = s(end) + norm(P(1,:) - P(end,:)); s = s/Ltot;
R = hypot(P(:,1), P(:,2)); z = P(:,3);
n_TF = opts.n_TF;
W = getf(opts, 'W', 0);
vault = false(M, 1); vi = getf(opts, 'vault_idx', []); vault(vi) = true;
pf = getf(opts, 'pf', zeros(0, 5)); clr = getf(opts, 'clearance', 0.1);
show = getf(opts, 'plot', has_display());
Z = struct('idx', {}, 't', {}, 's', {}, 'R', {}, 'z', {}, 'length', {}, 'gap', {});
used = false(M, 1);
for k = 1:size(spec, 1)
    s1 = spec(k,1); s2 = spec(k,2); t = spec(k,3);
    if s1 < 0 || s2 > 1 || s1 >= s2, error('tf3d_ois_zones:s', 'OIS %d: need 0 <= s_start < s_end <= 1 (got %g %g).', k, s1, s2); end
    if ~(t > 0), error('tf3d_ois_zones:t', 'OIS %d: thickness must be > 0.', k); end
    idx = find(s >= s1 & s <= s2);
    if numel(idx) < 2, error('tf3d_ois_zones:short', 'OIS %d: fewer than 2 centreline nodes in [%g %g]; refine n_seg or widen the zone.', k, s1, s2); end
    if any(vault(idx)), error('tf3d_ois_zones:vault', 'OIS %d overlaps the wedged inner-leg vault.', k); end
    if any(used(idx)), error('tf3d_ois_zones:overlap', 'OIS %d overlaps another OIS zone.', k); end
    used(idx) = true;
    gap = 2*R(idx)*sin(pi/n_TF) - W;
    Z(k).idx = idx(:)'; Z(k).t = t; Z(k).s = [s1 s2];
    Z(k).R = [min(R(idx)) max(R(idx))]; Z(k).z = [min(z(idx)) max(z(idx))];
    Z(k).length = (s(idx(end)) - s(idx(1)))*Ltot; Z(k).gap = [min(gap) max(gap)];
    if min(gap) <= 0
        warning('tf3d_ois_zones:no_gap', 'OIS %d: adjacent cases touch (free width %.3f m): this is a wedged region, not an OIS.', k, min(gap));
    elseif t > max(gap)
        warning('tf3d_ois_zones:thick', 'OIS %d: thickness %.3f m larger than the free width between cases (%.3f m).', k, t, max(gap));
    end
    for q = 1:size(pf, 1)
        dR = max(abs(R(idx) - pf(q,1)) - pf(q,3)/2, 0); dz = max(abs(z(idx) - pf(q,2)) - pf(q,4)/2, 0);
        dmin = min(hypot(dR, dz));
        if dmin < clr
            warning('tf3d_ois_zones:pf', 'OIS %d is %.3f m from PF/CS ring %d (clearance %.3f m).', k, dmin, q, clr);
        end
    end
end

fprintf('\nOIS zones (%d):\n', numel(Z));
for k = 1:numel(Z)
    fprintf('  %d: s %.3f-%.3f, R %.2f-%.2f m, z %+.2f..%+.2f m, length %.2f m, t %.3f m, free width %.3f-%.3f m\n', ...
        k, Z(k).s, Z(k).R, Z(k).z, Z(k).length, Z(k).t, Z(k).gap);
end

if show
    figure('Name', 'TF global model: OIS zones'); hold on; axis equal; grid on;
    plot(R, z, 'k-');
    tick = 0:0.05:1;
    for q = tick
        [~, j] = min(abs(s - q)); plot(R(j), z(j), 'k.');
        text(R(j), z(j), sprintf(' %.2f', q), 'FontSize', 7);
    end
    plot(R(vault), z(vault), 'b-', 'LineWidth', 4);
    col = lines(max(1, numel(Z)));
    for k = 1:numel(Z)
        plot(R(Z(k).idx), z(Z(k).idx), '-', 'Color', col(k,:), 'LineWidth', 3 + 40*Z(k).t);
        j = Z(k).idx(round(end/2));
        text(R(j), z(j), sprintf('  OIS %d, t = %.0f mm', k, 1e3*Z(k).t), 'Color', col(k,:));
    end
    for q = 1:size(pf, 1)
        rectangle('Position', [pf(q,1)-pf(q,3)/2, pf(q,2)-pf(q,4)/2, pf(q,3), pf(q,4)], 'EdgeColor', [0.6 0.6 0.6]);
    end
    xlabel('R [m]'); ylabel('z [m]'); title('Centreline s/s_{max} (black), vault (blue), OIS (colour)');
    drawnow;
    if getf(opts, 'confirm', true) && ~isempty(spec)
        a = input('Run the global model with these OIS zones? [y/n] ', 's');
        if ~strcmpi(strtrim(a), 'y'), error('tf3d_ois_zones:aborted', 'OIS zones not confirmed.'); end
    end
end
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function tf = has_display()
try
    tf = usejava('desktop') || (exist('OCTAVE_VERSION', 'builtin') && ~isempty(getenv('DISPLAY')));
catch
    tf = false;
end
end
