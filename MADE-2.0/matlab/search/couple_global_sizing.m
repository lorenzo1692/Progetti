function cpl = couple_global_sizing(row, p, opts)
%COUPLE_GLOBAL_SIZING Coupled sizing of a design point: 2D section sizing <-> 3D global model.
%
%   cpl = COUPLE_GLOBAL_SIZING(row, p) closes the loop between the 2D
%   sizing of the inner-leg section (jacket and case sized together:
%   jacket surrogate, case surrogate on the layered model, iterated in the
%   scan) and the 3D beam + shell model of the whole TF system
%   (TF3D_GLOBAL_FROM_DESIGN):
%
%     1. global model of the design: inner-leg vertical force F_inner and
%        the axial load factor k = F_inner / T_bf;
%     2. if k differs from the factor the design was sized with
%        (AXIAL_LOAD_FACTOR(p)) by more than tol: p.axial_load_factor = k
%        and the candidate of the design is sized again
%        (RESIZE_DESIGN_POINT: the full chain, jacket and case nose with
%        the new axial stress);
%     3. global model of the re-sized design (new section, new stiffness,
%        new current centroid), back to 2, until k no longer changes.
%
%   The loop converges in 1-2 re-sizings: k depends on the shape of the
%   coil and on the share of the vertical force between the inner and the
%   outer leg, little on the section.
%
%   opts: tol (0.01 on k), max_iter (4), cal, g, env (passed to
%   RESIZE_DESIGN_POINT), gm (global model already solved on row with p),
%   global (options of TF3D_GLOBAL_FROM_DESIGN), verbose (true).
%
%   cpl: row (coupled design point, same type as the input row), p (input
%   with axial_load_factor = k), k_axial, converged, rejected (the
%   candidate fails a check with the global load: the design is not
%   feasible, choose another), history (struct per iteration: k_sized,
%   k_global, radial_build, Nose, JT, case_Pm, jacket_Pm), gm (global model
%   of the last design).

if nargin < 3, opts = struct(); end
tol = getf(opts, 'tol', 0.01); max_iter = getf(opts, 'max_iter', 4);
gopt = getf(opts, 'global', struct()); gopt.verbose = false;
verbose = getf(opts, 'verbose', true);
ropts = struct();
for f = {'cal', 'g', 'env'}
    if isfield(opts, f{1}), ropts.(f{1}) = opts.(f{1}); end
end
gm = getf(opts, 'gm', []);
if isempty(gm), gm = tf3d_global_from_design(row, p, gopt); end
k = axial_load_factor(p);
hist = struct('k_sized', {}, 'k_global', {}, 'radial_build', {}, 'Nose', {}, 'JT', {}, ...
    'case_Pm', {}, 'jacket_Pm', {});
converged = false; rejected = false;
if verbose
    fprintf('\nCoupled sizing: 2D section (jacket + case) <-> 3D global model\n');
    fprintf('%5s %9s %9s %12s %9s %22s %9s %9s\n', 'iter', 'k sized', 'k global', 'radial [mm]', ...
        'nose[mm]', 'JT range [mm]', 'case Pm', 'jck Pm');
end
for it = 0:max_iter
    r = row; if ~isstruct(r), r = table2struct(r); end
    nl = r.n_layers;
    hist(end+1) = struct('k_sized', k, 'k_global', gm.k_axial, 'radial_build', r.radial_build, ...
        'Nose', r.Nose, 'JT', r.JT(1:nl), 'case_Pm', r.S_T_VT, 'jacket_Pm', getf(r, 'JT_Pm', NaN)); %#ok<AGROW>
    if verbose
        fprintf('%5d %9.4f %9.4f %12.1f %9.1f %10.2f ... %8.2f %9.0f %9.0f\n', it, k, gm.k_axial, ...
            1e3*r.radial_build, 1e3*r.Nose, 1e3*min(r.JT(1:nl)), 1e3*max(r.JT(1:nl)), r.S_T_VT, getf(r, 'JT_Pm', NaN));
    end
    if abs(gm.k_axial - k) <= tol*max(k, 1)
        converged = true; break
    end
    if it == max_iter, break, end
    k = gm.k_axial;
    p.axial_load_factor = k;
    [row2, info] = resize_design_point(row, p, ropts);
    if ~info.ok
        rejected = true;
        if verbose
            fprintf(['  with k = %.3f the candidate fails a check of the sizing (allowables, JT_max, ' ...
                'field): not feasible with the global axial load.\n'], k);
        end
        break
    end
    row = row2;
    gm = tf3d_global_from_design(row, p, gopt);
end
if verbose && converged
    fprintf('  converged: inner-leg axial force = %.3f x T_bf (set axial_load_factor = %.3f in the input to scan with it)\n', ...
        k, k);
elseif verbose && ~rejected
    fprintf('  not converged in %d re-sizings (last change %.3g)\n', max_iter, abs(gm.k_axial - k));
end
cpl = struct('row', row, 'p', p, 'k_axial', k, 'converged', converged, 'rejected', rejected, ...
    'history', hist, 'gm', gm);
end

function v = getf(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
