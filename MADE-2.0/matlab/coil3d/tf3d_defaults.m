function o = tf3d_defaults(p, overrides)
%TF3D_DEFAULTS Numeric settings: Excel fields, then explicit overrides.
% Mode 0 is the reproducible analytic baseline. Mode 1 is experimental.
o = struct('shape_mode',0,'n_seg',180,'shape_tol',1e-4,'max_iter',20, ...
    'shape_relax',0.35,'nose_ol_factor',0.5,'use_arcs',0,'export',1, ...
    'neighbor_gauss',2,'source_block',512,'point_block',128, ...
    'n_eval',48,'arc_tol',0.01,'ampere_tol',1e-3,'energy',0);
names = fieldnames(o);
for i=1:numel(names)
    key=['tf3d_' names{i}];
    if isfield(p,key), o.(names{i})=p.(key); end
end
if nargin>1 && ~isempty(overrides)
    names=fieldnames(overrides);
    for i=1:numel(names)
        if ~isfield(o,names{i}), error('tf3d:option','Unknown option: %s',names{i}); end
        o.(names{i})=overrides.(names{i});
    end
end
names=fieldnames(o);
for i=1:numel(names)
    v=o.(names{i});
    if ~(isnumeric(v)||islogical(v)) || ~isscalar(v) || ~isfinite(v)
        error('tf3d:option','%s must be a finite numeric scalar.',names{i});
    end
end
for key={'shape_mode','use_arcs','export','energy'}
    if ~ismember(o.(key{1}),[0 1]), error('tf3d:option','%s must be 0 or 1.',key{1}); end
end
for key={'n_seg','max_iter','source_block','point_block','n_eval'}
    v=o.(key{1});
    if v<1 || v~=round(v), error('tf3d:option','%s must be a positive integer.',key{1}); end
end
if o.n_seg<32 || o.n_eval<8 || o.shape_relax<=0 || o.shape_relax>1 || ...
        o.shape_tol<=0 || o.arc_tol<=0 || o.ampere_tol<=0 || o.nose_ol_factor<0 || ...
        ~ismember(o.neighbor_gauss,[0 1 2])
    error('tf3d:option','Invalid resolution, tolerance, relaxation or neighbour quadrature.');
end
end
