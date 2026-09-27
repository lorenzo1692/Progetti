function tg = wp_turn_geometry(row, p)
%WP_TURN_GEOMETRY Per-layer turn cross-section (cable, jacket, insulation).
%
%   tg = WP_TURN_GEOMETRY(row, p) is the single place where the geometry
%   of one turn is derived from a design point, for both conductor shapes
%   of SIZE_CICC_CABLE:
%
%     shape_cable = 201, 'Rect' - rectangular cable SC_w x SC_h with corner
%       fillet r_SC = clamp(JT, r_SC_min, r_SC_max); jacket of constant
%       thickness JT around it (outer corner radius r_SC + JT).
%     shape_cable = 200, 'RIS' (round-in-square) - circular cable of
%       diameter d = sqrt(4*A_cable/pi) centred in a square jacket of side
%       a = Cond_w - 2*t_ins; outer jacket corner radius
%       R_J = clamp(JT, r_SC_min, r_SC_max) (as SIZE_CICC_CABLE, whose
%       S_CICC uses the same R_J); JT = (a - d)/2 is the MINIMUM jacket
%       thickness (at the middle of the four sides), the jacket is thicker
%       towards the corners.
%
%   In both cases the turn insulation (t_ins) follows the jacket outer
%   contour (outer corner radius R_J + t_ins) and the corner filler fills
%   the rest of the square cell Cond_w x Cond_h.
%
%   The shape is taken from the SOLUTION (row.shape_cable, saved by
%   SCAN_WP_DESIGNS since this change) and only if the row does not carry
%   it (hand-written structs, results files written before this change)
%   from p.shape_cable. A shape stored in the row is never overridden by
%   the current p, so a RIS solution stays RIS when reloaded with a Rect
%   input file. Unknown codes raise an error (no silent fallback).
%
%   tg fields (1 x n_layers unless noted):
%     shape_code (scalar), shape_name (char), is_round (logical scalar)
%     cell_w, cell_h   - insulated cell
%     tins (scalar)    - turn insulation per face
%     jk_w, jk_h, jk_r - jacket outer rounded rectangle and corner radius
%     cab_w, cab_h, cab_r - cable rounded rectangle (RIS: w = h = d, r = d/2)
%     cab_d            - cable diameter (RIS), NaN for Rect
%     jt_min           - minimum jacket thickness
%     A_cable, A_jacket - cable and jacket steel areas
%     A_cable_row      - cable area implied by the row (for the check below)

code = resolve_shape_cable(row, p);
nl = row.n_layers;
Cond_w = row.Cond_w(1:nl); Cond_h = row.Cond_h(1:nl); JT = row.JT(1:nl);
tins = p.turn_insulation_nominal*p.Increm;
rmin = get_or(p, 'r_SC_min', 2e-3); rmax = get_or(p, 'r_SC_max', 6e-3);

tg.shape_code = code;
tg.cell_w = Cond_w; tg.cell_h = Cond_h; tg.tins = tins; tg.jt_min = JT;
tg.jk_w = Cond_w - 2*tins; tg.jk_h = Cond_h - 2*tins;
switch code
    case 201
        tg.shape_name = 'Rect'; tg.is_round = false;
        r = min(max(JT, rmin), rmax);
        tg.cab_w = Cond_w - 2*JT - 2*tins;
        tg.cab_h = Cond_h - 2*JT - 2*tins;
        tg.cab_r = min(r, 0.49*min(tg.cab_w, tg.cab_h));
        tg.cab_d = nan(1, nl);
        tg.jk_r = tg.cab_r + JT;
        tg.A_cable = tg.cab_w.*tg.cab_h - (4-pi)*tg.cab_r.^2;
    case 200
        tg.shape_name = 'RIS'; tg.is_round = true;
        if any(abs(Cond_w - Cond_h) > 1e-9*max(Cond_w))
            error('wp_turn_geometry:ris_not_square', ...
                ['RIS (shape_cable = 200) needs a square cell (Cond_w = Cond_h, as SIZE_CICC_CABLE ' ...
                 'builds it); this row has Cond_w ~= Cond_h: the design point is not a RIS one.']);
        end
        d = tg.jk_w - 2*JT;                         % diameter implied by the saved JT
        tg.cab_d = d; tg.cab_w = d; tg.cab_h = d; tg.cab_r = d/2;
        tg.jk_r = min(max(JT, rmin), rmax);
        tg.A_cable = pi*d.^2/4;
    otherwise
        error('wp_turn_geometry:shape', 'Unsupported shape_cable = %g (200 = RIS, 201 = Rect).', code);
end
if any(tg.cab_w <= 0) || any(tg.cab_h <= 0)
    error('wp_turn_geometry:size', 'Cable size <= 0: check Cond_w/Cond_h against JT and the turn insulation.');
end
if any(tg.jk_r >= 0.5*min(tg.jk_w, tg.jk_h))
    error('wp_turn_geometry:corner', 'Jacket outer corner radius >= half the jacket side.');
end
tg.A_jacket = tg.jk_w.*tg.jk_h - (4-pi)*tg.jk_r.^2 - tg.A_cable;
% consistency with the cable area the row was sized for (S_Cable), if saved
tg.A_cable_row = nan(1, nl);
if has_field(row, 'S_Cable')
    s = row.S_Cable; tg.A_cable_row = s(1:nl);
end
end

function code = resolve_shape_cable(row, p)
if has_field(row, 'shape_cable')
    code = row.shape_cable;
    if iscell(code), code = code{1}; end
    code = code(1);
elseif isfield(p, 'shape_cable') && ~isempty(p.shape_cable)
    code = p.shape_cable;
else
    error('wp_turn_geometry:no_shape', ...
        'The design point does not say whether it is RIS or Rect (no row.shape_cable nor p.shape_cable).');
end
if ~(isnumeric(code) && isscalar(code) && any(code == [200 201]))
    error('wp_turn_geometry:shape', 'Unsupported shape_cable value (200 = RIS, 201 = Rect).');
end
end

function tf = has_field(row, name)
if isstruct(row)
    tf = isfield(row, name);
else                                   % MATLAB table row
    tf = any(strcmp(row.Properties.VariableNames, name));
end
end

function v = get_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name)), v = s.(name); else, v = default; end
end
