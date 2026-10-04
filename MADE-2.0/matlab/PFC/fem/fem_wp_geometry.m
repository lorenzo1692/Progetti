function wp = fem_wp_geometry(row, p)
%FEM_WP_GEOMETRY Winding-pack description (single grade) from one row of the PFC results table.
%
%   wp = FEM_WP_GEOMETRY(row, p) rebuilds the geometry of a design point
%   found by SEARCH/SCAN_WP_DESIGNS (a one-row table, or a struct with the
%   same fields) as a flat struct used by every FEM module:
%     n_l, n_t       layers (radial) and turns (axial) of the winding pack
%     Cond_w, Cond_h conductor cell width (radial) and height (axial) [m]
%     JT, tins       jacket wall and turn insulation thickness [m]
%     S_Cable, S_JT  cable and steel areas of the cell [m^2]
%     Ri, Re         inner/outer radius of the turns region [m]
%     H              axial height of the turns region [m]
%     R_center       mid-radius of the coil [m]
%     Iop            operating current [A], N_spire total turns
%     type_cable     'LTS' or 'HTS'
%   Only single-grade designs (the PF_opt_VNS.m baseline) are supported.

n_l = rowget(row, 'n_layers'); n_t = rowget(row, 'n_turns');
if numel(n_l) > 1 || numel(n_t) > 1
    error('fem_wp_geometry:multi_grade', 'The FEM module supports single-grade winding packs only.');
end

wp.n_l = n_l; wp.n_t = n_t;
wp.N_spire = n_l*n_t;
wp.Cond_w = rowget(row, 'Cond_w'); wp.Cond_h = rowget(row, 'Cond_h');
wp.JT = rowget(row, 'JT');
wp.tins = rowget(row, 'tins_m');
wp.S_Cable = rowget(row, 'S_Cable_m2'); wp.S_JT = rowget(row, 'S_JT_m2');
wp.R_center = rowget(row, 'R_center');
wp.Iop = rowget(row, 'Iop_kA')*1e3;
tc = rowget(row, 'type_cable');
if iscell(tc), tc = tc{1}; end
wp.type_cable = tc;

wp.Ri = wp.R_center - wp.Cond_w*n_l/2;
wp.Re = wp.R_center + wp.Cond_w*n_l/2;
wp.H = wp.Cond_h*n_t;

wp.SC_w = wp.Cond_w - 2*wp.JT - 2*wp.tins;
wp.SC_h = wp.Cond_h - 2*wp.JT - 2*wp.tins;
wp.A_cell = wp.Cond_w*wp.Cond_h;
wp.S_CICC = wp.S_Cable + wp.S_JT;
wp.S_ins = wp.A_cell - wp.S_CICC;
wp.n_PF = rowget(row, 'n_PF');
wp.FZ_ext = rowget(row, 'FZ_N');
end

function v = rowget(row, name)
v = row.(name);
if iscell(v), v = v{1}; end
end
