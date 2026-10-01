function wp = fem_cs_geometry(row, p, g)
%FEM_CS_GEOMETRY Winding-pack description (single grade) from one row of the CS results table.
%
%   wp = FEM_CS_GEOMETRY(row, p, g) rebuilds the geometry of a design point
%   found by SEARCH/SCAN_WP_DESIGNS (a one-row table, or a struct with the
%   same fields) as a flat struct used by every CS FEM module:
%     n_l, n_t       layers (radial) and turns (axial) of one module
%     Cond_w, Cond_h conductor cell width (radial) and height (axial) [m]
%     JT, tins       jacket wall and turn insulation thickness [m]
%     S_Cable, S_JT, S_CICC, S_ins   cable, steel, cable+jacket, insulation areas of the cell [m^2]
%     Ri, Re         inner / outer radius of the turns region [m]
%     H              axial height of the turns region of one module [m]
%     Iop, N_mod     operating current [A], turns of one module
%     n_mod          number of modules of the stack
%     grins_w/h, spacer   ground insulation (radial / axial) and plate between modules [m]
%   Only single-grade designs are supported.
%
%   S_CICC and S_JT are not columns of the results table: they are rebuilt
%   with the same formulas as SIZE_CICC_CABLE (corner radius r_SC + JT).

n_l = rowget(row, 'n_layers'); n_t = rowget(row, 'n_turns');
if numel(n_l) > 1 || numel(n_t) > 1
    error('fem_cs_geometry:multi_grade', 'The FEM module supports single-grade winding packs only.');
end

wp.n_l = n_l; wp.n_t = n_t; wp.N_mod = n_l*n_t;
wp.n_mod = p.n_moduli;
wp.Cond_w = rowget(row, 'Cond_w'); wp.Cond_h = rowget(row, 'Cond_h');
wp.JT = rowget(row, 'JT'); wp.tins = g.tins;
wp.S_Cable = rowget(row, 'S_Cable');
R_J = g.r_SC + wp.JT;
wp.S_CICC = (wp.Cond_w - 2*wp.tins)*(wp.Cond_h - 2*wp.tins) - (4-pi)*R_J^2;
wp.S_JT = wp.S_CICC - wp.S_Cable;
wp.A_cell = wp.Cond_w*wp.Cond_h;
wp.S_ins = wp.A_cell - wp.S_CICC;
wp.SC_w = wp.Cond_w - 2*wp.JT - 2*wp.tins;
wp.SC_h = wp.Cond_h - 2*wp.JT - 2*wp.tins;
wp.Iop = rowget(row, 'Iop_kA')*1e3;
tc = rowget(row, 'type_cable');
if iscell(tc), tc = tc{1}; end
wp.type_cable = tc;

wp.Re = g.Re_WP_outer;                    % outer radius of the turns region
wp.Ri = wp.Re - wp.Cond_w*n_l;            % inner radius of the turns region
wp.H = wp.Cond_h*n_t;
wp.grins_w = g.grins_w; wp.grins_h = g.grins_h; wp.spacer = g.spacer;
wp.WP_h = wp.H + 2*wp.grins_h;            % module height incl. ground insulation
end

function v = rowget(row, name)
v = row.(name);
if iscell(v), v = v{1}; end
end
