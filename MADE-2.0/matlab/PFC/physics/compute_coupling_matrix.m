function L_matrix = compute_coupling_matrix(geom)
%COMPUTE_COUPLING_MATRIX Single-turn mutual/self inductance matrix of every CS+PFC+plasma conductor.
%
%   L_matrix = COMPUTE_COUPLING_MATRIX(geom) returns the square, symmetric
%   single-turn inductance matrix [H per (A-turn)^2, i.e. Vs/A for one turn
%   on each side] of every row in geom (one row per CS module, PF coil,
%   and the plasma filament), via XLM (shared with the TF/CS pipeline).
%   geom must have columns R, Z, dr, dz (as returned by
%   READ_COIL_GEOMETRY), one row per conductor, in the same row order used
%   throughout PFC (rows 1:p.cs_modules are the CS modules, rows
%   p.cs_modules+1 : p.cs_modules+6 are PF1..PF6, and any further row is
%   the plasma filament).
%
%   Relocated from the PF legacy archive (Poloidal_field.m, part of
%   MADE_PF.7z), which computed this matrix in a standalone script and
%   saved it to L_matrix_monoturn.xlsx for PF_opt_VNS.m to read back in a
%   separate step. Folded directly into the PFC pipeline here (called from
%   MAIN_WP_PFC_DESIGN.m before the per-coil scan) so the coupling matrix
%   no longer needs a manual, external two-script round trip - see manuale
%   PFC.

R = geom.R; Z = geom.Z; dr = geom.dr; dz = geom.dz;
L_matrix = xlm(R, Z, R, Z, dr, dz, 1);
end
