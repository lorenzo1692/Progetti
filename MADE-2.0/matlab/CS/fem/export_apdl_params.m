function files = export_apdl_params(row, p, g, out_dir, tag)
%EXPORT_APDL_PARAMS Write the chosen CS design point as input for the ANSYS CS model.
%
%   files = EXPORT_APDL_PARAMS(row, p, g, out_dir, tag) writes, for the
%   ANSYS model CS_model_2 (STR/MAG, units mm):
%
%   <tag>_Parametri_CS_design.lgw   overrides of the CS block of
%       Parametri_TCM_*.lgw (Re, WP_h, GR_INS_w/h, RENF_h, N_MOD, N_GRADES,
%       N_LAYERS_CS, N_TURNS_CS, EQV_Cable_Area, JT_, TINS_, R_SC_). Include
%       it right after the original Parametri file, e.g. by adding
%       /INPUT,'<tag>_Parametri_CS_design','lgw',../   below the line that
%       reads Parametri_TCM_*.lgw in Geom_CS.dat and STR_CS.dat (MAG too).
%   <tag>_DESIGN.csv   scenario table in the format of STR/input/SN.csv
%       (row 1 time, row 2 plasma, rows 3-8 CS modules, rows 9-14 PF, in
%       Ampere-turns): two time points with every CS module at the design
%       current Iop, no plasma and no PF - the load case assumed by the
%       analytic design models and by case 1 of FEM_CS_VERIFY. Add 'DESIGN'
%       to SCENARI (and NO_STP = 2, TIMEPOINT = 2) in the Parametri file.
%
%   The arrays (N_LAYERS_CS, N_TURNS_CS, ...) are only assigned, not
%   re-dimensioned: they already exist from the original Parametri file.
%   Not yet run in ANSYS.
%
%   Differences to keep in mind when comparing with the ANSYS results:
%   the original Parametri file sets R_SC_(1) = JT_(1) (cable corner
%   radius equal to the jacket thickness); here R_SC_ is the MATLAB r_SC,
%   so the cable-plus-jacket cross-section matches the design point.
%   RENF_h is the spacer plate between modules. N_MOD must stay 6 in the
%   ANSYS model (module names CSL3 ... CSU3 are fixed in Parametri).

if nargin < 4 || isempty(out_dir), out_dir = pwd; end
if nargin < 5, tag = 'CS'; end
if ~isfolder(out_dir), mkdir(out_dir); end
wp = fem_cs_geometry(row, p, g);
if wp.n_mod ~= 6
    warning('export_apdl_params:n_mod', 'The ANSYS model has 6 CS modules; the design has %d.', wp.n_mod);
end

f1 = fullfile(out_dir, [tag '_Parametri_CS_design.lgw']);
fid = fopen(f1, 'w');
fprintf(fid, '!!! CS design point exported by export_apdl_params.m - overrides the CS block of Parametri_TCM_*.lgw\n');
fprintf(fid, '!!! Iop = %.1f A, %d layers x %d turns per module, Fz (stack compression) = %.2f MN\n', ...
    wp.Iop, wp.n_l, wp.n_t, rowget(row, 'Fz_MN'));
fprintf(fid, 'N_GRADES = 1\n');
fprintf(fid, 'N_MOD = %d\n', wp.n_mod);
fprintf(fid, 'GR_INS_w = (%.6g)*Increm\n', wp.grins_w);
fprintf(fid, 'GR_INS_h = (%.6g)*Increm\n', wp.grins_h);
fprintf(fid, 'RENF_h = (%.6g)*Increm\n', wp.spacer);
fprintf(fid, 'Re = (%.6g)*Increm\n', g.Re);
fprintf(fid, 'WP_h = (%.6g)*Increm\n', wp.WP_h);
fprintf(fid, 'N_LAYERS_CS(1) = %d\n', wp.n_l);
fprintf(fid, 'N_TURNS_CS(1) = %d\n', wp.n_t);
fprintf(fid, '*VOPER,NSPR,N_LAYERS_CS,MULT,N_TURNS_CS\n');
fprintf(fid, '*VSCFUN,N_spire, SUM, NSPR\n');
fprintf(fid, 'EQV_Cable_Area(1) = (%.6g)*SCALEV**2\n', wp.S_Cable);
fprintf(fid, 'Cond_h_(1) = (WP_h-2*GR_INS_h)/N_TURNS_CS(1)\n');
fprintf(fid, 'JT_(1) = (%.6g)*Increm\n', wp.JT);
fprintf(fid, 'R_SC_(1) = (%.6g)*Increm\n', g.r_SC);
fprintf(fid, 'TINS_(1) = (%.6g)*Increm\n', wp.tins);
fprintf(fid, 'SC_h_(1) = (Cond_h_(1)-JT_(1)*2-TINS_(1)*2)\n');
fprintf(fid, 'SC_w_(1) = (EQV_Cable_Area(1)+(4-pi)*R_SC_(1)**2)/SC_h_(1)\n');
fprintf(fid, 'R_J_(1) = (R_SC_(1) + JT_(1))\n');
fprintf(fid, 'Cond_w_(1) = SC_w_(1)+JT_(1)*2+TINS_(1)*2\n');
fprintf(fid, 'Re_(1) = (Re-GR_INS_w)\n');
fprintf(fid, 'Ri_(1) = (Re_(1)-Cond_w_(1)*N_LAYERS_CS(1))\n');
fprintf(fid, 'Ri = (Ri_(N_GRADES)-GR_INS_w)\n');
fprintf(fid, 'WP_w = (Re-Ri)\n');
fprintf(fid, 'RENF_w = (Re-Ri)\n');
fprintf(fid, 'ELEM_SIZE(1) = EQV_Cable_Area(1)*2/SCALEV\n');
fclose(fid);

f2 = fullfile(out_dir, [tag '_DESIGN.csv']);
NI = wp.Iop*wp.N_mod;
M = zeros(14, 2);
M(1, :) = [0 1];
M(3:8, :) = NI;
fid = fopen(f2, 'w');
for k = 1:14
    fprintf(fid, '%.6E,%.6E\n', M(k, 1), M(k, 2));
end
fclose(fid);
files = {f1, f2};
fprintf('ANSYS input written: %s\n                     %s\n', f1, f2);
end

function v = rowget(row, name)
v = row.(name);
if iscell(v), v = v{1}; end
end
