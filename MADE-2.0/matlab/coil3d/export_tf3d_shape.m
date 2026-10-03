function files = export_tf3d_shape(out,outdir,tag)
%EXPORT_TF3D_SHAPE CAD polylines [mm], arc definitions and turn centrelines.
% Radial case envelopes are geometric assumptions, not full 3D CAD solids.
if nargin<2||isempty(outdir),outdir=pwd;end
if nargin<3||isempty(tag),tag='TF3D';end
if isempty(regexp(tag,'^[A-Za-z0-9_.-]+$','once')),error('tf3d:tag','Use letters, numbers, dot, hyphen or underscore.');end
c=out.checks;
if ~c.shape_converged||~c.geometry_offset_pass||~c.ampere_pass|| ...
        (out.options.use_arcs&&~c.arc_fit_pass)
    error('tf3d:export','Numerical/geometric checks failed: CAD export blocked. Inspect out.checks.');
end
if ~exist(outdir,'dir'),mkdir(outdir);end
d=out.design;L=out.loop;
u_outer=d.wp_u_max+d.ground_insulation+d.nose_il* ...
    (1+(out.options.nose_ol_factor-1)*(L.rz(:,1)-d.r1)/(d.r2-d.r1));
ci=L.rz+d.case_u_min*L.normal;
co=L.rz+u_outer.*L.normal;
wi=L.rz+d.wp_u_min*L.normal;wo=L.rz+d.wp_u_max*L.normal;
cl=(ci+co)/2;
files.shape=fullfile(outdir,[tag '_shape_mm.csv']);
writecsv(files.shape,'r_case_CL_mm,z_case_CL_mm,r_case_plasma_mm,z_case_plasma_mm,r_case_back_mm,z_case_back_mm,r_WP_plasma_mm,z_WP_plasma_mm,r_WP_back_mm,z_WP_back_mm,r_current_mm,z_current_mm', ...
    1e3*[cl ci co wi wo L.rz]);
files.arcs=fullfile(outdir,[tag '_arcs_mm.csv']);
a=out.arcs;
writecsv(files.arcs,'arc,centre_r_mm,centre_z_mm,radius_mm,sweep_start_rad,sweep_end_rad', ...
    [(1:3)' 1e3*a.centres 1e3*a.radii(:) a.angles(1:3)' a.angles(2:4)']);
files.turns=fullfile(outdir,[tag '_turn_centrelines_mm.csv']);
turns=zeros((L.n_seg+1)*d.n_turns,5);
for j=1:d.n_turns
    rz=L.rz+d.u(j)*L.normal;ix=(j-1)*(L.n_seg+1)+(1:L.n_seg+1);
    turns(ix,:)=[repmat(j,L.n_seg+1,1) (0:L.n_seg)' -1e3*d.v(j)*ones(L.n_seg+1,1) 1e3*rz];
end
writecsv(files.turns,'turn,node,X_mm,Y_mm,Z_mm',turns);
files.mat=fullfile(outdir,[tag '_results.mat']);save(files.mat,'out','-v7');
files.notes=fullfile(outdir,[tag '_README.txt']);
fid=fopen(files.notes,'w');if fid<0,error('tf3d:io','Cannot open export notes.');end
clean=onCleanup(@()fclose(fid));
fprintf(fid,['TF3D research export. Coordinates in mm; angles in radians.\n' ...
    'Coil 1 in YZ plane; local toroidal coordinate v=-X. Current goes up the inner leg.\n' ...
    'Turn paths are closed independent circuits; winding crossovers/leads omitted.\n' ...
    'Case back thickness interpolates in radius between inner and outer nose.\n' ...
    'These are section envelopes, not a wedged 3D case solid or clearance certification.\n' ...
    'Three-arc file describes the fitted reference even if use_arcs=0.\n' ...
    'Shape mode=%d; use_arcs=%d; convergence=%d; arc max error=%.6g mm.\n' ...
    '3D peak is a scalar 2D-uplift estimate. Forces require resolution studies.\n' ...
    'Fz_half/2 is a mean of two axial cut forces, not automatically T_inner.\n' ...
    'No validated 3D FEM comparison. No feedback to scan/discharge sizing.\n'], ...
    out.options.shape_mode,out.options.use_arcs,c.shape_converged,1e3*a.max_error);
end
function writecsv(path,header,data)
fid=fopen(path,'w');if fid<0,error('tf3d:io','Cannot open %s',path);end
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',header);
fmt=[repmat('%.12g,',1,size(data,2)-1) '%.12g\n'];fprintf(fid,fmt,data');
end
