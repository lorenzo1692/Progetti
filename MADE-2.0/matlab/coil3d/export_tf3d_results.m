function files = export_tf3d_results(out,outdir,tag)
%EXPORT_TF3D_RESULTS Numeric results, including failed/research diagnostics.
% Independent of CAD validity and plotting; preserves data if graphics fail.
if nargin<2||isempty(outdir),outdir=pwd;end
if nargin<3||isempty(tag),tag='TF3D';end
if isempty(regexp(tag,'^[A-Za-z0-9_.-]+$','once')),error('tf3d:tag','Invalid export tag.');end
if ~exist(outdir,'dir'),mkdir(outdir);end
files.mat=fullfile(outdir,[tag '_results.mat']);save(files.mat,'out','-v7');
f=out.field;
if isfield(f,'full'),full=f.full;else,full=tf3d_full_field(f,out.loop,out.design);end
if isfield(f,'summary'),s=f.summary;else,s=tf3d_field_summary(f,out.design,out.loop);end
files.field=fullfile(outdir,[tag '_field_force_coil1.csv']);
write_csv(files.field, ...
 'turn,layer,s_ref_m,X_m,Y_m,Z_m,BX_T,BY_T,BZ_T,Bcentre_T,FX_bin_N,FY_bin_N,FZ_bin_N,lower_half_reconstructed', ...
 [full.turn full.layer full.s full.points full.B full.Bmag full.force double(full.reconstructed)]);
files.turns=fullfile(outdir,[tag '_turn_summary.csv']);
write_csv(files.turns, ...
 'turn,layer,Bcentre_max_T,Xcentremax_m,Ycentremax_m,Zcentremax_m,Bpeak_ESTIMATE_T,Xpeak_m,Ypeak_m,Zpeak_m,FX_full_N,FY_full_N,FZ_full_N', ...
 [s.turn s.layer s.centre_max s.centre_max_xyz s.peak_estimate s.peak_estimate_xyz s.full_turn_force]);
c=out.checks;files.ripple=fullfile(outdir,[tag '_ripple.csv']);
write_csv(files.ripple,'phi_rad,Bphi_T,Bmag_T',[c.edge_phi c.edge_Bphi c.edge_Bmag]);
if isfield(c,'radial_R')
    files.radial=fullfile(outdir,[tag '_radial_profile.csv']);
    write_csv(files.radial,'R_m,Bphi_coil_plane_T,Bphi_intercoil_T,Bphi_Ampere_mean_T',[c.radial_R c.radial_Bphi c.radial_ampere_mean]);
end
files.status=fullfile(outdir,[tag '_numeric_status.txt']);
fid=fopen(files.status,'w');if fid<0,error('tf3d:io','Cannot write status.');end
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,['%s\nFEM validated: %d\nShape converged: %d\nAmpere pass: %d\nArc tolerance pass: %d\n' ...
 'Use arcs: %d\nOffset pass: %d\nRipple target met: %d\n' ...
 'Centre fields are evaluated. Peak columns are scalar 2D-uplift estimates.\n' ...
 'Lower-half rows are reconstructed by mid-plane symmetry.\n' ...
 'Forces are N per represented bin, not force per unit length.\n' ...
 'Z force cancellation follows enforced symmetry; it is not an independent check.\n' ...
 'This data export does not imply approval of geometry or structural design.\n'], ...
 out.status,out.validated,c.shape_converged,c.ampere_pass,c.arc_fit_pass,out.options.use_arcs,c.geometry_offset_pass,c.ripple_pass);
end
function write_csv(path,header,data)
fid=fopen(path,'w');if fid<0,error('tf3d:io','Cannot write %s',path);end
cleanup=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',header);
fmt=[repmat('%.12g,',1,size(data,2)-1) '%.12g\n'];fprintf(fid,fmt,data');
end
