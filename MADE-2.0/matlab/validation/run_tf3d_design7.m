function out=run_tf3d_design7(outdir)
%RUN_TF3D_DESIGN7 Reproducible downstream example, independent of a scan.
if nargin<1,outdir=fullfile(pwd,'TF3D_design7');end
here=fileparts(mfilename('fullpath'));addpath(genpath(fullfile(here,'..')));
[row,p]=tf3d_design7_fixture();
out=tf3d_from_design(row,p);
export_tf3d_results(out,outdir,'design7');
export_tf3d_shape(out,outdir,'design7');
plot_tf3d(out,p);
end
