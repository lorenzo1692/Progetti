function report = validate_tf3d(run_design)
%VALIDATE_TF3D Native MATLAB/Octave test gate; no Excel dependency.
% validate_tf3d(false): segment kernel, shape, offsets, arcs, current scaling.
% validate_tf3d(true): also run Design 7 at two resolutions and export check.
% These tests do not constitute 3D FEM validation or mechanical approval.
if nargin<1,run_design=false;end
here=fileparts(mfilename('fullpath'));addpath(genpath(fullfile(here,'..')));
report=struct();
A=[0 0 -1];C=[0 0 1];P=[1 0 0];expected=[0 2e-7/sqrt(2) 0];
B=biot_savart_segments(P,A,C,1);
assert(norm(B-expected)<1e-19,'Finite-segment analytic field failed.');
assert(norm(biot_savart_segments(P,C,A,1)+B)<1e-19,'Reversal failed.');
assert(norm(biot_savart_segments(P,A,C,3)-3*B)<1e-19,'Current scaling failed.');
assert(norm(biot_savart_segments([0 0 3],A,C,1))==0,'Collinear exterior field failed.');
assert(norm(biot_savart_segments(P,A+[2 4 5],C+[2 4 5],1)- ...
    biot_savart_segments(P-[2 4 5],A,C,1))<1e-19,'Translation failed.');
failed=false;
try,biot_savart_segments([0 0 0],A,C,1);catch err,failed=strcmp(err.identifier,'tf3d:singular');end
assert(failed,'On-source singularity must be rejected.');
errors=zeros(1,3);
for j=1:3
    n=32*2^j;t=linspace(0,2*pi,n+1)';Q=[cos(t) sin(t) zeros(size(t))];
    z=[0;.5;2];pts=[zeros(3,2) z];
    bb=biot_savart_segments(pts,Q(1:end-1,:),Q(2:end,:),1);
    exact=2*pi*1e-7./(1+z.^2).^(3/2);errors(j)=max(abs(bb(:,3)./exact-1));
end
assert(errors(3)<errors(2)/3 && errors(2)<errors(1)/3,'Circle refinement failed.');
report.circle_relative_errors=errors;
[row,p]=tf3d_design7_fixture();d=tf3d_design_geometry(row,p);
assert(d.n_turns==104 && abs(d.NI-6721312)<1e-6,'Turn/current mapping failed.');
assert(abs(d.r1-mean(d.y))<1e-14,'Current centroid failed.');
assert(abs((d.r2-d.d_c)-compute_outer_radius(p))<1e-12,'Mirrored outer radius failed.');
sh=tf_bending_free_shape(d.r1,d.r2,801);ar=tf_three_arc_fit(sh.upper);
report.arc_fit_10mm_gate_pass=ar.max_error<0.01;
% Keep the fit tolerance fixed and report failure; baseline uses the free span.
fprintf('Arc fit 10 mm gate: %d (error %.2f mm)\n',report.arc_fit_10mm_gate_pass,1e3*ar.max_error);
assert(norm(ar.joins(1,:)-sh.upper(end,:))<1e-10 && ...
    norm(ar.joins(end,:)-sh.upper(1,:))<1e-10,'Arc endpoints failed.');
for j=1:2
    n1=(ar.joins(j+1,:)-ar.centres(j,:))/ar.radii(j);
    n2=(ar.joins(j+1,:)-ar.centres(j+1,:))/ar.radii(j+1);
    assert(norm(n1-n2)<1e-10,'Arc tangency failed.');
end
L=tf3d_resample_loop(sh.upper,90);assert(norm(L.rz(1,:)-L.rz(end,:))<1e-12,'Loop closure failed.');
o=tf3d_defaults(p,struct('n_seg',90));src=tf3d_coil_filaments(L,d,p,o);
assert(norm(sum((src.C-src.A).*src.I,1))<1e-5,'Closed-circuit current balance failed.');
checks=tf3d_toroidal_checks(src,d,p,o);assert(checks.ampere_pass,'Ampere mean failed.');
report.ampere_relative_error=checks.ampere_rel_error;report.arc_max_error=ar.max_error;
% Direct kernel check of the mid-plane reconstruction used for full-coil plots.
probe=[-.06 1.15 .7;.08 2.5 2.0;.04 4.8 .9];
bup=biot_savart_segments(probe,src.A,src.C,src.I,o);
probe(:,3)=-probe(:,3);blo=biot_savart_segments(probe,src.A,src.C,src.I,o);
expected=bup;expected(:,3)=-expected(:,3);
report.reflection_error=max(abs(blo(:)-expected(:)))/max(abs(bup(:)));
assert(report.reflection_error<1e-10,'Mid-plane field parity failed.');

if run_design
    lo=tf3d_from_design(row,p,struct('n_seg',90,'n_eval',45,'export',0));
    hi=tf3d_from_design(row,p,struct('n_seg',180,'n_eval',90,'export',0));
    assert(lo.checks.ampere_pass&&hi.checks.ampere_pass,'Design 7 Ampere failed.');
    assert(hi.checks.geometry_offset_pass,'Offset orientation failed.');
    report.centre_refinement=max(abs(lo.checks.midplane_3D./hi.checks.midplane_3D-1));
    report.force_refinement=abs(lo.field.mean_axial_cut_force/hi.field.mean_axial_cut_force-1);
    % Development gates fixed before native execution; not adjusted to fit results.
    assert(report.centre_refinement<0.02,'Centre field not converged to 2%%.');
    assert(report.force_refinement<0.05,'Axial force not converged to 5%%.');
    tempdir=tempname;mkdir(tempdir);cleanup=onCleanup(@()rmdir(tempdir,'s'));
    files=export_tf3d_shape(hi,tempdir,'test_design7');assert(exist(files.shape,'file')==2,'Export failed.');
    numeric=export_tf3d_results(hi,tempdir,'test_design7');
    raw=dlmread(numeric.field,',',1,0);
    assert(size(raw,1)==2*hi.field.n_samples*hi.design.n_turns,'Incomplete full-coil field export.');
    assert(sum(raw(:,14))==hi.field.n_samples*hi.design.n_turns,'Reconstruction flags failed.');

    data=dlmread(files.shape,',',1,0);
    assert(max(abs(data(1,:)-data(end,:)))<1e-6,'Export not closed.');
    assert(abs(data(1,11)/1e3-d.r1)<1e-10,'CAD units/reference failed.');
    invalid=hi;invalid.checks.shape_converged=false;blocked=false;
    try,export_tf3d_shape(invalid,tempdir,'invalid');catch err,blocked=strcmp(err.identifier,'tf3d:export');end
    assert(blocked,'Invalid shape export must be blocked.');
end
report.native_tests_passed=true;report.fem_validated=false;
fprintf('PASSED: native TF3D numerical tests. Arc-fit gate=%d. 3D FEM validation pending.\n',report.arc_fit_10mm_gate_pass);
end
function r=compute_outer_radius(p)
r=(p.R0+p.R0/p.A)*p.ripple^(-1/p.n_TF);
end
