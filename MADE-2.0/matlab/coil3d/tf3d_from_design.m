function out = tf3d_from_design(row,p,opts)
%TF3D_FROM_DESIGN Downstream TF geometry and magnetic-field research module.
% out=tf3d_from_design(row,p,opts); overrides omit the tf3d_ prefix.
% No scan results, 2D mechanics or input currents are modified.
if nargin<3,opts=struct();end
clock=tic;o=tf3d_defaults(p,opts);d=tf3d_design_geometry(row,p);
shape=tf_bending_free_shape(d.r1,d.r2,max(401,2*o.n_seg+1));
if o.shape_mode==1,shape=tf3d_iterate_shape(shape,d,p,o);end
arcs=tf_three_arc_fit(shape.upper);
if o.use_arcs,upper=arcs.upper;else,upper=shape.upper;end
loop=tf3d_resample_loop(upper,o.n_seg);
% Reject folded inward offsets on the three arcs; geometric discretization
% checks below also cover the analytic curve.
if o.use_arcs && min(arcs.radii)+min(d.case_u_min,min(d.u-d.h/2))<=0
    error('tf3d:offset','An inward offset exceeds an arc radius.');
end
src=tf3d_coil_filaments(loop,d,p,o);
f=tf3d_evaluate_turns(loop,d,src,o);
if ~isstruct(row),row2=table2struct(row);else,row2=row;end
if ~isfield(row2,'B_TF'),g=compute_operating_params(p);row2.B_TF=g.B_PHI_TF;end
f2=compute_discrete_field_profile(row2,p);
% Scalar local uplift is an engineering ESTIMATE, not vector superposition
% or a rigorous 3D maximum. Its validity on curved legs remains unproven.
f.Bpeak_estimate=f.Bmag+repmat(f2.B_discrete-f2.B_center,f.n_samples,1);
f.full=tf3d_full_field(f,loop,d);
f.summary=tf3d_field_summary(f,d,loop);
P=[-d.v(:) d.y(:) zeros(d.n_turns,1)];
Bm=biot_savart_segments(P,src.A,src.C,src.I,o);
mid=sqrt(sum(Bm.^2,2))';
c=tf3d_toroidal_checks(src,d,p,o);
c.midplane_3D=mid;c.midplane_2D=f2.B_center;
c.midplane_max_relative_difference=max(abs(mid./f2.B_center-1));
c.arc_fit_pass=arcs.max_error<=o.arc_tol;c.shape_converged=shape.converged;
c.geometry_offset_pass=check_offsets(loop,d,o);
% Half vertical resultant equals the SUM of the two cut axial forces.
% T=Fz/2 is only their mean, not necessarily the inner-leg tension.
f.axial_cut_force_sum=f.half_force(3);
f.mean_axial_cut_force=f.half_force(3)/2;
f.centring_force_per_coil=-2*f.half_force(2);
g=compute_operating_params(p);
f.T_bf_MADE=0.5*g.k_bf*p.n_TF*d.NI^2*(4e-7*pi)/(2*pi);
f.axial_interpretation='Fz_half=T_inner+T_outer; half this is only their mean';
out=struct('design',d,'options',o,'shape',shape,'arcs',arcs,'loop',loop, ...
    'field',f,'field2D',f2,'checks',c,'validated',false, ...
    'status','RESEARCH: 3D FEM and discretization validation pending');
if o.energy,out.energy=tf3d_energy_estimate(loop,d,p);end
out.elapsed_seconds=toc(clock);
if ~shape.converged||~c.ampere_pass||~c.geometry_offset_pass||(o.use_arcs&&~c.arc_fit_pass)
    out.status='CHECK FAILED: inspect shape, offset, Ampere or arc-fit diagnostics';
end
fprintf('\nTF3D: %d turns, NI=%.6f MA.turn; r1=%.5f, r2=%.5f m\n',d.n_turns,d.NI/1e6,d.r1,d.r2);
fprintf('Arcs max deviation %.2f mm; Ampere error %.3g; edge ripple %.3f %%\n',arcs.max_error*1e3,c.ampere_rel_error,100*c.ripple_Bphi);
fprintf('Sampled centre max %.4f T; peak ESTIMATE %.4f T; mean cut force %.3f MN\n',max(f.Bmag(:)),max(f.Bpeak_estimate(:)),f.mean_axial_cut_force/1e6);
fprintf('%s (%.1f s)\n',out.status,out.elapsed_seconds);
end
function ok=check_offsets(loop,d,o)
% Positive segment direction detects local offset inversion. This is not
% a CAD solid-interference test of the finite-width cases of all coils.
back=d.wp_u_max+d.ground_insulation+max(d.nose_il,o.nose_ol_factor*d.nose_il);
back=max(back,d.case_u_max);
ok=true;
for u=[d.case_u_min back]
    Q=loop.rz+u*loop.normal;D=diff(Q);
    if any(Q(:,1)<=0)||any(sum(D.*loop.drz,2)<=0),ok=false;end
end
end
