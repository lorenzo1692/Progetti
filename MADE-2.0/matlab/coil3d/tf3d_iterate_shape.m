function shape = tf3d_iterate_shape(shape,d,p,o)
%TF3D_ITERATE_SHAPE Experimental fixed point on cross-section centre forces.
% This is NOT a structural equilibrium or a proof of zero bending moment.
% A failed iteration is retained with converged=false; never silently passed.
history=zeros(o.max_iter,2);shape.converged=false;
for it=1:o.max_iter
    loop=tf3d_resample_loop(shape.upper,o.n_seg);
    src=tf3d_coil_filaments(loop,d,p,o);
    f=tf3d_evaluate_turns(loop,d,src,o);
    rm=loop.mid(f.idx,1);
    free=rm>d.r1+1e-8*(d.r2-d.r1);
    [r,ix]=sort(rm(free));bf=f.normal_field(free);bf=bf(ix);
    [r,ia]=unique(r);bf=bf(ia);
    if numel(r)<8||any(bf<=0|~isfinite(bf))
        error('tf3d:shape_field','Insufficient positive normal-field samples for iteration.');
    end
    rr=[d.r1;r;d.r2];bb=[bf(1);bf;bf(end)];
    proposal=tf_bending_free_shape(d.r1,d.r2,size(shape.upper,1),@(x)interp1(rr,bb,x,'pchip'));
    residual=max(sqrt(sum((proposal.upper-shape.upper).^2,2)))/(d.r2-d.r1);
    history(it,:)=[it residual];
    fprintf('TF3D shape %d: unrelaxed relative update %.3g\n',it,residual);
    shape.upper=(1-o.shape_relax)*shape.upper+o.shape_relax*proposal.upper;
    shape.integral_B=proposal.integral_B;
    if residual<=o.shape_tol,shape.converged=true;break;end
end
shape.iterations=it;shape.relative_change=residual;shape.history=history(1:it,:);
shape.model='experimental finite-cable effective-field fixed point';
shape.tension_effective=d.NI*shape.integral_B/2;
end
