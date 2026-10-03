function summary = tf3d_field_summary(f,d,loop)
%TF3D_FIELD_SUMMARY Per-turn maxima and full-coil force integrals.
% Maxima refer to sampled centres and scalar peak estimates, respectively.
[bc,ic]=max(f.Bmag,[],1);[bp,ip]=max(f.Bpeak_estimate,[],1);
nt=d.n_turns;nb=f.n_samples;pc=zeros(nt,3);pp=pc;force=pc;
for j=1:nt
    ids=(j-1)*nb+(1:nb);
    pc(j,:)=f.points(ids(ic(j)),:);pp(j,:)=f.points(ids(ip(j)),:);
    hf=sum(f.force(ids,:),1);force(j,:)=[2*hf(1:2) 0];
end
summary=struct('turn',(1:nt)','layer',d.layer(:), ...
    'centre_max',bc(:),'centre_max_xyz',pc,'centre_max_s',f.s(ic), ...
    'peak_estimate',bp(:),'peak_estimate_xyz',pp,'peak_estimate_s',f.s(ip), ...
    'full_turn_force',force,'reference_length',loop.length);
summary.labels='centre_max: evaluated |B|; peak_estimate: unvalidated scalar 2D uplift';
end
