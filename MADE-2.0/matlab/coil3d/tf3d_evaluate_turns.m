function f = tf3d_evaluate_turns(loop,d,src,o)
%TF3D_EVALUATE_TURNS Centre fields and midpoint Lorentz forces, upper half.
% Sampling each bin once limits cost; ds weights cover the full half loop.
% Increasing n_eval to loop.n_half gives every upper segment midpoint.
nh=loop.n_half;nb=min(nh,o.n_eval);
edges=round(linspace(0,nh,nb+1));idx=zeros(nb,1);weights=zeros(nb,1);
for k=1:nb
    ids=(edges(k)+1):edges(k+1);idx(k)=ids(ceil(numel(ids)/2));
    weights(k)=sum(loop.ds(ids));
end
P=zeros(nb*d.n_turns,3);dl=zeros(size(P));
for j=1:d.n_turns
    rz=loop.rz+d.u(j)*loop.normal;
    Q=[-d.v(j)*ones(size(rz,1),1) rz(:,1) rz(:,2)];
    ids=(j-1)*nb+(1:nb);
    P(ids,:)=0.5*(Q(idx,:)+Q(idx+1,:));
    D=Q(idx+1,:)-Q(idx,:);
    dl(ids,:)=D.*(weights./loop.ds(idx));
end
B=biot_savart_segments(P,src.A,src.C,src.I,o);
force=d.Iop*cross(dl,B,2);
f=struct('points',P,'B',B,'Bmag',reshape(sqrt(sum(B.^2,2)),nb,d.n_turns), ...
    'Bphi',reshape(-B(:,1),nb,d.n_turns),'force',force,'idx',idx,'weights',weights, ...
    's',0.5*(loop.s(idx)+loop.s(idx+1)), ...
    'half_force',sum(force,1),'n_samples',nb);
f.normal_field=sum(reshape(-B(:,1).*sqrt(sum(dl.^2,2)),nb,d.n_turns),2)./(weights*d.n_turns);
end
