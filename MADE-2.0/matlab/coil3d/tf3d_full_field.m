function full = tf3d_full_field(f,loop,d)
%TF3D_FULL_FIELD Complete coil-1 samples using exact mid-plane symmetry.
% Source geometry is symmetric about Z=0, with reversed current circulation
% under reflection: Bx/By are even, Bz odd; Fx/Fy even, Fz odd.
% Lower samples are reconstructed, NOT additional Biot-Savart evaluations.
n=f.n_samples;nt=d.n_turns;
full.points=zeros(2*n*nt,3);full.B=full.points;full.force=full.points;
full.turn=zeros(2*n*nt,1);full.layer=full.turn;full.s=full.turn;
full.reconstructed=false(2*n*nt,1);
for j=1:nt
    src=(j-1)*n+(1:n);dst=(j-1)*2*n+(1:2*n);
    P=f.points(src,:);B=f.B(src,:);F=f.force(src,:);
    Plo=flipud(P);Plo(:,3)=-Plo(:,3);
    Blo=flipud(B);Blo(:,3)=-Blo(:,3);
    Flo=flipud(F);Flo(:,3)=-Flo(:,3);
    full.points(dst,:)=[P;Plo];full.B(dst,:)=[B;Blo];full.force(dst,:)=[F;Flo];
    full.turn(dst)=j;full.layer(dst)=d.layer(j);
    full.s(dst)=[f.s;loop.length-flipud(f.s)];
    full.reconstructed(dst(n+1:end))=true;
end
full.Bmag=sqrt(sum(full.B.^2,2));
full.force_units='N per represented integration bin, not N/m';
full.method='Upper half evaluated; lower half reconstructed by Z symmetry';
full.total_force=sum(full.force,1);
end
