function e = tf3d_energy_estimate(loop,d,p)
%TF3D_ENERGY_ESTIMATE OPTIONAL coarse current-centroid Neumann estimate.
% One loop per coil carrying NI. Same-coil kernel softened by WP GMD.
% This approximates the whole pack as a rectangular section; it is NOT the
% inductance of the individual turn circuits. No discharge sizing feedback.
ns=loop.n_seg;Lcoil=zeros(p.n_TF);q=loop.rz;
gmd=0.2235*(d.wp_width+d.wp_depth);
P=cell(p.n_TF,1);D=P;
for c=1:p.n_TF
    ph=pi/2+(c-1)*2*pi/p.n_TF;
    X=[q(:,1)*cos(ph) q(:,1)*sin(ph) q(:,2)];
    P{c}=(X(1:end-1,:)+X(2:end,:))/2;D{c}=diff(X);
end
for c=1:p.n_TF
    for b=c:p.n_TF
        val=0;
        for i=1:128:ns
            ix=i:min(i+127,ns);dx=P{c}(ix,1)-P{b}(:,1)';dy=P{c}(ix,2)-P{b}(:,2)';dz=P{c}(ix,3)-P{b}(:,3)';
            soft=0;if c==b,soft=gmd;end
            R=sqrt(dx.^2+dy.^2+dz.^2+soft^2);
            dotdl=D{c}(ix,:)*D{b}';term=dotdl./R;val=val+sum(term(:));
        end
        Lcoil(c,b)=1e-7*val;Lcoil(b,c)=Lcoil(c,b);
    end
end
W=0.5*d.NI^2*sum(Lcoil(:));
e=struct('W',W,'L_eq_per_coil',2*W/(d.Iop^2*p.n_TF), ...
    'centroid_loop_matrix',Lcoil,'gmd',gmd,'validated',false, ...
    'model','approximate centroid loops with softened WP GMD; not for discharge sizing');
end
