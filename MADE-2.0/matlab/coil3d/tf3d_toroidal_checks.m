function c = tf3d_toroidal_checks(src,d,p,o)
%TF3D_TOROIDAL_CHECKS Numerical Ampere mean and plasma-edge ripple.
% Ampere is exact in the continuum; sampled quadrature has a finite error.
n=16*p.n_TF;angles=(0:n-1)'*2*pi/n;
P=[p.R0*cos(angles) p.R0*sin(angles) zeros(n,1)];
B=biot_savart_segments(P,src.A,src.C,src.I,o);
bphi=-B(:,1).*sin(angles)+B(:,2).*cos(angles);
ref=2e-7*p.n_TF*d.NI/p.R0;
c.mean_Bphi=mean(bphi);c.ampere_reference=ref;c.input_B0=p.B0;
c.ampere_rel_error=abs(c.mean_Bphi/ref-1);
c.ampere_pass=c.ampere_rel_error<=o.ampere_tol;
% One periodic sector including its endpoints; report |B| and Bphi ripple.
a=linspace(pi/2,pi/2+2*pi/p.n_TF,129)';R=p.R0+p.R0/p.A;
B=biot_savart_segments([R*cos(a) R*sin(a) zeros(size(a))],src.A,src.C,src.I,o);
b=-B(:,1).*sin(a)+B(:,2).*cos(a);bm=sqrt(sum(B.^2,2));
c.edge_phi=a;c.edge_Bphi=b;c.edge_Bmag=bm;
c.ripple_Bphi=(max(b)-min(b))/(max(b)+min(b));
c.ripple_Bmag=(max(bm)-min(bm))/(max(bm)+min(bm));
c.ripple_target=p.ripple;c.ripple_pass=c.ripple_Bphi<=p.ripple;
% Radial profile through the plasma aperture at coil and inter-coil planes.
r=linspace(p.R0-p.R0/p.A,p.R0+p.R0/p.A,81)';
phis=[pi/2 pi/2+pi/p.n_TF];profile=zeros(numel(r),2);
for k=1:2
    ph=phis(k);pts=[r*cos(ph) r*sin(ph) zeros(size(r))];
    br=biot_savart_segments(pts,src.A,src.C,src.I,o);
    profile(:,k)=-br(:,1)*sin(ph)+br(:,2)*cos(ph);
end
c.radial_R=r;c.radial_Bphi=profile;c.radial_phi=phis;
c.radial_ampere_mean=2e-7*p.n_TF*d.NI./r;

end
