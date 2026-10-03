function shape = tf_bending_free_shape(r1,r2,n,Bfun)
%TF_BENDING_FREE_SHAPE Upper free span, ordered outboard -> inboard.
% Bfun(r) must give POSITIVE effective normal field [T]. Without Bfun,
% use B=1/r; only geometry is normalized, not a prediction of tension.
if r1<=0||r2<=r1||n<16,error('tf3d:shape','Invalid radii or resolution.');end
t=linspace(-pi/2,pi/2,n)';
if nargin<4||isempty(Bfun)
    k=0.5*log(r2/r1);r=sqrt(r1*r2)*exp(-k*sin(t));
    B=1./r;intB=log(r2/r1);model='analytic Princeton D';
else
    rr=linspace(r1,r2,max(2001,8*n))';bb=Bfun(rr);bb=bb(:);
    if numel(bb)~=numel(rr)||any(~isfinite(bb)|bb<=0)
        error('tf3d:shape','Effective field must be positive and finite.');
    end
    F=cumtrapz(rr,bb);k=F(end)/2;intB=F(end);
    r=interp1(F,rr,k*(1-sin(t)),'pchip');B=Bfun(r);B=B(:);
    model='effective-field iteration';
end
z=cumtrapz(t,-k*sin(t)./B);
if z(end)<=0||any(z< -1e-10),error('tf3d:shape','Invalid free span: positive inner-leg height required.');end
shape=struct('upper',[r z],'t',t,'integral_B',intB,'model',model, ...
    'converged',true,'iterations',0,'relative_change',0);
end
