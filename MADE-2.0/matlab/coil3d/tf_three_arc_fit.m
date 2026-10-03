function fit = tf_three_arc_fit(upper)
%TF_THREE_ARC_FIT Deterministic C1 three-arc approximation of the upper span.
% End positions and vertical tangents are exact. Sweep angles 0<a<b<pi
% and the middle radius parameterize the fit; the other radii solve two
% endpoint constraints. No random starts, logical residuals or open loop.
A=upper(end,:);B=upper(1,:);span=B(1)-A(1);
s=[0;cumsum(sqrt(sum(diff(upper).^2,2)))];
ref=interp1(s,upper,linspace(0,s(end),160)','pchip');
best=inf;bestx=[];
opts=optimset('Display','off','MaxIter',500,'MaxFunEvals',1200,'TolX',1e-8,'TolFun',1e-12);
for a=[0.5 0.9 1.2]
    for b=[1.9 2.3 2.7]
        for frac=[0.4 0.9 1.5]
            x0=[log(a/(pi-a)),log((b-a)/(pi-b)),log(frac)];
            [x,f]=fminsearch(@objective,x0,opts);
            if f<best,best=f;bestx=x;end
        end
    end
end
[ok,R,angles,centres,joins]=build(bestx);
if ~ok,error('tf3d:arcs','No admissible three-arc fit found.');end
pts=[];
for k=1:3
    al=linspace(angles(k),angles(k+1),101)';
    q=centres(k,:)+R(k)*[-cos(al) sin(al)];
    if k>1,q=q(2:end,:);end
    pts=[pts;q]; %#ok<AGROW>
end
% Bidirectional sampled point-to-curve check, not only the optimizer RMS.
d1=distance_to_arcs(ref,R,angles,centres,joins);
d2=point_polyline_distance(pts,upper);
fit=struct('upper',flipud(pts),'radii',R,'angles',angles, ...
    'centres',centres,'joins',joins,'rms_error',sqrt(mean(d1.^2)), ...
    'max_error',max([d1;d2]),'optimizer_objective',best);
    function f=objective(x)
        [valid,rr,aa,cc,jj]=build(x);
        if ~valid,f=1e4+sum(min(abs(x),100).^2);return;end
        dd=distance_to_arcs(ref,rr,aa,cc,jj);f=mean((dd/span).^2);
    end
    function [valid,rr,aa,cc,jj]=build(x)
        valid=false;rr=[];aa=[];cc=[];jj=[];
        if isempty(x)||any(~isfinite(x))||any(abs(x)>30),return;end
        a=pi/(1+exp(-x(1)));b=a+(pi-a)/(1+exp(-x(2)));
        aa=[0 a b pi];mid=span*exp(x(3));
        cr=cos(aa(1:3))-cos(aa(2:4));sz=sin(aa(2:4))-sin(aa(1:3));
        M=[cr([1 3]);sz([1 3])];
        if rcond(M)<1e-8,return;end
        ends=M\([span;-A(2)]-mid*[cr(2);sz(2)]);
        rr=[ends(1) mid ends(2)];
        if any(rr<=span*1e-4)||any(rr>100*span),return;end
        cc=zeros(3,2);jj=zeros(4,2);jj(1,:)=A;
        for k=1:3
            cc(k,:)=jj(k,:)+rr(k)*[cos(aa(k)) -sin(aa(k))];
            jj(k+1,:)=jj(k,:)+rr(k)*[cr(k) sz(k)];
        end
        valid=true;
    end
end
function d=distance_to_arcs(P,R,aa,C,J)
d=inf(size(P,1),1);
for k=1:3
    v=P-C(k,:);al=atan2(v(:,2),-v(:,1));al=mod(al,2*pi);
    inside=al>=aa(k)&al<=aa(k+1);
    dd=min(sqrt(sum((P-J(k,:)).^2,2)),sqrt(sum((P-J(k+1,:)).^2,2)));
    dd(inside)=abs(sqrt(sum(v(inside,:).^2,2))-R(k));d=min(d,dd);
end
end
function d=point_polyline_distance(P,Q)
d=inf(size(P,1),1);
for i=1:size(Q,1)-1
    v=Q(i+1,:)-Q(i,:);t=max(0,min(1,(P-Q(i,:))*v'/sum(v.^2)));
    d=min(d,sqrt(sum((P-Q(i,:)-t*v).^2,2)));
end
end
