function loop = tf3d_resample_loop(upper,nseg)
%TF3D_RESAMPLE_LOOP Closed clockwise R-Z path, starts at inner mid-plane.
% Keeps both mid-plane cuts and straight/free-span junctions as vertices.
% upper is outboard mid-plane -> top of the inner straight leg.
A=upper(end,:);B=upper(1,:);
if abs(B(2))>1e-10||A(2)<=0,error('tf3d:shape','Invalid upper span endpoints.');end
s=[0;cumsum(sqrt(sum(diff(upper).^2,2)))];
ns=max(8,round((nseg/2)*s(end)/(s(end)+A(2))));
nl=max(2,round(nseg/2)-ns);
U=interp1(s,upper,linspace(0,s(end),ns+1)','pchip');
U(1,:)=B;U(end,:)=A;
line=[repmat(A(1),nl+1,1) linspace(0,A(2),nl+1)'];
upperpath=[line;flipud(U(1:end-1,:))];
lowerpath=flipud(upperpath);lowerpath(:,2)=-lowerpath(:,2);
P=[upperpath;lowerpath(2:end,:)];
D=diff(P);L=sqrt(sum(D.^2,2));
if any(L<=0),error('tf3d:shape','Degenerate loop segment.');end
T=D./L;
Tv=T+[T(end,:);T(1:end-1,:)];Tv=Tv./sqrt(sum(Tv.^2,2));
N=[-Tv(:,2) Tv(:,1)];
loop=struct('rz',P,'normal',[N;N(1,:)],'s',[0;cumsum(L)], ...
    'length',sum(L),'n_seg',numel(L),'n_half',size(upperpath,1)-1);
loop.mid=0.5*(P(1:end-1,:)+P(2:end,:));
loop.drz=D;loop.ds=L;
end
