function [B,info] = biot_savart_segments(P,A,C,I,opts)
%BIOT_SAVART_SEGMENTS Exact finite straight segments, bounded working memory.
% P: Mx3, A/C: Sx3 segment starts/ends, I: scalar or Sx1 amperes.
% Points ON a source segment are singular and raise an error. Collinear
% points outside the segment have zero contribution. No self-field cutoff.
if nargin<5,opts=struct();end
pb=128;sb=512;
if isfield(opts,'point_block'),pb=opts.point_block;end
if isfield(opts,'source_block'),sb=opts.source_block;end
if size(P,2)~=3||size(A,2)~=3||~isequal(size(A),size(C))|| ...
        any(~isfinite([P(:);A(:);C(:)]))||pb<1||sb<1||pb~=round(pb)||sb~=round(sb)
    error('tf3d:segments','Invalid coordinates or block sizes.');
end
S=size(A,1);M=size(P,1);I=I(:);
if isscalar(I),I=repmat(I,S,1);end
if numel(I)~=S||any(~isfinite(I)),error('tf3d:segments','Invalid currents.');end
D=C-A;L=sqrt(sum(D.^2,2));
if any(L==0),error('tf3d:segments','Zero-length source segment.');end
E=D./L;B=zeros(M,3);info=struct('singular_count',0,'n_segments',S,'n_points',M);
for i=1:pb:M
    ip=i:min(i+pb-1,M);Q=P(ip,:);V=zeros(numel(ip),3);
    for j=1:sb:S
        js=j:min(j+sb-1,S);
        rx=Q(:,1)-A(js,1)';ry=Q(:,2)-A(js,2)';rz=Q(:,3)-A(js,3)';
        ex=E(js,1)';ey=E(js,2)';ez=E(js,3)';
        s=rx.*ex+ry.*ey+rz.*ez;
        % e cross (P-A), avoids cancellation in endpoint-vector formula.
        cx=ey.*rz-ez.*ry;cy=ez.*rx-ex.*rz;cz=ex.*ry-ey.*rx;
        rho2=cx.^2+cy.^2+cz.^2;
        ell=L(js)';q=s-ell;
        % A length-relative tolerance only identifies exact collinearity.
        near=rho2<=(64*eps*max(1,max(L(js))))^2;
        on=near & s>=0 & q<=0;
        if any(on(:)),error('tf3d:singular','Evaluation point lies on a source filament. Use a finite cross-section.');end
        ra=sqrt(rho2+s.^2);rb=sqrt(rho2+q.^2);
        den=rho2;den(near)=1;
        f=(s./ra-q./rb)./den;
        % Rationalized form is stable outside a segment, far along its axis.
        outside=s.*q>0 & ~near;
        stable=ell.*(s+q)./(ra.*rb.*(s.*rb+q.*ra));
        f(outside)=stable(outside);f(near)=0;
        f=1e-7*f.*I(js)';
        V=V+[sum(f.*cx,2) sum(f.*cy,2) sum(f.*cz,2)];
    end
    B(ip,:)=V;
end
end
