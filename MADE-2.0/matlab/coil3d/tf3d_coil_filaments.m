function src = tf3d_coil_filaments(loop,d,p,o)
%TF3D_COIL_FILAMENTS Closed planar source circuits from actual MADE turns.
% Coil 1 uses four cable quadrature filaments per turn, to evaluate a finite
% centre field including curved self-contributions (rectangular 2x2 Gauss
% approximation, not a validated rounded-cable 3D self-field model).
% Neighbours: 0 -> all turn quadratures; 1/2 -> 1/2 sources per layer.
% Two layer points preserve the actual turn-centre toroidal second moment.
A={};C={};I={};owner={};nc=0;
for coil=1:p.n_TF
    if coil==1||o.neighbor_gauss==0
        u=[];v=[];w=[];
        for a=[-1 1]
            for b=[-1 1]
                u=[u d.u+a*d.h/(2*sqrt(3))]; %#ok<AGROW>
                v=[v d.v+b*d.w/(2*sqrt(3))]; %#ok<AGROW>
                w=[w repmat(d.Iop/4,1,d.n_turns)]; %#ok<AGROW>
            end
        end
    else
        u=[];v=[];w=[];
        for k=1:max(d.layer)
            ids=find(d.layer==k);uk=mean(d.u(ids));
            if o.neighbor_gauss==1||numel(ids)==1
                u=[u uk];v=[v mean(d.v(ids))];w=[w d.Iop*numel(ids)]; %#ok<AGROW>
            else
                spread=sqrt(mean(d.v(ids).^2)+mean(d.w(ids).^2)/12);
                u=[u uk uk];v=[v -spread spread];w=[w d.Iop*numel(ids)/2 d.Iop*numel(ids)/2]; %#ok<AGROW>
            end
        end
    end
    phi=pi/2+(coil-1)*2*pi/p.n_TF;er=[cos(phi) sin(phi) 0];et=[-sin(phi) cos(phi) 0];
    for q=1:numel(u)
        rz=loop.rz+u(q)*loop.normal;
        if any(rz(:,1)<=0),error('tf3d:offset','A filament crosses the machine axis.');end
        X=rz(:,1)*er+v(q)*et+[zeros(size(rz,1),2) rz(:,2)];
        nc=nc+1;A{nc}=X(1:end-1,:);C{nc}=X(2:end,:); %#ok<AGROW>
        I{nc}=repmat(w(q),loop.n_seg,1);owner{nc}=repmat(coil,loop.n_seg,1); %#ok<AGROW>
    end
end
src=struct('A',vertcat(A{:}),'C',vertcat(C{:}),'I',vertcat(I{:}), ...
    'coil',vertcat(owner{:}),'n_filaments',nc,'model','finite cable quadrature / grouped neighbours');
end
