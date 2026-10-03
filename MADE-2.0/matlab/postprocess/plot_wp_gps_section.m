function fig = plot_wp_gps_section(out, fig_title, deformation_scale)
%PLOT_WP_GPS_SECTION Dedicated generalized-plane-strain FEM section figure.
% Uses the actual Q8/T6 mesh and solved nodal displacements, not sizing boxes.
% deformation_scale is optional; the automatic factor is reported explicitly.
if nargin<2||isempty(fig_title),fig_title='FEM surrogate: 2D GPS section';end
m=out.mesh;xy=m.xy;u=reshape(out.sol.U(:),2,[])';
if size(u,1)~=size(xy,1),error('plot_wp_gps_section:dofs','Displacement/mesh size mismatch.');end
umag=sqrt(sum(u.^2,2));
if nargin<3||isempty(deformation_scale)
    span=max(max(xy,[],1)-min(xy,[],1));
    deformation_scale=min(100,0.05*span/max(max(umag),eps));
end
if ~isscalar(deformation_scale)||~isfinite(deformation_scale)||deformation_scale<0
    error('plot_wp_gps_section:scale','Deformation scale must be finite and nonnegative.');
end
fig=figure('Name','FEM surrogate - 2D GPS section','NumberTitle','off', ...
    'Units','normalized','OuterPosition',[0.03 0.05 0.94 0.9]);
status='INVALID';if out.valid,status='Numerical checks passed';end
ax=subplot(2,2,1);hold(ax,'on');
paint(ax,m,xy,m.q8_mat(:),m.t6_mat(:),'flat',true);
colormap(ax,lines(7));set(ax,'CLim',[0.5 7.5]);colorbar(ax);
title(ax,{fig_title,'Material IDs: 1 LTS, 2 jacket, 3 insulation, 4 case, 5 wedge, 6 HTS, 7 filler'},'Interpreter','none');
ax=subplot(2,2,2);hold(ax,'on');
paint(ax,m,xy,out.elem_SINT.q8(:)/1e6,out.elem_SINT.t6(:)/1e6,'flat',false);
colormap(ax,jet(256));cb=colorbar(ax);ylabel(cb,'Tresca [MPa]');
title(ax,{sprintf('GPS total load case: eps_z = %.5g',out.sol.eps_z),status},'Interpreter','none');
ax=subplot(2,2,3);hold(ax,'on');
% Original mesh outline provides an explicit reference for the scaled shape.
wire(ax,m,xy,[.7 .7 .7]);
paint(ax,m,xy+deformation_scale*u,umag*1e3,umag*1e3,'interp',false);
colormap(ax,jet(256));cb=colorbar(ax);ylabel(cb,'|u| [mm], physical values');
title(ax,sprintf('Deformed section x %.3g; max |u| = %.4g mm',deformation_scale,max(umag)*1e3));
ax=subplot(2,2,4);hold(ax,'on');
paint(ax,m,xy,u(:,2)*1e3,u(:,2)*1e3,'interp',false);
colormap(ax,jet(256));cb=colorbar(ax);ylabel(cb,'u_y radial [mm]');
cav=out.geo.cav;
xlim(ax,[min(cav(:,1)) max(cav(:,1))]+[-.01 .01]);
ylim(ax,[min(cav(:,2)) max(cav(:,2))]+[-.01 .01]);
title(ax,'WP zoom: radial displacement on the undeformed mesh');
end
function paint(ax,m,xy,qval,tval,mode,edges)
edge='none';if edges,edge=[.3 .3 .3];end
if ~isempty(m.q8)
    patch(ax,'Faces',m.q8(:,[1 5 2 6 3 7 4 8]),'Vertices',xy, ...
        'FaceVertexCData',qval,'FaceColor',mode,'EdgeColor',edge,'LineWidth',.15);
end
if ~isempty(m.t6)
    patch(ax,'Faces',m.t6(:,[1 4 2 5 3 6]),'Vertices',xy, ...
        'FaceVertexCData',tval,'FaceColor',mode,'EdgeColor',edge,'LineWidth',.15);
end
axis(ax,'equal');box(ax,'on');xlabel(ax,'Toroidal x [m]');ylabel(ax,'Radial y [m]');
end
function wire(ax,m,xy,color)
if ~isempty(m.q8),patch(ax,'Faces',m.q8(:,[1 5 2 6 3 7 4 8]),'Vertices',xy,'FaceColor','none','EdgeColor',color,'LineWidth',.1);end
if ~isempty(m.t6),patch(ax,'Faces',m.t6(:,[1 4 2 5 3 6]),'Vertices',xy,'FaceColor','none','EdgeColor',color,'LineWidth',.1);end
end
