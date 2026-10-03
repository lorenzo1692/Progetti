function figs = plot_tf3d(out,p)
%PLOT_TF3D Full turn geometry and available EM outputs, with model labels.
L=out.loop;d=out.design;f=out.field;a=out.arcs;
figs(1)=figure('Name','TF3D geometry','NumberTitle','off');
plot(out.shape.upper(:,1),out.shape.upper(:,2),'k-',a.upper(:,1),a.upper(:,2),'r--');hold on
for u=[d.wp_u_min d.wp_u_max d.case_u_min]
    Q=L.rz+u*L.normal;plot(Q(:,1),Q(:,2));
end
% Use the same varying nose as the export, including the outer-leg factor.
u=d.wp_u_max+d.ground_insulation+d.nose_il* ...
    (1+(out.options.nose_ol_factor-1)*(L.rz(:,1)-d.r1)/(d.r2-d.r1));
Q=L.rz+u.*L.normal;plot(Q(:,1),Q(:,2));
axis equal;grid on;xlabel('R [m]');ylabel('Z [m]');
title(sprintf('Current-line fit: max deviation %.1f mm; use_arcs=%d',1e3*a.max_error,out.options.use_arcs));
legend('Free span','Three arcs','WP plasma','WP back','Case plasma','Case back (variable nose)');
figs(2)=figure('Name','TF3D full winding geometry','NumberTitle','off');hold on
% Draw every real turn of every coil. NaN-separated polylines keep the
% graphics-object count small. Crossovers and leads are not modelled.
for coil=1:p.n_TF
    ph=pi/2+(coil-1)*2*pi/p.n_TF;er=[cos(ph) sin(ph) 0];et=[-sin(ph) cos(ph) 0];
    X=nan((size(L.rz,1)+1)*d.n_turns,3);
    for j=1:d.n_turns
        rz=L.rz+d.u(j)*L.normal;
        ix=(j-1)*(size(L.rz,1)+1)+(1:size(L.rz,1));
        X(ix,:)=rz(:,1)*er+d.v(j)*et+[zeros(size(rz,1),2) rz(:,2)];
    end
    plot3(X(:,1),X(:,2),X(:,3),'Color',[.55 .62 .7],'LineWidth',.3);
end
if isfield(f,'full'),full=f.full;else,full=tf3d_full_field(f,L,d);end
scatter3(full.points(:,1),full.points(:,2),full.points(:,3),10,full.Bmag,'filled');
axis equal;grid on;view(3);xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]');
c=colorbar;ylabel(c,'Centre |B| [T], coil 1');
title({sprintf('%d coils x %d real turns',p.n_TF,d.n_turns), ...
    'Coil 1 lower-half field reconstructed by Z symmetry'});
figs(3)=figure('Name','TF3D field along complete turn','NumberTitle','off');
subplot(2,1,1);
plot([f.s;L.length-flipud(f.s)],[f.Bpeak_estimate;flipud(f.Bpeak_estimate)]);hold on
plot([0 L.length],max(out.field2D.B_discrete)*[1 1],'k--','LineWidth',1.5);
xlabel('Full reference-loop arc length [m]');ylabel('Estimated peak [T]');grid on
 title('Scalar 2D uplift estimate: not a resolved 3D surface maximum');
subplot(2,1,2);plot(1:d.n_turns,out.checks.midplane_2D,'k-',1:d.n_turns,out.checks.midplane_3D,'r--');
xlabel('Turn');ylabel('Centre |B| [T]');legend('2D','3D at z=0');grid on
figs(4)=figure('Name','TF3D ripple and radial profile','NumberTitle','off');
subplot(2,1,1);plot(out.checks.edge_phi*180/pi,out.checks.edge_Bphi,'LineWidth',1.4);grid on
xlabel('Toroidal angle [deg]');ylabel('B_phi [T]');
title(sprintf('Plasma-edge ripple %.3f%%; target %.3f%%',100*out.checks.ripple_Bphi,100*p.ripple));
if isfield(out.checks,'radial_R')
    subplot(2,1,2);plot(out.checks.radial_R,out.checks.radial_Bphi,'LineWidth',1.2);hold on
    plot(out.checks.radial_R,out.checks.radial_ampere_mean,'k--');grid on
    xlabel('R [m], plasma aperture');ylabel('B_phi [T]');legend('Coil plane','Between coils','Ampere toroidal mean');
end
figs(5)=figure('Name','TF3D Lorentz loads','NumberTitle','off');
subplot(1,2,1);quiver3(full.points(:,1),full.points(:,2),full.points(:,3), ...
    full.force(:,1),full.force(:,2),full.force(:,3),.6);
axis equal;grid on;view(3);xlabel('X [m]');ylabel('Y [m]');zlabel('Z [m]');
title('Lorentz force per integration bin; arrows autoscaled');
subplot(1,2,2);
if isfield(f,'summary'),s=f.summary;else,s=tf3d_field_summary(f,d,L);end
plot(s.turn,-s.full_turn_force(:,2)/1e3,'LineWidth',1.2);grid on
xlabel('Turn');ylabel('Net inward force [kN]');
title({sprintf('Coil centring force %.3f MN',f.centring_force_per_coil/1e6), ...
    sprintf('Mean axial cut force %.3f MN (not T_inner)',f.mean_axial_cut_force/1e6)},'Interpreter','none');
if isfield(out,'energy')
    figs(6)=figure('Name','TF3D approximate inductance','NumberTitle','off');
    imagesc(out.energy.centroid_loop_matrix);axis image;c=colorbar;ylabel(c,'H per equivalent NI loop');
    xlabel('Coil');ylabel('Coil');
    title({sprintf('Approximate W %.3f GJ; L_eq %.3f H',out.energy.W/1e9,out.energy.L_eq_per_coil), ...
        'Centroid/GMD approximation: not validated for discharge sizing'});
end
end
