clearvars; clc; close all;

%% D constant tension (guscio magnetico uniforme)

%%% AGGIORNARE LA SEZIONE DEL WP NEL FILE MAIN.dat %%%

printcond = 1;

TitleName0 = 'Three Circles with Tangent condition';
name_file  = sprintf('TF shape - VNS 07 2026 - %s.txt',today("datetime"));
TitleName  = sprintf('TF shape - VNS 07 2026 - %s',today("datetime"));
TitleName2 = sprintf('TF shape - VNS 07 2026 - %s ',today("datetime"));

%% VNS 2024 PF-in-TF
% R0 = 2.67;    
% B0 = 5.6; 
% A_V3 = 4.25;
% wb = 0.86;
% Ri = R0-R0/A_V3-wb; 
% Mu_0 = 4*pi*10^(-7);
% ripple = 0.01;
% N_TF = 12;
% Re = 5.15; % (R0+R0/A_V3)*(1+1/ripple)^(1/N_TF);  
% NI = (2*pi*R0*B0 /Mu_0)/N_TF;  

%% VNS 2026 PF-out-TF
R0 = 2.83;    
B0 = 5.6; 
A_V3 = 4.35;
wb = 0.92;
Ri = R0-R0/A_V3-wb; 
Mu_0 = 4*pi*10^(-7);
ripple = 0.012;
N_TF = 12;
Re =  (R0+R0/A_V3)*(1+1/ripple)^(1/N_TF);  
NI = (2*pi*R0*B0 /Mu_0)/N_TF;  

dz_0 = 0;
dr_plasma_side = 0.02;
WP_h = 0.38;
nose_il = 0.13;
nose_ol = nose_il/2;

%% DEMO LAR 2024 baseline
% R0 = 8.6;       
% B0 = 4.389;  
% A = 2.8;
% wb = 1.82;
% Ri = R0-R0/A-wb;  
% Mu_0 = 4*pi*10^-7;
% ripple = 0.006;
% N_TF = 16;
% Re = (R0+R0/A)*(1+1/ripple)^(1/N_TF);   
% NI = (2*pi*R0*B0 /Mu_0)/N_TF;       

%% VNS 2024 WF
% R0 = 2.67;    
% B0 = 5.6; 
% A_ = 4.25;
% wb = 0.86;
% Ri = R0-R0/A_-wb; 
% Mu_0 = 4*pi*10^(-7);
% ripple = 0.01;
% N_TF = 12;
% Re = 4.77;%(R0+R0/A_)*(1+1/ripple)^(1/N_TF);  
% ripple_ = 1/((Re/(R0+R0/A_))^N_TF-1);
% NI = (2*pi*R0*B0 /Mu_0)/N_TF; 

%% PROTO
% R0 = 6.7;       
% B0 = 4.7;  
% Ar = 2.8;
% N_TF = 16;    % TF number
% wb = 1.48;     % SB+VV+gap  
% rho = 0.01;  % ripple
% Mu_0 = 4*pi*10^(-7);
% a = R0/Ar;    % minor radius
% Ri = R0-a-wb;  % inner leg outer radius
% Re = (R0+a).*(1./rho).^(1/N_TF); % outer leg inner radius     
% NI = (2*pi*R0*B0 /Mu_0)/N_TF; 

%% CEFTR
% B0 = 6; % field on the axis
% R0 = 8; % major radius
% Ar  = 2.95; % aspect ratio
% N_TF = 16;% TF number
% d = 1.60;  % SB+VV+gap  
% Mu_0 = 4*pi*10^(-7);
% rho = 0.006; % ripple
% a = R0./Ar; % minor radius
% Ri = R0-a-d; % inner leg outer radius
% Re = (R0+a).*(1./rho).^(1/N_TF);
% % Re = 14.3;
% % ripple_ = 1/((Re/(R0+R0/Ar))^N_TF-1);
% NI = (2*pi*R0*B0 /Mu_0)/N_TF; 

% dz_0 = 0;
% dr_plasma_side = 0.06;
% WP_h = 0.6;
% nose_il = 0.3;
% nose_ol = nose_il/2;

%% PILOT B

% B0 = 5; % field on the axis
% R0 = 4; % major radius
% Ar  = 3.4; % aspect ratio
% N_TF = 16;% TF number
% d = 0.78;  % SB+VV+gap  
% Mu_0 = 4*pi*10^(-7);
% rho = 0.006; % ripple
% a = R0./Ar; % minor radius
% Ri = R0-a-d; % inner leg outer radius
% Re = (R0+a).*(1./rho).^(1/N_TF);
% ripple_ = 1/((Re/(R0+R0/Ar))^N_TF-1);
% NI = (2*pi*R0*B0 /Mu_0)/N_TF; 
% 
% dz_0 = 0;
% dr_plasma_side = 0.03;
% WP_h = 0.183;
% nose_il = 0.155;
% nose_ol = nose_il/2;

%% Bending Free TF shape
k = 0.5*log(Re/Ri);
r0 = Ri.*exp(k);
N = 100;
theta = linspace(pi/2, 3/2*pi, N);
dtheta = [0 diff(theta)];
r_BF = r0.*exp(k.*sin(theta));    
dz = (r0*k).*sin(theta).*exp(k.*sin(theta)).*dtheta;    
z = cumsum(dz);
dzeta = (min(z));
z_BF = z-dzeta;
z_il = linspace(-z(end),z(end),N);
r_il = zeros(1,N)+min(r_BF);
z_BF = [z_BF, -flip(z_BF), z_il]';
r_BF = [r_BF, flip(r_BF), r_il]';
% plot(r_BF,z_BF,'.')

%% Optimize bendingfree TF shape
delete(sprintf('B_.csv'))
delete(sprintf('radialMatlab_VNS.txt'))
delete(sprintf('verticalMatlab_VNS.txt'))
N = 600; % shape discretization - multipli di 3
[r, z, t] = bendingfree_opt(Ri, Re, N, N_TF, NI);

%% plot modified BF D-shape adding the innerleg
r = r(2:end,1);
z = z(2:end,1);
zcat = linspace(z(end),z(1),N/3);
zs = [z; zcat'; -flip(zcat'); -z];
rs = [r; zcat'.*0+r(1); flip(zcat'.*0+r(1)); r];

fileID = fopen('r_v.txt','w');
fprintf(fileID, '%f \n',rs);
fclose(fileID);
fileID = fopen('z_v.txt','w');
fprintf(fileID, '%f \n',zs);
fclose(fileID);

figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
plot(rs,zs,'.')
ylabel('z [m]','fontsize',20,'Interpreter','latex')
xlabel('r [m]','fontsize',20,'Interpreter','latex')
set(gca,'FontSize',20);
axis equal
grid on

%%
A = readmatrix('B_.csv');
Bx = (Mu_0*N_TF*NI./(2*pi.*A(:,2)));
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
plot(A(:,7))
hold on
plot(A(:,7).*0+Bx./2,'r')
set(gca,'FontSize',20);
ylabel('B [T]','fontsize',20,'Interpreter','latex')
xlabel('curv. abscissa','fontsize',20,'Interpreter','latex')
legend('FEM','ideal')
%
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
plot(A(:,2),A(:,7))
hold on
plot(A(:,2),A(:,7).*0+Bx/2,'r')
set(gca,'FontSize',20);
ylabel('B [T]','fontsize',20,'Interpreter','latex')
xlabel('m','fontsize',20,'Interpreter','latex')
legend('FEM','ideal')

%% Approximate the optimized BF shape with a set of tangent arcs
xA = min(r); yA = z(1);
xB = max(r); yB = 0;
max_distance = 1e9;

while max_distance >= 0.40

    % Funzione che restituisce le equazioni da risolvere
    fun = @(vars) [
        vars(4) - yA;
        vars(10) - yB;

        (vars(1) - vars(3))^2 + (vars(2) - vars(4))^2 - vars(5)^2;     % eq circonferenza 1 per T
        (vars(1) - vars(6))^2 + (vars(2) - vars(7))^2 - vars(8)^2;     % eq circonferenza 2 per T
    
        (vars(12) - vars(6))^2 + (vars(13) - vars(7))^2 - vars(8)^2;   % eq circonferenza 2 per T2
        (vars(12) - vars(9))^2 + (vars(13) - vars(10))^2 - vars(11)^2; % eq circonferenza 3 per T2
    
        (vars(1) - vars(3))^2 + (vars(2) - vars(4))^2 - (xA - vars(3))^2 + (yA - vars(4))^2;     % eq circonferenza 1 per A
        (vars(12) - vars(9))^2 + (vars(13) - vars(10))^2 - (xB - vars(9))^2 + (yB - vars(10))^2; % eq circonferenza 3 per B
    
        -(vars(1) - vars(3)) / (vars(2) - vars(4)) + (vars(1) - vars(6)) / (vars(2) - vars(7));    % circonferenza 1 tangente a circonferenza 2
        -(vars(12) - vars(6)) / (vars(13) - vars(7)) + (vars(12) - vars(9)) / (vars(13) - vars(10)); % circonferenza 2 tangente a circonferenza 3
        
        (xA - vars(3))^2 + (yA - vars(4))^2 - vars(5)^2;                % eq circonferenza 1 passa per punto A
        (xB - vars(9))^2 + (yB - vars(10))^2 - vars(11)^2;              % eq circonferenza 3 passa per punto B
    
        abs(sqrt((vars(6) - vars(3))^2 + (vars(7) - vars(4))^2) - abs(vars(8) - vars(5)));   % distanza tra i centri = differenza tra i raggi 1-2
        abs(sqrt((vars(9) - vars(6))^2 + (vars(10) - vars(7))^2) - abs(vars(11) - vars(8))); % distanza tra i centri = differenza tra i raggi 2-3       
    
        (atan2(vars(2) - vars(4), vars(1) - vars(3))    + ...
         atan2(vars(2) - vars(7), vars(1) - vars(6))    + ...
         atan2(vars(13) - vars(10), vars(12) - vars(9)))/2 < pi;
    
        vars(3) > xA; vars(3) < xB;
        vars(6) < xB;
        vars(9) < xB;
        vars(1) > xA; vars(1) < xB;
        vars(12)> xA;vars(12) < xB;
        vars(2) < yA; vars(2) > yB;
        vars(13) > yB;
        vars(12) > vars(1);    
        vars(13) < vars(2);  
        ];
    
    % Definizione delle opzioni di fsolve
    % options = optimoptions('fsolve', 'TolX', 1e-4, 'TolFun', 1e-4,'Display','iter','Algorithm','levenberg-marquardt');
    options = optimoptions('fsolve', 'TolX', 1e-3, 'TolFun', 1e-3,'Algorithm','levenberg-marquardt');

    % Risoluzione del sistema di equazioni
    xT = xA + (xB/2 - xA) * rand; yT = yA*1.5 + (yA - yA*1.5) * rand;
    xT2 = xT + (xB - xT) * rand; yT2 = yT + (yB - yT) * rand;
    x1 = xA + (xB - xA) * rand; y1 = yA + (yB - yA) * rand;
    r1 = xA * rand;
    x2 = xA + (xB - xA) * rand; y2 = yA + (yB - yA) * rand;
    r2 = (xB-xA) * rand;
    x3 = xA + (xB - xA) * rand; y3 = yA + (yB - yA) * rand;
    r3 = xB * rand;
    
    % Generazione casuale del punto di inizio
    x0 = [xT, yT, x1, y1, r1, x2, y2, r2, x3, y3, r3, xT2, yT2];     
    sol = fsolve(fun, x0, options);
    
    % Estrazione delle soluzioni
    solxT = sol(1);
    solyT = sol(2);
    solx1 = sol(3);
    soly1 = sol(4);
    solr1 = sol(5);
    solx2 = sol(6);
    soly2 = sol(7);
    solr2 = sol(8);
    solx3 = sol(9);
    soly3 = sol(10);
    solr3 = sol(11);
    solxT2 = sol(12);
    solyT2 = sol(13);
    
    % Centri delle circonferenze e i loro raggi
    center1 = [solx1, soly1];
    center2 = [solx2, soly2];
    center3 = [solx3, soly3];
    radius1 = solr1;
    radius2 = solr2;
    radius3 = solr3;
    
    % Punto di tangenza
    tangent_point  = [solxT, solyT];
    tangent_point2 = [solxT2, solyT2];
    
    % Calcolare gli angoli degli archi
    arc_angle_1 = atan2(solyT - soly1, solxT - solx1);
    arc_angle_2 = atan2(solyT - soly2, solxT - solx2);
    arc_angle_3 = atan2(solyT2 - soly3, solxT2 - solx3);
    
    % Angolo iniziale e finale per tracciare gli archi    
    theta1 = linspace(arc_angle_1,pi,N/3);
    theta2 = linspace(arc_angle_3,arc_angle_2,N/3);
    theta3 = linspace(0,arc_angle_3,N/3);

    % Coordinate degli archi delle circonferenze
    arc1_x = center1(1) + radius1 * cos(theta1);
    arc1_y = center1(2) + radius1 * sin(theta1);
    arc2_x = center2(1) + radius2 * cos(theta2);
    arc2_y = center2(2) + radius2 * sin(theta2);
    arc3_x = center3(1) + radius3 * cos(theta3);
    arc3_y = center3(2) + radius3 * sin(theta3);
    
    ARCS = [horzcat(flip(arc1_x),flip(arc2_x),flip(arc3_x));horzcat(flip(arc1_y),flip(arc2_y),flip(arc3_y))];

    % Calculate the distance between corresponding points
    distances = sqrt((r'-ARCS(1,:)).^2 + (z'-ARCS(2,:)).^2); 
    
    % Compute the mean distance
    max_distance = max(distances);
end

%% TF shape
ARCS = [horzcat(flip(arc1_x),flip(arc2_x),flip(arc3_x));horzcat(flip(arc1_y),flip(arc2_y),flip(arc3_y))]';
z_il = linspace(-ARCS(1,2),ARCS(1,2),N/3);
r_il = zeros(1,N/3)+min(r);
z_ARCS = [ARCS(:,2); -flip(ARCS(:,2)); z_il'];
r_ARCS = [ARCS(:,1); flip(ARCS(:,1)); r_il'];
% plot(r_ARCS,z_ARCS,'.')

C1r = center1(1);
C1z = center1(2);
C2r = center2(1);
C2z = center2(2);
C3r = center3(1);
C3z = center3(2);
theta1 = 180-arc_angle_1*180/pi;
theta3 = arc_angle_3*180/pi;
theta2 = arc_angle_2*180/pi - theta3;
tot_theta = theta1 + theta2 + theta3;
P_tan_1_2 = [arc1_x(1),arc1_y(1)];
P_tan_2_3 = [arc2_x(1),arc2_y(1)];

% plot TF shape
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
% Traccia gli archi delle circonferenze
plot(arc1_x, arc1_y, 'b'); hold on; axis equal
plot(arc2_x, arc2_y, 'r');
plot(arc3_x, arc3_y, 'g');

% Traccia i raggi verso il punto di tangenza
plot([center1(1), tangent_point(1)], [center1(2), tangent_point(2)], 'b--');
plot([center1(1), xA], [center1(2), yA], 'b--');
plot([center2(1), tangent_point(1)], [center2(2), tangent_point(2)], 'r--');
plot([center2(1), tangent_point2(1)], [center2(2), tangent_point2(2)], 'r--');
plot([center3(1), tangent_point2(1)], [center3(2), tangent_point2(2)], 'g--');
plot([center3(1), xB], [center3(2), yB], 'g--');

% Traccia i centri delle circonferenze
plot(center1(1), center1(2), 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
plot(center2(1), center2(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
plot(center3(1), center3(2), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');

% Traccia il punto di tangenza
plot(tangent_point(1), tangent_point(2), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
plot(tangent_point2(1), tangent_point2(2), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');

% % Traccia la retta verticale x = xA
% line([xA, xA], ylim, 'Color', 'k', 'LineStyle', '--');
% % Traccia la retta verticale x = xB
% line([xB, xB], ylim, 'Color', 'k', 'LineStyle', '--');

% Traccia i punti estremi della D
plot(xA, yA, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
plot(xB, yB, 'ks', 'MarkerSize', 8, 'MarkerFaceColor', 'k');

% Plot BF
plot(r,z, 'Color', 'k', 'LineStyle', '--');

% Etichette e titolo del grafico
xlabel('x');
ylabel('y');

title(TitleName0);

legend('Circonferenza 1', 'Circonferenza 2', 'Circonferenza 3');

axis equal;
% xlim([0 xB+xA])

if printcond == 1, print('-dpng','-r300',TitleName0), end

%% curves comparison
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
hold on
plot(rs,zs,'.')
plot(r_ARCS,z_ARCS,'.')
plot(r_BF,z_BF,'.')

legend('Optimized BF shape', 'Tangetn arcs', 'BF shape');

ylabel('z [m]','fontsize',20,'Interpreter','latex')
xlabel('r [m]','fontsize',20,'Interpreter','latex')
set(gca,'FontSize',20);
axis equal
grid on
hold off

if printcond == 1, print('-dpng','-r300',TitleName2), end
%%
% Open a text file to write
fileID = fopen('Shape Details.txt', 'w');

% Write outputs to the text file
fprintf(fileID, 'C1:\n');
fprintf(fileID, '%f %f\n', C1r, C1z);
fprintf(fileID, 'R1:\n');
fprintf(fileID, '%f\n', radius1);
fprintf(fileID, 'theta1:\n');
fprintf(fileID, '%f\n', theta1);

fprintf(fileID, 'C2:\n');
fprintf(fileID, '%f %f\n', C2r, C2z);
fprintf(fileID, 'R2:\n');
fprintf(fileID, '%f\n', radius2);
fprintf(fileID, 'theta2:\n');
fprintf(fileID, '%f\n', theta2);

fprintf(fileID, 'C3:\n');
fprintf(fileID, '%f %f\n', C3r, C3z);
fprintf(fileID, 'R3:\n');
fprintf(fileID, '%f\n', radius3);
fprintf(fileID, 'theta3:\n');
fprintf(fileID, '%f\n', theta3);

fprintf(fileID, 'P 1-2:\n');
fprintf(fileID, '%f %f\n', arc1_x(1), arc1_y(1));
fprintf(fileID, 'P 1-2:\n');
fprintf(fileID, '%f %f\n', arc2_x(1), arc2_y(1));

% Close the file
fclose(fileID);

%%
theta1 = linspace(arc_angle_1,pi,N/3);
theta2 = linspace(arc_angle_3,arc_angle_2,N/3);
theta3 = linspace(0,arc_angle_3,N/3);

% curva interna (lato plasma) del case
r_ARCS_casee = r_ARCS;
z_ARCS_casee = z_ARCS+dz_0;

% center line WP
arc1_x_WPi = center1(1) + (radius1+dr_plasma_side+WP_h/2) * cos(theta1);
arc1_y_WPi = center1(2) + (radius1+dr_plasma_side+WP_h/2) * sin(theta1);
arc2_x_WPi = center2(1) + (radius2+dr_plasma_side+WP_h/2) * cos(theta2);
arc2_y_WPi = center2(2) + (radius2+dr_plasma_side+WP_h/2) * sin(theta2);
arc3_x_WPi = center3(1) + (radius3+dr_plasma_side+WP_h/2) * cos(theta3);
arc3_y_WPi = center3(2) + (radius3+dr_plasma_side+WP_h/2) * sin(theta3);
ARCS_WPi = [horzcat(flip(arc1_x_WPi),flip(arc2_x_WPi),flip(arc3_x_WPi));horzcat(flip(arc1_y_WPi),flip(arc2_y_WPi),flip(arc3_y_WPi))];
z_il = linspace(-ARCS_WPi(2,1),ARCS_WPi(2,1),N/3);
r_il = zeros(1,N/3)+min(ARCS_WPi(1,:));
z_ARCS_CL = [ARCS_WPi(2,:)+dz_0, -flip(ARCS_WPi(2,:))+dz_0, z_il+dz_0]';
r_ARCS_CL = [ARCS_WPi(1,:), flip(ARCS_WPi(1,:)), r_il]';

% curva interna (lato plasma) WP
arc1_x_WPi = center1(1) + (radius1+dr_plasma_side) * cos(theta1);
arc1_y_WPi = center1(2) + (radius1+dr_plasma_side) * sin(theta1);
arc2_x_WPi = center2(1) + (radius2+dr_plasma_side) * cos(theta2);
arc2_y_WPi = center2(2) + (radius2+dr_plasma_side) * sin(theta2);
arc3_x_WPi = center3(1) + (radius3+dr_plasma_side) * cos(theta3);
arc3_y_WPi = center3(2) + (radius3+dr_plasma_side) * sin(theta3);
ARCS_WPi = [horzcat(flip(arc1_x_WPi),flip(arc2_x_WPi),flip(arc3_x_WPi));horzcat(flip(arc1_y_WPi),flip(arc2_y_WPi),flip(arc3_y_WPi))];
z_il = linspace(-ARCS_WPi(2,1),ARCS_WPi(2,1),N/3);
r_il = zeros(1,N/3)+min(ARCS_WPi(1,:));
z_ARCS_WPi = [ARCS_WPi(2,:)+dz_0, -flip(ARCS_WPi(2,:))+dz_0, z_il+dz_0]';
r_ARCS_WPi = [ARCS_WPi(1,:), flip(ARCS_WPi(1,:)), r_il]';

arc1_x_WPe = center1(1) + (radius1+dr_plasma_side+WP_h) * cos(theta1);
arc1_y_WPe = center1(2) + (radius1+dr_plasma_side+WP_h) * sin(theta1);
arc2_x_WPe = center2(1) + (radius2+dr_plasma_side+WP_h) * cos(theta2);
arc2_y_WPe = center2(2) + (radius2+dr_plasma_side+WP_h) * sin(theta2);
arc3_x_WPe = center3(1) + (radius3+dr_plasma_side+WP_h) * cos(theta3);
arc3_y_WPe = center3(2) + (radius3+dr_plasma_side+WP_h) * sin(theta3);
ARCS_WPe = [horzcat(flip(arc1_x_WPe),flip(arc2_x_WPe),flip(arc3_x_WPe));horzcat(flip(arc1_y_WPe),flip(arc2_y_WPe),flip(arc3_y_WPe))];
z_il = linspace(-ARCS_WPe(2,1),ARCS_WPe(2,1),N/3);
r_il = zeros(1,N/3)+min(ARCS_WPe(1,:));
z_ARCS_WPe = [ARCS_WPe(2,:)+dz_0, -flip(ARCS_WPe(2,:))+dz_0, z_il+dz_0]';
r_ARCS_WPe = [ARCS_WPe(1,:), flip(ARCS_WPe(1,:)), r_il]';

arc1_x_casei = center1(1) + (radius1+dr_plasma_side+WP_h+nose_il) * cos(theta1);
arc1_y_casei = center1(2) + (radius1+dr_plasma_side+WP_h+nose_il) * sin(theta1);
arc2_x_casei = center2(1) + (radius2+dr_plasma_side+WP_h+nose_ol) * cos(theta2);
arc2_y_casei = center2(2) + (radius2+dr_plasma_side+WP_h+nose_ol) * sin(theta2);
arc3_x_casei = center3(1) + (radius3+dr_plasma_side+WP_h+nose_ol) * cos(theta3);
arc3_y_casei = center3(2) + (radius3+dr_plasma_side+WP_h+nose_ol) * sin(theta3);
ARCS_casei = [horzcat(flip(arc1_x_casei),flip(arc2_x_casei),flip(arc3_x_casei));horzcat(flip(arc1_y_casei),flip(arc2_y_casei),flip(arc3_y_casei))];
z_il = linspace(-ARCS_casei(2,1),ARCS_casei(2,1),N/3);
r_il = zeros(1,N/3)+min(ARCS_casei(1,:));
z_ARCS_casei = [ARCS_casei(2,:)+dz_0, -flip(ARCS_casei(2,:))+dz_0, z_il+dz_0]';
r_ARCS_casei = [ARCS_casei(1,:), flip(ARCS_casei(1,:)), r_il]';

%
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9])
hold on
plot(r_ARCS_CL,z_ARCS_CL,'b--')
plot(r_ARCS_casee,z_ARCS_casee,'k')
plot(r_ARCS_WPe,z_ARCS_WPe,'r')
plot(r_ARCS_WPi,z_ARCS_WPi,'r')
plot(r_ARCS_casei,z_ARCS_casei,'k')

% Open file for writing
fileID = fopen(name_file, 'w');
% Define scaling factor
mm = 1e3;
% Add column headers
fprintf(fileID, 'r_case_CL z_case_CL r_case_i z_case_i r_case_e z_case_e r_WP_i z_WP_i r_WP_e z_WPe\n');
% Convert data to mm and ensure column-major order
data = [r_ARCS_CL, z_ARCS_CL, r_ARCS_casei, z_ARCS_casei, r_ARCS_casee, z_ARCS_casee, r_ARCS_WPi, z_ARCS_WPi, r_ARCS_WPe, z_ARCS_WPe] * mm;
% Write data row by row
fprintf(fileID, '%f %f %f %f %f %f %f %f %f %f\n', data.');
% Close file
fclose(fileID);

ylabel('z [m]','fontsize',20,'Interpreter','latex')
xlabel('r [m]','fontsize',20,'Interpreter','latex')
title(TitleName,'fontsize',20,'Interpreter','latex')
set(gca,'FontSize',20);
axis equal
grid minor
hold off

% plasma = readmatrix('plasma_shapes.txt');
% hold on
% plot(plasma(:,1)./1,plasma(:,2)./1,'m')
% 
% vessel = readmatrix('vessel_shape.txt');
% hold on
% plot(vessel(:,1)./1,vessel(:,2)./1,'sg')

if printcond == 1, print('-dpng','-r300',TitleName), end

 

%%

a =  readmatrix('TF shape - VNS 07 2026 - 09-Jul-2026.txt');

figure
hold on
plot(a(:,1),a(:,2))
plot(a(:,3),a(:,4))
plot(a(:,5),a(:,6))
plot(a(:,7),a(:,8))
plot(a(:,9),a(:,10))
axis equal
legend('1','2','3','4','5')

max(a(:,2))
max(a(:,4))
max(a(:,6))
max(a(:,8))
max(a(:,10))