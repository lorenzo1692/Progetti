% TF coils: keep user's D-shape block, 3D placement, Biot–Savart along coils
clearvars; clc; close all;

%% ------------------------ Physical constants
mu0 = 4*pi*1e-7;      % Vacuum permeability [T·m/A]

%% ------------------------ Tokamak / TF configuration (VNS-like)
B0      = 4.39;        % Field on magnetic axis [T]
R0      = 8.6;       % Major radius [m]
Ar      = 2.8;        % Aspect ratio
n_TF    = 16;         % Number of TF coils
d       = 1.82;       % Screen + VV gap [m]
rho     = 0.006;       % Plasma ripple
a       = R0/Ar;      % Minor radius [m]

Iop = 50e3;            % Operating current per conductor [A]
dx_WP   = 0.46;
dy_WP   = 0.33;
area_cross_section = dx_WP * dy_WP;  % Conductor pack area [m^2] (informative)
compute_energy  = false;           % Optional O(N^2) energy/inductance

Ri      = R0 - a - d;                  % Inner radius of the D coil [m]
Re      = (R0 + a) * (1/rho)^(1/n_TF); % Outer radius of the D coil [m]

Ric      = Ri-dy_WP/2; % WP centroid innerleg of the D coil [m]
Rec      = Re+dy_WP/2; % WP centroid outerleg of the D coil [m]

% Coil current per TF coil from an on-axis estimate
TF_current   = (2*pi*R0*B0/mu0) / n_TF;   % [A]

B_TF = mu0*n_TF*TF_current/(2*pi*Ri);

n_conductor  = TF_current / Iop;          % conductors per coil (info)
J            = TF_current / area_cross_section; % current density [A/m^2] (info)

const = mu0 * TF_current / (4*pi);

% --- Hurwitz regularization parameter ---
a_reg = dy_WP/2;                        % [m] characteristic smoothing length
a2_over_sqrte = (a_reg^2) / sqrt(exp(1));

%% ------------------------ Discretization
N = 300;                    % Number of points for each half (curve/straight)
% Total points per coil centerline = 2*N, segments = 2*N-1

%% ------------------------ D-shaped profile parameters
k = 0.5 * log(Rec / Ric);
r0 = Ric * exp(k);

% Calculate the curved D profile in cylindrical coordinates
theta1   = linspace(-pi/2, 0, N/4);
theta2   = linspace(0, pi, N);
theta3   = linspace(pi, 3/2*pi, N/4);
theta = [theta1 theta2 theta3];
dtheta  = [0 diff(theta)];
r_curve = r0 .* exp(k .* sin(theta));                                % radial profile
dz_curve= (r0 * k) .* sin(theta) .* exp(k .* sin(theta)) .* dtheta;  % vertical incr
z_curve = cumsum(dz_curve);
dzeta   = (min(z_curve) + max(z_curve)) / 2;
z_curve = z_curve - dzeta;

% Straight segment of the D shape
z_straight = linspace(z_curve(1), z_curve(end), N);
r_straight = z_straight * 0 + r_curve(1);

% Combine the curved and straight sections
r = [r_curve, fliplr(r_straight)]'; 
z = [z_curve, fliplr(z_straight)]';
BF = [r, z];   % [2N x 2], poloidal polyline (r,z)

% Quick 2D check (poloidal r–z)
figure('Name','Poloidal D-shape (r–z)');
plot(BF(:,1), BF(:,2), 'k-', 'LineWidth', 1.5); axis equal; grid on; hold on
xlabel('r [m]'); ylabel('z [m]'); hold off
title('Poloidal D-shape (kept from user block)');

%% ------------------------ Poloidal segments and tangents
dr   = diff(BF(:,1));
dz   = diff(BF(:,2));
dl2D = max(hypot(dr, dz), 1e-12);   % segment lengths, avoid zeros
t_r  = dr ./ dl2D;                  % unit tangent in r
t_z  = dz ./ dl2D;                  % unit tangent in z

% Segment midpoints (poloidal)
r_mid = 0.5*(BF(1:end-1,1) + BF(2:end,1));
z_mid = 0.5*(BF(1:end-1,2) + BF(2:end,2));
% figure; plot(r_mid,z_mid,'.'); hold on; plot(BF(:,1), BF(:,2),'s'); axis equal; grid on;

% Counts
Npts = size(BF,1);         % N points
Nsegm   = Npts - 1;        % segments (and midpoints) per coil

%% ------------------------ Build 3D geometry for each coil (fixed toroidal φ per coil)
phi_all = (0:n_TF-1).' * (2*pi/n_TF);  % [n_TF x 1]

% Preallocate midpoints and segment vectors in 3D
x_mid  = zeros(n_TF, Nsegm);
y_mid  = zeros(n_TF, Nsegm);
z_mid3 = repmat(z_mid.', n_TF, 1);

dlx = zeros(n_TF, Nsegm);
dly = zeros(n_TF, Nsegm);
dlz = repmat((t_z .* dl2D).', n_TF, 1);   % z component same for all coils

for ii = 1:n_TF
    phi = phi_all(ii);
    c = cos(phi); s = sin(phi);

    % Midpoints in 3D
    x_mid(ii,:) = r_mid.' .* c;
    y_mid(ii,:) = r_mid.' .* s;

    % Segment vectors in 3D: dl = (dr * e_r) + (dz * e_z),
    % with e_r = (cosφ, sinφ, 0) at that coil
    % % % dlx(ii,:) = (t_r .* dl2D).' * c;
    % % % dly(ii,:) = (t_r .* dl2D).' * s;

    dlx(ii,:) = dr.' * c;
    dly(ii,:) = dr.' * s;
end

%% ------------------------ 3D visualization: coil polylines + current direction
figure('Name','TF coils (3D) and current directions'); hold on; axis equal; grid on; view(3);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('D-shaped TF coils and current directions');

for ii = 1:n_TF
    phi = phi_all(ii); c = cos(phi); s = sin(phi);
    x_pts = BF(:,1) .* c;
    y_pts = BF(:,1) .* s;
    z_pts = BF(:,2);
    % plot3(x_pts, y_pts, z_pts, 'k-', 'LineWidth', 1.0);
    quiver3(x_mid(ii,:), y_mid(ii,:), z_mid3(ii,:), ...
            dlx(ii,:),     dly(ii,:),     dlz(ii,:), 0.5, 'r');
end
hold off;

%% ------------------------ Biot–Savart: B at coil midpoints (Hurwitz regularized)
% B(P) = μ0/(4π) * I * Σ [ dl × R ] / (|R|^2 + a^2/√e)^(3/2)
% Evaluates the regularized field at each segment midpoint of each coil.

% Flatten all coil segments into 1D arrays (length = n_TF * Nsegm)
Xs  = x_mid(:);         % source segment midpoints (X)
Ys  = y_mid(:);         % source segment midpoints (Y)
Zs  = z_mid3(:);        % source segment midpoints (Z)
DLx = dlx(:);           % segment vector components (X)
DLy = dly(:);           % segment vector components (Y)
DLz = dlz(:);           % segment vector components (Z)
nSegTotal = numel(Xs);


%% ===================== Valutazione su TF #1 (theta=0): matrice NL×NT =====================
% Parametri matrice conduttori nella sezione del coil
NL = 12;                         % strati lungo la normale poloidale (spessore dx_WP)
NT = 9;                         % file lungo la direzione toroidale (spessore dy_WP)
u_list = linspace(-dx_WP/2, +dx_WP/2, NL);   % [m] offset normale poloidale
v_list = linspace(-dy_WP/2, +dy_WP/2, NT);   % [m] offset toroidale (spessore)
[UU,VV] = meshgrid(u_list, v_list);
u_vec = UU(:);                          % [Ncond x 1]
v_vec = VV(:);                          % [Ncond x 1]
Ncond = numel(u_vec);

% Normale poloidale lungo la polilinea (r–z)
n_r = -t_z;                 % [Nsegm x 1]
n_z =  t_r;                 % [Nsegm x 1]

% Centerline TF #1 è a phi=0
phi0 = 0;

% Midpoints poloidali come riga
r0 = r_mid.';               % [1 x Nsegm]
z0 = z_mid.';               % [1 x Nsegm]

% Offset u (normale poloidale): r_u, z_u sono [Ncond x Nsegm]
r_u = bsxfun(@plus, r0, u_vec * (n_r.'));      % r + u*n_r
z_u = bsxfun(@plus, z0, u_vec * (n_z.'));      % z + u*n_z

% Offset toroidale come spostamento d’arco: phi' = phi0 + v/r_u
phi_u = bsxfun(@plus, phi0, bsxfun(@rdivide, v_vec, r_u));  % [Ncond x Nsegm]

% Coordinate 3D dei punti di valutazione (midpoints di ogni filamento della TF #1)
Xe = r_u .* cos(phi_u);          % [Ncond x Nsegm]
Ye = r_u .* sin(phi_u);          % [Ncond x Nsegm]
Ze = z_u;                        % [Ncond x Nsegm]

% Appiattisci i punti di valutazione
Xe = Xe(:);  Ye = Ye(:);  Ze = Ze(:);
Ne = numel(Xe);

% ===================== Biot–Savart su questi punti (sorgenti = 12 spire) =================
Bx_eval = zeros(Ne,1);
By_eval = zeros(Ne,1);
Bz_eval = zeros(Ne,1);

% Stima blocco per le SORGENTI (per RAM): dimensione dominante ~ Ne × Nb
K = 9;   % ~n° di grandi temporanei
try
    m = memory;                                  % solo Windows
    M_free = 0.5 * m.MaxPossibleArrayBytes;      % ~50% di margine
    block  = floor( M_free / (8 * K * Ne) );
    block  = max(256, block - mod(block,256));
catch
    block = 4000;                                 % fallback
end
block = min(block, nSegTotal);

for startIdx = 1:block:nSegTotal
    stopIdx = min(startIdx+block-1, nSegTotal);

    % sorgenti del blocco (riga 1×Nb)
    Xsb = Xs(startIdx:stopIdx).';
    Ysb = Ys(startIdx:stopIdx).';
    Zsb = Zs(startIdx:stopIdx).';
    DLxb= DLx(startIdx:stopIdx).';
    DLyb= DLy(startIdx:stopIdx).';
    DLzb= DLz(startIdx:stopIdx).';

    % vettori R = r_eval - r_src  --> [Ne × Nb]
    Rx = Xe - Xsb;
    Ry = Ye - Ysb;
    Rz = Ze - Zsb;

    % distanza regolarizzata
    R2 = Rx.^2 + Ry.^2 + Rz.^2;
    R3 = (R2 + a2_over_sqrte).^(3/2);

    % dl × R
    cx =  DLyb.*Rz - DLzb.*Ry;
    cy =  DLzb.*Rx - DLxb.*Rz;
    cz =  DLxb.*Ry - DLyb.*Rx;

    % accumula contributo
    Bx_eval = Bx_eval + const * sum(cx ./ R3,  2);
    By_eval = By_eval + const * sum(cy ./ R3,  2);
    Bz_eval = Bz_eval + const * sum(cz ./ R3,  2);
end

% Ricomponi nella griglia [NT × NL × Nsegm] (ordine come meshgrid v,u)
Bx_pack = reshape(Bx_eval, [NT, NL, Nsegm]);
By_pack = reshape(By_eval, [NT, NL, Nsegm]);
Bz_pack = reshape(Bz_eval, [NT, NL, Nsegm]);
Bmag_pack = sqrt(Bx_pack.^2 + By_pack.^2 + Bz_pack.^2);

% ===================== Plot rapidi (facoltativi) =====================
% 1) |B| lungo il percorso poloidale per il "filamento centrale" (u=0,v=0)
[~,iu0] = min(abs(u_list));   % indice più vicino a u=0
[~,iv0] = min(abs(v_list));   % indice più vicino a v=0
B_central = squeeze(Bmag_pack(iv0, iu0, :));
s_end = [0; cumsum(dl2D(:))];                 
s_mid = 0.5*(s_end(1:end-1) + s_end(2:end));
figure('Name','|B| lungo la TF #1 (filamento centrale)'); 
plot(s_mid, B_central, 'LineWidth', 1.8);
grid on; xlabel('s poloidale [m]'); ylabel('|B| [T]');
title('|B| lungo TF #1 (u=0, v=0) dovuto a 12 spire sorgenti');

% 3) (opzione) quiver di B su sottoinsieme di filamenti e segmenti
iv = iv0; iu = iu0;    % centrale
% ricostruisci coordinate 3D corrispondenti per quel filamento
r_u1  = r0 + u_list(iu) * (n_r.');
z_u1  = z0 + u_list(iu) * (n_z.');
phi_u1= phi0 + v_list(iv) ./ r_u1;
xq = r_u1 .* cos(phi_u1);  yq = r_u1 .* sin(phi_u1);  zq = z_u1;
uq = squeeze(Bx_pack(iv,iu,:)); vq = squeeze(By_pack(iv,iu,:)); wq = squeeze(Bz_pack(iv,iu,:));
figure('Name','B su TF #1 (filamento centrale)'); hold on; grid on; axis equal; view(3);
plot3(xq, yq, zq, 'k.', 'MarkerSize', 8);
quiver3(xq, yq, zq, uq', vq', wq', 0.6, 'b'); 
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Vettori B su TF #1 (u=0, v=0)'); hold off;

%% =============== QUIVER COMPLETO SUL COIL #1 (theta=0) =================
% Requisiti: BF, r_mid, z_mid, t_r, t_z, NL, NT, dx_WP, dy_WP, n_TF, phi_all
%            Bx_pack, By_pack, Bz_pack  [NT x NL x Nsegm]  già calcolati

% --- ricostruisci le coordinate 3D dei punti (NT x NL x Nsegm)
u_list = linspace(-dx_WP/2, +dx_WP/2, NL);   % offset normale poloidale
v_list = linspace(-dy_WP/2, +dy_WP/2, NT);   % offset toroidale
[UU, VV] = meshgrid(u_list, v_list);         % UU: NT×NL, VV: NT×NL

n_r = -t_z(:);                               % [Nsegm x 1]
n_z =  t_r(:);                               % [Nsegm x 1]
r0  = r_mid(:).';                            % [1 x Nsegm]
z0  = z_mid(:).';                            % [1 x Nsegm]
phi0 = 0;

% offset u lungo la normale poloidale
r_u_all = r0 + reshape(UU, [],1) * (n_r.');  % [Ncond x Nsegm], Ncond=NT*NL
z_u_all = z0 + reshape(UU, [],1) * (n_z.');

% offset v toroidale come spostamento d'arco: phi' = phi0 + v/r
phi_u_all = phi0 + (reshape(VV, [],1) ./ r_u_all);

% coordinate 3D
X_all = r_u_all .* cos(phi_u_all);
Y_all = r_u_all .* sin(phi_u_all);
Z_all = z_u_all;

% rimetti in [NT x NL x Nsegm] (ordine come meshgrid v,u)
Nsegm = size(Bx_pack,3);
X_pack = reshape(X_all, [NT, NL, Nsegm]);
Y_pack = reshape(Y_all, [NT, NL, Nsegm]);
Z_pack = reshape(Z_all, [NT, NL, Nsegm]);

% --- selezione (decimazione) per la resa grafica
decimV   = 1;                     % 1 = TUTTI gli strati toroidali
decimU   = 1;                     % 1 = TUTTI gli strati poloidali (normale)
decimSeg = 1;                     % 1 = TUTTI i segmenti lungo la polilinea
iv = 1:decimV:NT;
iu = 1:decimU:NL;
is = 1:decimSeg:Nsegm;

% estrai sottoinsieme (o tutto) e appiattisci in vettori
xq = X_pack(iv,iu,is);   xq = xq(:);
yq = Y_pack(iv,iu,is);   yq = yq(:);
zq = Z_pack(iv,iu,is);   zq = zq(:);

uq = Bx_pack(iv,iu,is);  uq = uq(:);
vq = By_pack(iv,iu,is);  vq = vq(:);
wq = Bz_pack(iv,iu,is);  wq = wq(:);

Bmag_q = sqrt(uq.^2 + vq.^2 + wq.^2);

%% --- figura quiver completa sul coil #1
figure('Name','Quiver B su coil #1 (theta=0)','units','normalized','outerposition',[0.25 0.1 0.5 0.9]); 
clf; hold on; grid on; axis equal vis3d;
view(35,20);
xlabel('[m]','fontsize',20,'Interpreter','latex')
ylabel('[m]','fontsize',20,'Interpreter','latex')
zlabel('[m]','fontsize',20,'Interpreter','latex')
title(sprintf('Vettori B su coil #1 (NL=%d, NT=%d) — sorgenti: 12 TF coils', NL, NT));
set(gca,'FontSize',20);
scale = 5;                          % fattore di scala delle frecce
q = quiver3(xq, yq, zq, uq, vq, wq, scale, 'b');
[q] = ColorQuiver(q);
q.Color = 'b';
q.LineWidth = 1.1;
q.ShowArrowHead = 'on';
q.MaxHeadSize = 2;
q.AutoScale
axis equal
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'B [T]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([min(Bmag_q) max(Bmag_q)]);
MAXMIN=get(c,'Limits');
T = linspace(MAXMIN(1),MAXMIN(2),9);
set(c,'Ticks',T)

% colorazione del modulo |B| con punti alla base delle frecce
% sc = scatter3(xq, yq, zq, 8, Bmag_q, 'filled'); 
% colormap(jet); c = colorbar; c.Label.String = '|B| [T]';

% sovrapponi le 12 centerline delle spire a D (riferimento)
for ii = 1:n_TF
    phi = phi_all(ii); cph = cos(phi); sph = sin(phi);
    plot3(BF(:,1).*cph, BF(:,1).*sph, BF(:,2), 'k-', 'LineWidth', 0.6);
end
hold off;








