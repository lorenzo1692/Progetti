% TF coils: keep user's D-shape block, 3D placement, Biot–Savart along coils
clearvars; clc; close all;

%% ------------------------ Physical constants
mu0 = 4*pi*1e-7;      % Vacuum permeability [T·m/A]

%% ------------------------ Tokamak / TF configuration (VNS-like)
B0      = 5.6;        % Field on magnetic axis [T]
R0      = 2.67;       % Major radius [m]
Ar      = 4.25;       % Aspect ratio
n_TF    = 12;         % Number of TF coils
d       = 0.86;       % Screen + VV gap [m]
rho     = 0.01;       % Plasma ripple
a       = R0/Ar;      % Minor radius [m]
Iop     = 50e3;       % Operating current per conductor [A]
dx_WP   = 0.46;
dy_WP   = 0.33;
NL = 12; 
NT = 9;
area_cross_section = dx_WP * dy_WP;  % Conductor pack area [m^2] (informative)
area_conductor = area_cross_section/NT/NL;
Ri  = R0 - a - d;                  % Inner radius of the D coil [m]
Re  = (R0 + a) * (1/rho)^(1/n_TF); % Outer radius of the D coil [m]
Ric = Ri-dy_WP/2; % WP centroid innerleg of the D coil [m]
Rec = Re+dy_WP/2; % WP centroid outerleg of the D coil [m]

% Coil current per TF coil from an on-axis estimate
TF_current   = (2*pi*R0*B0/mu0) / n_TF;   % [A]
B_TF = mu0*n_TF*TF_current/(2*pi*Ri);
n_conductor  = TF_current / Iop;          % conductors per coil (info)
J            = TF_current / area_cross_section; % current density [A/m^2] (info)
const = mu0 * TF_current / (4*pi);

% --- Hurwitz regularization parameter ---
a_reg = sqrt(area_conductor/pi); % [m] characteristic smoothing length
a2_over_sqrte = (a_reg^2) / sqrt(exp(1));

%% ------------------------ Discretization
N = 50;                    % Number of points for each half (curve/straight)
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
z_straight = linspace(z_curve(1), z_curve(end), N/2);
r_straight = z_straight * 0 + r_curve(1);

% Combine the curved and straight sections
r = [r_curve, fliplr(r_straight)]'; 
z = [z_curve, fliplr(z_straight)]';
BF = [r, z];   % [2N x 2], poloidal polyline (r,z)

% Quick 2D check (poloidal r–z)
figure('Name','Poloidal D-shape (r–z)');
plot(BF(:,1), BF(:,2), 'k.', 'LineWidth', 1.5); axis equal; grid on; hold on
xlabel('r [m]'); ylabel('z [m]'); hold off
title('Poloidal D-shape (kept from user block)');

%% ------------------------ Poloidal segments and tangents
dr   = (BF(:,1)-BF([2:end,1],1));
dz   = (BF(:,2)-BF([2:end,1],2));
dl2D = max(hypot(dr, dz), 1e-12);   % segment lengths, avoid zeros
t_r  = dr ./ dl2D;                  % unit tangent in r
t_z  = dz ./ dl2D;                  % unit tangent in z

% Segment midpoints (poloidal)
r_mid = (BF(:,1)+BF([2:end,1],1))*0.5; % 0.5*(BF(1:end-1,1) + BF(2:end,1));
z_mid = (BF(:,2)+BF([2:end,1],2))*0.5; % 0.5*(BF(1:end-1,2) + BF(2:end,2));
% figure; plot(r_mid,z_mid,'.'); hold on; plot(BF(:,1), BF(:,2),'s'); axis equal; grid on;

% Counts
Npts = size(BF,1);         % N points
Nsegm   = Npts;        % segments (and midpoints) per coil

%% ------------------------ Build 3D-line geometry for each coil (fixed toroidal φ per coil)
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

% ------------------------ 3D visualization: coil polylines + current direction
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

%% ===================== SORGENTI: 12 TF ciascuna con pacchetto NL×NT =====================
% Griglia nella sezione del coil: u = normale poloidale (spessore dx_WP), v = toroidale (dy_WP)

u_list = linspace(-dx_WP/2, +dx_WP/2, NL);    % [m] normale poloidale
v_list = linspace(-dy_WP/2, +dy_WP/2, NT);    % [m] toroidale (spessore)
[UU,VV] = meshgrid(u_list, v_list);
u_vec = UU(:);                    % [Ncond x 1]
v_vec = VV(:);                    % [Ncond x 1]
Ncond = numel(u_vec);

% Normale poloidale lungo la polilinea (r–z)
n_r = -t_z(:);                    % [Nsegm x 1]
n_z =  t_r(:);                    % [Nsegm x 1]

% Centerline poloidale come righe
r0 = r_mid(:).';                  % [1 x Nsegm]
z0 = z_mid(:).';                  % [1 x Nsegm]

% Offsets u (uguali per tutte le TF)
r_u_cond = bsxfun(@plus, r0, u_vec * (n_r.'));   % [Ncond x Nsegm]
z_u_cond = bsxfun(@plus, z0, u_vec * (n_z.'));   % [Ncond x Nsegm]

% Prealloc sorgenti (midpoints e dl)   [n_TF x Ncond x Nsegm]
x_mid_src  = zeros(n_TF, Ncond, Nsegm);
y_mid_src  = zeros(n_TF, Ncond, Nsegm);
z_mid_src  = zeros(n_TF, Ncond, Nsegm);
dlx_src    = zeros(n_TF, Ncond, Nsegm);
dly_src    = zeros(n_TF, Ncond, Nsegm);
dlz_src    = zeros(n_TF, Ncond, Nsegm);

% ------------------------ 3D visualization: coil polylines + current direction
figure('Name','TF coils (3D) and current directions'); hold on; axis equal; grid on; view(3);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('D-shaped TF coils and current directions');

for ii = 1:n_TF
    phi0 = phi_all(ii);
    % phi' depends on v and local r
    phi_u = phi0 + bsxfun(@rdivide, v_vec, r_u_cond);  % [Ncond x Nsegm]

    % midpoints 3D
    x_mid_src(ii,:,:) = r_u_cond .* cos(phi_u);
    y_mid_src(ii,:,:) = r_u_cond .* sin(phi_u);
    z_mid_src(ii,:,:) = z_u_cond;

    % dl 3D (parallelo alla centerline, ruotato a phi')
    dlx_src(ii,:,:) = bsxfun(@times, dr.', cos(phi_u));
    dly_src(ii,:,:) = bsxfun(@times, dr.', sin(phi_u));
    dlz_src(ii,:,:) = repmat(dz.', Ncond, 1);

    plot3(squeeze(x_mid_src(ii,:,:)), squeeze(y_mid_src(ii,:,:)), squeeze(z_mid_src(ii,:,:)), 'k.', 'LineWidth', 1.0);
end
hold off;

%% ===================== PUNTI DI VALUTAZIONE: coil #1 (theta=0), NL×NT completi ==========
phi_eval = 0;
% costruiamo i midpoints dei filamenti su coil #1 (stesso reticolo u,v)
phi_u_eval = phi_eval + bsxfun(@rdivide, v_vec, r_u_cond);   % [Ncond x Nsegm]
Xe_mat = r_u_cond .* cos(phi_u_eval);    % [Ncond x Nsegm]
Ye_mat = r_u_cond .* sin(phi_u_eval);
Ze_mat = z_u_cond;

Xe = Xe_mat(:); Ye = Ye_mat(:); Ze = Ze_mat(:);
Ne = numel(Xe);

u_list = linspace(-dx_WP/2, +dx_WP/2, NL);
v_list = linspace(-dy_WP/2, +dy_WP/2, NT);
[UU,VV] = meshgrid(u_list, v_list);
u_vec = UU(:); v_vec = VV(:);
n_r = -t_z(:); n_z = t_r(:);
r0 = r_mid(:).'; z0 = z_mid(:).';
r_u_cond = bsxfun(@plus, r0, u_vec * (n_r.'));  % [Ncond x Nsegm]
z_u_cond = bsxfun(@plus, z0, u_vec * (n_z.'));  % [Ncond x Nsegm]

% --- Costruisci i PUNTI DI VALUTAZIONE per tutte le TF (rotando φ) ---
X_eval_all = zeros(n_TF, Ncond, Nsegm);
Y_eval_all = zeros(n_TF, Ncond, Nsegm);
Z_eval_all = zeros(n_TF, Ncond, Nsegm);
for ii = 1:n_TF
    phi_eval = phi_all(ii);
    phi_u_eval = phi_eval + bsxfun(@rdivide, v_vec, r_u_cond);  % [Ncond x Nsegm]
    X_eval_all(ii,:,:) = r_u_cond .* cos(phi_u_eval);
    Y_eval_all(ii,:,:) = r_u_cond .* sin(phi_u_eval);
    Z_eval_all(ii,:,:) = z_u_cond;
end

% --- Appiattisci EVAL e SORGENTI ---
Xe_all = X_eval_all(:);  Ye_all = Y_eval_all(:);  Ze_all = Z_eval_all(:);
Ne_all = numel(Xe_all);

Xs = x_mid_src(:); Ys = y_mid_src(:); Zs = z_mid_src(:);
DLx = dlx_src(:);  DLy = dly_src(:);  DLz = dlz_src(:);
nSegTotal = numel(Xs);

% Corrente per filamento (corrente totale per coil mantenuta)
I_per_fil = TF_current / (NL*NT);
const_fil = mu0 * I_per_fil / (4*pi);

%% ------------------------ Biot–Savart: B at coil midpoints (Hurwitz regularized)
% B(P) = μ0/(4π) * I * Σ [ dl × R ] / (|R|^2 + a^2/√e)^(3/2)
% Evaluates the regularized field at each segment midpoint of each coil.

Bx_all = zeros(Ne_all,1);
By_all = zeros(Ne_all,1);
Bz_all = zeros(Ne_all,1);

K = 9; frac = 0.5; bytes = 8;   % stima memoria temporanei
try
    m = memory; M_free = frac * m.MaxPossibleArrayBytes;
catch
    M_free = 1.0 * 1024^3;  % fallback 1 GiB
end
target = max(1, floor(M_free/(bytes*K)));
eval_chunk = min(Ne_all, max(256, floor(sqrt(target))));
src_block  = min(nSegTotal, max(256, floor(target / eval_chunk)));
eval_chunk = eval_chunk - mod(eval_chunk,64); eval_chunk = max(eval_chunk,256);
src_block  = src_block  - mod(src_block,64);  src_block  = max(src_block,256);

for eStart = 1:eval_chunk:Ne_all
    eStop = min(eStart+eval_chunk-1, Ne_all);
    Xe_blk = Xe_all(eStart:eStop);
    Ye_blk = Ye_all(eStart:eStop);
    Ze_blk = Ze_all(eStart:eStop);

    for sStart = 1:src_block:nSegTotal
        sStop = min(sStart+src_block-1, nSegTotal);
        Xsb = Xs(sStart:sStop).';   Ysb = Ys(sStart:sStop).';   Zsb = Zs(sStart:sStop).';
        DLxb= DLx(sStart:sStop).';  DLyb= DLy(sStart:sStop).';  DLzb= DLz(sStart:sStop).';

        Rx = Xe_blk - Xsb;   Ry = Ye_blk - Ysb;   Rz = Ze_blk - Zsb;
        R2 = Rx.^2 + Ry.^2 + Rz.^2;
        R3 = (R2 + a2_over_sqrte).^(3/2);

        cx =  DLyb.*Rz - DLzb.*Ry;
        cy =  DLzb.*Rx - DLxb.*Rz;
        cz =  DLxb.*Ry - DLyb.*Rx;

        Bx_all(eStart:eStop) = Bx_all(eStart:eStop) + const_fil * sum(cx ./ R3, 2);
        By_all(eStart:eStop) = By_all(eStart:eStop) + const_fil * sum(cy ./ R3, 2);
        Bz_all(eStart:eStop) = Bz_all(eStart:eStop) + const_fil * sum(cz ./ R3, 2);
    end
end

% --- Rimetti in shape [n_TF x NT x NL x Nsegm] per plotting comodo ---
Bx_full = reshape(Bx_all, [n_TF, Ncond, Nsegm]);   % Ncond = NT*NL
By_full = reshape(By_all, [n_TF, Ncond, Nsegm]);
Bz_full = reshape(Bz_all, [n_TF, Ncond, Nsegm]);

% ricostruisci coordinate corrispondenti
X_full  = reshape(Xe_all, [n_TF, Ncond, Nsegm]);
Y_full  = reshape(Ye_all, [n_TF, Ncond, Nsegm]);
Z_full  = reshape(Ze_all, [n_TF, Ncond, Nsegm]);

% risagoma in [n_TF x NT x NL x Nsegm] (ordine come meshgrid v,u)
Bx_full = reshape(Bx_full, [n_TF, NT, NL, Nsegm]);
By_full = reshape(By_full, [n_TF, NT, NL, Nsegm]);
Bz_full = reshape(Bz_full, [n_TF, NT, NL, Nsegm]);
X_full  = reshape(X_full , [n_TF, NT, NL, Nsegm]);
Y_full  = reshape(Y_full , [n_TF, NT, NL, Nsegm]);
Z_full  = reshape(Z_full , [n_TF, NT, NL, Nsegm]);

Bmag_full = sqrt(Bx_full.^2 + By_full.^2 + Bz_full.^2);

%% ===================== PLOT: campo su tutti i punti di tutte le TF (quiver) =================
% Decimazione per resa grafica (metti =1 per TUTTI i punti)
decimCoil = 1;                   % seleziona 1 ogni 'decimCoil' coils
decimV    = 1;                   % tra le file toroidali NT
decimU    = 1;                   % tra gli strati NL
decimS    = 1;                   % tra i segmenti lungo la polilinea (aumenta se fitto)

ii_coils = 1:decimCoil:n_TF;
iv = 1:decimV:NT;  iu = 1:decimU:NL;  is = 1:decimS:Nsegm;

% prepara vettori per quiver/scatter
xq = X_full(ii_coils,iv,iu,is); xq = xq(:);
yq = Y_full(ii_coils,iv,iu,is); yq = yq(:);
zq = Z_full(ii_coils,iv,iu,is); zq = zq(:);

uq = Bx_full(ii_coils,iv,iu,is); uq = uq(:);
vq = By_full(ii_coils,iv,iu,is); vq = vq(:);
wq = Bz_full(ii_coils,iv,iu,is); wq = wq(:);

Bmag_q = sqrt(uq.^2 + vq.^2 + wq.^2);

%
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

%
figure('Name','Campo su tutti i punti di tutte le TF'); clf; hold on; grid on; axis equal vis3d;
view(35,20);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('Vettori B su %d TF (NL=%d, NT=%d) — sorgenti: %d TF × %d cond.', ...
    numel(ii_coils), NL, NT, n_TF, NL*NT));
sc = scatter3(xq, yq, zq, 8, Bmag_q, 'filled'); 
colormap(jet); c = colorbar; c.Label.String = '|B| [T]';
for ii = 1:n_TF
    phi = phi_all(ii); cph = cos(phi); sph = sin(phi);
    plot3(BF(:,1).*cph, BF(:,1).*sph, BF(:,2), 'k-', 'LineWidth', 0.6);
end
hold off;

%% Mappa |B| sulla sezione [4D: n_TF, NT, NL, Nsegm]
figure('Name','|B| sezione del coil TF #1 (piano locale), segmento ~metà','units','normalized','outerposition',[0.25 0.1 0.5 0.9]); 
plot(BF(:,1), BF(:,2), 'k-', 'LineWidth', 1.0); axis equal; grid on; hold on;
hold on

for sez = 20
    B_slice = Bmag_full(1,:,:,sez); B_slice = reshape(B_slice, [NT, NL]); 
    xq = X_full(1,:,:,sez); xq = reshape(xq, [NT, NL]); % n_TF, NT, NL, Nsegm
    yq = Y_full(1,:,:,sez); yq = reshape(yq, [NT, NL]);
    zq = Z_full(1,:,:,sez); zq = reshape(zq, [NT, NL]);
    surf(xq, zq, yq, B_slice);
end
% view(2)
set(gca,'FontSize',20);
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'B [T]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
MAXMIN=get(c,'Limits');
T = linspace(MAXMIN(1),MAXMIN(2),9);
set(c,'Ticks',T)

xlabel('[m]','fontsize',20,'Interpreter','latex')
ylabel('[m]','fontsize',20,'Interpreter','latex')
zlabel('[m]','fontsize',20,'Interpreter','latex')
title('|B| sulla sezione NL×NT alla metà della polilinea (TF #1)');

%% BUILD_INDUCTANCE_MATRIX

% [n_TF x Ncond x Nsegm]

n=0;
for i = 1:n_TF
    for j = 1: Ncond
        for k = 1: Nsegm
            n=n+1;
            % midpoints 3D
            Xe_all(n,1) = x_mid_src(i,j,k); 
            Ye_all(n,1) = y_mid_src(i,j,k); 
            Ze_all(n,1) = z_mid_src(i,j,k); 
            
            % dl 3D (parallelo alla centerline, ruotato a phi')
            DLx(n,1) = dlx_src(i,j,k); 
            DLy(n,1) = dly_src(i,j,k); 
            DLz(n,1) = dlz_src(i,j,k); 
        end 
    end 
end 

a = x_mid_src(1,1,:);
b = y_mid_src(1,1,:);
c = z_mid_src(1,1,:);

figure; plot3(Xe_all(1:Nsegm),Ye_all(1:Nsegm),Ze_all(1:Nsegm),'.'); axis equal

figure; plot3(a(:),b(:),c(:),'.'); axis equal

I_spire = Iop*ones(n_TF,1);
const_A = mu0/(4*pi) * 1;

nSegPerCoil = NT * NL * Nsegm;

% Build segment index lists if not provided (assume coil-major packing)
segIndexCoil = cell(n_TF,1);
for q = 1:n_TF
    i0 = (q-1)*nSegPerCoil + 1;
    i1 = q*nSegPerCoil;
    segIndexCoil{q} = i0:i1;
end

L = zeros(n_TF, n_TF);

% === Loop over source coils q: build A_q everywhere, then fill column q ===
for q = 1:1 %n_TF
    idx_source = segIndexCoil{q}(:).';

    % Pre-extract sources
    Xs_q  = Xe_all(idx_source).';  Ys_q  = Ye_all(idx_source).';  Zs_q  = Ze_all(idx_source).';
    DLx_q = DLx(idx_source).';     DLy_q = DLy(idx_source).';     DLz_q = DLz(idx_source).';

    

    % Da tutti i punti campo sottraggo tutti i punti sorgente
    % Ottengo una matrice

    Rx = Xe_all - Xs_q;   Ry = Ye_all - Ys_q;   Rz = Ze_all - Zs_q;
    R2 = Rx.^2 + Ry.^2 + Rz.^2;
    R  = sqrt(R2 + a2_over_sqrte);   % regularized |R|

    for p = 1:n_TF
        idx_p = segIndexCoil{p};
        
        % Ottengo le componenti di A in tutti i punti campo dovute 
        % a tutti i punti sorgente in parallelo

        Ax_p = const_A * sum(DLx(idx_p) ./ R(idx_p,:), 2) ;
        Ay_p = const_A * sum(DLy(idx_p) ./ R(idx_p,:), 2) ;
        Az_p = const_A * sum(DLz(idx_p) ./ R(idx_p,:), 2) ;
        
        % circuitazione di A (prodotto scalare)
        
        sdot = DLx(idx_p).*Ax_p + DLy(idx_p).*Ay_p + DLz(idx_p).*Az_p;
        L(p,q) = sum(sdot); % [H/A]
    end
end

Lt = L*(N_spire)^2;
E = (1/2*I_spira'*Lt*I_spira)*1e-9; % [GJ]

