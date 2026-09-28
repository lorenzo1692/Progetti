% TF coils: keep user's D-shape block, 3D placement, Biot–Savart along coils
clearvars; clc; close all;

%% ------------------------ Physical constants
mu0 = 4*pi*1e-7;      % Vacuum permeability [T·m/A]

%% ------------------------ Tokamak / TF configuration (VNS-like)
B0      = 5.6;        % Field on magnetic axis [T]
R0      = 2.67;       % Major radius [m]
Ar      = 4.25;        % Aspect ratio
n_TF    = 12;         % Number of TF coils
d       = 0.86;       % Screen + VV gap [m]
rho     = 0.01;       % Plasma ripple
a       = R0/Ar;      % Minor radius [m]

compute_energy  = false;           % Optional O(N^2) energy/inductance
Ri      = R0 - a - d;                  % Inner radius of the D coil [m]
Re      = (R0 + a) * (1/rho)^(1/n_TF); % Outer radius of the D coil [m]

% Coil current per TF coil from an on-axis estimate
TF_current   = (2*pi*R0*B0/mu0) / n_TF;   % [A]
B_TF = mu0*n_TF*TF_current/(2*pi*Ri)*1.08;
const = mu0 * TF_current / (4*pi);

% Jeng_user = 50*1e6; % [A/m2]
% A_user = TF_current; % [m2]
% dx_WP   = 0.05; % [m]
% dy_WP   = A_user/dx_WP;

dx_WP   = 0.46;
dy_WP   = 0.33;
area_cross_section = dx_WP * dy_WP;  % Conductor pack area [m^2] (informative)
Ric      = Ri-dy_WP/2; % WP centroid innerleg of the D coil [m]
Rec      = Re+dy_WP/2; % WP centroid outerleg of the D coil [m]

Iop = 6e4;
N_spire = ceil(TF_current/Iop);
k_bf = 0.5*log(Rec/Ric); % k bending free
T_bf = 0.5*(k_bf*n_TF*(N_spire*Iop)^2*mu0/(2*pi))*1e-6; % [MN] Hoop tension along TF longitudinal axis
Sv = T_bf/area_cross_section;

% --- Hurwitz regularization parameter ---
a_reg = dy_WP/2;                        % [m] characteristic smoothing length
a2_over_sqrte = (a_reg^2) / sqrt(exp(1));

%% ------------------------ Discretization
N = 300;                    % Number of points for each half (curve/straight)

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

dlxm = zeros(n_TF, Nsegm);
dlym = zeros(n_TF, Nsegm);
dlzm = repmat((t_z .* dl2D).', n_TF, 1);   % z component same for all coils

for ii = 1:n_TF
    phi = phi_all(ii);
    c = cos(phi); s = sin(phi);

    % Midpoints in 3D
    x_mid(ii,:) = r_mid.' .* c;
    y_mid(ii,:) = r_mid.' .* s;

    dlxm(ii,:) = dr.' * c;
    dlym(ii,:) = dr.' * s;
end

% Preallocate source and segment vectors in 3D
x_sour  = zeros(n_TF, Nsegm);
y_sour  = zeros(n_TF, Nsegm);
z_sour3 = repmat(BF(1:end-1,2).', n_TF, 1);

dlxs = zeros(n_TF, Nsegm);
dlys = zeros(n_TF, Nsegm);
dlzs = repmat((t_z .* dl2D).', n_TF, 1);   % z component same for all coils

for ii = 1:n_TF
    phi = phi_all(ii);
    c = cos(phi); s = sin(phi);

    % Midpoints in 3D
    x_sour(ii,:) = BF(1:end-1,1).' .* c;
    y_sour(ii,:) = BF(1:end-1,1).' .* s;

    dlxs(ii,:) = dr.' * c;
    dlys(ii,:) = dr.' * s;
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
            dlxm(ii,:),     dlym(ii,:),     dlzm(ii,:), 0.5, 'r');
end
hold off;
%% ------------------------ Inductance matrix via vector potential (Hurwitz regularized)
% L(p,q) = (1/I_q) * sum_{i in coil p} dl_i · A_q(r_i)
% A_q(r) = mu0/(4*pi) * I_q * sum_{s in coil q} dl_s / sqrt(|r - r_s|^2 + a^2/sqrt(e))

% Flatten all coil segments into 1D arrays (length = n_TF * Nsegm)
Xs  = x_sour(:);         % source segment midpoints (X)
Ys  = y_sour(:);         % source segment midpoints (Y)
Zs  = z_sour3(:);        % source segment midpoints (Z)
DLx = dlxm(:);           % segment vector components (X)
DLy = dlym(:);           % segment vector components (Y)
DLz = dlzm(:);           % segment vector components (Z)
nSegTotal = numel(Xs);

Ne_all = nSegTotal;                        % = n_TF * Nsegm
Xe_all = Xs;  Ye_all = Ys;  Ze_all = Zs;   % eval points = all segment midpoints

% ATTENZIONE: Xs = x_mid(:) usa il layout colonna-major [n_TF x Nsegm] -> indicizzazione segment-major.
% Costruiamo quindi gli indici dei segmenti per ciascuna TF con passo n_TF.
segIndexByCoil = cell(n_TF,1);
for q = 1:n_TF
    segIndexByCoil{q} = q:n_TF:Ne_all;    % tutti i segmenti appartenenti alla TF q
end

% Corrente per ciascuna TF (qui tutte uguali, ma puoi personalizzare)
I_spira = Iop * ones(n_TF,1);

% Heuristica memoria (chunking 2D: blocchi di evaluation x blocchi di sorgenti)
Ktmp = 9; bytes = 8; frac = 0.5;
try
    m = memory; M_free = frac * m.MaxPossibleArrayBytes;
catch
    M_free = 1.0 * 1024^3;   % fallback 1 GiB
end
target     = max(1, floor(M_free/(bytes*Ktmp)));
eval_chunk = min(Ne_all, max(256, floor(sqrt(target))));
eval_chunk = eval_chunk - mod(eval_chunk,64); eval_chunk = max(eval_chunk,256);

L = zeros(n_TF, n_TF);      % output inductance matrix [H]

for q = 1:n_TF
    idx_src = segIndexByCoil{q};
    nSrc_q  = numel(idx_src);
    const_A = mu0/(4*pi) * 1; % su 1 A

    % Sorgenti (solo coil q) come vettori riga per broadcasting
    Xs_q  = Xe_all(idx_src).';
    Ys_q  = Ye_all(idx_src).';
    Zs_q  = Ze_all(idx_src).';
    DLx_q = DLx   (idx_src).';
    DLy_q = DLy   (idx_src).';
    DLz_q = DLz   (idx_src).';

    % Dimensione blocco sorgenti coerente con il budget
    src_block = min(nSrc_q, max(256, floor(target / max(eval_chunk,1))));
    src_block = src_block - mod(src_block,64); src_block = max(src_block,256);

    % Accumulatori del potenziale vettore A_q su TUTTI i punti di valutazione
    Ax_all = zeros(Ne_all,1);
    Ay_all = zeros(Ne_all,1);
    Az_all = zeros(Ne_all,1);

    for eStart = 1:eval_chunk:Ne_all
        eStop = min(eStart+eval_chunk-1, Ne_all);

        Xe_blk = Xe_all(eStart:eStop);
        Ye_blk = Ye_all(eStart:eStop);
        Ze_blk = Ze_all(eStart:eStop);

        for sStart = 1:src_block:nSrc_q
            sStop = min(sStart+src_block-1, nSrc_q);

            Xsb  = Xs_q (sStart:sStop);
            Ysb  = Ys_q (sStart:sStop);
            Zsb  = Zs_q (sStart:sStop);
            DLxb = DLx_q(sStart:sStop);
            DLyb = DLy_q(sStart:sStop);
            DLzb = DLz_q(sStart:sStop);

            Rx = Xe_blk - Xsb;   Ry = Ye_blk - Ysb;   Rz = Ze_blk - Zsb;
            R2 = Rx.^2 + Ry.^2 + Rz.^2;
            R  = sqrt(R2 + a2_over_sqrte);  % distanza regolarizzata

            Ax_all(eStart:eStop) = Ax_all(eStart:eStop) + const_A * sum( DLxb ./ R , 2);
            Ay_all(eStart:eStop) = Ay_all(eStart:eStop) + const_A * sum( DLyb ./ R , 2);
            Az_all(eStart:eStop) = Az_all(eStart:eStop) + const_A * sum( DLzb ./ R , 2);
        end
    end

    % Colonna q della matrice: per ogni coil p somma dl·A_q sui suoi segmenti
    for p = 1:n_TF
        idx_p = segIndexByCoil{p};
        sdot  = DLx(idx_p).*Ax_all(idx_p) + DLy(idx_p).*Ay_all(idx_p) + DLz(idx_p).*Az_all(idx_p);
        L(p,q) = sum(sdot);      % [H/A]
    end
end

% Opzionale: forza simmetria per ridurre rumore numerico
Lt = L*(N_spire)^2;
E = (1/2*I_spira'*Lt*I_spira)*1e-9; % [GJ]

% Stampa rapida
fprintf('\nInductance matrix L [%dx%d] in H (symmetrized):\n', n_TF, n_TF);
disp(L);
disp(E);

% Opzionale: energia magnetica per le correnti attuali
Ivec = I_spira;               % [n_TF x 1]
Wmag = 0.5 * Ivec.' * L * Ivec * 1e-9;  % [GJ]
fprintf('Magnetic energy with current pattern I: W = %.6g GJ\n', Wmag);

%% ------------------------ Biot–Savart: B at coil midpoints (Hurwitz regularized)
% B(P) = μ0/(4π) * I * Σ [ dl × R ] / (|R|^2 + a^2/√e)^(3/2)
% Evaluates the regularized field at each segment midpoint of each coil.

% Flatten all coil segments into 1D arrays (length = n_TF * Nsegm)
Xs  = x_sour(:);         % source segment midpoints (X)
Ys  = y_sour(:);         % source segment midpoints (Y)
Zs  = z_sour3(:);        % source segment midpoints (Z)
DLx = dlxm(:);           % segment vector components (X)
DLy = dlym(:);           % segment vector components (Y)
DLz = dlzm(:);           % segment vector components (Z)
nSegTotal = numel(Xs);

% Output arrays for B at each midpoint of each coil
Bx = zeros(n_TF, Nsegm);
By = zeros(n_TF, Nsegm);
Bz = zeros(n_TF, Nsegm);

% Blocked summation to manage memory
K = 9;                      % # matrici temporanee "grandi" in contemporanea (~8–10)
try
  m = memory;                                % solo Windows
  M_free = 0.5 * m.MaxPossibleArrayBytes;    % prendi ~50% come budget
  block = floor(M_free / (8 * K * Nsegm));
  block = max(256, block - mod(block,256));
catch
  % fallback al valore fisso/empirico
  block = 4000;
end

for ii = 1:n_TF
    % Evaluation midpoints for this coil
    Xe = x_mid(ii,:).';  
    Ye = y_mid(ii,:).';  
    Ze = z_mid3(ii,:).';  

    % Accumulators for B components
    Bx_i = zeros(Nsegm,1); 
    By_i = zeros(Nsegm,1); 
    Bz_i = zeros(Nsegm,1);

    for startIdx = 1:block:nSegTotal
        stopIdx = min(startIdx+block-1, nSegTotal);

        % Source segment positions and dl vectors for this block
        Xsb = Xs(startIdx:stopIdx).';   
        Ysb = Ys(startIdx:stopIdx).';
        Zsb = Zs(startIdx:stopIdx).';
        DLxb= DLx(startIdx:stopIdx).';
        DLyb= DLy(startIdx:stopIdx).';
        DLzb= DLz(startIdx:stopIdx).';

        % Vector from source midpoint to evaluation point
        Rx = Xe - Xsb;   
        Ry = Ye - Ysb;
        Rz = Ze - Zsb;

        % Regularized distance cubed
        R2 = Rx.^2 + Ry.^2 + Rz.^2;
        R3 = (R2 + a2_over_sqrte).^(3/2);

        % dl × R  (source dl crossed with displacement vector)
        cx =  DLyb.*Rz - DLzb.*Ry;
        cy =  DLzb.*Rx - DLxb.*Rz;
        cz =  DLxb.*Ry - DLyb.*Rx;

        % Accumulate contribution from this block
        Bx_i = Bx_i + const * sum(cx ./ R3, 2);
        By_i = By_i + const * sum(cy ./ R3, 2);
        Bz_i = Bz_i + const * sum(cz ./ R3, 2);
    end

    % Store results for this coil
    Bx(ii,:) = Bx_i; 
    By(ii,:) = By_i; 
    Bz(ii,:) = Bz_i;
end

% Magnitude of B at each midpoint along each coil
Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);  % [n_TF x Nx]

%% ------------------------ Plot |B| along one coil vs poloidal arclength
coil_to_plot = 6;   % choose coil index

% Poloidal arclength along coil centerline (midpoints)
s_end = [0; cumsum(dl2D(:))];                 
s_mid = 0.5*(s_end(1:end-1) + s_end(2:end));  

figure('Name','|B| along a coil midpoint path');
plot(s_mid, Bmag(coil_to_plot,:),'LineWidth',1.5);
xlabel('Poloidal arclength s along coil [m]');
ylabel('|B| at coil midpoints [T]');
title(sprintf('Hurwitz-regularized |B| along coil %d', coil_to_plot));
grid on;

%% ------------------------ Quiver 3D of B along the D (selected coil)
coil_to_plot = 1:12;  % scegli la bobina
skip = 2;             % dirada le frecce per leggibilità (aumenta se è fitto)
idx  = 1:skip:Nsegm;

xq = x_mid(coil_to_plot, idx);
yq = y_mid(coil_to_plot, idx);
zq = z_mid3(coil_to_plot, idx);

u = Bx(coil_to_plot, idx);
v = By(coil_to_plot, idx);
w = Bz(coil_to_plot, idx);

Bsum = sqrt(u.^2+v.^2+w.^2);

figure('Name','B field quiver along D (3D)','units','normalized','outerposition',[0.25 0.1 0.5 0.9])
hold on; axis equal; grid on; view(3);
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title(sprintf('Magnetic field direction'));
plot3(xq, yq, zq, '.k', 'LineWidth', 0.2);
q = quiver3(xq, yq, zq, u, v, w, 1, 'b');
[q]  = ColorQuiver(q);
 set(gca,'FontSize',20);
colormap(jet);    
c = colorbar;
c.Label.FontSize = 20;
c.Label.Color = 'k';
c.Label.Rotation = 90;
c.Label.String = 'B [T]';
c.Label.Interpreter = 'latex';
c.TickLabelInterpreter = 'latex';
clim([min(min(Bsum)) max(max(Bsum))]);
MAXMIN=get(c,'Limits');
T = linspace(MAXMIN(1),MAXMIN(2),6);
set(c,'Ticks',T,'Fontsize',20)
grid minor  
hold off

%% ===== Campo magnetico lungo la linea r(t) = r1 + t*(r2-r1), t in [0,1] =====
% INPUT richiesti:
r1 = [Ric 0 0];
r2 = [Rec 0 0];
Nline = 200;

t  = linspace(0,1,Nline).';                % [Nline x 1]
dr = (r2 - r1);                            % vettore direzione linea
Xe = r1(1) + t * dr(1);                    % [Nline x 1]
Ye = r1(2) + t * dr(2);
Ze = r1(3) + t * dr(3);

% Metti le sorgenti come righe per il broadcasting
Xsr = Xs.';  Ysr = Ys.';  Zsr = Zs.';
DLxr = DLx.'; DLyr = DLy.'; DLzr = DLz.';

% R = r_eval - r_src
Rx = Xe - Xsr;           % [Nline x nSegTotal]
Ry = Ye - Ysr;
Rz = Ze - Zsr;

% Distanza regolarizzata
R2 = Rx.^2 + Ry.^2 + Rz.^2;
R3 = (R2 + a2_over_sqrte).^(3/2);

% dl × R  (prodotto vettoriale)
cx =  DLyr.*Rz - DLzr.*Ry;
cy =  DLzr.*Rx - DLxr.*Rz;
cz =  DLxr.*Ry - DLyr.*Rx;

% Somma su tutte le sorgenti (sulle colonne)
Bx_line = const * sum(cx ./ R3, 2);   % [Nline x 1]
By_line = const * sum(cy ./ R3, 2);
Bz_line = const * sum(cz ./ R3, 2);

% Modulo del campo
Bmag_line = sqrt(Bx_line.^2 + By_line.^2 + Bz_line.^2);

% Esempio di plot rapido:
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9]); 
plot(Xe, Bmag_line,'LineWidth',3); xlabel('s [m]'); ylabel('|B| [T]'); grid on; hold on
xline(R0, '--', 'R0', 'LineWidth',1.2, 'Color',[0 0 0]); 
xline(Ri, '--', 'Ri', 'LineWidth',1.2, 'Color',[0 0 0]); 
xline(Ric, '--', 'Ric', 'LineWidth',1.2, 'Color',[0 0 0]); 
xline(Re, '--', 'Re', 'LineWidth',1.2, 'Color',[0 0 0]); 
xline(Rec, '--', 'Rec', 'LineWidth',1.2, 'Color',[0 0 0]); 
% plot(x_mid(1, idx), z_mid3(1, idx)+max(Bmag_line)/2, 'k', 'LineWidth', 0.2);
hold off; 

%% ===== Campo magnetico sezione poloidale =====
% INPUT richiesti:

x_list = linspace(0, 6, 200);     % [m] normale poloidale
z_list = linspace(-4, +4, 200);   % [m] toroidale (spessore)
[X, Z] = meshgrid(x_list, z_list);
Xe = X(:);                                % [Ncond x 1]
Ye = Z(:).*0;                             % [Ncond x 1]
Ze = Z(:);

% Metti le sorgenti come righe per il broadcasting
Xsr = Xs.';  Ysr = Ys.';  Zsr = Zs.';
DLxr = DLx.'; DLyr = DLy.'; DLzr = DLz.';

% R = r_eval - r_src
Rx = Xe - Xsr;           % [Nline x nSegTotal]
Ry = Ye - Ysr;
Rz = Ze - Zsr;

% Distanza regolarizzata
R2 = Rx.^2 + Ry.^2 + Rz.^2;
R3 = (R2 + a2_over_sqrte).^(3/2);

% dl × R  (prodotto vettoriale)
cx =  DLyr.*Rz - DLzr.*Ry;
cy =  DLzr.*Rx - DLxr.*Rz;
cz =  DLxr.*Ry - DLyr.*Rx;

% Somma su tutte le sorgenti (sulle colonne)
Bx_line = const * sum(cx ./ R3, 2);   % [Nline x 1]
By_line = const * sum(cy ./ R3, 2);
Bz_line = const * sum(cz ./ R3, 2);

% Modulo del campo
Bmag_line = sqrt(Bx_line.^2 + By_line.^2 + Bz_line.^2);

%%
figure('units','normalized','outerposition',[0.25 0.1 0.5 0.9]); 
pb = reshape(Bmag_line,size(Z));
surf(X,Z,pb)
view(2)
set(gca,'FontSize',20);
colormap(jet)
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
axis equal
xlabel('[m]','fontsize',20,'Interpreter','latex')
ylabel('[m]','fontsize',20,'Interpreter','latex')
zlabel('[m]','fontsize',20,'Interpreter','latex')
title('|B| sulla sezione ');



