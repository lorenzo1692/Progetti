function field = tf3d_coil_filaments(row, p, geo, path_rz, n_fil_layer, field2d) %#ok<INUSL>
%   row is accepted for interface symmetry with the rest of MADE (every
%   coil3d/* driver-level function takes row,p) though not used directly
%   here (geo and field2d already carry everything derived from it).
%TF3D_COIL_FILAMENTS 3D field checks from a single-filament-per-coil model.
%
%   field = TF3D_COIL_FILAMENTS(row, p, geo, path_rz, n_fil_layer, field2d)
%   represents each of the n_TF coils by ONE filament following path_rz
%   (public convention: index 1 = inner leg top, index end = outer-leg
%   midplane; r,z of the current centroid, from the bending-free shape or
%   its three-arc fit), closed into a full loop exactly as
%   TF_BENDING_FREE_SHAPE's internal path, and:
%
%     - checks Ampere's law at (R0, z=0): a circular path there threads
%       the coils' bore (r1 < R0 < r2 at z=0), so plain (uncalibrated)
%       Biot-Savart from all n_TF coils must reproduce
%       mu0*n_TF*I_coil/(2*pi*R0) to high accuracy - this is an
%       independent, in-context check of BIOT_SAVART_SEGMENTS itself;
%     - the field ripple at the outboard midplane (R0+a, z=0), compared to
%       p.ripple (MADE's own target, used to size RTFo);
%     - the axial hoop tension T = Fz(upper half)/2 from the Lorentz force
%       on the (self, ideal filament + 2x other-coils Biot-Savart -
%       see TF_BENDING_FREE_SHAPE for why the factor of 2) field acting
%       along the coil's own path, compared to MADE's shell-model T_bf;
%     - a single scalar correction rho = B_3D/B_2D_ideal at the inner-leg
%       cross-section (same construction as the ripple/tension field),
%       applied to the ALREADY VALIDATED detailed 2D peak-per-turn field
%       (COMPUTE_DISCRETE_FIELD_PROFILE) to estimate the 3D peak per turn:
%       this avoids a full multi-turn x multi-coil 3D Biot-Savart (tens of
%       thousands of segments squared, not tractable for an interactive
%       tool) while still capturing the genuine 3D/finite-n_TF correction,
%       assumed uniform across the small WP cross-section relative to the
%       machine scale (checked by keeping the WP off the coil's own
%       curved-path source resolution, see n_seg in TF_BENDING_FREE_SHAPE).
%       n_fil_layer is accepted for a future per-layer refinement but not
%       used by this scalar-correction model (kept at whole-coil
%       resolution for speed); field2d is accepted for the same reason
%       even though only its peak-per-turn array is scaled here, by the
%       caller (see TF3D_FROM_DESIGN).
%
%   field fields:
%     B_phi_R0   - toroidal field at (R0,0,0) from all n_TF coils [T]
%     ripple     - (Bmax-Bmin)/(Bmax+Bmin) at the outboard midplane
%     T_axial    - hoop tension from the 3D Lorentz force [N]
%     rho_mid    - B_3D/B_2D_ideal at the inner-leg cross-section
%     peak_mid, peak2d_mid - rho_mid-scaled and original 2D peak field at
%                  that cross-section (informational; the per-turn scaled
%                  array is built by the caller)

Mu_0 = 4e-7*pi;
n_TF = geo.n_TF; I_coil = geo.I_coil;

[S1c, S2c] = closed_path(path_rz);                  % coil 1 (phi=0) closed path segments

%% Ampere's law check at (R0, 0, 0)
Pamp = [p.R0, 0, 0];
Bamp = [0 0 0];
for c = 1:n_TF
    a = (c-1)*2*pi/n_TF;
    [S1r, S2r] = rot_path(S1c, S2c, a);
    Bamp = Bamp + biot_savart_segments(Pamp, S1r, S2r, I_coil);
end
field.B_phi_R0 = Bamp(2);                            % toroidal component = y at phi=0

%% Ripple at the outboard midplane
a_minor = p.R0/p.A;
Prip0   = [p.R0 + a_minor, 0, 0];                    % phi = 0 (under a coil)
Prip1   = [(p.R0 + a_minor)*cos(pi/n_TF), (p.R0 + a_minor)*sin(pi/n_TF), 0]; % phi = pi/n_TF (between coils)
Bmax = [0 0 0]; Bmin = [0 0 0];
for c = 1:n_TF
    a = (c-1)*2*pi/n_TF;
    [S1r, S2r] = rot_path(S1c, S2c, a);
    Bmax = Bmax + biot_savart_segments(Prip0, S1r, S2r, I_coil);
    Bmin = Bmin + biot_savart_segments(Prip1, S1r, S2r, I_coil);
end
bmax = norm(Bmax); bmin = norm(Bmin);
field.ripple = (bmax - bmin)/(bmax + bmin);

%% Inner-leg cross-section correction factor rho = B_3D / B_2D_ideal
r1 = path_rz(1,1); z1 = path_rz(1,2);
Pmid = [r1, 0, z1];
Bmid = Mu_0*I_coil/(2*pi*r1);                        % self, ideal-filament term (see TF_BENDING_FREE_SHAPE)
for c = 2:n_TF
    a = (c-1)*2*pi/n_TF;
    [S1r, S2r] = rot_path(S1c, S2c, a);
    Bc = biot_savart_segments(Pmid, S1r, S2r, I_coil);
    Bmid = Bmid + 2*Bc(2);
end
B2d_ideal = Mu_0*n_TF*I_coil/(2*pi*r1);
field.rho_mid = Bmid/B2d_ideal;
field.peak2d_mid = max(field2d.B_discrete);
field.peak_mid = field.rho_mid*field.peak2d_mid;

%% Axial hoop tension: T = Fz(upper half)/2
mid_idx = find(path_rz(:,2) <= 0, 1);                % first upper-half point reaching z<=0 (outer leg midplane)
if isempty(mid_idx), mid_idx = size(path_rz,1); end
Su1 = S1c(1:mid_idx-1,:); Su2 = S2c(1:mid_idx-1,:);   % upper-half segments of coil 1
Pm = 0.5*(Su1 + Su2);
% self term at each segment midpoint: ideal filament, direction = local
% e_phi at phi=0 = global y (see TF_BENDING_FREE_SHAPE for the same idealization)
r_seg = sqrt(Pm(:,1).^2 + Pm(:,2).^2);
Btot = zeros(size(Pm));
Btot(:,2) = Mu_0*I_coil./(2*pi*r_seg);
for c = 2:n_TF
    a = (c-1)*2*pi/n_TF;
    [S1r, S2r] = rot_path(S1c, S2c, a);
    Bc = biot_savart_segments(Pm, S1r, S2r, I_coil);
    Btot = Btot + 2*Bc;
end
dl = Su2 - Su1;
F = I_coil*cross(dl, Btot, 2);                        % Lorentz force per segment [N]
field.T_axial = sum(F(:,3))/2;                        % Fz(upper half)/2
end

function [S1, S2] = closed_path(rz)
% Public-convention half-shape (index1=r1,z_top .. index end=r2,z=0) into
% a full closed loop: upper half, mirrored lower half, straight inner leg.
ru = rz(:,1); zu = rz(:,2);
rl = flipud(ru); zl = -flipud(zu);
rs = [ru; rl(2:end)]; zs = [zu; zl(2:end)];
rs = [rs; rs(1)]; zs = [zs; zs(1)];
x = rs; y = zeros(size(rs)); zc = zs;
S1 = [x(1:end-1), y(1:end-1), zc(1:end-1)];
S2 = [x(2:end),   y(2:end),   zc(2:end)];
end

function [S1r, S2r] = rot_path(S1, S2, a)
c = cos(a); s = sin(a);
S1r = [S1(:,1)*c - S1(:,2)*s, S1(:,1)*s + S1(:,2)*c, S1(:,3)];
S2r = [S2(:,1)*c - S2(:,2)*s, S2(:,1)*s + S2(:,2)*c, S2(:,3)];
end
