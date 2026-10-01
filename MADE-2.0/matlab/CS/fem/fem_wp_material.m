function mat = fem_wp_material(wp, p)
%FEM_WP_MATERIAL Homogenized orthotropic properties of one CICC winding-pack cell.
%
%   mat = FEM_WP_MATERIAL(wp, p) smears the jacket, cable and turn
%   insulation of a conductor cell into one orthotropic material for the
%   axisymmetric model (axes r, z, theta). Moduli in GPa come from p
%   (E_jckt, E_cbl_LTS, E_cbl_HTS, E_ins), Poisson ratio from p.fem_nu_wp.
%
%   Hoop (theta): the three materials share the same strain, so the
%   stiffness is the area-weighted sum (parallel model) - the same
%   k_jkt/k_tot ratio used by EQV_STRESS_COIL_RING_CICC.
%   Radial / axial: the cell is a network of two parallel paths, the wall
%   pair along the load direction and the inner chain
%   (insulation - jacket - cable - jacket - insulation) in series, the
%   same structure as the Ke_cavo stiffness of the analytic model, but
%   with each path's own cross-section.
%
%   Output mat:
%     D           4x4 constitutive matrix [Pa], stress order [sr sz sth trz]
%     Er, Ez, Eth, Grz   homogenized moduli [Pa]
%     f_hoop      steel-to-average stress ratio, hoop (= E_jckt/Eth)
%     f_z, f_r    jacket-wall stress per unit average stress, vertical and
%                 radial (load share of the walls x cell width / wall
%                 thickness; used to recover jacket stress from the
%                 smeared FE stress)
%     E_j, E_c, E_i  constituent moduli [Pa]

Ej = p.E_jckt*1e9; Ei = p.E_ins*1e9;
if strcmp(wp.type_cable, 'HTS')
    Ec = p.E_cbl_HTS*1e9;
else
    Ec = p.E_cbl_LTS*1e9;
end

w = wp.Cond_w; h = wp.Cond_h; t = wp.JT; ti = wp.tins;
scw = wp.SC_w; sch = wp.SC_h;

% hoop: parallel
Eth = (wp.S_JT*Ej + wp.S_Cable*Ec + wp.S_ins*Ei)/wp.A_cell;

% radial (load along w): walls = top/bottom pair (thickness t+ti, length w),
% chain = side walls + cable across the inner height
hi = h - 2*t - 2*ti;
Kw_r = (Ej*2*t + Ei*2*ti)/w;
Kc_r = hi/(2*ti/Ei + 2*t/Ej + scw/Ec);
Er = (Kw_r + Kc_r)*w/h;
share_r = (Ej*2*t/w)/(Kw_r + Kc_r);

% axial (load along h)
wi = w - 2*t - 2*ti;
Kw_z = (Ej*2*t + Ei*2*ti)/h;
Kc_z = wi/(2*ti/Ei + 2*t/Ej + sch/Ec);
Ez = (Kw_z + Kc_z)*h/w;
share_z = (Ej*2*t/h)/(Kw_z + Kc_z);

nu = p.fem_nu_wp;
Grz = min(Er, Ez)/(2*(1+nu));

% orthotropic compliance, order [r z theta]; nu_ij/E_i = nu_ji/E_j
S = zeros(3);
S(1,1) = 1/Er; S(2,2) = 1/Ez; S(3,3) = 1/Eth;
S(1,2) = -nu/Ez; S(2,1) = S(1,2);
S(1,3) = -nu/Eth; S(3,1) = S(1,3);
S(2,3) = -nu/Eth; S(3,2) = S(2,3);
if any(eig(S) <= 0)
    error('fem_wp_material:not_positive_definite', ...
        'The homogenized compliance is not positive definite: check fem_nu_wp and the constituent moduli.');
end
D3 = inv(S);
D = zeros(4);
D(1:3, 1:3) = D3;
D(4,4) = Grz;

mat.D = D;
mat.Er = Er; mat.Ez = Ez; mat.Eth = Eth; mat.Grz = Grz;
mat.f_hoop = Ej/Eth;
mat.f_z = share_z*w/(2*t);
mat.f_r = share_r*h/(2*t);
mat.E_j = Ej; mat.E_c = Ec; mat.E_i = Ei;
end
