function k = axial_load_factor(p)
%AXIAL_LOAD_FACTOR Factor on the vertical force T_bf of the inner leg (p.axial_load_factor, default 1).
%
%   k = AXIAL_LOAD_FACTOR(p) is the ratio between the vertical force that
%   the inner leg of one coil really carries and the bending-free formula
%   T_bf = 0.5 k_bf n_TF (N I)^2 mu0/(2 pi) of the 2D sizing, which splits
%   the vertical Lorentz force of the upper half equally between the inner
%   and the outer leg. The 3D beam + shell global model of the whole TF
%   system (TF3D_GLOBAL_FROM_DESIGN) gives it: the inner leg, stiffer than
%   the outer leg and wedged in the vault, takes more than half.
%
%   Every sizing and verification step multiplies T_bf by k: the scan
%   (jacket axial stress, SIZE_CASE_VAULT), the layered model and the case
%   surrogate, the jacket surrogate, REFINE_JACKET_FE and the 2D FE
%   (WP_MECH_SURROGATE). Input files without the parameter: k = 1, the
%   results of the previous versions.

k = 1;
if isfield(p, 'axial_load_factor') && ~isempty(p.axial_load_factor) && isfinite(p.axial_load_factor)
    k = p.axial_load_factor;
end
if ~(isscalar(k) && k > 0)
    error('axial_load_factor:value', 'p.axial_load_factor must be a positive scalar.');
end
end
