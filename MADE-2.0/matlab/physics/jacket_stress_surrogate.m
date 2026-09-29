function [Pm, PmPb] = jacket_stress_surrogate(sigma_nom, n_turns, B_layer, is_transition, S_z) %#ok<INUSD>
%JACKET_STRESS_SURROGATE Fast jacket primary Pm and Pm+Pb per layer (scf_model = 1).
%
%   Calibration on the 2D FE (WP_MECH_SURROGATE) in progress: this
%   function is a placeholder and stops with an explicit error, so that
%   scf_model = 1 cannot be used before it is calibrated and validated.
%   scf_model = 0 (default) does not call it.
error('jacket_stress_surrogate:not_calibrated', ...
    'The jacket stress surrogate is not calibrated yet: use scf_model = 0.');
end
