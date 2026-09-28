function [B_layer, k] = wp_field_calibrated(cal, W, Iop, n_spire)
%WP_FIELD_CALIBRATED Per-layer peak field from the start-of-scan calibration.
%
%   [B_layer, k] = WP_FIELD_CALIBRATED(cal, W, Iop, n_spire) returns the
%   field to size each layer with, B_k = k(W, Iop) * B_Ampere * f(s_k),
%   from the tables of WP_FIELD_CALIBRATION: W toroidal WP width [m], Iop
%   operating current [A], n_spire(k) = turns from layer k inward (as in
%   SCAN_WP_DESIGNS), s_k = n_spire(k)/n_spire(1). Linear interpolation in
%   W and Iop, clamped to the calibrated range (no extrapolation).

Wc = min(max(W, min(cal.W)), max(cal.W));
Ic = min(max(Iop, min(cal.Iop)), max(cal.Iop));
[Wg, ix] = sort(cal.W);
kw = cal.k(ix,:); fw = cal.f(ix,:,:);
k = interp2_lin(cal.Iop, Wg, kw, Ic, Wc);
ns = numel(cal.s); fs = zeros(1, ns);
for q = 1:ns
    fs(q) = interp2_lin(cal.Iop, Wg, fw(:,:,q), Ic, Wc);
end
s = n_spire/n_spire(1);
B_layer = k*cal.B_amp*interp1(cal.s, fs, min(max(s, 0), 1), 'linear');
end

function v = interp2_lin(x, y, Z, xq, yq)
% bilinear on a grid that may have a single point in either direction
if numel(x) == 1 && numel(y) == 1, v = Z(1); return, end
if numel(x) == 1, v = interp1(y, Z(:,1), yq, 'linear'); return, end
if numel(y) == 1, v = interp1(x, Z(1,:), xq, 'linear'); return, end
v = interp2(x, y, Z, xq, yq, 'linear');
end
