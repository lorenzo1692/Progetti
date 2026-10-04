function [Bin, Bout] = system_field_eval(sys, Ri, Re, h)
%SYSTEM_FIELD_EVAL Background Bz [T] on the inner and outer radius of a PF winding pack, per scenario.
%
%   [Bin, Bout] = SYSTEM_FIELD_EVAL(sys, Ri, Re, h) returns row vectors (one
%   value per scenario) with the axial field of the rest of the machine
%   (see SYSTEM_FIELD_SETUP), averaged over the winding-pack height h
%   (5 axial points over 90% of h, centred on the coil), at radius Ri
%   (Bin) and Re (Bout). Positive Bz is the sign of the self field of a
%   coil carrying a positive current.

z = sys.Zrow + linspace(-0.45, 0.45, 5)'*h;
Sc = 1e-32*ones(numel(sys.Rf), 1);
Bin = mean(xbz(Ri*ones(5,1), z, sys.Rf, sys.Zf, Sc)*sys.Ifil, 1);
Bout = mean(xbz(Re*ones(5,1), z, sys.Rf, sys.Zf, Sc)*sys.Ifil, 1);
end
