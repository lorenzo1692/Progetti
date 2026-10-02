function [type_cable,N_Cu,N_Sc,N_tot,S_Cable,S_REBCO,S_Cu_HTS,THS,mat,Ic_sc] = cicc(B,Iop,Tau_discharge,WP_SC_type,THS_max_LTS,THS_max_HTS)
%CICC Size one CICC grade (SC strands/tapes + segregated Cu) at field B.
%
%   The hot-spot limit is type-specific: THS_max_LTS for LTS (default
%   250 K) and THS_max_HTS for HTS (default 150 K). Pass them explicitly
%   (e.g. p.THS_max_LTS/p.THS_max_HTS from the input Excel) to override
%   the defaults in cicc_params.m. All other constants: cicc_params.m.

cp = cicc_params();
if nargin < 5 || isempty(THS_max_LTS), THS_max_LTS = cp.THS_max_LTS; end
if nargin < 6 || isempty(THS_max_HTS), THS_max_HTS = cp.THS_max_HTS; end

CunonCu         = cp.CunonCu;
S_tapes         = cp.S_tapes;        % [m2]
d_fili          = cp.d_fili;         % [m]
theta           = cp.theta;
T_dim           = cp.T_dim;
d_cc            = cp.d_cc;           % cooling channel diameter [m]
cos_theta       = cp.cos_theta;      % twisted strands: effective strand cross-section is larger
VF              = cp.VF;             % void fraction
N_fili          = 1:cp.N_fili_max;
N_Cu            = 0;
N_tot           = 0;
type_cable      = {'X'};
S_REBCO         = 0;
S_Cu_HTS        = 0;
S_Cable         = 0;
S_Cable_0       = 0;
THS             = 0;

%% Superconductor selection

switch WP_SC_type

    case 101
        % Pure HTS / REBCO option
        [Ic_sc, type_cable, mat, Tlim] = select_HTS(B,T_dim,theta);

    case 100
        % Pure LTS option
        [Ic_sc, type_cable, mat, Tlim] = select_LTS(B,d_fili);

    case 102
        % Hybrid option:
        % HTS above 15 T, LTS below 15 T
        if B > cp.B_hybrid_HTS
            [Ic_sc, type_cable, mat, Tlim] = select_HTS(B,T_dim,theta);
        else
            [Ic_sc, type_cable, mat, Tlim] = select_LTS(B,d_fili);
        end

    otherwise
        error('Unknown WP_SC_type: %d', WP_SC_type);

end


%% Number of required strands / tapes

N_Sc = ceil(Iop/Ic_sc);


%% Check maximum allowed number of strands/tapes

if N_Sc > max(N_fili)

    switch WP_SC_type

        case 102
            % In hybrid mode, if the LTS solution requires too many strands,
            % force the use of HTS.
            [Ic_sc, type_cable, mat, Tlim] = select_HTS(B,T_dim,theta);

            N_Sc = ceil(Iop/Ic_sc);

        otherwise
            % For non-hybrid options, the design point is not feasible.
            return

    end

end


%% Local functions

function [Ic_sc,type_cable,mat,Tlim] = select_HTS(B,T,theta)
% REBCO / HTS critical current
    type_cable = {'HTS'};
    % Ic_sst33 signature: Ic_sst33(B,T,theta,opt)
    Ic_sc = Ic_sst33(B,T,theta,[3,4]);
    mat = 2;     % 0 = Nb3Sn, 1 = NbTi, 2 = REBCO
    Tlim = THS_max_HTS;
end


function [Ic_sc,type_cable,mat,Tlim] = select_LTS(B,d_fili)
% LTS selection:
%   B < cp.B_NbTi_max  -> NbTi
%   B >= cp.B_NbTi_max -> Nb3Sn
    type_cable = {'LTS'};
    Tlim = THS_max_LTS;
    if B < cp.B_NbTi_max
        Ic_sc = Ic_NbTi(B,d_fili*1e3);
        mat = 1;     % NbTi
    else
        [Ic_sc,~] = Ic_Nb3Sn(B,d_fili*1e3,'DTT_TF_KAT');
        mat = 0;     % Nb3Sn
    end
end   
             
% if mat ~= 2 
%     a_ = N_fili > N_Sc;
%     z_ = N_fili(a_);
%     N_Sc = min(z_);  
% end

%% Heat Balance
% Smallest copper strand count N_Cu with the hot spot THS(N_Cu) <= Tlim,
% within 5 K below the limit (THS falls monotonically with N_Cu). Bracketed
% search on log(N_Cu) in [1, N_Cu_max]: geometric bisection until both ends
% are known, then regula falsi (Illinois) on log(N_Cu); stops when the safe
% end is within 5 K of the limit or the bracket is one strand. About 5-8
% transient solutions instead of the 30-40 of the former 10-point grid
% refinement, which also stopped on the coarse grid (up to ~1100 strands of
% extra copper) once the bracket passed 980 strands.
ths = @(n) heat_balance_cicc_ode(N_Sc,n,d_fili,CunonCu,Iop,B,Tau_discharge,mat,d_cc,VF,cos_theta,S_tapes,cp.Tau_delay,2*Tlim);
N_Cu_max = 10000;
hi = N_Cu_max; T_hi = ths(hi);
if T_hi > Tlim
    % not even the largest copper keeps the hot spot below the limit: the
    % scan rejects the candidate on THS
    N_Cu = hi; THS = T_hi;
else
    lo = 1; T_lo = Inf;              % N_Cu = 1: runaway, above the limit
    side = 0;
    while hi - lo > 1 && T_hi < Tlim - 5
        if isfinite(T_lo)
            % regula falsi on log(N) with the Illinois correction
            w_lo = T_lo - Tlim; w_hi = T_hi - Tlim;
            if side == 1, w_lo = w_lo/2; elseif side == -1, w_hi = w_hi/2; end
            x = log(lo) + (log(hi) - log(lo))*w_lo/(w_lo - w_hi);
            n = round(exp(x));
        else
            n = round(sqrt(lo*hi));
        end
        n = min(max(n, lo + 1), hi - 1);
        T_n = ths(n);
        if T_n <= Tlim
            hi = n; T_hi = T_n; side = 1;
        else
            lo = n; T_lo = T_n; side = -1;
        end
    end
    N_Cu = hi; THS = T_hi;
end
%%
N_tot = N_Sc+N_Cu;

if mat == 2 
    S_REBCO   = N_Sc*S_tapes;
    S_Cu_HTS  =  N_Cu*pi*d_fili^2/4;
    S_Cable_0 = (S_REBCO + S_Cu_HTS + (pi*d_cc^2/4)); 
elseif N_Sc < max(N_fili) 
    S_Cable_0 = (N_tot*pi*d_fili^2/4/cos_theta/VF+(pi*d_cc^2/4)); 
end

D_eqv = (S_Cable_0*4/pi)^0.5;           % Diametro cavo equivalente da area eqv
A_w  = ((D_eqv+0.0004)^2-D_eqv^2)*pi/4; % Considereo l'area del wrapping in acciaio di spessore 0.4 mm
S_Cable =  pi/4*D_eqv^2+A_w;            % Correggo area equivalente del cavo LTS                                     
end


