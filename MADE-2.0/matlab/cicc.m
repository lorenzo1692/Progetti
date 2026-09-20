function [type_cable,N_Cu,N_Sc,N_tot,S_Cable,S_REBCO,S_Cu_HTS,THS,mat,Ic_sc] = cicc(B,Iop,Tau_discharge,WP_SC_type)

CunonCu         = 1;
S_tapes         = (4*1e-7);        % [m2]
d_fili          = 0.001 + 0*0.00082;         % [m] 
theta           = 30;
T_dim           = 20;
d_cc            = 0.005;           % Diametro del cooling channel [m] 
cos_theta       = 0.97;            % cos(theta) tiene conto del fatto che i fili sono twistati, quindi la sezione effettiva di fili è maggiore
VF              = 0.7;             % Void Fraction nel conduttore 
% N_fili          = [1500 1440 1350 1296 1200 1152 1080 972 960 900 864 810 768 720 675 648 540 486 360 324 300 216 180 162 144]; 
N_fili          = 1:1500;
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
        if B > 15
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
            [Ic_sc, type_cable, mat, Tlim] = select_HTS(B,20,theta);

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
    Tlim = 150;
end


function [Ic_sc,type_cable,mat,Tlim] = select_LTS(B,d_fili)
% LTS selection:
%   B < 5 T  -> NbTi
%   B >= 5 T -> Nb3Sn
    type_cable = {'LTS'};
    Tlim = 250;
    if B < 6
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
div = 10;
N_Cu0 = linspace(1,10000,div); 
while true
    THS = zeros(size(N_Cu0,2),1);
    for j=1:size(N_Cu0,2)         
        THS(j) = heat_balance_cicc_ode(N_Sc,N_Cu0(j),d_fili,CunonCu,Iop,B,Tau_discharge,mat,d_cc,VF,cos_theta,S_tapes);           
    end    
    [~,indx] = min(abs(THS-Tlim)); 
    %
    if THS(indx) > Tlim
        b = round(N_Cu0(min(max(1,indx+1),div))); a = round(N_Cu0(indx));
    else
        a = round(N_Cu0(max(1,indx-1))); b = round(N_Cu0(indx));
    end          
    %
    if abs(Tlim-THS(indx)) <= 5 || b-a <= 1 || a > 980 
        N_Cu = ceil(N_Cu0(indx));
        THS = THS(indx);
        break
    end
    N_Cu0 = linspace(a,b,div);
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


