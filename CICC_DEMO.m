function [type_cable,N_Cu,N_Sc,N_tot,S_Cable,S_REBCO,S_Cu_HTS,THS,mat] = CICC_DEMO(B,Iop,Tau_discharge,WP_SC_type)

CunonCu         = 1;
S_tapes         = (4*1e-7);        % [m2]
d_fili          = 0.001;         % [m] 
theta           = 30;
T_dim           = 20;
d_cc            = 0.01;           % Diametro del cooling channel [m] 
cos_theta       = 0.97;            % cos(theta) tiene conto del fatto che i fili sono twistati, quindi la sezione effettiva di fili è maggiore
VF              = 0.7;             % Void Fraction nel conduttore 
N_fili          = [1500 1440 1350 1296 1200 1152 1080 972 960 900 864 810 768 720 675 648 540 486 360 324 300 216 180 162 144]; 
N_Cu            = 0;
N_tot           = 0;
type_cable      = {'X'};
S_REBCO         = 0;
S_Cu_HTS        = 0;
S_Cable         = 0;
S_Cable_0       = 0;
THS             = 0;
                
if WP_SC_type == 101
    type_cable = {'HTS'};
    % [y,Ic_sc]  = Jc_REBCO(B,S_tapes); % Critical current REBCO 
    Ic_sc = Ic_sst33(T_dim,B,theta,[3,4]);
    mat = 2; % '0 = Nb3Sn' - '1 = NbTi' - '2 = REBCO';
    Tlim = 150;
elseif WP_SC_type == 100
    type_cable = {'LTS'};
    Ic_sc = Ic_Nb3Sn_WST(B,d_fili*1e3); % Critical current Nb3Sn
    mat = 0; % '0 = Nb3Sn' - '1 = NbTi';
    Tlim = 250;
    if B < 0
        Ic_sc = Ic_NbTi_TF(B,d_fili*1e3); % Critical current NbTi
        mat = 1; % '0 = Nb3Sn' - '1 = NbTi';
        Tlim = 250;
    end
elseif WP_SC_type == 102
    if B > 15
        type_cable = {'HTS'};
        % [~,Ic_sc]  = Jc_REBCO(B,S_tapes); % Critical current REBCO 
        Ic_sc = Ic_sst33(T_dim,B,theta,[3,4]);
        mat = 2; % '0 = Nb3Sn' - '1 = NbTi' - '2 = REBCO';
        Tlim = 150;
    else
        type_cable = {'LTS'};
        Ic_sc = Ic_Nb3Sn_WST(B,d_fili*1e3); % Critical current Nb3Sn
        mat = 0; % '0 = Nb3Sn' - '1 = NbTi';
        Tlim = 250;
        if B < 0
            Ic_sc = Ic_NbTi_TF(B,d_fili*1e3); % Critical current NbTi
            mat = 1; % '0 = Nb3Sn' - '1 = NbTi';
            Tlim = 250;
        end
    end
end

N_Sc = ceil(Iop/Ic_sc); % Needed wires or tapes

if N_Sc > max(N_fili) && WP_SC_type == 102
    type_cable = {'HTS'};
    % [~,Ic_sc]  = Jc_REBCO(B);              % Critical current REBCO 
    Ic_sc = Ic_sst33(20,B,theta,[3,4]);
    mat = 2; % '0 = Nb3Sn' - '1 = NbTi' - '2 = REBCO';
    N_Sc = ceil(Iop/Ic_sc); % Needed wires or tapes
    Tlim = 250;
elseif N_Sc > max(N_fili) && WP_SC_type ~= 102  
    return
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
        THS(j) = heat_balance_CICC_ODE(N_Sc,N_Cu0(j),d_fili,CunonCu,Iop,B,Tau_discharge,2,d_cc,VF,cos_theta,S_tapes);           
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


