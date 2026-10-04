function [Ka, Kc] = ParametriY(a_, c_, w, t, Sa)
%PARAMETRIY Newman-Raju stress intensity shape factors for a surface flaw.
%
%   [Ka, Kc] = PARAMETRIY(a_, c_, w, t, Sa) returns the stress intensity
%   factors [MPa*sqrt(m)] at the deepest point (Ka, psi=pi/2) and at the
%   surface (Kc, psi=0) of a semi-elliptical surface crack (depth a_,
%   half-length c_) in a plate of width w and thickness t under remote
%   stress Sa.
%
%   Relocated unchanged from the CS legacy archive (ParametriY.m). Used by
%   FATIGUE_FCGR.

%% Ka
psi = pi/2;
if (a_/c_) <= 1
    Q=1+1.464*((a_/c_))^1.65;
    M1=1.13-0.09*(a_/c_);
    M2=-0.54+0.89/(0.2+(a_/c_));
    M3=0.5-1/(0.65+(a_/c_))+14*(1-(a_/c_))^24;
    g=1+(0.1+0.35*(a_/t)^2)*(1-sin(psi))^2;
    fphi=(((a_/c_))^2*(cos(psi))^2+(sin(psi))^2)^0.25;
    fw=sqrt(sec(pi*c_/w*sqrt(a_/t)));
    p=0.2+(a_/c_)+0.6*a_/t;
    G21=-1.22-0.12*(a_/c_);
    G22=0.55-1.05*((a_/c_))^0.75+0.47*((a_/c_))^1.5;
    H1=1-0.34*a_/t-0.11*((a_/c_))*(a_/t);
    H2=1+G21*(a_/t)+G22*(a_/t)^2;
    Hs=H1+(H2-H1)*(sin(psi)^p);
elseif (a_/c_) > 1
    Q=1+1.464*(c_/a_)^1.65;
    M1=sqrt(c_/a_)*(1+0.04*c_/a_);
    M2=0.2*(c_/a_)^4;
    M3=-0.11*(c_/a_)^4;
    g=1+(0.1+0.35*(c_/a_)*(a_/t)^2)*(1-sin(psi))^2;
    fphi=((c_/a_)^2*(sin(psi))^2+(cos(psi))^2)^0.25;
    fw=sqrt(sec(pi*c_/w*sqrt(a_/t)));
    p=0.2+c_/a_+0.6*a_/t;
    G11=-0.04-0.41*(c_/a_);
    G12=0.55-1.93*(c_/a_)^0.75+1.38*(c_/a_)^1.5;
    G21=-2.11+0.77*(c_/a_);
    G22=0.55-0.72*(c_/a_)^0.75+0.14*(c_/a_)^1.5;
    H1=1+G11*(a_/t)+G12*(a_/t)^2;
    H2=1+G21*(a_/t)+G22*(a_/t)^2;
    Hs=H1+(H2-H1)*(sin(psi)^p);
end
Fs = (M1+M2*(2*a_/t)^2+M3*(2*a_/t)^4)*g*fphi*fw;
Ka = Fs*(Sa+Hs*0)*sqrt(pi*a_/Q);

%% Kc
psi = 0;
if (a_/c_) <= 1
    Q=1+1.464*((a_/c_))^1.65;
    M1=1.13-0.09*(a_/c_);
    M2=-0.54+0.89/(0.2+(a_/c_));
    M3=0.5-1/(0.65+(a_/c_))+14*(1-(a_/c_))^24;
    g=1+(0.1+0.35*(a_/t)^2)*(1-sin(psi))^2;
    fphi=(((a_/c_))^2*(cos(psi))^2+(sin(psi))^2)^0.25;
    fw=sqrt(sec(pi*c_/w*sqrt(a_/t)));
    p=0.2+(a_/c_)+0.6*a_/t;
    G21=-1.22-0.12*(a_/c_);
    G22=0.55-1.05*((a_/c_))^0.75+0.47*((a_/c_))^1.5;
    H1=1-0.34*a_/t-0.11*((a_/c_))*(a_/t);
    H2=1+G21*(a_/t)+G22*(a_/t)^2;
    Hs=H1+(H2-H1)*(sin(psi)^p);
elseif (a_/c_) > 1
    Q=1+1.464*(c_/a_)^1.65;
    M1=sqrt(c_/a_)*(1+0.04*c_/a_);
    M2=0.2*(c_/a_)^4;
    M3=-0.11*(c_/a_)^4;
    g=1+(0.1+0.35*(c_/a_)*(a_/t)^2)*(1-sin(psi))^2;
    fphi=((c_/a_)^2*(sin(psi))^2+(cos(psi))^2)^0.25;
    fw=sqrt(sec(pi*c_/w*sqrt(a_/t)));
    p=0.2+c_/a_+0.6*a_/t;
    G11=-0.04-0.41*(c_/a_);
    G12=0.55-1.93*(c_/a_)^0.75+1.38*(c_/a_)^1.5;
    G21=-2.11+0.77*(c_/a_);
    G22=0.55-0.72*(c_/a_)^0.75+0.14*(c_/a_)^1.5;
    H1=1+G11*(a_/t)+G12*(a_/t)^2;
    H2=1+G21*(a_/t)+G22*(a_/t)^2;
    Hs=H1+(H2-H1)*(sin(psi)^p);
end
Fs = (M1+M2*(2*a_/t)^2+M3*(2*a_/t)^4)*g*fphi*fw;
Kc = Fs*(Sa+Hs*0)*sqrt(pi*a_/Q);
end
