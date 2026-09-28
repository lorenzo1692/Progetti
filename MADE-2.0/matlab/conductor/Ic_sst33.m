function [Ic,FA]=Ic_sst33(B,T,theta,opt)
%theta=0 perp field, theta=90 parr field

    if nargin<4
        opt=[3,4];
    end
    
    switch opt(1) %B-T dependence
        case 1 %first approach
            FF=jc_ybco(B,T,'sst2')*3.3e-9;
        case 2 %more accurate
            FF=jc_ybco(B,T,'sst')/jc_ybco(0,77,'sst')*43*3.3;
        case 3 %SO'23 4 mm wide
            FF=jc_ybco(B,T,'sst')/jc_ybco(0,77,'sst')*43*4*1.92;
    end
    
    switch opt(2) %angular dependence
        case 1 %interp 4.2 K curve
            data=load('SP_angular_4.2_DU.mat');
            FA=interp1(data.angle,data.alpha,theta/180*pi);
        case 2 %analytical
            FA=1+3.222*exp((theta-90)/8.234)+1.278*exp((theta-90)/0.5901);
%             a=0.6084; c=0.2315; r=5.4041;
%             FA=((cosd(theta).^2+a*sind(theta).^2)./(cosd(theta).^2+a/r^(1/c)*sind(theta).^2)).^c;
        case 3 %interp var T curves (B influence neglected)
            data=load('SST_angular_low_T_Wimbush+DU.mat');
            FA=griddata(data.angle,data.temp,data.Icnorm,theta,T);
        case 4 %SO'23 analytical
            FA=1+4.482*exp((theta-90)/10.03);
    end
    
    Ic=FF.*FA;