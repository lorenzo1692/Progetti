!=======================================================================
! postpro/Plot.f - plots of the structural solution (PNG)
!
! Called by Solution_STR.f after SOLVE. Last load step: Tresca (S,INT),
! components and displacements of the whole section, of the jacket
! (material 2) and of the case (material 10), shear and normal stresses
! of the insulation (materials 3, 5, 7-8, 9), and the insulation fatigue
! index LHD_xy = S_X/sigma0 + (S_XY/tau0)^2. Contour scales end at
! 0.867 GPa (jacket) and 1 GPa (case). PNG files go to the working
! directory (.\output when run through RUN.dat).
!=======================================================================
/POST1
SET,LAST
ALLSEL,ALL,ALL
! /DEV,FONT,LEGEND,MENU 
/dev,font,1,Calibri,700,0,-21,0,0,,, 
       
/RGB,INDEX,100,100,100, 0   
/RGB,INDEX, 80, 80, 80,13   
/RGB,INDEX, 60, 60, 60,14   
/RGB,INDEX, 0, 0, 0,15  

/GFILE,1600 
/SHOW,PNG,,,9
 
/PLOPTS,INFO,2  
/PLOPTS,LEG1,1  
/PLOPTS,LEG2,0  
/PLOPTS,LEG3,1  
/PLOPTS,FRAME,0 
/PLOPTS,TITLE,0 
/PLOPTS,MINM,1  
/PLOPTS,FILE,0  
/PLOPTS,SPNO,0  
/PLOPTS,WINS,1  
/PLOPTS,WP,0
/PLOPTS,DATE,0  
/TRIAD,OFF 

/PSF,DEFA, ,1,0,0   
/PBF,DEFA, ,1   
/PSYMB,CS,0 
/PSYMB,NDIR,0   
/PSYMB,ESYS,0   
/PSYMB,LDIV,0   
/PSYMB,LDIR,0   
/PSYMB,ADIR,0   
/PSYMB,ECON,0   
/PSYMB,XNODE,0  
/PSYMB,DOT,1
/PSYMB,PCONV,   
/PSYMB,LAYR,0   
/PSYMB,FBCS,0 

/EDGE,1,0,45
/GLINE,1,-1 
 
/GRAPH,FULL
ALLSEL,ALL,ALL
!!!!!!!!!!!!!!!!!!!!!!!!!!!
RSYS,0
CSYS,0

PLNSOL,S,INT
*GET, SINT_min, PLNSOL, 0,MIN,
/CONTOUR,,9,SINT_min,(0.867e9-SINT_min)/9,0.867e9
PLNSOL,S,INT
/CONTOUR,,9,
PLNSOL,U,X
PLNSOL,U,Y
PLNSOL,U,SUM

ESEL,S,MAT,,2
NSLE,S,1
PLNSOL,S,INT
*GET, SINT_min, PLNSOL, 0,MIN,
/CONTOUR,,9,SINT_min,(0.867e9-SINT_min)/9,0.867e9
PLNSOL,S,INT
/CONTOUR,,9,
PLNSOL,S,X
PLNSOL,S,Y
PLNSOL,S,Z
PLNSOL,U,X
PLNSOL,U,Y
PLNSOL,U,SUM

! ESEL,S,MAT,,2,4
! NSLE,S,1
! PLNSOL,S,INT
! *GET, SINT_min, PLNSOL, 0,MIN,
! /CONTOUR,,9,SINT_min,(1e9-SINT_min)/9,1e9
! PLNSOL,S,INT
! /CONTOUR,,9,
! PLNSOL,S,X
! PLNSOL,S,Y
! PLNSOL,S,Z

ESEL,S,MAT,,10
NSLE,S,1
PLNSOL,S,INT
*GET, SINT_min, PLNSOL, 0,MIN,
/CONTOUR,,9,SINT_min,(1e9-SINT_min)/9,1e9
PLNSOL,S,INT
/CONTOUR,,9,
PLNSOL,S,X
PLNSOL,S,Y
PLNSOL,S,Z
PLNSOL,U,X
PLNSOL,U,Y
PLNSOL,U,SUM

!!!!!!!!!!!!!!!!!!!!!!!!!!!
RSYS,SOLU
CSYS,0

ESEL,S,MAT,,3 
NSLE,S,1
PLNSOL,S,XY
PLNSOL,S,X

ESEL,S,MAT,,5
NSLE,S,1
PLNSOL,S,XY
PLNSOL,S,X 

ESEL,S,MAT,,7,8
NSLE,S,1
PLNSOL,S,XY
PLNSOL,S,X 

ESEL,S,MAT,,9
NSLE,S,1
PLNSOL,S,XY
PLNSOL,S,X 

!--------------------------------!
!!!	INSULATION FATIGUE INDEX LHD_xy = S_X/sigma0 + (S_XY/tau0)^2 (failure at 1)
!--------------------------------!

sigma0 = 38e6	! fatigue limit, normal stress [Pa]
tau0 = 27e6   ! fatigue limit, shear stress [Pa]

! INSULATION
ETABLE,ERAS
/AUTO

RSYS,SOLU
ESEL,S,MAT,,3 
NSLE,S,1
ETABLE,norm_stress,S,X
ETABLE,shear_stress_xy,S,XY
SEXP,shear_stress_xy2,shear_stress_xy,,2
SADD,LHD_xy,norm_stress,shear_stress_xy2,1/sigma0,1/(tau0**2)
/TITLE, Insulation - Coulomb/Mohr Sxy   
/GRAPH,FULL
PLETAB,LHD_xy,NOAV

ALLSEL
ESEL,S,MAT,,5 
ESEL,R,ELEM,,CSU1_E
NSLE,S,1
PLETAB,LHD_xy,NOAV
*GET,LHD_xy_MAX,PLNSOL,0,MAX
/CONTOUR,1,8,1,,LHD_xy_MAX
PLETAB,LHD_xy,NOAV

ALLSEL
ESEL,S,MAT,,7,8
ESEL,R,ELEM,,CSU1_E
NSLE,S,1
PLETAB,LHD_xy,NOAV
*GET,LHD_xy_MAX,PLNSOL,0,MAX
/CONTOUR,1,8,1,,LHD_xy_MAX
PLETAB,LHD_xy,NOAV


ALLSEL
ESEL,S,MAT,,9 
ESEL,R,ELEM,,CSU1_E
NSLE,S,1
PLETAB,LHD_xy,NOAV
*GET,LHD_xy_MAX,PLNSOL,0,MAX
/CONTOUR,1,8,1,,LHD_xy_MAX
PLETAB,LHD_xy,NOAV





