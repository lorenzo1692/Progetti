!=======================================================================
! Solution_EM.f - electromagnetic solution and plots
!
! Called by RUN.dat after EM_MODEL.f. Re-opens EM_2D.db, switches the
! MESH200 elements to PLANE233 (conductors and air), sets the flux
! boundary condition and the current of every turn (TF_CURRENT, coupled
! VOLT DOFs per turn), solves (EM_2D.rst), writes TFBM_em_plane233.txt
! (EXPORT_TF_EM_PLANE233.mac) and the field plots (PNG). The Lorentz
! forces in EM_2D.rst are read by Solution_STR.f (LDREAD).
!=======================================================================
RESUME,EM_2D,db

/PREP7

! element types: PLANE233 (2D electromagnetic) for conductors (1) and air (20),
! null elements for the other materials (not present in the EM model)
ET,1,233,1
ET,2,0
ET,3,0
ET,4,0
ET,5,0
ET,6,0
ET,7,0
ET,8,0
ET,9,0
ET,10,0
ET,20,233,1

KEYOPT,1,8,1                                ! forces computed as Lorentz forces J x B (read by LDREAD)
KEYOPT,1,3,0

KEYOPT,20,8,1
KEYOPT,20,3,0

ALLSEL,ALL,ALL

!-----------------------------------------!
! boundary condition: flux parallel to the outer air boundary (AZ = 0 at r = AIR_HEIGHT)
CSYS,1
KSEL,S,LOC,X,AIR_HEIGHT
LSLK,S,1
NSLL,S,1
D,ALL,AZ,0
ALLSEL

R,THIS_REAL,SL_LENGTH                       ! depth SL_LENGTH (no N_TF symmetry factor)
!-----------------------------------------!
! current load: one real constant per turn, VOLT coupled, current TF_CURRENT
ALLSEL
THIS_REAL = 100
n=0
*DO,i,1,GRADES	
	*DO,j,1,NL(i)	
		*DO,k,1,NT(i)
			n=n+1	
			THIS_REAL = THIS_REAL+1			
			CSYS,REF_CS
			ASEL,S,REAL,,100+n
			R,THIS_REAL,SL_LENGTH,!,N_TF,,,,,,1/N_TF
			ALLSEL,BELOW,AREA
			CP,NEXT,VOLT,ALL
			F,NDNEXT(0),AMPS,TF_CURRENT
			ALLSEL  
		*ENDDO
	*ENDDO
*ENDDO
FINISH
!----------------------------!
/SOLU
ALLSEL,ALL,ALL
SOLVE
FINISH
!----------------------------!

! /INPUT,'%MDIR%EXPORT_TF_EM_BENCHMARK','mac'        ! alternative conductor-only export
/INPUT,'%MDIR%EXPORT_TF_EM_PLANE233','mac'

!---- plots (PNG files in the working directory) ----!
/POST1
SET,LAST
ALLSEL,ALL,ALL

PLNSOL,B,SUM
!*  
PLVECT,B, , , ,VECT,ELEM,ON,0   
 
EMFT

! /DEV,FONT,LEGEND,MENU 
/dev,font,1,Calibri,500,0,72,0,0,,, 
/dev,font,2,Calibri,500,0,72,0,0,,, 
/dev,font,3,Calibri,500,0,72,0,0,,, 
/TSPEC, 4, 24 , 2, 
 
/RGB,INDEX,100,100,100, 0   
/RGB,INDEX, 80, 80, 80,13   
/RGB,INDEX, 60, 60, 60,14   
/RGB,INDEX, 0, 0, 0,15  

/GFILE,2400 
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
PLNSOL,B,SUM

ESEL,U,MAT,,20
NSLE,S,1
PLNSOL,B,SUM

! /CONT,,15,0,,20.3
! /REP
! ALLSEL
! /REP
 
/CONT,,15,0,,12.5
/REP

ALLSEL

! |B| along the radial path from the axis to AIR_HEIGHT
CSYS,0
PATH,FLUX,2,30,50,  
PPATH,1,,0,0,0,0, 
PPATH,2,,0,AIR_HEIGHT,0,0,   
PDEF,B__,B,SUM	
PLPATH,B__




