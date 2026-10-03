!=======================================================================
! Solution_STR.f - structural solution, checks and post-processing
!
! Called by RUN.dat after STR_MODEL.f. Re-opens STR_2D.db and applies:
! - symmetry of the wedged sector: UY = 0 (local cylindrical CS 2011) on
!   the outer face of the wedge insulation at both flanks;
! - vertical force of the inner leg, T_bf/2 (T_bf = vertical tension of
!   the whole coil, bending-free D-shape formula), applied through the
!   generalized plane strain (GSGDATA/GSBDATA);
! - Lorentz forces on the conductors, read from EM_2D.rst (LDREAD);
! - uniform temperature 4.2 K from the reference 293 K (cool-down).
! Writes CHECK_TF.txt [MN]: resultant of the Lorentz nodal forces Fx
! (~0 by symmetry) and Fy (centering force of the modelled sector), and
! Fz + T_bf/2 (= T_bf/2, no nodal FZ in 2D). Then solves (STR_2D.rst),
! plots (postpro\Plot.f) and writes the TFBM_*.txt diagnostic export
! (EXPORT_TF_BENCHMARK.mac).
!=======================================================================
FINISH
RESUME,STR_2D,db

ALLSEL,ALL,ALL
! optional: different contact stiffness for the conductor-jacket pairs
! STIFF = 10
! PENTR = 0.0
! PNBLL = 0.0
! RMODIF, CONT,3,STIFF,PENTR

! SAVE,PRE_SOL,db,,all                     ! optional: save before the solution

!------------------------------------!
!!! BOUNDARY CONDITIONS: UY = 0 on both flanks (cylindrical CS 2011)
!------------------------------------! 
/PREP7
ALLSEL
CSYS,2011
ALLS
ALLSEL,BELOW,AREA
*GET,MAX_Y,KP,,MXLOC,Y
*GET,MIN_Y,KP,,MNLOC,Y

ang = MAX_Y-MIN_Y                           ! sector angle seen by the model

CSYS,2011
ASEL,S,MAT,,9
ALLSEL,BELOW,AREA
NSEL,R,LOC,Y,0
NROTAT,ALL
D,ALL,UY,0
!*
ASEL,S,MAT,,9
ALLSEL,BELOW,AREA
NSEL,R,LOC,Y,-ang
NROTAT,ALL
D,ALL,UY,0

ALLSEL
CSYS,REF_CS 
!------------------------------------!
!!! CYCLIC SYMMETRY (alternative to the flank constraints, not used)
!------------------------------------!
! /PREP7
! ALLSEL,ALL,ALL
! CSYS,2011
! ASEL,S,MAT,,9
! ALLSEL,BELOW,AREA
! NSEL,R,LOC,Y,0
! CM,TF_M01H,NODE
! !*
! ALLSEL,BELOW,AREA
! NSEL,R,LOC,Y,-ang
! CM,TF_M01L,NODE

! CMSEL,S,TF_M01H,NODE
! CMSEL,A,TF_M01L,NODE
! CPCYC,ALL,,2011,,ang
! ALLSEL,ALL,ALL
! CYCLIC,16,,2011,TF,1
!------------------------------------!
!!! SOLUTION  
!------------------------------------!
! cool-down: uniform temperature set in /SOLU below (TUNIF,4.20)
/PREP7
ALLSEL,ALL,ALL
! BF,ALL,TEMP,4.2
! FINISH
!------------------------------------!
!!! VERTICAL LOAD: T_bf/2 on the inner leg, generalized plane strain
!------------------------------------!
CSYS,0
ASEL,,TYPE,,CAB_MAT,CASE_MAT-2
ASUM
*GET,REFX,AREA,0,CENT,X
*GET,REFY,AREA,0,CENT,Y
!*
k_bf = 0.5*log((OUTERLEG_INTR)/(INNLEG_OUTR))
Mu_0 = 4e-7*pi
T_bf = (k_bf*N_TF*(N_WIRES*TF_CURRENT)**2)*Mu_0/(2*PI)
! r0 = sqrt(Rm_OUTERLEG*Rm_INNERLEG)
! theta = PI/2
! dz = (r0*k_bf)*sin(theta)*exp(k_bf*sin(theta))*dtheta
ALLSEL
GSGDATA,SL_LENGTH,REFX,REFY,,
GSBDATA,F,T_bf/2,ROTX,0.0,ROTY,0.0
!------------------------------------!
!!! LORENTZ LOADS: nodal forces of the EM solution (EM_2D.rst)
!------------------------------------!
/PREP7
ALLSEL,ALL,ALL
ASEL,S,TYPE,,CAB_MAT
NSLA,S,1
LDREAD,FORC,,,,,EM_2D,rst

! check: resultant of the Lorentz forces and vertical load [MN] -> CHECK_TF.txt
*DIM,CHECK_TF,ARRAY,1,3
*GET, NODI_TF, NODE, 0, COUNT
NODO=0
Fx_TF = 0
Fz_TF = 0
Fy_TF = 0
*DO,h,1,NODI_TF
	NODO = NDNEXT(NODO)		
	*GET, temp_fx,  NODE, NODO,F,FX,
	*GET, temp_fy,  NODE, NODO,F,FY,
	*GET, temp_fz,  NODE, NODO,F,FZ,
	Fx_TF = Fx_TF + temp_fx
	Fy_TF = Fy_TF + temp_fy
	Fz_TF = Fz_TF + temp_fz
*ENDDO
CHECK_TF(1,1) = Fx_TF*1e-6 
CHECK_TF(1,2) = Fy_TF*1e-6 
CHECK_TF(1,3) = Fz_TF*1e-6 + T_bf/2*1e-6
!------------------------------------
*MWRITE,CHECK_TF,CHECK_TF,txt
%G %G %G 
!------------------------------------	
ALLSEL,ALL,ALL
FINISH
!------------------------------------
/SOLU
ANTYPE,0,
CSYS,REF_CS
ALLSEL,ALL,ALL
TUNIF,4.20
KBC,0
PRED,ON
AUTOTS,AUTO
SOLVE
FINISH
!------------------------------------!
/POST1
/INPUT,'%MDIR%postpro\Plot','f'			! plots (PNG files in the working directory)
!------------------------------------!
FINISH
RESUME,STR_2D,db                            ! export from the solved model
/INPUT,'%MDIR%EXPORT_TF_BENCHMARK','mac'
FINISH