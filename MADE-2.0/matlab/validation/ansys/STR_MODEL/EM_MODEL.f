!=======================================================================
! EM_MODEL.f - geometry and mesh of the 2D electromagnetic model
!
! Called by RUN.dat (after /FILNAME,EM_2D). Reads the design parameters
! (input\ParametriTF.f, written by MADE export_ansys_input) and the
! materials (input\Materiali.f), builds the conductors of every turn and
! the surrounding air, and saves EM_2D.db. Only the conductors carry
! current: jacket, insulation and case are not modelled here.
!=======================================================================
FINISH
/PREP7
! Working directory: chosen by RUN.dat (the .\output folder). The former
! hard-coded /CWD to the network folder was removed, so the model runs
! wherever it is launched; model files are read through MDIR.
*GET,MDIR_T,PARM,MDIR,TYPE                  ! MDIR = model folder, set by RUN.dat ('..\' from .\output)
*IF,MDIR_T,EQ,-1,THEN                       ! not defined: run by hand from the model folder
	MDIR = '.\'
*ENDIF

/PREP7
/INPUT,'%MDIR%input\ParametriTF','f' 
/INPUT,'%MDIR%input\Materiali','f' 

! Local coordinate systems
LOCAL,2001,0,,,,ANG_DIV/2 				! Cartesian local CS used for duplicating half-WP section
LOCAL,2011,1,,,,ANG_DIV/2+90 			! Cylindrical local CS used for assigning boundary conditions
LOCAL,2002,0,,,,ANG_DIV					! Cartesian local CS to assign conditions to second WP
REF_CS = 0                                 ! reference (global Cartesian) CS
REF_CS_SYMM = 0
SYMM_PLANE = 0

! material and element-type numbers of the EM model (properties in input\Materiali.f)
CAB_MAT    = 1
CAB_ET     = 1
JCK_MAT    = 2
JCK_ET     = 2
TINS_MAT   = 3
TINS_ET    = 3
FILL_MAT   = 4
FILL_ET    = 4
INTLINS_MAT = 5
INTLINS_ET  = 5
DLINS_MAT  = 6
DLINS_ET   = 6
GINS_MAT   = 7
GINS_ET    = 7
BUFF_MAT   = 8
BUFF_ET    = 8
WINS_MAT   = 9
WINS_ET    = 9
CASE_MAT   = 10
CASE_ET    = 10
AIR_MAT    = 20
AIR_ET     = 20

*IF,WP_TOPOLOGY,EQ,'RIS_only',THEN
	/INPUT,'%MDIR%geom\TF_cab_RIS','f',
*ELSE
	! /INPUT,'%MDIR%geom\TF_cab','f',			! alternative: equal-width layers (not used)
	/INPUT,'%MDIR%geom\TF_cab_trpz','f',
*ENDIF
/INPUT,'%MDIR%geom\TF_air','f',			! air region around the WP up to AIR_HEIGHT

/PREP7
ALLSEL,ALL,ALL

FINISH
SAVE,EM_2D,db,,all

