!=======================================================================
! STR_MODEL.f - geometry, mesh and contacts of the 2D structural model
!
! Called by RUN.dat (after /FILNAME,STR_2D). Reads the design parameters
! (input\ParametriTF.f, written by MADE export_ansys_input) and the
! materials (input\Materiali.f), builds turn by turn the conductors,
! jackets, turn insulation and corner fillers, then the inter-layer and
! ground insulation, the case and the wedge insulation, defines the
! contacts (BC\Contacts.f) and saves STR_2D.db.
!=======================================================================

FINISH
/PREP7
SHPP,OFF                                    ! no shape checks during meshing (thin insulation layers)
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
LOCAL,2003,0,,,,90

REF_CS = 0
REF_CS_SYMM = 0
SYMM_PLANE = 0

! material and element-type numbers (properties in input\Materiali.f)
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
STEEL_FIL_MAT = 11
STEEL_FIL_ET = 10
AIR_MAT    = 20
AIR_ET	   = 20
! Winding pack. WP_TOPOLOGY (input\ParametriTF.f): 'RIS_only' = round-in-square
! conductors; otherwise rectangular conductors, layers of different width
! (trapezoidal/stepped WP, *_trpz macros)
*IF,WP_TOPOLOGY,EQ,'RIS_only',THEN
	/INPUT,'%MDIR%geom\TF_cab_RIS','f', 
	/INPUT,'%MDIR%geom\TF_jck_RIS','f', 
	/INPUT,'%MDIR%geom\TF_tins_RIS','f', 
	/INPUT,'%MDIR%geom\TF_fill_RIS','f', 
	! /INPUT,'%MDIR%geom\TF_fill_TRP','f', 	! alternative (not used)
*ELSE
	! /INPUT,'%MDIR%geom\TF_cab','f',		! alternative: equal-width layers (not used) 
	/INPUT,'%MDIR%geom\TF_cab_trpz','f',
	! /INPUT,'%MDIR%geom\TF_jck','f',		! alternative: equal-width layers (not used)
	/INPUT,'%MDIR%geom\TF_jck_trpz','f',
	! /INPUT,'%MDIR%geom\TF_tins','f',		! alternative: equal-width layers (not used)
	/INPUT,'%MDIR%geom\TF_tins_trpz','f',
	! /INPUT,'%MDIR%geom\TF_fill','f',		! alternative: equal-width layers (not used)
	/INPUT,'%MDIR%geom\TF_fill_trpz','f',
*ENDIF
ALLSEL,ALL,ALL
*IF,WP_TOPOLOGY,EQ,'RIS_only',THEN
	/INPUT,'%MDIR%geom\TF_lwins','f',
	! /INPUT,'%MDIR%geom\TF_cfill','f',
	/INPUT,'%MDIR%geom\TF_gins_panck','f',	! NOTE: TF_gins_panck.f is missing from geom\: restore it before using 'RIS_only'
*ELSE
	/INPUT,'%MDIR%geom\TF_lwins_trpz','f',
	
	! /INPUT,'%MDIR%geom\TF_gins_trpz','f',	! alternative ground insulation (not used)
	/INPUT,'%MDIR%geom\TF_gins_trpz_g10','f',	! ground insulation around the stepped WP
	
	! /INPUT,'%MDIR%geom\TF_buff','f',		! buffer insulation (not used)
	! /INPUT,'%MDIR%geom\TF_intlins','f',	! inter-layer insulation, old version (not used)
	! /INPUT,'%MDIR%geom\TF_dlins','f',		! double-layer insulation (not used)
	! external boundary lines of the WP (jackets and steel fillers excluded):
	! line attribute BUFF_MAT, the side of the WP / ground insulation contact
	ASEL,,MAT,,JCK_MAT,BUFF_MAT
	ALLSEL,BELOW,AREA
	LSEL,R,EXT
	LSEL,U,MAT,,JCK_MAT
	ASEL,S,MAT,,STEEL_FIL_MAT
	LSLA,U
	LATT,BUFF_MAT
*ENDIF
ALLSEL,ALL,ALL
! SAVE,WP,db                                 ! optional: save the WP alone
/INPUT,'%MDIR%geom\TF_case','f',			! case (arc bore at WEDGE_INTR_INS, flat plasma side)
/INPUT,'%MDIR%geom\TF_wins','f',			! wedge insulation on the case flanks
/PREP7
LSEL,,LCCAT
LDELE,ALL
FINISH

!------------------------------------!
/INPUT,'%MDIR%BC\Contacts','f'				! structural elements and contact pairs
!------------------------------------!

SAVE,STR_2D,db,,all
