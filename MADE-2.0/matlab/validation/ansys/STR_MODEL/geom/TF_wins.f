!=======================================================================
! TF_wins.f - Wedge insulation (thickness CASE_INS) on the case flanks, mirrored to
! both sides, with its material axes oriented normal to the flank.
! Called by STR_MODEL.f.
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! WEDGE INSULATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = WINS_MAT
THIS_ET = WINS_ET
THIS_REAL = WINS_MAT
ELM_SZ = CASE_INS

ASEL,,MAT,,CASE_MAT
ALLSEL,BELOW,AREA
CSYS,REF_CS
*GET,REF_Xmax,KP,,MXLOC,X
*GET,REF_Xmin,KP,,MNLOC,X
*GET,REF_Ymax,KP,,MXLOC,Y
*GET,REF_Ymin,KP,,MNLOC,Y

CSYS,2011
KSEL,ALL
*GET,ULT_KP,KP,0,NUM,MAX
KSEL,NONE
 
! K,ULT_KP+1,WEDGE_OUTR,THETA_INS-ANG_DIV 
! K,ULT_KP+2,WEDGE_OUTR_INS,-ANG_DIV,

CSYS,2003
K,ULT_KP+1,INNLEG_OUTR,-CASE_WIDTH/2
K,ULT_KP+2,INNLEG_OUTR,-CASEW_UP/2

LSEL,NONE
L,ULT_KP+1,ULT_KP+2
CM,WINS_EXT_LN,LINE

CSYS,REF_CS
ASEL,,MAT,,CASE_MAT
ALLSEL,BELOW,AREA
LSEL,R,EXT
LSEL,U,LOC,Y,REF_Ymax
LSEL,U,LOC,Y,REF_Ymin
LSEL,U,MAT,,CASE_MAT
LSEL,U,LOC,X,-1e3,0

CM,CASE_SIDE_FACE,LINE
LGEN,2,ALL
CMSEL,U,CASE_SIDE_FACE,LINE

ASEL,NONE
ADRAG,ALL,,,,,,WINS_EXT_LN
AATT,THIS_MAT,,THIS_ET,WEDGE_CS
ESIZE,ELM_SZ
MSHKEY,1
MSHAPE,0
AMESH,ALL

CSYS,2003
KSEL,ALL
*GET,ULT_KP,KP,0,NUM,MAX
K,ULT_KP+3,INNLEG_OUTR,CASE_WIDTH/2
K,ULT_KP+4,INNLEG_OUTR,CASEW_UP/2

LSEL,NONE
L,ULT_KP+3,ULT_KP+4
CM,WINS_EXT_LN,LINE

CSYS,REF_CS
ASEL,,MAT,,CASE_MAT
ALLSEL,BELOW,AREA
LSEL,R,EXT
LSEL,U,LOC,Y,REF_Ymax
LSEL,U,LOC,Y,REF_Ymin
LSEL,U,MAT,,CASE_MAT
LSEL,U,LOC,X,0,1e3

CM,CASE_SIDE_FACE,LINE
LGEN,2,ALL
CMSEL,U,CASE_SIDE_FACE,LINE

ASEL,NONE
ADRAG,ALL,,,,,,WINS_EXT_LN
AATT,THIS_MAT,,THIS_ET,WEDGE_CS
ESIZE,ELM_SZ
MSHKEY,1
MSHAPE,0
AMESH,ALL

!------------------------------------------------------!
! !!! Mirror the area	
! CSYS,REF_CS
! ASEL,,MAT,,THIS_MAT,,,1
! ALLSEL,BELOW,AREA
! ARSYM,X,ALL
!------------------------------------------------------!
!!! Orient the material axes
CSYS,REF_CS_SYMM
*GET,MAX_CS,CDSY,,NUM,MAX
WEDGE_CS_SYMMd = MAX_CS+1
LOCAL,WEDGE_CS_SYMMd,0,WEDGE_OUTR_INS*sin(ANG_DIV/2*pi/180),WEDGE_OUTR_INS*cos(ANG_DIV/2*pi/180),,-INCL
CSYS,REF_CS_SYMM

ASEL,,MAT,,THIS_MAT,,,1
ALLSEL,BELOW,AREA
ESLA,,1
ESEL,R,CENT,X,0,1e4! INNLEG_OUTR*sin(ANG_DIV/2*pi/180)
EMODIF,ALL,ESYS,WEDGE_CS_SYMMd

CSYS,REF_CS_SYMM
*GET,MAX_CS,CDSY,,NUM,MAX
WEDGE_CS_SYMMs = MAX_CS+1
LOCAL,WEDGE_CS_SYMMs,0,WEDGE_OUTR_INS*sin(-ANG_DIV/2*pi/180),WEDGE_OUTR_INS*cos(-ANG_DIV/2*pi/180),,INCL
CSYS,REF_CS_SYMM

ASEL,,MAT,,THIS_MAT,,,1
ALLSEL,BELOW,AREA
ESLA,,1
ESEL,R,CENT,X,-1e4,0!INNLEG_OUTR*sin(-ANG_DIV/2*pi/180),0,
EMODIF,ALL,ESYS,WEDGE_CS_SYMMs

ASEL,,MAT,,THIS_MAT,,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE
NUMMRG,KP
