!=======================================================================
! TF_case.f - Case of one coil: flat plasma side at INNLEG_OUTR, flanks at
! +-ANG_DIV/2, bore arc at WEDGE_INTR_INS = Rk_/cos(ANG_DIV/2), cavity
! around the WP. Called by STR_MODEL.f.
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! CASING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = CASE_MAT
THIS_ET = CASE_ET
THIS_REAL = CASE_MAT
ELM_SZ = CASE_THICK*(1-CASE_FLAT)/3

! ASEL,,MAT,,BUFF_MAT
ALLSEL,BELOW,AREA
CSYS,REF_CS
*GET,REF_Xmax,KP,,MXLOC,X
*GET,REF_Xmin,KP,,MNLOC,X
*GET,REF_Ymax,KP,,MXLOC,Y
*GET,REF_Ymin,KP,,MNLOC,Y

!!! TOP !!!
CSYS,2003
KSEL,ALL
*GET,ULT_KP,KP,0,NUM,MAX
KSEL,NONE
K,ULT_KP+1,INNLEG_OUTR,-CASE_WIDTH/2
K,ULT_KP+2,INNLEG_OUTR,CASE_WIDTH/2 
K,ULT_KP+3,INNLEG_INTR,-CASE_WIDTH_LOW/2 
K,ULT_KP+4,INNLEG_INTR,CASE_WIDTH_LOW/2 
K,ULT_KP+5,

CSYS,0
K,ULT_KP+5,REF_Xmin*4,REF_Ymin
K,ULT_KP+6,REF_Xmax*4,REF_Ymin
					
K,ULT_KP+7,REF_Xmin*4,REF_Ymax
K,ULT_KP+8,REF_Xmax*4,REF_Ymax

LSEL,NONE
L,ULT_KP+1,ULT_KP+2
CM,LN1,LINE
LSEL,NONE
L,ULT_KP+1,ULT_KP+3
CM,LN2,LINE
LSEL,NONE
L,ULT_KP+4,ULT_KP+2
CM,LN3,LINE
LSEL,NONE
! L,ULT_KP+3,ULT_KP+4
LARC,ULT_KP+3,ULT_KP+4,ULT_KP+5,-Rk_
CM,LN4,LINE		   
		   
LSEL,S,MAT,,BUFF_MAT
LSEL,R,EXT
CM,TEMP1,LINE
LGEN,2,ALL,,,,,,,1
CMSEL,U,TEMP1,LINE
LATT,CASE_MAT
CM,CASEINT,LINE
		   
CMSEL,A,LN1,LINE
CMSEL,A,LN2,LINE
CMSEL,A,LN3,LINE	   
CMSEL,A,LN4,LINE	  

KSLL,S
NUMMRG,KP

ASEL,NONE
AL,ALL 
CM,CASETEMP,AREA

KSEL,ALL

LSEL,NONE
L,ULT_KP+5,ULT_KP+6
L,ULT_KP+7,ULT_KP+8
CM,LTEMP,LINE
ASBL,CASETEMP,LTEMP
CM,CASE,AREA
ALLSEL,BELOW,AREA

LSEL,R,LOC,Y,REF_Ymin
LCCAT,ALL

ALLSEL,BELOW,AREA
LSEL,R,LOC,Y,REF_Ymax
LCCAT,ALL

ASEL,,AREA,,CASE
ALLSEL,BELOW,AREA
AATT,THIS_MAT,,THIS_ET
MSHKEY,0
MSHAPE,0
SMRTSIZE,OFF
ESIZE,0.01
AMESH,ALL

ASEL,,MAT,,THIS_MAT,,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE,1e-5
NUMMRG,KP,1e-5

LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

ALLSEL,ALL,ALL