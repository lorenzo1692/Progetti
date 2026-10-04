!=======================================================================
! TF_air.f - Air region of the electromagnetic model, up to AIR_HEIGHT.
! Called by EM_MODEL.f.
!=======================================================================
/PREP7

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! AIR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = AIR_MAT
THAT_MAT = CAB_MAT
THIS_ET = AIR_ET
ELM_SZ = 0.01!4*(HCH_N-CFR)/2

THIS_REAL = THIS_REAL+1

CSYS,REF_CS
KSEL,ALL
*GET,ULT_KP,KP,,NUM,MAX
K,ULT_KP+1,,,
K,ULT_KP+2,,AIR_HEIGHT
K,ULT_KP+3,-AIR_HEIGHT*sin(INCL*pi/180),AIR_HEIGHT*cos(INCL*pi/180)

LSEL,NONE
L,ULT_KP+1,ULT_KP+2
LARC,ULT_KP+2,ULT_KP+3,ULT_KP+1,AIR_HEIGHT
L,ULT_KP+3,ULT_KP+1
ASEL,NONE
AL,ALL
CM,AIR_AREA,AREA
ASEL,,MAT,,CAB_MAT,,,1
ALLSEL,BELOW,AREA
CM,CAB_AREA,AREA
AGEN,2,ALL,,,,,,,1
CMSEL,U,CAB_AREA,AREA
CM,CAB_AREA,AREA
CMSEL,A,AIR_AREA,AREA
ASBA,AIR_AREA,CAB_AREA
AATT,THIS_MAT,THIS_REAL,THIS_ET
ALLSEL,ALL,ALL
NUMMRG,NODE
NUMMRG,KP
ASEL,S,MAT,,THIS_MAT,,,1
SMRTSIZE,5
ESIZE,ELM_SZ!/2
MSHKEY,0
MSHAPE,0	
AMESH,ALL
		
CSYS,SYMM_PLANE
ASEL,,MAT,,THIS_MAT,,,1
ALLSEL,BELOW,AREA
ARSYM,X,ALL

ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
ALLSEL,BELOW,AREA
! Tight tolerance (1 micron): only the truly coincident nodes are merged (air
! on the symmetry axis, mirrored air on the cables of the other half). The
! default 1e-4 also merged two distinct nodes of a small air triangle (< 0.1 mm)
! and collapsed it (zero Jacobian -> run stopped). Same tolerance as TF_fill_trpz.
NUMMRG,NODE,1e-6
NUMMRG,KP,1e-6
