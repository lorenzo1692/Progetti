!=======================================================================
! TF_gins_trpz_g10.f - Ground insulation around the stepped WP (convex hull of the layers,
! thickness GIT), meshed with the buffer material number.
! Called by STR_MODEL.f (WP_TOPOLOGY other than RIS_only).
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

THIS_MAT = BUFF_MAT
THIS_ET = BUFF_ET
THAT_MAT = JCK_MAT

CSYS,0
ASUM
*GET,REFX,AREA,0,CENT,X
*GET,REFY,AREA,0,CENT,Y
*GET, GI_CSYS ,CDSY, 0, NUM,MAX
LOCAL,GI_CSYS+1,0,REFX,REFY

LSEL,R,EXT
LSEL,U,MAT,,JCK_MAT
LSEL,U,MAT,,CAB_MAT

CM,LINE_GI_INT,LINE

!----------------------------------------------------!
n = 1
k = 1 
KSEL,ALL
*GET,ULT_KP,KP,,NUM,MAX

CSYS,CSYS_CABLE(n,k)+2
KSEL,S,LOC,Y,HJH(n,k)+ICIT  
LSLK,S,1
CM,THIS_INTLINSd,LINE
*GET,Yd,KP,,MNLOC,Y
*GET,Xs,KP,,MNLOC,X
*GET,Xd,KP,,MXLOC,X
K,ULT_KP+1,Xs-GIT,Yd+GIT
K,ULT_KP+2,Xd+GIT,Yd+GIT
	
n=0
nn=2
*DO,i,1,GRADES-1
	n = n+NL(i)	
	k=1+i
	
	CSYS,CSYS_CABLE(n,k)+2
	KSEL,,LOC,Y,-HJH(n,k)-ICIT 
	LSLK,,1
	*GET,Yd,KP,,MNLOC,Y
	*GET,Xs,KP,,MNLOC,X
	*GET,Xd,KP,,MXLOC,X
	K,ULT_KP+nn+1,Xs-GIT,Yd-GIT
	K,ULT_KP+nn+2,Xd+GIT,Yd-GIT

	nn = nn+2
*ENDDO

i=GRADES
n = n+NL(i)
k=1+i

CSYS,CSYS_CABLE(n,k)+2
KSEL,,LOC,Y,-HJH(n,k)-ICIT 
LSLK,,1
*GET,Yd,KP,,MNLOC,Y
*GET,Xs,KP,,MNLOC,X
*GET,Xd,KP,,MXLOC,X
K,ULT_KP+nn+1,Xs-GIT,Yd-GIT
K,ULT_KP+nn+2,Xd+GIT,Yd-GIT

LSEL,NONE
KSEL,ALL
L,ULT_KP+1,ULT_KP+2
L,ULT_KP+nn+1,ULT_KP+nn+2

*DO,i,1,(nn+2)/2-1
	k=2
	L,ULT_KP+1+k*(i-1),ULT_KP+1+k*i
	L,ULT_KP+2+k*(i-1),ULT_KP+2+k*i
*ENDDO

CM,LINE_GI_EXT,LINE

LSEL,S,LINE,,LINE_GI_INT
LSEL,A,LINE,,LINE_GI_EXT

AL,ALL
AATT,THIS_MAT,,THIS_ET, 
MSHKEY,0
MSHAPE,1
ESIZE,0.005
AMESH,ALL

LSEL,U,LINE,,LINE_GI_INT
LATT,BUFF_MAT
CM,TEMP_CONT,LINE
































