!=======================================================================
! TF_lwins_trpz.f - Inter-layer insulation of the layer-wound WP (thickness LWIT) between
! consecutive layers. Stepped WP.
! Called by STR_MODEL.f (WP_TOPOLOGY other than RIS_only).
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTER-LAYER INSULATION (layer-wound WP, thickness LWIT)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = INTLINS_MAT
THIS_ET = INTLINS_ET
THAT_MAT = JCK_MAT
ELM_SZ = LWIT!*2

SELTOL,1e-4
THIS_REAL = 100
n=0
*DO,i,1,GRADES	
	*DO,j,1,NL(i)-1
		n = n+1
		k=1
		
		CSYS,CSYS_CABLE(n,k)+2
		KSEL,,LOC,Y,-HJH(n,k)-ICIT
		LSLK,,1
		CM,THIS_INTLINSu,LINE
		*GET,Xs,KP,,MNLOC,X
		*GET,Xd,KP,,MXLOC,X
		*GET,Yu,KP,,MNLOC,Y
		KSEL,ALL
		KSEL,S,KP,,KP(Xs,Yu,0)
		*GET,KP1s,KP,,NUM,MAX
		KSEL,ALL
		KSEL,S,KP,,KP(Xd,Yu,0)
		*GET,KP1d,KP,,NUM,MAX		
		
		KSEL,,LOC,Y,-HJH(n,k)-ICIT-LWIT 
		LSLK,,1
		CM,THIS_INTLINSd,LINE
		*GET,Yd,KP,,MNLOC,Y
		KSEL,ALL
		KSEL,S,KP,,KP(Xs,Yd,0)
		*GET,KP2s,KP,,NUM,MAX
		KSEL,ALL
		KSEL,S,KP,,KP(Xd,Yd,0)
		*GET,KP2d,KP,,NUM,MAX
		
		ALLSEL,ALL,ALL
		
		LSEL,NONE
		L,KP1s,KP2s
		L,KP1d,KP2d
		CMSEL,A,THIS_INTLINSu,LINE
		CMSEL,A,THIS_INTLINSd,LINE
		ASEL,NONE
		AL,ALL
		ALLSEL,BELOW,AREA
		NUMMRG,NODE
		NUMMRG,KP
		CMSEL,S,THIS_INTLINSu,LINE
		LCCAT,ALL
		CMSEL,S,THIS_INTLINSd,LINE
		LCCAT,ALL
		AATT,THIS_MAT,,THIS_ET,TOP_CS
		MSHKEY,0
		MSHAPE,1
		ESIZE,ELM_SZ
		AMESH,ALL
		
		LSEL,,LCCAT
		LDELE,ALL		
	*ENDDO
	n = n+1
*ENDDO

*ENDDO

n=0
*DO,i,1,GRADES-1
	n = n+NL(i)
	k=1+i
	
	CSYS,CSYS_CABLE(n,k)+2
	KSEL,,LOC,Y,-HJH(n,k)-ICIT-LWIT 
	LSLK,,1
	CM,THIS_INTLINSd,LINE
	*GET,Yd,KP,,MNLOC,Y
	*GET,Xs,KP,,MNLOC,X
	*GET,Xd,KP,,MXLOC,X
	KSEL,ALL
	KSEL,S,KP,,KP(Xs,Yd,0)
	*GET,KP2s,KP,,NUM,MAX
	KSEL,ALL
	KSEL,S,KP,,KP(Xd,Yd,0)
	*GET,KP2d,KP,,NUM,MAX		
	
	KSEL,,LOC,Y,-HJH(n,k)-ICIT
	LSLK,,1
	LSEL,R,LOC,X,Xs,Xd
	KSLL,S,1
	CM,THIS_INTLINSu,LINE
	*GET,Yu,KP,,MNLOC,Y
	KSEL,ALL
	KSEL,S,KP,,KP(Xs,Yu,0)
	*GET,KP1s,KP,,NUM,MAX
	KSEL,ALL
	KSEL,S,KP,,KP(Xd,Yu,0)
	*GET,KP1d,KP,,NUM,MAX		

	ALLSEL,ALL,ALL
	
	LSEL,NONE
	L,KP1s,KP2s
	L,KP1d,KP2d
	CMSEL,A,THIS_INTLINSu,LINE
	CMSEL,A,THIS_INTLINSd,LINE
	ASEL,NONE
	AL,ALL
	ALLSEL,BELOW,AREA
	NUMMRG,NODE
	NUMMRG,KP
	CMSEL,S,THIS_INTLINSu,LINE
	LCCAT,ALL
	CMSEL,S,THIS_INTLINSd,LINE
	LCCAT,ALL
	AATT,THIS_MAT,,THIS_ET,TOP_CS
	MSHKEY,0
	MSHAPE,1
	ESIZE,ELM_SZ
	AMESH,ALL
	
	LSEL,,LCCAT
	LDELE,ALL	
*ENDDO

ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE,1e-5
NUMMRG,KP,1e-5


