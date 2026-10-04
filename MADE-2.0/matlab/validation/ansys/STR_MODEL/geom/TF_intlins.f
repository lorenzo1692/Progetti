!=======================================================================
! TF_intlins.f - Inter-layer insulation, older version (see TF_lwins_trpz.f).
! Not called by the current STR_MODEL.f.
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTER-LAYER INSULATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = INTLINS_MAT
THIS_ET = INTLINS_ET
THAT_MAT = JCK_MAT
ELM_SZ = INTLT*2

VER_CS = MAX_CS

*DO,i,1,GRADES,2
	*IF,i,LE,NL_HIGHF,THEN
		k = 1
		CSYS,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+2
		KSEL,,LOC,Y,-HJH(i,k)-ICIT
		LSLK,,1
		CM,THIS_INTLINS,LINE 
		*GET,XREF,KP,,MNLOC,X
		
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX
		KSEL,NONE
		K,ULT_KP+1,XREF,-HJH(i,k)-ICIT
		K,ULT_KP+2,XREF,-HJH(i,k)-ICIT-INTLT*2
			LSEL,NONE
			L,ULT_KP+1,ULT_KP+2
			INTLINS_EXT_LN = LSNEXT(0)
			CMSEL,A,THIS_INTLINS,LINE
			NSLL,,1
				ASEL,NONE
				ADRAG,ALL,,,,,,INTLINS_EXT_LN
				KSLL,,1
				NSLL,,1
				NUMMRG,NODE
				NUMMRG,KP
				ALLSEL,BELOW,AREA
				AATT,THIS_MAT,,THIS_ET,TOP_CS
				MSHKEY,1
				MSHAPE,0
				ESIZE,ELM_SZ
				AMESH,ALL
	*ELSE
		k = 1
		CSYS,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+2
		KSEL,,LOC,Y,-HJH(i,k)-ICIT
		LSLK,,1
		CM,THIS_INTLINS,LINE 
		*GET,XREF,KP,,MNLOC,X
		
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX
		KSEL,NONE
		K,ULT_KP+1,XREF,-HJH(i,k)-ICIT
		K,ULT_KP+2,XREF,-HJH(i,k)-ICIT-INTLT*2
			LSEL,NONE
			L,ULT_KP+1,ULT_KP+2
			INTLINS_EXT_LN = LSNEXT(0)
			CMSEL,A,THIS_INTLINS,LINE
			NSLL,,1
				ASEL,NONE
				ADRAG,ALL,,,,,,INTLINS_EXT_LN
				KSLL,,1
				NSLL,,1
				NUMMRG,NODE
				NUMMRG,KP
				ALLSEL,BELOW,AREA
				AATT,THIS_MAT,,THIS_ET,TOP_CS
				MSHKEY,1
				MSHAPE,0
				ESIZE,ELM_SZ
				AMESH,ALL
	*ENDIF
*ENDDO

ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE
NUMMRG,KP


