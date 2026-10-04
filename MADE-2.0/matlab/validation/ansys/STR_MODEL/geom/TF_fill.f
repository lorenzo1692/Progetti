!=======================================================================
! TF_fill.f - Corner fillers of a WP with layers of equal width (older version of
! TF_fill_trpz.f). Not called by the current STR_MODEL.f.
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

SHPP,OFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FILLER MATERIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_COMP = 'FILL'
THIS_MAT = FILL_MAT
THAT_MAT = JCK_MAT
THIS_ET = FILL_ET

*DO,i,1,GRADES
	*IF,i,LE,NL_HIGHF,THEN
		*DO,k,1,NT_HIGHF
		KSEL,ALL
		LSEL,ALL
		NUMCMP,KP
		NUMCMP,LINE
		*GET,ULT_KP,KP,,NUM,MAX				 
		! Geometric parameters
		THICK = ICIT
		! R(i,k) = JT_h(i,k)/2			!Jacket internal radius
		R_N = R(i,k)	
		IN_RAD = R_N+JT_h(i,k)       	!Jacket external radius		
		EXT_RAD = IN_RAD+THICK
		
		JH_N = JH(i,k)			!Half Jacket Height
		JW_N = 2*HJW(i,k)       !Half jacket width
		HCW_N = HCW(i,k)        !Half Conductor width
		HCH_N = HCH(i,k)        !Half Conductor Height
		JT_Nh = JT_h(i,k)
		JT_Nw = JT_w(i,k)

		LHI_N = LHI(i,k)
		LWI_N = LWI(i,k)

		! Local CS to identify cable centre and later define orthotropic material properties
		CSYS,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+2

		THIS_REAL = 100+NT_HIGHF*(i-1)+k
		THIS_CS = MAX_CS+8*((k-1)+NT_HIGHF*(i-1))

     	! Create filler material areas
		!TP-RGT
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX	
		K,ULT_KP+1,HCW_N-R_N,HCH_N+JT_Nh+THICK
		K,ULT_KP+2,HCW_N+JT_Nw+THICK,HCH_N+JT_Nh+THICK
		K,ULT_KP+3,HCW_N+JT_Nw+THICK,HCH_N+JT_Nh-IN_RAD
		K,ULT_KP+4,HCW_N-R_N,HCH_N+JT_Nh-IN_RAD
		LSEL,NONE
		L,ULT_KP+1,ULT_KP+2
		L,ULT_KP+3,ULT_KP+2
		LARC,ULT_KP+1,ULT_KP+3,ULT_KP+4,EXT_RAD
		ASEL,NONE
		AL,ALL
		AATT,THIS_MAT,THIS_REAL,THIS_ET
		!BT_RGT
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX	
		K,ULT_KP+1,HCW_N-R_N,-HCH_N-JT_Nh-THICK
		K,ULT_KP+2,HCW_N+JT_Nw+THICK,-HCH_N-JT_Nh-THICK
		K,ULT_KP+3,HCW_N+JT_Nw+THICK,-HCH_N-JT_Nh+IN_RAD
		K,ULT_KP+4,HCW_N-R_N,-HCH_N-JT_Nh+IN_RAD
		LSEL,NONE
		L,ULT_KP+1,ULT_KP+2
		L,ULT_KP+3,ULT_KP+2
		LARC,ULT_KP+1,ULT_KP+3,ULT_KP+4,EXT_RAD
		ASEL,NONE
		AL,ALL
		AATT,THIS_MAT,THIS_REAL,THIS_ET

		LSEL,,LCCAT
		LDELE,ALL

		ASEL,S,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT-1,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		ALLSEL,BELOW,AREA
		! NUMMRG,NODE
		! NUMMRG,KP
		
		LSEL,R,RADIUS,,EXT_RAD
		LSEL,R,LINE,,LSNEXT(0)
		*GET,NUMDIVL,LINE,LSNEXT(0),ATTR,NDIV
		
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ALLSEL,BELOW,AREA
		ALLSEL,BELOW,AREA
		LSEL,U,RADIUS,,EXT_RAD
		LESIZE,ALL,,,NUMDIVL/2,,1
	
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ALLSEL,BELOW,AREA
		ESIZE,,NUMDIVL/2
	
		MSHKEY,0
		MSHAPE,1
		! AMESH,ALL		
		ARSYM,X,ALL
		*ENDDO
	*ELSE
		*DO,k,1,NT_LOWF
		KSEL,ALL
		LSEL,ALL
		NUMCMP,KP
		NUMCMP,LINE
		*GET,ULT_KP,KP,,NUM,MAX	
			
		! Geometric parameters
		THICK = ICIT
		! R(i,k) = JT_h(i,k)/2			!Jacket internal radius
		R_N = R(i,k)	
		IN_RAD = R_N+JT_h(i,k)       	!Jacket external radius		
		EXT_RAD = IN_RAD+THICK
		
		JH_N = JH(i,k)			!Half Jacket Height
		JW_N = 2*HJW(i,k)       !Half jacket width
		HCW_N = HCW(i,k)        !Half Conductor width
		HCH_N = HCH(i,k)        !Half Conductor Height
		JT_Nh = JT_h(i,k)
		JT_Nw = JT_w(i,k)

		LHI_N = LHI(i,k)
		LWI_N = LWI(i,k)

		! Local CS to identify cable centre and later define orthotropic material properties
		CSYS,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+2

		THIS_REAL = 100+NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+k
		THIS_CS = MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))

		! Create filler material areas
		!TP-RGT
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX	
		K,ULT_KP+1,HCW_N-R_N,HCH_N+JT_Nh+THICK
		K,ULT_KP+2,HCW_N+JT_Nw+THICK,HCH_N+JT_Nh+THICK
		K,ULT_KP+3,HCW_N+JT_Nw+THICK,HCH_N+JT_Nh-IN_RAD
		K,ULT_KP+4,HCW_N-R_N,HCH_N+JT_Nh-IN_RAD
		LSEL,NONE
		L,ULT_KP+1,ULT_KP+2
		L,ULT_KP+3,ULT_KP+2
		LARC,ULT_KP+1,ULT_KP+3,ULT_KP+4,EXT_RAD
		ASEL,NONE
		AL,ALL
		AATT,THIS_MAT,THIS_REAL,THIS_ET
		!BT_RGT
		KSEL,ALL
		*GET,ULT_KP,KP,,NUM,MAX	
		K,ULT_KP+1,HCW_N-R_N,-HCH_N-JT_Nh-THICK
		K,ULT_KP+2,HCW_N+JT_Nw+THICK,-HCH_N-JT_Nh-THICK
		K,ULT_KP+3,HCW_N+JT_Nw+THICK,-HCH_N-JT_Nh+IN_RAD
		K,ULT_KP+4,HCW_N-R_N,-HCH_N-JT_Nh+IN_RAD
		LSEL,NONE
		L,ULT_KP+1,ULT_KP+2
		L,ULT_KP+3,ULT_KP+2
		LARC,ULT_KP+1,ULT_KP+3,ULT_KP+4,EXT_RAD
		ASEL,NONE
		AL,ALL
		AATT,THIS_MAT,THIS_REAL,THIS_ET

		LSEL,,LCCAT
		LDELE,ALL

		ASEL,S,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT-1,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		ALLSEL,BELOW,AREA
		! NUMMRG,NODE
		! NUMMRG,KP
		
		LSEL,R,RADIUS,,EXT_RAD
		LSEL,R,LINE,,LSNEXT(0)
		*GET,NUMDIVL,LINE,LSNEXT(0),ATTR,NDIV
		
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ALLSEL,BELOW,AREA
		ALLSEL,BELOW,AREA
		LSEL,U,RADIUS,,EXT_RAD
		LESIZE,ALL,,,NUMDIVL/2,,1
	
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ALLSEL,BELOW,AREA
		ESIZE,,NUMDIVL/2
	
		MSHKEY,0
		MSHAPE,1
		! AMESH,ALL		
		ARSYM,X,ALL
		*ENDDO
	*ENDIF
*ENDDO

ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE,1e-6 
NUMMRG,KP,1e-6  

 AMESH,ALL		