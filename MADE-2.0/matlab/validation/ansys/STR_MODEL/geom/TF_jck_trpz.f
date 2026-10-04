!=======================================================================
! TF_jck_trpz.f - Jackets of every turn (constant thickness JT_w/JT_h of the layer,
! internal radius R, external radius R+JT). Stepped WP.
! Called by STR_MODEL.f (WP_TOPOLOGY other than RIS_only).
!=======================================================================
/PREP7
*GET,CFRG_TYPE,PARM,CFRG,TYPE		! 1 = per-grade jacket internal radius array defined in ParametriTF
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! JACKETS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_COMP = 'JCK'
THIS_MAT = JCK_MAT
THAT_MAT = JCK_MAT
THIS_ET = JCK_ET
ELM_SZ = 0.001 ! JT_h(1,1)/4
ARC_DIV = 6

THIS_REAL = 100
n=0
*DO,i,1,GRADES	
	*DO,j,1,NL(i)
		n=n+1
		*DO,k,1,NT(i)
			*GET,ULT_KP,KP,,NUM,MAX
			*IF,CFRG_TYPE,EQ,1,THEN
				R(n,k) = CFRG(i,1)			! Jacket internal radius of the grade (clamp(JT, r_SC_min, r_SC_max), as MADE)
			*ELSE
				R(n,k) = CFR			! parameter files without CFRG: one radius for the whole WP
			*ENDIF
			R_N      = R(n,k)                       
			EXT_RAD  = R_N+JT_h(n,k)       	    ! Jacket external radius		
			! ! ! JER(n,k) = 3.0*1.E-3**SI_UNIT       
			HJH(n,k) = JH(n,k)*0.5              ! Half jacket Height
			HJW(n,k) = LTOT*0.5-ICIT  			! Half jacket width
			HCW(n,k) = LTOT*0.5-JT_w(n,k)-ICIT  ! Half conductor width (jacket of THIS layer; was JT_w(1,1): with a variable jacket the turn overflowed its cell)			
			HCH(n,k) = JH(n,k)*0.5-JT_h(n,k)    ! Half conductor Height

			HCH_N = HCH(n,k)
			HCW_N = HCW(n,k)	
			JW_N  = 2*HJW(n,k)	
			JH_N  = JH(n,k)
			JT_Nh = JT_h(n,k)
			JT_Nw = JT_w(n,k)
		
			CSYS,CSYS_CABLE(n,k)+2
			THIS_REAL = THIS_REAL+1
			
			THIS_CS = CSYS_CABLE(n,k)
			
			! create jacket areas
		*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
		
			! TOP
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,,HCH_N
			K,ULT_KP+2,HCH_N*cos(pi/4),HCH_N*sin(pi/4)
			K,ULT_KP+3,HCW_N-R_N,HCH_N+JT_Nh
			K,ULT_KP+4,,HCH_N+JT_Nh
			K,ULT_KP+5,, 
			LSEL,NONE				
			LARC,ULT_KP+1,ULT_KP+2,ULT_KP+5,HCH_N
			*GET,ULT_LN,LINE,,NUM,MAX
			CM,TOP_%THIS_COMP%_%THIS_REAL%,LINE
			L,ULT_KP+3,ULT_KP+4
			L,ULT_KP+4,ULT_KP+1
			L,ULT_KP+2,ULT_KP+3		
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+1
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			!TP-RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCH_N*cos(pi/4),HCH_N*sin(pi/4)
			K,ULT_KP+2,HCW_N-R_N,HCH_N+JT_Nh
			K,ULT_KP+3,HCW_N+JT_Nw,HCH_N-R_N
			K,ULT_KP+4,HCW_N+JT_Nw-EXT_RAD,HCH_N+JT_Nh-EXT_RAD			
			LSEL,NONE
			L,ULT_KP+1,ULT_KP+2
			L,ULT_KP+1,ULT_KP+3
			LARC,ULT_KP+2,ULT_KP+3,ULT_KP+4,EXT_RAD
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+5
			! LESIZE,ALL,,,ARC_DIV
			! RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCH_N*cos(pi/4),HCH_N*sin(pi/4)
			K,ULT_KP+2,HCW_N+JT_Nw,HCH_N+JT_Nh-EXT_RAD
			K,ULT_KP+3,HCW_N+JT_Nw,-(HCH_N+JT_Nh-EXT_RAD)
			K,ULT_KP+4,HCH_N*cos(pi/4),-HCH_N*sin(pi/4)
			K,ULT_KP+5,,,
			LSEL,NONE
			LARC,ULT_KP+1,ULT_KP+4,ULT_KP+5,HCH_N
			CM,RGT_%THIS_COMP%_%THIS_REAL%,LINE
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+1,ULT_KP+2
			L,ULT_KP+2,ULT_KP+3
			L,ULT_KP+4,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+2
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			!BT-RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCH_N*cos(pi/4),-HCH_N*sin(pi/4)
			K,ULT_KP+2,HCW_N-R_N,-HCH_N-JT_Nh
			K,ULT_KP+3,HCW_N+JT_Nw,-HCH_N+R_N
			K,ULT_KP+4,HCW_N+JT_Nw-EXT_RAD,-HCH_N-JT_Nh+EXT_RAD			
			LSEL,NONE
			L,ULT_KP+1,ULT_KP+2
			L,ULT_KP+1,ULT_KP+3
			LARC,ULT_KP+2,ULT_KP+3,ULT_KP+4,EXT_RAD
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+5
			! BOT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,,-HCH_N
			K,ULT_KP+2,HCH_N*cos(pi/4),-HCH_N*sin(pi/4)
			K,ULT_KP+3,HCW_N-R_N,-HCH_N-JT_Nh
			K,ULT_KP+4,,-HCH_N-JT_Nh
			K,ULT_KP+5,, 
			LSEL,NONE										
			LARC,ULT_KP+1,ULT_KP+2,ULT_KP+5,HCH_N
			*GET,ULT_LN,LINE,,NUM,MAX
			CM,TOP_%THIS_COMP%_%THIS_REAL%,LINE
			L,ULT_KP+4,ULT_KP+1
			L,ULT_KP+2,ULT_KP+3
			L,ULT_KP+3,ULT_KP+4
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+3
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
		*ELSE
			! TOP
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,,HCH_N
			K,ULT_KP+2,HCW_N-R_N,HCH_N
			K,ULT_KP+3,HCW_N-R_N,HCH_N+JT_Nh
			K,ULT_KP+4,,HCH_N+JT_Nh
			LSEL,NONE
			L,ULT_KP+3,ULT_KP+4
			CM,TOP_%THIS_COMP%_%THIS_REAL%,LINE
			L,ULT_KP+4,ULT_KP+1
			L,ULT_KP+1,ULT_KP+2
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+2,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+1
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			!TP-RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCW_N-R_N,HCH_N
			K,ULT_KP+2,HCW_N,HCH_N-R_N
			K,ULT_KP+3,HCW_N+JT_Nw,HCH_N-R_N
			!K,ULT_KP+4,HCW_N+JT_Nw,HCH_N+JT_Nh-EXT_RAD
			K,ULT_KP+5,HCW_N-R_N,HCH_N+JT_Nh
			!K,ULT_KP+6,HCW_N-R_N,HCH_N+JT_Nh
			K,ULT_KP+7,HCW_N-R_N,HCH_N-R_N
			!K,ULT_KP+8,HCW_N+JT_Nw-EXT_RAD,HCH_N+JT_Nh-EXT_RAD
			LSEL,NONE
			! LARC,ULT_KP+4,ULT_KP+5,ULT_KP+8,EXT_RAD
			LARC,ULT_KP+3,ULT_KP+5,ULT_KP+7,EXT_RAD
			!L,ULT_KP+3,ULT_KP+4
			!L,ULT_KP+5,ULT_KP+6
			LESIZE,ALL,,,ARC_DIV
			CM,TRG_%THIS_COMP%_%THIS_REAL%,LINE
			! L,ULT_KP+6,ULT_KP+1
			
			LARC,ULT_KP+1,ULT_KP+2,ULT_KP+7,R_N
			LESIZE,ALL,,,ARC_DIV
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+5,ULT_KP+1
			L,ULT_KP+2,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+5
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			! LESIZE,ALL,,,ARC_DIV
			! RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCW_N,HCH_N-R_N
			K,ULT_KP+2,HCW_N,-HCH_N+R_N
			K,ULT_KP+3,HCW_N+JT_Nw,-HCH_N+R_N
			K,ULT_KP+4,HCW_N+JT_Nw,HCH_N-R_N
			LSEL,NONE
			L,ULT_KP+3,ULT_KP+4
			CM,RGT_%THIS_COMP%_%THIS_REAL%,LINE
			L,ULT_KP+4,ULT_KP+1
			L,ULT_KP+1,ULT_KP+2
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+2,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+2
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			!BT_RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,HCW_N-R_N,-HCH_N
			K,ULT_KP+2,HCW_N,-HCH_N+R_N
			K,ULT_KP+3,HCW_N+JT_Nw,-HCH_N+R_N	
			!K,ULT_KP+4,HCW_N+JT_Nw,-HCH_N-JT_Nh+EXT_RAD		
			K,ULT_KP+5,HCW_N-R_N,-HCH_N-JT_Nh		
			!K,ULT_KP+6,HCW_N-R_N,-HCH_N-JT_Nh	
			K,ULT_KP+7,HCW_N-R_N,-HCH_N+R_N		
			!K,ULT_KP+8,HCW_N+JT_Nw-EXT_RAD,-HCH_N-JT_Nh+EXT_RAD
			LSEL,NONE
			!LARC,ULT_KP+4,ULT_KP+5,ULT_KP+8,EXT_RAD
			LARC,ULT_KP+3,ULT_KP+5,ULT_KP+7,EXT_RAD
			!L,ULT_KP+3,ULT_KP+4
			!L,ULT_KP+5,ULT_KP+6
			LESIZE,ALL,,,ARC_DIV
			CM,BRG_%THIS_COMP%_%THIS_REAL%,LINE
			!L,ULT_KP+6,ULT_KP+1
			
			LARC,ULT_KP+1,ULT_KP+2,ULT_KP+7,R_N
			LESIZE,ALL,,,ARC_DIV
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+5,ULT_KP+1
			L,ULT_KP+2,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+5
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
			! LESIZE,ALL,,,ARC_DIV
			! BOT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			KSEL,NONE
			K,ULT_KP+1,,-HCH_N
			K,ULT_KP+2,HCW_N-R_N,-HCH_N
			K,ULT_KP+3,HCW_N-R_N,-HCH_N-JT_Nh
			K,ULT_KP+4,,-HCH_N-JT_Nh
			LSEL,NONE
			L,ULT_KP+3,ULT_KP+4
			CM,BOT_%THIS_COMP%_%THIS_REAL%,LINE
			L,ULT_KP+4,ULT_KP+1
			L,ULT_KP+1,ULT_KP+2
			*GET,ULT_LN,LINE,,NUM,MAX	
			L,ULT_KP+2,ULT_KP+3
			ASEL,NONE
			AL,ALL
			AATT,THIS_MAT,THIS_REAL,THIS_ET,THIS_CS+3
			LSEL,,LINE,,ULT_LN
			LATT,THIS_MAT,THIS_REAL
		*ENDIF	

		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		! CMSEL,,TRG_%THIS_COMP%_%THIS_REAL%,LINE
		! LCCAT,ALL
		! CMSEL,,BRG_%THIS_COMP%_%THIS_REAL%,LINE
		! LCCAT,ALL
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THAT_MAT,THIS_MAT
		ALLSEL,BELOW,AREA
		KSLL,,1
		NSLL,,1
		NUMMRG,NODE
		NUMMRG,KP

		ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
		ALLSEL,BELOW,AREA
		NUMMRG,NODE
		NUMMRG,KP
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		MSHKEY,1
		MSHAPE,0
		*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
			MSHKEY,0
			ESIZE,ELM_SZ*2 
		*ELSE
			ESIZE,ELM_SZ
		*ENDIF
		AMESH,ALL	
		CM,AREA_%THIS_COMP%_RGT,AREA

		LSEL,,LCCAT
		LDELE,ALL

		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		ALLSEL,BELOW,AREA
		NUMMRG,NODE,1e-4
		NUMMRG,KP,1e-4
		ARSYM,X,ALL
		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THAT_MAT,THIS_MAT
		ALLSEL,BELOW,AREA
		NUMMRG,NODE
		NUMMRG,KP

		ASEL,,REAL,,THIS_REAL
		ASEL,R,MAT,,THIS_MAT
		ASEL,R,LOC,X,-LTOT/2,LTOT/2
		CMSEL,U,AREA_%THIS_COMP%_RGT
		ESLA,,1
		ESEL,R,ESYS,,THIS_CS+5
		EMODIF,ALL,ESYS,THIS_CS+8
		ESLA,,1
		ESEL,R,ESYS,,THIS_CS+6
		EMODIF,ALL,ESYS,THIS_CS+7
		ESLA,,1
		ESEL,R,ESYS,,THIS_CS+2
		EMODIF,ALL,ESYS,THIS_CS+4
		*ENDDO			
	*ENDDO
*ENDDO

ALLSEL,ALL,ALL
NUMCMP,NODE
NUMCMP,KP
NUMCMP,ELEM
NUMCMP,LINE