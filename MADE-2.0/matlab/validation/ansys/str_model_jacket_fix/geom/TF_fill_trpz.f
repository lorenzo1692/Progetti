/PREP7
*GET,CFRG_TYPE,PARM,CFRG,TYPE		! 1 = per-grade jacket internal radius array defined in ParametriTF
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
			JER(n,k) = 3.0*1.E-3**SI_UNIT       
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
			
			IN_RAD = R_N+JT_h(n,k)       	!Jacket external radius		
			EXT_RAD = IN_RAD+ICIT
		
			CSYS,CSYS_CABLE(n,k)+2
			THIS_REAL = THIS_REAL+1
			
			THIS_CS = CSYS_CABLE(n,k)
			
			! Create filler material areas
			!TP-RGT
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX	
			K,ULT_KP+1,HCW_N-R_N,HCH_N+JT_Nh+ICIT
			K,ULT_KP+2,HCW_N+JT_Nw+ICIT,HCH_N+JT_Nh+ICIT
			K,ULT_KP+3,HCW_N+JT_Nw+ICIT,HCH_N+JT_Nh-IN_RAD
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
			K,ULT_KP+1,HCW_N-R_N,-HCH_N-JT_Nh-ICIT
			K,ULT_KP+2,HCW_N+JT_Nw+ICIT,-HCH_N-JT_Nh-ICIT
			K,ULT_KP+3,HCW_N+JT_Nw+ICIT,-HCH_N-JT_Nh+IN_RAD
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
	*ENDDO
*ENDDO

ASEL,,MAT,,THAT_MAT,THIS_MAT,,1
ALLSEL,BELOW,AREA
NUMMRG,NODE,1e-6 
NUMMRG,KP,1e-6  

AMESH,ALL	