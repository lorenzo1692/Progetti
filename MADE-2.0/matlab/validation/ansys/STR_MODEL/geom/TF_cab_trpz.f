!=======================================================================
! TF_cab_trpz.f - Conductors (cable areas) of every turn, rectangular cable with rounded
! corners (fillet R = jacket internal radius). Layers of different width
! (stepped WP). Also defines the local CS of every cable (CSYS_CABLE),
! used by the following macros and for the orthotropic materials.
! Called by STR_MODEL.f and EM_MODEL.f (WP_TOPOLOGY other than RIS_only).
!=======================================================================
/PREP7
*GET,CFRG_TYPE,PARM,CFRG,TYPE		! 1 = per-grade jacket internal radius array defined in ParametriTF
/PSYMB,CSYS,0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
THIS_MAT = CAB_MAT
THIS_ET = CAB_ET
THIS_ET_3D = CAB_ET_3D

LHI_0 = Ri_-CASE_THICK-GIT-(JH(1,1)*0.5+ICIT) 
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

			*IF,n,NE,1,THEN,
				LHI(n,k) = LHI_0-(HJH(n,1)+HJH(n-1,1)+2*ICIT+LWIT) 
			*ELSE
				LHI(n,k) = LHI_0 
			*ENDIF

			LWI(n,k) = -W_GRADES(i)/2+LTOT/2*(2*k-1) 

			HCH_N = HCH(n,k)
			HCW_N = HCW(n,k)	
			JW_N  = 2*HJW(n,k)	
			JH_N  = JH(n,k)
			
			LHI_N = LHI(n,k)
			LWI_N = LWI(n,k)

			! Local CS to identify cable centre and later define orthotropic material properties
			*GET,CSYS_CABLE(n,k),CDSY,0,NUM,MAX
			CSYS,REF_CS
			CLOCAL,CSYS_CABLE(n,k)+1,0,LWI_N,LHI_N,0.0,90
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+2,0,LWI_N,LHI_N,0.0,
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+3,0,LWI_N,LHI_N,0.0,-90
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+4,0,LWI_N,LHI_N,0.0,-180
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+5,1,LWI_N+HCW_N-R_N,LHI_N+HCH_N-R_N,0.0,45
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+6,1,LWI_N+HCW_N-R_N,LHI_N-HCH_N+R_N,0.0,-45
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+7,1,LWI_N-HCW_N+R_N,LHI_N-HCH_N+R_N,0.0,-135
			CSYS,REF_CS   
			CLOCAL,CSYS_CABLE(n,k)+8,1,LWI_N-HCW_N+R_N,LHI_N+HCH_N-R_N,0.0,135
						
			CSYS,CSYS_CABLE(n,k)+2
			THIS_REAL = THIS_REAL+1
			
			! Create cable areas
			*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
				K,ULT_KP+1,,HCH_N
				K,ULT_KP+2,HCH_N,0
				K,ULT_KP+3,0
				K,ULT_KP+4,HCH_N*cos(pi/4),HCH_N*sin(pi/4)
				LSEL,NONE		
				LARC,ULT_KP+1,ULT_KP+4,ULT_KP+3,HCH_N
				LARC,ULT_KP+2,ULT_KP+4,ULT_KP+3,HCH_N
				CM,EXT_CAB_%THIS_REAL%,LINE
				L,ULT_KP+3,ULT_KP+2
				L,ULT_KP+1,ULT_KP+3				
				KSLL,,1
				NSLL,,1
				NUMMRG,NODE
				NUMMRG,KP
				ASEL,NONE
				AL,ALL
				AATT,THIS_MAT,THIS_REAL,THIS_ET
				ESIZE,HCH_N/4
				MSHKEY,0
				MSHAPE,0
				AMESH,ALL
				ASEL,,REAL,,THIS_REAL
				ASEL,R,MAT,,1
				ASEL,R,LOC,X,-LTOT/2,LTOT/2
				ALLSEL,BELOW,AREA
				NUMMRG,NODE
				NUMMRG,KP
				ARSYM,X,ALL
				ARSYM,Y,ALL
				ALLSEL,BELOW,AREA
				NUMMRG,NODE
				NUMMRG,KP
				LSEL,,EXT
				LATT,THIS_MAT,THIS_REAL
				LSLA,S
				LSEL,R,EXT
				CM,EXT_CAB_%THIS_REAL%,LINE
			*ELSE
				K,1+ULT_KP,,HCH_N-R_N
				K,2+ULT_KP,,-HCH_N+R_N
				K,3+ULT_KP,HCW_N-R_N,-HCH_N+R_N
				K,4+ULT_KP,HCW_N-R_N,+HCH_N-R_N
				LSEL,NONE
				L,ULT_KP+2,ULT_KP+3
				L,ULT_KP+3,ULT_KP+4
				L,ULT_KP+4,ULT_KP+1
				CM,INT_CAB_%THIS_REAL%,LINE
				L,ULT_KP+1,ULT_KP+2
				ASEL,NONE
				AL,ALL
				AATT,THIS_MAT,THIS_REAL,THIS_ET
				ESIZE,(HCH_N-R_N)/4
				MSHKEY,1
				MSHAPE,0
				AMESH,ALL
				KSEL,ALL
				*GET,ULT_KP,KP,,NUM,MAX
				K,ULT_KP+1,,-HCH_N
				K,ULT_KP+2,HCW_N-R_N,-HCH_N
				K,ULT_KP+3,HCW_N,-HCH_N+R_N
				K,ULT_KP+4,HCW_N,HCH_N-R_N
				K,ULT_KP+5,HCW_N-R_N,HCH_N
				K,ULT_KP+6,,HCH_N
				K,ULT_KP+7,HCW_N-R_N,-HCH_N+R_N
				K,ULT_KP+8,HCW_N-R_N,HCH_N-R_N
				K,ULT_KP+9,,-HCH_N+R_N
				K,ULT_KP+10,,HCH_N-R_N
				LSEL,NONE
				L,ULT_KP+1,ULT_KP+2
				LARC,ULT_KP+2,ULT_KP+3,ULT_KP+7,R_N
				L,ULT_KP+3,ULT_KP+4
				LARC,ULT_KP+4,ULT_KP+5,ULT_KP+8,R_N
				L,ULT_KP+5,ULT_KP+6
				CM,EXT_CAB_%THIS_REAL%,LINE
				L,ULT_KP+6,ULT_KP+10
				L,ULT_KP+9,ULT_KP+1
				CMSEL,A,INT_CAB_%THIS_REAL%,LINE
				KSLL,,1
				NSLL,,1
				NUMMRG,NODE
				NUMMRG,KP
				ASEL,NONE
				AL,ALL
				AATT,THIS_MAT,THIS_REAL,THIS_ET
				CMSEL,,INT_CAB_%THIS_REAL%,LINE
				LCCAT,ALL
				CMSEL,,EXT_CAB_%THIS_REAL%,LINE
				LCCAT,ALL
				LSLA,A,1
				ESIZE,R_N/4
				MSHKEY,1
				MSHAPE,0
				AMESH,ALL
				LSEL,,LCCAT
				LDELE,ALL
				ASEL,,REAL,,THIS_REAL
				ASEL,R,MAT,,1
				ASEL,R,LOC,X,-LTOT/2,LTOT/2
				ALLSEL,BELOW,AREA
				NUMMRG,NODE
				NUMMRG,KP
				ARSYM,X,ALL
				ALLSEL,BELOW,AREA
				NUMMRG,NODE
				NUMMRG,KP
				LSEL,,EXT
				LATT,THIS_MAT,THIS_REAL	
			*ENDIF
		*ENDDO
		LHI_0 =	LHI(n,k)	
	*ENDDO
*ENDDO
 	
ALLSEL,ALL,ALL
NUMCMP,NODE
NUMCMP,KP
NUMCMP,ELEM