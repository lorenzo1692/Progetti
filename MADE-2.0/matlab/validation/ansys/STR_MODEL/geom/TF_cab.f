!=======================================================================
! TF_cab.f - Conductors of a WP with layers of equal width (older version of
! TF_cab_trpz.f). Not called by the current STR_MODEL.f / EM_MODEL.f.
!=======================================================================
/PREP7
LSEL,,LCCAT
LDELE,ALL
ALLSEL,ALL,ALL
NUMCMP,LINE
/PSYMB,CSYS,0

!!!!!!!!!!!!!! 
!!!!!!!!!!!!!! 
!!! CABLES !!! 
!!!!!!!!!!!!!! 
!!!!!!!!!!!!!! 

THIS_MAT = CAB_MAT
THIS_ET = CAB_ET

LHI_0 = INNLEG_OUTR-CASE_THICK-CASE_GIT-BUFF

*GET,MAX_CS,CDSY,0,NUM,MAX
! TOP_CS = MAX_CS+1
! RGT_CS = MAX_CS+2
! BOT_CS = MAX_CS+3
! LFT_CS = MAX_CS+4

! ESIZE,CFR

*DO,i,1,GRADES
	*IF,i,LE,NL_HIGHF,THEN
		*DO,k,1,NT_HIGHF
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX
			R(i,k) = CFR!JT_h(i,k)		    !Jacket internal radius
			R_N =  R(i,k)                       
			EXT_RAD = R_N+JT_h(i,k)       	    !Jacket external radius		
			JER(i,k) = 3.0*1.E-3**SI_UNIT       
			HJH(i,k) = JH(i,k)*0.5              !Half jacket Height
			HJW(i,k) = LTOT*0.5-ICIT  			!Half jacket width
			HCW(i,k) = LTOT*0.5-JT_w(i,k)-ICIT  !Half conductor width			
			HCH(i,k) = JH(i,k)*0.5-JT_h(i,k)    !Half conductor Height		 
											    
			*IF,i,EQ,1,THEN
					LHI(i,k) = LHI_0-(HJH(i,k)+ICIT) 
				*ELSE
					LHI(i,k) = LHI_0-(HJH(i,k)+ICIT+LWIT)-(HJH(i-1,k)+ICIT+LWIT) 
			*ENDIF

			LWI(i,k) = -WPW_HIGHF/2+LTOT/2*(2*k-1) 

			HCH_N = HCH(i,k)
			HCW_N = HCW(i,k)	
			JW_N = 2*HJW(i,k)	
			JH_N = JH(i,k)
			
			LHI_N = LHI(i,k)
			LWI_N = LWI(i,k)

			! Local CS to identify cable centre and later define orthotropic material properties
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+1,0,LWI_N,LHI_N,0.0,90
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+2,0,LWI_N,LHI_N,0.0,
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+3,0,LWI_N,LHI_N,0.0,-90
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+4,0,LWI_N,LHI_N,0.0,-180
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+5,1,LWI_N+HCW_N-R_N,LHI_N+HCH_N-R_N,0.0,45
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+6,1,LWI_N+HCW_N-R_N,LHI_N-HCH_N+R_N,0.0,-45
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+7,1,LWI_N-HCW_N+R_N,LHI_N-HCH_N+R_N,0.0,-135
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+8,1,LWI_N-HCW_N+R_N,LHI_N+HCH_N-R_N,0.0,135
			CSYS,REF_CS
			CSYS,MAX_CS+8*((k-1)+NT_HIGHF*(i-1))+2

			THIS_REAL = 100+NT_HIGHF*(i-1)+k

			! Create cable areas
			*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
				*GET,ULT_KP,KP,,NUM,MAX
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
				! ESIZE,R_N/4
				MSHKEY,0!1
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
LHI_0 =	LHI(i,k)		
	*ELSE		
		*DO,k,1,NT_LOWF
			
			KSEL,ALL
			*GET,ULT_KP,KP,,NUM,MAX
			R(i,k) = CFR! JT_h(i,k)/2			 !Jacket internal radius
			R_N = R(i,k)
			EXT_RAD = R_N+JT_h(i,k)       	 !Jacket external radius		
			
			HJH(i,k) = JH(i,k)*0.5               !Half Jacket Height
			HJW(i,k) = LTOT*0.5-ICIT  			 !Half jacket width
			HCW(i,k) = LTOT*0.5-JT_w(i,k)-ICIT   !Half Conductor width			
			HCH(i,k) = JH(i,k)*0.5-JT_h(i,k)     !Half Conductor Height		

			*IF,i,EQ,1,THEN
				LHI(i,k) = LHI_0-(HJH(i,k)+ICIT) 
				*ELSE
				LHI(i,k) = LHI_0-(HJH(i,k)+ICIT+LWIT)-(HJH(i-1,k)+ICIT+LWIT) 
			*ENDIF

			LWI(i,k) = -WPW_LOWF/2+LTOT/2*(2*k-1)

			HCH_N = HCH(i,k)
			HCW_N = HCW(i,k)	
			JW_N = 2*HJW(i,k)	
			JH_N = JH(i,k)

			LHI_N = LHI(i,k)
			LWI_N = LWI(i,k)

			! Local CS to identify cable centre and later define orthotropic material properties
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+1,0,LWI_N,LHI_N,0.0,90
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+2,0,LWI_N,LHI_N,0.0,
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+3,0,LWI_N,LHI_N,0.0,-90
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+4,0,LWI_N,LHI_N,0.0,-180
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+5,1,LWI_N+HCW_N-R_N,LHI_N+HCH_N-R_N,0.0,45
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+6,1,LWI_N+HCW_N-R_N,LHI_N-HCH_N+R_N,0.0,-45
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+7,1,LWI_N-HCW_N+R_N,LHI_N-HCH_N+R_N,0.0,-135
			CSYS,REF_CS
			CLOCAL,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+8,1,LWI_N-HCW_N+R_N,LHI_N+HCH_N-R_N,0.0,135
			CSYS,REF_CS
			CSYS,MAX_CS+8*(NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+(k-1))+2

			THIS_REAL = 100+NL_HIGHF*NT_HIGHF+NT_LOWF*(i-NL_HIGHF-1)+k

			! Create cable areas
			*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
				*GET,ULT_KP,KP,,NUM,MAX
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
				! ESIZE,R_N/4
				MSHKEY,0!1
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
	*ENDIF
LHI_0 = LHI(i,k)  
*ENDDO
