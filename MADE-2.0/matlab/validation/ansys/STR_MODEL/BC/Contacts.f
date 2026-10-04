!=======================================================================
! Contacts.f - Structural element types (PLANE183, generalized plane strain) and the
! contact pairs (CONTA172/TARGE169): conductor-jacket of every turn,
! ground insulation-case and case-wedge insulation.
! Called by STR_MODEL.f.
!=======================================================================
!---------------------------------------------------!
/PREP7
ET,1, 183,,,5
ET,2, 183,,,5
ET,3, 183,,,5 
ET,4, 183,,,5
ET,5, 183,,,5
ET,6, 183,,,5
ET,7, 183,,,5
ET,8, 183,,,5
ET,9, 183,,,5
ET,10,183,,,5
!---------------------------------------------------!
!!! Modeling structural contacts
!---------------------------------------------------!
ET,11,169
ET,12,172
KEYOPT,12,1,0
KEYOPT,12,2,0
KEYOPT,12,4,0
KEYOPT,12,5,0
KEYOPT,12,6,2
KEYOPT,12,9,3
KEYOPT,12,10,2
KEYOPT,12,12,0

ET,14,172
KEYOPT,14,1,0
KEYOPT,14,2,0
KEYOPT,14,4,0
KEYOPT,14,5,0
KEYOPT,14,6,2
KEYOPT,14,9,3
KEYOPT,14,10,2
KEYOPT,14,12,0

STIFF = 10.0
PENTR = 0.0
PNBLL = 0.0

! Conductor-jacket contact
MP,MU,12,0.20
MAT,12

ALLSEL,ALL,ALL
*GET,MAX_REAL,ELEM,,RELM
CSYS,REF_CS
CAB_CONT = 1.0
THIS_REAL = 100

n=0
*DO,i,1,GRADES	
	*DO,j,1,NL(i)
		n=n+1
		*DO,k,1,NT(i)
		THIS_REAL = THIS_REAL+1
		R,THIS_REAL
		REAL,THIS_REAL
		R,THIS_REAL,,,STIFF,PENTR,0,PNBLL 
		! Generate the target surface  
		LSEL,,MAT,,JCK_MAT
		LSEL,R,REAL,,THIS_REAL
		TYPE,11  
		NSLL,S,1
		ESLN,S
		ESURF   
		! Generate the contact surface  
		LSEL,,MAT,,CAB_MAT
		LSEL,R,REAL,,THIS_REAL
		*IF,CABLE_TYPE(i,1),EQ,'RIS',THEN
			TYPE,14 
		*ELSE
			TYPE,12 
		*ENDIF
		NSLL,S,1
		ESLN,S,0
		ESURF   
		ALLSEL  
		*ENDDO	
	*ENDDO
*ENDDO
!---------------------------------------------------!
! ground_insulation-case
ALLSEL,ALL,ALL

CSYS,REF_CS
ET,15,172
KEYOPT,15,1,0
KEYOPT,15,2,0
KEYOPT,15,4,0
KEYOPT,15,5,0
KEYOPT,15,6,2
KEYOPT,15,9,5
KEYOPT,15,10,2
KEYOPT,15,12,0

MP,MU,13,0.2
MAT,13

STIFF = 1.0

LSEL,S,MAT,,CASE_MAT
NSLL,S,1
*GET,Y1,NODE,,MNLOC,Y
*GET,Y2,NODE,,MXLOC,Y
*GET,X1,NODE,,MNLOC,X
*GET,X2,NODE,,MXLOC,X

LSEL,S,MAT,,CASE_MAT
LSEL,R,LOC,Y,Y1
NSLL,S,1
CM,CASE_MAT_LOW,NODE
!*
LSEL,S,MAT,,BUFF_MAT
LSEL,R,EXT
LSEL,R,LOC,Y,Y1
NSLL,S,1
CM,BUFF_MAT_LOW,NODE

LSEL,S,MAT,,CASE_MAT
LSEL,R,LOC,Y,Y2
NSLL,S,1
CM,CASE_MAT_UP,NODE
!*
LSEL,S,MAT,,BUFF_MAT
LSEL,R,EXT
LSEL,R,LOC,Y,Y2
NSLL,S,1
CM,BUFF_MAT_UP,NODE

LSEL,S,MAT,,CASE_MAT
LSEL,U,LOC,Y,Y1
LSEL,U,LOC,Y,Y2
LSEL,U,LOC,X,0,1E9
NSLL,S,1
CM,CASE_MAT_LEFT,NODE
!*
LSEL,S,MAT,,BUFF_MAT
LSEL,R,EXT
LSEL,U,LOC,Y,Y1
LSEL,U,LOC,Y,Y2
LSEL,U,LOC,X,0,1E9
NSLL,S,1
CM,BUFF_MAT_LEFT,NODE

LSEL,S,MAT,,CASE_MAT
LSEL,U,LOC,Y,Y1
LSEL,U,LOC,Y,Y2
LSEL,U,LOC,X,-1E9,0
NSLL,S,1
CM,CASE_MAT_RIGHT,NODE
!*
LSEL,S,MAT,,BUFF_MAT
LSEL,R,EXT
LSEL,U,LOC,Y,Y1
LSEL,U,LOC,Y,Y2
LSEL,U,LOC,X,-1E9,0
NSLL,S,1
CM,BUFF_MAT_RIGHT,NODE

! *GET,MAX_REAL,ELEM,,RELM
! R,MAX_REAL+1
! REAL,MAX_REAL+1
! R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL 
! ! Generate the target surface  
! NSEL,S,NODE,,CASE_MAT_UP
! TYPE,11  
! ESLN,S 
! ESURF   
! ! Generate the contact surface  
! NSEL,S,NODE,,BUFF_MAT_UP
! TYPE,15  
! ESLN,S 
! ESURF   
! ALLSEL 

*GET,MAX_REAL,ELEM,,RELM
R,MAX_REAL+1
REAL,MAX_REAL+1
R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL 
! Generate the target surface  
NSEL,S,NODE,,CASE_MAT_LOW
TYPE,11  
ESLN,S 
ESURF   
! Generate the contact surface  
NSEL,S,NODE,,BUFF_MAT_LOW
TYPE,15  
ESLN,S 
ESURF   
ALLSEL 

*GET,MAX_REAL,ELEM,,RELM
R,MAX_REAL+1
REAL,MAX_REAL+1
R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL 
! Generate the target surface  
NSEL,S,NODE,,CASE_MAT_LEFT
TYPE,11  
ESLN,S 
ESURF   
! Generate the contact surface  
NSEL,S,NODE,,BUFF_MAT_LEFT
TYPE,15  
ESLN,S 
ESURF   
ALLSEL 

*GET,MAX_REAL,ELEM,,RELM
R,MAX_REAL+1
REAL,MAX_REAL+1
R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL 
! Generate the target surface  
NSEL,S,NODE,,CASE_MAT_RIGHT
TYPE,11  
ESLN,S 
ESURF   
! Generate the contact surface  
NSEL,S,NODE,,BUFF_MAT_RIGHT
TYPE,15  
ESLN,S 
ESURF   
ALLSEL 

! *GET,MAX_REAL,ELEM,,RELM
! R,MAX_REAL+1
! REAL,MAX_REAL+1
! R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL 
! ! Generate the target surface  
! LSEL,S,MAT,,CASE_MAT
! TYPE,11  
! NSLL,S,1
! ESLN,S,0
! ESURF   
! ! Generate the contact surface  
! LSEL,S,MAT,,BUFF_MAT
! LSEL,R,EXT
! TYPE,15  
! NSLL,S,1
! ESLN,S,0
! ESURF   
! ALLSEL  
!---------------------------------------------------!
! case wedge contacts
ALLSEL,ALL,ALL
*GET,MAX_REAL,ELEM,,RELM
CSYS,REF_CS

ET,13,172
KEYOPT,13,1,0
KEYOPT,13,2,0
KEYOPT,13,4,0
KEYOPT,13,5,0
KEYOPT,13,6,2
KEYOPT,13,9,5
KEYOPT,13,10,2
KEYOPT,13,12,0

MP,MU,14,0.20
MAT,14

STIFF = 10.0
PENTR = 0.1
PNBLL = 1

CSYS,REF_CS
R,MAX_REAL+1
REAL,MAX_REAL+1
R,MAX_REAL+1,,,STIFF,PENTR,0,PNBLL
! Generate the target surface  
ASEL,,MAT,,CASE_MAT
ALLSEL,BELOW,AREA
CSYS,WEDGE_CS_SYMMs
KSEL,R,LOC,X,CASE_INS
LSLK,,1
TYPE,11  
NSLL,S,1
ESLN,S,0
ESURF   
! Generate the contact surface  
ASEL,,MAT,,WINS_MAT
ALLSEL,BELOW,AREA
CSYS,WEDGE_CS_SYMMs
KSEL,R,LOC,X,CASE_INS
LSLK,,1
TYPE,13 
NSLL,S,1
ESLN,S,0
ESURF   
ALLSEL 
 
!*!*!*!*!

CSYS,REF_CS_SYMM
R,MAX_REAL+2
REAL,MAX_REAL+2
R,MAX_REAL+2,,,STIFF,PENTR,0,PNBLL
! Generate the target surface  
ASEL,,MAT,,CASE_MAT
ALLSEL,BELOW,AREA
CSYS,WEDGE_CS_SYMMd
KSEL,R,LOC,X,-CASE_INS
LSLK,,1
TYPE,11  
NSLL,S,1
ESLN,S,0
ESURF   
! Generate the contact surface  
ASEL,,MAT,,WINS_MAT
ALLSEL,BELOW,AREA
CSYS,WEDGE_CS_SYMMd
KSEL,R,LOC,X,-CASE_INS
LSLK,,1
TYPE,13  
NSLL,S,1
ESLN,S,0
ESURF   
ALLSEL  
CSYS,REF_CS
