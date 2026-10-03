!=======================================================================
! Materiali.f - Element types (MESH200 placeholders, switched to PLANE183/PLANE233 by
! BC\Contacts.f and Solution_EM.f) and material properties at 4.2 K.
! Called by STR_MODEL.f and EM_MODEL.f.
!=======================================================================
/PREP7

!!! MATERIAL PROPERTIES !!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! 
!!! Components of orthotropic materials have been defined
!!! with respect to local cylindrical coordinate systems
!!! X = THROUGH-THICKNESS DIRECTION (i.e. radial direction)
!!! Y = WRAPPING DIRECTION (i.e. azimuthal direction)
!!! Z = CABLE-AXIS DIRECTION (i.e. axial direction)
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! ELEMENT TYPE !!!

ET,1, 200,7 ! Conductor
ET,2, 200,6 ! SS Jacket
ET,3, 200,6 ! Turn-to-turn insulation
ET,4, 200,6 ! Filler material
ET,5, 200,6 ! Inter-layer insulation
ET,6, 200,6 ! Double layer insulation
ET,7, 200,6 ! Ground insulation
ET,8, 200,6 ! Buffer insulation
ET,9, 200,6 ! Wedge insulation
ET,10,200,6 ! Case
ET,11,200,6 ! Steel Filler
ET,20,200,7 ! Air

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

TOFFST,0              !DEFINE KELVIN AS TEMPERATURE UNITS
TREF,293.0

!!! MATERIALS DEFINITION !!!
!!! CONDUCTOR #1
MPTEMP,1,4.2 
MPDATA,EX,1,,10E+009          !UNITS ARE [N/m^2]
MPDATA,PRXY,1,,0.3  
MPDATA,DENS,1,,6100     
MPDATA,ALPX,1,,5.54E-6 ! 10.4E-6         
MP,MURX,1,1.0
MP,RSVX,1,1.0

!!! STEEL JACKET #2
MPTEMP,1,4.2
MPDATA,EX,2,,205000E+006       !UNITS ARE [N/m^2]
MPDATA,PRXY,2,,0.29		
MPDATA,DENS,2,,7900     
MPDATA,ALPX,2,,10.38E-6
MP,MURX,2,1

!!! TURN INSULATION #3
MPTEMP,1,4.2 
MPDATA,EX,3,,12000E+006       
MPDATA,EY,3,,20000E+006       
MPDATA,EZ,3,,20000E+006       
MPDATA,PRXY,3,,0.198
MPDATA,PRYZ,3,,0.17
MPDATA,PRXZ,3,,0.198
MPDATA,GXY,3,,6.E+9
MPDATA,GYZ,3,,6.E+9
MPDATA,GXZ,3,,6.E+9 
MPDATA,DENS,3,,1948    
MPDATA,ALPX,3,,24.22E-6         
MPDATA,ALPY,3,,8.65E-6        
MPDATA,ALPZ,3,,8.65E-6       
MP,MURX,3,1

!!! FILLER MATERIAL #4 
MPTEMP,1,4.2
MPDATA,EX,4,,7000E+006         
MPDATA,PRXY,4,,0.3
MPDATA,DENS,4,,1948   
MPDATA,ALPX,4,,17.3E-6        
MP,MURX,4,1

!!! INTER LAYER INSULATION #5
MPTEMP,1,4.2
MPDATA,EX,5,,12000E+006         
MPDATA,EY,5,,20000E+006         
MPDATA,EZ,5,,20000E+006         
MPDATA,PRXY,5,,0.198
MPDATA,PRYZ,5,,0.17
MPDATA,PRXZ,5,,0.198
MPDATA,GXY,5,,6.E+9
MPDATA,GYZ,5,,6.E+9
MPDATA,GXZ,5,,6.E+9
MPDATA,DENS,5,,1948   
MPDATA,ALPX,5,,24.22E-6         
MPDATA,ALPY,5,,8.65E-6        
MPDATA,ALPZ,5,,8.65E-6         
MP,MURX,5,1

!!! DOUBLE LAYER INSULATION #6
MPTEMP,1,4.2
MPDATA,EX,6,,12000E+006         
MPDATA,EY,6,,20000E+006         
MPDATA,EZ,6,,20000E+006         
MPDATA,PRXY,6,,0.198
MPDATA,PRYZ,6,,0.17
MPDATA,PRXZ,6,,0.198
MPDATA,GXY,6,,6.E+9
MPDATA,GYZ,6,,6.E+9
MPDATA,GXZ,6,,6.E+9
MPDATA,DENS,6,,1948   
MPDATA,ALPX,6,,24.22E-6         
MPDATA,ALPY,6,,8.65E-6        
MPDATA,ALPZ,6,,8.65E-6         
MP,MURX,6,1

!!! GROUND INSULATION #7
MPTEMP,1,4.2
MPDATA,EX,7,,12000E+006         
MPDATA,EY,7,,20000E+006         
MPDATA,EZ,7,,20000E+006         
MPDATA,PRXY,7,,0.198
MPDATA,PRYZ,7,,0.17
MPDATA,PRXZ,7,,0.198
MPDATA,GXY,7,,6.E+9
MPDATA,GYZ,7,,6.E+9
MPDATA,GXZ,7,,6.E+9
MPDATA,DENS,7,,1948   
MPDATA,ALPX,7,,24.22E-6         
MPDATA,ALPY,7,,8.65E-6        
MPDATA,ALPZ,7,,8.65E-6         
MP,MURX,7,1

!!! GROUND INSULATION #8
MPTEMP,1,4.2
MPDATA,EX,8,,12000E+006         
MPDATA,EY,8,,20000E+006         
MPDATA,EZ,8,,20000E+006         
MPDATA,PRXY,8,,0.198
MPDATA,PRYZ,8,,0.17
MPDATA,PRXZ,8,,0.198
MPDATA,GXY, 8,,6.E+9
MPDATA,GYZ, 8,,6.E+9
MPDATA,GXZ, 8,,6.E+9
MPDATA,DENS,8,,1948   
MPDATA,ALPX,8,,24.22E-6         
MPDATA,ALPY,8,,8.65E-6        
MPDATA,ALPZ,8,,8.65E-6         
MP,MURX,8,1

!!! WEDGE INSULATION #9
MPTEMP,1,4.2
MPDATA,EX,9,,12000E+006         
MPDATA,EY,9,,20000E+006         
MPDATA,EZ,9,,20000E+006         
MPDATA,PRXY,9,,0.198
MPDATA,PRYZ,9,,0.17
MPDATA,PRXZ,9,,0.198
MPDATA,GXY,9,,6.E+9
MPDATA,GYZ,9,,6.E+9
MPDATA,GXZ,9,,6.E+9
MPDATA,DENS,9,,1948   
MPDATA,ALPX,9,,24.22E-6         
MPDATA,ALPY,9,,8.65E-6        
MPDATA,ALPZ,9,,8.65E-6         
MP,MURX,9,1

!!! CASING #10
MPTEMP,1,4.2
MPDATA,EX,10,,205000E+006       
MPDATA,PRXY,10,,0.29		
MPDATA,DENS,10,,7900      
MPDATA,ALPX,10,,10.38E-6
MP,MURX,10,1

!!! WP STEEL FILLER #11
MPTEMP,1,4.2
MPDATA,EX,11,,205000E+006       
MPDATA,PRXY,11,,0.29		
MPDATA,DENS,11,,7900      
MPDATA,ALPX,11,,10.38E-6
MP,MURX,11,1

! MPTEMP,1,4.2
! MPDATA,EX,11,,7000E+006         
! MPDATA,PRXY,11,,0.3
! MPDATA,DENS,11,,1948   
! MPDATA,ALPX,11,,17.3E-6        
! MP,MURX,11,1

! MPDATA,EX,11,,205000E+006       
! MPDATA,PRXY,11,,0.29		
! MPDATA,DENS,11,,7900      
! MPDATA,ALPX,11,,10.38E-6
! MP,MURX,11,1


!!! AIR #20
MP,MURX,20,1
MP,RSVX,20,0