!**************************************************************
!         *************************************************
!         *                                               *
!         *     twinsol.f90                               *
!         *                                               *
!         *     UM-UND Solenoid Ray Tracing Program       *
!         *                                               *
!         *     adapted by J. Meissner (from original     *
!         *     program by Wenzhong Liu) for UNIX systems *
!         *     and ACE/gr X-Window output.               *
!         *     modified by M. Lee using fortran90 for    *
!         *          double solenoid system (10/95)       *
!         *          radial electrode       (1/96)        *
!         *          more accurate BigSol specs           *
!         *          using Cryomagnetics Spec Info        *
!         *          use smaller step size: DZ =0.002     *
!         *                                 (2/21/96)     * 
!         *          include Ziegler's STOP routine       *
!         *             to calculate stopping for delta-E *
!         *               absorbers        (3/28/96)      *  
!         *          include defocusing solenoid          *
!         *             doesn't work!!     (4/21/96)      *
!         *          include gas filled mode              *
!         *                                (4/25/96)      *
!         *          include x-y plot      (4/30/96)      *
!         *          change twinsol solenoid parameters   *
!         *             using coild dimensions provided   *
!         *             by Cryomagnetics   (6/24/96)      *
!         *                               1/96            *
!         *************************************************

!**********************************************************************
!                                                      
!  the FORTRAN units used in this program are:        
!                                     
!    the usual ...        
!    5     standard input  (keyboard)             
!    6     standard output (terminal screen)          
!                                                          
!    and also ...                                       
!    9     FILENAME.DAT    (input parameter file)          
!    10    RAYS.DAT        (output file for rho vs. z view of trajectory)
!    11    RAYS2.DAT       (output file for axial view of trajectory)
!    12    RAYS3.DAT       (output file xyz)
!**********************************************************************

       PROGRAM twinsol_f90
       
       IMPLICIT NONE
       
!            ************
!            dummy vars
!            ************
       INTEGER XYZ,npth,iraz
       REAL ZYX, TEMP, ZZCM
       
!            ************
!            input file
!            ************
       CHARACTER *30 FILENAME
       REAL    ZO(700),AMU(700),EO(700),THETA(700)
       REAL    RN
       INTEGER NR

!            ************
!            double or single?
!            ************
       CHARACTER *1 MODE1
       REAL    PAR_SOLENOID
       INTEGER SOLNAME1,SOLNAME2

!            ************
!            solenoid(s)
!            ************
       REAL    S1,LL,LLS,ER,ERS
       REAL    MA1,MA2,I1,I2  
       REAL    ZTGT,Z1,Z2,SEPAR
       REAL    CC,CCC,BM1,BM2
       REAL    LENGTH,BORE

!            ************
!            radial electrode
!            ************
       CHARACTER *1 MODE2  
       INTEGER     E_RAD
       REAL    V_INNER,V_OUTER
       REAL    R_INNER,R_OUTER,L_ELECTRODE,Z_E_1,Z_E_2,phy

!            ************
!            mid-plane absorber
!            ************
       CHARACTER *1 MODE3
       INTEGER ION, ABSORBER, POS_ABSORBER
       REAL THICKNESS_ABSORBER,POSITION_ABSORBER,AMU_ABSORBER
       REAL M1,M2,RHO,ATRHO,VFERMI,LFCTR
       REAL, DIMENSION(8) :: PCOEF
       INTEGER MM1
       REAL EPSIL, SE, SN, VELSQR
       
!            ************
!            defocusing solenoid
!            ************
       CHARACTER *1 MODE4
       REAL  DS_S1, DS_LL, DS_ER, DS_ERS
       REAL  DS_MA, DS_Z1, DS_CC, DS_BM, BM3
       REAL  DS_LENGTH, DS_BORE
!            note that for the defocusing element S1=LL=LENGTH and ER=BORE
       
!            ************
!            gas absorber
!            ************
       CHARACTER *1 MODE5
!      use the same variables as in the midplane absorber
       REAL TORR
            
!            ************
!            particle vars
!            ************
       REAL    N,M,E,T,ANG3,V0,VX,VY,VZ,V1,V2,V3              
       REAL    RR,RRS,INRR,INRRS
       REAL    X0,Y0,Z0,X,Y,Z,GX,GY,GZ,HGX,HGY,HGZ
       REAL    F7,F8,F9,HF7,HF8,HF9
       REAL    XX,YY,ZZ,XXS,YYS,xxr,yyr
       REAL    AX,AY,AZ

!            ************
!            calculational
!            ************
       REAL    ZND,DZ,HDZ,DT,DTS,CQM,XP,YP

!            ************
!            field variables
!            ************
       REAL    BZSL,BRSL,BR,BZ,BX,BY
       REAL    E_FIELD,EX,EY,EZ
       
!            ************
!            loop indices
!            ************
 
       INTEGER I,STEP,STEPNUM,J, J_MAX
       REAL    STEPNOM

!**********************************************************************
! initialize variables for  use in program      
!**********************************************************************
       
       V0=0.0
       HDZ=0.0
       EX=0.0
       EY=0.0
       EZ=0.0
       E_FIELD = 0.0

!**********************************************************************
! read in the input parameters                               
!**********************************************************************
!      WRITE(6,*)' Enter the input parameter file (FILENAME.DAT) '
      READ (5,'(A)') FILENAME
      OPEN (UNIT=9,FILE=FILENAME)

!**********************************************************************
!     example of a *.dat file is shown below:                       *
!                                  *                              *
!         5                        *  5 is read as RN, the number *
!         2,4.002603,16.4863,1.    *  of rays. Then charge, mass, *
!         2,4.002603,16.4835,2.    *  energy and angle are read   *
!         2,4.002603,16.4789,3.    *  (ZO,AMU,EO,THETA) in units  *
!         2,4.002603,16.4725,4.    *  n,AMU,MeV,degrees. RN can   *
!         2,4.002603,16.4642,5.    *  be up to 20 rays. Make sure *
!                                  *  to include a <return> after *
!                                  *  the last data point.        *
!**********************************************************************

      READ(9,*)RN
      NR=INT(RN)

      DO I=1,NR
!      WRITE(6,*)'Now loading input parameters ...'
      READ(9,*)ZO(I),AMU(I),EO(I),THETA(I)
      END DO




!*************************************
!            rotation phy
!        included by  OCBSantos
!*************************************       
               READ (5,*) phy

!**********************************************************************
!      EDITED 8/9/90 FOR TRAJ. PLOT                         
!      NOTRE DAME L'il SOLENOID SETUP                          
!      OBJECT DISTANCE 52.3 CM                           
!      IMAGE DISTANCE 157.2 CM                         
!      SUM OF THE DISTANCES 209.5 CM             
!      I=125 amps BM=3.5 TELSA AT SOLENOID CENTER  MA1=2.125   
!      PARAMETERS FOR THE PARTICLE             
!      UNITS :                                            
!               LENGTH[METERS] TIME[SECONDS]                
!               CHARGE[COULOMBS] VOLTAGE[VOLTS]            
!               M[AMU] E[MeV] N[CHARGE OF PROTON]            
!**********************************************************************

!**********************************************************************
!  DZ is the step size for plotting
!**********************************************************************

      DZ=.0020000
      HDZ=.0010000

!**********************************************************************
!  initialize the values for the different solenoids
!     S1 is the effective length of the solenoid
!     ER is the effective radius of the solenoid
!     MA1 is a scaling factor that will be used later to calculate
!           the magnetic field -- see the variables BM and MAERS1 later 
!           in this program.  See also the appendix at end.
!     Z1 is the entrance position of the solenoid
!**********************************************************************

!      WRITE(6,*) 'Distances'
!      WRITE(6,*) 'ZTGT = dist from tgt to (first) Solenoid center '
!      WRITE(6,*) 'Enter target distance ZTGT (meters):           '
      READ (5,*) ZTGT
!      WRITE(6,*) 'Enter max. distance to plot rays ZND (meters): '
      READ (5,*) ZND
      STEPNOM = ZND/DZ
      STEPNUM = NINT(STEPNOM)
      iraz=stepnum/75+1
      stepnum=iraz*75
      znd=stepnum*dz

      !Lendo valores de X0(m), Y0(m), Z0(m), modificado por OCBSantos
      READ (5,*) X0
      READ (5,*) Y0
      READ (5,*) Z0

!      WRITE(6,*) 'Single or Double Solenoid? (enter "1" or "2")  '
      READ (5,*) SOLNAME1
      IF ((SOLNAME1 .NE. 1) .AND. (SOLNAME1 .NE. 2)) THEN
         WRITE(6,*) 'You have selected an invalid choice.'
         WRITE(6,*) 'I will assume you want only one solenoid.'
         SOLNAME1 = 1
      ENDIF
      CALL SOL_INPUT_VAL(SOLNAME1,I1,I2,S1,ER,MA1,MA2,LENGTH,BORE,Z1,Z2,ZTGT,SEPAR,MODE1,PAR_SOLENOID)

      
!*************************************
!  For Radial Electric field
!*************************************      
!      WRITE(6,*) 'Do you want to include cylindrical electrode? (y/n)'
      READ (5,*) MODE2
      IF ((ICHAR(MODE2) .EQ. 89) .OR. (ICHAR(MODE2) .EQ. 121)) THEN
         E_RAD = 1
      ELSE
         E_RAD = 0
      ENDIF
      IF (E_RAD .EQ. 1) THEN
         WRITE(6,*) 'Enter the following parameters in meters and volts as appropriate.'
         WRITE(6,*) 'Length of the electrode: '
         READ (5,*) L_ELECTRODE
         WRITE(6,*) 'Inner electrode radius: '
         READ (5,*) R_INNER
         WRITE(6,*) 'Outer electrode radius: '
         READ (5,*) R_OUTER
         WRITE(6,*) 'Inner electrode potential: '
         READ (5,*) V_INNER
         WRITE(6,*) 'Outer electrode potential: '
         READ (5,*) V_OUTER
         WRITE(6,*) 'Center to center distance from (first) solenoid to electrode: ' 
         READ (5,*) Z_E_1
         
         Z_E_1 = ZTGT + Z_E_1 - 0.5*L_ELECTRODE
         Z_E_2 = Z_E_1 + L_ELECTRODE
      ENDIF
     
!*************************************
!  For Mid-Plane Absorber
!************************************* 
!      WRITE(6,*)'Do you want to include a mid-plane absorber?'
      READ(5,*) MODE3
      MODE5='n'  !*** don't allow gas absorber (see below) ***!
      IF ((ICHAR(MODE3) .EQ. 89) .OR. (ICHAR(MODE3) .EQ. 121)) THEN
         WRITE(6,*)'Enter the atomic number of the absorber.'
         READ(5,*)ABSORBER
         WRITE(6,*)'Enter the distance from the target to the absorber (m).'
         READ(5,*)POSITION_ABSORBER
         POS_ABSORBER = POSITION_ABSORBER/DZ
         WRITE(6,*)'Enter the thickness the target (mm).'
         READ(5,*)THICKNESS_ABSORBER
      ENDIF
      
!*************************************
!  For defocusing solenoid
!************************************* 
!      WRITE(6,*)'Do you want to include a defocusing element?'
      READ(5,*) MODE4
      IF ((ICHAR(MODE4) .EQ. 89) .OR. (ICHAR(MODE4) .EQ. 121)) THEN
         WRITE(6,*)'Enter the radius of the defocusing element.'
         READ(5,*)DS_ER
         DS_ERS = DS_ER*DS_ER
         WRITE(6,*)'Enter the length of the defocusing element.'
         READ(5,*)DS_LENGTH         
         WRITE(6,*)'Enter the distance from the target to the defocusing element (m).'
         READ(5,*)DS_Z1
         DS_CC = DS_Z1 - 0.5*DS_LENGTH 
         WRITE(6,*)'Enter the maximum magnetic field of the defocusing element.'
         READ(5,*)DS_MA        
         DS_MA = DS_MA * SQRT(DS_ERS + DS_LENGTH*DS_LENGTH/4) / DS_LENGTH
         DS_S1 = DS_LENGTH
         DS_LL = DS_LENGTH
         DS_BORE = DS_ER
      ENDIF

!*************************************
!  Gas Absorber
!************************************* 
      IF ((ICHAR(MODE3) .NE. 89) .OR. (ICHAR(MODE3) .NE. 121)) THEN    
!         WRITE(6,*)'Do you want a gas absorber in the first magnet?'
         READ(5,*) MODE5
         IF ((ICHAR(MODE5) .EQ. 89) .OR. (ICHAR(MODE5) .EQ. 121)) THEN
            WRITE(6,*)'Enter the atomic number of the gas.'
            READ(5,*)ABSORBER
            POS_ABSORBER = (ZTGT - LENGTH/2) / DZ
            WRITE(6,*)'Enter the pressure of the gas (torr).'
            READ(5,*) TORR
         ENDIF
      ENDIF
      
      
!**********************************************************************
!  SOME FIXED FIELD RELATIONS
!**********************************************************************
      
!     ***********************
!     * effective length of solenoid
!     ***********************  
      LL = S1                

!     ***********************
!     * entrance position of solenoid
!     ***********************
      CC = Z1                
      CCC = Z2

!     ************************
!     * generally an "S" indicates a squared variable
!     ************************
      LLS=LL*LL
      ERS=ER*ER

!     ***********************
!     * BM = Max B at the center of the solenoid
!     ***********************
      BM1=MA1*LL/SQRT(ERS+LLS/4)
      BM2=MA2*LL/SQRT(ERS+LLS/4)
      BM3=DS_MA
      
!**********************************************************************
!                 START PLOTTING THE TRAJECTORIES
!**********************************************************************
      

      OPEN (UNIT=10,FILE='rays.dat')
!     WRITE(10,'(a,f5.1,a,f5.1,a,a,a)') '@    subtitle "I1:',I1,',I2:',I2,'   Amps. File: ',FILENAME,' "'
!      WRITE(10,'(a)') '@G0 on'
!      WRITE(10,'(a)') '@with g0'
!      OPEN (UNIT=11,FILE='rays2.dat') 
!      WRITE(11,'(a,f5.1,a,f5.1,a,a,a)') '@    subtitle "I1:',I1,',I2:',I2,'   Amps. File: ',FILENAME,' "'
!      WRITE(10,'(a)') '@G0 on'
!      WRITE(10,'(a)') '@with g0'
!      OPEN (UNIT=12,FILE='rays3.dat') 
      OPEN (UNIT=56,FILE='output.dat')
      
!***********************************
! draw outlines of the solenoids
!***********************************

!      WRITE(10,'(a)') '@type xy'
!      WRITE(10,*) Z1,0.0
!      WRITE(10,*) Z1,ER
!      WRITE(10,*) Z1+S1,ER
!      WRITE(10,*) Z1+S1,0.0
!      WRITE(10,'(a)') '&'
!      WRITE(10,'(a)') '@type xy'
!      WRITE(10,*) ZTGT-0.5*LENGTH,0.0
!      WRITE(10,*) ZTGT-0.5*LENGTH,BORE
!      WRITE(10,*) ZTGT+0.5*LENGTH,BORE
!      WRITE(10,*) ZTGT+0.5*LENGTH,0
!      WRITE(10,'(a)') '&'      

      IF (SOLNAME1 .EQ. 2) THEN
!         WRITE(10,'(a)')'@type xy'
!         WRITE(10,*) Z2,0
!         WRITE(10,*) Z2,ER
!         WRITE(10,*) Z2+S1,ER
!         WRITE(10,*) Z2+S1,0
!         WRITE(10,'(a)')'&'
!         WRITE(10,*)'@type xy'
!         WRITE(10,*) ZTGT+SEPAR-0.5*LENGTH,0.0
!         WRITE(10,*) ZTGT+SEPAR-0.5*LENGTH,BORE
!         WRITE(10,*) ZTGT+SEPAR+0.5*LENGTH,BORE
!         WRITE(10,*) ZTGT+SEPAR+0.5*LENGTH,0
!         WRITE(10,'(a)')'&'
      ENDIF
      
!***********************************
! draw outline of electrode
!***********************************
      IF (E_RAD .EQ. 1) THEN
!         WRITE(10,'(a)')'@type xy'
!         WRITE(10,*) Z_E_1,R_INNER
!         WRITE(10,*) Z_E_2,R_INNER
!         WRITE(10,'(a)')'&'
!         WRITE(10,'(a)')'@type xy'
!         WRITE(10,*) Z_E_1,R_OUTER
!         WRITE(10,*) Z_E_2,R_OUTER         
!         WRITE(10,'(a)')'&'
      ENDIF
      
!***********************************
! draw outline of mid-plane absorber
!***********************************
      IF ((ICHAR(MODE3) .EQ. 89) .OR. (ICHAR(MODE3) .EQ. 121)) THEN
!         WRITE(10,'(a)')'@type xy'
!         WRITE(10,*) POSITION_ABSORBER,0
!         WRITE(10,*) POSITION_ABSORBER,ER
!         WRITE(10,'(a)')'&'
      ENDIF

!***********************************
! draw outline of defocusing solenoid
!***********************************
      IF ((ICHAR(MODE4) .EQ. 89) .OR. (ICHAR(MODE4) .EQ. 121)) THEN
!         WRITE(10,'(a)')'@type xy'
!         WRITE(10,*) DS_Z1,0.0
!         WRITE(10,*) DS_Z1,DS_ER
!         WRITE(10,*) DS_Z1+DS_S1,DS_ER
!         WRITE(10,*) DS_Z1+DS_S1,0.0
!         WRITE(10,'(a)')'&' 
      ENDIF
      
!***********************************
! MAIN LOOP - done for each particle 
!***********************************
     
      DO I=1,NR
!         WRITE(10,'(a)')'@type xy'
!         WRITE(11,'(a)')'@type xy'
         IF(I.GT.20) THEN
            CLOSE(UNIT=9)
            STOP
         ENDIF
         N=ZO(I)
         M=AMU(I)
         ANG3=THETA(I)
!         THETA=THETA(I)
         E=EO(I)
         T=ANG3/57.3  
!   EMERGING ANGLE
      
         V0=1.389E+07*SQRT(E/M)
!****************************************
!  see the appendix for calculation of V0
!****************************************

         
!***********************************************************
!  initialize some variables - this is done for each new ray
!***********************************************************
!         X0=0.0
!         Y0=0.0
!         Z0=0.0
         X=X0
         Y=Y0
         Z=Z0
         VX=V0*SIN(T)
         VY=0.0
         VZ=V0*COS(T)
         GX=0.0
         GY=0.0
         GZ=0.0
         HGX=0.5*GX
         HGY=0.5*GY
         HGZ=0.5*GZ
         F7=0.0
         F8=0.0
         F9=0.0
         HF7=0.5*F7
         HF8=0.5*F8
         HF9=0.5*F9
          
!**********************************************************************
!   CALCULATIONS OF POSITIONS AND FIELDS STARTS HERE
!   THIS HALF STEP TECHNIQUE ALLOWS LARGE STEPS - 
!   Note that all the field calculations are done at
!   the halfway point through the interval.
!**********************************************************************
!         WRITE(6,'(a,I1)') 'Now calculating the trajectory of ray #',I
!         WRITE(6,*) 'Now calculating the trajectory of ray #'
         IF ((ICHAR(MODE3) .EQ. 89) .OR. (ICHAR(MODE3) .EQ. 121)) THEN
            ION = INT(N) 
            M1=AMU(I)
            CALL GETCOEF(ION,MM1,M1,ABSORBER,M2,RHO,ATRHO,VFERMI,LFCTR,PCOEF)
         ENDIF
         IF ((ICHAR(MODE5) .EQ. 89) .OR. (ICHAR(MODE5) .EQ. 121)) THEN
            ION = INT(N) 
            M1=AMU(I)
            CALL GETCOEF(ION,MM1,M1,ABSORBER,M2,RHO,ATRHO,VFERMI,LFCTR,PCOEF)
         ENDIF
         DO STEP = 1,STEPNUM
            V1=VX+HF7
            V2=VY+HF8
            V3=VZ+HF9
            DT=DZ/V3
            DTS=DT*DT
            XX=X+HGX
            YY=Y+HGY
            ZZ=Z+HDZ
            XXS=XX*XX
            YYS=YY*YY
            RRS=XXS+YYS
            RR=SQRT(RRS)
            INRRS=1/(RRS+0.00000001)
            INRR=1/(RR+0.00000001) 
           
!***********************************************
!  This next section uses rationalized B fields.
!***********************************************
            BZSL=0.
            BRSL=0.
            TEMP=1.
            CALL B_CALC(MA1,ERS,ZZ,CC,LL,RR,BZSL,BRSL,TEMP)
            IF (SOLNAME1 .EQ. 2) THEN   
               CALL B_CALC(MA2,ERS,ZZ,CCC,LL,RR,BZSL,BRSL,PAR_SOLENOID)
               !*** for second solenoid  ***!
            ENDIF
            IF ((ICHAR(MODE4) .EQ. 89) .OR. (ICHAR(MODE4) .EQ. 121)) THEN    
               CALL B_CALC(DS_MA,DS_ERS,ZZ,DS_CC,DS_LENGTH,RR,BZSL,BRSL,TEMP) 
               !*** for defocusing element  ***!        
            ENDIF
            BR=BRSL
            BZ=BZSL
            BX=BR*INRR*XX
            BY=BR*INRR*YY
            zzcm=zz*100.
            write(7,*)zzcm,brsl,bzsl                                       
!*********************************
!   Next section only for radial electrode
!*********************************
            IF (E_RAD .EQ. 1) THEN 
            	IF ((Z .GE. Z_E_1) .AND. (Z .LE. Z_E_2)) THEN
            		IF ((RR .LE. R_INNER) .OR. (RR .GE. R_OUTER)) THEN
            		  WRITE(6,'(a,I2,a)')'RAY #',I,' HIT THE ELECTRODE !!!'
            			EXIT
            		ENDIF
                    E_FIELD = (V_INNER-V_OUTER) / Log(R_OUTER/R_INNER) * INRR
                    EX = E_FIELD*INRR*XX
                    EY = E_FIELD*INRR*YY
                    EZ = 0
                ELSE 
                	EX = 0
                	EY = 0
                	EZ = 0
                ENDIF
            ENDIF
! *****************************
 

!**************************
! see appendix regarding CQM
!***************************
            CQM=9.648E+07*(N/M)


!**********************************************************************
!  calculate the cross products V x B 
!  and also include the electrostatic field
!**********************************************************************
            AX=CQM*(V2*BZ-V3*BY+EX) 
            AY=CQM*(V3*BX-V1*BZ+EY)
            AZ=CQM*(V1*BY-V2*BX+EZ)
      
!**********************************************************************
!  calculate the X and Y displacements and new X and Y positions
!**********************************************************************
            GX=(V1*DT+.5*AX*DTS)
            GY=(V2*DT+.5*AY*DTS)
            X=X+GX
            Y=Y+GY
            Z=DZ*FLOAT(STEP)         
!**********************************************************************
!  calculate the increments of the velocity 
!  calculate the half steps to be used in next iteration of loop
!**********************************************************************

            F7=AX*DT
            F8=AY*DT
            F9=AZ*DT
            HF7=.5*F7
            HF8=.5*F8
            HF9=.5*F9
            VX=VX+F7
            VY=VY+F8
            VZ=VZ+F9
             npth=stepnum/75
            IF (MOD(STEP,npth) .EQ. 0) THEN
!                         Write only every 20th point to disk    
!**********************************************************************
!  The projections on R (radial position) vs Z (axial position)
!  will be used to actually display the rays on the graphical output.
!  The horizontal axis is the axial position and the vertical axis is 
!  the radial position.  Thus:
!**********************************************************************
               XP=Z
               YP=SQRT(X**2+Y**2) 
               XX=(YP*COS(T))*100
               YY=(YP*SIN(T))*100
               xxr=X*COS(phy)-Y*Sin(phy)
               yyr=X*SIN(phy)+Y*COS(phy)
!               WRITE(10,*)XP,XX,YY
!               WRITE(10,*) XP,YP
!               WRITE(11,*) GX,GY
               WRITE(10,*) XP,xxr,yyr
               WRITE(56,*) E 
            ENDIF
 
            
 !*********************************************************************
 !  this next section is for the mid-plane absorber
 !  the algorithm was adapted from Ziegler's code found in
 !  the back of range tables series, vol 1.
 !*********************************************************************        
            IF ((ICHAR(MODE3) .EQ. 89) .OR. (ICHAR(MODE3) .EQ. 121)) THEN
            IF (STEP .EQ. POS_ABSORBER) THEN
               TEMP = (SQRT(VX*VX+VY*VY+VZ*VZ)/VZ) * THICKNESS_ABSORBER
               J_MAX = TEMP/.001
               !***************************************************
               ! The actual thickness encountered by the ion depends
               ! on the angle of incidence. Dividing by a .001 implies
               ! we are taking steps of microns of thickness.
               !***************************************************
               DO J = 1,J_MAX
                  VELSQR = E*1000/M1
                  !******
                  !calculate universal nuclear stopping power: SN
                  ! in units of MeV/micron
                  !******
                  EPSIL=32.53*M2*E*1000/(N*ABSORBER*(M1+M2)*(N**.23+ABSORBER**.23))
                  IF (EPSIL .GE. 30.) THEN
                     SN = LOG(EPSIL)/(2*EPSIL)
                  ELSE 
                     TEMP = (.01321*EPSIL**.21226) + (.19593*EPSIL**.5)
                     SN = .5*LOG(1+1.1383*EPSIL)/(EPSIL+TEMP)
                  ENDIF
                  SN = SN*N*ABSORBER*M1*8.462/((M1+M2)*(N**.23+ABSORBER**.23))
                  SN = SN * ATRHO * .001  !** convert to MeV/micron  **!
                  IF (ION .EQ. 1) CALL PSTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
                  IF (ION .EQ. 2) CALL HESTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
                  IF (ION .GT. 2) CALL HISTOP(ION,M1,ABSORBER,M2,VELSQR,E,VFERMI,LFCTR,PCOEF,SE)
                  SE = SE * ATRHO * .001  !** convert to MeV/micron  **!
                  E = E - SE - SN
                  IF (E .LT. .000002) E = .000002
               END DO
               TEMP = E/EO(I)  !** Final Energy / Initial Energy **!
               TEMP = SQRT(TEMP)
               VX = TEMP*VX
               VY = TEMP*VY
               VZ = TEMP*VZ  
!               N = N + 1.0 
               WRITE(*,*) "Energy =", E      
            ENDIF
            ENDIF
            
!           ***********
!           gas mode
!           unlike the mid-plane absorber which is a one time calculation done at the
!           location of the absorber this section is calculated every loop while
!           we are inside the first solenoid.
!           ***********
            IF ((ICHAR(MODE5) .EQ. 89) .OR. (ICHAR(MODE5) .EQ. 121)) THEN
            IF ( (Z .GE. (ZTGT-LENGTH/2)) .AND. (Z .LE. (ZTGT+LENGTH/2)) ) THEN
               TEMP = (SQRT(VX*VX+VY*VY+VZ*VZ)/VZ) * DZ/.000001
               !***************************************************
               ! The actual thickness of gas encountered by the ion depends
               ! on the angle of incidence. Dividing by  .000001 
               ! gives us the gas thickness in microns.
               !***************************************************
               ATRHO = ATRHO * TORR / 760
               VELSQR = E*1000/M1
               !******
               !calculate universal nuclear stopping power: SN
               ! in units of MeV/micron
               !******
               EPSIL=32.53*M2*E*1000/(N*ABSORBER*(M1+M2)*(N**.23+ABSORBER**.23))
               IF (EPSIL .GE. 30.) THEN
                  SN = LOG(EPSIL)/(2*EPSIL)
               ELSE 
                  TEMP = (.01321*EPSIL**.21226) + (.19593*EPSIL**.5)
                  SN = .5*LOG(1+1.1383*EPSIL)/(EPSIL+TEMP)
               ENDIF
               SN = SN*N*ABSORBER*M1*8.462/((M1+M2)*(N**.23+ABSORBER**.23))
               SN = SN * ATRHO * .001  !** convert to MeV/micron  **!
               SN = SN * TEMP          !** multiply by gas thickness **!
               IF (ION .EQ. 1) CALL PSTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
               IF (ION .EQ. 2) CALL HESTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
               IF (ION .GT. 2) CALL HISTOP(ION,M1,ABSORBER,M2,VELSQR,E,VFERMI,LFCTR,PCOEF,SE)
               SE = SE * ATRHO * .001  !** convert to MeV/micron  **!
               SE = SE * TEMP          !** multiply by gas thickness **!
               E = E - SE - SN
               IF (E .LT. .000002) E = .000002
               TEMP = E/EO(I)  !** Final Energy / Initial Energy **!
               TEMP = SQRT(TEMP)
               VX = TEMP*VX
               VY = TEMP*VY
               VZ = TEMP*VZ  
!               WRITE(*,*) "Energy =", E          
            ENDIF
            ENDIF
            
         END DO   !*** calculate next step ***!
!      WRITE(10,'(a)')'&'
!      WRITE(11,'(a)')'&'
      END DO      !*** calculate next ion ***!
      CLOSE(UNIT=10)
!      CLOSE(UNIT=11)
!      CLOSE(UNIT=12)
      
      CONTAINS
            
      SUBROUTINE SOL_INPUT_VAL(SOLNAME1,I1,I2,S1,ER,MA1,MA2,LENGTH,BORE,Z1,Z2,ZTGT,SEPAR,MODE1,PAR_SOLENOID)
      !******************************************************************
      !  this subroutine initializes the various solenoid parameters after
      !  requesting from the user which solenoids and settings to use
      !******************************************************************
      REAL I1,I2,S1,ER,MA1,MA2,LENGTH,BORE,Z1,Z2,ZTGT,SEPAR,PAR_SOLENOID
      INTEGER SOLNAME1,SOLNAME2
      CHARACTER *1 MODE1
         IF (SOLNAME1 .EQ. 1) THEN
        !   WRITE(6,*) 'Which solenoid do you want?'
        !   WRITE(6,*) 'UND Lil SOL,RIBRAS, MSU BIGSOL (enter "1", "2=RIBRAS" or "3")'
        !   WRITE(6,*) 'If you dont choose 1, 2 or 3, I will default to 2.' 
           READ (5,*) SOLNAME2
           if(solname2.ne.2)then
           if(solname2.ne.1)then
           if(solname2.ne.3)then
            solname2=2 
           endif
           endif
           endif      
        !   WRITE (6,*) 'Max currents: Lil SOL[125A], RIBRAS[91.88], NuSOL[100A]'
        !   WRITE (6,*) 'Enter the current in amps  (eg. 80.0) = '
           READ  (5,*) I1
           write(11,*)I1
         !************************************
         ! values for L'il SOL @ UND
         !************************************
            IF ( SOLNAME2 .EQ. 1) THEN 
               S1=0.35
               ER=0.1205
               MA1=2.125*I1/125.
               WRITE(6,*) 'Solenoid center to back position:    0.523  m'
               WRITE(6,*) 'Back to center (disk radius):        0.0967 m'
               LENGTH = 1.04
               BORE = 0.097

         !************************************
         ! values for RIBRAS
         !************************************
            ELSE IF ( SOLNAME2 .EQ. 2) THEN 
                S1=0.679958
               ER=0.192034
!               MA1=3.636*I1/100
               MA1=3.74636*I1/91.86
          !     WRITE(6,*) 'Coil Length        0.679958 m'
          !     WRITE(6,*) 'Bore Radius        0.15 m'
               LENGTH = 1.00
               BORE = 0.15

         !************************************
         ! values for BIGSOL
         ! Inner Diameter (coil) = 19.50 inches --> R-inner = 0.24765 m 
         ! Outer Diameter (coil) = 25.75 inches --> R-outer = 0.327025 m
         ! RMS Radius = Sqr[(R-inner^2 + R-outer^2)/2] = 0.2900654 m
         ! Warm Bore Diameter = 17.75 inches --> Bore Radius = 0.225425 m
         ! Coil Length = 37.3333 inches --> S1  = 0.948266666 m
         ! Bore Length = 60.32 inches --> LENGTH = 1.532128 m
         ! Max Specs: 6.5 Tesla at 171.25 amps
         !************************************
            ELSE IF ( SOLNAME2 .EQ. 3) THEN  
               S1=0.948
               ER=0.29
               MA1=3.80995*I1/171.25
               WRITE(6,*) 'Coil Length        0.94827 m'
               WRITE(6,*) 'Bore Radius        0.22543 m'
               LENGTH = 1.532
               BORE = 0.225
         !************************************
         ! if you don't choose 1,2,3
         !************************************
             ELSE
         WRITE(6,*) 'Your choice is invalid.  I will assume you want RIBRAS.'
!               S1=.35
!               ER=.1205
!               MA1=2.125*I1/125.          
               S1=0.679958
               ER=0.1920342
               MA1=3.74636*I1/91.86
               WRITE(6,*) 'Coil Length        0.679958 m'
               WRITE(6,*) 'Bore Radius        0.15 m'
               LENGTH = 1.00
               BORE = 0.15
!               WRITE(6,*) 'Solenoid center to back position:    0.523  m'
!               WRITE(6,*) 'Back to center (disk radius):        0.0967 m'
            ENDIF
         Z1 = ZTGT - 0.5*S1   
         ENDIF  
         !*************************************
         !  For Two Solenoids - Values Corrected 6/25/96!!!
         !*************************************
         IF (SOLNAME1 .EQ. 2) THEN
            WRITE (6,*)'Enter current in amps for RIBRAS1: '
            READ  (5,*) I1
            WRITE (6,*)'Enter current in amps for RIBRAS2: '
            READ  (5,*) I2
            WRITE (6,*) 'Enter the center to center separation of NuSOLS (2.0m min.): '
            READ  (5,*) SEPAR
!            S1 = 0.80
!            ER = 0.2123
!            LENGTH = 1.0
!            BORE = .15
!            MA1 = 3.636*I1/100 
!            MA2 = 3.636*I2/100
               S1=0.680008
               ER=0.19204235
               Bore=0.15
               MA1=3.74636*I1/91.86
               MA2=3.74610*I2/91.88
               WRITE(6,*) 'Coil Length        0.679958 m'
               WRITE(6,*) 'Bore Radius        0.15 m'
               LENGTH = 1.00
               BORE = 0.15
            WRITE(6,*) 'Solenoids in parallel or antiparallel mode? '
            WRITE(6,*) 'Enter "p" or "a" '
            READ (5,*) MODE1
            IF (ICHAR(MODE1) .LT. ICHAR('n')) THEN 
               PAR_SOLENOID = -1.0
            ELSE
               PAR_SOLENOID = 1.0
         ENDIF
         Z1 = ZTGT - 0.5*S1
         Z2 = ZTGT + SEPAR - 0.5*S1 
      ENDIF
      END SUBROUTINE
      
      SUBROUTINE B_CALC(MA,ERS,ZZ,CC,LL,RR,BZSL,BRSL,PAR_SOLENOID)
          !***********************
          ! calculate the radial and axial components of the magnetic
          ! field using rationalized expressions from Liu's original
          ! code.
          !***********************
          REAL MA,ERS,ZZ,CC,LL,RR,BZSL,BRSL,PAR_SOLENOID

          REAL MAERS,ZC,ZCS,ZL,ZLS,SCS,SC,SLS,SL,SC3,SC5,SC7
          REAL SL3,SL5,SL7,HR,HRS,HR3,BZZ,B1,B3,BRR
         
          MAERS=MA*ERS
          ZC=ZZ-CC
          ZCS=ZC*ZC
          ZL=ZZ-CC-LL
          ZLS=ZL*ZL
          SCS=ERS+ZCS
          SC=SQRT(SCS)
          SLS=ERS+ZLS
          SL=SQRT(SLS)
          SC3=SCS*SC
          SC5=SC3*SCS
          SC7=SC5*SCS
          SL3=SLS*SL
          SL5=SL3*SLS
          SL7=SL5*SLS
          HR=.5*RR
          HRS=HR*HR
          HR3=HRS*HR
          BZZ=MA*(ZC/SC-ZL/SL) - HRS*3*MAERS*(-ZC/SC5+ZL/SL5)
          B1=-HR*MAERS*(1/SC3-1/SL3)
          B3=HR3/2*3*MAERS*((4*ZCS-ERS)/SC7-(4*ZLS-ERS)/SL7)
          BRR=B1+B3
          BZSL=BZSL + PAR_SOLENOID*BZZ
          BRSL=BRSL + PAR_SOLENOID*BRR
!         *******************************************************
!           Differences between the first and second solenoids when 
!           calling this subroutine.
!           The variable on the left is the one found inside the 
!           subroutine.  The two variables on the right are the ones 
!           that are sent to the subroutine from the main program.
!         
!                              1st Solenoid        2nd Solenoid
!          
!            MA            =       MA1                 MA2
!            CC            =       CC                  CCC
!            PAR_SOLENOID  =        1              PAR_SOLENOID
!         
!            Before the B_CALC subroutine is called for the first 
!            solenoid, make sure that BZSL and BRSL have been 
!            initialzed to zero. Also when calling the subroutine 
!            for the first solenoid, the variable PAR_SOLENOID 
!            receives zero as its value.  PAR_SOLENOID specifies
!            whether the two solenoids are in parallel or antiparallel
!            configurations.  Since solenoids are always focusing, 
!            the absolute orientations of the two magnetic fields 
!            is unimportant, only their relative orientations. Note 
!            that having the two solenoids parallel will reduce the 
!            strength of the fringe fields in the space between the 
!            them.  Also, there is slightly greater focusing power 
!            when the two solenoids are parallel.
!         *******************************************************
       END SUBROUTINE
             
         SUBROUTINE GETCOEF(ION,MM1,M1,ABSORBER,M2,RHO,ATRHO,VFERMI,LFCTR,PCOEF)
            INTEGER ION,ABSORBER,MM1,AAAA,BBBB,IIII
            REAL M1,M2,RHO,ATRHO,VFERMI,LFCTR
            REAL CCCC,DDDD,EEEE,FFFF,GGGG,HHHH
            REAL, DIMENSION(8) :: XXXX,PCOEF
!            INTEGER :: io     
            OPEN(20, FILE = 'scoef_unix.dat')
            DO IIII = 1,92
!               READ(20,*,IOSTAT=io)AAAA,BBBB,CCCC,DDDD,EEEE,FFFF,GGGG,HHHH      
            READ(20,*)AAAA,BBBB,CCCC,DDDD,EEEE,FFFF,GGGG,HHHH          
               IF (IIII .EQ. ION) THEN
                  LFCTR=HHHH
               ENDIF
               IF (IIII .EQ. ABSORBER) THEN
                  M2=DDDD
                  RHO=EEEE
                  ATRHO=FFFF
                  VFERMI=GGGG
               ENDIF
            END DO     
            DO IIII = 1,92 
!               READ(20,*,IOSTAT=io)AAAA,XXXX
            READ(20,*)AAAA,XXXX
               IF (IIII .EQ. ABSORBER) THEN
                 PCOEF = XXXX
               END IF
            END DO
            CLOSE(20)
            !*****************************
            ! M1, MM1 and LFCTR refer to the ion
            ! M2, RHO, ATRHO, VFERMI refer to the absorber
            !*****************************
         END SUBROUTINE
      
         SUBROUTINE PSTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
            INTEGER ION,ABSORBER
            REAL M1,M2,VELSQR,SE
            REAL, DIMENSION(8) :: PCOEF
            REAL PE, SL, SH, VELPWR
            PE = MAX(25.,VELSQR)
            SL=( PCOEF(1) * PE**PCOEF(2) ) + ( PCOEF(3) * PE**PCOEF(4) )
            SH=(PCOEF(5)/(PE**PCOEF(6))) * LOG((PCOEF(7)/PE) + PCOEF(8)*PE)
            SE = SL*SH/(SL+SH)
            IF (VELSQR .LT. 25.) THEN
               VELPWR = 0.45
               IF (ABSORBER .LE. 6) VELPWR=0.25
               SE = SE*(VELSQR/25.)**VELPWR
            ENDIF              
         END SUBROUTINE
         
         SUBROUTINE HESTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SE)
            INTEGER ION,ABSORBER
            REAL M1,M2,VELSQR,SE
            REAL, DIMENSION(8) :: PCOEF         
            REAL A, B, HE, HEH, SP
            HE = MAX(1.,VELSQR)
            B = LOG(HE)
            A = .2865+.1266*B-.001429*B*B+.02402*B*B*B-.01135*B**4+.001475*B**5
            HEH = 1. - EXP(-MIN(30.,A))
            A = (1.+(.007+.00005*ABSORBER)*EXP(-(7.6-MAX(0.,LOG(HE)))**2))
            HEH = A*A*HEH
            CALL PSTOP(ION,M1,ABSORBER,M2,HE,PCOEF,SP)
            SE = SP * HEH * ION * ION
            IF (VELSQR .LT. 1.) SE = SE*SQRT(VELSQR/1.)
         END SUBROUTINE
         
         SUBROUTINE HISTOP(ION,M1,ABSORBER,M2,VELSQR,E,VFERMI,LFCTR,PCOEF,SE)
            INTEGER ION,ABSORBER
            REAL M1,M2,VELSQR,E,VFERMI,LFCTR,SE,SP
            REAL, DIMENSION(8) :: PCOEF        
            REAL VRMIN, YRMIN, VVVV, VR, YR,EEE,POWER,VMIN
            REAL AAAA,BBBB,QQQQ,QQQ1,LLL1,LLL0,LLLL,ZETA
            YRMIN = .13
            VRMIN = 1.0
            VVVV = SQRT(VELSQR/25)/VFERMI
            IF (VVVV .GE. 1.) THEN
               VR = VVVV*VFERMI*(1+1/(5*VVVV*VVVV))
            ELSE
               VR = (3*VFERMI/4)*(1+(2*VVVV*VVVV/3)-VVVV**4/15)
            ENDIF
            YR = MAX(YRMIN,VR/ION**.6667,VRMIN/ION**.6667)
            AAAA=-.803*YR**0.3+1.3167*YR**0.6+.38157*YR+.008983*YR*YR
            QQQQ=MIN(1.,MAX(0.,1.-EXP(-MIN(AAAA,50.))))
            BBBB=(MIN(0.43,MAX(.32,.12+.025*ION)))/ION**.3333
            LLL0=(.8-QQQQ*(MIN(1.2,0.6+ION/30.)))/ION**.3333
            LLL1=BBBB*(1.-QQQQ)/(.025*MIN(16.,1.*ION))
            IF (QQQQ .LT. 0.2) LLL1 = 0.
            IF (QQQQ .LT. (MAX(0.,.9-.025*ION))) THEN
               QQQ1=0.2
               LLL1=BBBB*(QQQQ-.2)/ABS(MAX(0.,.9-.025*ION)-.2000001)
            ENDIF
            IF (QQQQ .LT. (MAX(0.,1.-.025*MIN(16.,1.*ION)))) LLL1 = BBBB
            LLLL = MAX(LLL1,LLL0*LFCTR)
            ZETA = QQQQ+(1./(2.*VFERMI**2))*(1.-QQQQ)*LOG(1+(4*LLLL*VFERMI/1.919)**2)
            AAAA=-(7.6-MAX(0.,LOG(VELSQR)))**2
            ZETA=ZETA*(1.+(1./ION**2)*(.18+.0015*ABSORBER)*EXP(AAAA))
            IF (YR .LE. MAX(YRMIN,VRMIN/ION**.6667)) THEN
               VRMIN=MAX(VRMIN,YRMIN*ION**.6667)
               VMIN=.5*(VRMIN+SQRT(MAX(0.,VRMIN**2-0.8*VFERMI**2)))
               EEE = 25*VMIN**2
               CALL PSTOP(ION,M1,ABSORBER,M2,EEE,PCOEF,SP)
               POWER=.5
               IF ((ABSORBER.EQ.6).OR.(((ABSORBER.EQ.14).OR.(ABSORBER.EQ.32)).AND.(ION.LE.19))) POWER = .375
               SE=(SP*(ZETA*ION)**2)*(VELSQR/EEE)**POWER               
            ELSE
               CALL PSTOP(ION,M1,ABSORBER,M2,VELSQR,PCOEF,SP)
               SE=SP*(ZETA*ION)**2
            ENDIF           
         END SUBROUTINE         

      END



!**********************************************************************
!    ********    APPENDIX     *********          
!**********************************************************************
!
!*****************************
!   MA is derived empirically
!*****************************
!   For example,
!           BIGSOL has a 6.5T field at 171.25 Amps.  MA has the
!           form:   MA = constant x (current / max current) for the 
!                   the superconducting solenoids.  For a defocusing
!                   permanent magnet it has the form MA = constant
!           The constant can be determined as follows:
!          
!         MA     x    LL       where LL is the effective length of
! BM = ---------------------   the solenoid and ER is the effective
!      SQRT[ER^2 + (LL/2)^2]   radius (RMS radius). BM is the max
!                              field of the solenoid along the axis
!                              (in Tesla).
!
! therefore MA = BM x SQRT[ER^2 + (LL/2)^2] / LL
! units of MA are Tesla 
!
!*****************************************************************
!    for a finite single sheet solenoid the magnetic field along the axis is:
!
!       B = 0.5 * mu_zero * n * I * [ cos(theta1) - cos(180-theta2) ]   
!
!    where:
!          n is the number of turns per unit length
!          I is in Amperes
!          theta1 and theta2 are the half of the vertex angles
!                 subtended by the two openings of the solenoid
!
!                  ..............
!                   \        _-
!                    \     _-
!            theta2   \  _-     theta1
!              --------()------------------   z-axis
!
!
!
!                  xxxxxxxxxxxxxx
!
!     In our program the various physical constants and numerical
!     factors are included in the 'constant' found in the variable MA.
!

!******************************
!  V0 is calculated as follows                                
!******************************                                                                  
!  1 AMU = 1.66054 E-27 kg   or   931.494 MeV/c^2                  
!  1 MeV/nucleon = (c^2)/931.494                                   
!                = (2.99792E+08 m/s)^2 / 931.494                   
!                = 9.64850E+13 (m/s)^2                             
!  V = sqrt (2E/M)  ---> V = 1.389E+07*SQRT(E/M)                   
!
!****************************
!  calculation CQM
!****************************
!   NOTE: Liu's and ND update by J. Meissner had 9.651E-07 
!         for the constant in CQM. I calculate CQM as follows:
!  e    = 1.60219E-19 coulombs
!  u    = 1.66057E-27 kg
!  e/u  = 9.64843E+07
!
!******************************
!   Check out the postscript file solray.f90ps
!   for a diagram of how the fields are calculated
!******************************
!
!   A thin sheet approximation was used for the field calculations.
!	An expansion is used to calculate the magnetic field off axis.
!
!****************************
!  approximation of Electric Field   (01/11/96)
!****************************
!	For this preliminary calculation we will neglect fringe
!	effects and use the electric field for an infinitely long
!	cylindrical electrode and only take into effect the electric
!	field as the particle passes directly over the electrode.
!
!	For infinitely long concentric cylindrical electrode where
!	R1 = R_INNER, V1 = V_INNER, R2 = R_OUTER. V2 = V_OUTER
!	
!	V(r) = V1 + {  [ (V2 - V1) / ln(R2/R1) ] * ln(r/R1)  }
!	
!	E(r) == - dV/dr = [ (V1-V2) / ln(R2/R1)]  *  1/r   
!
!
!****************************
!  straggling
!****************************
!   in Enge's Nuclear Physics Book (Addison-Wesley, Reading MA)
!   Ch 7 (Stopping and Detecting Nuclear Radiations) he refers to
!   an article by N. Bohr in Mat Fys Medd Dan Vid Selsk vol18, no 8 (1948)
!   that give the standard deviation of the straggling as
!   8.85 * z * sqr(2Z/A) * sqr(areal density) --> in units of keV
!   z is the projectile z and Z is the target Z.  For now (4/14) I'll
!   use this simple gaussian distribution for the straggling.  However,
!   given that it's an article from 1948 I suspect it's not very
!   universally applicable.
!******************************************

