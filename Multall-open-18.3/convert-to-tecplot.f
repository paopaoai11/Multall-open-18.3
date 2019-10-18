C
C     THIS PROGRAM COMVERTS A MULTALL PLOTTING OUTPUT FILE TO A "Tecplot" INPUT FILE.
C     IT READS IN THE FILES  "flow_out"  and "grid_out"  written by MULTALL AND WRITES 
C     A FILE "tecplot-input.dat" which should be readable by "tecplot" .
C
      PARAMETER(NI=95,NJ=2500,NK=95)
C
      COMMON/BLKINT/  IM,JM,KM,NSTEP
C
      COMMON/BLKPRIM/ RO(NI,NJ,NK),ROVX(NI,NJ,NK),ROVR(NI,NJ,NK),
     &                ROVT(NI,NJ,NK),ROE(NI,NJ,NK),DUM(NI,NJ,NK),
     &                Q(NI,NJ,NK),QQ(NI,NJ,NK)
C
      COMMON/BLKFLOW/ VX(NI,NJ,NK),VR(NI,NJ,NK),VT(NI,NJ,NK),
     &                WT(NI,NJ,NK),PSTAT(NI,NJ,NK),MACH(NI,NJ,NK),
     &                TSTAT(NI,NJ,NK),TSTAG(NI,NJ,NK),PSTAG(NI,NJ,NK),
     &                ENT_FUN(NI,NJ,NK),VMER(NI,NJ,NK),YAW(NI,NJ,NK),
     &                PITCH(NI,NJ,NK)
C
      COMMON/BLKGEOM/ IND(NJ),W(NJ),NBLADES(NJ),CP,CV,RGAS,GA,PI,
     &                X(NJ,NK),R(NJ,NK),RTHETA(NI,NJ,NK),THETA(NI,NJ,NK)
C
      REAL  MACH
C
      CHARACTER*1  ANS
      CHARACTER*72 FLOWFILE, GRIDFILE, TITLE
C
C
       PI     = 3.14159265
C
C
      OPEN  (UNIT=1, FILE= '/dev/tty')
      OPEN  (UNIT=10,FILE= 'tecplot-input.dat')
C
      TITLE = ' TECPLOT INPUT DATA FROM MULTALL'
C
C     ANNOUNCE PROGRAM
C
      WRITE (6,*) ' PROGRAM TO CONVERT MULTALL OUTPUT DATA TO A TECPLOT
     & INPUT FILE'
C
C     Find out the input file name
C
      FLOWFILE   = 'flow_out'
      WRITE(6,*) 'The flow file is usually named flow_out, is this OK ?'
      WRITE(6,*) 'Answer Y or N '
      READ(1,*) ANS
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
           WRITE(6,*) ' INPUT THE NAME OF THE INPUT FILE'
           READ(1,10) FLOWFILE
      END IF
   10 FORMAT(A80)
C
      WRITE(6,*) ' Input file = ' ,FLOWFILE
C
      OPEN(UNIT=7,FILE= FLOWFILE,   FORM='unformatted')
C
C
      GRIDFILE = 'grid_out'
      WRITE(6,*) 'The grid file is usually named grid_out, is this OK ?'
      WRITE(6,*) 'Answer Y or N '
      READ(1,*) ANS
      IF(ANS.EQ.'N'.OR.ANS.EQ.'n') THEN
           WRITE(6,*) ' INPUT THE NAME OF THE INPUT FILE'
           READ(1,10) GRIDFILE
      END IF
C
      WRITE(6,*)  ' The grid geometry file = ', GRIDFILE
C
      OPEN(UNIT=21,FILE= GRIDFILE, FORM= 'unformatted' )
C
C****************************************************************************
C****************************************************************************
C     CALL GRIDIN  TO INPUT THE GRID GEOMETRY FROM UNIT 21 .
C
      WRITE(6,*) ' calling GRIDIN '
      CALL GRIDIN 
      WRITE(6,*) ' done GRIDIN '
C
C****************************************************************************
C****************************************************************************
C
C   CALL FLOWIN  TO INPUT A NEW SET OF FLOW VARIABLES.
C
      WRITE(6,*) ' calling sub FLOWIN'
        CALL FLOWIN
      WRITE(6,*) ' done sub FLOWIN '
C
C****************************************************************************
C****************************************************************************
C
C     CALCULATE MORE FLOW PROPERTIES
C
      IMID  = IM/2
      KMID  = KM/2
      WRITE(6,*) ' IM   = ', IM  , 'JM   = ', JM,   'KM  = ', KM
      WRITE(6,*) ' IMID = ', IMID, 'KMID = ', KMID
C
      CV    = CP/GA
      RGAS  = CP - CV
      FGA   = GA/(GA-1.0)
      WRITE(6,*) ' CP = ', CP, 'CV= ', CV, 'GAMMA= ', GA, ' RGAS= ',RGAS
C
      DO 100 K=1,KM
      DO 100 J=1,JM
      DO 100 I=1,IM
      VX(I,J,K) = ROVX(I,J,K)/RO(I,J,K)
      VT(I,J,K) = ROVT(I,J,K)/RO(I,J,K)
      WT(I,J,K) = ROVT(I,J,K)/RO(I,J,K) - W(J)*R(J,K)
      VR(I,J,K) = ROVR(I,J,K)/RO(I,J,K)
      EINT      = ROE(I,J,K)/RO(I,J,K)
      VTABS     = VT(I,J,K)
      WTREL     = WT(I,J,K)
      VMSQ      = VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K)
      VMER(I,J,K)  = SQRT(VMSQ)
      YAW(I,J,K)   = ATAN(WT(I,J,K)/VMER(I,J,K))*180./PI
      PITCH(I,J,K) = ATAN(VR(I,J,K)/VX(I,J,K))*180./PI
      EKE       =  VMSQ + VTABS*VTABS 
      EKEREL    =  VMSQ + WTREL*WTREL 
      TSTAT(I,J,K)=  (EINT - 0.5*EKE)/CV
      TSTAG(I,J,K)= TSTAT(I,J,K) + 0.5*EKE/CP
      VSSQ        =  GA*RGAS*TSTAT(I,J,K)
      MACH(I,J,K) = SQRT(EKEREL/VSSQ)
      PSTAT(I,J,K)= RGAS*TSTAT(I,J,K)*RO(I,J,K)
      PSTAG(I,J,K)= PSTAT(I,J,K)*(TSTAG(I,J,K)/TSTAT(I,J,K))**FGA
  100 CONTINUE
C
      PREF = PSTAG(IMID,1,KMID)
      TREF = TSTAG(IMID,1,KMID)
C
      WRITE(6,*) ' PREF = ', PREF, ' TREF = ', TREF
C
C     CALCULATE THE ENTROPY FUNCTION
C
      DO 200 K= 1,KM
      DO 200 J= 1,JM
      DO 200 I= 1,IM
      DELTA_ENTPY = CP*ALOG(TSTAG(I,J,K)/TREF)
     &            - RGAS*ALOG(PSTAG(I,J,K)/PREF)
      ENT_FUN(I,J,K) = EXP(-DELTA_ENTPY/RGAS)
  200 CONTINUE
C
C    WRITE TECPLOT FILE
C
      WRITE(6,*)
      WRITE(6,*) ' THE FLOW IS THE SAME IN ALL PASSAGES BUT THE PLOT'
      WRITE(6,*) ' LOOKS MUCH BETTER WITH TWO OR MORE PASSAGES'
      WRITE(6,*)
      WRITE(6,*) ' INPUT THE NUMBER OF BLADE PASSAGES TO BE PLOTTED'
      WRITE(6,*)
C
      READ(1,*) NPASS    
C
C
      WRITE(6,*)
      WRITE(6,*)' WRITING THE FILE "tecplot-input.dat". '
      WRITE(6,*)' THIS MAY TAKE A FEW MINUTES FOR LARGE FILES.'
      WRITE(6,*)
C
      DO 500 NPAS = 1,NPASS 
      WRITE(10,11)  TITLE 
   11 FORMAT(' TITLE = "',A40, ' " ')
      WRITE(10,12)
   12 FORMAT('Variables = "X","Y","Z","VX","VY","WT","VZ","RHO","PSTAT",
     &"TSTAT","PSTAG","TSTAG","MACH","YAW","PITCH","ENT_FUN" ')
      WRITE(10,13) NPAS, IM, JM, KM
   13 FORMAT(' ZONE T = "PASSAGE NUMBER',I5,' " ',' I=',I5,' J=',I6,
     &        ' K=',I5)
C 
      DO 300 K=1,KM
      DO 300 J=1,JM
      DO 300 I=1,IM
      ANGL = RTHETA(I,J,K)/R(J,K) + 2*PI*(NPAS-1)/NBLADES(J)
      XVAL = X(J,K)
      YVAL = R(J,K)*SIN(ANGL)
      ZVAL = R(J,K)*COS(ANGL)
C
      WRITE(10,*)  XVAL,YVAL,ZVAL,VX(I,J,K),VT(I,J,K),WT(I,J,K),
     &             VR(I,J,K),RO(I,J,K),PSTAT(I,J,K),TSTAT(I,J,K),
     &             PSTAG(I,J,K),TSTAG(I,J,K),MACH(I,J,K),YAW(I,J,K),
     &             PITCH(I,J,K),ENT_FUN(I,J,K)
C
  300 CONTINUE
C
  500 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*) ' FILE "tecplot-input.dat"  SUCCESSFULLY WRITTEN. '
      WRITE(6,*)
C
      STOP
      END
C
C****************************************************************************
C****************************************************************************
C*******************************************************************C
C
C            SUBROUTINE TO INPUT THE FLOW VARIABLES
C
C******************************************************************C
C******************************************************************C
C
      SUBROUTINE FLOWIN
C
      PARAMETER(NI=95,NJ=2500,NK=95)
C
C
      COMMON/BLKINT/  IM,JM,KM,NSTEP
C
      COMMON/BLKPRIM/ RO(NI,NJ,NK),ROVX(NI,NJ,NK),ROVR(NI,NJ,NK),
     &                ROVT(NI,NJ,NK),ROE(NI,NJ,NK),DUM(NI,NJ,NK),
     &                Q(NI,NJ,NK),QQ(NI,NJ,NK)
C
      COMMON/BLKFLOW/ VX(NI,NJ,NK),VR(NI,NJ,NK),VT(NI,NJ,NK),
     &                WT(NI,NJ,NK),PSTAT(NI,NJ,NK),MACH(NI,NJ,NK),
     &                TSTAT(NI,NJ,NK),TSTAG(NI,NJ,NK),PSTAG(NI,NJ,NK),
     &                ENT_FUN(NI,NJ,NK),VMER(NI,NJ,NK),YAW(NI,NJ,NK),
     &                PITCH(NI,NJ,NK)
C
      COMMON/BLKGEOM/ IND(NJ),W(NJ),NBLADES(NJ),CP,CV,RGAS,GA,PI,
     &                X(NJ,NK),R(NJ,NK),RTHETA(NI,NJ,NK),THETA(NI,NJ,NK)
C

C
C*****************************************************************************
C*****************************************************************************
C      READ IN THE PLOTTING/RESTART FILE FROM UNIT 7 .
C
      WRITE(6,*)
      WRITE(6,*) ' STARTING TO READ IN THE FLOW FIELD DATA FROM UNIT 7.'
      WRITE(6,*)
      READ(7)  NSTEP   
      WRITE(6,*)  ' READ NSTEP'  
      WRITE(6,*)  ' NSTEP =', NSTEP
      READ(7) (((RO(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ RO'  
      READ(7) (((ROVX(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ ROVX' 
      READ(7) (((ROVR(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ ROVR'  
      READ(7) (((ROVT(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ ROVT' 
      READ(7) (((ROE(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ ROE'  
      READ(7) (((DUM(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ DUM'  
      READ(7) (((  Q(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ Q' 
      READ(7) ((( QQ(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(6,*) ' READ QQ' 
C
      WRITE(6,*) ' LEAVING SUBROUTINE FLOWIN, FLOW DATA INPUT OK. '
      RETURN
      END
C
C****************************************************************************
C****************************************************************************
C******************************************************************C
C
C         SUBROUTINE TO INPUT GRID COORDINATES 
C
C******************************************************************C
C******************************************************************C
C
      SUBROUTINE GRIDIN 
C
      PARAMETER(NI=95,NJ=2500,NK=95) 
C
      COMMON/BLKINT/  IM,JM,KM,NSTEP
C
      COMMON/BLKPRIM/ RO(NI,NJ,NK),ROVX(NI,NJ,NK),ROVR(NI,NJ,NK),
     &                ROVT(NI,NJ,NK),ROE(NI,NJ,NK),DUM(NI,NJ,NK),
     &                Q(NI,NJ,NK),QQ(NI,NJ,NK)
C
      COMMON/BLKFLOW/ VX(NI,NJ,NK),VR(NI,NJ,NK),VT(NI,NJ,NK),
     &                WT(NI,NJ,NK),PSTAT(NI,NJ,NK),MACH(NI,NJ,NK),
     &                TSTAT(NI,NJ,NK),TSTAG(NI,NJ,NK),PSTAG(NI,NJ,NK),
     &                ENT_FUN(NI,NJ,NK),VMER(NI,NJ,NK),YAW(NI,NJ,NK),
     &                PITCH(NI,NJ,NK)
C
C
      COMMON/BLKGEOM/ IND(NJ),W(NJ),NBLADES(NJ),CP,CV,RGAS,GA,PI,
     &                X(NJ,NK),R(NJ,NK),RTHETA(NI,NJ,NK),THETA(NI,NJ,NK)
C
C****************************************************************************
C****************************************************************************
      WRITE(6,*) ' IN SUBROUTINE  GRIDIN '
C****************************************************************************
C****************************************************************************
C
      WRITE(6,*) ' STARTING TO READ IN GRID GEOMETRY DATA FROM UNIT 21 '
C
      READ(21) NSTEPS
      WRITE(6,*) ' MAXIMUM NUMBER OF TIME STEPS = ', NSTEPS
      READ(21) IM,JM,KM
      WRITE(6,*) ' GRID POINT NUMBERS, IM, JM, KM = ',IM,JM,KM
      READ(21) CP,GA
      READ(21) (IND(J),J=1,JM)
      READ(21) (W(J),J=1,JM)
      READ(21) (NBLADES(J),J=1,JM)
C 
      DO 250 J=1,JM
      DO 250 K=1,KM 
      READ(21) X(J,K),R(J,K),(RTHETA(I,J,K),I=1,IM)
      
  250 CONTINUE
C
      WRITE(6,*) ' LEAVING SUBROUTINE GRIDIN, GRID GEOMETRY INPUT OK'
C
      CLOSE(21)
C
      RETURN
      END


