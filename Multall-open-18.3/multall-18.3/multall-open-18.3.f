C                                                                      C
C                                                                      C
C             THIS PROGRAM CALCULATES STEADY 3D FLOW THROUGH           C
C             MULTIPLE TURBOMACHINE BLADE ROWS                         C
C             BY SOLVING THE UNSTEADY CONTINUITY,MOMENTUM,             C
C             AND ENERGY EQUATIONS UNTIL A STEADY STATE IS REACHED.    C
C                                                                      C
C           THIS VERSION ALLOWS FOR MULTIPLE BLADE ROWS. THERE         C
C           IS NO INTRINSIC LIMIT ON THE NUMBER OF ROWS WHICH          C
C           IS DETERMINED BY THE PARAMETER  "NRS" .                    C
C           INCREASE THIS IF REQUIRED.                                 C
C                                                                      C
c                      DEVELOPMENT HISTORY                             C
C                                                                      C
C   The original develoment started about 1973 at CEGB Marchwood       C
C   This was for inviscid flow through a single blade row with         C
C   extremely coarse grids, typically 10x25x10 points.                 C
C                                                                      C
C   The early versions used the "opposed difference" numerical scheme  C
C   with cell centre storage on overlapping cells.                     C
C                                                                      C
C   It was extended to deal with two or more blade rows using a simple C 
C   plane model in 1979, oon after the author moved to Cambridge.      C
C                                                                      C
C   Around 1980 the scheme was changed to use cell corner storage of   C
C   the variables with non-overlapping cells and multigrid was         C
C   introduced. Both pf these gave major improvements in performance.  C
C                                                                      C
C   The first nmoves to include viscous effects (1985) was by a separate
C   boundary layer calculation with transpiration through the blade    C
C   surfaces to allow for the boundary layer blockage. This was        C
C   surprisingly useful for compressors.                               C
C                                                                      C
C   Around 1986 the first approach to including viscous forces         C
C   in the calculation used a skin friction coefficient and an         C
C   empircal  distribution of the viscous viscous stresses through the C
C   flow. This was soon extended to obtain the skin friction from wall C
C   functions and a simple mixing length model for the turbulent       C
C   viscosity. The viscous forces are included via a body force which  C
C   only needs to be upadted about every 5 time steps, giving          C
C   significant savings in run time.                                   C
C                                                                      C
C   Numerous minor improvements were included throughout the 1990's    C
C   especially to the mixing plane treatment, and the program was      C
C   widely used for multistage compressors and turbines.               C
C   turbines. Cooling flow and bleed flows were added together with    C
C   the pinched tip model for tip leakage and a model for shroud       C
C   leakage flows.                                                     C
C                                                                      C
C   A major development around 1998 was the change from the "opposed   C
C   difference" scheme to the "scree" scheme. This was simpler and     C
C   more accurate but the mixing plane treatment needed changing.      C
C                                                                      C
c   A major "tidying up" was performed in 2005 and some of the older   C
C   loss routines were removed or combined into a single subroutine    C
C   "LOSS" . An option to use real gas properties was added in 2006    C
c   and an allowance for surface roughness was added in 2008.          C
C                                                                      C
C   An option to use the Spalart-Allmaras turbulence model was added   C
C  in 2010 and at the same time an updated mixing length model "NEWLOS"C
C  was added. Both of these are full Navier-Stokes models whilst       C
C  the original routine "LOSS" is a thin shear layer approximation.    C
C                                                                      C
C   In 2014 the option to perform Q3D calculations on a blade-to blade C
C   stream surface was added together with a new solution algorith,    C
C   the "SSS" scheme, which permits larger CFL numbers.                C
C                                                                      C
C   The mixing plane model was further improved in 2015 to allow betterC
C   interation with shock waves and to pernmit reverse flows across it.C
C                                                                      C
C   Also in 2015 a major tidying up was performed to obtain the currentC
C   version. This included considerable rearranging of the input data  C
C   so that data from previous versions is no longer compatiable.      C
C                                                                      C
C======================================================================C
C                                                                      C
C   TO CHANGE THE DIMENSIONS OF THE ARRAYS PERFORM A GLOBAL CHANGE     C
C   OF THE  PARAMETER STATEMENTS IN "commall-open-18.3"                C
C                                                                      C
C   i.e of 'ID, JD, KD, MAXKI, NRS, IG1, JG1, KG1, IG2, JG2, KG2, JG3. C
C                                                                      C
C************************************************************************
C
C
C        THIS IS THE MAIN PROGRAM WHICH IS JUST USED TO CALL THE
C                            SUBROUTINES.
C   
C
      CHARACTER*1 ANSIN
C
      OPEN(UNIT=1,FILE='/dev/tty')
      OPEN(UNIT=4,FILE ='stage.log')
      OPEN(UNIT=7,FILE ='flow_out',   FORM = 'unformatted')
      OPEN(UNIT=11,FILE='global.plt', FORM = 'unformatted')
      OPEN(UNIT=3,FILE ='results.out')
      OPEN(UNIT=12,FILE='stopit')
C
      IFSTOP = 0
      WRITE(12,*) IFSTOP
      CLOSE(12)
C
C******************************************************************************
C     DECIDE WHICH INPUT STYLE TO USE
C
      OPEN(UNIT=13,FILE='intype')
            READ(13,*,END=10,ERR=10)   ANSIN
      GO TO 20
   10 WRITE(6,*) 'STOPPING BECAUSE FILE  "intype" DOES NOT EXIST '
      STOP 
C
   20 CONTINUE
C
      CLOSE(13)
C
      IF(ANSIN.EQ.'N'.OR.ANSIN.EQ.'n') THEN
             WRITE(6,*) ' NEW_READIN DATA FORMAT SPECIFIED.'
             CALL NEW_READIN
      ELSE
C
      IF(ANSIN.EQ.'O'.OR.ANSIN.EQ.'o')THEN
             WRITE(6,*)' OLD_READIN DATA FORMAT SPECIFIED.'
             CALL OLD_READIN
      ELSE
C
      WRITE(6,*) ' STOPPING BECAUSE FILE "intype" DOES NOT CONTAIN  "O" 
     & OR  "N" '
      STOP
C
      END IF
C
      END IF
C******************************************************************************
C 
      CALL SETUP(ANSIN)
C
C   LOOP calls  many other subroutines, especially TSTEP .
      CALL LOOP
C
C*********************************************************************************
C
      STOP
      END
C
C**************************************************************************************
C**************************************************************************************
C**************************************************************************************
C
      SUBROUTINE NEW_READIN
C
C       THIS SUBROUTINE READS IN THE DATA IN
C       ====================

C
      INCLUDE 'commall-open-18.3'
C
      COMMON/BKODDS/
     &           XINT(JD),YINT(JD),RINT(JD),SDIST(JD),
     &           ICUSP(NRS),LCUSP(NRS),LCUSPUP(NRS),
     &           FRACNEW(JD),BETANEW(JD),SLOPE(JD),
     &           THICKUP(JD),THICKLOW(JD),
     &           XINT1(KD),XINT2(KD),XINT3(KD),XINT4(KD)
C
      DIMENSION XHUB(JD),RHUB(JD),XTIP(JD),RTIP(JD),XQO(MAXKI),
     &          RQO(MAXKI),XNEWHUB(JD),RNEWHUB(JD),XNEWTIP(JD),
     &          RNEWTIP(JD)
C
C       START OF INPUT SECTION
C       THROUGHOUT THE INPUT SECTION THE VARIABLES ARE AS FOLLOWS
C       =====================
C
C       XSURF(J,K)    IS THE STORE FOR THE INPUT BLADE AXIAL COORDINATES
C       RT_UPP(J,K)   IS THE STORE FOR THE INPUT BLADE SUCTION SURFACE CO-ORD'S
C       RT_THICK(J,K) IS THE STORE FOR THE INPUT THE BLADE TANGENTIAL THICKNESS
C       RSURF(J,K)    IS THE STORE FOR THE INPUT BLADE RADIAL COORDINATES
C
C       =====================
 1200 FORMAT(A72)
 1700 FORMAT(8F12.6)
 1730 FORMAT(8F10.3)
 1800 FORMAT(8F10.1)
 1610 FORMAT(40I2)
 1000 FORMAT(10I5)
C
      PI     = 3.14159265
      DEGRAD = PI/180.
      RADDEG = 180./PI
C
C**************************************************************************************
C**************************************************************************************
C
C    			 START TO INPUT DATA.
C     FIRST FOR QUANTITIES WHICH ARE NOT DEPENDENT ON THE BLADE ROW
C
C**************************************************************************************
C**************************************************************************************
C
C     INPUT A TITLE FOR THE RUN. ANY CHARACTERS IN ROWS 1 to 72.
C
      WRITE(6,*)
      WRITE(6,*)
      READ(5,1200)  TITLE
      WRITE(6,1200) TITLE
      WRITE(6,*)
      WRITE(6,*)
C*******************************************************************************
C     INPUT THE GAS PROPERTIES, DEFAULTS TO CP =1005, GAMMA = 1.4
C
      IFGAS = 0
      CP = 1005.0
      GA = 1.4
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=7000) CP,GA
 7000 CONTINUE
      WRITE(6,*) ' GAS PROPERTIES:    CP = ', CP, ' GAMMA = ',GA
      WRITE(6,*)
C
C    USE REAL GAS PROPERTIES IF CP IS INPUT AS NEGATIVE.
C    TYPICAL VALUES FOR COMBUSTION PRODUCTS  ARE:CP1 = 1272.5, CP2 = 0.2125,
C    CP3 = 0.000015625, RGAS = 287.15 AT  TREF = 1400 K.
C
      IF(CP.LT.0.0) THEN
      CP1  = 1272.5
      CP2  = 0.2125
      CP3  = 0.000015625
      TREF = 1400.0
      RGAS = 287.15
      READ(5,*,ERR=7001) CP1, CP2, CP3, TREF, RGAS
 7001 CONTINUE
C
      WRITE(6,*) ' IDEAL GAS PROPERTIES READ IN '
      WRITE(6,*) ' CP1, CP2, CP3, TREF, RGAS =',CP1,CP2,CP3,TREF,RGAS
      WRITE(6,*)
C
      CPGAS  = CP1
      GAGAS  = CP1/(CP1 - RGAS)
      CP     = CP1
      GA     = GAGAS
      CV     = CP/GA
      IFGAS  = 1
      CALL SET_COEFFS
      END IF
C
C******************************************************************************
C    INPUT THE TIME STEPPING OPTION, 3 , 4 , -4,  5 or 6 . DEFAULT = 3 .
C
      ITIMST = 3
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR= 7002) ITIMST
 7002 CONTINUE
      WRITE(6,*) ' TIME STEP TYPE , ITIMST = ', ITIMST
      WRITE(6,*)
C
C    SET THE COEFFICIENTS FOR THE SSS SCHEME
C
      IF(ITIMST.EQ.3.OR.ITIMST.EQ.5.OR.ITIMST.EQ.6) THEN
              F1          =  2.0000
              F2          = -1.000
              F3          =  0.00
              F2EFF       = -1.0
              NRSMTH      =  0
              RSMTH       =  0.40
      END IF
C
      IF(ITIMST.EQ.-3.OR.ITIMST.EQ.-5.OR.ITIMST.EQ.-6) THEN
              F1          =  2.0000
              F2          = -1.65
              F3          = -0.65
              F2EFF       = -1.0
              NRSMTH      =  1
              RSMTH       =  0.40
              ITIMST      =  ABS(ITIMST)
      END IF
C
      IF(ITIMST.EQ.4.OR.ITIMST.EQ.-4) THEN
              READ(5,*) F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*) ' F1, F2EFF, F3 , RSMTH, NRSMTH ',
     &                     F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*)
              IF(F2EFF.GT.0.0) THEN 
              WRITE(6,*)
              WRITE(6,*) ' ERROR,  F2EFF  MUST BE NEGATIVE.'
              WRITE(6,*) ' THE SIGN OF THE INPUT VALUE WILL BE CHANGED.'
              WRITE(6,*)
              F2EFF = -F2EFF
              END IF
              IF(F3.GT.0.0) THEN 
              WRITE(6,*)
              WRITE(6,*) ' ERROR,  F3  MUST BE NEGATIVE.'
              WRITE(6,*) ' THE SIGN OF THE INPUT VALUE WILL BE CHANGED.'
              WRITE(6,*)
              F3  = -F3
              END IF
              F2     = F2EFF*(1.0 - F3)
              ITIMST = 3
      END IF
C
              WRITE(6,*) ' F1, F2EFF, F3 , RSMTH, NRSMTH ',
     &                     F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*)
C
C**********************************************************************************
C
C    INPUT THE ARTIFICIAL SPEED OF SOUND IF ITIMST = 5.
C    THIS SHOULD BE ABOUT HALF THE MAXIMUM RELATIVE VELOCITY IN THE FLOW.
C
      IF(ITIMST.GE.5) THEN
      VSOUND    = 150.
      RF_PTRU   = 0.01
      RF_VSOUND = 0.002
      DENSTY    = 1.20
      VS_VMAX   = 2.0
      IF(ITIMST.EQ.5)  READ(5,*,END= 2350)
     &  VSOUND, RF_PTRU, RF_VSOUND, VS_VMAX
      IF(ITIMST.EQ.6)  READ(5,*,END= 2350)
     &  VSOUND, RF_PTRU, RF_VSOUND, VS_VMAX, DENSTY
 2350 CONTINUE
         RF_PTRU1    = 1.0 - RF_PTRU
         RF_VSOUND1  = 1.0 - RF_VSOUND
         WRITE(6,*)
         WRITE(6,*) ' CALCULATION USING ARTIFICIAL COMPRESSIBILITY '
         WRITE(6,*) ' ARTIFICIAL SPEED OF SOUND = ', VSOUND
         WRITE(6,*) ' DENSITY RELAXATION FACTOR = ', RF_PTRU
         WRITE(6,*) ' SOUND SPEED RELAXATION FACTOR = ', RF_VSOUND
         WRITE(6,*) ' RATIO OF SOUND SPEED TO MAXIMUM SPEED = ',VS_VMAX
         IF(ITIMST.EQ.6) WRITE(6,*) 
     &   ' INCOMPRESSIBLE FLOW WITH DENSITY = ',DENSTY
         WRITE(6,*)
      END IF
C
C******************************************************************************
C   INPUT THE CFL NUMBER, DEFAULT VALUE = 0.4, BUT INCREASE TO 0.7 IF USING THE SSS SCHEME.
C   ALSO THE DAMPING FACTOR, MACH NUMBER LIMITER AND FRACTION OF PRESSURE DOWNWINDING.
C
      CFL     = 0.4
      DAMPIN  = 10.0
      MACHLIM = 2.0
      F_PDOWN = 0.0
      READ(5,*) DUMMY_INPUT
      READ(5,*,END=7003) CFL, DAMPIN, MACHLIM, F_PDOWN
 7003 CONTINUE
      WRITE(6,*) ' CFL NUMBER = ', CFL,' DAMPING FACTOR = ',DAMPIN,
     & ' MACHLIM = ',MACHLIM, ' F_PDOWN = ', F_PDOWN
      WRITE(6,*)
C
C******************************************************************************
C    READ IN THE RESTART AND OUTPUT FILE OPTIONS.
C    SET  "IF_RESTART"  = 1 TO START FROM A RESTART FILE.  
C    THE COMBINED RESTART AND PLOTTING FILE,  "flow_out"   IS ALWAYS WRITTEN.
C
      IF_RESTART  = 0
      READ(5,*)            DUMMY_INPUT
      READ(5,*,ERR = 7017) IF_RESTART
 7017 CONTINUE
      WRITE(6,*) 'RESTART FILE OPTIONS, IF_RESTART= ',IF_RESTART
      WRITE(6,*)
C
C**************************************************************************
C   INPUT THE MAXIMUM NUMBER OF TIME STEPS AND CONVERGENCE LIMIT.
C
      NSTEPS_MAX = 5000
      CONLIM     = 0.005
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR = 7004) NSTEPS_MAX, CONLIM
 7004 CONTINUE
      WRITE(6,*) ' MAXIMUM TIME STEPS= ',NSTEPS_MAX,
     &           ' CONVERGENCE LIMIT = ',CONLIM
      WRITE(6,*)
C
C****************************************************************************
C      INPUT THE SMOOTHING FACTORS, FRACTION OF FOURTH ORDER SMOOTHING AND STEPS OVER
C      WHICH THEY ARE GRADUALLY DECREASED
C      DEFAULT VALUES ARE:   0.005,  0.005 , 0.8 AND NSTEPS/4  .
C
      SFXIN   = 0.005
      SFTIN   = 0.005
      FAC_4TH = 0.8
      NCHANGE = NSTEPS_MAX/4
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=7005) SFXIN,SFTIN,FAC_4TH,NCHANGE
 7005 CONTINUE
C    8/4/2017.  SET NCHANGE = 100 IF STARTING FROM A RESTART FILE.
      IF(IF_RESTART.NE.0) NCHANGE = 100
      WRITE(6,*) ' SFXIN,SFTIN,FAC_4TH,NCHANGE ',
     &             SFXIN,SFTIN,FAC_4TH,NCHANGE
      WRITE(6,*)
C
C******************************************************************************
C******************************************************************************
C   INPUT THE NUMBER OF BLADE ROWS TO BE CALCULATED
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    NROWS
      WRITE(6,*) ' NUMBER OF BLADE ROWS TO BE CALCULATED = ',NROWS
      WRITE(6,*)
C
C**************************************************************************
C   INPUT THE NUMBER OF GRID POINTS IN THE PITCHWISE (I) AND SPANWISE (K) DIRECTIONS.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    IM,KM
      WRITE(6,*) ' NUMBER OF PITCHWISE GRID POINTS = ',IM
      WRITE(6,*) ' NUMBER OF SPANWISE GRID POINTS  = ',KM
      IMM1 = IM-1
      IMM2 = IM-2
      IF(IM.EQ.2) IMM2=1
      KMM1 = KM-1
      KMM2 = KM-2
      WRITE(6,*)
C
C
C   INPUT THE RELATIVE SPACING OF THE GRID POINTS IN THE PITCHWISE DIRECTION. 
      READ(5,*)     DUMMY_INPUT
      WRITE(6,*)    DUMMY_INPUT
      READ(5,*)    (FP(I),I=1,IMM1)
      WRITE(6,*) ' RELATIVE SPACING OF THE GRID POINTS IN THE PITCHWISE
     & DIRECTION .'
      WRITE(6,*) ' THIS IS THE SAME FOR ALL BLADE ROWS.'
      WRITE(6,1700) (FP(I),I=1,IMM1)
      WRITE(6,*)
C
C   INPUT THE RELATIVE SPACING OF THE GRID POINTS IN THE SPANWISE DIRECTION.'
      READ(5,*)      DUMMY_INPUT
      WRITE(6,*)     DUMMY_INPUT
      READ(5,*)     (FR(K),K=1,KMM1)
      WRITE(6,*) ' RELATIVE SPACING OF THE GRID POINTS IN THE SPANWISE
     & DIRECTION .'
      WRITE(6,*) ' THIS IS THE SAME FOR ALL BLADE ROWS.'
      WRITE(6,1700) (FR(K),K=1,KMM1)
      WRITE(6,*)
C**************************************************************************      
C    INPUT THE MULTIGRID BLOCK SIZES. DEFAULTS = 3  and 9
C
      IR   = 3
      JR   = 3
      KR   = 3
      IRBB = 9
      JRBB = 9
      KRBB = 9
      READ(5,*) DUMMY_INPUT
      READ (5,*,ERR=7007) IR,JR,KR,IRBB,JRBB,KRBB
 7007 CONTINUE
      IF(KM.EQ.2) THEN
           KR   = 1
           KRBB = 1
      END IF
      IF(IM.EQ.2) THEN
           IR   = 1
           IRBB = 1
      END IF
C
      WRITE(6,*) ' MULTIGRID BLOCK SIZES = ', IR,JR,KR,IRBB,JRBB,KRBB
      WRITE(6,*)
C
C************************************************************************
C   INPUT THE MULTIGRID TIME STEP FACTORS, FBLK1,  FBLK2, FBLK3
C
      FBLK1 = 0.4
      FBLK2 = 0.2
      FBLK3 = 0.1
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=8007) FBLK1,FBLK2,FBLK3
      WRITE(6,*) ' MULTIGRID TIME STEP FACTORS= ',FBLK1,FBLK2,FBLK3
 8007 CONTINUE
      WRITE(6,*)
C
C******************************************************************************
C    READ IN THE MIXING PLANE PARAMETERS
C
      IFMIX = 1
      READ(5,*)    DUMMY_INPUT
      READ(5,*,ERR=7009) IFMIX
 7009 CONTINUE
      WRITE(6,*) ' MIXING PLANE MARKER, IFMIX = ', IFMIX
      WRITE(6,*)
C
      IF(IFMIX.NE.0) THEN
      RFMIX    = 0.025
      FSMTHB   = 1.0
      FEXTRAP  = 0.80
      FANGLE   = 0.80
      READ(5,*) DUMMY_INPUT
      READ(5,*,END=7010,ERR=7010)  RFMIX,FEXTRAP,FSMTHB,
     &                             FANGLE
 7010 CONTINUE
      WRITE(6,*)' THE MIXING PLANE PARAMETERS ARE '
      WRITE(6,*)' RFMIX, FEXTRAP, FSMTHB ,FANGLE',
     &            RFMIX, FEXTRAP, FSMTHB, FANGLE
      WRITE(6,*)
      END IF
C
C******************************************************************************
C   READ IN MARKERS FOR BLEED FLOWS, COOLING FLOWS AND SURFACE ROUGHNESS
C
      IFCOOL   = 0
      IFBLEED  = 0
      IF_ROUGH = 0
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR = 7011) IFCOOL,IFBLEED,IF_ROUGH
 7011 CONTINUE
      WRITE(6,*) ' IFCOOL,IFBLEED,IF_ROUGH ', IFCOOL,IFBLEED,IF_ROUGH
      WRITE(6,*)
C
C******************************************************************************
C   READ IN THE NUMBER OF BLADE SECTIONS ON WHICH THE GEOMETRY IS TO BE GENERATED.
C 
      READ(5,*)   DUMMY_INPUT 
      WRITE(6,*)  DUMMY_INPUT 
      READ(5,*)   NSECS_IN
      WRITE(6,*)
     &'NUMBER OF SECTIONS USED FOR BLADE GEOMETRY GENERATION= ',NSECS_IN
      WRITE(6,*)
C
C******************************************************************************
C      READ IN THE INLET BOUNDARY CONDITION OPTIONS
C
      IN_PRESS   = 0
      IN_VTAN    = 0
      IN_VR      = 1
      IN_FLOW     = 0
      IF_REPEAT  = 0
      RFIN       = 0.1
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=7012) IN_PRESS,IN_VTAN,IN_VR,IN_FLOW,IF_REPEAT,RFIN
 7012 CONTINUE
      WRITE(6,*) ' IN_PRESS,IN_VTAN,IN_VR,IN_FLOW,IF_REPEAT ',
     &             IN_PRESS,IN_VTAN,IN_VR,IN_FLOW,IF_REPEAT
      WRITE(6,*) ' RFIN = ', RFIN
      WRITE(6,*)
C
C******************************************************************************
C
C      READ IN THE EXIT BOUNDARY CONDITION OPTIONS
      IPOUT          = 1
      SFEXIT         = 0.0
      NSFEXIT        = 0
      READ(5,*)          DUMMY_INPUT
      READ(5,*,ERR=8012) IPOUT,SFEXIT,NSFEXIT
 8012 CONTINUE
      WRITE(6,*) ' IPOUT, SFEXIT ', IPOUT, SFEXIT, NSFEXIT
      WRITE(6,*) ' SFEXIT= ',SFEXIT, ' NSFEXIT= ', NSFEXIT
C
C****************************************************************************** 
C
      PLATE_LOSS    = 0.0
      THROTTLE_EXIT = 0.0
      READ(5,*)   DUMMY_INPUT
      READ(5,*,ERR=8013)   PLATE_LOSS, THROTTLE_EXIT
 8013 CONTINUE
      WRITE(6,*)' LOSS COEFFICIENT OF A PERFORATED PLATE AT EXIT = ',
     &            PLATE_LOSS
      WRITE(6,*)' MARKER FOR SETTING A THROTTLE AT THE EXIT BOUNDARY= ',
     &            THROTTLE_EXIT 
C
      IF(THROTTLE_EXIT.GT.0.001) THEN
            READ(5,*)  THROTTLE_PRES,THROTTLE_MAS,RFTHROTL
            WRITE(6,*)'A THROTTLE EXIT BOUNDARY CONDITION HAS BEEN SET.'
            WRITE(6,*)'THROTTLE_PRES, THROTTLE_MAS, RFTHROTL=',
     &                 THROTTLE_PRES,THROTTLE_MAS,RFTHROTL 
           RFTHROTL1 = 1.0 - RFTHROTL
      END IF
      WRITE(6,*)
C
C******************************************************************************
C     READ IN IN THE SPECIFIED INLET FLOW AND RELAXATION FACTOR IF IN_FLOW
C     IS NOT = ZERO
C
      IF(IN_FLOW.NE.0) THEN
           READ(5,*)    DUMMY_INPUT
           READ(5,*)    FLOWIN, RFLOW
           WRITE(6,*) ' INLET FLOW FORCING, FLOWIN, RFLOW ',FLOWIN,RFLOW
           WRITE(6,*)
      END IF
C
C******************************************************************************
C     READ IN THE FACTORS USED FOR REPEATING FLOW CONDITIONS
C
      IF(IF_REPEAT.NE.0) THEN
           NINMOD = 10
           RFINBC = 0.025
           READ(5,*) DUMMY_INPUT
           READ(5,*,ERR=7013)  NINMOD, RFINBC
 7013      CONTINUE
           WRITE(6,*) ' REPEATING STAGE SPECIFIED, NINMOD = ', NINMOD,
     &                ' RFINBC = ', RFINBC
           WRITE(6,*)
      END IF
C
C******************************************************************************
C    READ IN THE CHOICE OF VISCOUS MODEL
C
      ILOS   = 100
      NLOS   = 5
      IBOUND = 0
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=7014) ILOS, NLOS, IBOUND
 7014 CONTINUE
      IF(KM.EQ.2) IBOUND = 2
      WRITE(6,*) ' VISCOUS MODEL SET BY ILOS= ',ILOS,' NLOS= ',NLOS,
     &           ' IBOUND= ', IBOUND
      WRITE(6,*)
C
C******************************************************************************
C    READ IN THE PARAMETERS FOR THE VISCOUS MODEL.

      REYNO       = 500000.0
      RF_VIS      = 0.5
      FTRANS      = 0.0001
      PRANDTL     = 1.0
      YPLUSWALL   = 0.0
      TURBVIS_LIM = 3000.0
C
      IF(ILOS.NE.0) THEN
      READ(5,*)          DUMMY_INPUT
      READ(5,*,ERR=7015) REYNO,RF_VIS,FTRANS,TURBVIS_LIM,
     &                   PRANDTL,YPLUSWALL
 7015 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*) ' REYNOLDS NUMBER                       = ',REYNO
      WRITE(6,*) ' VISCOUS TERM RELAXATION FACTOR        = ',RF_VIS
      WRITE(6,*) ' TRANSITION FACTOR, FTRANS             = ',FTRANS
      WRITE(6,*) ' LIMIT ON TURBULENT/LAMINAR VISCOSITY  = ',TURBVIS_LIM
      WRITE(6,*) ' PRANDTL NUMBER                        = ',PRANDTL
      WRITE(6,*) ' YPLUSWALL - IF USED (NOT OFTEN USED)  = ',YPLUSWALL
      WRITE(6,*)
C
      END IF
C
C   IF USING THE SA TURBULENCE MODEL
C
      IF(ILOS.GE.200)THEN
           FAC_STMIX = 0.0
           FAC_ST0   = 1.0
           FAC_ST1   = 1.0
           FAC_ST2   = 1.0
           FAC_ST3   = 1.0
           FAC_SFVIS = 2.0
           FAC_VORT  = 0.0
           FAC_PGRAD = 0.0
           READ(5,*) DUMMY_INPUT
           READ(5,*,ERR= 7016) FAC_STMIX, FAC_ST0, FAC_ST1,
     &                         FAC_ST2  , FAC_ST3, FAC_SFVIS,
     &                         FAC_VORT, FAC_PGRAD
 7016 CONTINUE
C
           WRITE(6,*)
           WRITE(6,*) ' SPALART-ALLMARAS TURBULENCE MODEL IS BEING USED'
           WRITE(6,*) ' THE S_A SOURCE TERM MULTIPLIERS ARE ',
     &                  FAC_STMIX, FAC_ST0, FAC_ST1, FAC_ST2, FAC_ST3,
     &                  FAC_VORT, FAC_PGRAD 
           WRITE(6,*) ' THE TURBULENT VISCOSITY SMOOTHING FACTOR IS ', 
     &                  FAC_SFVIS
           WRITE(6,*)
      END IF
C
C**********************************************************************************
C    INPUT THE RANGE OF YPLUS VALUES OVER WHICH THE TURBULENT VISCOSITY WILL BE REDUCED.
C
      IF(ILOS.NE.0) THEN
           YPLAM    = 5.0
           YPTURB   = 25.0
           READ(5,*) DUMMY_INPUT
           READ(5,*, ERR= 590) YPLAM,YPTURB
  590 CONTINUE
           WRITE(6,*)' THE TURBULENT VISCOSITY IS REDUCED OVER THE RANGE 
     &  YPLUS = ',YPLAM,' TO ',YPTURB 
      END IF
C
C******************************************************************************
C   INPUT THE FORCING FACTOR AND SMOOTHING FACTOR IF DOING A THROUGHFLOW  CALCULATION
C
      IF(IM.EQ.2) THEN
      Q3DFORCE = 1.0
      SFPBLD   = 0.1
      NSFPBLD  = 2
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=591) Q3DFORCE, SFPBLD, NSFPBLD
  591 CONTINUE
      SFPBLD1 = 1.0 - SFPBLD
      WRITE(6,*)
      WRITE(6,*) 'THROUGHFLOW CALCULATION REQUESTED, Q3DFORCE= ',
     & Q3DFORCE, 'SMOOTHING FACTOR=',SFPBLD,' No OF SMOOTHINGS=',NSFPBLD
      WRITE(6,*)
      END IF       
C
C******************************************************************************
C    READ IN THE BLADE ROW AND GRID ALIGNMENT OPTIONS
C
      ISHIFT      = 2
      NEXTRAP_LE  = 10
      NEXTRAP_TE  = 10
      READ(5,*)          DUMMY_INPUT
      READ(5,*,ERR=7008) ISHIFT,NEXTRAP_LE,NEXTRAP_TE
 7008 CONTINUE
      WRITE(6,*) ' ISHIFT ,NEXTRAP_LE, NEXTRAP_TE = ',
     &             ISHIFT, NEXTRAP_LE, NEXTRAP_TE
      WRITE(6,*)
C
C**************************************************************************************
C
C     READ IN THE STAGE NUMBERS AND SORT OUT THE START AND END OF EACH STAGE.
C
C     FIRST SET DEFAULTS.
      DO  N = 1,NROWS
      NSTAGE(N)    = 1 + (N-1)/2
      END DO
C
C     NOW READ IN THE ACTUAL STAGE NUMBER FOR EACH BLADE ROW.
C
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=588)(NSTAGE(N),N=1,NROWS)
  588 CONTINUE
C
C**************************************************************************
C**************************************************************************
C   INPUT 5 TIME STEPS WHERE OUTPUT IS REQUESTED.
C   SET TO A NUMBER ABOVE NSTEPS_MAX IF NONE REQUIRED.
C
      DO N=1,5
           NOUT(N) = NSTEPS_MAX+10 
      END DO
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=7006) (NOUT(N),N=1,5)
 7006 CONTINUE
      WRITE(6,*) ' OUTPUT REQUESTED AT TIME STEPS = ', (NOUT(N),N=1,5)
      WRITE(6,*)
C******************************************************************************
C     INPUT A LIST OF VARIABLES TO BE SENT TO THE OUTPUT FILE ' RESULTS.OUT'.
      DO I=1,20
           IOUT(I) = 0
      END DO
      READ(5,*)      DUMMY_INPUT
      READ(5,*,ERR=7027) (IOUT(I),I=1,13)
 7027 CONTINUE
      WRITE(6,*)'THE VARIABLES TO BE OUTPUT TO THE FILE RESULTS.OUT ARE'
      WRITE(6,1610) (IOUT(I),I=1,13)
C***************************************************************************************
C     CHOOSE WHICH K VALUES (STREAM SURFACES) ARE TO BE OUTPUT TO THE FILE 'RESULTS.OUT'.
      DO K=1,KM
           KOUT(K) = 0
      END DO
      READ(5,*)      DUMMY_INPUT
      READ(5,*,ERR=2028) (KOUT(K),K=1,KM)
 2028 CONTINUE
      WRITE(6,*)  'THE STREAM SURFACES ON WHICH RESULTS ARE SENT TO THE
     &OUTPUT FILE ARE:'
      WRITE(6,1610) (KOUT(K),K=1,KM)
C
C******************************************************************************
C**************************************************************************************
C    END OF DATA INPUT WHICH IS NOT DEPENDENT ON THE BLADE ROW.
C**************************************************************************************
C**************************************************************************************
C
C     CHECK THAT THE DIMENSIONS ARE NOT TOO LARGE.
C
      IF(MAXKI.LT.KD)  WRITE(6,*)
     &            ' STOPPING BECAUSE MAXKI IS LESS THAN KD',
     &            ' MAXKI = ',MAXKI, ' KD = ',KD
      IF(MAXKI.LT.ID)  WRITE(6,*)
     &            ' STOPPING BECAUSE MAXKI IS LESS THAN ID',
     &            ' MAXKI = ',MAXKI, ' ID = ',ID
      IF(IM.GT.ID)  WRITE(6,*) ' STOPPING BECAUSE IM TOO LARGE.',
     &            ' IM= ',IM,  ' DIMENSION LIMIT = ',ID
      IF(KM.GT.KD)  WRITE(6,*) ' STOPPING BECAUSE KM TOO LARGE.',
     &            ' KM= ',KM,  ' DIMENSION LIMIT = ',KD
C
      IF(MAXKI.LT.KD.OR.MAXKI.LT.ID.OR.IM.GT.ID.OR.KM.GT.KD) STOP
C
      IF(NROWS.GT.NRS)  WRITE(6,*) ' STOPPING BECAUSE NROWS TOO LARGE.',
     &            ' NROWS= ',NROWS,' DIMENSION LIMIT = ',NRS
      IF(NROWS.GT.NRS) STOP
C
C**************************************************************************************
C**************************************************************************************
C**************************************************************************************
C**************************************************************************************
C
C            NOW INPUT THE CONTROL PARAMETERS FOR EACH BLADE ROW
C            AND SET THE VALUES OF SOME BLADE ROW VARIABLES
C
      J1        = 1
      IFSHROUD  = 0
C
      WRITE(6,*)'******************************************************'
      WRITE(6,*)'******************************************************'
C
      WRITE(6,*)
      WRITE(6,*) ' STARTING INPUT FOR EACH BLADE ROW '
      WRITE(6,*)
C
      DO 1550 NR = 1,NROWS
C
C    READ IN TWO BLANK LINES TO HELP SPACE OUT THE DATA.
      READ(5,*) DUMMY_INPUT
      READ(5,*) DUMMY_INPUT
      WRITE(6,*)'******************************************************'
      WRITE(6,*)'******************************************************'
C
C******************************************************************************
C     READ IN THE ROW TITLE. THIS IS NEVER USED BUT HELPS TO LAY OUT THE DATA
      READ(5,1200) ROWTYP
C
      WRITE(6,*) ' INPUTTING DATA FOR BLADE ROW NUMBER ',NR
      WRITE(6,*) ' BLADE ROW TITLE = ', ROWTYP
      WRITE(6,*)'******************************************************'
      WRITE(6,*)'******************************************************'
C
C******************************************************************************
C     READ IN THE NUMBER OF BLADES IN THIS ROW  
      READ(5,*)    DUMMY_INPUT    
      READ(5,*)    NBLADES_IN_ROW
      WRITE(6,*)
      WRITE(6,*) ' NUMBER OF BLADES IN ROW No.',NR,' =',NBLADES_IN_ROW
      WRITE(6,*)
C
C******************************************************************************
C    JMROW = NUMBER OF J GRID POINTS ON THIS ROW.
C    JLEROW AND JTEROW ARE MEASURED FROM THE START OF THIS ROW. 
C    THEY ARE NOT THE FINAL OVERALL VALUES.
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    JMROW,JLEROW,JTEROW
      WRITE(6,*) ' JM = ',JMROW,'JLE =',JLEROW,'JTE = ',JTEROW
      WRITE(6,*)
C     SET THE OVERALL J VALUES OF THE LEADING AND TRAILING EDGE POINTS ON THIS ROW.   
      JLEE = J1+JLEROW-1
      JTEE = J1+JTEROW-1
C
C******************************************************************************
C
C      KTIPS IS THE K VALUE OF THE POINT WHERE THE TIP GAP STARTS
C      KTIPE IS THE K VALUE OF THE POINT WHERE THE TIP GAP ENDS.
C      SET KTIPS(NR) = 0  FOR NO TIP GAP ON THIS ROW.    
C
C      SET KTIPS(NR) NEGATIVE TO USE THE SHROUD LEAKAGE MODEL ON THIS
C      ROW. EXTRA INPUT DATA IS THEN NEEDED AT THE END OF THE DATA FILE.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    KTIPS(NR),KTIPE(NR)
      WRITE(6,*) ' K TIP START= ',KTIPS(NR),' K TIP END= ',KTIPE(NR)
      WRITE(6,*)
C     FOR Q3D
      IF(KM.EQ.2) THEN
           KTIPS(NR) = 0
           KTIPE(NR) = 0
      END IF
C     END Q3D
C     SET MARKER  "IFSHROUD" FOR A TIP SHROUUD IF KTIPS IS NEGATIVE ON ANY BLADE ROW.
C
      IF(KTIPS(NR).LE.0) KTIPE(NR) = 0
C     SET   "IFSHROUD" = 1   IF ANY BLADE ROW HAs A SHROUD.
C     IF NO BLADE ROW HAS A SHROUD   "IFSHROUD"   REMAINS AT 0 .
      IF(KTIPS(NR).LT.0) IFSHROUD  = 1
C
C*******************************************************************************
C      INPUT THE TIP CLEARANCE AS A FRACTION OF THE BLADE SPAN AT THE LEADING EDGE
C      AND TRAILING EDGE.
       IF(KTIPS(NR).GT.0) THEN
            READ(5,*)   DUMMY_INPUT
            READ(5,*)   FRACTIP1(NR),FRACTIP2(NR)
            WRITE(6,*)' FRACTIP1, FRACTIP2 = ',FRACTIP1(NR),FRACTIP2(NR)
            WRITE(6,*)
      END IF
C
C*******************************************************************************
C      FTHICK(NR,K) IS THE MULTIPLYING FACTOR ON THE BLADE THICKNESS SO
C      THAT IT CAN BE REDUCED AT AND BEYOND THE TIP.
C
      DO 2507 K=1,KM
 2507 FTHICK(NR,K)=1.0
C
C      FTHICK  IS ASSUMED TO BE 1.0 UNLESS INPUT HERE.
C      SET FTHICK = 0.0 FOR POINTS IN THE TIP GAP, INCLUDING THE POINT AT THE BLADE TIP.
      IF(KTIPS(NR).GT.0.AND.KTIPE(NR).GT.0) THEN
           WRITE(6,*)   ' BLADE TIP THINNING FACTORS '
           READ(5,*)      DUMMY_INPUT
           READ(5,*)     (FTHICK(NR,K),K=1,KM)
           WRITE(6,1700) (FTHICK(NR,K),K=1,KM)
           WRITE(6,*) 
      END IF    
C
C******************************************************************************
C   BOUNDARY LAYER TRANSITION IS SPECIFIED SEPARATELY ON EACH SURFACE.
C   BUT IS ALSO DETERMINED BY FTRANS WHICH HAS ALREADY BEEN INPUT . 
C   JTRAN_I1, etc,  ARE MEASRED FROM THE START OF THIS ROW.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    JTRAN_I1(NR),JTRAN_IM(NR),JTRAN_K1(NR),JTRAN_KM(NR)
      WRITE(6,*) ' BOUNDARY LAYER TRANSITION IS SPECIFIED AS FOLLOWS '
      WRITE(6,*) ' TRANSITION ON I= 1 SURFACE FIXED AT J= ',JTRAN_I1(NR)
      WRITE(6,*) ' TRANSITION ON I=IM SURFACE FIXED AT J= ',JTRAN_IM(NR)
      WRITE(6,*) ' TRANSITION ON K= 1 SURFACE FIXED AT J= ',JTRAN_K1(NR)
      WRITE(6,*) ' TRANSITION ON K=KM SURFACE FIXED AT J= ',JTRAN_KM(NR)
      WRITE(6,*)
C
C******************************************************************************
C     SET NEW_GRID = 1  TO GENERATE A NEW GRID IN THE "J" (STREAMWISE) DIRECTION
C     DATA FOR THIS IS READ IN LATER. SET NEW_GRID = 0 TO USE THE GRID INPUT HERE.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    NEW_GRID
      WRITE(6,*) ' NEWGRID = ', NEW_GRID
      IF(NEW_GRID.NE.0) WRITE(6,*)
     & ' A NEW GRID WILL BE GENERATED USING DATA TO BE INPUT LATER. '
      WRITE(6,*)
C
C
C******************************************************************************
C     READ IN THE RPM AND HUB ROTATION OF THIS ROW. 
C     THE CASING ROTATION IS TAKEN TO BE RPMROW BETWEEN JROTTS AND JROTTE AND ZERO ELSEWHERE.
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    RPMROW, RPMHUB  
      WRITE(6,*) ' BLADE ROW RPM = ',RPMROW,' HUB RPM = ',RPMHUB
      WRITE(6,*)
C
C******************************************************************************
C     THE HUB IS TAKEN TO BE ROTATING AT RPMHUB BETWEEN JROTHS AND JROTHE .
C     THE CASING IS TAKEN TO BE ROTATING AT RPMROW BETWEEN JROTTS AND JROTTE.
C     OUTSIDE THESE LIMITS THE HUB AND CASING ARE NOT ROTATING.
C     JROTHS, JROTHE, etc  ARE MEASURED FROM THE J VALUE AT START OF THIS BLADE ROW
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    JROTHS,JROTHE,JROTTS,JROTTE
C
C   SET ALL ROTATION POINTS TO JMROW IF THEY ARE OUTSIDE THE LIMITS FOR THIS ROW.
      IF(JROTHS.EQ.0.OR.JROTHS.GT.JMROW) JROTHS = JMROW
      IF(JROTHE.EQ.0.OR.JROTHE.GT.JMROW) JROTHE = JMROW
      IF(JROTTS.EQ.0.OR.JROTTS.GT.JMROW) JROTTS = JMROW
      IF(JROTTE.EQ.0.OR.JROTTE.GT.JMROW) JROTTE = JMROW
C
      WRITE(6,*) ' HUB ROTATION STARTS AT J = ',JROTHS
      WRITE(6,*) ' HUB ROTATION ENDS   AT J = ',JROTHE
      WRITE(6,*) ' TIP ROTATION STARTS AT J = ',JROTTS
      WRITE(6,*) ' TIP ROTATION ENDS   AT J = ',JROTTE
      WRITE(6,*)
C
C****************************************************************************** 
C     INPUT AN INITIAL GUESS OF THE UPSTREAM, LEADING EDGE, TRAILING EDGE AND 
C     DOWNSTREAM STATIC PRESSURES FOR THIS ROW. 
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    PUPROW,PLEROW,PTEROW,PDNROW
      WRITE(6,*) ' INITIAL GUESS OF PRESSURES = ',
     &             PUPROW,PLEROW,PTEROW,PDNROW 
      WRITE(6,*)
C
C******************************************************************************
C     INPUT THE NUMBER OF INPUT SECTIONS FOR THIS BLADE ROW.
      READ(5,*)   DUMMY_INPUT
      NSECS_ROW  = NSECS_IN
      INSURF     = 0
      READ(5,*,ERR=7019) NSECS_ROW,INSURF
 7019 CONTINUE
      WRITE(6,*)'NUMBER OF INPUT SECTIONS FOR THIS BLADE ROW=',NSECS_ROW
      WRITE(6,*)'     NEW ENDWALL GENERATION OPTION, INSURF =', INSURF
      IF(NSECS_IN.EQ.1)  NSECS_ROW = 1
      IF(NSECS_IN.EQ.1)  NSECS_IN  = 2
      IF(NSECS_ROW.EQ.1) NSECS_ROW = 2
C
C******************************************************************************
C******************************************************************************
C     Set the cusp generation parameters if IF_CUSP is not 0.
C     Maintain the original grid and no cusp is generated if IF_CUSP = 0
C     The cusp is centred on the blade centre line if ICUSP = 0.
C     The cusp makes the I=1 surface continuous on the cusp if  ICUSP =  1.
C     The cusp makes the I=IM surface continuous on the cusp if ICUSP = -1.
C     The cusp is of length  LCUSP and starts LCUSPUP points before the
C     trailing edge.
C
      IF_CUSP(NR)    = 0
      IF_ANGLES(NR)  = 0
      READ(5,*)  DUMMY_INPUT
      READ(5,*,ERR=7020) IF_CUSP(NR), IF_ANGLES(NR)
 7020 CONTINUE
      WRITE(6,*)'           CUSP GENERATION OPTION, IF_CUSP =',
     &           IF_CUSP(NR)
      WRITE(6,*)'UPSTREAM AND DOWNSTREAM ANGLE INPUT OPTION =',
     &           IF_ANGLES(NR)
      WRITE(6,*)
C
      IF(IF_CUSP(NR).EQ.1) THEN
	   READ(5,*)    ICUSP(NR),LCUSP(NR),LCUSPUP(NR)
           WRITE(6,*) ' ICUSP, LCUSP , LCUSPUP = ',
     &                  ICUSP(NR),LCUSP(NR),LCUSPUP(NR)
           WRITE(6,*)
      END IF
C
C    IF  IFCUSP = 2 a body force is used to force separation at the trailing edge.
C    The force starts NUP_I1 grid points upstream of the TE  on the I=1 blade surface, 
C    and at NUP_IM points upstream on the I=IM blade surface. It extends N_WAKE points
C    downstream of the TE. The thickness of the affected layer is determined by SEP_THIK,
C    typical value  0.01, and the strength of the body force by SEP_DRAG, typical value 0.99.
C
      IF(IF_CUSP(NR).EQ.2) THEN
           ICUSP(NR)   = 0
           LCUSP(NR)   = 0
           LCUSPUP(NR) = 0
           READ(5,*)  NUP_I1(NR),NUP_IM(NR),N_WAKE(NR),
     &                  SEP_THIK(NR),SEP_DRAG(NR)
           WRITE(6,*)' BODY FORCE USED TO FORCE TRAILING EDGE SEPARTION'
           WRITE(6,*)' CUSP BODY FORCE POINTS, NUP_I1, NUP_IM, N_WAKE, 
     & SEP_THICK, SEP_DRAG'
             WRITE(6,*) NUP_I1(NR),NUP_IM(NR),N_WAKE(NR),
     &                  SEP_THIK(NR),SEP_DRAG(NR)
             WRITE(6,*)
      END IF
C******************************************************************************
C******************************************************************************
C
C   SET THE J VALUE RELATIVE TO THE OVERALL FIRST GRID POINT
      J2 = J1 + JMROW - 1
C
C*******************************************************************************
C*******************************************************************************
C*******************************************************************************
C
C      NOW  READ IN MAIN GEOMETRICAL DATA FOR THE CURRENT BLADE ROW 
C      ON  "NSECS_ROW"  BLADE SECTIONS.
C
C*******************************************************************************
C*******************************************************************************
C*******************************************************************************
C
      WRITE(6,*)'*******************************************************
     &********************'
      WRITE(6,*)
     &      ' STARTING TO INPUT THE BLADE COORDINATES FOR ROW NUMBER',NR
      WRITE(6,*)'*******************************************************
     &********************'
      WRITE(6,*)
C
      DO 1600 K=1,NSECS_ROW
C
      WRITE(6,*)'*******************************************************
     &********************'
      WRITE(6,*)'*******************************************************
     &********************'
C
      IF_DESIGN  = 0
      IF_STAGGER = 0
      IF_LEAN    = 0
      READ(5,*)  DUMMY_INPUT
      READ(5,*)  DUMMY_INPUT 
      READ(5,*,  ERR=1601) IF_DESIGN,IF_STAGGER,IF_LEAN
 1601 CONTINUE 
C
      WRITE(6,*) ' IF_DESIGN = ', IF_DESIGN
      IF(IF_DESIGN.NE.0) WRITE(6,*)
     & ' THIS BLADE SECTION IS TO BE DESIGNED'
C
      WRITE(6,*) ' IF_STAGGER = ', IF_STAGGER
      IF(IF_STAGGER.NE.0) WRITE(6,*)
     & ' THIS BLADE SECTION IS TO BE RE-STAGGERED'
C
      WRITE(6,*) ' IF_LEAN = ', IF_LEAN
      IF(IF_LEAN.NE.0) WRITE(6,*)
     & ' THIS BLADE SECTION IS TO BE LEANED'
C
C******************************************************************************
C     DO NOT READ IN ANY GEOMETRY IF "IF_DESIGN" IS NON ZERO. 
C     INSTEAD CALL SUBROUTINE RE_DESIGN TO GENERATE A NEW SECTION.
C     THEN JUMP TO 1595 .
C
      IF(IF_DESIGN.NE.0) THEN
          WRITE(6,*)'CALLING RE_DESIGN TO GENERATE A NEW BLADE GEOMETRY'
          CALL RE_DESIGN(NR,K,J1,J2,JLEROW,JTEROW) 
          GO TO 1595
      END IF
C
C******************************************************************************
C******************************************************************************
C
C   NOW INPUT THE BLADE GEOMETRY ON THE QUASI STREAM SURFACE .
C******************************************************************************
C******************************************************************************
C
C******************************************************************************
C***********XSURF(J,K) IS AXIAL COORDINATE OF POINTS ON THE STREAMWISE SURFACE.
C
      WRITE(6,*)  ' ROW NUMBER ',NR,' SECTION NUMBER ',K,
     &            ' X  COORDINATES '      
      READ(5,*)     FAC1, XSHIFT
      WRITE(6,*)  ' FAC1, XSHIFT = ', FAC1,XSHIFT
      WRITE(6,*)
      READ(5,*)     (XSURF(J,K),J=J1,J2)
      WRITE(6,1700) (XSURF(J,K),J=J1,J2)
      WRITE(6,*)
      XRANGE = XSURF(J2,K) - XSURF(J1,K)      
C
C***********RT_UPP(J,K) IS THE R-THETA COORDINATE OF POINTS ON THE STREAMWISE
C           SURFACE ON THE BLADE SURFACE WITH LARGEST VALUE OF THETA.
C           i.e. THE UPPER SURFACE OF THE BLADE AND THE LOWER SURFACE OF THE
C           BLADE TO BLADE PASSAGE.
C
      WRITE(6,*)  ' ROW NUMBER ',NR,' SECTION NUMBER ', K,
     &            ' R_THETA OF UPPER BLADE SURFACE'
      READ(5,*)     FAC2, TSHIFT
      WRITE(6,*)  ' FAC2, TSHIFT = ', FAC2,TSHIFT
      WRITE(6,*) 
      READ(5,*)     (RT_UPP(J,K),J=J1,J2)
      WRITE(6,1700) (RT_UPP(J,K),J=J1,J2)
      WRITE(6,*)
C
C***********RT_THICK(J,K) IS BLADE THICKNESS DELTA R-THETA MEASURED
C          IN THE TANGENTIAL DIRECTION.
C 
      WRITE(6,*)  ' ROW NUMBER ',NR,' SECTION NUMBER ',K,
     &            ' BLADE TANGENTIAL THICKNESS. '
      READ(5,*)     FAC3
      WRITE(6,*)  ' FAC3= ', FAC3
      WRITE(6,*)
      READ(5,*)     (RT_THICK(J,K),J=J1,J2)
      WRITE(6,1700) (RT_THICK(J,K),J=J1,J2)
      WRITE(6,*)
C
C      RSURF(J,K) IS THE RADII OF POINTS ON THE STREAMWISE SURFACES ON WHICH DATA IS INPUT.
C
      WRITE(6,*)   ' ROW NUMBER ',NR,' SECTION NUMBER ',K,
     &             ' STREAM SURFACE RADIUS. '
      READ(5,*)      FAC4, RSHIFT
      WRITE(6,*)   ' FAC4, RSHIFT  = ', FAC4,RSHIFT
      WRITE(6,*)
      READ(5,*)     (RSURF(J,K),J=J1,J2)
      WRITE(6,1700) (RSURF(J,K),J=J1,J2)
      WRITE(6,*)
C    END OF ALL GEOMETRICAL INPUT FOR THIS BLADE SECTION
C
C******************************************************************************
C******************************************************************************
C     CHECK THAT THE INPUT MERIDIONAL GRID SPACINGS ARE NOT TOO SMALL
C     AND ADJUST THEM IF THEY ARE.
C
      SDIST(J1) = 0.0
      DO 1450 J= J1+1,J2
           XD = XSURF(J,K)  - XSURF(J-1,K)
           RD = RSURF(J,K)  - RSURF(J-1,K)
           SDIST(J) = SDIST(J-1) + SQRT(XD*XD+RD*RD)
 1450 CONTINUE
      SRANGE = SDIST(J2) - SDIST(J1)
      SLIM   = 0.0001*SRANGE
      DO 1460 J = J1+2,J2
        IF((SDIST(J)- SDIST(J-1)).LT.SLIM) THEN
           WRITE(6,*) 
     &   ' WARNING, THE INPUT GRID SPACINGS ARE TOO SMALL AT J= ',
     &     J,'K= ',K
           WRITE(6,*) ' ADJUSTING THE SPACINGS AND CONTINUING. '
           XD = SLIM*(XSURF(J-2,K)-XSURF(J-1,K))/(SDIST(J-2)-SDIST(J-1))
           RD = SLIM*(RSURF(J-2,K)-RSURF(J-1,K))/(SDIST(J-2)-SDIST(J-1))
           XSURF(J,K) = XSURF(J-1,K) + XD
           RSURF(J,K) = RSURF(J-1,K) + RD
        END IF
 1460 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     SCALE AND SHIFT THE INPUT BLADE SECTIONS BY XSHIFT, TSHIFT or RSHIFT.
C
      DO 1500 J=J1,J2
      XSURF(J,K)    = FAC1*(XSHIFT + XSURF(J,K))
      RT_THICK(J,K) = FAC3*RT_THICK(J,K)
      RSURF(J,K)    = FAC4*(RSHIFT + RSURF(J,K))
      RT_UPP(J,K)   = FAC2*(TSHIFT + RT_UPP(J,K))
 1500 CONTINUE
C
C******************************************************************************
C****************************************************************************** 
C   RE ENTER HERE IF A NEW SECTION WAS DESIGNED,  i.e  IF  IF_DESIGN WAS NON ZERO.
C
 1595 CONTINUE
C
C******************************************************************************
C******************************************************************************
C    NOW START TO ROTATE OR RESTAGGER THIS BLADE SECTION .
C
C******************************************************************************
C******************************************************************************
C     RESTAGGER THE BLADE SECTION IF IF_STAGGER IS NON ZERO.
C     POSITIVE  "ROTATE"  MEANS CLOCKWISE ROTATION.
C
      IF(IF_STAGGER.NE.0 ) THEN
           WRITE(6,*) ' RESTAGGERING SECTION NUMBER ',K
           CALL RESTAGGER(K,J1,J2,JLEROW,JTEROW,ROTATE,FRACX_ROT)
      END IF
C
C******************************************************************************
C******************************************************************************
C     LEAN THE BLADE BY ANGLEAN IF "IF_LEAN" IS GREATER THAN ZERO.
C     IF "ANGLEAN" IS POSITIVE THE HUB IS HELD FIXED AND THE OTHER SECTIONS ARE 
C     MOVED IN THE POSITIVE THETA DIRECTION .
C
      IF(IF_LEAN.NE.0) THEN
      WRITE(6,*) '  LEANING SECTION NUMBER ', K
      CALL LEAN(K,J1,J2,ANGLEAN)
      END IF
C
C************************** ****************************************************
C******************************************************************************
C   Q3D
C     CALL SET SSTHICK TO SET THE STREAM SURFACE THICKNESS IF DOING A Q3D
C     BLADE TO BLADE CALCULATION. THEN JUMP TO 1605 TO END THE GEOMETRY INPUT
C     FOR THIS BLADE ROW AS ONLY A SINGLE STREAM SURFACE IS USED FOR Q3D.
C
      IF(KM.EQ.2) THEN
           CALL SET_SSTHICK(J1,J2)
C    JUMP OUT OF THE DO 1600 LOOP AS ONLY ONE STREAM SURFACE IS NEEDED.
           GO TO 1605
      END IF
C   END Q3D
C
C******************************************************************************
C******************************************************************************
C     END OF ALL INPUT FOR THIS BLADE SECTION
 1600 CONTINUE
C
C  RE-ENTER HERE IF  KM = 2
 1605 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     THE BLADE GEOMETRY HAS NOW BEEN INPUT FOR ALL SECTIONS.
C     ALL THE FOLLOWING INPUT APPLIES TO THE WHOLE ROW.
C*******************************************************************************
C*******************************************************************************
C
C     SET UP NEW HUB AND CASING COORDINATES AND BLADE SECTIONS IF INSURF > 0 .
C     USE THE INPUT COORDINATES AND SECTIONS IF INSURF = 0 .
C
C*******************************************************************************
C*******************************************************************************
C
      IF(INSURF.EQ.1.OR.INSURF.GT.2) THEN
           READ(5,*) NHUB
           READ(5,*) (XHUB(N),N=1,NHUB)
           WRITE(6,*) ' XHUB = ', (XHUB(N),N=1,NHUB)
           READ(5,*) (RHUB(N),N=1,NHUB)
           WRITE(6,*) ' RHUB = ', (RHUB(N),N=1,NHUB)
           DO N = 1,NHUB
               XHUB(N) = (XSHIFT + XHUB(N) ) * FAC1
               RHUB(N) = (RSHIFT + RHUB(N) ) * FAC4
           END DO
      END IF
C
      IF(INSURF.GE.2) THEN
           READ(5,*) NTIP
           READ(5,*) (XTIP(N),N=1,NTIP)
           WRITE(6,*) ' XTIP =', (XTIP(N),N=1,NTIP)
           READ(5,*) (RTIP(N),N=1,NTIP)
           WRITE(6,*) ' RTIP =', (RTIP(N),N=1,NTIP)
           DO N = 1,NTIP
               XTIP(N) = (XSHIFT + XTIP(N) ) * FAC1
               RTIP(N) = (RSHIFT + RTIP(N) ) * FAC4
           END DO
      END IF
C
C******************************************************************************
C******************************************************************************
      DO 1650 J = J1,J2
C
      DO 1660 K = 1,NSECS_ROW
           XQO(K)  = XSURF(J,K)
           RQO(K)  = RSURF(J,K)
 1660 CONTINUE
C
C      CALLING  INSECT  TO FIND THE INTERSECTION OF THE HUB CURVE WITH
C      THE BLADE QO CURVE.
C
      IF(INSURF.EQ.1.OR.INSURF.GT.2) THEN
            CALL INSECT(JD,MAXKI,NHUB,XHUB,RHUB,NSECS_ROW,XQO,RQO,
     &           XNEWHUB(J),RNEWHUB(J))
      ELSE
            XNEWHUB(J) = XSURF(J,1)
            RNEWHUB(J) = RSURF(J,1)
      END IF
C
C      CALLING  INSECT  TO FIND THE INTERSECTION OF THE  CASING CURVE WITH
C      THE BLADE QO CURVE.
C
      IF(INSURF.GE.2) THEN
            CALL INSECT(JD,MAXKI,NTIP,XTIP,RTIP,NSECS_ROW,XQO,RQO,
     &           XNEWTIP(J),RNEWTIP(J))
      ELSE
            XNEWTIP(J) = XSURF(J,NSECS_ROW)
            RNEWTIP(J) = RSURF(J,NSECS_ROW)
      END IF
C
 1650 CONTINUE
C
C   NEW HUB AND CASING CORDINATES HAVE BEEN FOUND.
C
C******************************************************************************
C******************************************************************************
C  INTERPOLATE TO OBTAIN  "NSECS_IN"  NEW BLADE SECTIONS EQUALLY SPACED BETWEEN 
C  THE HUB AND CASING.
C
      DO 1670 J = J1,J2
C
      XSPAN  = XNEWTIP(J) - XNEWHUB(J)
      RSPAN  = RNEWTIP(J) - RNEWHUB(J)
C
      DO 1671 K = 1,NSECS_ROW
           XINT1(K)  = XSURF(J,K)
           XINT2(K)  = RSURF(J,K)
           XINT3(K)  = RT_UPP(J,K)
           XINT4(K)  = RT_THICK(J,K)
 1671 CONTINUE
C
C  FIRST IF THE SPANWISE DIRECTION IS PREDOMINANTLY AXIAL.
C
      IF(ABS(XSPAN).GT.ABS(RSPAN) ) THEN
C
      DO 1672 K = 1, NSECS_IN
      F_SPAN     = FLOAT(K-1)/FLOAT(NSECS_IN -1)
      XARG       = XNEWHUB(J) + F_SPAN*(XNEWTIP(J) - XNEWHUB(J) )
      XSURF(J,K) = XARG 
      CALL INTP(NSECS_ROW,XINT1,XINT2,XARG,RSURF(J,K) )
      CALL INTP(NSECS_ROW,XINT1,XINT3,XARG,RT_UPP(J,K) )
      CALL INTP(NSECS_ROW,XINT1,XINT4,XARG,RT_THICK(J,K) )
 1672 CONTINUE
C
      ELSE    
C
C     NEXT IF THE SPANWISE DIRECTION IS PREDOMINANTLY RADIAL
C
      DO 1673 K = 1, NSECS_IN
      F_SPAN     = FLOAT(K-1)/FLOAT(NSECS_IN -1)
      RARG       = RNEWHUB(J) + F_SPAN*(RNEWTIP(J) - RNEWHUB(J) )
      RSURF(J,K) = RARG 
      CALL INTP(NSECS_ROW,XINT2,XINT1,RARG,XSURF(J,K) )
      CALL INTP(NSECS_ROW,XINT2,XINT3,RARG,RT_UPP(J,K) )
      CALL INTP(NSECS_ROW,XINT2,XINT4,RARG,RT_THICK(J,K) )
 1673 CONTINUE
C
      END IF 
C 
 1670 CONTINUE
C
C   FINISHED INTERPOLATING IN INPUT SECTIONS TO FORM NSECS_IN SECTIONS.
C
C************************************************** ****************************
C******************************************************************************
C     INTERPOLATE IN THE UPSTREAM AND DOWNSTREAM GRID ANGLES IF IF_ANGLES > 0
C
      IF(IF_ANGLES(NR).GT.0) THEN
C
      WRITE(6,*)  ' INPUTTING THE UPSTREAM AND DOWNSTREAM GRID ANGLES'
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    N_ANGLES
      READ(5,*)   (FRACN_SPAN(NK),NK=1,N_ANGLES)
      READ(5,*)   (ANGL_UP(NK),   NK=1,N_ANGLES)
      READ(5,*)   (ANGL_DWN1(NK), NK=1,N_ANGLES)
      READ(5,*)   (ANGL_DWN2(NK), NK=1,N_ANGLES)
C
      QSPAN(1) = 0.0
      J_MID     = J1 + (JLEROW + JTEROW)/2
      DO 1680 K = 2,NSECS_IN
           XDIF = XSURF(J_MID,K) - XSURF(J_MID,K-1)
           RDIF = RSURF(J_MID,K) - RSURF(J_MID,K-1)
           QSPAN(K) = QSPAN(K-1) + SQRT(XDIF*XDIF + RDIF*RDIF)
 1680 CONTINUE
      DO 1681 K = 1,NSECS_IN
      QSPAN(K) = QSPAN(K)/QSPAN(NSECS_IN)
 1681 CONTINUE 
C
      DO 1685 K = 1,NSECS_IN
      ARG = QSPAN(K)
      CALL INTP(N_ANGLES,FRACN_SPAN,ANGL_UP,ARG,BETAUP(NR,K))
      CALL INTP(N_ANGLES,FRACN_SPAN,ANGL_DWN1,ARG,BETADWN1(NR,K))
      CALL INTP(N_ANGLES,FRACN_SPAN,ANGL_DWN2,ARG,BETADWN2(NR,K))
C
      WRITE(6,*)'INTERPOLATED GRID ANGLES ON THE NEW STREAM SURFACES',
     &           BETAUP(NR,K),BETADWN1(NR,K),BETADWN2(NR,K)
 1685 CONTINUE
      WRITE(6,*)
C
      END IF
C
C******************************************************************************
C******************************************************************************
C  TFLOW
C******************************************************************************
C******************************************************************************
C  IF A THROUGHFLOW CALCULATION READ IN THE REQUIRED EXIT FLOW ANGLE OR DEVIATION 
C  ANGLE FOR THIS BLADE ROW. ANGL_TYP = 'A' FOR EXIT ANGLE, = 'D' FOR DEVIATION ANGLE.
      IF(IM.EQ.2) THEN
      READ(5,*)    DUMMY_INPUT
      WRITE(6,*) ' INPUTTING THE THROUGHFLOW DATA '
      READ(5,*)    ANGL_TYP(NR), NANGLES(NR)
      READ(5,*)  ( FRAC_SPAN(NR,NK),NK=1,NANGLES(NR))
      READ(5,*)  ( EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
      WRITE(6,*) ' NANGLES = ',NANGLES(NR),' ANGL_TYP = ',ANGL_TYP(NR)
      WRITE(6,9048) (FRAC_SPAN(NR,NK),NK=1,NANGLES(NR))
      IF(ANGL_TYP(NR).EQ.'A')
     &     WRITE(6,9045) (EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
      IF(ANGL_TYP(NR).EQ.'D')
     &     WRITE(6,9046) (EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
 9048 FORMAT(' FRACTION OF SPAN',/,(10F10.3))
 9045 FORMAT(' EXIT FLOW ANGLE IN DEGREES ',/,(10F10.3))
 9046 FORMAT(' EXIT FLOW ANGLE DEVIATION FROM GRID ANGLE, IN DEGREES.',
     & /,(10F10.3))
      WRITE(6,*) ' END THROUGHFLOW DATA '
      END IF
C ******************************************************************************
C  END TFLOW
C******************************************************************************
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
C      CALL NEWGRID FOR THE CURRENT BLADE ROW IF REQUESTED
C      THIS READS IN MORE DATA AND GENERATES NEW STREAMWISE (J) GRID POINTS.
C
      IF(NEW_GRID.NE.0) THEN
            CALL NEWGRID(JMROW,J1,J2,
     &                   JLEROW,JTEROW,JROTHS,JROTHE,JROTTS,JROTTE)
      END IF
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
C      SET THE OVERALL  "J"  MARKERS TO BE RELATIVE TO THE FIRST POINT ON THE FIRST ROW.
C
      JLE(NR) = J1+JLEROW-1
      JTE(NR) = J1+JTEROW-1
      J2      = J1 -1 + JMROW
      JROW    = 0
C
      DO 100 J = J1,J2
      JM1      = J-1
      IF(J.EQ.1) JM1=1
      JROW     = JROW+1
C
C********************************************************************************
C      SET THE ROTATION OF THE HUB AND CASING.
C      THE DEFAULT IF JROTHS AND JROTTS ARE ZERO
C      IS THAT THE WHOLE HUB IS ROTATING AT THE SAME SPEED AS
C      THE BLADE ROW AND THE CASING IS STATIONARY.
C
      WRAD(J) = RPMROW*PI/30.
      WHUB(J) = 0.0
      IF(JROW.GT.JROTHS.AND.JROW.LT.JROTHE) WHUB(J) = RPMHUB*PI/30.
      WTIP(J) = 0.0
      IF(JROW.GT.JROTTS.AND.JROW.LT.JROTTE) WTIP(J) = WRAD(J)
C
C     SET LIMITS OF SOLID PART OF THE BLADE ELEMENTS, IE WHERE NO FLOW.
C     THESE ARE FROM 1 TO KTIPS-1  OR FROM KTIPE TO KM-1
C
      IF(KTIPS(NR).LE.0) THEN
           KS1(NR) = 1
           KS2(NR) = KMM1
      ENDIF
      IF(KTIPS(NR).EQ.1) THEN
           KS1(NR) = KTIPE(NR)
           KS2(NR) = KMM1
      ENDIF
      IF(KTIPS(NR).GT.1) THEN
           KS1(NR)   = 1
           KS2(NR)   = KTIPS(NR)-1
           KTIPE(NR) = KM
      ENDIF
C
C     SET SOME MARKERS
C
      IND(J)    = 0
      INDLE(J)  = 0
      INDTE(J)  = 0
      INDMIX(J) = 0
      INDMID(J) = 0
      NBLADE(J) = NBLADES_IN_ROW
      NROW(J)   = NR
      JMIDRW    = 0.5*(JLEROW+JTEROW)
      IF(JROW.GT.JLEROW.AND.JROW.LE.JTEROW) IND(J) = 1
      IF(JROW.EQ.JLEROW) INDLE(J)    = 1
      IF(JROW.EQ.JLEROW) JLED(NR)    = J
      IF(JROW.EQ.JTEROW) INDTE(J)    = 1
      IF(JROW.EQ.JLEROW) INDLE(J-1) = -1
      IF(JROW.EQ.JLEROW) INDLE(J-2) = -2
      IF(JROW.EQ.JMIDRW) INDMID(J)   = 1
C
C     SET THE INITIAL GUESS OF STATIC PRESSURE
C
      IF(JROW.LE.JLEROW) PGUESS(J) = PUPROW+(PLEROW-PUPROW)*
     & (JROW-1)/(JLEROW-1)
      IF(JROW.GT.JLEROW.AND.JROW.LE.JTEROW) PGUESS(J) = 
     &          PLEROW + (PTEROW-PLEROW)*(JROW-JLEROW)/(JTEROW-JLEROW)
      IF(JROW.GT.JTEROW) PGUESS(J) =
     &          PTEROW + (PDNROW-PTEROW)*(JROW-JTEROW)/(JMROW-JTEROW)
  100 CONTINUE
C
       INDMIX(J2) = 1
       JMIX(NR)   = J2
       JSTART(NR) = J1
C
C     RESET THE J INDEX TO CONTINUE THROUGH THE NEXT BLADE ROW, IF ANY.
C
      J1 = J2 + 1
C
C**********************************************************************************
C        END OF DATA INPUT ON THIS BLADE ROW. RETURN TO INPUT DATA ON THE
C        NEXT ROW UNLESS THIS IS THE LAST, IE UNLESS NR = NROWS.
C**********************************************************************************
C**********************************************************************************
 1550 CONTINUE
C**********************************************************************************
C**********************************************************************************
C
C   SET THE TOTAL NUMBER OF "J" POINTS = JM.
C
      JM         = J2
      JMM1       = JM-1
      JMM2       = JM-2
      INDMIX(JM) = 0
C
C     CHECK THAT JM IS NOT TOO LARGE.
C
      IF(JM.GT.JD)  WRITE(6,*) 'STOPPING BECAUSE JM TOO LARGE.',
     &            ' JM= ',JM, ' DIMENSION LIMIT = ',JD
      IF(JM.GT.JD) STOP
C**********************************************************************************
C     SET THE START AND END POINTS OF EACH STAGE
C
      JSTG_START(1)   = 1
      NSTAGES         = NSTAGE(NROWS)
      NSTAGE(NROWS+1) = NSTAGES + 1
      JSTART(NROWS+1) = JM +1
C
      DO N = 1,NROWS
      NSTG   = NSTAGE(N)
      NSTGP1 = NSTAGE(N+1)
      IF(NSTG.NE.NSTGP1) THEN
             JSTG_END(NSTG)     = JMIX(N)
             JSTG_START(NSTGP1) = JSTART(N+1)
      END IF
      WRITE(6,*) 'ROW NUMBER ',N,'STAGE NUMBER ',NSTG,'JSTART ',
     &             JSTART(N),'JEND ',JMIX(N)
      END DO
C
      DO N = 1,NSTAGES
      WRITE(6,*) ' STAGE NUMBER ',N,' JSTART= ',JSTG_START(N),
     &           ' JEND = ',JSTG_END(N)
      END DO
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C      AUTOMATICALLY SHIFT THE BLADES TO MAKE THE  X  SPACING CONTINUOUS
C      BETWEEN THE LAST GRID POINT ON ONE BLADE AND THE FIRST GRID
C      POINT ON THE NEXT ONE ON THE HUB STREAMLINE. IF ISHIFT = 1 OR > 2.
C**********************************************************************************
C**********************************************************************************
C    START TO SET THE GRID BETWEEN BLADE ROWS AND MAKE THE MIXING PLANES CONTIGUOUS
C**********************************************************************************
C
C     SKIP THE NEXT PART AND DO NOT MOVE THE BLADES OR GRID IF ISHIFT = 0 .
      IF(ISHIFT.EQ.0) GO TO 140
C     
      IF(ISHIFT.GE.2) GO TO 133
C
C      NEXT SECTION ONLY FOR ISHIFT = 1
C      IF  ISHIFT = 1  AUTOMATICALLY SHIFT THE BLADES TO MAKE THE X SPACING CONTINUOUS
C      BETWEEN THE LAST GRID POINT ON ONE BLADE AND THE FIRST GRID POINT ON THE NEXT
C      ONE ON THE HUB STREAMLINE BUT DO NOT CHANGE THE GRID SPACINGS OR THE RADII .
C
      XSHIFT = 0.0
      DO 130 J = 2,JM
      IF(INDMIX(J-1).NE.1) GO TO 129
      DIFF1  = XSURF(J-1,1) - XSURF(J-2,1)
      DIFF2  = XSURF(J+1,1) - XSURF(J,1)
      XSHIFT = XSURF(J-1,1) + 0.01*(DIFF1+DIFF2) - XSURF(J,1)
  129 CONTINUE
      DO 131 K=1,KM
      XSURF(J,K) = XSURF(J,K) + XSHIFT
      X(J,K)     = X(J,K)     + XSHIFT
  131 CONTINUE
  130 CONTINUE
C
C     MAKE NO MORE CHANGES TO THE GRID AND JUMP TO 140 IF ISHIFT = 1.
      IF(ISHIFT.EQ.1) GO TO 140
C
C**********************************************************************************
C     RE ENTER HERE IF ISHIFT = 2
  133 CONTINUE
C
C******************************************************************************C
C******************************************************************************C
C******************************************************************************C
C     START THE LOOP TO GRID THE GAPS BETWEEN BLADE ROWS AND ALSO THE GRIDS
C     UPSTREAM OF THE FIRST ROW AND DOWNSTREAM OF THE LAST ROW.
C     NROW_GAP  IS THE NUMBER OF REGIONS TO BE GRIDDED, INCLUDING UPSTREAM OF THE
C     FIRST ROW, ALL INTER ROW GAPS AND DOWNSTREAM OF THE LAST ROW. 
C     SO THERE ARE NROWS + 1  GAPS .
C
      DO 6666 NROW_GAP = 1,NROWS+1
C
      IF(NROW_GAP.EQ.1) THEN
           NRW  = 1
           JST  = 1
           JEND = JLE(1)
      END IF
      IF(NROW_GAP.GT.1.AND.NROW_GAP.NE.NROWS+1) THEN
          NRW  = NROW_GAP
          JST  = JTE(NRW-1)
          JEND = JLE(NRW)
          JMID = JMIX(NRW-1)
      END IF
      IF(NROW_GAP.EQ.NROWS+1) THEN
           NRW  = NROWS + 1
           JST  = JTE(NROWS)
           JEND = JM
      END IF
C
      NRWM1 = NRW-1
C******************************************************************************C
C     NOTE THAT THIS SECTION SETS THE GRID IN THE GAP BETWEEN BLADE ROWS ALSO
C     UPSTREAM OF THE FIRST ROW AND DOWNSTREAM OF THE LAST ROW.
C     "NRW"  IS THE ROW NUMBER OF THE ROW DOWNSTREAM OF THE GAP BEING GRIDDED.
C******************************************************************************C
      DO 8000 K=1,NSECS_IN
C
C    SET THE MERIDIONAL DISTANCE SMERID
C
      SMERID(1,K) = 0.0
      DO 213 J=2,JM
      XD = XSURF(J,K)  - XSURF(J-1,K)
      RD = RSURF(J,K)  - RSURF(J-1,K)
      SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)
  213 CONTINUE
C
C    IF ISHIFT = 2 MAINTAIN THE MERIDIONAL CURVE IN THE GAP BETWEEN BLADE ROWS.
C    IF ISHIFT = 3 MAKE THE MERIDIONAL VIEW OF THE GRID LINEAR IN THE GAP BETWEEN BLADE ROWS.
C    IF ISHIFT = 4 SAME AS 3 BUT DO NOT CHANGE THE HUB AND CASING PROFILES
C
      DLAST= -1.0
      NPOINTS   = 1
      XGAP = XSURF(JEND,K)  - XSURF(JST,K)
      RGAP = RSURF(JEND,K)  - RSURF(JST,K)
      GAP  = SQRT(XGAP*XGAP + RGAP*RGAP)
C
      DO 110 J = JST,JEND
      XD   = XSURF(J,K)  - XSURF(JST,K)
      RD   = RSURF(J,K)  - RSURF(JST,K)
      PROJ = (XD*XGAP + RD*RGAP)
      IF(ISHIFT.GE.3) PROJ = PROJ/GAP
      DIST = SQRT(XD*XD+RD*RD)
      IF(PROJ.LT.0.0) GO TO 110
C
      IF(DIST.GT.(1.00001*GAP)) THEN
           XINT(NPOINTS) = XSURF(JEND,K)
           RINT(NPOINTS) = RSURF(JEND,K)
           NPOINTS = NPOINTS+1
      GO TO 111
      ENDIF
C
      IF(DIST.LT.DLAST) GO TO 110
C
      DLAST     = DIST
C
      IF(ISHIFT.EQ.3) THEN
           XINT(NPOINTS)  = XSURF(JST,K)   + PROJ*XGAP/GAP
           RINT(NPOINTS)  = RSURF(JST,K)   + PROJ*RGAP/GAP
      ENDIF
C
      IF(ISHIFT.EQ.2) THEN
           XINT(NPOINTS)  = XSURF(J,K)
           RINT(NPOINTS)  = RSURF(J,K)
      ENDIF
C
      IF(ISHIFT.EQ.4) THEN
           IF(K.EQ.1.OR.K.EQ.NSECS_IN) THEN
                XINT(NPOINTS)  = XSURF(J,K)
                RINT(NPOINTS)  = RSURF(J,K)
           ELSE
                XINT(NPOINTS)  = XSURF(JST,K)   + PROJ*XGAP/GAP
                RINT(NPOINTS)  = RSURF(JST,K)   + PROJ*RGAP/GAP
           END IF
      END IF
C
      NPOINTS = NPOINTS+1
C
  110 CONTINUE
C
  111 CONTINUE
C******************************************************************************C
C    SMOOTH THE X AND R COORDINATES IN THE BLADE TO BLADE GAP.
C
      NGAPS   = NPOINTS - 1
      NSMOOTH = 10
      SFGAP   = 0.2
      DO 113 NS = 1,NSMOOTH
      DO 114 NN = 2,NGAPS-1
      XINT(NN) = (1.-SFGAP)*XINT(NN) + SFGAP*0.5*(XINT(NN-1)+XINT(NN+1))
      RINT(NN) = (1.-SFGAP)*RINT(NN) + SFGAP*0.5*(RINT(NN-1)+RINT(NN+1))
  114 CONTINUE
  113 CONTINUE
C******************************************************************************C
C     RE CALCULATE THE MERIDIONAL DISTANCE  "SDIST "  AS IT WILL HAVE BEEN CHANGED
C     BY THE SMOOTHING.
C
      SDIST(1)  = SMERID(JST,K)
      DO 112 NN = 2,NGAPS
      XD = XINT(NN) - XINT(NN-1)
      RD = RINT(NN) - RINT(NN-1)
      SDIST(NN) = SDIST(NN-1) + SQRT(XD*XD+RD*RD)
  112 CONTINUE
C
C
      SMID     = 0.5*(SDIST(1) + SDIST(NGAPS))
      STE      = SDIST(1)
      SLE      = SDIST(NGAPS)
      ANGLUP   = BETAUP(NRW,K)
      IF(NRW.EQ.1) THEN
           ANGLDWN1 = BETADWN1(1,K)
           ANGLDWN2 = BETADWN2(1,K)
      ELSE
           ANGLDWN1 = BETADWN1(NRWM1,K)
           ANGLDWN2 = BETADWN2(NRWM1,K)
      END IF
      IF(NRW.EQ.NROWS+1) THEN
           ANGLDWN1 = BETADWN1(NROWS,K)
           ANGLDWN2 = BETADWN2(NROWS,K)
      END IF
C
C******************************************************************************
C    START TO CALL GRID_UP AND GRID_DOWN TO SET THE GRID UPSTREAM AND DOWNSTREAM OF
C    ALL BLADE ROWS. ALSO SETS THE TRAILING EDGE CUSPS.
C******************************************************************************
C    CALL GRID_UP TO FORM THE GRID UPSTREAM OF THE FIRST BLADE ROW.
C
            WRITE(6,*)
      IF(NRW.EQ.1.AND.NRW.NE.NROWS+1.AND.K.EQ.1) THEN
            WRITE(6,*)
     &'****************************************************************'
      WRITE(6,*) ' CALLING GRID_UP AND GRID_DOWN TO SET THE GRID AND ANY
     & CUSPS IN THE GAPS BETWEEN BLADE ROWS.'
            WRITE(6,*)
     &'****************************************************************'
            WRITE(6,*)
            WRITE(6,*)  ' STARTING TO CALL GRID_UP AND GRID_DOWN FOR ROW
     & NUMBER', 1
            WRITE(6,*)
      END IF
C
      IF(NRW.EQ.1) THEN
           WRITE(6,*)
     &    'CALLING GRID_UP FOR ROW 1, K = ',K, ' JSTART,JLE1= ',JST,JEND
           CALL GRID_UP(K,JST,JEND,STE,SLE,NGAPS,SDIST,XINT,
     &                  RINT,NEXTRAP_LE,ANGLUP,IF_ANGLES(1) )
      END IF
C
      IF(NRW.EQ.NROWS+1) GO TO 7900
C
C     CALL GRID_UP AND GRID_DOWN TO FORM THE GRIDS UPSTREAM AND DOWNSTREAM OF THE
C     MIXING PLANE FOR INTERIOR BLADE ROWS.
C
      IF( (NRW.NE.1.AND.NRW.NE.NROWS+1).AND.K.EQ.1) THEN
            WRITE(6,*)
            WRITE(6,*)
     &'****************************************************************'
            WRITE(6,*)  ' STARTING TO CALL GRID_UP AND GRID_DOWN FOR ROW
     & NUMBER', NRW
            WRITE(6,*)
      END IF
C
      IF (NRW.NE.1)  THEN
C
       WRITE(6,*)
     &' CALLING GRID_DOWN FOR ROW No.',NRW,'K = ',K,' JTE,JMID= ',
     &  JST,JMID, 'IF_CUSP=', IF_CUSP(NRWM1)
      CALL GRID_DOWN(K,JST,JMID,STE,SMID,NGAPS,SDIST,XINT,RINT,
     &        NEXTRAP_TE,ANGLDWN1,ANGLDWN2,IF_CUSP(NRWM1),
     &        ICUSP(NRWM1),LCUSP(NRWM1),LCUSPUP(NRWM1),IF_ANGLES(NRWM1))
C
      WRITE(6,*)
     &' CALLING GRID_UP FOR ROW No.  ',NRW,'K = ',K,'JMID, JLE= ',
     &  JMID+1,JEND
      CALL GRID_UP(K,JMID+1,JEND,SMID,SLE,NGAPS,SDIST,XINT,RINT,
     &             NEXTRAP_LE,ANGLUP,IF_ANGLES(NRW))
C
      ENDIF
C
      GO TO 8000
C
 7900 CONTINUE
C
C     CALL GRID_DOWN TO FORM THE GRID DOWNSTREAM OF THE LAST BLADE ROW.
C
      IF(NRW.EQ.NROWS+1.AND.K.EQ.1) THEN
            WRITE(6,*)
     &'****************************************************************'
            WRITE(6,*)'CALLING GRID_DOWN THE LAST BLADE ROW, ROW NUMBER'
     &                 ,NROWS
            WRITE(6,*)
      END IF
C
      IF(NRW.EQ.NROWS+1)  THEN
C
      WRITE(6,*)
     & ' CALLING GRID_DOWN FOR THE LAST BLADE ROW, K = ',K,'JTE,JEND= ',
     &   JST,JEND, 'IF_CUSP=', IF_CUSP(NROWS)
      CALL GRID_DOWN(K,JST,JEND,STE,SLE,NGAPS,SDIST,
     &        XINT,RINT,NEXTRAP_TE,ANGLDWN1,ANGLDWN2,IF_CUSP(NROWS),
     &        ICUSP(NROWS),LCUSP(NROWS),LCUSPUP(NROWS),IF_ANGLES(NROWS))
      END IF
C
 8000 CONTINUE
C
      WRITE(6,*)
C
C*******************************************************************************
C*******************************************************************************
C
C     IF  JEND = JM GAP GRID GENERATION COMPLETED
C     
      IF(JEND.EQ.JM) THEN
             WRITE(6,*) 'END OF SETTING THE GRID AND CUSPS UPSTREAM AND 
     & DOWNSTREAM OF ALL BLADE ROWS'
             WRITE(6,*)
     &'****************************************************************'
             WRITE(6,*)
     &'****************************************************************'
             WRITE(6,*)
      END IF
C
C******************************************************************************
C     JUMP BACK TO 6666 TO START ON THE NEXT ROW UNLESS JEND = JM.
C
 6666 CONTINUE
C
C******************************************************************************C
C******************************************************************************C
C      RESET SMERID  WHICH WILL HAVE BEEN CHANGED BY THE NEW GRID.
C
      DO 216 K=1,NSECS_IN
      DO 216 J=2,JM
      XD = XSURF(J,K) - XSURF(J-1,K)
      RD = RSURF(J,K) - RSURF(J-1,K)
  216 SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)
C
C******************************************************************************
C******************************************************************************
C   END OF ALL GEOMETRY INPUT AND MANIPULATION.
C
  140 CONTINUE
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
C    READ IN THE SPANWISE VARIATION IN INLET AND EXIT BOUNDARY CONDITIONS.
C
C******************************************************************************
C******************************************************************************     
C******************************************************************************
C   INPUT THE NUMBER OF DATA POINTS USED FOR THE INLET BOUNDARY CONDITIONS.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    KIN
      WRITE(6,*) ' NUMBER OF DATA POINTS USED FOR THE INLET BOUNDARY'
      WRITE(6,*) ' CONDITIONS = ', KIN
      WRITE(6,*)
C
      IF(KIN.GT.MAXKI) WRITE(6,*) ' STOPPING BECAUSE KIN TOO LARGE.',
     &            ' KIN= ',KIN,' DIMENSION LIMIT = ',MAXKI
      IF(KIN.GT.MAXKI) STOP
C
      KIN_MID=IFIX(0.5*KIN)
C  Q3D
      IF(KM.EQ.2) KIN_MID = 1
C  END Q3D
C
C
C*******************************************************************************
C   INPUT THE RELATIVE SPANWISE SPACING OF THE POINTS WHERE THE INLET BOUNDARY
C   CONDITIONS ARE GIVEN'
C
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (FR_IN(K),K=1,KIN-1)
      WRITE(6,*)  ' RELATIVE SPACING OF THE POINTS WHERE THE INLET BOUND
     &ARY CONDITIONS ARE GIVEN'
      WRITE(6,1700) (FR_IN(K),K=1,KIN-1)
      WRITE(6,*)
C
C    INPUT THE INLET ABSOLUTE STAGNATION PRESSURE VARIATION WITH SPAN
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (PO1(K),K=1,KIN)
      WRITE(6,*)  ' INLET STAGNATION PRESSURE VARIATION WITH SPAN.'
      WRITE(6,1800) (PO1(K),K=1,KIN)
      PO_IN_MID = PO1(KIN_MID)
C
C    INPUT THE EXIT STATIC PRESSURE VARIATION WITH SPAN. ONLY IF  "IPOUT" = 3 .
      IF(IPOUT.EQ.3) THEN
           READ(5,*)     DUMMY_INPUT
           READ(5,*)    (PD(K),K=1,KIN)
           WRITE(6,*)  ' EXIT STATIC PRESSURE VARIATION WITH SPAN.'
           WRITE(6,1800)(PD(K),K=1,KIN)
      END IF
C
C     INPUT THE INLET ABSOLUTE STAGNATION TEMPERATURE VARIATION WITH SPAN.
      READ(5,*)       DUMMY_INPUT
      READ(5,*)      (TO1(K),K=1,KIN)
      WRITE(6,*)    ' INLET STAGNATION TEMPERATURE VARIATION WITH SPAN.'
      WRITE(6,1730)  (TO1(K),K=1,KIN)
      TO_IN_MID = TO1(KIN_MID)
C
C   INPUT THE INLET TANGENTIAL VELOCITY VARIATION WITH SPAN.
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (VTIN(K),K=1,KIN)
      WRITE(6,*)   ' INLET TANGENTIAL VELOCITY VARIATION WITH SPAN.'
      WRITE(6,1730) (VTIN(K),K=1,KIN)
      VT_IN_MID  = VTIN(KIN_MID)
C
C    INPUT THE INLET MERIDIONAL VELOCITY VARIATION WITH SPAN.
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (VM1(K),K=1,KIN)
      WRITE(6,*)   ' INLET MERIDIONAL VELOCITY VARIATION WITH SPAN.'
      WRITE(6,1730) (VM1(K),K=1,KIN)
C
C    INPUT THE INLET YAW ANGLE VARIATION WITH SPAN.
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (BS(K),K=1,KIN)
      WRITE(6,*)  ' INLET YAW ANGLE VARIATION WITH SPAN.'
      WRITE(6,1730) (BS(K),K=1,KIN)
      YAW_IN_MID = BS(KIN_MID)
C
C    INPUT THE INLET MERIDIONAL PITCH ANGLE VARIATION WITH SPAN.
      READ(5,*)      DUMMY_INPUT
      READ(5,*)     (BR(K),K=1,KIN)
      WRITE(6,*)   ' INLET MERIDIONAL PITCH ANGLE VARIATION WITH SPAN.'
      WRITE(6,1730) (BR(K),K=1,KIN)
      PITCH_IN_MID  = BR(KIN_MID)
C
C    END OF INPUTTING THE SPANWISE VARIATION OF BOUNDARY CONDITIONS.
C****************************************************************************** 
C****************************************************************************** 
C     INPUT THE EXIT STATIC PRESSURES ON THE HUB AND CASING.
C     THIS IS THE MAIN EXIT BOUNDARY CONDITION.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    PDOWN_HUB,PDOWN_TIP
      WRITE(6,*) ' SPECIFIED DOWNSTREAM PRESSURES ON HUB AND CASING=' ,
     &             PDOWN_HUB,PDOWN_TIP
      WRITE(6,*)
C
C*****************************************************************************
C******************************************************************************
C******************************************************************************
C     FORM THE FP(I) AND FR(K) AND FR_IN(K) INTO A GEOMETRIC SERIES IF REQUESTED.
C     BY SETTING FP(3) OR FR(3) OR FR_IN(3) TO ZERO.
C
      NUM3 = 3
C     
C   Q3D
      IF(KM.EQ.2) THEN
           FR(1)    = 1.0
           FR(2)    = 0.0
      ELSE
C   END Q3D
C
C   FIRST FOR THE SPANWISE GRID SPACINGS,  RF(K) .
      IF(FR(NUM3).GT.0.00001) GO TO 555
      FRRAT = FR(1)
      FRMAX = FR(2)
      FR(1) = 1.0
      FR(KM)  = 0.0
      DO 556 K = 2,KMM1
  556 FR(K) = FR(K-1)*FRRAT
      DO 557 K = 1,KMM1
      FREV  = FR(KM-K)
      IF(FREV.LT.FR(K))  FR(K) = FREV
      IF(FR(K).GT.FRMAX) FR(K) = FRMAX
  557 CONTINUE
C
  555 CONTINUE
C
      END IF
C
C    NEXT FOR THE BOUNDARY CONDITION SPACINGS, FR_IN(K) .
      IF(KIN.EQ.2) THEN
           FR_IN(1) = 1.0
           FR_IN(2) = 0.0
      ELSE
C
      IF(FR_IN(NUM3).GT.0.00001) GO TO 560
      FRRAT = FR_IN(1)
      FRMAX = FR_IN(2)
      FR_IN(1)   = 1.0
      FR_IN(KIN) = 0.0
      DO 558 K = 2,KIN-1
      FR_IN(K) = FR_IN(K-1)*FRRAT
  558 CONTINUE
      DO 559 K = 1,KIN-1
      FREV  = FR_IN(KIN-K)
      IF(FREV.LT.FR_IN(K))  FR_IN(K) = FREV
      IF(FR_IN(K).GT.FRMAX) FR_IN(K) = FRMAX
  559 CONTINUE
C
  560 CONTINUE
C
      END IF
C
C   NEXT FOR THE PITCHWISE GRID SPACINGS, FP(I) .
C   THROUGHFLOW
      IF(IM.EQ.2) THEN
           FP(1) = 1.0
           FP(2) = 0.0
      ELSE
C    END THROUGHFLOW
C
      IF(FP(NUM3).GT.0.00001) GO TO 563
      FPRAT = FP(1)
      FPMAX = FP(2)
      FP(1) = 1.0
      FP(IM)= 0.0
      DO 564 I=2,IMM1
      FP(I) = FP(I-1)*FPRAT
  564 CONTINUE
      DO 565 I=1,IMM1
      FREV  = FP(IM-I)
      IF(FREV.LT.FP(I))  FP(I) = FREV
      IF(FP(I).GT.FPMAX) FP(I) = FPMAX
  565 CONTINUE
C
  563 CONTINUE
C
      END IF
C
C***************************************************************************************
C***************************************************************************************
C   CALL "INPINT" TO INTERPOLATE IN THE SPANWISE VARIATION OF INFLOW PROPERTIES.
C   AND SET UP VALUES ON THE GRID POINTS.
C
      WRITE(6,*) ' CALLING INPINT TO SET UP THE INLET FLOW. '
           CALL INPINT
      WRITE(6,*) ' LEAVING INPINT.'
C
C*************************************************************************************** 
C******************************************************************************
C      READ IN THE MIXING LENGTH LIMITS IF ILOS IS NOT ZERO.
C
      DO 571 N = 1, NROWS
           XLLIM_I1(N)  = 0.03
           XLLIM_IM(N)  = 0.03
           XLLIM_K1(N)  = 0.03
           XLLIM_KM(N)  = 0.03
           XLLIM_DWN(N) = 0.03
           XLLIM_UP(N)  = 0.02
           XLLIM_IN(N)  = 0.02
           XLLIM_LE(N)  = 0.03
           XLLIM_TE(N)  = 0.04
           XLLIM_DN(N)  = 0.05
           FSTURB(N)    = 1.0
           TURBVIS_DAMP(N) = 0.5
  571 CONTINUE
C
C      READ IN THE MIXING LENGTH LIMITS IF ILOS IS NOT ZERO.
C
      IF(ILOS.NE.0) THEN
C
      READ(5,*) DUMMY_INPUT
C
      DO 570 N = 1,NROWS

      IF(ILOS.EQ.10) THEN
C
           READ(5,*,ERR=580)  XLLIM_I1(N),XLLIM_IM(N),
     &                 XLLIM_K1(N),XLLIM_KM(N),XLLIM_DWN(N),XLLIM_UP(N)
      END IF
C
      IF(ILOS.GE.100) THEN
C
           READ(5,*,ERR=580)  XLLIM_IN(N),XLLIM_LE(N),
     &                XLLIM_TE(N),XLLIM_DN(N),FSTURB(N),TURBVIS_DAMP(N)
      END IF
C           
  580 CONTINUE
C
           WRITE(6,*) ' ROW NUMBER ', N
           IF(ILOS.EQ.10) WRITE(6,*)
     &     'ILOS = 10 MIXING LENGTH LIMITS= ',
     &     XLLIM_I1(N),XLLIM_IM(N),XLLIM_K1(N),XLLIM_KM(N),XLLIM_DWN(N),
     &     XLLIM_UP(N)
           IF(ILOS.GE.100) WRITE(6,*) 'ILOS > 100 MIXING LENGTH LIMITS',
     &     XLLIM_IN(N),XLLIM_LE(N),XLLIM_TE(N),XLLIM_DN(N)
           IF(ILOS.GE.100) WRITE(6,*) 
     &     ' FREE STREAM TURBULENT VISCOSITY RATIO',FSTURB(N),
     &     ' MIXING PLANE TURBULENCE DECAY',TURBVIS_DAMP(N)
C
  570 CONTINUE
C
C     READ IN A FACTOR TO INCREASE THE TURBULENT VISCOSITY FOR THE FIRST NMIXUP STEPS.
C
           FACMIXUP = 2.0
           NMIXUP   = 1000
           READ(5,*) DUMMY_INPUT
           READ(5,*, ERR= 585)     FACMIXUP, NMIXUP
           IF(FACMIXUP.LT.1.0)     FACMIXUP = 1.0
           IF(IF_RESTART.NE.0)     FACMIXUP = 1.0
  585 CONTINUE
           WRITE(6,*) ' FACMIXUP = ', FACMIXUP,' NMIXUP = ',NMIXUP
C
      ENDIF
C
C******************************************************************************
C     WRITE OUT THE VALUE OF YPLUS AT THE WALL IF YPLUSWALL > 5.0
C
      IF(YPLUSWALL.GT.5.0)  THEN
           WRITE(6,*)
           WRITE(6,*)
     &    ' WALL SHEAR STRESSES CALCULATED USING THE YPLUSWALL MODEL. '
           WRITE(6,*) ' YPLUS AT THE WALL IS TAKEN AS, ', YPLUSWALL
           CFWALL = 1/(YPLUSWALL*YPLUSWALL)
      ENDIF
C******************************************************************************
C  READ IN THE SURFACE ROUGHNESSES IN MICRONS IF  IF_ROUGH  >  0  .
C
      IF(IF_ROUGH.GT.0) THEN
C
	WRITE(6,*) ' Non-hydraulic smooth surfaces specified'
	WRITE(6,*) ' Input Surface Roughness in microns for all 4'
	WRITE(6,*) ' surfaces of each row in turn.'
           READ(5,*) DUMMY_INPUT
      DO 586 N = 1,NROWS
           READ(5,*)  ROUGH_H(N),ROUGH_T(N),ROUGH_L(N),ROUGH_U(N)
C
           WRITE(6,*) ' ROW NUMBER ', N
           WRITE(6,*) ' Surface Roughnesses in microns ',
     &     ROUGH_H(N),ROUGH_T(N),ROUGH_L(N),ROUGH_U(N)
C
C   Change to physical roughness in metres . 
           ROUGH_H(N) = ROUGH_H(N)*1.0E-06
           ROUGH_T(N) = ROUGH_T(N)*1.0E-06
           ROUGH_L(N) = ROUGH_L(N)*1.0E-06
           ROUGH_U(N) = ROUGH_U(N)*1.0E-06
C
  586 CONTINUE
C
      ELSE
C
      DO 587 N = 1,NROWS
           ROUGH_H(N)  = 0.0
           ROUGH_T(N)  = 0.0
           ROUGH_L(N)  = 0.0
           ROUGH_U(N)  = 0.0
  587 CONTINUE
C
      ENDIF
C
C******************************************************************************
C******************************************************************************
C
      WRITE(6,*)
      WRITE(6,*)  ' Subroutine NEW_READIN completed, input data OK.'
      WRITE(6,*)
C
C
      RETURN
      END
C
C*************************************************************************
C
      SUBROUTINE OLD_READIN
C
C       THIS SUBROUTINE READS IN THE DATA IN
C       ====================

C
      INCLUDE 'commall-open-18.3'
C
      COMMON/BKODDS/
     &           XINT(JD),YINT(JD),RINT(JD),SDIST(JD),
     &           ICUSP(NRS),LCUSP(NRS),LCUSPUP(NRS),
     &           FRACNEW(JD),BETANEW(JD),SLOPE(JD),
     &           THICKUP(JD),THICKLOW(JD),
     &           XINT1(KD),XINT2(KD),XINT3(KD),XINT4(KD)
C
C       START OF INPUT SECTION
C       THROUGHOUT THE INPUT SECTION THE VARIABLES ARE AS FOLLOWS
C       =====================
C
C       RT_UPP(J,K)  IS TEMPORARY STORE FOR BLADE SUCTION SURFACE CO-ORD'S
C       RT_THICK(J,K) IS THE BLADE TANGENTIAL THICKNESS
C       RCYL(K)    IS RADIUS OF K'TH CYLINDRICAL SURFACE IF INPUT IS ON
C                  CYLINDRICAL SURFACES.
C       R(J,1)     IS RADII OF HUB
C       R(J,KM)    IS RADII OF CASING
C
C       =====================
 1000 FORMAT(16I5)
 1100 FORMAT(8F10.6)
 1200 FORMAT(A72)
 1600 FORMAT(16I5)
 1700 FORMAT(8F10.6)
 1701 FORMAT(8F12.6)
 1720 FORMAT(8F10.1)
 1730 FORMAT(F10.3,4F10.1,2F10.3)
 1800 FORMAT(18A4)
C
      PI     = 3.14159265
      DEGRAD = PI/180.
      RADDEG = 180./PI
C
      READ(5,1200)  TITLE
      WRITE(6,1800) TITLE
C
C************FIRST READ IN MAIN INTEGER CONTROL VARIABLES ************
C
      READ(5,1000) IM,JDUM,KM,IF_ROUGH,NSTEPS_MAX,IFCOOL,IFBLEED,NOSECT,
     &             NROWS,IFMIX,ISHIFT,KIN,NEXTRAP_LE,NEXTRAP_TE,NCHANGE
      WRITE(6,1600)IM,JDUM,KM,IF_ROUGH,NSTEPS_MAX,IFCOOL,IFBLEED,NOSECT,
     &             NROWS,IFMIX,ISHIFT,KIN,NEXTRAP_LE,NEXTRAP_TE,NCHANGE
C
      READ(5,1000) IN_VTAN,INSURF,IN_VR,ITIMST,IDUMMY,IPOUT,IN_FLOW,
     &             ILOS,NLOS,IF_RESTART,IDUMY,IBOUND,IF_REPEAT
      WRITE(6,1600)IN_VTAN,INSURF,IN_VR,ITIMST,IDUMMY,IPOUT,IN_FLOW,
     &             ILOS,NLOS,IF_RESTART,IDUMY,IBOUND,IF_REPEAT
C
      IN_PRESS = IN_VR
      
C
C******************************************************************************
C    SET THE COEFFICIENTS   F1, F2,  F3  FOR THE TIME STEPPING SCHEME
C
      WRITE(6,*)
      WRITE(6,*) ' TIME STEP TYPE , ITIMST = ', ITIMST
      WRITE(6,*)
C
C    SET THE COEFFICIENTS FOR THE SSS SCHEME
C
      IF(ITIMST.EQ.3.OR.ITIMST.EQ.5.OR.ITIMST.EQ.6) THEN
              F1          =  2.0000
              F2          = -1.000
              F3          =  0.00
              F2EFF       = -1.0
              NRSMTH      =  0
              RSMTH       =  0.40
      END IF
C
      IF(ITIMST.EQ.-3.OR.ITIMST.EQ.-5.OR.ITIMST.EQ.-6) THEN
              F1          =  2.0000
              F2          = -1.65
              F3          = -0.65
              F2EFF       = -1.0
              NRSMTH      =  1
              RSMTH       =  0.40
              ITIMST      =  ABS(ITIMST)
      END IF
C
      IF(ITIMST.EQ.4.OR.ITIMST.EQ.-4) THEN
              READ(5,*) F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*) ' F1, F2EFF, F3 , RSMTH, NRSMTH ',
     &                     F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*)
              IF(F2EFF.GT.0.0) THEN 
              WRITE(6,*)
              WRITE(6,*) ' ERROR,  F2EFF  MUST BE NEGATIVE.'
              WRITE(6,*) ' THE SIGN OF THE INPUT VALUE WILL BE CHANGED.'
              WRITE(6,*)
              F2EFF = -F2EFF
              END IF
              IF(F3.GT.0.0) THEN 
              WRITE(6,*)
              WRITE(6,*) ' ERROR,  F3  MUST BE NEGATIVE.'
              WRITE(6,*) ' THE SIGN OF THE INPUT VALUE WILL BE CHANGED.'
              WRITE(6,*)
              F3  = -F3
              END IF
              F2     = F2EFF*(1.0 - F3)
              ITIMST = 3
      END IF
C
              WRITE(6,*) ' F1, F2EFF, F3 , RSMTH, NRSMTH ',
     &                     F1, F2EFF, F3 , RSMTH, NRSMTH
              WRITE(6,*)
C
C**********************************************************************************
C      
C  Q3D
      IF(KM.EQ.2) IBOUND = 2
C  END Q3D
C
      READ (5,1000) IR,JR,KR,IRBB,JRBB,KRBB
      WRITE(6,1600) IR,JR,KR,IRBB,JRBB,KRBB
C
C  Q3D
      IF(KM.EQ.2) THEN
           KR   = 1
           KRBB = 1
      END IF
      IF(IM.EQ.2) THEN
           IR   = 1
           IRBB = 1
      END IF
C  END Q3D
C
      IF_KINT = 0
      IF(KIN.LT.0)     IF_KINT = 1
      KIN = ABS(KIN)      
      IF(KIN.EQ.0)     KIN  = KM
      IF(KIN.NE.KM)    IF_KINT = 1
C
      IF(NEXTRAP_LE.EQ.0)     NEXTRAP_LE  = 5
      IF(NEXTRAP_TE.EQ.0)     NEXTRAP_TE  = 5
      IF(NCHANGE.EQ.0) NCHANGE = NSTEPS_MAX/4
      IF(IF_RESTART.NE.0) NCHANGE = 100
      IF(NLOS.EQ.0)    NLOS = 5
      IMM1 = IM-1
      IMM2 = IM-2
      IF(IM.EQ.2) IMM2 =1
      KMM1 = KM-1
      KMM2 = KM-2
C
C
C     CHECK THAT THE DIMENSIONS ARE NOT TOO LARGE.
C
      IF(MAXKI.LT.KD)  WRITE(6,*)
     &      ' STOPPING BECAUSE MAXKI IS LESS THAN KD',
     &      ' MAXKI = ',MAXKI, ' KD = ',KD
      IF(MAXKI.LT.ID)  WRITE(6,*)
     &      ' STOPPING BECAUSE MAXKI IS LESS THAN ID',
     &      ' MAXKI = ',MAXKI, ' ID = ',ID
      IF(IM.GT.ID)  WRITE(6,*) ' STOPPING BECAUSE IM TOO LARGE.',
     &            ' IM= ',IM,  ' DIMENSION LIMIT = ',ID
      IF(KM.GT.KD)  WRITE(6,*) ' STOPPING BECAUSE KM TOO LARGE.',
     &            ' KM= ',KM,  ' DIMENSION LIMIT = ',KD
      IF(KIN.GT.MAXKI) WRITE(6,*) ' STOPPING BECAUSE KIN TOO LARGE.',
     &            ' KIN= ',KIN,' DIMENSION LIMIT = ',MAXKI
C
      IF(IM.GT.ID.OR.KM.GT.KD)  STOP
      IF(KIN.GT.MAXKI.OR.MAXKI.LT.KD.OR.MAXKI.LT.ID) STOP
C
      IF(NROWS.GT.NRS)  WRITE(6,*) 'STOPPING BECAUSE NROWS TOO LARGE.',
     &            ' NROWS= ',NROWS,' DIMENSION LIMIT = ',NRS
      IF(NROWS.GT.NRS) STOP
C
C     CHECK IF THERE ARE VARIABLE NUMBERS OF INPUT SECTIONS, i.e.  IF NOSECT IS NEGATIVE
C
      IF_SECTS = 0
      IF(NOSECT.LT.0) THEN
           IF_SECTS = 1
           NOSECT   = ABS(NOSECT)
      END IF
C
C
C
C***********INPUT THE BLADE GEOMETRY ON NOSECT STREAMWISE SURFACES******
C*********** THE GEOMETRY IS INPUT FOR EACH BLADE ROW SEPARATELY******
C
C            FIRST INPUT THE CONTROL PARAMETERS FOR THE BLADE ROW
C            AND SET THE VALUES OF SOME BLADE ROW VARIABLES
C
      J1        = 1
      IFSHROUD  = 0
C
      DO 1550 NR = 1,NROWS
C
C      READ IN THE ROW TITLE. THIS IS NEVER USED BUT HELPS TO LAY OUT THE DATA
C
      READ(5,1200) ROWTYP
C
C
C      KTIPS IS THE K VALUE OF THE POINT WHERE THE TIP GAP STARTS
C      KTIPE IS THE K VALUE OF THE POINT WHERE THE TIP GAP ENDS.
C      SET KTIPS(NR) = 0  FOR NO TIP GAP ON THIS ROW.    
C
C      SET KTIPS(NR) NEGATIVE TO USE THE SHROUD LEAKAGE MODEL ON THIS
C      ROW. EXTRA INPUT DATA IS THEN NEEDED AT THE END OF THE DATA FILE.
C
      READ(5,1000) JMROW,JLEROW,JTEROW,NBLADES_IN_ROW,
     &   KTIPS(NR),KTIPE(NR),JROTHS,JROTHE,JROTTS,JROTTE,
     &   NEW_GRID,JTRAN_I1(NR),JTRAN_IM(NR),JTRAN_K1(NR),
     &   JTRAN_KM(NR),IF_CUSP(NR)
C
C     FOR Q3D
      IF(KM.EQ.2) THEN
           KTIPS(NR) = 0
           KTIPE(NR) = 0
      END IF
C     END Q3D
C
C     INPUT THE NUMBER OF INPUT SECTIONS IF THIS IS NOT CONSTANT.
C
      IF(IF_SECTS.EQ.1) THEN
            READ(5,*) NSECS_ROW
      ELSE
            NSECS_ROW = NOSECT
      END IF
C    NSECS_IN  MUST BE SET AS IT IS USED IN INTPOL
            NSECS_IN = NSECS_ROW
C
C     END VARIABLE INPUT SECTIONS
C   
      JLEE = J1+JLEROW-1
      JTEE = J1+JTEROW-1
C
      IF(KTIPS(NR).LE.0) KTIPE(NR) = 0
      IF(KTIPS(NR).LT.0) IFSHROUD  = 1
C
      IF(JROTHS.EQ.0.OR.JROTHS.GT.JMROW) JROTHS = JMROW
      IF(JROTHE.EQ.0.OR.JROTHE.GT.JMROW) JROTHE = JMROW
      IF(JROTTS.EQ.0.OR.JROTTS.GT.JMROW) JROTTS = JMROW
      IF(JROTTE.EQ.0.OR.JROTTE.GT.JMROW) JROTTE = JMROW
C
C
      WRITE(6,1600) JMROW,JLEROW,JTEROW,NBLADES_IN_ROW,
     &   KTIPS(NR),KTIPE(NR),JROTHS,JROTHE,JROTTS,JROTTE,NEW_GRID,
     &   JTRAN_I1(NR),JTRAN_IM(NR),JTRAN_K1(NR),JTRAN_KM(NR),IF_CUSP(NR)
C
C
C     Set the cusp generation parameters if IF_CUSP is not 0.
C     Maintain the original grid and no cusp is generated if IF_CUSP = 0
C     The cusp is centred on the blade centre line if ICUSP = 0.
C     The cusp makes the I=1 surface continuous on the cusp if  ICUSP =  1.
C     The cusp makes the I=IM surface continuous on the cusp if ICUSP = -1.
C     The cusp is of length  LCUSP and starts LCUSPUP points before the
C     trailing edge.
C
C
      IF(IF_CUSP(NR).EQ.0) WRITE(6,*) '  NO CUSP SPECIFIED '
C
      IF(IF_CUSP(NR).EQ.1) THEN
	   READ(5,*)    ICUSP(NR),LCUSP(NR),LCUSPUP(NR)
           WRITE(6,*) ' ICUSP, LCUSP , LCUSPUP = ',
     &                  ICUSP(NR),LCUSP(NR),LCUSPUP(NR)
           WRITE(6,*)
      END IF
C
C    IF  IF_CUSP = 2 a body force is used to force separation at the trailing edge.
C    The force starts NUP_I1 grid points upstream of the TE  om the I=1 blade surface, 
C    and at NUP_IM points upstream on the I=IM blade surface. It extends N_WAKE points
C    downstream of the TE. The thickness of the affected layer is determined bt SEP_THIK
C    , typical value  0.01, and the strength of the body force by SEP_DRAG, typical value 0.99.
C
      IF(IF_CUSP(NR).EQ.2) THEN
           ICUSP(NR)   = 0
           LCUSP(NR)   = 0
           LCUSPUP(NR) = 0
           READ(5,*)  NUP_I1(NR),NUP_IM(NR),N_WAKE(NR),
     &                  SEP_THIK(NR),SEP_DRAG(NR)
           WRITE(6,*)' BODY FORCE USED TO FORCE TRAILING EDGE SEPARTION'
           WRITE(6,*)' CUSP BODY FORCE POINTS, NUP_I1, NUP_IM, N_WAKE, 
     & SEP_THICK, SEP_DRAG'
             WRITE(6,*) NUP_I1(NR),NUP_IM(NR),N_WAKE(NR),
     &                  SEP_THIK(NR),SEP_DRAG(NR)
             WRITE(6,*)
      END IF
C
C
C     READ IN THE RPM, TIP CLEARANCE AND INITIAL GUESS OF INLET AND EXIT 
C     PRESSURES FOR THIS BLADE ROW.
C
      READ(5,1100) RPMROW,PUPROW,PLEROW,PTEROW,PDNROW,FRACTIP(NR),RPMHUB
      WRITE(6,1730)RPMROW,PUPROW,PLEROW,PTEROW,PDNROW,FRACTIP(NR),RPMHUB
C
      IF(FRACTIP(NR).LT.0.0) THEN 
           READ(5,*) FRACTIP1(NR),FRACTIP2(NR)
           WRITE(6,*)' FRACTIP1,  FRACTIP2 = ',FRACTIP1(NR),FRACTIP2(NR)
      ELSE
           FRACTIP1(NR) = FRACTIP(NR)
           FRACTIP2(NR) = FRACTIP(NR)
      END IF
C
C
C      FTHICK(NR,K) IS THE MULTIPLYING FACTOR ON THE BLADE THICKNESS SO
C      THAT IT CAN BE REDUCED AT THE TIP. IT IS ASSUMED TO BE 1.0 UNLESS
C      INPUT HERE.
C
      DO 2507 K=1,KM
 2507 FTHICK(NR,K)=1.0
      IF(KTIPS(NR).GT.0)  READ(5,1100) (FTHICK(NR,K),K=1,KM)
      IF(KTIPS(NR).GT.0) WRITE(6,1700) (FTHICK(NR,K),K=1,KM)
C
C*********************************************************************************
C*********************************************************************************
C*********************************************************************************
C
      J2 = J1 + JMROW - 1
      ANGLEAN = 0.0
      IF_ANGLES(NR) = 0
C
C*********************************************************************************
C*********************************************************************************
C      NOW  READ IN MAIN GEOMETRICAL DATA FOR THE BLADE ROW 
C      ON  NOSECT  BLADE SECTIONS OF THE CURRENT BLADE ROW
C*********************************************************************************
C**********************************************************************************
C
      DO 1555 K=1,NSECS_ROW
C
C
      READ(5,1111)  FAC1,XSHIFT,IF_DESIGN,IF_RESTAGGER,IF_LEAN
      WRITE(6,1111) FAC1,XSHIFT,IF_DESIGN,IF_RESTAGGER,IF_LEAN
 1111 FORMAT(2F10.5,3I10)
C 
C******************************************************************************
C******************************************************************************
C    DO NOT READ IN AN EXISTING BLADE SECTION IF  "IF_DESIGN"  IS NON-ZERO. 
C    JUMP TO 1510 TO DESIGN A NEW SECTION
C 
      IF(IF_DESIGN.NE.0) GO TO 1510
C
C******************************************************************************
C******************************************************************************
C***********XSURF(J,K) IS AXIAL COORDINATE OF POINTS ON THE STREAMWISE SURFACE.
C
      READ(5,1100)  (XSURF(J,K),J=J1,J2)
      WRITE(6,1700) (XSURF(J,K),J=J1,J2)
C
C
C     LEAN THE WHOLE BLADE BY AN ANGLE ANGLEAN IF ANGLEAN IS GREATER THAN ZERO)
C
           READ(5,1100)  FAC2,TSHIFT
           WRITE(6,1700) FAC2,TSHIFT   
C
C***********RT_UPP(J,K) IS THE R-THETA COORDINATE OF POINTS ON THE STREAMWISE
C           SURFACE ON THE BLADE SURFACE WITH LARGEST VALUE OF THETA.
C           i.e. THE UPPER SURFACE OF THE BLADE AND THE LOWER SURFACE OF THE
C                BLADE TO BLADE PASSAGE.
C
      READ(5,1100)  (RT_UPP(J,K),J=J1,J2)
      WRITE(6,1700) (RT_UPP(J,K),J=J1,J2)
C
C
      READ(5,1100)  FAC3,BETAUP(NR,K),BETADWN1(NR,K),BETADWN2(NR,K)
      WRITE(6,1700) FAC3,BETAUP(NR,K),BETADWN1(NR,K),BETADWN2(NR,K)
C
C
      IF(ABS(BETADWN2(NR,K)).LT.0.0001) BETADWN2(NR,K) = BETADWN1(NR,K)
      IF(ABS(BETADWN2(NR,K)).GT.0.0001) IF_ANGLES(NR) = 1
C
C
C***********RT_THICK(J,K) IS BLADE THICKNESS DELTA R-THETA MEASURED
C           IN THE TANGENTIAL DIRECTION.
C
      READ(5,1100)  (RT_THICK(J,K),J=J1,J2)
      WRITE(6,1700) (RT_THICK(J,K),J=J1,J2)
C
C
C      IF INSURF =1 OR 2 RSURF(J,K) IS RADII OF POINTS ON THE STREAMWISE
C      SURFACES ON WHICH DATA IS INPUT.
C      IF INSURF =0 RCYL(J) IS THE RADIUS OF THE CYLINDRICAL SURFACES ON
C       DATA IS INPUT
C
      READ(5,1100)   FAC4,RSHIFT
      WRITE(6,1700)  FAC4,RSHIFT
C
C
      IF(INSURF.NE.0) READ(5,1100)  (RSURF(J,K),J=J1,J2)
      IF(INSURF.NE.0) WRITE(6,1700) (RSURF(J,K),J=J1,J2)
      IF(INSURF.EQ.0) READ(5,1100)  RCYL(K)
      IF(INSURF.EQ.0) WRITE(6,1700) RCYL(K)
C
C******************************************************************************
C******************************************************************************
C     SHIFT THE BLADE SECTION BY XSHIFT OR TSHIFT
C
      DO 1500 J=J1,J2
      IF(INSURF.EQ.0)  RSURF(J,K) = RCYL(K)
      XSURF(J,K)    = FAC1*(XSHIFT + XSURF(J,K))
      RT_THICK(J,K) = FAC3*RT_THICK(J,K)
      RSURF(J,K)    = FAC4*(RSHIFT + RSURF(J,K))
      RT_UPP(J,K)   = FAC2*(TSHIFT + RT_UPP(J,K))
 1500 CONTINUE
C 
C******************************************************************************
C******************************************************************************
 1510 CONTINUE
C******************************************************************************
C******************************************************************************
C     JDD ADDITION TO ALLOW REDESIGN OF THIS BLADE SECTION.
C     CHANGED TO USE SUBROUTINE RE_DESIGN. AUGUST 2016.
C
      IF(IF_DESIGN.NE.0) THEN
          WRITE(6,*)'CALLING RE_DESIGN TO GENERATE A NEW BLADE GEOMETRY'
          CALL RE_DESIGN(NR,K,J1,J2,JLEROW,JTEROW) 
          BETAUP(NR,K)   = 0.0
          BETADWN1(NR,K) = 0.0
          BETADWN2(NR,K) = 0.0
          IF_ANGLES(NR)  = 0
      END IF
C   
C    END OF OPTION TO REDESIGN THE BLADE SECTION
C******************************************************************************
C******************************************************************************
C
C     RESTAGGER THE BLADE SECTION IF "IF_RESTAGGER" IS NON-ZERO.
C
      IF(IF_RESTAGGER.NE.0) THEN
      WRITE(6,*) ' CALLING SUBROUTINE RESTAGGER TO ROTATE THE BLADE '
      CALL RESTAGGER(K,J1,J2,JLEROW,JTEROW,ROTATE,FRACX_ROT)
      END IF
C
C******************************************************************************
C******************************************************************************
C 
C    LEAN THE BLADE SECTION IF "IF_LEAN"  IS NOT ZERO.
C
      IF(IF_LEAN.NE.0) THEN
      WRITE(6,*) ' CALLING SUBROUTINE LEAN TO LEAN THE BLADE SECTION '
      CALL LEAN(K,J1,J2,ANGLEAN)
      END IF
C
C******************************************************************************
C******************************************************************************
C   Q3D
C     CALL SET SSTHICK TO SET THE STREAM SURFACE THICKNESS IF DOING A Q3D
C     BLADE TO BLADE CALCULATION. THEN END THE GEOMETRY INPUT FOR THIS BLADE ROW
C     AS ONLY A SINGLE STREAM SURFACE IS USED .
C
      IF(KM.EQ.2) THEN
           CALL SET_SSTHICK(J1,J2)
C    JUMP OUT OF THE DO 1555 LOOP AS ONLY ONE STREAM SURFACE IS NEEDED.
C    SKIP THE SEPARATE HUB AND CASING GEOMETRY INPUT IF DOING A Q3D CALCULATION.
           GO TO 1504
      END IF
C   END Q3D
C
C******************************************************************************
C******************************************************************************
C     END IF INPUT FOR THIS BLADE SECTION
C
 1555 CONTINUE
C
C
C******************************************************************************
C******************************************************************************
C  TFLOW
C******************************************************************************
C******************************************************************************
C  IF A THROUGHFLOW CALCULATION READ IN THE REQUIRED EXIT FLOW ANGLE OR DEVIATION 
C  ANGLE FOR THIS BLADE ROW. ANGL_TYP = 'A' FOR EXIT ANGLE, = 'D' FOR DEVIATION ANGLE.
C
      IF(IM.EQ.2) THEN  
      READ(5,*)    DUMMY_INPUT
      WRITE(6,*) ' INPUTTING THE THROUGHFLOW DATA '
      READ(5,*)    ANGL_TYP(NR), NANGLES(NR)
      READ(5,*)  ( FRAC_SPAN(NR,NK),NK=1,NANGLES(NR))
      READ(5,*)  ( EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
      WRITE(6,*) ' NANGLES = ',NANGLES(NR),' ANGL_TYP = ',ANGL_TYP(NR)
      WRITE(6,9048) (FRAC_SPAN(NR,NK),NK=1,NANGLES(NR))
      IF(ANGL_TYP(NR).EQ.'A')
     &     WRITE(6,9045) (EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
      IF(ANGL_TYP(NR).EQ.'D')
     &     WRITE(6,9046) (EXIT_ANGL(NR,NK),NK=1,NANGLES(NR))
 9048 FORMAT(' FRACTION OF SPAN',/,(10F10.3))
 9045 FORMAT(' EXIT FLOW ANGLE IN DEGREES ',/,(10F10.3))
 9046 FORMAT(' EXIT FLOW ANGLE DEVIATION FROM GRID ANGLE, IN DEGREES.',
     & /,(10F10.3))
      WRITE(6,*) ' END THROUGHFLOW DATA '
      END IF
C ******************************************************************************
C  END TFLOW
C******************************************************************************
C******************************************************************************
C
c          INPUT THE HUB AND CASING GEOMETRY.
C          IF INSURF = 2 THE HUB AND CASING  ARE TAKEN AS THE
C          FIRST AND LAST STREAMWISE SURFACES. OTHERWISE READ IN
C          HUB AND CASING COORDINATES AT THE ENDS OF THE QUASI ORTHOGONALS.
C          THIS CANNOT BE USED IF KM = 2 for a Q3D CALCULATION.
C
      IF(INSURF.NE.2) THEN
C
      READ(5,1100)  (X(J,1),J=J1,J2)
      READ(5,1100)  (R(J,1),J=J1,J2)
      WRITE(6,1700) (X(J,1),J=J1,J2)
      WRITE(6,1700) (R(J,1),J=J1,J2)
      READ(5,1100)  (X(J,KM),J=J1,J2)
      READ(5,1100)  (R(J,KM),J=J1,J2)
      WRITE(6,1700) (X(J,KM),J=J1,J2)
      WRITE(6,1700) (R(J,KM),J=J1,J2)
      DO 1505 J=J1,J2
      X(J,1)  = (XSHIFT + X(J,1))*FAC1
      R(J,1)  = (RSHIFT + R(J,1))*FAC4
      X(J,KM) = (XSHIFT + X(J,KM))*FAC1
 1505 R(J,KM) = (RSHIFT + R(J,KM))*FAC4
C
      ELSE
C
      DO 1503 J = J1,J2
      R(J,1)    = RSURF(J,1)
      R(J,KM)   = RSURF(J,NSECS_ROW)
      X(J,1)    = XSURF(J,1)
      X(J,KM)   = XSURF(J,NSECS_ROW)
 1503 CONTINUE
C
      END IF
C
C*********************************************************************************
C*********************************************************************************
c     interpolate extra input sections if  NSECS_ROW  NOT EQUAL TO   NOSECT
C     THIS ROUTINE PROVIDED BY S GALLIMORE.
c
      IF(NSECS_ROW.NE.NOSECT)THEN
        DO J = J1,J2
	  FSPAN(1) = 0.0
          DO K = 1,NSECS_ROW
            XINT1(K) = xsurf(J,K)
            XINT2(K) = rsurf(J,K)
            XINT3(K) = rt_upp(J,K)
            XINT4(K) = rt_thick(J,K)
            IF (K.GT.1) THEN
              RD = rsurf(J,K)-rsurf(J,K-1)
              XD = xsurf(J,K)-xsurf(J,K-1)
	    FSPAN(K) = FSPAN(K-1) + SQRT(XD*XD+RD*RD)
            END IF
          END DO
c
          xarg = 0.0
          DO K = 1,NOSECT
            IF(k.gt.1) xarg = xarg+(FSPAN(NSECS_ROW)/(nosect-1))
            CALL INTP(NSECS_ROW,FSPAN,XINT1,XARG,xsurf(J,K))
            CALL INTP(NSECS_ROW,FSPAN,XINT2,XARG,rsurf(J,K))
            CALL INTP(NSECS_ROW,FSPAN,XINT3,XARG,rt_upp(J,K))
            CALL INTP(NSECS_ROW,FSPAN,XINT4,XARG,rt_thick(J,K))
          END DO
        END DO
      END IF
C
C****************************************************************************
C  RE-ENTER HERE IF KM = 2
 1504 CONTINUE
C
C****************************************************************************
C      CALL NEWGRID FOR THE CURRENT BLADE ROW IF REQUESTED
C      THIS READS IN MORE DATA AND GENERATES NEW STREAMWISE (J) GRID POINTS.
C
      IF(NEW_GRID.NE.0) CALL NEWGRID(JMROW,J1,J2,
     &                JLEROW,JTEROW,JROTHS,JROTHE,JROTTS,JROTTE)
C
C*****************************************************************************
C
C      SET THE J MARKERS, INITIAL GUESS OF PRESSURE, ETC
C
      JLE(NR) = J1+JLEROW-1
      JTE(NR) = J1+JTEROW-1
      J2      = J1 -1 + JMROW
      JROW    = 0
C
      DO 100 J = J1,J2
      JM1      = J-1
      IF(J.EQ.1) JM1=1
      JROW     = JROW+1
C
C      SET THE ROTATION OF THE HUB AND CASING.
C      THE DEFAULT IF JROTHS AND JROTTS ARE ZERO
C      IS THAT THE WHOLE HUB IS ROTATING AT THE SAME SPEED AS
C      THE BLADE ROW AND THE CASING IS STATIONARY.
C
      WRAD(J) = RPMROW*3.14159/30.
      WHUB(J) = 0.0
      IF(JROW.GT.JROTHS.AND.JROW.LT.JROTHE) WHUB(J) = RPMHUB*3.14159/30.
      WTIP(J) = 0.0
      IF(JROW.GT.JROTTS.AND.JROW.LT.JROTTE) WTIP(J) = WRAD(J)
C
C     SET LIMITS OF SOLID PART OF THE BLADE ELEMENTS, IE WHERE NO FLOW.
C     THESE ARE FROM 1 TO KTIPS-1  OR FROM KTIPE TO KM-1
C
      IF(KTIPS(NR).LE.0) THEN
      KS1(NR) = 1
      KS2(NR) = KMM1
      ENDIF
      IF(KTIPS(NR).EQ.1) THEN
      KS1(NR) = KTIPE(NR)
      KS2(NR) = KMM1
      ENDIF
      IF(KTIPS(NR).GT.1) THEN
      KS1(NR)   = 1
      KS2(NR)   = KTIPS(NR)-1
      KTIPE(NR) = KM
      ENDIF
C
      IND(J)    = 0
      INDLE(J)  = 0
      INDTE(J)  = 0
      INDMIX(J) = 0
      INDMID(J) = 0
      NBLADE(J) = NBLADES_IN_ROW
      NROW(J)   = NR
      JMIDRW    = 0.5*(JLEROW+JTEROW)
      IF(JROW.GT.JLEROW.AND.JROW.LE.JTEROW) IND(J)=1
      IF(JROW.EQ.JLEROW) INDLE(J)    = 1
      IF(JROW.EQ.JLEROW) JLED(NR)    = J
      IF(JROW.EQ.JTEROW) INDTE(J)    = 1
      IF(JROW.EQ.JLEROW) INDLE(J-1) = -1
      IF(JROW.EQ.JLEROW) INDLE(J-2) = -2
      IF(JROW.EQ.JMIDRW) INDMID(J)   = 1
      IF(JROW.LE.JLEROW) PGUESS(J) = PUPROW+(PLEROW-PUPROW)*
     & (JROW-1)/(JLEROW-1)
      IF(JROW.GT.JLEROW.AND.JROW.LE.JTEROW) PGUESS(J) = 
     &          PLEROW + (PTEROW-PLEROW)*(JROW-JLEROW)/(JTEROW-JLEROW)
      IF(JROW.GT.JTEROW) PGUESS(J) =
     &          PTEROW + (PDNROW-PTEROW)*(JROW-JTEROW)/(JMROW-JTEROW)
  100 CONTINUE
C
       INDMIX(J2) = 1
       JMIX(NR)   = J2
       JSTART(NR) = J1
C
C     RESET THE J INDEX TO CONTINUE THROUGH THE NEXT BLADE ROW, IF ANY.
C
      J1 = J2 + 1
C
C        END OF DATA INPUT ON THIS BLADE ROW. RETURN TO INPUT DATA ON THE
C        NEXT ROW UNLESS THIS IS THE LAST, IE UNLESS NR=NROWS.
C
 1550 CONTINUE
C
C       END OF INPUT OF BLADE GEOMETRY DATA
C
      JM         = J2
      JMM1       = JM-1
      JMM2       = JM-2
      INDMIX(JM) = 0
C
C
C     CHECK THAT JM IS NOT TOO LARGE.
C
      IF(JM.GT.JD)  WRITE(6,*) 'STOPPING BECAUSE JM TOO LARGE.',
     &            ' JM= ',JM,' DIMENSION LIMIT = ',JD
      IF(JM.GT.JD) STOP
C
C**********************************************************************************
C**********************************************************************************
C    START TO SET THE GRID BETWEEN BLADE ROWS AND MAKE THE MIXING PLANES CONTIGUOUS
C**********************************************************************************
C
C     SKIP THE NEXT PART AND DO NOT MOVE THE BLADES OR GRID IF ISHIFT = 0 .
      IF(ISHIFT.EQ.0) GO TO 140
C     
      IF(ISHIFT.GE.2) GO TO 133
C
C      NEXT SECTION ONLY FOR ISHIFT = 1
C      IF  ISHIFT = 1  AUTOMATICALLY SHIFT THE BLADES TO MAKE THE X SPACING CONTINUOUS
C      BETWEEN THE LAST GRID POINT ON ONE BLADE AND THE FIRST GRID POINT ON THE NEXT
C      ONE ON THE HUB STREAMLINE BUT DO NOT CHANGE THE GRID SPACINGS.
C
      XSHIFT = 0.0
      DO 130 J = 2,JM
      IF(INDMIX(J-1).NE.1) GO TO 129
      DIFF1  = XSURF(J-1,1) - XSURF(J-2,1)
      DIFF2  = XSURF(J+1,1) - XSURF(J,1)
      XSHIFT = XSURF(J-1,1) + 0.01*(DIFF1+DIFF2) - XSURF(J,1)
  129 CONTINUE
      DO 131 K=1,KM
      XSURF(J,K) = XSURF(J,K) + XSHIFT
      X(J,K)     = X(J,K)     + XSHIFT
  131 CONTINUE
  130 CONTINUE
C
C     MAKE NO MORE CHANGES TO THE GRID IF ISHIFT = 1.
      IF(ISHIFT.EQ.1) GO TO 140
C
C**********************************************************************************
C     RE ENTER HERE IF ISHIFT = 2
  133 CONTINUE
C
C********************************************************************************
C      START TO MAKE THE GRID SPACING VARY GEOMETRICALLY BETWEEN BLADE ROWS 
C      IF ISHIFT >= 2.
C*******************************************************************************
C
      NR   = 1
      JST  = 1
      JEND = JLE(1)
C
C******************************************************************************
C
 6666 CONTINUE
C
C******************************************************************************
      DO 210 K=1,NOSECT
C
      SMERID(1,K) = 0.0
      DO 213 J=2,JM
      XD = XSURF(J,K)  - XSURF(J-1,K)
      RD = RSURF(J,K)  - RSURF(J-1,K)
  213 SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)
C
C    IF ISHIFT = 2 MAINTAIN THE MERIDIONAL CURVE IN THE GAP BETWEEN BLADE ROWS.
C    IF ISHIFT = 3 MAKE THE MERIDIONAL VIEW OF THE GRID LINEAR IN THE GAP BETWEEN BLADE ROWS.
C    IF ISHIFT = 4 SAME AS 3 BUT DO NOT CHANGE THE HUB AND CASING PROFILES
C
      DLAST = -1.0
      NP = 1
      XGAP = XSURF(JEND,K)  - XSURF(JST,K)
      RGAP = RSURF(JEND,K)  - RSURF(JST,K)
      GAP  = SQRT(XGAP*XGAP + RGAP*RGAP)
C
      DO 110 J = JST,JEND
      XD   = XSURF(J,K)  - XSURF(JST,K)
      RD   = RSURF(J,K)  - RSURF(JST,K)
      PROJ = (XD*XGAP + RD*RGAP)

      IF(ISHIFT.GE.3) PROJ = PROJ/GAP

      DIST = SQRT(XD*XD+RD*RD)
      IF(PROJ.LT.0.0) GO TO 110
      IF(DIST.GT.(1.00001*GAP)) THEN
      XINT(NP) = XSURF(JEND,K)
      RINT(NP) = RSURF(JEND,K)
      NP = NP+1
      GO TO 111
      ENDIF
      IF(DIST.LT.DLAST) GO TO 110
      DLAST     = DIST

      IF(ISHIFT.EQ.3) THEN
           XINT(NP)  = XSURF(JST,K)   + PROJ*XGAP/GAP
           RINT(NP)  = RSURF(JST,K)   + PROJ*RGAP/GAP
      ENDIF

      IF(ISHIFT.EQ.2) THEN
           XINT(NP)  = XSURF(J,K)
           RINT(NP)  = RSURF(J,K)
      ENDIF

      IF(ISHIFT.EQ.4) THEN
           IF(K.EQ.1.OR.K.EQ.NOSECT) THEN
                XINT(NP)  = XSURF(J,K)
                RINT(NP)  = RSURF(J,K)
           ELSE
                XINT(NP)  = XSURF(JST,K)   + PROJ*XGAP/GAP
                RINT(NP)  = RSURF(JST,K)   + PROJ*RGAP/GAP
           END IF
      END IF
C
      NP = NP+1
C
  110 CONTINUE
C
  111 CONTINUE
C
C    SMOOTH THE X AND R COORDINATES IN THE BLADE TO BLADE GAP
C
      NGAP    = NP - 1
      NSMOOTH = 10
      SFGAP   = 0.2
      DO 113 NS = 1,NSMOOTH
      DO 114 NN = 2,NGAP-1
      XINT(NN) = (1.-SFGAP)*XINT(NN) + SFGAP*0.5*(XINT(NN-1)+XINT(NN+1))
      RINT(NN) = (1.-SFGAP)*RINT(NN) + SFGAP*0.5*(RINT(NN-1)+RINT(NN+1))
  114 CONTINUE
  113 CONTINUE
C
C
      SDIST(1)  = SMERID(JST,K)
      DO 112 NN = 2,NGAP
      XD = XINT(NN) - XINT(NN-1)
      RD = RINT(NN) - RINT(NN-1)
  112 SDIST(NN) = SDIST(NN-1) + SQRT(XD*XD+RD*RD)
C
C
      SMID   = 0.5*(SDIST(1) + SDIST(NGAP))
      STE    = SDIST(1)
      SLE    = SDIST(NGAP)
      IF(NR.LE.NROWS) ANGLUP = BETAUP(NR,K)
C
      IF(NR.NE.1)  THEN
                   ANGLDWN1 = BETADWN1(NR-1,K)
                   ANGLDWN2 = BETADWN2(NR-1,K)
      END IF
C
C******************************************************************************
C    START TO CALL GRID_UP AND GRID_DOWN TO SET THE GRID UPSTREAM AND DOWNSTREAM OF
C    ALL BLADE ROWS.
C******************************************************************************
C    CALL GRID_UP TO FORM THE GRID UPSTREAM OF THE FIRST BLADE ROW.
C
      IF(NR.EQ.1.AND.K.EQ.1) THEN
            WRITE(6,*)
            WRITE(6,*) ' STARTING NEW BLADE ROW, ROW NUMBER', 1
            WRITE(6,*)
      END IF
      IF(NR.EQ.1) WRITE(6,*)
     & ' CALLING GRID_UP FOR ROW 1, K = ',K, ' JSTART,JLE1= ',JST,JEND
C
      IF(NR.EQ.1)  CALL GRID_UP(K,JST,JEND,STE,SLE,NGAP,SDIST,XINT,
     &             RINT,NEXTRAP_LE,ANGLUP,IF_ANGLES(NR))
C
C     CALL GRID_UP AND GRID_DOWN TO FORM THE GRIDS UPSTREAM AND DOWNSTREAM OF THE
C     MIXING PLANE FOR INTERIOR BLADE ROWS.
C
      IF( (NR.NE.1.AND.NR.NE.NROWS+1).AND.K.EQ.1) THEN
            WRITE(6,*)
            WRITE(6,*) ' STARTING NEW BLADE ROW, ROW NUMBER',NR
            WRITE(6,*)
      END IF
      IF(NR.NE.1.AND.NR.NE.NROWS+1) THEN
C
      WRITE(6,*)
     &' CALLING GRID_DOWN FOR ROW No.',NR,'K = ',K,' JTE,JMID= ',
     &  JST,JMID
C
      CALL GRID_DOWN(K,JST,JMID,STE,SMID,NGAP,SDIST,XINT,RINT,
     &         NEXTRAP_TE,ANGLDWN1,ANGLDWN2,IF_CUSP(NR-1),ICUSP(NR-1),
     &         LCUSP(NR-1),LCUSPUP(NR-1),IF_ANGLES(NR-1) )
      WRITE(6,*)
     &' CALLING GRID_UP FOR ROW No.  ',NR,'K = ',K,'JMID, JLE= ',
     &  JMID+1,JEND
      CALL GRID_UP(K,JMID+1,JEND,SMID,SLE,NGAP,SDIST,XINT,RINT,
     &             NEXTRAP_LE,ANGLUP,IF_ANGLES(NR) )
C
      ENDIF
C
C     CALL GRID_DOWN TO FORM THE GRID DOWNSREAM OF THE LAST BLADE ROW.
C
      IF(NR.EQ.NROWS+1.AND.K.EQ.1) THEN
            WRITE(6,*)
            WRITE(6,*) ' STARTING NEW BLADE ROW, ROW NUMBER',NROWS
            WRITE(6,*)
      END IF
C
      IF(NR.EQ.NROWS+1)  THEN
           ANGLDWN1 = BETADWN1(NROWS,K)
           ANGLDWN2 = BETADWN2(NROWS,K)
      WRITE(6,*)
     & ' CALLING GRID_DOWN FOR THE LAST BLADE ROW, K = ',K,'JTE,JEND= ',
     &   JST,JEND
C
      CALL GRID_DOWN(K,JST,JEND,STE,SLE,NGAP,SDIST,
     &        XINT,RINT,NEXTRAP_TE,ANGLDWN1,ANGLDWN2,IF_CUSP(NR),
     &        ICUSP(NROWS),LCUSP(NROWS),LCUSPUP(NROWS),IF_ANGLES(NROWS))
C
      END IF
C
      IF(ABS(ANGLUP).GT.0.0.OR.ABS(ANGLDWN1).GT.0.0) THEN
           WRITE(6,*)
           WRITE(6,*) 'ROW No ',NR,'SECTION No ',K,'ANGLUP= ',ANGLUP,
     &                'ANGLDWN1 & 2= ',ANGLDWN1,ANGLDWN2,' DEGREES.'
           WRITE(6,*)
      END IF
C
C     END OF SETTING THE GRID UPSTREAM AND DOWNSREAM OF ALL BLADE ROWS.
C
  210 CONTINUE
C
      WRITE(6,*)
C
C********************************************************************************
C
      IF(JEND.EQ.JM) GO TO 200
C
      NR = NR + 1
      IF(NR.NE.NROWS+1) THEN
      JST    = JTE(NR-1)
      JEND   = JLE(NR)
      JMID   = JMIX(NR-1)
      GO TO 6666
      ENDIF
      IF(NR.EQ.NROWS+1) THEN
      JST  = JTE(NROWS)
      JEND = JM
      GO TO 6666
      ENDIF
C
  200 CONTINUE
C
C      RESET SMERID
C
      DO 216 K=1,NOSECT
      DO 216 J=2,JM
      XD = XSURF(J,K) - XSURF(J-1,K)
      RD = RSURF(J,K) - RSURF(J-1,K)
  216 SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)

C      RESET THE HUB AND CASING COORDINATES TO NEW STREAM SURFACE VALUES IF INSURF = 2.
C
      IF(INSURF.EQ.2) THEN
C
       DO 214 J=1,JM
       X(J,1)  = XSURF(J,1)
       R(J,1)  = RSURF(J,1)
       X(J,KM) = XSURF(J,NOSECT)
       R(J,KM) = RSURF(J,NOSECT)
  214 CONTINUE
C
      ELSE
C
C    RESET THE HUB AND CASING COORDINATES IF INSURF IS NOT = 2.
C
      DO 218       NR = 1,NROWS+1
      IF(NR.NE.1)  JT = JTE(NR-1)
      IF(NR.EQ.1)  JT = 1
      IF(NR.NE.NROWS+1) JL = JLE(NR)
      IF(NR.EQ.NROWS+1) JL = JM
      DO 215 K = 1,KM,KMM1
      M = 1
      IF(K.EQ.KM) M = NOSECT
      DO 333      J = JT,JL
      NJ            = J - JT + 1
      IF(NJ.EQ.1) THEN
      SDIST(NJ) = 0.0
      ELSE
      XD = X(J,K) - X(J-1,K)
      RD = R(J,K) - R(J-1,K)
      SDIST(NJ) = SDIST(NJ-1) + SQRT(XD*XD+RD*RD)
      ENDIF
      XINT(NJ) = X(J,K)
      RINT(NJ) = R(J,K)
  333 CONTINUE
C
      GAP = SDIST(NJ) - SDIST(1)
      DO 334 NN = 2,NJ
  334 SDIST(NN) = (SDIST(NN)-SDIST(1))/GAP
C
      SGAP = SMERID(JL,M) - SMERID(JT,M)
      DO 217 J = JT,JL
      ARG = (SMERID(J,M) - SMERID(JT,M))/SGAP
      CALL LININT(NJ,SDIST,XINT,ARG,X(J,K))
      CALL LININT(NJ,SDIST,RINT,ARG,R(J,K))
  217 CONTINUE
  215 CONTINUE
  218 CONTINUE
C
      ENDIF
C
  140 CONTINUE
C
C   END OF GEOMETRY INPUT AND MANIPULATION.
C******************************************************************************
C******************************************************************************
C
C
C******************************************************************************
C******************************************************************************
C**********READ IN GAS CONSTANTS ,TIME STEP LENGTH,SMOOTHING FACTOR,ETC.
C
      WRITE(6,*)'READING IN THE GAS CONSTANTS,TIME STEP LENGTH,SMOOTHING 
     &FACTOR,ETC. '
C
      IFGAS = 0
      READ(5,1100)  CP,GA,CFL,SFTIN,SFXIN,MACHLIM
      WRITE(6,1701) CP,GA,CFL,SFTIN,SFXIN,MACHLIM
      IF(MACHLIM.LT.1.0) MACHLIM = 2.0
C
C    USE REAL GAS PROPERTIES IF CP IS INPUT AS NEGATIVE.
C    TYPICAL VALUES FOR COMBUSTION PRODUCTS  ARE:CP1 = 1272.5, CP2 = 0.2125, CP3 = 0.000015625, RGAS = 287.15
C    AT  TREF = 1400 K.
C
      IF(CP.LT.0.0) THEN
      READ(5,*) CP1, CP2, CP3, TREF, RGAS
C
      WRITE(6,*) ' IDEAL GAS PROPERTIES READ IN '
      WRITE(6,*) ' CP1, CP2, CP3, TREF, RGAS =',CP1,CP2,CP3,TREF,RGAS
C
      CPGAS  = CP1
      GAGAS  = CP1/(CP1 - RGAS)
      CP     = CP1
      GA     = GAGAS
      CV     = CP/GA
      IFGAS  = 1
      CALL SET_COEFFS
      END IF
C
C******************************************************************************
C******************************************************************************
C
      READ(5,1100)  DAMPIN,DUMM,FBLK1,FBLK2,FBLK3,SFEXIT,CONLIM,RFIN
      IF(SFEXIT.GT.0.0001) READ(5,*) NSFEXIT
      WRITE(6,1700) DAMPIN,DUMM,FBLK1,FBLK2,FBLK3,SFEXIT,CONLIM,RFIN
      IF(SFEXIT.GT.0.0001) WRITE(6,*) NSFEXIT
C
      IF(DAMPIN.LT.0.001)    DAMPIN = 10.
      IF(CONLIM.LT.0.000001) CONLIM = 0.005
      IF(RFIN.LT.0.000001)   RFIN = 0.10
      RFTHROTL = 0.0
C
C******************************************************************************
C******************************************************************************
C      READ IN INITIAL GUESS OF UPSTREAM AND DOWNSTREAM PRESSURES
C      ,STAGNATION TEMPERATURE , VELOCITIES,FLOW DIRECTIONS,ETC.
C
      READ(5,1100)  PUPHUB,PUPTIP,PDOWN_HUB,PDOWN_TIP,PLATE_LOSS,
     &              THROTTLE_EXIT,DUMM,F_PDOWN
      IF(THROTTLE_EXIT.GT.0.001) READ(5,*) THROTTLE_PRES, THROTTLE_MAS,
     &                                      RFTHROTL
      WRITE(6,1720) PUPHUB,PUPTIP,PDOWN_HUB,PDOWN_TIP,PLATE_LOSS,
     &              THROTTLE_EXIT,DUMM,F_PDOWN
      IF(THROTTLE_EXIT.GT.0.001) WRITE(6,*) THROTTLE_PRES,THROTTLE_MAS,
     &                                       RFTHROTL
C
      RFTHROTL1 = 1.0 - RFTHROTL
C
      KMID=IFIX(0.5*KIN)
C  Q3D
      IF(KM.EQ.2) KMID = 1
C  END Q3D
C
C******************************************************************************
C******************************************************************************
      WRITE(6,*)
      WRITE(6,*) ' READING IN THE INLET BOUNDARY CONDITIONS '
C
      READ(5,1100)  (PO1(K),K=1,KIN)
      WRITE(6,1720) (PO1(K),K=1,KIN)
      PO_IN_MID = PO1(KMID)
C
      IF(IPOUT.EQ.3) READ(5,1100)(PD(K),K=1,KIN)
      IF(IPOUT.EQ.3) WRITE(6,1720)(PD(K),K=1,KIN)
C
      READ(5,1100)  (TO1(K),K=1,KIN)
      WRITE(6,1700) (TO1(K),K=1,KIN)
      TO_IN_MID = TO1(KMID)
C
      READ(5,1100)  (VTIN(K),K=1,KIN)
      WRITE(6,1700) (VTIN(K),K=1,KIN)
      VT_IN_MID  = VTIN(KMID)
C
      READ(5,1100)  (VM1(K),K=1,KIN)
      WRITE(6,1700) (VM1(K),K=1,KIN)
C
      READ(5,1100)  (BS(K),K=1,KIN)
      WRITE(6,1700) (BS(K),K=1,KIN)
      YAW_IN_MID = BS(KMID)
C
      READ(5,1100)  (BR(K),K=1,KIN)
      WRITE(6,1700) (BR(K),K=1,KIN)
      PITCH_IN_MID  = BR(KMID)
C
      READ(5,1100)  (FR_IN(K),K=1,KIN-1)
      WRITE(6,1700) (FR_IN(K),K=1,KIN-1)
C
      READ(5,1100)  (FP(I),I=1,IMM1)
      WRITE(6,1700) (FP(I),I=1,IMM1)
C
      READ(5,1000)  (NOUT(L),L=1,5)
      WRITE(6,1600) (NOUT(L),L=1,5)
C
      READ(5,1610)  (IOUT(I),I=1,13)
      WRITE(6,1610) (IOUT(I),I=1,13)
C
      READ(5,1610)  (KOUT(K),K=1,KM)
      WRITE(6,1610) (KOUT(K),K=1,KM)
 1610 FORMAT(40I2)
C
C***************************************************************************************
C***************************************************************************************
C     FORM THE FP(I) AND FR_IN(K) INTO A GEOMETRIC SERIES IF REQUESTED.
C     BY SETTING FP(3) OR FR_IN(3) TO ZERO.
C
      NUM3 = 3
C
C   Q3D
      IF(KM.EQ.2) THEN
           FR_IN(1) = 1.0
           FR_IN(2) = 0.0
      ELSE
C   END Q3D
C
      IF(FR_IN(NUM3).GT.0.00001) GO TO 555
      FRRAT = FR_IN(1)
      FRMAX = FR_IN(2)
      FR_IN(1) = 1.0
      FR_IN(KM)  = 0.0
      FR_IN(KIN) = 0.0
      DO 556 K = 2,KIN-1
           FR_IN(K) = FR_IN(K-1)*FRRAT
  556 CONTINUE
      DO 557 K = 1,KIN-1
           FREV  = FR_IN(KIN-K)
           IF(FREV.LT.FR_IN(K))  FR_IN(K) = FREV
           IF(FR_IN(K).GT.FRMAX) FR_IN(K) = FRMAX
  557 CONTINUE
C
  555 CONTINUE
C
      END IF
C
C   NEXT FOR THE PITCHWISE GRID SPACINGS, FP(I) .
C   THROUGHFLOW
      IF(IM.EQ.2) THEN
            FP(1) = 1.0
            FP(2) = 0.0
      ELSE
C  END THROUGHFLOW
C
      IF(FP(NUM3).GT.0.00001) GO TO 563
      FPRAT = FP(1)
      FPMAX = FP(2)
      FP(1) = 1.0
      FP(IM)= 0.0
      DO 564 I=2,IMM1
      FP(I) = FP(I-1)*FPRAT
  564 CONTINUE
      DO 565 I=1,IMM1
      FREV  = FP(IM-I)
      IF(FREV.LT.FP(I))  FP(I) = FREV
      IF(FP(I).GT.FPMAX) FP(I) = FPMAX
  565 CONTINUE
C
  563 CONTINUE
C
      END IF
C
C***************************************************************************************
C***************************************************************************************
C
C     READ IN IN THE SPECIFIED INLET FLOW AND RELAXATION FACTOR
C     IF IN_FLOW NOT = ZERO
C
      IF(IN_FLOW.NE.0) THEN
           READ(5,1100)  FLOWIN, RFLOW
           WRITE(6,*)  ' FLOWIN, RFLOW ', FLOWIN,RFLOW
           FLOWIN = FLOWIN/NBLADE(1)
      END IF
C
C******************************************************************************
C******************************************************************************
C      READ IN DATA FOR VISCOUS FLOW MODELLING.
C
      REYNO       = 500000.0
      RF_VIS      = 0.5
      FTRANS      = 0.0001
      FAC_4TH     = 0.8
      PRANDTL     = 1.0
      YPLUSWALL   = 0.0
      TURBVIS_LIM = 3000.0
C
      IF(ILOS.NE.0) THEN
      READ(5,*,END=448,ERR=448)REYNO,RF_VIS,FTRANS,FAC_4TH,TURBVIS_LIM,
     &                         PRANDTL,YPLUSWALL
  448 CONTINUE
C
      IF(TURBVIS_LIM.LT.0.1) TURBVIS_LIM = 3000.0
C
      WRITE(6,*)
      WRITE(6,*) ' REYNOLDS NUMBER                       = ',REYNO
      WRITE(6,*) ' VISCOUS TERM RELAXATION FACTOR        = ',RF_VIS
      WRITE(6,*) ' TRANSITION FACTOR, FTRANS             = ',FTRANS
      WRITE(6,*) ' PROPORTION OF FOURTH ORDER SMOOTHING  = ', FAC_4TH
      WRITE(6,*) ' LIMIT ON TURBULENT/LAMINAR VISCOSITY  = ',TURBVIS_LIM
      WRITE(6,*) ' PRANDTL NUMBER                        = ',PRANDTL
      WRITE(6,*) ' YPLUSWALL - IF USED (NOT OFTEN USED)  = ',YPLUSWALL
      WRITE(6,*) ' MACH NUMBER  LIMITER, USUALLY 2.0     = ',MACHLIM
      WRITE(6,*)
C
      IF(ILOS.GE.200)THEN
           FAC_STMIX = 0.0
           FAC_ST0   = 1.0
           FAC_ST1   = 1.0
           FAC_ST2   = 1.0
           FAC_ST3   = 1.0
           FAC_SFVIS = 2.0
           READ(5,*,END = 450,ERR= 450) FAC_STMIX, FAC_ST0, FAC_ST1,
     &                                  FAC_ST2, FAC_ST3, FAC_SFVIS,
     &                                  FAC_VORT,FAC_PGRAD 
  450 CONTINUE
C
           WRITE(6,*)
           WRITE(6,*) ' SPALART-ALLMARAS TURBULENCE MODEL IS BEING USED'
           WRITE(6,*) ' THE S_A SOURCE TERM MULTIPLIERS ARE ',
     &                  FAC_STMIX, FAC_ST0, FAC_ST1, FAC_ST2, FAC_ST3,
     &                  FAC_VORT,FAC_PGRAD 
           WRITE(6,*) ' THE TURBULENT VISCOSITY SMOOTHING FACTOR IS ', 
     &                  FAC_SFVIS
           WRITE(6,*)
      END IF
C
      END IF 
C
      WRITE(6,*)
C
C***************************************************************************************
C**************************************************************************************
C     READ IN THE MIXING PLANE PARAMETERS
C
      RFMIX    = 0.025
      FSMTHB   = 1.0
      FEXTRAP  = 0.95
      FANGLE   = 0.95
C
      IF(IFMIX.NE.0) THEN
      READ(5,*,END=449,ERR=449)  RFMIX,FEXTRAP,FSMTHB,
     &                           FANGLE
  449 CONTINUE
      WRITE(6,*)' THE MIXING PLANE PARAMETERS ARE '
      WRITE(6,*)' RFMIX, FEXTRAP, FSMTHB ,FANGLE',
     &            RFMIX, FEXTRAP, FSMTHB, FANGLE
      WRITE(6,*)
      END IF
C
C***************************************************************************************
C***************************************************************************************
C    SET THE SPANWISE GRID SPACINGS,  "FR(K)" .
C   THESE ARE THE SAME AS THE INLET BC SPACINGS, "FR_IN(K)" ,  IF "IF_KINT" = 0.
C
      IF(IF_KINT.EQ.1) THEN
C
C     READ IN THE RELATIVE SPANWISE SPACING OF THE GRID POINTS.
C     NOTE THAT THIS IS FORMATTED INPUT.
      READ(5,10)    (FR(K),K=1,KMM1)
   10 FORMAT(8F10.5)
C
      WRITE(6,*)
      WRITE(6,*)       ' THE RELATIVE SPACINGS OF THE GRID POINTS IN THE
     &  SPANWISE DIRECTION ARE: '
      WRITE(6,1700)(FR(K),K=1,KMM1)
      WRITE(6,*)
C
C   MAKE THE SPANWISE SPACINGS INTO A GEOMETRICAL PROGRESSION IF FR(3) = 0 .
      NUM3 = 3
C
      IF(KM.EQ.2) THEN
           FR(1) = 1.0
           FR(2) = 0.0
      ELSE
C
      IF(FR(NUM3).GT.0.00001) GO TO 652
      FRRAT = FR(1)
      FRMAX = FR(2)
      FR(1)=1.0
      DO 650 K=2,KM-1
      FR(K) = FR(K-1)*FRRAT
  650 CONTINUE
      DO 651 K=1,KMM1
      FREV = FR(KM-K)
      IF(FREV.LT.FR(K))  FR(K) = FREV
      IF(FR(K).GT.FRMAX) FR(K) = FRMAX
  651 CONTINUE
C
  652 CONTINUE
C
      END IF
C
C
      WRITE(6,*) ' CALLING INPINT TO SET UP A NEW INLET FLOW. '
C
            CALL INPINT
C
      ELSE
C
      DO 655 K=1,KMM1
           FR(K) = FR_IN(K)
  655 CONTINUE
C  END OF  IF IF_KINT = 1 LOOP
      ENDIF
C
C***************************************************************************************
C***************************************************************************************
C      READ IN THE MIXING LENGTH LIMITS IF ILOS IS NOT ZERO.
C
      DO 2344 N = 1,NROWS
           XLLIM_I1(N)  = 0.03
           XLLIM_IM(N)  = 0.03
           XLLIM_K1(N)  = 0.03
           XLLIM_KM(N)  = 0.03
           XLLIM_DWN(N) = 0.03
           XLLIM_UP(N)  = 0.02
           XLLIM_IN(N)  = 0.02
           XLLIM_LE(N)  = 0.03
           XLLIM_TE(N)  = 0.04
           XLLIM_DN(N)  = 0.05
           FSTURB(N)    = 1.0
           TURBVIS_DAMP(N) = 0.5
 2344 CONTINUE
C
C
      IF(ILOS.NE.0) THEN
C
      DO 2345 N = 1,NROWS

           IF(ILOS.LT.100) READ(5,*,ERR=2346)  XLLIM_I1(N),XLLIM_IM(N),
     &                 XLLIM_K1(N),XLLIM_KM(N),XLLIM_DWN(N),XLLIM_UP(N)
           IF(ILOS.GE.100) READ(5,*,ERR=2346)  XLLIM_IN(N),XLLIM_LE(N),
     &                XLLIM_TE(N),XLLIM_DN(N),FSTURB(N),TURBVIS_DAMP(N)
C           
 2346 CONTINUE
C
           WRITE(6,*) ' ROW NUMBER ', N
           IF(ILOS.LT.100) WRITE(6,*)'ILOS = 9/10 MIXING LENGTH LIMITS',
     &     XLLIM_I1(N),XLLIM_IM(N),XLLIM_K1(N),XLLIM_KM(N),XLLIM_DWN(N),
     &     XLLIM_UP(N)
           IF(ILOS.GE.100) WRITE(6,*) 'ILOS > 100 MIXING LENGTH LIMITS',
     &     XLLIM_IN(N),XLLIM_LE(N),XLLIM_TE(N),XLLIM_DN(N)
           IF(ILOS.GE.100) WRITE(6,*) 
     &     ' FREE STREAM TURBULENT VISCOSITY RATIO',FSTURB(N),
     &     ' MIXING PLANE TURBULENCE DECAY', TURBVIS_DAMP(N)
C
 2345 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     READ IN A FACTOR TO INCREASE THE TURBULENT VISCOSITY FOR THE FIRST NMIXUP STEPS.
C
           FACMIXUP = 2.0
           NMIXUP   = 1000
           READ(5,*, ERR= 2341)    FACMIXUP, NMIXUP
           IF(FACMIXUP.LT.1.0) FACMIXUP = 1.0
           IF(IF_RESTART.NE.0)     FACMIXUP = 1.0
 2341 CONTINUE
           WRITE(6,*) ' FACMIXUP = ', FACMIXUP,' NMIXUP = ',NMIXUP
C
      ENDIF
C
C******************************************************************************
C******************************************************************************
C     WRITE OUT THE VALUE OF YPLUS AT THE WALL IF YPLUSWALL > 5.0.
C
      IF(YPLUSWALL.GT.5.0)  THEN
         WRITE(6,*)'YPLUSWALL IS BEING USED TO OBTAIN THE SKIN FRICTION'
         WRITE(6,*)'YPLUS AT THE WALL TAKEN AS, ', YPLUSWALL
           CFWALL = 1/(YPLUSWALL*YPLUSWALL)
      ENDIF
C
C******************************************************************************
C******************************************************************************
C  READ IN THE SURFACE ROUGHNESSES IN MICRONS IF  IF_ROUGH  >  0  .
C
      IF(IF_ROUGH.GT.0) THEN
C
	WRITE(6,*) ' Non-hydraulic smooth surfaces specified'
	WRITE(6,*) ' Input Surface Roughness in microns for ALL'
	WRITE(6,*) ' surfaces.'
      DO 2347 N = 1,NROWS
           READ(5,*)  ROUGH_H(N),ROUGH_T(N),ROUGH_L(N),ROUGH_U(N)
C
           WRITE(6,*) ' ROW NUMBER ', N
           WRITE(6,*) ' Surface Roughnesses in microns ',
     &     ROUGH_H(N),ROUGH_T(N),ROUGH_L(N),ROUGH_U(N)
C
C   Change to physical roughness in metres . 
           ROUGH_H(N) = ROUGH_H(N)*1.0E-06
           ROUGH_T(N) = ROUGH_T(N)*1.0E-06
           ROUGH_L(N) = ROUGH_L(N)*1.0E-06
           ROUGH_U(N) = ROUGH_U(N)*1.0E-06
C
 2347 CONTINUE
C
      ELSE
C
      DO 2348 N = 1,NROWS
           ROUGH_H(N)  = 0.0
           ROUGH_T(N)  = 0.0
           ROUGH_L(N)  = 0.0
           ROUGH_U(N)  = 0.0
 2348 CONTINUE
C
      END IF
C
C**********************************************************************************
C**********************************************************************************
C    INPUT THE ARTIFICIAL SPEED OF SOUND IF ITIMST = 5.
C    THIS SHOULD BE ABOUT HALF THE MAXIMUM RELATIVE VELOCITY IN THE FLOW.
C
      IF(ITIMST.GE.5) THEN
      VSOUND    = 150.
      RF_PTRU   = 0.01
      RF_VSOUND = 0.002
      DENSTY    = 1.20
      VS_VMAX   = 2.0
      IF(ITIMST.EQ.5)  READ(5,*,END= 2350)
     &  VSOUND, RF_PTRU, RF_VSOUND, VS_VMAX
      IF(ITIMST.EQ.6)  READ(5,*,END= 2350)
     &  VSOUND, RF_PTRU, RF_VSOUND, VS_VMAX, DENSTY
 2350 CONTINUE
         RF_PTRU1    = 1.0 - RF_PTRU
         RF_VSOUND1  = 1.0 - RF_VSOUND
         WRITE(6,*)
         WRITE(6,*) ' CALCULATION USING ARTIFICIAL COMPRESSIBILITY '
         WRITE(6,*) ' ARTIFICIAL SPEED OF SOUND = ', VSOUND
         WRITE(6,*) ' DENSITY RELAXATION FACTOR = ', RF_PTRU
         WRITE(6,*) ' SOUND SPEED RELAXATION FACTOR = ', RF_VSOUND
         WRITE(6,*) ' RATIO OF SOUND SPEED TO MAXIMUM SPEED = ',VS_VMAX
         IF(ITIMST.EQ.6) WRITE(6,*) 
     &   ' INCOMPRESSIBLE FLOW WITH DENSITY = ',DENSTY
         WRITE(6,*)
      END IF
C
C******************************************************************************
C******************************************************************************
C     READ IN THE OPTION TO USE REPEATING FLOW CONDITIONS
      IF(IF_REPEAT.NE.0) THEN
           READ(5,*)    NINMOD, RFINBC
           WRITE(6,*) ' REPEATING STAGE SPECIFIED, NINMOD, RFINBC = ',
     &                  NINMOD, RFINBC 
      END IF
C
C******************************************************************************
C******************************************************************************
C
C     READ IN THE STAGE NUMBERS AND SORT OUT THE START AND END OF EACH STAGE.
C     FIRST SET DEFAULTS.
      DO  N = 1,NROWS
      NSTAGE(N)    = 1 + (N-1)/2
      END DO
C
C     NOW READ IN THE ACTUAL STAGE NUMBER FOR EACH BLADE ROW.
C
      READ(5,*,END=3456,ERR=3456)(NSTAGE(N),N=1,NROWS)
 3456 CONTINUE
C
      WRITE(6,*) ' READ IN THE STAGE NUMBERS FOR EACH BLADE ROW'
      DO N= 1,NROWS
           WRITE(6,*) ' ROW NUMBER ',N,'IS IN STAGE NUMBER ',NSTAGE(N)
      END DO
C
C  NOW SET THE START AND END POINTS OF EACH STAGE
C
      JSTG_START(1)   = 1
      NSTAGES         = NSTAGE(NROWS)
      NSTAGE(NROWS+1) = NSTAGES + 1 
      JSTART(NROWS+1) = JM +1
C
      DO N = 1,NROWS
      NSTG   = NSTAGE(N)
      NSTGP1 = NSTAGE(N+1)
      IF(NSTG.NE.NSTGP1) THEN
             JSTG_END(NSTG)     = JMIX(N)
             JSTG_START(NSTGP1) = JSTART(N+1)
      END IF
      WRITE(6,*) 'ROW NUMBER ',N,'STAGE NUMBER ',NSTG,'JSTART ',
     &             JSTART(N),'JEND ',JMIX(N)
      END DO
C
      DO N = 1,NSTAGES
      WRITE(6,*) ' STAGE NUMBER ',N,' JSTART= ',JSTG_START(N),
     &           ' JEND = ',JSTG_END(N)
      END DO 
C
C**********************************************************************************
C**********************************************************************************
C   INPUT THE FORCING FACTOR AND SMOOTHING FACTOR IF DOING A THROGHFLOW  CALCULATION
C
      IF(IM.EQ.2) THEN
      Q3DFORCE = 1.0
      SFPBLD   = 0.1
      NSFPBLD  = 2
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=591) Q3DFORCE, SFPBLD, NSFPBLD
  591 CONTINUE
      SFPBLD1 = 1.0 - SFPBLD
      WRITE(6,*) 'THROUGHFLOW CALCULATION REQUESTED, Q3DFORCE= ',
     & Q3DFORCE, 'SMOOTHING FACTOR=',SFPBLD,' No OF SMOOTHINGS=',NSFPBLD
      WRITE(6,*)
      END IF       
C
C******************************************************************************
C******************************************************************************
C    INPUT THE RANGE OF YPLUS VALUES OVER WHICH THE TURBULENT VISCOSITY WILL BE REDUCED.
C
      IF(ILOS.GT.0) THEN
           YPLAM    = 5.0
           YPTURB   = 25.0
           READ(5,*, ERR= 3458,END = 3458) YPLAM,YPTURB
 3458 CONTINUE
           WRITE(6,*) ' TURBULENT VISCOSITY REDUCED OVER THE RANGE YPLUS
     & = ',YPLAM,' TO',YPTURB
      END IF
C
C******************************************************************************
C******************************************************************************
      WRITE(6,*)  ' Subroutine OLD_READIN completed, input data OK'
C
C
      RETURN
      END
C
C*************************************************************************
C*************************************************************************
C*************************************************************************
C*************************************************************************
C
      SUBROUTINE LOOP
C
C      THIS IS THE MAIN TIME STEPPING LOOP. IT IS EXECUTED MANY HUNDREDS OF
C      TIMES AND USES MOST OF THE CPU TIME. ITS MAIN SUBROUTINE IS 'TSTEP'
C      WHICH ALSO USES MUCH OF THE TIME.
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION CHECK_FLOW(JD),BLADE_FLOW(JD),ROVMSQ(KD),ANG_INC(KD)
      SAVE START,FINI,CHECK_FLOW,BLADE_FLOW
C
C****************************************************************************** 
C    SET SOME REFERENCE VALUES
C
      RO_REF = PO1(KMID)/(RGAS*TO1(KMID))
      IF(ITIMST.EQ.6) RO_REF = DENSTY
      ROMIN  = AMIN1(RO(IMID,1,KMID),RO(IMID,JM,KMID))
      ROLIM  = 0.25*ROMIN
      RIJKM  = 2./(IM*JM*KM)
C
      DO 4001 K=1,KM
      PLIM(K) = 0.9*PO1(K) + 0.1*P(1,1,K)
 4001 CONTINUE
C
C     SET THE REFERENCE PRESSURE. THIS IS SUBTRACTED FRON THE TRUE PRESSURE TO
C     MINIMISE ROUNDING ERRORS.
      PREFF = 0.5*(P(IMID,1,KMID) + P(IMID,JM,KMID))
      WRITE(6,*) ' REFERENCE PRESSURE IN BAR = ', PREFF/1.0E05
C
C****************************************************************************** 
C****************************************************************************** 
C    SET MORE CONSTANTS
C
C     MAKE THE SMOOTHING PROPORTIONAL TO THE CFL NUMBER
      SFXIN  = SFXIN*CFL/0.5
      SFTIN  = SFTIN*CFL/0.5
      RFIN1  = 1.0 - RFIN
      RFIN2  = 2.*RFIN
C
      ECONT    = 1.
      EAVG     = 1.0E06
      EMAX     = 1.0E06
      HBLEDTOT = 0.0
      TCOOL_MIN = 1.0E6
C
      IFEND    = 0
      NSTEP    = 0
C
      IBLEED = 0
      KBLEED = 0
      NCOOLB = 0
      NCOOLW = 0
C
C******************************************************************************
C     CALL COOLIN TO SET THE COOLANT SOURCE TERMS
C
      IF(IFCOOL.NE.0) THEN
           CALL COOL_INPUT
           WRITE(6,*) ' CALLED COOL_INPUT '
           IF(IFCOOL.EQ.1) CALL COOLIN_1
           IF(IFCOOL.EQ.2) CALL COOLIN_2
           WRITE(6,*) ' CALLED COOLIN '
      ENDIF
C
C******************************************************************************
C     CALL BLEEDOUT TO SET THE BLEED FLOW TERMS
C
      IF(IFBLEED.NE.0) THEN
           CALL BLEEDOUT(IBLEED,KBLEED)
           WRITE(6,*) ' CALLED BLEEDOUT '
      ENDIF
C
C
C******************************************************************************
C      CALL THE TIMING ROUTINES. THESE ARE MACHINE SPECIFIC AND MAY NEED CHANGING.
C      CALL CLOCK@(START)
       START = FLOAT(MCLOCK())
C
C*************************************************************************
C**********************START OF THE MAIN TIME STEPPING LOOP ********************
C*************************************************************************
C
      WRITE(6,*)
      WRITE(6,*)'*******************************************************
     &**********'
      WRITE(6,*) ' STARTING THE MAIN TIME STEPPING LOOP '
C
C       RETURN HERE AFTER EVERY TIME STEP
C
 5000 CONTINUE
C
C******************************************************************************
C     CHECK IF CONVERGED OR IF THE MAXIMUM NUMBER OF ITERATIONS HAS BEEN REACHED.
C     SET  "IFEND" = 1 TO STOP THE CALCULATION .
C
      IFPRINT = 0
      NSTEP  = NSTEP + 1
      N      = NSTEP
      IF(NSTEP.EQ.NSTEPS_MAX) IFEND = 1
      IF(EAVG.LT.CONLIM.AND.ECONT.LT.0.01)          IFEND = 1
      IF(ABS(EMAX).LT.(4*CONLIM).AND.ECONT.LT.0.01) IFEND = 1
C
C     CHECK IF PRINTOUT REQUIRED ON THIS TIME STEP
C
      DO 5010 I=1,5
 5010 IF(NSTEP.EQ.NOUT(I)) IFPRINT = 1
C
C    CHECK IF A REQUEST TO STOP HAS BEEN RECEIVED BY "stopit" BEING SET
C    GREATER THAN ZERO.
C
      IF(MOD(N,10).EQ.0) THEN
           OPEN(UNIT=12,FILE='stopit')
           READ(12,*) IFSTOP
           IF(IFSTOP.GE.1) IFEND = 1
           CLOSE(12)
      ENDIF 
C
C******************************************************************************
C      VARY THE SMOOTHING AND DAMPING OVER THE FIRST  "NCHANGE"  STEPS.
C      IF STARTING FROM AN INITIAL GUESS.
C
C     8/4/2017  MODIFY THE NEXT SECTION SO THAT DAMPING AND SMOOTHING ARE INCREASED
C     OVER THE FIRST NCHANGE STEPS EVEN WHEN STARTING FROM A RESTART FILE.
      IF(IF_RESTART.NE.0.AND.NSTEP.GT.100) GO TO 5555
C
      FCHANGE     = 1.0 -  FLOAT(NSTEP-1)/NCHANGE
C
      IF(FCHANGE.GT.0.0)   THEN
C     DELETE NEXT LINE 8/4/2017
C           IF(IF_RESTART.NE.0) FCHANGE = 0.0
           SFT  = SFTIN   + 0.02*FCHANGE
           SFX  = SFXIN   + 0.02*FCHANGE
           DAMP = DAMPIN*(1. - 0.75*FCHANGE)
      ENDIF
C
      SFT1    = 1.0-SFT
      SFTH    = 0.5*SFT
      SFX1    = 1.-SFX
      SFXH    = 0.5*SFX
      SFX15   = 1.0-1.5*SFX
      SFX14   = 0.25*SFX
      SFXB    = FSMTHB*SFX
      SFXB1   = 1.0 - SFXB
      SFXBH   = 0.5*SFXB
      SFXHM1  = 1.-SFXH
C
C    SET FMIXUP TO INCREASE THE TURBULENT VISCOSITY OVER THE FIRST NMIXUP STEPS
C
      FMIXUP    = 1.0 + (FACMIXUP-1.)*FLOAT(NMIXUP-NSTEP)/FLOAT(NMIXUP)
      IF(FMIXUP.LT.1.0) FMIXUP = 1.0
      FMIXLEN   = FMIXUP
      FMIXUP    = FMIXUP*FMIXUP
C
C    JUMP TO HERE IF STARTING FROM A RESTART FILE
C
 5555 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C
C       CALCULATE VX,VR,VT FROM ROVX,ROVR,ROVT
C
      DO 7450 K=1,KM
      DO 7450 J=2,JM
           RRNOW  = 1/R(J,K)
      DO 7450 I=1,IM
           RECIP       = 1.0/RO(I,J,K)
           VX(I,J,K)   = ROVX(I,J,K)*RECIP
           VR(I,J,K)   = ROVR(I,J,K)*RECIP
           ROVT(I,J,K) = RORVT(I,J,K)*RRNOW
           VT(I,J,K)   = ROVT(I,J,K)*RECIP
 7450 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C
      IF(MOD(NSTEP,5).EQ.0.OR.IFEND.EQ.1) THEN
C
      IF(IFCOOL.EQ.2) CALL COOLIN_2
C
C      EVERY 5 STEPS SUM THE RELATIVE VELOCITIES TO FIND A REFERENCE VELOCITY
C      AND THE MAXIMUM RELATIVE VELOCITY.
C
      VRMS=0.0
      VMAX=0.0
      DO 5910 K=1,KM
      DO 5910 J=2,JM
      DO 5910 I=1,IM
      EKE = VR(I,J,K)*VR(I,J,K) + WT(I,J,K)*WT(I,J,K) 
     &       + VX(I,J,K)*VX(I,J,K)
      IF(EKE.GT.VMAX) THEN
	 VMAX  = EKE
	 IVMAX = I
	 JVMAX = J
	 KVMAX = K
      ENDIF
      VRMS = VRMS + EKE
 5910 CONTINUE
      VMAX = SQRT(VMAX)
C
      ENDIF
C
C*******************************************************************************
C*******************************************************************************
C
C   CALCULATE NEW PRESSURES AND TEMPERATURES.THE METHOD DEPENDS ON THE FLOW MODEL, ITIMST.
C
      IF(ITIMST.LT.5) THEN
C
      DO 5900 K=1,KM
      DO 5900 J=2,JMM1
      DO 5900 I=1,IM
           RONOW = RO(I,J,K)
           EKE   =.5*(VR(I,J,K)*VR(I,J,K) + VT(I,J,K)*VT(I,J,K) 
     &           + VX(I,J,K)*VX(I,J,K))
C
      IF(IFGAS.EQ.0) THEN
           TSTATIC      = (ROE(I,J,K)/RONOW - EKE)/CV
           P(I,J,K)     = RONOW*RGAS*TSTATIC
           HO(I,J,K)    = CP*TSTATIC + EKE 
      ELSE
           ESTAT        = ROE(I,J,K)/RONOW - EKE
           TSTATIC      = TFROME(ESTAT,TREF,EREF,ET1,ET2,ET3,ET4)
C
           IF(TSTATIC.LT.0.1*TREF) TSTATIC = 0.1*TREF
           P(I,J,K)     = RONOW*RGAS*TSTATIC
           HO(I,J,K)    = ESTAT + RGAS*TSTATIC + EKE
      END IF
C
           T_STATIC(I,J,K) = TSTATIC
C
 5900 CONTINUE
C
      END IF
C
C
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) THEN
C
C     VARY THE ARTIFICIAL SPEED OF SOUND SO IT VERY GRADUALLY TENDS TO A VALUE
C      = 2 x THE MAXIMUM RELATIVE VELOCITY.
C
      VS_NEW   = VS_VMAX*VMAX
C      VSOUND   = RF_VSOUND1*VSOUND + RF_VSOUND*VS_NEW
      FRACDIF  =  VS_NEW/VSOUND - 1.0
      DVSOUND  =  VSOUND*FRACDIF*(1.0 + 3*ABS(FRACDIF) )
      IF(ABS(FRACDIF).LT.0.1) DVSOUND = 0.0
      VSOUND   =  VSOUND + RF_VSOUND*DVSOUND
      DP_DRO   =  VSOUND*VSOUND
C
      DO 5950 K=1,KM
      DO 5950 J=2,JM
      DO 5950 I=1,IM
      EKE          =.5*(VR(I,J,K)*VR(I,J,K) + VT(I,J,K)*VT(I,J,K) 
     &               + VX(I,J,K)*VX(I,J,K))
      P(I,J,K)     = PO_REF  - DP_DRO*(RO_REF - ROSUB(I,J,K))
      TSTATIC      = (ROE(I,J,K)/ROSUB(I,J,K) - EKE)/CV
      RO(I,J,K)    = RF_PTRU*P(I,J,K)/RGAS/TSTATIC 
     &             + RF_PTRU1*RO(I,J,K)
      IF(ITIMST.EQ.6) RO(I,J,K) = DENSTY 
      HO(I,J,K)    = CP*TSTATIC  + EKE 
      T_STATIC(I,J,K) = TSTATIC
 5950 CONTINUE
C     END OF ITIMST = 5  OR 6 LOOP
      END IF
C
C*******************************************************************************
C*******************************************************************************
C
C         CALL LOSS ROUTINES TO UPDATE BODY FORCES EVERY NLOS ITERATIONS
C
C     IF ILOS = 10  USE THE ORIGINAL MIXING LENGTH MODEL.
      IF((ILOS.EQ.10).AND.MOD(NSTEP,NLOS).EQ.0)  CALL LOSS
C
C     IF ILOS > 100 and < 200  USE THE NEW MIXING LENGTH MODEL.
      IF((ILOS.GE.100.AND.ILOS.LT.200).AND.MOD(NSTEP,NLOS).EQ.0)
     &    CALL NEW_LOSS
C
C     IF(ILOS > 200  USE THE SPALART ALLMARAS MODEL.  
      IF(ILOS.GE.200.AND.MOD(NSTEP,NLOS).EQ.0)   CALL SPAL_LOSS
C
C*******************************************************************************
C*******************************************************************************
C      CALL SHROUDFLOW TO SET SHROUD LEAKAGE FLUXES IF 'IFSHROUD' NE ZERO.
C
      IF((NSTEP.EQ.1.OR.MOD(NSTEP,NLOS).EQ.0).AND.(IFSHROUD.EQ.1))
     &  CALL SHROUDFLOW
C
C*******************************************************************************
C*******************************************************************************
C
C          WORK OUT THE MASS FLUXES AND STORE AS FLOWX, FLOWT, FLOWR.
C
C       FLOWR .... MASS FLOW RATE THROUGH THE STREAM-SURFACE (TOP) FACE
C       FLOWT .... MASS FLOW RATE THROUGH THE BLADEWISE (SIDE) FACE
C       FLOWX .... MASS FLOW RATE THROUGH THE QUASI-ORTHOGONAL (UPSTREAM) FACE
C
      DO 5020 K=1,KMM1
      DO 5020 J=1,JM
      DO 5020 I=1,IMM1
      AVGRVX = ROVX(I,J,K)+ROVX(I,J,K+1)+ROVX(I+1,J,K)+ROVX(I+1,J,K+1)
      AVGRVR = ROVR(I,J,K)+ROVR(I,J,K+1)+ROVR(I+1,J,K)+ROVR(I+1,J,K+1)
      FLOWX(I,J,K)   = 0.25*(AVGRVX*AQX(I,J,K) + AVGRVR*AQR(I,J,K))
      SOURCE(I,J,K ) = 0.0
      XFLUX(I,J,K)   = FLOWX(I,J,K)
 5020 CONTINUE
C
C     CALCULATE AND STORE  ROWT(I,J,K) AND WT(I,J,K) 
C
      DO 5025 K=1,KM
      DO 5025 J=1,JM
      DO 5025 I=1,IM
      ROWT(I,J,K) = ROVT(I,J,K) - UBLADE(J,K)*RO(I,J,K)
      WT(I,J,K)   = ROWT(I,J,K)/RO(I,J,K)
 5025 CONTINUE
C
      DO 5030 K=1,KMM1
      DO 5030 J=2,JM
      DO 5030 I=1,IM
      AVGRVX = ROVX(I,J,K)+ROVX(I,J-1,K)+ROVX(I,J-1,K+1)+ROVX(I,J,K+1)
      AVGRVT = ROWT(I,J,K)+ROWT(I,J-1,K)+ROWT(I,J-1,K+1)+
     &         ROWT(I,J,K+1)
      AVGRVR = ROVR(I,J,K)+ROVR(I,J-1,K)+ROVR(I,J-1,K+1)+ROVR(I,J,K+1)
      AVGRO  = RO(I,J,K)+RO(I,J-1,K)+RO(I,J-1,K+1)+RO(I,J,K+1)
      FLOWT(I,J,K) = (AVGRVX*ABX(I,J,K) + AVGRVT*ABT(J,K)
     &             +  AVGRVR*ABR(I,J,K))*0.25 
      TFLUX(I,J,K) = FLOWT(I,J,K) 
 5030 CONTINUE
C
C      BALANCE MASS FLUXES INTO AND OUT OF CUSPS
C
      DO 5050 J=2,JM
      IF(IND(J).EQ.1) GO TO 5050
      DO 5051 K=1,KMM1
      AVGFLUX        = 0.5*(FLOWT(1,J,K) + FLOWT(IM,J,K))
      FLOWT(1,J,K)   = AVGFLUX 
      FLOWT(IM,J,K)  = AVGFLUX
      TFLUX(1,J,K)   = AVGFLUX
      TFLUX(IM,J,K)  = AVGFLUX
 5051 CONTINUE
 5050 CONTINUE
C
C**********************************************************************
C    TFLOW ADDITION
C**********************************************************************
      IF(IM.NE.2) GO TO 5655
C
      CSQ = GA*RGAS*TO1(KMID) * Q3DFORCE
C
C  SET THE FLOW ANGLE WITHIN THE BLADE PASSAGE USING THE DEVIATION FROM THE GRID ANGLE.
C  THE DEVIATION VARIES LINEARLY WITH THE J INDEX.  
C  GRADUALLY CHANGE THE INCIDENCE OVER "FRAC_INC" OF THE J VALUES FROM THE LEADING EDGE.
C  ADD THE DEVIATION OVER "FRAC_DEV" OF THE J VALUES FROM THE TRAILING EDGE.
C
      FRAC_INC = 0.3333
      FRAC_DEV = 0.5
C
      DO 5552 K=1,KMM1
C
      DO 5553 J= 2,JM
C
      NR    = NROW(J) 
      JLEE  = JLE(NR)
      JLEM1 = JLEE-1
      JLEP1 = JLEE+1
      JTEE  = JTE(NR) 
C
      FDEVN = FLOAT(J-JLEE)/FLOAT(JTEE - JLEE)
      FDEVN = 1.0 - (1.0 - FDEVN)/FRAC_DEV
      IF(FDEVN.GT.1.00001)    FDEVN = 1.0 
      IF(FDEVN.LT.1.0E-4)     FDEVN = 0.0
C
      FINC = FLOAT(J - JLEE)/FLOAT(JTEE - JLEE)/FRAC_INC
      FINC =  1.0 - FINC
      IF(FINC.GT.1.00001)    FINC = 1.0
      IF(FINC.LT.1.0E-4)     FINC = 0.0
C 
      FLOAD = 1.0
      IF(J.LE.JLEM1)  FLOAD = 0.0
      IF(J.GT.JTEE+1) FLOAD = 0.0
C 
C   CALCULATE THE FLOW ANGLE AT THE LEADING EDGE      
C
           ROWTAVG    = 0.25*(ROWT(1,J,K)   + ROWT(1,J,K+1)
     &                + ROWT(1,J-1,K) + ROWT(1,J-1,K+1) )
           XD = DX(J,K)
           RD = DR(J,K)
           SD = DS(J,K)
           ROVMK = (XD*(ROVX(1,J,K) + ROVX(1,J-1,K))
     &           +  RD*(ROVR(1,J,K) + ROVR(1,J-1,K)))/SD
           XD = DX(J,K+1)
           RD = DR(J,K+1)
           SD = DS(J,K+1)
           ROVMKP1  = (XD*(ROVX(1,J,K+1) + ROVX(1,J-1,K+1))
     &              +  RD*(ROVR(1,J,K+1) + ROVR(1,J-1,K+1)))/SD
           ROVMAVG  = 0.25*(ROVMK + ROVMKP1) 
C
           FLO_ANGL    = ATAN(ROWTAVG/ROVMAVG)
C 
C   SET THE INCIDENCE ANGLE BASED ON THE CENTRE LINE ANGLE AT THE LEADING EDGE
C    
      IF(J.EQ.JLEE)   ANG_INC(K)  = FLO_ANGL  - ALPHA_CENT(J,K)
C
      IF(J.LT.JLEE.OR.J.GT.JTEE) ANG_INC(K)  = 0.0 

C    SET THE ANGLE TO WHICH THE FLOW IS FORCED AS A COMBINATION OF THE CENTRE
C    LINE ANGLE, THE INCIDENCE ANGLE AND THE DEVIATION ANGLE.
C
      GRID_ANGL    = ALPHA_CENT(J,K)
      DEV_ANGL     = DEVN_ANGL(NR,K)*FDEVN  
      ANGL_INC     = ANG_INC(K)*FINC
      ALPHA_FORCED = ALPHA_CENT(J,K) - DEV_ANGL + ANGL_INC
C  
C     CALCULATE THE CHANGE IN BLADE SURFACE PRESSURE DUE TO ANY FLOW ACROSS
C     THE FORCED ANGLE SURFACE.
C 
      AVGTFLUX1  = 0.5*(TFLUX(1,J,K) + TFLUX(2,J,K) ) 
      AVGTFLUX2  = AVGTFLUX1 + ROVMAVG*ABT(J,K)*(TAN(ALPHA_CENT(J,K))
     &           - TAN(ALPHA_FORCED))
      ROVARN        = AVGTFLUX2*STEP(1,J,K)
      PVARN         = CSQ*(2*ROVARN - ROVAR_M1(J,K) )
      ROVAR_M1(J,K) = ROVARN
      PBLADE(J,K)   = PVARN + PBLADE(J,K)
      PBLADE(J,K)   = PBLADE(J,K)*FLOAD
      IF(K.GE.KTIPS(NR).AND.K.LT.KTIPE(NR) ) PBLADE(J,K) = 0.0
 5553 CONTINUE
C
 5552 CONTINUE
C
C    SMOOTH THE BLADE SURFACE PRESSURE DISTRIBUTION.
C
      DO 5651 N = 1,NSFPBLD
      DO 5652 K=1,KMM1
      DO 5653 J=3,JM-3
      SLEFT  = SMERID(J,K)   - SMERID(J-1,K)
      SRIGHT = SMERID(J+1,K) - SMERID(J,K)
      AVGPBLADE(J) = SFPBLD1*PBLADE(J,K)
     &             + SFPBLD*(SRIGHT*PBLADE(J-1,K) + SLEFT*PBLADE(J+1,K)) 
     &             /(SLEFT + SRIGHT)
 5653 CONTINUE
      DO 5654 J=3,JM-3
      PBLADE(J,K) = AVGPBLADE(J)    
 5654 CONTINUE
 5652 CONTINUE
 5651 CONTINUE
C
C
 5655 CONTINUE
C
C**********************************************************************
C   END TFLOW ADDITION
C**********************************************************************
C
C      SET ZERO MASS FLOW THROUGH THE BLADE SURFACES - WHERE NO BLEED.
C
      DO 5060 J=2,JM
      IF(IND(J).NE.1) GO TO 5060
      NR= NROW(J)      
      DO 5061 K=KS1(NR),KS2(NR)
      TFLUX(1,J,K)  = 0.0
      TFLUX(IM,J,K) = 0.0
      FLOWT(1,J,K)  = 0.0
      FLOWT(IM,J,K) = 0.0
 5061 CONTINUE
 5060 CONTINUE
C
C    REMOVE ANY BLEED FLOWS THROUGH THE BLADE SURFACES
C
      IF(IBLEED.NE.0) THEN
      DO 5065 J = 2,JM
      IF(IND(J).NE.1) GO TO 5065
      DO 5064 K = 1,KMM1
      TFLUX(1,J,K)  = -BLFLOWI1(J,K)
      TFLUX(IM,J,K) =  BLFLOWIM(J,K)
      FLOWT(1,J,K)  = -BLFLOWI1(J,K)
      FLOWT(IM,J,K) =  BLFLOWIM(J,K)
 5064 CONTINUE
 5065 CONTINUE
      ENDIF
C          
C    NOW SET RFLUX TO THE MASS FLOW THROUGH THE STREAMWISE FACES.
C  Q3D  
      IF(KM.NE.2) THEN
           K1 = 2
           K2 = KMM1
      ELSE
C   END Q3D
           K1 = 1
           K2 = 2
      END IF
C
      DO 5040 K=K1,K2
      DO 5040 J=2,JM
      DO 5040 I=1,IMM1
      AVGRVX = ROVX(I,J,K)+ROVX(I,J-1,K)+ROVX(I+1,J-1,K)+ROVX(I+1,J,K)
      AVGRVR = ROVR(I,J,K)+ROVR(I,J-1,K)+ROVR(I+1,J-1,K)+ROVR(I+1,J,K)
      FLOWR(I,J,K)    = 0.25*(AVGRVX*ASX(I,J,K)+AVGRVR*ASR(I,J,K))
      RFLUX(I,J,K)    = FLOWR(I,J,K)
 5040 CONTINUE
C
C***************************************************************************************
C     SET THE BODY FORCE TO KEEP THE FLOW ON THE STREAM SURFACE IF DOING A Q3D CALCULATION.
C   Q3D
      IF(KM.NE.2) GO TO 5042
C
      CSQ = GA*RGAS*TO1(1)*Q3DFORCE
      DO 5041 J=2,JM
      DO 5041 I=1,IMM1
      AVRFLUX     = RFLUX(I,J,1) + RFLUX(I,J,2)
      ROVARN      = AVRFLUX*STEP(I,J,1)
      PVARN       = CSQ*ROVARN
      XFORCE(I,J,1) = XFORCE(I,J,1) + PVARN*(ASX(I,J,1)+ASX(I,J,2))
      RFORCE(I,J,1) = RFORCE(I,J,1) + PVARN*(ASR(I,J,1)+ASR(I,J,2))
 5041 CONTINUE
C
 5042 CONTINUE
C      
C  END Q3D
C***************************************************************************************
C      SET ZERO MASS FLOW THROUGH THE HUB AND CASING - WHERE NO BLEED.
C
      DO 5045 J=2,JM
      DO 5045 I=1,IMM1
      FLOWR(I,J,1)  = 0.0
      FLOWR(I,J,KM) = 0.0
      RFLUX(I,J,1)  = 0.0
      RFLUX(I,J,KM) = 0.0
 5045 CONTINUE
C
C   Q3D
      IF(KM.EQ.2) GO TO 5081
C   END Q3D
C
C     REMOVE ANY BLEED FLOWS THROUGH THE HUB OR CASING
C
      IF(KBLEED.NE.0) THEN
      DO 5080 J=1,JM
      DO 5080 I=1,IMM1
      FLOWR(I,J,1)     = -BLFLOWK1(I,J)
      FLOWR(I,J,KM)    =  BLFLOWKM(I,J)
      RFLUX(I,J,1)     =  FLOWR(I,J,1)
      RFLUX(I,J,KM)    =  FLOWR(I,J,KM)
 5080 CONTINUE
      ENDIF
C
 5081 CONTINUE
C
C     ADD COOLANT MASS FLUXES THROUGH THE BLADE SURFACES.
C
      IF(NCOOLB.NE.0) THEN
      DO 5550 NC=1,NCOOLB
      DO 5556 J=JCBS(NC)+1,JCBE(NC)
      DO 5556 K=KCBS(NC),KCBE(NC)-1
      IF(IC(NC).EQ.1)  SOURCE(1,J,K)     =   - CFLOWI1(J,K)
      IF(IC(NC).EQ.IM) SOURCE(IMM1,J,K)  =   - CFLOWIM(J,K)
 5556 CONTINUE
 5550 CONTINUE
      ENDIF
C
C
C  Q3D
      IF(KM.EQ.2) GO TO 5557
C  END Q3D
C
C     ADD COOLANT MASS FLOWS THROUGH THE HUB AND CASING
C
      IF(NCOOLW.NE.0) THEN
      DO 5558 NC=1,NCOOLW
      DO 5559 J=JCWS(NC)+1,JCWE(NC)
      DO 5559 I=ICWS(NC),ICWE(NC)-1     
      IF(KC(NC).EQ.1) SOURCE(I,J,1)    = SOURCE(I,J,1)    - CFLOWK1(I,J)
      IF(KC(NC).EQ.KM)SOURCE(I,J,KMM1) = SOURCE(I,J,KMM1) - CFLOWKM(I,J)
 5559 CONTINUE
 5558 CONTINUE
      ENDIF
C
 5557 CONTINUE
C
C***********************************************************************
C      SET THE SHROUD LEAKAGE MASS FLUXES
C
      IF(IFSHROUD.EQ.1) CALL SHROUDFLUX(SHROUDGR)
C
C*******************************************************************
C*******************************************************************
C       CALL SUBROUTINE TSTEP TO UPDATE THE DENSITY AT ALL POINTS
C
      IF(ITIMST.GE.5) THEN
           CALL TSTEP(ROSUB,DRO,1)
      ELSE
           CALL TSTEP(RO,DRO,1)
      END IF
C
C***********************************************************************
C***********************************************************************
C   FORM THE CELL AVERAGE DENSITY IF USING NEWLOSS OR SPAL_LOSS.
C
      IF(ILOS.GE.100) THEN
C
      DO 5600 K=1,KMM1
      DO 5600 J=2,JM
      DO 5600 I=1,IMM1
      ROAVG = 0.125*(RO(I,J,K)+RO(I+1,J,K)+RO(I+1,J,K+1)+RO(I,J,K+1)
     &      + RO(I,J-1,K)+RO(I+1,J-1,K)+RO(I+1,J-1,K+1)+RO(I,J-1,K+1))
      ROAVG_CELL(I,J,K) = ROAVG
 5600 CONTINUE
C
      END IF
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C     SET THE FLUXES OF TURBULENT VISCOSITY FOR THE SA MODEL ONLY
C     NOTE THAT THESE AND THE SOURCE TERM ARE DOUBLED TO GIVE AN EFFECTIVE
C     LARGER TIME STEP FOR THE TURBULENT VISCOSITY EQUATION.
C
      IF(ILOS.GE.200) THEN
C
      DO 6001 K = 1,KMM1
      DO 6001 J = 1,JM
      JP1 = J+1
      IF(J.EQ.JM) JP1 = JM
      IF(J.EQ.1)  JP1  = 1
      IF(INDMIX(J).EQ.1) JP1 = J
      DO 6001 I = 1,IMM1
      XFLUX(I,J,K) = FLOWX(I,J,K)*(TRANS_KVIS(I,J,K)
     &             + TRANS_KVIS(I,JP1,K))
 6001 CONTINUE
C
      DO 6008 K=1,KMM1
      DO 6008 J=2,JM
      DO 6008 I=1,IMM1
      SOURCE(I,J,K) = -2.0*T_SOURCE(I,J,K)
 6008 CONTINUE
C
C
      DO 6002 K=1,KMM1
      DO 6002 J=2,JM
      DO 6002 I=2,IMM1
      TFLUX(I,J,K) = FLOWT(I,J,K)*
     &               (TRANS_KVIS(I-1,J,K)+TRANS_KVIS(I,J,K))
 6002 CONTINUE
      DO 6003 K=1,KMM1
      DO 6003 J= 2,JM
      TFLUX(1,J,K) = FLOWT(1,J,K)*
     &               (TRANS_KVIS(IMM1,J,K)+TRANS_KVIS(1,J,K))
      TFLUX(IM,J,K) = TFLUX(1,J,K)
 6003 CONTINUE
C
C  Q3D
      IF(KM.EQ.2) GO TO 6005
C  END Q3D
C
      DO 6004 K=2,KMM1
      DO 6004 J=2,JM
      DO 6004 I=1,IMM1
      RFLUX(I,J,K) = FLOWR(I,J,K)*
     &               (TRANS_KVIS(I,J,K-1) + TRANS_KVIS(I,J,K))
 6004 CONTINUE
C
 6005 CONTINUE
C
      DO 6006 J=2,JM
      DO 6006 I=1,IMM1
      RFLUX(I,J,1) = 0.0
      RFLUX(I,J,KM)= 0.0
 6006 CONTINUE
C
C********************************************************************************
C********************************************************************************
C     CALL TSTEP TO UPDATE THE TURBULENT VISCOSITY FOR THE SA MODEL ONLY.
C
      CALL TSTEP(TRANS_DYN_VIS, DEL_DYNVIS,2)
C
C********************************************************************************
C********************************************************************************
C  
      TRANSVISMIN = 0.1*VISC_LAM(IMID,2,KMID)
      DO 6007 K=1,KMM1
      DO 6007 J=2,JM
      DO 6007 I=1,IMM1
      IF(TRANS_DYN_VIS(I,J,K).LT.TRANSVISMIN) 
     &  TRANS_DYN_VIS(I,J,K) = TRANSVISMIN
      TRANS_KVIS(I,J,K) = TRANS_DYN_VIS(I,J,K)/ROAVG_CELL(I,J,K)
 6007 CONTINUE
C
C     END OF SETTING THE TURBULENT VISCOSITY FLUXES FOR THE SA MODEL.
C
      END IF
C
C******************************************************************************
C********************************************************************************
C   NOW WORK OUT THE FLUXES OF ENERGY AND STORE AS XFLUX, TFLUX AND RFLUX.
C
      DO 5501 K=1,KMM1
      DO 5501 J=1,JM
      DO 5501 I=1,IMM1
      AVGHO = HO(I,J,K)+HO(I,J,K+1)+HO(I+1,J,K)+HO(I+1,J,K+1)
      XFLUX(I,J,K) = FLOWX(I,J,K)*AVGHO*0.25
 5501 CONTINUE
C
      DO 5505 J=2,JM
      DO 5505 K=1,KMM1
      DO 5505 I=1,IMM1
      SOURCE(I,J,K) = QSOURCE(I,J,K)
 5505 CONTINUE
C
C  TFLOW ADDITION
      IF(IM.EQ.2) THEN
      DO 5506 K=1,KMM1
      DO 5506 J=2,JM
      SOURCE(1,J,K) = SOURCE(1,J,K) + BFORCE_Q(J,K)
 5506 CONTINUE
      END IF
C   END TFLOW ADDITION
C
      DO 5502 K=1,KM
      DO 5502 J=2,JM
      DO 5502 I=1,IMM1
      AVGHO = 0.25*(HO(I,J,K)+HO(I,J-1,K)+HO(I+1,J-1,K)+HO(I+1,J,K))
      RFLUX(I,J,K) = FLOWR(I,J,K)*AVGHO
 5502 CONTINUE
C
C         
      DO 5503 K=1,KMM1
      DO 5503 J=2,JM
      DO 5503 I=1,IM
      AVGHO = HO(I,J,K)+HO(I,J-1,K)+HO(I,J-1,K+1)+HO(I,J,K+1)
      AVGPB = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I,J-1,K+1)+PEFF(I,J,K+1) 
      TFLUX(I,J,K) = (FLOWT(I,J,K)*AVGHO + WRABT(J,K)*AVGPB)*0.25
 5503 CONTINUE
C
C
C   CALCULATE THE ENTHALPY BLED FROM THE MAIN FLOW - EVERY 200 STEPS.
C
      IF((MOD(NSTEP,200).EQ.0.OR.IFPRINT.EQ.1.OR.IFEND.EQ.1)
     &    .AND.IFBLEED.NE.0)  THEN
C
      HBLEDTOT = 0.0
      IF(KBLEED.NE.0) THEN
      DO 5507 J=2,JM
      DO 5507 I=1,IMM1
      HBLEDTOT  = HBLEDTOT + NBLADE(J)*(RFLUX(I,J,KM) - RFLUX(I,J,1))
 5507 CONTINUE
      ENDIF
C
      IF(IBLEED.NE.0) THEN
      DO 5508 J=2,JM
      DO 5508 K=1,KMM1
      HBLEDTOT = HBLEDTOT + NBLADE(J)*(TFLUX(IM,J,K)  - TFLUX(1,J,K))
 5508 CONTINUE
      ENDIF
C
      ENDIF
C
C     ADD ANY COOLANT ENERGY FLUXES THROUGH THE BLADE SURFACES.
C
      IF(NCOOLB.NE.0) THEN
      DO 5650 NC=1,NCOOLB
      DO 5656 J = JCBS(NC)+1,JCBE(NC)
      DO 5656 K = KCBS(NC),KCBE(NC)-1
      IF(IC(NC).EQ.1) SOURCE(1,J,K)    = SOURCE(1,J,K) -  HOCWLI1(J,K)
      IF(IC(NC).EQ.IM)SOURCE(IMM1,J,K) = SOURCE(IMM1,J,K)-HOCWLIM(J,K)
 5656 CONTINUE
 5650 CONTINUE
      ENDIF
C
C     ADD ANY COOLANT ENERGY FLOWS THROUGH THE HUB AND CASING.
C
C   Q3D
      IF(KM.EQ.2) GO TO 5657
C   END Q3D
C
      IF(NCOOLW.NE.0) THEN
      DO 5658 NC=1,NCOOLW
      DO 5659 J=JCWS(NC)+1,JCWE(NC)
      DO 5659 I=ICWS(NC),ICWE(NC)-1     
      IF(KC(NC).EQ.1)  SOURCE(I,J,1)    = SOURCE(I,J,1) -  HOCWLK1(I,J)
      IF(KC(NC).EQ.KM) SOURCE(I,J,KMM1) = SOURCE(I,J,KMM1)-HOCWLKM(I,J)
 5659 CONTINUE
 5658 CONTINUE
      ENDIF
C
 5657 CONTINUE
C
C      BALANCE THE FLUXES OF ENERGY ACROSS PERIODIC BOUNDARIES
C
      DO 900 J=2,JM
      IF(IND(J).EQ.1) GO TO 900
      DO 901 K=1,KMM1
      TFLUX(1,J,K)  = 0.5*(TFLUX(1,J,K)+TFLUX(IM,J,K))
      TFLUX(IM,J,K) = TFLUX(1,J,K)
  901 CONTINUE
  900 CONTINUE
C
C      SET SHROUD LEAKAGE ENERGY FLUXES
C
      IF(IFSHROUD.EQ.1) CALL SHROUDFLUX(SHROUDHO)
C
C*******************************************************************************
C*******************************************************************************   
C         CALL TSTEP TO UPDATE THE INTERNAL ENERGY AT ALL POINTS
C
            CALL TSTEP(ROE,DROE,3)   
C
C*******************************************************************************
C*******************************************************************************
C
C   NOW START TO UPDATE THE INFLOW AND OUTFLOW BOUNDARY CONDITIONS.
C
C*******************************************************************************
C*******************************************************************************
C
      IF(PLATE_LOSS.GT.0.001) THEN
           DO 5999 K=1,KM
           ROVMSQ(K) = 0.0
           DO 5999 I=1,IM
           ROVMSQ(K) = ROVMSQ(K) +
     &            ROVX(I,JM,K)*VX(I,JM,K) + ROVR(I,JM,K)*VR(I,JM,K)
 5999      CONTINUE
           ROVMSQ(K) = 0.5*ROVMSQ(K)/IM
      ENDIF
C
C     VARY THE FIXED HUB OR CASING PRESSURES IF USING THE THROTTLE EXIT CONDITION.
C
      PTHROTTLE = 0.0
      PDMID     = 0.0
      IF(NSTEP.GT.100.AND.THROTTLE_EXIT.GT.0.001) THEN
           FFLOW      = CHECK_FLOW(JM)/THROTTLE_MAS 
           PTHROTTLE  = THROTTLE_PRES*FFLOW*FFLOW
           PDOWN_TIP  = RFTHROTL1*PDOWN_TIP + RFTHROTL*PTHROTTLE
           PDOWN_HUB  = RFTHROTL1*PDOWN_HUB + RFTHROTL*PTHROTTLE
           PDMID      = PTHROTTLE - PD(KMID)
C   8/4/2017.  GRADUALLY RELAX THROTTLE_PRESS TO THE CURRENT VALUE OF PTHROTTLE.
           RF_P       = 0.02*RFTHROTL
           THROTTLE_PRES  = (1.-RF_P)*THROTTLE_PRES + RF_P*PTHROTTLE
           IF(MOD(NSTEP,50).EQ.0) THEN
           WRITE(6,*) 'THROTTLE SET ',THROTTLE_PRES,
     &                'PTHROTTLE= ',PTHROTTLE,' FLOW RATIO=', FFLOW
           END IF
      END IF
C
      PDIFF = 0.0
      IF(IPOUT.EQ.-1) PDIFF = PDOWN_TIP - PD(KM)
C
C************************************************************************
C************************************************************************
C   START A LOOP OVER ALL SPANWISE POINTS
      DO 6000 K=1,KM
C
C  FORM SOME PITCHWISE AVERAGED PRESSURES
      PAVG2   = 0.0
      PAVG3   = 0.0
      PAVGOUT = 0.0
      SUM     = 0.0
      DO 6010 I=1,IMM1
      PAVG2   = PAVG2   + 0.5*(P(I,2,K)+P(I+1,2,K))*FP(I)
      PAVG3   = PAVG3   + 0.5*(P(I,3,K)+P(I+1,3,K))*FP(I)
      PAVGOUT = PAVGOUT + 0.5*(P(I,JMM1,K)+P(I+1,JMM1,K))*FP(I)
      SUM     = SUM     + 0.5*FP(I)*
     &(ROVT(I,JM,K)*VT(I,JM,K)+ROVT(I+1,JM,K)*VT(I+1,JM,K))/R(JM,K)
 6010 CONTINUE
C
C***********************************************************************
C************************************************************************
C       RESET THE PRESSURE  "PD(K)" AT THE DOWNSTREAM (EXIT) BOUNDARY
C       USING RADIAL EQUILIBRIUM AND ASSUMING NO OTHER RADIAL ACCELERATION
C       IF  IPOUT = 0   or  = -1 .
C       OR MAINTAINING CONSTANT INPUT VALUES IF IPOUT = 1, OR = 3.
C
      IF(IPOUT.LE.0) THEN
           PD(1) = PDOWN_HUB
           IF(K.EQ.1)     GO TO 6017
           DP    = (SUM + SUMP)*0.5*(R(JM,K)-R(JM,K-1))
           PD(K) = 0.9*PD(K)+ 0.1*(PD(K-1) + DP)
 6017 CONTINUE
           SUMP  = SUM
           IF(IPOUT.EQ.-1) PD(K) = PD(K) + PDIFF
      END IF       
C
C     VARY THE FIXED EXIT PRESSURE WITH THE EXIT MASS FLOW RATE IF THROTTLE_EXIT GT 0.
C
      IF(IPOUT.GE.1.AND.THROTTLE_EXIT.GT.0.001) THEN
             PD(K) = PD(K) + RFTHROTL*PDMID
      END IF
C
C   USE A SPANWISE PRESSURE LOSS PROPORTIONAL TO ROVM**2 TO MAKE THE EXIT FLOW
C   MORE UNIFORM IF PLATE_LOSS > ZERO.
C
      IF(PLATE_LOSS.GT.0.001) PD(K) = PD(K)
     &  + 0.01*PLATE_LOSS*(ROVMSQ(K) - ROVMSQ(KMID))
C
C************************************************************************
C************************************************************************
C
C     START A LOOP OVER ALL PITCHWISE GRID POINTS
C
      DO 6030 I=1,IM
C
C         UPDATE FLOW CONDITIONS AT INLET BY EXTRAPOLATING THE PRESSURE TO THE
C         INLET BOUNDARY AND TAKING ISENTROPIC FLOW FROM THE INLET
C         STAGNATION CONDITIONS
C
      IF(ITIMST.GE.5) THEN
           PNEW     = PO_REF - (RO_REF - ROSUB(I,1,K))*DP_DRO
           P(I,1,K) = RFIN1*P(I,1,K) + RFIN*PNEW
      ELSE
           RO_STAG = PO1(K)/RGAS/TO1(K)
           RO_IN   = RO(I,1,K)
           IF(RO_IN.GT.0.9999*RO_STAG) RO_IN = 0.9999*RO_STAG
           PNEW    = PO1(K) * (RO_IN/RO_STAG)**GA
C
           IF(IN_PRESS.EQ.0) P(I,1,K) = RFIN1*P(I,1,K) + RFIN*PNEW
C
           IF(IABS(IN_PRESS).EQ.1) P(I,1,K) =
     &     RFIN1*P(I,1,K) + RFIN2*P(I,2,K)- RFIN*P(I,3,K)
C
           IF(IABS(IN_PRESS).EQ.3) P(I,1,K) = 
     &     RFIN1*P(I,1,K) + RFIN2*PAVG2 - RFIN*PAVG3
C
           IF(IABS(IN_PRESS).EQ.4) P(I,1,K) =
     &     RFIN1*P(I,1,K) + RFIN2*P(IMID,2,KMID) - RFIN*P(IMID,3,KMID)
C
      END IF
C
      IF(P(I,1,K).GT.PLIM(K))    P(I,1,K)  = 0.9*P(I,1,K) + 0.1*PLIM(K)
      IF(P(I,1,K).GT.0.99999*PO1(K)) P(I,1,K) = PO1(K)*0.99999
C
C  NEXT SET THE VELOCITIES AND TEMPERATURES AT INLET
C
      IF(ITIMST.NE.6) THEN
             IF(IFGAS.EQ.0) THEN
                  TSTATIC = TO1(K)*(P(I,1,K)*RPO1(K))**FGA
                  VABSQ   = 2*CP*(TO1(K) - TSTATIC)
             ELSE
                  PRAT    = P(I,1,K)*RPO1(K)
                  TRAT    = TRAT_FROM_PRAT(TO1(K),PRAT,FGAGAS,R_ALPHA,
     &                      BALPHA1,BALPHA2)
                  IF(TRAT.GT.1.0) TRAT = 0.999999
                  TSTATIC = TO1(K)*TRAT
                  HSTAG   = HFROMT(TO1(K),TREF,HREF,CP1,CP2,CP3) 
                  HSTAT   = HFROMT(TSTATIC,TREF,HREF,CP1,CP2,CP3)
                  IF(HSTAT.GT.HSTAG) HSTAT = 0.99999*HSTAG 
                  VABSQ   = 2.0*(HSTAG - HSTAT)
             END IF  
      ELSE
             VABSQ   = 2.0*(PO1(K) - P(I,1,K))/DENSTY
             TSTATIC = TO1(K) - 0.5*VABSQ/CP
      END IF
C
      T_STATIC(I,1,K) = TSTATIC
      VABS    = SQRT(VABSQ)         
C
      IF(IN_VTAN.EQ.0) THEN
           VMER      = VABS*BSCOS(K)
           VT(I,1,K) = VABS*BSSIN(K)
      END IF
C
      IF(IN_VTAN.EQ.1) THEN
           VMERSQ = VABSQ - VTIN(K)*VTIN(K)
           IF(VMERSQ.LT.0.01) VMERSQ = 0.01
           VMER      = SQRT(VMERSQ)
           VT(I,1,K) = VTIN(K)
      END IF
C
      IF(IN_VTAN.EQ.2) THEN
           WR        = UBLADE(1,K)
           VRELSQ    = VABSQ - WR*WR*BSCOS(K)*BSCOS(K)
           IF(VRELSQ.LT.0.01) VRELSQ = 0.01
           VREL      = SQRT(VRELSQ) - WR*BSSIN(K)
           VMER      = VREL*BSCOS(K)
           VT(I,1,K) = VREL*BSSIN(K) + WR
      END IF
C
C
      IF(IN_VR.GT.0) VX(I,1,K) = VMER*BRCOS(K)
      VXSQ = VMER*VMER - VR(I,2,K)*VR(I,2,K)
      IF(VXSQ.LT.0.0000001*TO1(K)*CP) VXSQ = 0.0000001*TO1(K)*CP
      IF(IN_VR.LE.0) VX(I,1,K) = SQRT(VXSQ)
      IF(IN_VR.GT.0) VR(I,1,K) = VMER*BRSIN(K)
      IF(IN_VR.LE.0) VR(I,1,K) = VR(I,2,K)
C
      IF(ITIMST.GE.5) THEN
C      
           IF(ITIMST.EQ.5) RO(I,1,K) = RF_PTRU*P(I,1,K)/RGAS/TSTATIC 
     &                               + RF_PTRU1*RO(I,1,K) 

           IF(ITIMST.EQ.6) RO(I,1,K) = DENSTY
C	   
           ROVT(I,1,K)  = RO(I,1,K)*VT(I,1,K)
           ROVR(I,1,K)  = RO(I,1,K)*VR(I,1,K)
           RORVT(I,1,K) = ROVT(I,1,K)*R(1,K)
           ROVX(I,1,K)  = RO(I,1,K)*VX(I,1,K)
	   ROE(I,1,K)   = ROSUB(I,1,K)*(CV*TSTATIC + 0.5*VABSQ)
	   HO(I,1,K)    = CP*TO1(K)
      ELSE
           RO(I,1,K)    = P(I,1,K)/(RGAS*TSTATIC)
           ROVT(I,1,K)  = RO(I,1,K)*VT(I,1,K)
           ROVR(I,1,K)  = RO(I,1,K)*VR(I,1,K)
           RORVT(I,1,K) = ROVT(I,1,K)*R(1,K)
           ROVX(I,1,K)  = RO(I,1,K)*VX(I,1,K)
C
           IF(IFGAS.EQ.0) THEN
                ROE(I,1,K)   = P(I,1,K)/GA1 + RO(I,1,K)*0.5*VABSQ
   	        HO(I,1,K)    = CP*TO1(K)
           ELSE
                ROE(I,1,K)   = RO(I,1,K)*(HSTAT-RGAS*TSTATIC+0.5*VABSQ)
                HO(I,1,K)    = HSTAG
           END IF
C
      END IF
C
C******************************************************************************
C******************************************************************************
C  UPDATE THE TURBULENT VISCOSITY AT INLET IF USING THE SA MODEL.
C
      IF(ILOS.GE.200) THEN
           IF(FSTURB(1).LT.1.0) THEN 
           TRANSVISIN  = VISC_LAM(I,2,K)*4.35*FSTURB(1)**0.25
           ELSE
           TRANSVISIN  = VISC_LAM(I,2,K)*(3.5 + FSTURB(1)*0.85)
           END IF
           TRANS_DYN_VIS(I,1,K) = TRANSVISIN
           TRANS_DYN_VIS(I,2,K) = TRANS_DYN_VIS(I,1,K)
           TRANS_KVIS(I,1,K)    = TRANS_DYN_VIS(I,1,K)/RO(IMID,1,KMID)
           TRANS_KVIS(I,2,K)    = TRANS_KVIS(I,1,K)
      END IF
C
C     END OF SETTING THE INLET BOUNDARY CONDITIONS.
C
C*********************************************************************
C*********************************************************************
C 
C     THE NEXT SECTION ALLOWS THE EXIT PRESSURE TO VARY ALONG THE PITCH BY
C     EXTRAPOLATING THE VARIATION FROM UPSTREAM. THE STAGNATION ENTHALPY 
C     IS ALSO EXTRAPOLATED
C
C
       EKE  =.5*(VR(I,JM,K)*VR(I,JM,K) + VT(I,JM,K)*VT(I,JM,K) 
     &      + VX(I,JM,K)*VX(I,JM,K))  
C
      IF(ITIMST.GE.5) THEN
            P(I,JM,K)    = PD(K)  + FEXTRAP*(P(I,JMM1,K) - PAVGOUT)     
            ROSUB(I,JM,K)= RO_REF - (PO_REF - P(I,JM,K))/DP_DRO
            HO(I,JM,K)   = 2.0*HO(I,JMM1,K) - HO(I,JM-2,K)
            TEXIT        = (HO(I,JM,K) - EKE)/CP
            ROE(I,JM,K)  = ROSUB(I,JM,K)*(CV*TEXIT + EKE)  
	    RO(I,JM,K)   = P(I,JM,K)/RGAS/TEXIT 
            IF(ITIMST.EQ.6) RO(I,JM,K) = DENSTY   
      ELSE
C
            P(I,JM,K)    = PD(K)  + FEXTRAP*(P(I,JMM1,K) - PAVGOUT)
            HO(I,JM,K)   = 2.0*HO(I,JMM1,K) - HO(I,JM-2,K)
C
            IF(IFGAS.EQ.0) THEN
                  TEXIT        = (HO(I,JM,K) - EKE)/CP
	          RO(I,JM,K)   = P(I,JM,K)/RGAS/TEXIT
                  ROE(I,JM,K)  = RO(I,JM,K)*(CV*TEXIT + EKE)
            ELSE
                  HEXIT        = HO(I,JM,K) - EKE
                  TEXIT        = TFROMH(HEXIT,TREF,HREF,HT1,HT2,HT3,HT4)
	          RO(I,JM,K)   = P(I,JM,K)/RGAS/TEXIT
                  ROE(I,JM,K)  = RO(I,JM,K)*(HEXIT - RGAS*TEXIT + EKE)               
            END IF
C
      END IF
C
C   END OF THE " DO I "  LOOP
 6030 CONTINUE
C
C     END OF THE  " DO  K " LOOP
 6000 CONTINUE
C
C******************************************************************************
C******************************************************************************
C    ALL INLET AND EXIT BOUNDARY CONDITIONS ARE NOW UPDATED
C
C******************************************************************************
C******************************************************************************
C       CALL SETFLO TO ADJUST THE ROV'S TO TRY TO ACHIEVE A SPECIFIED MASS
C       FLOW RATE IF IN_FLOW=3. OR TO RELAX TO AN AVERAGE FLOW IF IN_FLOW=2
C       IN_FLOW=2  OFTEN GIVES IMPROVED CONVERGENCE.
C
      IF(IN_FLOW.NE.0) CALL SETFLO
C
C    FORM PEFF TO MINIMISE ROUNDING ERRORS.
C    ALSO ALLOW DOWNWINDING OF PRESSURE TO SMEAR SHOCKS USING  FP_DOWN .
C
      DO 4321 K=1,KM
      DO 4321 J=1,JM
      JP1 = J+1
      IF(J.EQ.JM) JP1 = JM
      DO 4321 I=1,IM
      PEFF(I,J,K) = P(I,J,K) - PREFF + F_PDOWN*(P(I,JP1,K) - P(I,J,K))
 4321 CONTINUE
C
C*******************************************************************************
C   TFLOW ADDITION
C*******************************************************************************
      IF(IM.EQ.2) THEN
C
C     RESOLVE THE BLADE FORCE PERPENDICULAR TO THE STREAMLINE SO IT DOES NOT
C     GENERATE ANY LOSS.
C
      DO 6050 K=1,KMM1
      KP1 = K+1
      DO 6050 J=2,JM
      JM1 = J-1
      P_DIFF = 0.5*PBLADE(J,K)
      PEFF(1,J,K)    = PEFF(1,J,K)     - P_DIFF
      PEFF(1,J,KP1)  = PEFF(1,J,KP1)   - P_DIFF
      PEFF(1,JM1,K)  = PEFF(1,JM1,K)   - P_DIFF
      PEFF(1,JM1,KP1)= PEFF(1,JM1,KP1) - P_DIFF
      PEFF(2,J,K)    = PEFF(2,J,K)     + P_DIFF
      PEFF(2,J,KP1)  = PEFF(2,J,KP1)   + P_DIFF
      PEFF(2,JM1,K)  = PEFF(2,JM1,K)   + P_DIFF
      PEFF(2,JM1,KP1)= PEFF(2,JM1,KP1) + P_DIFF
 6050 CONTINUE
C
      DO 6051 J=2,JM

      DO 6052 K=1,KMM1
      AVGP1 =(PEFF(1,J,K)+PEFF(1,J,K+1)+PEFF(1,J-1,K)+PEFF(1,J-1,K+1))/4
      AVGP2 =(PEFF(2,J,K)+PEFF(2,J,K+1)+PEFF(2,J-1,K)+PEFF(2,J-1,K+1))/4
      ABXAVG = (ABX(1,J,K) + ABX(2,J,K))/2
      ABRAVG = (ABR(1,J,K) + ABR(2,J,K))/2
      XFOR   = (AVGP1 - AVGP2)*ABXAVG
      RFOR   = (AVGP1 - AVGP2)*ABRAVG
      TFOR   = (AVGP1 - AVGP2)*ABT(J,K)
C      
      VXAVG = VX(1,J,K)+VX(1,J,K+1)+VX(1,J-1,K)+VX(1,J-1,K+1)
      VRAVG = VR(1,J,K)+VR(1,J,K+1)+VR(1,J-1,K)+VR(1,J-1,K+1)
      WTAVG = WT(1,J,K)+WT(1,J,K+1)+WT(1,J-1,K)+WT(1,J-1,K+1)
      WAVG = SQRT(VXAVG*VXAVG + VRAVG*VRAVG + WTAVG*WTAVG)
      SFOR = (XFOR*VXAVG + RFOR*VRAVG + TFOR*WTAVG)/WAVG
C        
      IF(IND(J).EQ.0) GO TO 6052
C
      BFORCE_X(J,K) = SFOR*VXAVG/WAVG
      BFORCE_T(J,K) = SFOR*WTAVG/WAVG*RAVG_CELL(J,K)
      BFORCE_R(J,K) = SFOR*VRAVG/WAVG
      BFORCE_Q(J,K) = BFORCE_T(J,K)*WRAD(J)
C
 6052 CONTINUE
 6051 CONTINUE
C      
      END IF
C*******************************************************************************
C   END OF TFLOW ADDITION
C*******************************************************************************
C
C*******************************************************************************
C       CALCULATE FLUXES FOR THE AXIAL-MOMENTUM EQUATION.
C*******************************************************************************
C******************************************************************************
      DO 6100 K=1,KMM1
      DO 6100 J=1,JM
      DO 6100 I=1,IMM1
      SOURCE(I,J,K) = XFORCE(I,J,K)
      AVGP  = PEFF(I,J,K)+PEFF(I,J,K+1)+PEFF(I+1,J,K+1)+PEFF(I+1,J,K)
      AVGVX = VX(I,J,K)+VX(I,J,K+1)+VX(I+1,J,K+1)+VX(I+1,J,K)
      XFLUX(I,J,K) = 0.25*(AVGP*AQX(I,J,K) + FLOWX(I,J,K)*AVGVX)
 6100 CONTINUE
C
C  TFLOW ADDITION
      IF(IM.EQ.2) THEN
      DO 6101 K=1,KMM1
      DO 6101 J=2,JM
      SOURCE(1,J,K) = SOURCE(1,J,K) + BFORCE_X(J,K)
 6101 CONTINUE
      END IF
C   END TFLOW ADDITION
C
      DO 6110 K=1,KM
      DO 6110 J=2,JM
      DO 6110 I=1,IMM1
      AVGP  = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I+1,J-1,K)+PEFF(I+1,J,K)
      AVGVX = VX(I,J,K)+VX(I,J-1,K)+VX(I+1,J-1,K)+VX(I+1,J,K)
      RFLUX(I,J,K) = 0.25*(AVGP*ASX(I,J,K) + FLOWR(I,J,K)*AVGVX)
 6110 CONTINUE
C
C
      DO 6120 K=1,KMM1
      DO 6120 J=2,JM
      DO 6120 I=1,IM
      AVGVX  = VX(I,J,K)+VX(I,J-1,K)+VX(I,J-1,K+1)+VX(I,J,K+1)
      AVGP   = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I,J-1,K+1)+PEFF(I,J,K+1)
      TFLUX(I,J,K) = 0.25*(FLOWT(I,J,K)*AVGVX + ABX(I,J,K)*AVGP)
 6120 CONTINUE
C
C     ADD COOLANT AXIAL MOMENTUM FLUXES THROUGH BLADES
C
      IF(NCOOLB.NE.0) THEN
      DO 6650 NC=1,NCOOLB
      DO 6656 J = JCBS(NC)+1,JCBE(NC)
      DO 6656 K = KCBS(NC),KCBE(NC)-1
      IF(IC(NC).EQ.1) SOURCE(1,J,K)   = SOURCE(1,J,K)    - VXCWLI1(J,K)
      IF(IC(NC).EQ.IM)SOURCE(IMM1,J,K)= SOURCE(IMM1,J,K) - VXCWLIM(J,K)
 6656 CONTINUE
 6650 CONTINUE
      ENDIF
C
C  Q3D
      IF(KM.EQ.2) GO TO 6657
C  END Q3D
C     ADD COOLANT AXIAL MOMENTUM FLOWS THROUGH THE HUB AND CASING
C
      IF(NCOOLW.NE.0) THEN
      DO 6658 NC=1,NCOOLW
      DO 6659 J=JCWS(NC)+1,JCWE(NC)
      DO 6659 I=ICWS(NC),ICWE(NC)-1
      IF(KC(NC).EQ.1)  SOURCE(I,J,1)   = SOURCE(I,J,1)    - VXCWLK1(I,J)
      IF(KC(NC).EQ.KM) SOURCE(I,J,KMM1)= SOURCE(I,J,KMM1) - VXCWLKM(I,J)
 6659 CONTINUE
 6658 CONTINUE
      ENDIF
C
 6657 CONTINUE
C
C          BALANCE FLUXES ON PERIODIC BOUNDARIES FOR ROVX
C
      DO 6026 J=2,JM
      IF(IND(J).EQ.1) GO TO 6026
      DO 6027 K=1,KMM1
      AVGPCUSP =  PEFF(1,J,K)     + PEFF(1,J-1,K) + PEFF(1,J,K+1)
     &          + PEFF(1,J-1,K+1) + PEFF(IM,J,K)  + PEFF(IM,J-1,K)
     &          + PEFF(IM,J,K+1)  + PEFF(IM,J-1,K+1)
      AVGVXCUSP = VX(1,J,K)       + VX(1,J-1,K) + VX(1,J,K+1)
     &          + VX(1,J-1,K+1)   + VX(IM,J,K)  + VX(IM,J-1,K)
     &          + VX(IM,J,K+1)    + VX(IM,J-1,K+1)      
      TFLUX(1,J,K) =0.125*(FLOWT(1,J,K)*AVGVXCUSP +AVGPCUSP*ABX(1,J,K))
      TFLUX(IM,J,K)=0.125*(FLOWT(IM,J,K)*AVGVXCUSP+AVGPCUSP*ABX(IM,J,K))
 6027 CONTINUE
 6026 CONTINUE
C
C      SET SHROUD LEAKAGE AXIAL MOMENTUM FLUXES
C
      IF(IFSHROUD.EQ.1)  CALL SHROUDFLUX(SHROUDVX)
C
C******************************************************************************
C******************************************************************************
C      CALL SUBROUTINE TSTEP TO UPDATE THE AXIAL-MOMENTUM
C
      CALL TSTEP(ROVX,DROVX,4)
C
C******************************************************************************
C******************************************************************************
C       CALCULATE FLUXES FOR MOMENT OF MOMENTUM EQUATION.
C******************************************************************************
C
      DO 6200 K=1,KMM1
      DO 6200 J=1,JM
      AVGR  =  0.125*(R(J,K)+R(J,K+1))
      DO 6200 I=1,IMM1
      AVGVT = VT(I,J,K)+VT(I,J,K+1)+VT(I+1,J,K)+VT(I+1,J,K+1)
      SOURCE(I,J,K) = TFORCE(I,J,K)
      XFLUX(I,J,K)  = FLOWX(I,J,K)*AVGR*AVGVT
 6200 CONTINUE
C
C  TFLOW ADDITION
      IF(IM.EQ.2) THEN
      DO 6201 K=1,KMM1
      DO 6201 J=2,JM
      SOURCE(1,J,K) = SOURCE(1,J,K) + BFORCE_T(J,K)
 6201 CONTINUE
      END IF
C   END TFLOW ADDITION
C
      DO 6210 K=1,KMM1
      DO 6210 J=2,JM
      AVGR  = 0.25*RAVG_CELL(J,K)
      DO 6210 I=1,IM
      AVGVT = VT(I,J,K)+VT(I,J-1,K)+VT(I,J-1,K+1)+VT(I,J,K+1)
      AVGP  = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I,J-1,K+1)+PEFF(I,J,K+1)
      TFLUX(I,J,K) = (FLOWT(I,J,K)*AVGVT  +  ABT(J,K)*AVGP)*AVGR
 6210 CONTINUE
C
      DO 6230 K=1,KM
      DO 6230 J=2,JM
      AVGR=(R(J,K)+R(J-1,K))*0.125
      DO 6230 I=1,IMM1
      AVGVT = VT(I,J,K)+VT(I,J-1,K)+VT(I+1,J-1,K)+VT(I+1,J,K)
      RFLUX(I,J,K) = FLOWR(I,J,K)*AVGVT*AVGR
 6230 CONTINUE
C
C      BALANCE FLUXES OF ANGULAR MOMENTUM ACROSS PERIODIC BOUNDARIES
C
      DO 902 J=2,JM
      IF(IND(J).EQ.1) GO TO 902
      DO 903 K=1,KMM1
      TFLUX(1,J,K)  = 0.5*(TFLUX(1,J,K) + TFLUX(IM,J,K))
      TFLUX(IM,J,K) = TFLUX(1,J,K)
  903 CONTINUE
  902 CONTINUE
C
C     ADD COOLANT TANGENTIAL MOMENTUM FLUXES THROUGH BLADES
C
      IF(NCOOLB.NE.0) THEN
      DO 7750 NC = 1,NCOOLB
      DO 7756 J  = JCBS(NC)+1,JCBE(NC)
      DO 7756 K  = KCBS(NC),KCBE(NC)-1
      IF(IC(NC).EQ.1)  SOURCE(1,J,K)  = SOURCE(1,J,K)    - RVTCWLI1(J,K)
      IF(IC(NC).EQ.IM) SOURCE(IMM1,J,K)=SOURCE(IMM1,J,K) - RVTCWLIM(J,K)
 7756 CONTINUE
 7750 CONTINUE
      ENDIF
C
C  Q3D
      IF(KM.EQ.2) GO TO 7757
C  END Q3D
C
C     ADD COOLANT TANGENTIAL MOMENTUM FLOWS THROUGH HUB AND CASING
C
      IF(NCOOLW.NE.0) THEN
      DO 7758 NC = 1,NCOOLW
      DO 7759 J  = JCWS(NC)+1,JCWE(NC)
      DO 7759 I=ICWS(NC),ICWE(NC)-1
      IF(KC(NC).EQ.1)  SOURCE(I,J,1)  = SOURCE(I,J,1)    - RVTCWLK1(I,J)
      IF(KC(NC).EQ.KM) SOURCE(I,J,KMM1)=SOURCE(I,J,KMM1) - RVTCWLKM(I,J)
 7759 CONTINUE
 7758 CONTINUE
      ENDIF
C
 7757 CONTINUE
C
C     SET SHROUD LEAKAGE ANGULAR MOMENTUM FLUXES
C
      IF(IFSHROUD.EQ.1) CALL SHROUDFLUX(SHROUDRVT)
C
C*******************************************************************************
C******************************************************************************
C       CALL SUBROUTINE TSTEP TO UPDATE THE MOMENT OF MOMENTUM
C
      CALL TSTEP(RORVT,DRORVT,5)
C
C*******************************************************************************
C******************************************************************************
C       CALCULATE FLUXES FOR THE RADIAL MOMENTUM EQUATION
C******************************************************************************
C
      DO 6300 K=1,KMM1
      DO 6300 J=1,JM
      DO 6300 I=1,IMM1
      AVGVR = VR(I,J,K)+VR(I,J,K+1)+VR(I+1,J,K)+VR(I+1,J,K+1)
      AVGP  = PEFF(I,J,K)+PEFF(I,J,K+1)+PEFF(I+1,J,K+1)+PEFF(I+1,J,K)
      XFLUX(I,J,K) = 0.25*(FLOWX(I,J,K)*AVGVR + AQR(I,J,K)*AVGP)
 6300 CONTINUE

      DO 6310 K=1,KMM1
      DO 6310 J=2,JM
      DO 6310 I=1,IM
      AVGVR = VR(I,J,K)+VR(I,J-1,K)+VR(I,J-1,K+1)+VR(I,J,K+1)
      AVGP  = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I,J-1,K+1)+PEFF(I,J,K+1)
      TFLUX(I,J,K) =.25*(FLOWT(I,J,K)*AVGVR + AVGP*ABR(I,J,K))
 6310 CONTINUE
C
      DO 6320 K=1,KM
      DO 6320 J=2,JM
      DO 6320 I=1,IMM1
      AVGVR = VR(I,J,K)+VR(I,J-1,K)+VR(I+1,J-1,K)+VR(I+1,J,K)
      AVGP  = PEFF(I,J,K)+PEFF(I,J-1,K)+PEFF(I+1,J-1,K)+PEFF(I+1,J,K)
      RFLUX(I,J,K) = 0.25*(FLOWR(I,J,K)*AVGVR + ASR(I,J,K)*AVGP)
 6320 CONTINUE
C
C     CALCULATE THE SOURCE TERM. THIS IS DUE TO THE CENTRIFUGAL FORCE AND ANY 
C     IMBALANCE IN THE RADIAL PROJECTED AREAS.
C
      DO 6330 K=1,KMM1
      DO 6330 J=2,JM
      DO 6330 I=1,IMM1
      AVGRVT = ROVT(I,J,K)+ROVT(I,J,K+1)+ROVT(I,J-1,K)+ROVT(I,J-1,K+1)
     & + ROVT(I+1,J,K)+ROVT(I+1,J,K+1)+ROVT(I+1,J-1,K)+ROVT(I+1,J-1,K+1)
      AVGVT = 0.125*(VT(I,J,K)+VT(I,J-1,K)+VT(I,J-1,K+1)+VT(I,J,K+1)
     & + VT(I+1,J,K)+VT(I+1,J-1,K)+VT(I+1,J-1,K+1)+VT(I+1,J,K+1))
      AVGP = PEFF(I,J,K)+ PEFF(I,J-1,K)+ PEFF(I,J-1,K+1)+ PEFF(I,J,K+1)
     & + PEFF(I+1,J,K)+PEFF(I+1,J-1,K)+PEFF(I+1,J-1,K+1)+PEFF(I+1,J,K+1)
      SOURCE(I,J,K) = RFORCE(I,J,K) - (AVGRVT*AVGVT + AVGP)*VOLOR(I,J,K)
 6330 CONTINUE
C
C  TFLOW ADDITION
      IF(IM.EQ.2) THEN
      DO 6331 K=1,KMM1
      DO 6331 J=2,JM
      SOURCE(1,J,K) = SOURCE(1,J,K) + BFORCE_R(J,K)
 6331 CONTINUE
      END IF
C   END TFLOW ADDITION
C
C     ADD COOLANT RADIAL MOMENTUM FLUXES THROUGH THE BLADE SURFACES.
C
      IF(NCOOLB.NE.0) THEN
      DO 9650 NC=1,NCOOLB
      DO 9656 J = JCBS(NC)+1,JCBE(NC)
      DO 9656 K = KCBS(NC),KCBE(NC)-1
      IF(IC(NC).EQ.1)  SOURCE(1,J,K)   = SOURCE(1,J,K)    - VRCWLI1(J,K)
      IF(IC(NC).EQ.IM) SOURCE(IMM1,J,K)= SOURCE(IMM1,J,K) - VRCWLIM(J,K)
 9656 CONTINUE
 9650 CONTINUE
      ENDIF
C
C     ADD COOLANT RADIAL MOMENTUM FLOWS THROUGH THE HUB AND CASING
C
C  Q3D
      IF(KM.EQ.2) GO TO 8657
C  END Q3D
C
      IF(NCOOLW.NE.0) THEN
      DO 8658 NC=1,NCOOLW
      DO 8659 J=JCWS(NC)+1,JCWE(NC)
      DO 8659 I=ICWS(NC),ICWE(NC)-1
      IF(KC(NC).EQ.1)  SOURCE(I,J,1)   = SOURCE(I,J,1)    - VRCWLK1(I,J)
      IF(KC(NC).EQ.KM) SOURCE(I,J,KMM1)= SOURCE(I,J,KMM1) - VRCWLKM(I,J)
 8659 CONTINUE
 8658 CONTINUE
      ENDIF
C
 8657 CONTINUE
C
C          BALANCE FLUXES ON PERIODIC BOUNDARIES FOR ROVR
C
      DO 6260 J=2,JM
      IF(IND(J).EQ.1) GO TO 6260
      DO 6261 K=1,KMM1
      AVGP1  =  0.25*(PEFF(1,J,K)+PEFF(1,J-1,K)+PEFF(1,J,K+1)
     &              + PEFF(1,J-1,K+1))
      AVGPM  =  0.25*(PEFF(IM,J,K)+PEFF(IM,J-1,K)+PEFF(IM,J,K+1)
     &              + PEFF(IM,J-1,K+1))
      T1 = TFLUX(1,J,K) - AVGP1*ABR(1,J,K)
      T2 = TFLUX(IM,J,K)- AVGPM*ABR(IM,J,K)
      TFLUX(1,J,K)   = 0.5*((T1+T2) + (AVGP1+AVGPM)*ABR(1,J,K))
      TFLUX(IM,J,K)  = 0.5*((T1+T2) + (AVGP1+AVGPM)*ABR(IM,J,K))
 6261 CONTINUE
 6260 CONTINUE
C
C      SET SHROUD LEAKAGE RADIAL MOMENTUM FLUXES
C
      IF(IFSHROUD.EQ.1) CALL SHROUDFLUX(SHROUDVR)
C
C******************************************************************************
C******************************************************************************
C       CALL SUBROUTINE TSTEP TO UPDATE THE RADIAL-MOMENTUM
C
      CALL TSTEP(ROVR,DROVR,6)
C
C******************************************************************************
C******************************************************************************
C      ALL THE CONSERVATION EQUATIONS HAVE NOW BEEN UPDATED
C******************************************************************************
C******************************************************************************
C
C      CHECK FOR ANY VERY STEEP STREAMWISE DENSITY GRADIENTS OR ANY VERY LOW DENSITIES..
C
      DO 7200 K=1,KM
      DO 7200 J=3,JM-3
      DO 7200 I=1,IM
      ROAVG = 0.25*(RO(I,J-2,K)+RO(I,J-1,K)+RO(I,J+1,K)+RO(I,J+2,K))
      IF(RO(I,J,K).LT.0.6*ROAVG) RO(I,J,K) = 0.5*(RO(I,J,K)+0.6*ROAVG)
      IF(RO(I,J,K).GT.1.5*ROAVG) RO(I,J,K) = 0.5*(RO(I,J,K)+1.5*ROAVG)
      IF(RO(I,J,K).LT.ROLIM) RO(I,J,K) = 0.5*(RO(I,J,K) + ROLIM)
 7200 CONTINUE
C
C******************************************************************************
C  CHECK FOR ANY VERY HIGH MACH NUMBERS - LIMIT = MACHLIM
C
      GM1      = GA - 1.0
      TRATIO   = 1.0 + MACHLIM*MACHLIM*0.5*GM1
      RFGM1    = 1.0/GM1
      RORATIO  = TRATIO**RFGM1
      DO 7250 NR = 1,NROWS
      JREF   = JLE(NR)
      TREF   = T_STATIC(IMID,JREF,KMID)
      WREFSQ = VX(IMID,JREF,KMID)*VX(IMID,JREF,KMID) +
     &         VR(IMID,JREF,KMID)*VR(IMID,JREF,KMID) +
     &         WT(IMID,JREF,KMID)*WT(IMID,JREF,KMID)
      TOREL_REF   = TREF + 0.5*WREFSQ/CP
      TLIMIT      = TOREL_REF/TRATIO
      IF(TLIMIT.GT.TCOOL_MIN) TLIMIT = TCOOL_MIN
      VLIMIT(NR)  = SQRT(2*CP*(TOREL_REF - TLIMIT) )
      ROSTAG_REL  = RO(IMID,JREF,KMID) *(TOREL_REF/TREF)**RFGM1
      ROLIMIT(NR) = ROSTAG_REL/RORATIO
 7250 CONTINUE
C
C    CHECK AND APPLY LIMITS TO BOTH THE RELATIVE MACH NUMBER AND RELATIVE MASS FLUX .
C
      DO 7100 K=1,KM
      DO 7100 J=1,JM
      NR       = NROW(J)
      V_LIMIT  = VLIMIT(NR)
      RO_LIMIT = ROLIMIT(NR)
      DO 7100 I=1,IM
      IF(RO(I,J,k).LT.RO_LIMIT) RO(I,J,K) = RO_LIMIT
      RONOW       = RO(I,J,K)
      WRELSQ      = VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K)
     &            + WT(I,J,K)*WT(I,J,K)
      WREL        = SQRT(WRELSQ)
      FACSAFE     = V_LIMIT/WREL
      IF(FACSAFE.LT.1.0) THEN
           ROVX(I,J,K) = RONOW*VX(I,J,K) * FACSAFE
           ROVR(I,J,K) = RONOW*VR(I,J,K) * FACSAFE
           ROVTREL     = RONOW*WT(I,J,K) * FACSAFE
           ROWT(I,J,K) = ROVTREL
           ROVT(I,J,K) = ROVTREL + RONOW*UBLADE(J,K) 
           RORVT(I,J,K)= ROVT(I,J,K)*R(J,K)
      END IF
C
 7100 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     CALL SUBROUTINE NEW_MIXPLAN TO TRANSFER THE FLOW ACROSS THE MIXING PLANES
C
      IF(IM.GT.2) THEN
           IF(IFMIX.GT.0.AND.RFMIX.GT.1.0E-3) CALL NEW_MIXPLAN
      END IF
C
C******************************************************************************
C******************************************************************************
C      FORCE THE TRAILING EDGE SEPARATION POINTS IF  "IF_CUSP" = 2
C******************************************************************************
C******************************************************************************
C
           DO 7355 NR = 1,NROWS
C
       IF(IF_CUSP(NR).EQ.2) THEN
                JTEDGE   = JTE(NR)
                JSEP_I1  = JTEDGE - NUP_I1(NR)
                JSEP_IM  = JTEDGE - NUP_IM(NR)
                JTHIK    = MAX0(JSEP_I1,JSEP_IM) 
                FDRAG    = SEP_DRAG(NR)
                NWKE     = N_WAKE(NR)               
C
           DO 7350 K=1,KM
           THIK_LIM = SEP_THIK(NR)*RT_THICK(JTHIK,K)
           SSSLOPE = (RTHETA(1,JSEP_I1+1,K) - RTHETA(1,JSEP_I1-1,K))
     &              /(X(JSEP_I1+1,K) - X(JSEP_I1-1,K))
           PSSLOPE = (RTHETA(IM,JSEP_IM+1,K) - RTHETA(IM,JSEP_IM-1,K))
     &              /(X(JSEP_IM+1,K) - X(JSEP_IM-1,K))
C
           DO 7351 J = JSEP_I1,JTEDGE+NWKE
           RTEXTRAP = RTHETA(1,JSEP_I1,K)+ SSSLOPE*(X(J,K)-X(JSEP_I1,K))
           DO 7353 I = 1,6
           IF((RTEXTRAP-RTHETA(I,J,K)).GT.THIK_LIM) THEN
                ROVX(I,J,K)  = FDRAG*ROVX(I,J,K)
                ROVR(I,J,K)  = FDRAG*ROVR(I,J,K)
                RO_VBLADE  = RO(I,J,K)*UBLADE(J,K)
                RO_WT      = ROVT(I,J,K) - RO_VBLADE
                RO_VT      = FDRAG*RO_WT + RO_VBLADE                 
                RORVT(I,J,K) = RO_VT*R(J,K)
           END IF
 7353      CONTINUE
 7351      CONTINUE
C
           DO 7352 J = JSEP_IM,JTEDGE+NWKE
           RTEXTRAP = RTHETA(IM,JSEP_IM,K)+PSSLOPE*(X(J,K)-X(JSEP_IM,K))
           DO 7354 I = IM-6,IM
           IF((RTHETA(I,J,K)-RTEXTRAP).GT.THIK_LIM) THEN
                ROVX(I,J,K)  = FDRAG*ROVX(I,J,K)
                ROVR(I,J,K)  = FDRAG*ROVR(I,J,K)
                RO_VBLADE  = RO(I,J,K)*UBLADE(J,K)
                RO_WT      = ROVT(I,J,K) - RO_VBLADE
                RO_VT      = FDRAG*RO_WT + RO_VBLADE                 
                RORVT(I,J,K) = RO_VT*R(J,K)
           END IF
 7354      CONTINUE
 7352      CONTINUE
C
 7350      CONTINUE
C
       END IF
C
 7355      CONTINUE
C
C******************************************************************************
C*********************  END OF THE MAIN TIME STEPPING LOOP ********************
C
C     GO TO 8000 TO CHECK THE CONVERGENCE EVERY 5 TIME STEPS.
C
      IF(MOD(NSTEP,5).EQ.0.OR.(IFEND.EQ.1)) GO TO 8000
C
C******************************************************************************
C******************************************************************************
C      RETURN TO THE START OF THE MAIN LOOP FOR THE NEXT TIMESTEP
C
      GO TO 5000
C
C****************************************************************************
C      THE REMAINDER OF THIS SUB PROGRAM IS ONLY EXECUTED EVERY 5 STEPS
C****************************************************************************
C
 8000 CONTINUE
C
C       EVERY FIVE TIME STEPS CHECK CONVERGENCE AND PRINT OUT SUMMARY TO UNITS 4 AND 6
C
C        EVALUATE THE INLET AND LOCAL MASS FLOWS.
C
C        BLADE_FLOW  is the total mass flow through the blade passages,
C        including and cooling, leakage  or bleed flows.
C
C        CHECK_FLOW  Is the local mass flow excluding cooling flows, leakage and bleed flows.
C        This should be conserved and is used as a check for global continuity.  
C
      DO 5675 J=1,JM
      SUMAS = 0
      NB = NBLADE(J)
      DO 5665 I=1,IMM1
      DO 5665 K=1,KMM1
      SUMAS = SUMAS  - FLOWX(I,J,K)
 5665 CONTINUE
      BLADE_FLOW(J) =  SUMAS*NB
      CHECK_FLOW(J)  = (SUMAS + SHRDFLOW(J))*NB+ SUMBLEED(J)- SUMCWL(J)
 5675 CONTINUE
      ECONT    = 0.0
      DO 5677 J=1,JM
      RATIO = CHECK_FLOW(J)/CHECK_FLOW(1)
      EMASS = ABS(1.-RATIO)
      IF(EMASS.GT.ECONT) THEN
          ECONT = EMASS
          JCONT = J
      END IF
 5677 CONTINUE
      FLOWRAT =  BLADE_FLOW(1)
C
C******************************************************************************
C******************************************************************************
C     CALCULATE THE MAXIMUM PERCENTAGE CHANGE IN MERIDIONAL VELOCITY AND SAVE
C     IT AS STORE(I,J,K) FOR THE CONVERGENCE CHECK.
C     THIS IS ONLY DONE EVERY 5 ITERATIONS.
C
      EMAX=0.0
      EAVG=0.0
      VREF = SQRT(VRMS*RIJKM)
      RVEF = 100./VREF
      DO 8500 K=1,KM
      DO 8510 J=2,JM
      XD = DX(J,K)
      RD = DR(J,K)
      SD = DS(J,K)
      DO 8520 I=1,IM
      VM_START = (XD*VX(I,J,K)  + RD*VR(I,J,K))/SD
      VM_END   = (XD*ROVX(I,J,K)
     &         +  RD*ROVR(I,J,K))/SD/RO(I,J,K)
      DVMER  = (VM_END - VM_START)*RVEF
      STORE(I,J,K) = DVMER
      EAVG   = EAVG + ABS(DVMER)
      IF(ABS(DVMER).GT.ABS(EMAX)) THEN
           EMAX = DVMER
           IMAX = I
           JMAX = J
           KMAX = K
      ENDIF
 8520 CONTINUE
 8510 CONTINUE
 8500 CONTINUE
C
      EAVG=EAVG/(IM*JMM1*KM)
C
C******************************************************************************
C******************************************************************************
C    CALCULATE THE PITCHWISE AVERAGE VALUES AT EXIT FOR USE IN SETTING UP A REPEATING STAGE CONDITION.
C
      IF(IF_REPEAT.GT.0)     THEN
      IF(MOD(N,NINMOD).EQ.0) THEN
C
      CALL MIX_BCONDS(1)
C
      CALL NEWBCONDS                 
C
      END IF
      END IF
C
C******************************************************************************
C******************************************************************************
C    CALL OUTPUT TO PRINT OUT MAIN PRINTED OUTPUT IF REQUESTED, OR IF CONVERGED.
C    CALL THE LOSS ROUTINES TO UPDATE THE VISCOUS FORCES BEFORE PRINTING OUT.
C
      IF(IFPRINT.EQ.1.OR.IFEND.EQ.1) THEN
      IF(ILOS.EQ.10)                  CALL LOSS
      IF(ILOS.GE.100.AND.ILOS.LT.200) CALL NEW_LOSS
      IF(ILOS.GE.200)                 CALL SPAL_LOSS
      CALL OUTPUT
      END IF
C
C******************************************************************************
C******************************************************************************
C      CALL STEPUP EVERY TO UPDATE THE TIMESTEP IF ITIMST=3
C
      IF(ITIMST.GE.3) CALL STEPUP(0.25)
C
C******************************************************************************
C******************************************************************************
C
C      CALL EFICOOL TO PRINT OUT MASS AVERAGED QUANTITIES
C      AND MACHINE EFFICIENCY AND SHROUD LEAKAGE FLOWS EVERY 200 STEPS.
C
      IF(MOD(NSTEP,200).EQ.0.OR.IFPRINT.EQ.1.OR.IFEND.EQ.1) THEN
      CALL EFICOOL(HBLEDTOT)
C
C******************************************************************************
C******************************************************************************
C
C     WRITE OUT THE SHROUD LEAKAGE FLOWS AND FRICTIONAL WORK. ALSO EVERY 200 STEPS.
C
      IF(IFSHROUD.EQ.1) THEN
C
      DO 9099 NR = 1,NROWS
      IF(KTIPS(NR).LT.0) THEN
           WRITE(6,*)
           WRITE(6,*) 'ROW NO ',NR,'PERCENTAGE SHROUD LEAKAGE FLOW = ',
     &                 100.*SLEAK(NR)/BLADE_FLOW(1),' % '
	   WRITE(6,*) 'WORK DONE ON THE SHROUD AND HUB AND CASING BY THE
     & LEAKAGE FLOW = ',SWORK(NR),'WATTS'
      ENDIF
 9099 CONTINUE
C
       WRITE(6,*)
       WRITE(6,*) ' TOTAL FRICTIONAL WORK DONE ON ALL THE SHROUDS = ',
     &              SWORKTOT,' WATTS'
C
C    END OF OUTPUT FOR SHROUD LEAKAGE
      ENDIF
C
C   END OF OUTPUT EVERY 200 STEPS OR WHEN FINISHED OR CONVERGED.
      ENDIF
C
C******************************************************************************
C******************************************************************************
C      WRITE OUT A SHORT OUTPUT SUMMARY EVERY 5 STEPS
C
      IF(NSTEP.EQ.5.OR.MOD(NSTEP,50).EQ.0) WRITE(6,9002)ITIMST,CFL,
     &    DAMP,SFT,SFX,FEXTRAP,F1,F2EFF,F3,NRSMTH
 9002 FORMAT(/,' ISTEP=',I2,' CFL=',F5.2,' DAMP= ',F5.2,' SFT,SFX= ',
     &2F6.3,' FEXTRAP= ',F5.3,' F1=',F5.2,' F2EFF=',F5.2,' F3=',F5.2,
     &' NRSMTH=',I3,/)
C
      IF(ITIMST.GE.5.AND.MOD(NSTEP,50).EQ.0) WRITE(6,*)
     &   ' ARTIFICIAL SPEED OF SOUND   = ', VSOUND
C
C******************************************************************************
C******************************************************************************
C     THE FOLLOWING TIMING CALL WILL DIFFER FOR DIFFERENT COMPUTERS.
C
      IF(MOD(NSTEP,50).EQ.0) THEN
C      CALL CLOCK@(FINI)
      FINI = FLOAT(MCLOCK())
      RUNTIME   = (FINI - START)/50*1.0E-06
      POINTIME  = RUNTIME/(IM*JM*KM)
      WRITE(6,*) ' CPU TIME PER STEP=', RUNTIME, 'SECONDS.',
     & ' CPU TIME PER POINT PER STEP=', POINTIME,'SECONDS.'
      WRITE(6,*)
C      CALL CLOCK@(START)
      START = FLOAT(MCLOCK())
      ENDIF
C
C******************************************************************************
C     END OF TIMING CALL
C******************************************************************************
C******************************************************************************
C     WRITE OUT A SUMMARY TO THE SCREEN EVERY 5 STEPS
C
      IF(NSTEP.EQ.5.OR.MOD(NSTEP,50).EQ.0) WRITE(6,9001)
 9001 FORMAT(' STEP   EMAX AT I   J   K  EAVG    ECONT AT J=    VREF   V 
     &MAX  AT I   J   K     IN FLOW     OUT FLOW    RAT FLOW')
C
      RATFLOW = BLADE_FLOW(JM)/FLOWRAT
      WRITE(6,9000) NSTEP,EMAX,IMAX,JMAX,KMAX,EAVG,ECONT,JCONT,VREF,
     &              VMAX,IVMAX,JVMAX,KVMAX,FLOWRAT,BLADE_FLOW(JM),
     &              RATFLOW
 9000 FORMAT(I5,F8.4,3I4,2F8.4,I6,2F8.2,3I4,3F12.4)
C
C******************************************************************************
C******************************************************************************
C   WRITE A CONVERGENCE HISTORY TO UNIT 4 OUTPUT SUMMARY EVERY 5 STEPS
C
C
      WRITE(4,9009) EMAX,EAVG,ECONT,FLOWRAT,NSTEP,IMAX,JMAX,KMAX
 9009 FORMAT('EMAX,EAVG,ECONT',3E12.5,' FLOW',E12.5,'STEP',I5,'AT',3I5)
C
C******************************************************************************
C******************************************************************************
C    CALL MIX_BCONDS TO WRITE OUT THE PITCHWISE AVERAGE VALUES AT THE MIXING PLANES.
C
      IF(IFEND.EQ.1) THEN
      CALL MIX_BCONDS(1)
      END IF
C
C*******************************************************************************
C*********STOP IF CONVERGED OR MAXIMUM ITERATIONS REACHED**********
C
      IF(IFEND.EQ.1) STOP
C
C********************** RETURN TO START OF MAIN LOOP ************
C**************  IF NOT YET CONVERGED OR AT LAST ITERATION ******
C
      GO TO 5000
C
C******************************************************************************
C******************************************************************************
C
C     END OF THE SUBROUTINE LOOP
C
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE TSTEP(D,DIFF,NCALL)
C       ====================
C
C       THIS ROUTINE UPDATES THE GIVEN VARIABLE EVERY
C       TIMESTEP (I.E. RO ROVR ROVT ROVX ROE)
C       AND ALSO PERFORMS THE PITCHWISE SMOOTHING
C       MOST OF THE COMPUTATIONAL TIME IS USED BY THIS SUBROUTINE
C
C       ====================
      INCLUDE  'commall-open-18.3'
C
      DIMENSION D(ID,JD,KD), DIFF(ID,JD,KD),AVG(KD),
     &          B1CHG(IG1,JG1,KG1),B2CHG(IG2,JG2,KG2),SBCHG(JD)
C
C******************************************************************************
C     SET THE MULTIGRID CHANGES TO ZERO
C
      DO 110 K=1,NKB1
      DO 110 J=1,NJB1+1
      DO 110 I=1,NIB1
      B1CHG(I,J,K) = 0.0
  110 CONTINUE
      DO 210 K=1,NKB2
      DO 210 J=1,NJB2+1
      DO 210 I=1,NIB2
      B2CHG(I,J,K) = 0.0
  210 CONTINUE
      DO 310 J = 1,NSBLK
      SBCHG(J) = 0.0
  310 CONTINUE
C*******************************************************************************
C     BALANCE FLUXES ACROSS THE TIP GAP
C
      DO 520 J=2,JM
      NR = NROW(J)
      IF(KTIPS(NR).LE.0) GO TO 520
      DO 521      K = KTIPS(NR),KTIPE(NR)-1
      TFLUX(1,J,K)  = 0.5*(TFLUX(1,J,K)+TFLUX(IM,J,K))
      TFLUX(IM,J,K) = TFLUX(1,J,K)
  521 CONTINUE
  520 CONTINUE
C
C***********************************************************
C     EXTRAPOLATE THE FLUXES ON THE UPSTREAM FACE OF THE MIXING PLANES. 
C     UNLESS FEXTRAP = 0.0
C 
      IF(FEXTRAP.LT.0.001) GO TO 4141

      DO 4140 NR = 1,NRWSM1
      J = JMIX(NR)
      JP1 = J+1
      JM1 = J-1
      JP2 = J+2
C
      DO 4100 K = 1,KMM1
C
      FLOWDIRN = -FLOWX(IMID,J,K)
C
      AVFLXJM1 = 0.0
      AVFLXJP2 = 0.0
      DO 4120 I = 1,IMM1
      AVFLXJM1   = AVFLXJM1  + XFLUX(I,JM1,K)
      AVFLXJP2   = AVFLXJP2  + XFLUX(I,JP2,K)
 4120 CONTINUE
     
      DO 4130 I=1,IMM1
C
      IF(FLOWDIRN.GT.0.0) THEN
           DFLUXJM1 = XFLUX(I,JM1,K) - AVFLXJM1*FP(I)
           DFLUX  = DFLUXJM1*FEXTRAP
           SOURCE(I,J,K) = SOURCE(I,J,K) - DFLUX
      ELSE
           DFLUXJP2 = XFLUX(I,JP2,K) - AVFLXJP2*FP(I)
           DFLUX = DFLUXJP2*FEXTRAP
           SOURCE(I,JP2,K) = SOURCE(I,JP2,K) + DFLUX
      END IF
C
      SOURCE(I,JP1,K) = 0.0
C
 4130 CONTINUE
C
 4100 CONTINUE
C
 4140 CONTINUE

 4141 CONTINUE
C*******************************************************************************
C*******************************************************************************
C      SUM THE FLUXES TO FORM CHANGE AND SAVE IT IN STORE(I,J,K) 
C
      DO 1000 K=1,KMM1
      DO 1000 J=2,JM
      RATPITCH = FLOAT(NBLADE(J-1))/NBLADE(J)
      DO 1000 I=1,IMM1
      DELTA        = XFLUX(I,J,K) - XFLUX(I,J-1,K)*RATPITCH
     &             + RFLUX(I,J,K) - RFLUX(I,J,K+1)
     &             + TFLUX(I,J,K) - TFLUX(I+1,J,K)
     &             - SOURCE(I,J,K)
      STORE(I,J,K)  = F1*DELTA + F2*DIFF(I,J,K)
      DIFF(I,J,K)   = DELTA    + F3*DIFF(I,J,K)
 1000 CONTINUE
C
C**********************************************************************************
C**********************************************************************************
C     PITCHWISE AVERAGE THE CHANGES AT ANY MIXING PLANES
C     THIS IS ONLY DONE ON THE DOWNSTREAM FACE OF THE MIXING PLANE UNLESS FEXTRAP = 0
C     IN WHICH CASE IT IS DONE ON BOTH SIDES AS IN TBLOCK-13.
C
      DO 1750 NR = 1,NRWSM1
      J   = JMIX(NR)
      DO 1790 K = 1,KMM1
C
      FLOWDIRN = -FLOWX(IMID,J,K)
C    FIRSTLY AT J = JMIX + 2
      IF((FLOWDIRN.GT.0.0).OR.(FEXTRAP.LT.0.001)) THEN
      JMX = J + 2
      SUM_STORE = 0.0
      DO 1770 I=1,IMM1
      SUM_STORE  = SUM_STORE  + STORE(I,JMX,K)
 1770 CONTINUE
      DO 1780 I=1,IMM1
         STORE(I,JMX,K) = SUM_STORE*FP(I)
 1780 CONTINUE
      END IF
C  NEXT IF FEXTRAP = 0  AT J = JMIX
      IF((FLOWDIRN.LT.0.0).OR.(FEXTRAP.LT.0.001)) THEN
      SUM_STORE = 0.0
      DO 1771 I=1,IMM1
      SUM_STORE  = SUM_STORE  + STORE(I,J,K)
 1771 CONTINUE
      DO 1781 I=1,IMM1
         STORE(I,J,K) = SUM_STORE*FP(I)
 1781 CONTINUE
      END IF
C
 1790 CONTINUE
C
 1750 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C    JUMP TO 1020  IF NO MULTIGRID. THIS IS VERY UNUSUAL.
C
      IF(IR.LE.1.AND.JR.LE.1.AND.KR.LE.1) GO TO 1020
C
C*******************************************************************************
C*******************************************************************************
C      SUM THE ELEMENT CHANGES TO FORM THE CHANGES FOR THE MULTIGRID
C      BLOCKS. STORE(I,J,K) IS NOW USED AS A STORE FOR THE CHANGE IN
C      THE ELEMENTS.
C
      DO 700 K=1,KMM1
      K1 = KB1(K)
      K2 = KB2(K)
      DO 700 J=2,JM
      JSB= JSBLK(J)
      J1 = JB1(J)
      J2 = JB2(J)
      DO 700 I=1,IMM1
      I1 = IB1(I)
      I2 = IB2(I)
      DELTA = STORE(I,J,K)
      B1CHG(I1,J1,K1) = B1CHG(I1,J1,K1) + DELTA
      B2CHG(I2,J2,K2) = B2CHG(I2,J2,K2) + DELTA
      SBCHG(JSB)      = SBCHG(JSB)      + DELTA
  700 CONTINUE
C
 1020 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C     ADD THE BLOCK CHANGES TO THE ELEMENT CHANGES, MULTIPLY BY THE TIME STEP
C     JDD REMOVED THE CALCULATION OF THE AVERAGE CHANGE FROM THIS APRIL 2018.
C     IT IS NOW CALCULATED IN THE DO 1501 LOOP.
C
      DO 1500 K=1,KMM1
      K1 = KB1(K)
      K2 = KB2(K)
      DO 1500 J=2,JM
      JSB= JSBLK(J)
      J1 = JB1(J)
      J2 = JB2(J)
      DO 1500 I=1,IMM1
      I1 = IB1(I)
      I2 = IB2(I)
      DELTA = STORE(I,J,K)*STEP(I,J,K)
     &      + (B1CHG(I1,J1,K1)*STEP1(I1,J1,K1)
     &      + B2CHG(I2,J2,K2)*STEP2(I2,J2,K2)
     &      + SBCHG(JSB)*STEPSBK(JSB))*RSTEP(I,J,K)
      STORE(I,J,K) = DELTA
 1500 CONTINUE
C
c      END OF CALCULATING THE MULTIGRID CHANGES
C*******************************************************************************
C*******************************************************************************
C    USE RESIDUAL SMOOTHING IF NRSMTH > 0
C
      IF(NRSMTH.GT.0) THEN
      CALL SMOOTH_RESID(STORE,RSMTH,NRSMTH)
      END IF
C
C   
C*******************************************************************************
C*******************************************************************************
C
C     APPLY THE NEGATIVE FEEDBACK. SKIP IT IF DAMP IS SMALL OR LARGE
C
      IF(DAMP.LT.2.0.OR.DAMP.GT.100) GO TO 1550
C
C     CHANGED BY JDD TO USE THE AVERAGE CHANGE PER ROW - AVG_BLK(NR)- 14/02/2018
C     CALCULATE THE AVERAGE CHANGES FOR EACH ROW,  AVG_CHG(NR).
C     THIS IS NEW BY JDD  14/02/2018.
C
      DO 1501 NR = 1,NROWS
      SUMCHG = 0.0
      JST = JSTART(NR) + 1 
      JEN = JMIX(NR)   - 1
      JCHANGE = JEN - JST  + 1
      NSUM = IMM1*KMM1*JCHANGE
      DO 1502 K=1,KMM1
      DO 1502 J = JST, JEN
      DO 1502 I=1,IMM1
      SUMCHG = SUMCHG + ABS(STORE(I,J,K))
 1502 CONTINUE
      AVG_CHG(NR) = SUMCHG/NSUM
 1501 CONTINUE
C
C     SMOOTH THE BLADE ROW CHANGES
C
      IF(NROWS.EQ.1) AVG_BLK(1) = AVG_CHG(1)
      IF(NROWS.EQ.2) THEN
            AVG_BLK(2) = 0.5*(AVG_CHG(1) + AVG_CHG(2))
            AVG_BLK(1) = AVG_BLK(2)
      END IF
      IF(NROWS.GT.2) THEN
      DO 1503 NR = 2,NROWS-1
      AVG_BLK(NR) = 0.25*(AVG_CHG(NR-1)+AVG_CHG(NR+1)) + 0.5*AVG_CHG(NR)
 1503 CONTINUE
      AVG_BLK(NROWS) = 0.5*(AVG_CHG(NROWS-1) + AVG_CHG(NROWS))
      AVG_BLK(1)     = 0.5*(AVG_CHG(1) + AVG_CHG(2)) 
      END IF 
C
C*******************************************************************************
C      APPLY THE NEGATIVE FEEDBACK TO LIMIT THE MAXIMUM CHANGE.
C*******************************************************************************
C
      DO 1525 K=1,KMM1
      DO 1525 J=2,JM
      NR     = NROW(J)
      DO 1525 I=1,IMM1
      DELTA  = STORE(I,J,K)
      ABSCHG = ABS(DELTA)
      FDAMP  = ABSCHG/AVG_BLK(NR)
      STORE(I,J,K) = DELTA/(1. + FDAMP/DAMP )
 1525 CONTINUE
C
C   END OF JDD  14/02/2018  CHANGES.
C*******************************************************************************
C*******************************************************************************
C
 1550 CONTINUE
C   END OF JDD 14/02/2018 CHANGES
C*****************************************************************************
C******************************************************************************
C   NEXT IS SPECIAL TREATMENT FOR THE TURBULENT VISCOSITY IF USING THE SA MODEL
C
      IF(NCALL.EQ.2) THEN
C  
C     UPDATE THE TURBULENT VISCOSITY - REMEMBERNG THAT IT IS CELL CENTRED.
C    
      DO 1600 K=1,KMM1
      DO 1600 J=2,JM
      DO 1600 I=1,IMM1
      D(I,J,K) = D(I,J,K) + STORE(I,J,K)
 1600 CONTINUE
C
C   TRANSFER THE AVERAGE TURBULENT VISCOSITY ACROSS THE MIXING PLANE
C
      DO 1560 NR = 1,NROWS-1
      J = JMIX(NR) 
      FMULT = TURBVIS_DAMP(NR)
      DO 1570 K=1,KMM1
      AVG(K) = 0.0
      DO 1580 I=1,IMM1
      AVG(K) = AVG(K) + FP(I)*D(I,J,K)
 1580 CONTINUE
      JS = J + 1
      JE = J + 2
      DO 1591 JAV = JS,JE
      DO 1590 I=1,IMM1
 1590 D(I,JAV,K) = FMULT*AVG(K)
 1591 CONTINUE
 1570 CONTINUE
 1560 CONTINUE      
C
C******************************************************************************
C   SMOOTH THE TURBULENT VISCOSITY
C
      NSMTH = 1
      SFTVIS  = FAC_SFVIS*SFT
      CALL SMOOTH_RESID(D,SFTVIS,NSMTH)
C
C     JUMP TO 8700 FOR THE TURBULENT VISCOSITY ONLY
      GO TO 8700
C
C   END OF SPECIAL TREATMENT FOR THE TURBULENT VISCOSITY
C
      END IF
C
C**************************************************************************************
C**************************************************************************************
C     AGAIN PITCHWISE AVERAGE THE CHANGES AT ANY MIXING PLANES AFTER INCLUDING 
C     MULTIGRID AND MULTIPLYING BY THE TIMESTEP. 
C     THIS IS ONLY DONE ON THE DOWNSTREAM FACE OF THE MIXING PLANE UNLESS  FEXTRAP = 0
C     IN WHICH CASE IT IS DONE ON BOTH SIDES AS IN TBLOCK-13.
C
      DO 1850 NR = 1,NRWSM1
      J   = JMIX(NR)
      JP1 = J+1
      DO 1890 K = 1,KMM1
C
      FLOWDIRN = -FLOWX(IMID,J,K)
C
      IF((FLOWDIRN.GT.0.0).OR.(FEXTRAP.LT.0.001)) THEN
      JMX = J + 2
      SUM_STORE = 0.0
      DO 1870 I=1,IMM1
      SUM_STORE  = SUM_STORE  + STORE(I,JMX,K)
 1870 CONTINUE
      SUM_STORE = SUM_STORE/IMM1
      DO 1880 I=1,IMM1
         STORE(I,JMX,K) = SUM_STORE
 1880 CONTINUE
      END IF
C    IF FEXTRAP = 0  ALSO AT THE CELLS UPSTREAM OF THE MIXING PLANE, J = JMIX  
      IF((FLOWDIRN.LT.0.0).OR.(FEXTRAP.LT.0.001)) THEN
      SUM_STORE = 0.0
      DO 1871 I=1,IMM1
      SUM_STORE  = SUM_STORE  + STORE(I,J,K)
 1871 CONTINUE
      SUM_STORE = SUM_STORE/IMM1
      DO 1881 I=1,IMM1
         STORE(I,J,K) = SUM_STORE
 1881 CONTINUE
      END IF
C
 1890 CONTINUE
C
 1850 CONTINUE
C
C****************************************************************************
C****************************************************************************
C          ADD THE CHANGES TO THE OLD VALUE OF THE VARIABLE D
C              DISTRIBUTING THE CHANGES TO THE FOUR CORNERS
C              WITH DOUBLE WEIGHTING AT THE BOUNDARIES.
C              THE FACTOR OF 1/8 IS INCLUDED IN THE 'FMI' TERMS.
C****************************************************************************
C
      DO 1100 K=1,KMM1
      DO 1100 J=2,JM
      DO 1100 I=1,IMM1
C
      ADD = STORE(I,J,K)
C
      D(I,J,K)      =  D(I,J,K)       + ADD*FBL(I,K)*FACDWN(J)
      D(I+1,J,K)    =  D(I+1,J,K)     + ADD*FBR(I,K)*FACDWN(J)
      D(I,J,K+1)    =  D(I,J,K+1)     + ADD*FTL(I,K)*FACDWN(J)
      D(I+1,J,K+1)  =  D(I+1,J,K+1)   + ADD*FTR(I,K)*FACDWN(J)
      D(I,J-1,K)    =  D(I,J-1,K)     + ADD*FBL(I,K)*FACUP(J)
      D(I+1,J-1,K)  =  D(I+1,J-1,K)   + ADD*FBR(I,K)*FACUP(J)
      D(I,J-1,K+1)  =  D(I,J-1,K+1)   + ADD*FTL(I,K)*FACUP(J)
      D(I+1,J-1,K+1)=  D(I+1,J-1,K+1) + ADD*FTR(I,K)*FACUP(J)
 1100 CONTINUE
C
C************************************************************************************
C************************************************************************************
C    CALLL  "SMOOTHVAR" TO APPLY THE SMOOTHING (ARTIFICIAL VISCOSITY) TO THE VARIABLE "D". 
C
      CALL SMOOTH_VAR(D)
C
C*******************************************************************************
C*******************************************************************************
C    RE ENTER HERE IF NCALL = 2 FOR THE TURBULENT VISCOSITY.
 8700 CONTINUE
C*******************************************************************************
C*******************************************************************************
C**********APPLY THE PERIODIC BOUNDARY CONDITIONS AT I=1 AND IM ***************
C  
      ILAST = IM
      IF(NCALL.EQ.2) ILAST = IMM1
C
      DO 8000 J=1,JM
      IF((IND(J).EQ.1).OR.(INDLE(J).EQ.1)) GO TO 8000
      DO 7000 K=1,KM
      D(1,J,K)      = 0.5*(D(1,J,K)+D(ILAST,J,K))
      D(ILAST,J,K) = D(1,J,K)
 7000 CONTINUE
 8000 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C  Q3D
      IF(KM.EQ.2) GO TO 509
C  END Q3D
C*******************************************************************************
C*******************************************************************************
C
C     APPLY PERIODICITY ACROSS THE TIP GAP.
C     NOTE THAT THE TIP POINT ITSELF,  KTIP,  IS NOT PERIODIC.
C
      DO 510 J=1,JM
      NR = NROW(J)
      IF(KTIPS(NR).LE.0) GO TO 510
      K1 = KTIPS(NR)
      K2 = KTIPE(NR)
      IF(K1.EQ.1)  K2 = K2-1
      IF(K2.EQ.KM) K1 = K1+1
      DO 511 K=K1,K2
      D(1,J,K)     = 0.5*(D(1,J,K)+D(ILAST,J,K))
      D(ILAST,J,K) = D(1,J,K)
  511 CONTINUE
  510 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
  509 CONTINUE
C*******************************************************************************
C*******************************************************************************
C
C    PITCHWISE AVERAGE THE VARIABLES ON BOTH FACES OF THE MIXING PLANE
C
      IF(IFMIX.EQ.0) GO TO 1400
C
      DO 100 NR = 1,NRWSM1
      J   = JMIX(NR)
      JP1 = J+1
C
      DO 201 K=1,KM
      AVGVARJ    = 0.0
      AVGVARJP1  = 0.0
      DO 301 I=1,IMM1
      AVGVARJ   = AVGVARJ   + FP(I)*(D(I,J,K)+D(I+1,J,K)) 
      AVGVARJP1 = AVGVARJP1 + FP(I)*(D(I,JP1,K)+D(I+1,JP1,K))  
  301 CONTINUE
      AVGVARJ   = 0.5*AVGVARJ
      AVGVARJP1 = 0.5*AVGVARJP1
      DO 401 I=1,IM
      D(I,J,K)   = AVGVARJ
      D(I,JP1,K) = AVGVARJP1
  401 CONTINUE 
  201 CONTINUE
C 
  100 CONTINUE
C
C     END OF  MIXING PLANE TREATMENT
C
 1400 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C
      RETURN
      END
C*******************************************************************************
C*******************************************************************************
C
      SUBROUTINE OUTPUT
C       ====================
C
C       THIS ROUTINE CALCULATES FLOW PROPERTIES AND PRINTS HEADINGS
C       BEFORE CALLING PRINT TO OUTPUT THE RESULT ARRAYS
C
      INCLUDE  'commall-open-18.3'
C
C
      WRITE(6,*) 'IN OUTPUT, WRITING PLOTTING/RESTART FILE "flow_out" ' 
C
C        WRITE  THE COMBINED RESTART/PLOTTING  FILE   "flow_out"  IF  "IFPRINT"  = 1
C        OR  IF  "IFEND " = 1 WHEN THE PROGRAM HAS CONVERGED.
C
      PREFF = 0.5*(P(IMID,1,KMID) + P(IMID,JM,KMID))
      WRITE(7) NSTEP
      WRITE(7)(((RO(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(7)(((ROVX(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(7)(((ROVR(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      WRITE(7)(((ROVT(I,J,K),I=1,IM),J=1,JM),K=1,KM)
C     CORRECT THE VALUE OF ROE  TO THE TRUE VALUE BEFORE SENDING IT
C     TO OUTPUT FILE IF USING ARETIFICIAL COMPRESSIBILITY.
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) THEN
                WRITE(7)(((ROE(I,J,K)*RO(I,J,K)/ROSUB(I,J,K)
     &         ,I=1,IM),J=1,JM),K=1,KM)
      ELSE
                WRITE(7)(((ROE(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      END IF
C
      WRITE(7)(((ROSUB(I,J,K),I=1,IM),J=1,JM),K=1,KM)
C
C     WRITE OUT THE TRANSFORMED VISCOSITY IF USING THE SA MODEL
C     OTHERWISE WRITE OUT THE TURBULENT/LAMINAR VISCOSITY RATIO FOR PLOTTING ONLY.
      IF(ILOS.GE.200) THEN
          WRITE(7)(((TRANS_DYN_VIS(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      ELSE
          WRITE(7)(((DPDS_CELL(I,J,K)*1.0E-6,I=1,IM),J=1,JM),K=1,KM)
      END IF
C
C      WRITE OUT THE PEFF FOR THE SURFACE PRESSURES IN A THROUGHFLOW CALCULATION.
       WRITE(7)((( (PEFF(I,J,K)+PREFF),I=1,IM),J=1,JM),K=1,KM)
C      WRITE OUT SPAREVAR  WHICH CAN BE SET TO ANY REQUIRED 3D VARIBLE.
C      WRITE(7) (((SPAREVAR(I,J,K),I=1,IM),J=1,JM),K=1,KM)
C
      WRITE(6,*) 'IN OUTPUT, PLOTTING/RESTART FILE "flow_out" WRITTEN '
C
C********************************************************************************
C********************************************************************************
C********************************************************************************
C
C    NOW  WRITE OUT THE FORMATTED OUTPUT FILE  "results.out" TO UNIT 3 
C
C********************************************************************************
C********************************************************************************
C********************************************************************************
C
C     CALCULATE SOME FLOW PARAMETERS FOR PRINTING OUT 
C
C      TEMP1 = RELATIVE SWIRL VELOCITY
C      TEMP2 = RELATIVE MACH NUMBER
C      TEMP3 = PRESSURE IN BAR
C      TEMP4 = ABSOLUTE STAGNATION TEMPERATURE
C     
      DO 1100 I=1,IM
      DO 1100 K=1,KM
      DO 1100 J=1,JM
      EKE = 0.5*(VX(I,J,K)*VX(I,J,K) + VT(I,J,K)*VT(I,J,K)
     &         + VR(I,J,K)*VR(I,J,K))
C
      IF(IFGAS.EQ.0) THEN
            TSTATIC = (HO(I,J,K) - EKE)*RCP
            IF(TSTATIC.LT.1.) TSTATIC = 1.
            V_SONIC = SQRT(GA*RGAS*TSTATIC)
            TEMP4(I,J,K)  = HO(I,J,K)*RCP
      ELSE
            HSTAG   = HO(I,J,K)
            HSTAT   = HSTAG - EKE
            TSTAG   = TFROMH(HSTAG,TREF,HREF,HT1,HT2,HT3,HT4)
            TSTAT   = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
            TEMP4(I,J,K) = TSTAG
            IF(TSTAT.LT.1.) TSTAT = 1.
            CPNOW   = CP1 + CP2*(TSTAT-TREF) 
     &              + CP3*(TSTAT-TREF)*(TSTAT-TREF)
            GAMNOW  = CPNOW/(CPNOW-RGAS)
            V_SONIC = SQRT(GAMNOW*RGAS*TSTAT)
      END IF
C
      TEMP1(I,J,K)  = VT(I,J,K) - UBLADE(J,K)
      EKEREL  = EKE - 0.5*(VT(I,J,K)*VT(I,J,K)
     &              - TEMP1(I,J,K)*TEMP1(I,J,K))
      TEMP3(I,J,K)  = 0.00001*P(I,J,K)
      TEMP2(I,J,K)  = SQRT (2.0*EKEREL)/V_SONIC
 1100 CONTINUE
C
C***************************************************************************
C***************************************************************************
C    STORE IS THE PERCENTAGE CHANGE IN MERIDIONAL VELOCITY IN THE LAST TIME STEP.
C
      WRITE(3,5) TITLE
    5 FORMAT(1H ,24X,18A4/)
      IF(IOUT(1).EQ.0) GO TO 10
      IOPT=IOUT(1)
      WRITE(3,5400) NSTEP
 5400 FORMAT(1H ,'** PERCENTAGE CHANGE IN VM ** TIMESTEP =',I5,' **')
      CALL PRINT (IOPT,STORE)
C
   10 CONTINUE
C
      IF(IOUT(2).EQ.0) GO TO 20
      IOPT=IOUT(2)
      WRITE(3,5100) NSTEP
 5100 FORMAT(1H ,'*** AXIAL VELOCITY  M/SEC TIMESTEP =',I5,' **')
      CALL PRINT(IOPT,VX)
C
   20 CONTINUE
C
      IF(IOUT(3).EQ.0) GO TO 30
      IOPT=IOUT(3)
      WRITE(3,5200) NSTEP
 5200 FORMAT(1H ,' RADIAL VELOCITY M/SEC TIMESTEP =',I5,' **')
      CALL PRINT(IOPT,VR)
C
   30 CONTINUE
C
      IF(IOUT(4).EQ.0) GO TO 40
      IOPT=IOUT(4)
      WRITE(3,5300) NSTEP
 5300 FORMAT(1H0,' RELATIVE SWIRL  VELOCITY M/SEC TIMESTEP =',I5,' **')
      CALL PRINT(IOPT,TEMP1)
C
   40 CONTINUE
C
      IF(IOUT(5).EQ.0) GO TO 50
      IOPT=IOUT(5)
      WRITE(3,5500) NSTEP
 5500 FORMAT(1H0,'********** PRESSURE IN BAR *** TIMESTEP =',I5,' **')
      CALL PRINT(IOPT,TEMP3(I,J,K))
C
   50 CONTINUE
C
      IF(IOUT(6).EQ.0) GO TO 60
      IOPT=IOUT(6)
      WRITE(3,5600) NSTEP
 5600 FORMAT(1H ,' ** RELATIVE  MACH  NUMBER  ** TIMESTEP =',I5,' **')
      CALL PRINT(IOPT,TEMP2)
C
   60 CONTINUE
C
      IF(IOUT(7).EQ.0) GO TO 70
      WRITE(3,5700) NSTEP
      IOPT=IOUT(7)
 5700 FORMAT(1H0,'ABSOLUTE STAG. TEMPERATURE DEG K TIMESTEP= ',I5,'**')
      CALL PRINT(IOPT,TEMP4)
C
C     CALCULATE THE ENTROPY FUNCTION
C
      DO 1200 I=1,IM
      DO 1200 K=1,KM
      DO 1200 J=1,JM
      D=.5*(VX(I,J,K)*VX(I,J,K)+VT(I,J,K)*VT(I,J,K)+VR(I,J,K)*VR(I,J,K))
C
      IF(IFGAS.EQ.0) THEN
            TSTATIC = (HO(I,J,K) - D)*RCP
            IF(TSTATIC.LT.1.) TSTATIC = 1.
            STORE(I,J,K) = P(I,J,K)/PO1(KMID)*
     &      (TO1(KMID)/TSTATIC)**(GA/(GA-1))
      ELSE
           TSTAG  = TEMP4(I,J,K)
           HSTAT  = HSTAG - D
           TSTAT  = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
           TRAT   = TSTAT/TSTAG
           PRAT   = PRAT_FROM_TRAT(TO1(KMID),TRAT,ALPHA,BETA1,BETA2)
           STORE(I,J,K) = P(I,J,K)/(PO1(KMID)*PRAT)
      END IF
C
 1200 CONTINUE
C
   70 CONTINUE
C
      IF(IOUT(12).EQ.0) GO TO 75
      IOPT=IOUT(12)
      WRITE(3,5710) NSTEP
 5710 FORMAT(1H0,' P/(T**(GA/(GA-1)))/ INLET VALUE ON MID K GRID LINE,
     &  TIMESTEP NO ',I5)
      CALL PRINT(IOPT,STORE)
C
   75 CONTINUE
C*******************************************************************************
C    CALCULATE MORE PROPERTIES FOR PRINTING OUT
C
C          STORE IS MERIDIONAL VELOCITY
C          TEMP1 IS RELATIVE SWIRL ANGLE, ATAN(VTREL/VM)
C          TEMP2 IS MERIDIONAL PITCH ANGLE, ATAN(VR/VX)
C
      PIN = 0.0
      DO 2000 K=1,KM
      DO 2010 I=1,IM
      DO 2020 J=1,JM
      STORE(I,J,K) = SQRT(VX(I,J,K)*VX(I,J,K)+VR(I,J,K)*VR(I,J,K))
      TEMP1(I,J,K) = ATAN(TEMP1(I,J,K)/STORE(I,J,K)) * 57.296
      TEMP2(I,J,K) = ASIN(VR(I,J,K)/STORE(I,J,K)) * 57.296
 2020 CONTINUE
      PIN = PIN + P(I,1,K)
 2010 CONTINUE
 2000 CONTINUE
C**********************************************************************************
C
      IF(IOUT(8).EQ.0) GO TO 80
      IOPT=IOUT(8)
      WRITE(3,2030) NSTEP
 2030 FORMAT(1H0,' MERIDIONAL VELOCITY, M/SEC, TIMESTEP=',I5 )
      CALL PRINT(IOPT,STORE)
C
   80 CONTINUE
C
      IF(IOUT(9).EQ.0) GO TO 90
      IOPT=IOUT(9)
      WRITE(3,2040) NSTEP
 2040 FORMAT(1H0,' SWIRL ANGLE ON STREAM SURFACE, ATAN(WTREL/VM), DEG,
     1 TIMESTEP=',I5)
      CALL PRINT(IOPT,TEMP1)
C
   90 CONTINUE
C
      IF(IOUT(10).EQ.0) GO TO 100
      IOPT=IOUT(10)
      WRITE(3,2050) NSTEP
 2050 FORMAT(1H0,' MERIDIONAL PITCH ANGLE, ATAN(VR/VX), DEG, TIMESTEP=
     1 ',I5)
      CALL PRINT(IOPT,TEMP2)
C
  100 CONTINUE
C
      IF(IOUT(11).EQ.0) GO TO 110
      IOPT=IOUT(11)
      WRITE(3,2060) NSTEP
 2060 FORMAT(1H0,' DENSITY KG/M**3, TIMESTEP=   ',I5)
      CALL PRINT(IOPT,RO)
C
  110 CONTINUE
C
      IF(IOUT(12).EQ.0) GO TO 120
      IOPT=IOUT(12)
      WRITE(3,2070) NSTEP
 2070 FORMAT(1H0,' ARTIFICIAL DENSITY KG/M**3, TIMESTEP=   ',I5)
      CALL PRINT(IOPT,ROSUB)
C
  120 CONTINUE
C
      IF(IOUT(13).EQ.0) GO TO 130
      IOPT=IOUT(13)
      PIN=PIN/(IM*KM)
      DO 113 K=1,KM
      DO 113 I=1,IM
      DO 113 J=1,JM
      STORE(I,J,K)=(P(I,J,K)-PIN)/(PO1(KMID)-PIN)
  113 CONTINUE
      WRITE(3,2080) NSTEP
 2080 FORMAT(1H0, ' PRESSURE COEFFICIENT BASED ON  AVERAGE INLET
     & STATIC PRESSURE AND STAGNATION PRESSURE AT KMID',I5)
      CALL PRINT(IOPT,STORE)
C
  130 CONTINUE
C
C    END OF PRINTED OUTPUT
C
C************************************************************************************
C************************************************************************************
C
      RETURN
      END
C
C************************************************************************************
C************************************************************************************
C
      SUBROUTINE PRINT(IOPT,F)
C
C       ====================
C       THIS ROUTINE PRINT'S OUT 3-D ARRAYS OR PITCHWISE AVERAGE VALUES.
C       ====================
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION F(ID,JD,KD),SUMFUN(MAXKI)
C
      IF(IOPT.EQ.1) GO TO 4000
      IF(IOPT.EQ.3) GO TO 4000
      IF(IOPT.EQ.2) GO TO 500
      RETURN
C**********************************************************************
C      PRINT OUT WHOLE FLOW FIELD OF THE INPUT VARIABLE  "F(I,J,K)" .
C
 4000 CONTINUE
C
      DO 1000 K=1,KM
      IF(KOUT(K).EQ.0 )GO TO 1000
      WRITE(3,5200) K
 5200 FORMAT(/' STREAM SURFACE NUMBER  ',I4,/)
      DO 1100 J=1,JM
      WRITE(3,5300) J,R(J,K),X(J,K),SMERID(J,K),(F(I,J,K),I=1,IM)
 1100 CONTINUE
 1000 CONTINUE
 5300 FORMAT(1H ,'J= ',I4,'R= ',F10.5,' X= ',F10.5,' MERIDIONAL DISTANCE
     & =',F10.5,' VALUE I=1->IM =',/,(10F11.3))
C
      RETURN
C
C************************************************************************
C      PRINT OUT PITCHWISE MASS AVERAGED VALUES IF REQUESTED. THE VALUES ARE THE 
C      AVERAGED VALUES IN EACH CELL SO ONLY KM-1 VALUES
C
  500 CONTINUE
C
      WRITE(3,2000)
 2000 FORMAT( /,'   FRACTIONAL SPAN OF EACH K LINE ' )
      WRITE(3,2001) (FSPAN(K),K=1,KMM1)
 2001 FORMAT(10F12.4)
      WRITE(3,*)
      WRITE(3,*)
C
      DO 100 J=1,JM
      DO 110 K=1,KMM1
      SUMAS     = 0.0
      SUMFUN(K) = 0.0
      DO 120 I=1,IMM1
      SUMAS     = SUMAS + FLOWX(I,J,K)
      SUMFUN(K) = SUMFUN(K)+FLOWX(I,J,K)*(F(I,J,K)+F(I+1,J,K)
     &          + F(I+1,J,K+1)+F(I,J,K+1))*0.25
  120 CONTINUE
      SUMFUN(K) = SUMFUN(K)/SUMAS
  110 CONTINUE
C
      WRITE(3,200) J,R(J,1),R(J,KM),(SUMFUN(K),K=1,KMM1)
  200 FORMAT( /,' J=',I5,' RHUB=',F10.4,' RTIP=',F10.4,' PITCHWISE MASS 
     & AVERAGE, K=1,KMM1,= ',/, (10F12.3))
C
  100 CONTINUE
C
      RETURN
      END
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE SETUP(ANSIN)
C
C**********THIS SUBROUTINE SETS UP THE GRID AND INITIALISES ALL THE
C          VARIABLES USED IN THE MAIN PROGRAM.
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION JSTART1(JD),JEND1(JD),JSTART2(JD),JEND2(JD),PDOWN(KD),
     &          INDLETE(JD)
C
      CHARACTER*1 ANSW,ANSIN
C
C      SET VARIOUS CONSTANTS AND INTEGERS NEEDED THROUGHOUT THE CALCULATION
C
      WRITE(6,1111)
 1111 FORMAT('  ENTERING SUBROUTINE SETUP SO DATA INPUT WAS OK ')
C
      JLEP1  = JLE(1)+1
      NRWSM1 = NROWS-1
      IMID   = IFIX(0.5*IM)
      IF(IM.EQ.2) IMID = 1
      KMID   = IFIX(0.5*KM)
      IF(KM.EQ.2) KMID = 1
      IFPRINT= 0
      IFEND  = 0
C
      DO J=1,JM
           INDLETE(J) = 0
           IF(INDLE(J).EQ.1) INDLETE(J) = 1
           IF(INDTE(J).EQ.1) INDLETE(J) = 1
      END DO 
C  Q3D
      IF(KM.EQ.2) THEN
           KMID = 1
           KMM2 = 1
      END IF
C  END Q3D

      DO 1000 K=1,KM
      BSSIN(K) = SIN(BS(K)*DEGRAD)
      BSCOS(K) = COS(BS(K)*DEGRAD)
      BRSIN(K) = SIN(BR(K)*DEGRAD)
 1000 BRCOS(K) = COS(BR(K)*DEGRAD)
C
C******************************************************************************
C     MAKE THE FR(K) SUM TO UNITY
C
      SUM = 0.0
      DO 2000 K=1,KMM1
 2000 SUM = SUM+FR(K)
      DO 2001 K=1,KMM1
 2001 FR(K) = FR(K)/SUM
C
C    FSPAN IS THE FRACTIONAL SPAN AT THE CENTRE OF THE ELEMENTS.
      FSPAN(1) = 0.5*FR(1)
      DO 1999 K=2,KMM1
 1999 FSPAN(K) = FSPAN(K-1)+0.5*(FR(K)+FR(K-1))
C
C    MAKE THE FP(I) SUM TO UNITY
      SUM=0.0
      DO 2002 I=1,IMM1
 2002 SUM = SUM+FP(I)
      DO 2003 I=1,IMM1
 2003 FP(I) = FP(I)/SUM
C 
C  THE FU(I)  AND FD(I)  ARE USED IN THE SMOOTHING ROUTINE.
C
      IF(IM.GT.2) THEN
C
      DO 2004 I=2,IMM1
      FU(I) = FP(I)/(FP(I)+FP(I-1))
 2004 FD(I) = FP(I-1)/(FP(I)+FP(I-1))
      FU(1) = FP(1)/FP(2)
      FU(IM)= FP(IMM1)/FP(IMM2)  
C
      ELSE
C   FOR TFLOW
      FU(1) = 0.5
      FD(1) = 0.5
      FU(2) = 0.5
      FD(2) = 0.5
C
      END IF    
C
C      SET FKU(K)  AND  FKD(K) WHICH  ARE USED IN THE SMOOTHING ROUTINE.
C
      IF(KM.GT.2) THEN
           DO 2005 K=2,KMM1
           FKU(K) = FR(K)/(FR(K)+FR(K-1))
 2005      FKD(K) = FR(K-1)/(FR(K)+FR(K-1))
           FKU(1) = FR(1)/FR(2)
           FKU(KM)= FR(KMM1)/FR(KMM2)
      ELSE
C   Q3D
           FKU(1) = 0.5
           FKU(2) = 0.5
      END IF
C   END Q3D
C******************************************************************************
C
C     SET THE DISTRIBUTION FUNCTIONS FOR USE IN SUBROUTINE TSTEP
C
      DO 2006 K=1,KMM1
      DO 2006 I=1,IMM1
           FTL(I,K)  = 0.125
           FBL(I,K)  = 0.125
           FTR(I,K)  = 0.125
           FBR(I,K)  = 0.125
 2006 CONTINUE
      DO 2007 K=1,KMM1
           FTL(1,K)    = 0.25
           FBL(1,K)    = 0.25
           FTR(IMM1,K) = 0.25
           FBR(IMM1,K) = 0.25
 2007 CONTINUE
      DO 2008 I=1,IMM1
           FBL(I,1)    = 0.25
           FBR(I,1)    = 0.25
           FTL(I,KMM1) = 0.25
           FTR(I,KMM1) = 0.25
 2008 CONTINUE
           FTL(1,KMM1)    = 0.5
           FBL(1,1)       = 0.5
           FTR(IMM1,KMM1) = 0.5
           FBR(IMM1,1)    = 0.5
      DO J = 2,JM
           FACUP(J)  = 1.0
           FACDWN(J) = 1.0
      END DO
           FACUP(2)    = 2.0
           FACDWN(JM)  = 2.0
C
C********************************************************************************
C********************************************************************************
C
C      CALL INTPOL TO INTERPOLATE IN INPUT DATA TO SET UP GRID NODES ON THE
C      BLADE SURFACES.
C      IF KM = 2 THEN THE TWO INPUT STREAM SURFACES ARE THE HUB AND CASING
C      AND THERE IS NO INTERPOLATION.
C  Q3D
      IF(KM.EQ.2) THEN
           DO K=1,2
           DO J=1,JM
                X(J,K) = XSURF(J,K)
                R(J,K) = RSURF(J,K)
           END DO
           END DO
      ELSE
C
C  END Q3D
C
      WRITE(6,*) ' INPUT STYLE, ANSIN = ', ANSIN
C
      IF(ANSIN.EQ.'n'.OR.ANSIN.EQ.'N')THEN
          WRITE(6,*) ' CALLING NEW_INTPOL FROM SETUP'
          CALL NEW_INTPOL
      ELSE
          WRITE(6,*) ' CALLING OLD_INTPOL FROM SETUP'          
          CALL OLD_INTPOL
      END IF
C
      WRITE(6,1888)
 1888 FORMAT( '  INTERPOLATION IN INPUT BLADE SECTIONS COMPLETED OK ')
C
      ENDIF
C
C******************************************************************************
C******************************************************************************
C     SHIFT THE BLADES SO THEY ARE ALIGNED CIRCUMFERENTIALLY AT I=1 AT MID-SPAN
C
      SHIFT = 0.0
      DO 3022 J=3,JM
      IF(INDMIX(J-1).EQ.1) SHIFT = (RT_UPP(J-1,KMID)
     &                           -  RT_UPP(J,KMID))/R(J,KMID)
      DO 3022 K=1,KM
      RT_UPP(J,K) = RT_UPP(J,K) + SHIFT*R(J,K)
 3022 CONTINUE
C
C******************************************************************************
C     SCALE THE BLADE THICKNESS FOR THE PINCHED TIP MODEL OF TIP CLEARANCE.
C
      DO 2506 J=1,JM
      NR = NROW(J)
      DO 2505 K=1,KM
      RTMID = RT_UPP(J,K) - 0.5*RT_THICK(J,K)
      THICK = FTHICK(NR,K)*RT_THICK(J,K)
      RT_UPP(J,K)    = RTMID + 0.5*THICK
      RT_THICK(J,K)  = THICK
 2505 CONTINUE
 2506 CONTINUE
C
C*******************************************************************************
C      SET VARIOUS CONSTANTS
C
      GA1   = GA - 1.0
      FGA   = GA1/GA
      RFGA  = 1.0/FGA
      RCP   = 1.0/CP
      CV    = CP/GA
      RCV   = 1.0/CV
      RGAS  = CP - CV
      SFEX1 = 1. - SFEXIT
      SFEXH = 0.5*SFEXIT
C
C*******************************************************************************
C          SET TEMPORARY GRID COORDINATES FOR USE IN WORKING OUT AREAS.
C
C       THETA .......STORE FOR THETA CO-ORD  OF GRID NODES
C       RTHETA ......STORE FOR RTHETA CO-ORD OF GRID NODES
C       DX ..........AXIAL EXTENT OF ELEMENTS
C       DS ... ......MERIDIONAL EXTENT OF ELEMENTS.
C       RT_PITCH ... LOCAL BLADE PITCH.
C
      DO 3000 J=1,JM
      PITCH = 2.*PI/NBLADE(J)
      DO 3021 K=1,KM
      GAP           = PITCH*R(J,K) - RT_THICK(J,K)
      RTHETA(1,J,K) = RT_UPP(J,K)
      THETA(1,J,K)  = RT_UPP(J,K)/R(J,K)
      RT_PITCH(J,K) = PITCH*R(J,K)
      DO 3020 I=2,IM
      RTHETA(I,J,K) = RTHETA(I-1,J,K) + FP(I-1)*GAP
      THETA(I,J,K)  = RTHETA(I,J,K)/R(J,K)
 3020 CONTINUE
 
 3021 CONTINUE
 3000 CONTINUE
C
C*******************************************************************************
C      EVALUATE THE MERIDIONAL DISTANCE, SMERID.
C
      DSMIN=1.0E06
      DO 3530 K=1,KM
      SMERID(1,K) = 0.0
      DO 3540 J=2,JM
      DR(J,K)     = (R(J,K)-R(J-1,K))
      DX(J,K)     = (X(J,K)-X(J-1,K))
      DS(J,K)     = SQRT(DR(J,K)*DR(J,K) + DX(J,K)*DX(J,K))
      SMERID(J,K) = SMERID(J-1,K) + DS(J,K)
      IF(DS(J,K).LT.DSMIN) DSMIN = DS(J,K)
 3540 CONTINUE
      DX(1,K) = DX(2,K)
      DS(1,K) = DS(2,K)
      DR(1,K) = DR(2,K)
 3530 CONTINUE
C
C   JDD ADDITION AUGUST 2017
C   FIND THE BLADE CHORDS FOR EVERY ROW
      DO 3550 NR = 1,NROWS
      JLEDGE    = JLE(NR)
      JTEDGE    = JTE(NR)
      SCHORD    = SMERID(JTEDGE,KMID)-SMERID(JLEDGE,KMID)
      TCHORD    = RTHETA(1,JLEDGE,KMID) - RTHETA(1,JTEDGE,KMID)
      CHORD(NR) = SQRT(SCHORD*SCHORD + TCHORD*TCHORD)
 3550 CONTINUE
C
C      WRITE "SMERID" TO THE FILE "GLOBPLOT"
C
      WRITE(11)    JM,1,1
      WRITE(11)    (SMERID(J,KMID), J=1,JM)
C
C*******************************************************************************
C   SET RAVG_CELL(J,K) = AVERAGE RADIUS OF AN ELEMENT.
C
      DO 4000 J=2,JM
      DO 4000 K=1,KMM1 
      RAVG_CELL(J,K) = 0.25*(R(J,K)+R(J-1,K)+R(J-1,K+1)+R(J,K+1))
 4000 CONTINUE
C
C     SET THE BLADE SPEED  UBLADE(J,K)
C
      DO 4005 J=1,JM
      DO 4005 K=1,KM
      UBLADE(J,K) = WRAD(J)*R(J,K)
 4005 CONTINUE
C
C*******************************************************************************
C   SET FMUP ,FMDN  FOR USE IN THE SMOOTHING ROUTINE.
C
      DO 4010 K=1,KM
      DO 4010 J=2,JMM1
      DSUP      = SMERID(J,K)   - SMERID(J-1,K)
      DSDWN     = SMERID(J+1,K) - SMERID(J,K)
      FMDN(J,K) = DSDWN/(DSUP + DSDWN)
      FMUP(J,K) =  DSUP/(DSUP + DSDWN)
 4010 CONTINUE      
C
C*******************************************************************************
C      STORE EXTRAPOLATED VALUES OF THETA AND R*THETA ON THE MIXING PLANE
C      TEMPORARILY AS TEMP1 AND TEMP2
C
      DO 3025 J=1,JM
      IF(INDMIX(J).NE.1) GO TO  3025
      DO 3024 K=1,KM
      FAC = DS(J+1,K)/DS(J+2,K)
      DO 3024 I=1,IM
      TEMP1(I,J,K) = THETA(I,J+1,K)+ FAC*(THETA(I,J+1,K)-THETA(I,J+2,K))
      TEMP2(I,J,K) = R(J,K)*TEMP1(I,J,K)
 3024 CONTINUE
 3025 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     WRITE THE GRID GEOMETRY TO UNIT 21 FOR USE WHEN PLOTTING
C
      OPEN(UNIT=21,FILE='grid_out',form= 'unformatted' )
C
      WRITE(21) NSTEPS_MAX
      WRITE(21) IM,JM,KM
      WRITE(21) CP,GA
      WRITE(21) (INDLETE(J),J=1,JM)
      WRITE(21) (WRAD(J),J=1,JM)
      WRITE(21) (NBLADE(J),J=1,JM)
C 
      DO 20 J=1,JM
      DO 20 K=1,KM 
      WRITE(21) X(J,K),R(J,K),(RTHETA(I,J,K),I=1,IM)
   20 CONTINUE
C
      CLOSE(21)
C
C******************************************************************************
C******************************************************************************
C           WORK OUT THE PROJECTED AREAS OF THE BLADE FACES
C
      WRITE(6,*) ' STARTING TO WORK OUT AREAS OF FACES OF ELEMENTS'
C
C           QUASI ORTHOGONAL FACE FIRST.
C
      DO 3030 J=1,JM
      DO 3030 I=1,IMM1
      DO 3015 K=1,KMM1
      X1=X(J,K+1)-X(J,K)
      X2=X1
      R1=R(J,K+1)-R(J,K)
      R2=R1
      T1=R(J,K+1)*(THETA(I+1,J,K+1)-THETA(I,J,K))
      T2=R(J,K+1)*(THETA(I,J,K+1)-THETA(I,J,K))
     &    -R(J,K)*(THETA(I+1,J,K)-THETA(I,J,K))
      AQX(I,J,K) = -0.5*(R2*T1-R1*T2)
      AQR(I,J,K) = -0.5*(X1*T2-X2*T1)
      AQTOT(I,J,K) = SQRT(AQX(I,J,K)*AQX(I,J,K)+AQR(I,J,K)*AQR(I,J,K))
 3015 CONTINUE
      AQX(I,J,KM)=AQX(I,J,KMM1)
 3030 AQR(I,J,KM)=AQR(I,J,KMM1)
C
      DO 3031 J=1,JM
      DO 3031 K=1,KM
      AQX(IM,J,K) = AQX(IMM1,J,K)
      AQR(IM,J,K) = AQR(IMM1,J,K)
 3031 CONTINUE
C
C***********NEXT WORK OUT THE AREAS OF THE BLADEWISE FACE*********
C
      DO 3040 J=2,JM
      DO 3040 K=1,KMM1
      X1 = X(J-1,K+1)-X(J,K)
      X2 = X(J,K+1)-X(J-1,K)
      R1 = R(J-1,K+1)-R(J,K)
      R2 = R(J,K+1)-R(J-1,K)
      ABT(J,K) = -0.5*(X1*R2-X2*R1)
C
      WRABT(J,K) = WRAD(J)*0.25*(R(J,K)+R(J-1,K)+R(J-1,K+1)+R(J,K+1))
     &            *ABT(J,K)
C
      DO 3040 I=1,IM
      T1 = R(J-1,K+1)*(THETA(I,J-1,K+1)-THETA(I,J,K))
      T2 =   R(J,K+1)*(THETA(I,J,K+1)-THETA(I,J,K))
     &   -   R(J-1,K)*(THETA(I,J-1,K)-THETA(I,J,K))
C
      IF(INDMIX(J-1).EQ.1) THEN
      T1 = R(J-1,K+1)*(TEMP1(I,J-1,K+1)-THETA(I,J,K))
      T2 =   R(J,K+1)*(THETA(I,J,K+1)-THETA(I,J,K))
     &   -   R(J-1,K)*(TEMP1(I,J-1,K)-THETA(I,J,K))
      ENDIF
C
      ABX(I,J,K) = -0.5*(R1*T2-R2*T1)
      ABR(I,J,K) = -0.5*(X2*T1-X1*T2)
      ABTOT(I,J,K) = SQRT(ABX(I,J,K)**2 + ABR(I,J,K)**2 + ABT(J,K)**2)
C
      IF(I.EQ.IM) GO TO 3040
C
C   WORK OUT THE VOLUME OF THE ELEMENTS
C
      VOL(I,J,K) = ABT(J,K)*0.25*(RTHETA(I+1,J,K) + RTHETA(I+1,J-1,K)
     &  + RTHETA(I+1,J-1,K+1) + RTHETA(I+1,J,K+1) - RTHETA(I,J,K)
     &  - RTHETA(I,J-1,K)     - RTHETA(I,J-1,K+1) - RTHETA(I,J,K+1))
C
      IF(INDMIX(J-1).EQ.1)
     &    VOL(I,J,K) = ABT(J,K)*0.25*(RTHETA(I+1,J,K) + TEMP2(I+1,J-1,K)
     &  + TEMP2(I+1,J-1,K+1) + RTHETA(I+1,J,K+1) - RTHETA(I,J,K)
     &  - TEMP2(I,J-1,K)     - TEMP2(I,J-1,K+1)  - RTHETA(I,J,K+1))
C
 3040 CONTINUE
C
C**********WORK OUT AREAS OF STREAMWISE FACE.
C
      DO 3050 J=2,JM
      DO 3050 I=1,IMM1
      DO 3050 K=1,KM
      X1=X(J-1,K)-X(J,K)
      X2=X1
      R1=R(J-1,K)-R(J,K)
      R2=R1
      T1 = R(J-1,K)*(THETA(I,J-1,K)  - THETA(I,J,K))
     &   -   R(J,K)*(THETA(I+1,J,K)  - THETA(I,J,K))
      T2 = R(J-1,K)*(THETA(I+1,J-1,K)- THETA(I,J,K))
C
      IF(INDMIX(J-1).EQ.1) THEN
      T1 = R(J-1,K)*(TEMP1(I,J-1,K)  - THETA(I,J,K))
     &   -   R(J,K)*(THETA(I+1,J,K)  - THETA(I,J,K))
      T2 = R(J-1,K)*(TEMP1(I+1,J-1,K)- THETA(I,J,K))
      ENDIF
C
      ASX(I,J,K)=0.5*(R1*T2-R2*T1)
      ASR(I,J,K)=0.5*(X2*T1-X1*T2)
C
 3050 CONTINUE
C
C**********************************************************************************
C     SET THE VOLUME/RADIUS TERM NEEDED FOR THE RADIAL MOMENTUM EQUATION.
C
      NNEG = 0
      DO 4045 K=1,KMM1
      DO 4045 J=2,JM
      DO 4045 I=1,IMM1
      IF(VOL(I,J,K).LT.0.0) NNEG = NNEG + 1
      RAREA  = ABR(I,J,K) - ABR(I+1,J,K)
     &       + AQR(I,J,K) - AQR(I,J-1,K)
     &       + ASR(I,J,K) - ASR(I,J,K+1)
      VOLOR(I,J,K) = - 0.125*RAREA
C      RATIO = ABS(RAREA)/ABS(ASR(I,J,K))
C      CHECK CLOSURE OF THE VOLUMES IN THE RADIAL DIRECTION
C      IF( RATIO.GT.0.00001 )THEN
C      WRITE(6,*) ' RAREA NOT WELL  CLOSED AT  J , K =', J,K
C      WRITE(6,*) ' ASR, RAREA',  ASR(I,J,K), RAREA, RATIO
C      END IF

      XAREA  = AQX(I,J,K) - AQX(I,J-1,K)
     &       + ABX(I,J,K) - ABX(I+1,J,K)
     &       + ASX(I,J,K) - ASX(I,J,K+1)
C      RATIO = ABS(XAREA)/ABS(AQX(I,J,K) )
C      CHECK CLOSURE OF THE VOLUMES IN THE AXIAL DIRECTION
C      IF( RATIO.GT.0.00001 )THEN
C      WRITE(6,*) ' XAREA NOT WELL  CLOSED AT  J , K =', J,K
C      WRITE(6,*) ' AQX, XAREA',  AQX(I,J,K), XAREA, RATIO
C      END IF

 4045 CONTINUE
C
      WRITE(6,*)  ' AREAS OF ALL FACES OF THE ELEMENTS EVALUATED '
C
C**********************************************************************************
C**********************************************************************************
C
C************SET THE INITIAL GUESS OF P AND RO AT ALL GRID POINTS********
C
      WRITE(6,*)' STARTING TO MAKE THE INITIAL GUESS OF THE FLOW FIELD.'
C
C   SET THE SPANWISE VARIATION OF THE  EXIT PRESSURE IF IPOUT = 3.
      IF(IPOUT.EQ.3) GO TO 3212
      PD(1) = PDOWN_HUB
      DO 3211 K=2,KM
      PD(K) = PD(K-1) + FR(K-1)*(PDOWN_TIP - PDOWN_HUB)
 3211 CONTINUE
 3212 CONTINUE
C
C**********************************************************************************
C     SET VALUES NEEDED IF USING ARTIFICIAL COMPRESSIBILITY
C
      IF(ITIMST.GE.5) THEN
           PO_REF  = PO1(KMID)
           RO_REF  = PO1(KMID)/RGAS/TO1(KMID)
           IF(ITIMST.EQ.6) RO_REF = DENSTY
           DP_DRO  = VSOUND*VSOUND
      END IF
C
C*************************************************************************  
C  SET THE INITIAL GUESS OF PRESSURE, DENSITY AND TEMPERATURE
C
      DO 3210 K=1,KM
           RPO1(K)   = 1./PO1(K)
           VR(1,1,K) = VM1(K)*(R(2,K)-R(1,K))/DS(1,K)
           VINSQ     = VM1(K)*VM1(K)+VTIN(K)*VTIN(K)
C
      IF(IFGAS.EQ.0) THEN
            T1  = TO1(K) - 0.5*VINSQ/CP
            P1  = PO1(K)*(T1/TO1(K))**RFGA
      ELSE
            HSTAG    = HFROMT(TO1(K),TREF,HREF,CP1,CP2,CP3)
            HSTAT    = HSTAG - 0.5*VINSQ
            T1       = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
            TRAT     = T1/TO1(K)
            PRAT     = PRAT_FROM_TRAT(TO1(K),TRAT,ALPHA,BETA1,BETA2)
            P1       = PO1(K)*PRAT
      END IF
C          
      DO 3210 J=1,JM
           PS   = PGUESS(J)
           PRAT = PS/PO1(K)
C
      IF(IFGAS.EQ.0) THEN
           TS = (PRAT**FGA)*TO1(K)
      ELSE
           TRAT = TRAT_FROM_PRAT(TO1(K),PRAT,FGAGAS,R_ALPHA,
     &            BALPHA1,BALPHA2)
           TS   = TO1(K)*TRAT
      END IF
C
           ROS = PS/TS/RGAS
           IF(ITIMST.EQ.6) ROS  = DENSTY
C
      DO 3220 I=1,IM
           P(I,J,K)     = PS
           RO(I,J,K)    = ROS
           IF(ITIMST.EQ.5.OR.ITIMST.EQ.6)
     &     ROSUB(I,J,K) = RO_REF - (PO_REF - P(I,J,K))/DP_DRO
           ROAVG_CELL(I,J,K) = ROS
           T_STATIC(I,J,K)   = TS
 3220 CONTINUE
 3210 CONTINUE
C
C
C**********************************************************************************
C
      WRITE(6,*) ' INITAL GUESS OF P & RO COMPLETED '
C
C**********************************************************************************
C   SET THE INITIAL GUESS OF VELOCITY COMPONENTS.  MAKE THE TANGENTIAL VELOCITY
C   VARY SMOOTHLY BETWEEN BLADED AND UNBLADED REGIONS USING  "FACJ" .
C
      KAVG = 0.75*KM
      DO 3300 K=1,KM
      DO 3310 I=1,IM
      VT(I,1,K) = VTIN(K)
      WT(I,1,K) = VTIN(K) - UBLADE(1,K)
C
      IF(IFGAS.EQ.0) THEN
           HO(I,1,K) = CP*TO1(K)
      ELSE
           HO(I,1,K) = HFROMT(TO1(K),TREF,HREF,CP1,CP2,CP3)
      END IF
C
      FACJ = 0.9
      DO 3320 J = 2,JM
C
      IF(J.GT.2.AND.INDTE(J-2).EQ.1) FACJ = 1.0
C
      AREA=SQRT(AQX(IMID,J,KAVG)*AQX(IMID,J,KAVG)+AQR(IMID,J,KAVG)
     &         *AQR(IMID,J,KAVG))*NBLADE(J)
      IF(J.EQ.2) AIN=AREA
      VS        = VM1(KAVG)*AIN*RO(IMID,1,KAVG)/(AREA*RO(IMID,J,KAVG))
      VTGRID = VS*R(J,K)*(THETA(I,J,K) - THETA(I,J-1,K))/DS(J,K)
     &          + UBLADE(J,K)
C  RELAX  vtheta   BETWEEN THE GRID ANGLE VAUE AND THE UPSTREAM POINT VALUE
      VT(I,J,K) =  (1.- FACJ)*VTGRID + FACJ*VT(I,J-1,K)*R(J-1,K)/R(J,K)
C     OMIT THE VALUE AT THE MIXING PLANE WHERE THE GRID ANGLE VALUE IS WRONG
      IF(INDMIX(J-1).EQ.1) VT(I,J,K) = VT(I,J-1,K)
C     USE THE INLET BOUNDARY VALUE UPSTREAM OF THE FIRST LEADING EDGE.
      IF(J.LT.JLE(1))  VT(I,J,K) = VTIN(K)*R(1,K)/R(J,K)
C
      WT(I,J,K) = VT(I,J,K) - UBLADE(J,K)
      VX(I,J,K) = VS*DX(J,K)/DS(J,K)
      VR(I,J,K) = VS*DR(J,K)/DS(J,K)
      HO(I,J,K) = HO(I,J-1,K)+WRAD(J)*(R(J,K)*VT(I,J,K)-R(J-1,K)
     &            *VT(I,J-1,K))
C   RELAX THE UPSTREAM VALUE BY FACJ = 0.9.
      IF(INDLE(J).EQ.1) FACJ = 0.9
C
      IF(K.EQ.KMID.AND.I.EQ.IMID) THEN 
          WRITE(6,3321) J, T_STATIC(I,J,K), P(I,J,K), RO(I,J,K),
     &                  VX(I,J,K), VR(I,J,K), VT(I,J,K) 
      END IF
 3321 FORMAT('J=  ',I5,'  INITIAL GUESS OF: T, P ,RO, VX, VR, VT',
     &      F10.2, F12.1, 4F10.3)
C
 3320 CONTINUE
      VR(I,1,K) = VR(I,2,K)
      VX(I,1,K) = VX(I,2,K)
      WT(I,1,K) = WT(I,2,K)
 3310 CONTINUE
 3300 CONTINUE
C
C TFLOW 
      IF(IM.EQ.2) THEN
           DO 3333 K=1,KM
           DO 3333 J=1,JM
                VTAVG = 0.5*(VT(1,J,K) + VT(2,J,K))
                VT(1,J,K) = VTAVG
                VT(2,J,K) = VTAVG
 3333      CONTINUE
      END IF
C  END TFLOW 
C*********************************************************************************
C     SET THE MASS FLUXES, ROVX, ROVR, RORVT  AND ROE '
C
      PREFF = 0.5*(P(IMID,1,KMID) + P(IMID,JM,KMID))
C
      DO 3500 J=1,JM
      DO 3500 K=1,KM
      DO 3500 I=1,IM
      EKE  = 0.5*(VX(I,J,K)**2+VT(I,J,K)**2+VR(I,J,K)**2)
C
      IF(IFGAS.EQ.0) THEN
           ROE(I,J,K)    = P(I,J,K)/(GA-1) + EKE*RO(I,J,K)
      ELSE
           ESTAT         = HO(I,J,K) - EKE - P(I,J,K)/RO(I,J,K) 
           ROE(I,J,K)    = RO(I,J,K)*(ESTAT + EKE)
      END IF
C
      IF(ITIMST.GE.5) ROE(I,J,K)=ROSUB(I,J,K)*((HO(I,J,K)-EKE)/GA + EKE)
C
      ROS           = RO(I,J,K)
      ROVR(I,J,K)   = ROS*VR(I,J,K)
      ROVT(I,J,K)   = ROS*VT(I,J,K)
      ROWT(I,J,K)   = ROS*WT(I,J,K)
      RORVT(I,J,K)  = ROVT(I,J,K)*R(J,K)
      ROVX(I,J,K)   = ROS*VX(I,J,K)
      PEFF(I,J,K)   = P(I,J,K) - PREFF
      DRO(I,J,K)    = 0.0
      DROE(I,J,K)   = 0.0
      DROVX(I,J,K)  = 0.0
      DROVR(I,J,K)  = 0.0
      DRORVT(I,J,K) = 0.0
      DPDS_CELL(I,J,K)  = 0.0
      DEL_DYNVIS(I,J,K) = 0.0
      Y_PLUS(I,J,K) = 25.0
      VISC_RAT(I,J,K) = 1.0
      TRANS_DYN_VIS(I,J,K) = 0.0
      YPLUS_K1(I,J) = 25.0
      YPLUS_KM(I,J) = 25.0
      YPLUS_I1(J,K) = 25.0
      YPLUS_IM(J,K) = 25.0
 3500 CONTINUE
C
      WRITE(6,*)  ' DONE INITIAL GUESS OF ALL FLOW VELOCITIES ETC '
C
C**********************************************************************************
C**********************************************************************************
C      READ IN FROM RESTART FILE IF "IF_RESTART" = 1 TO OVERWRITE INITIAL GUESS
C
      IF(IF_RESTART.EQ.0) GO TO 3700
C
C******************************************************************************
      WRITE(6,*)
      WRITE(6,*)   ' READING IN THE RESTART FILE FROM UNIT 7.'
C
      READ(7) NSTEP
      READ(7)(((RO(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      READ(7)(((ROVX(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      READ(7)(((ROVR(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      READ(7)(((ROVT(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      READ(7)(((ROE(I,J,K),I=1,IM),J=1,JM),K=1,KM)
      READ(7)(((ROSUB(I,J,K),I=1,IM),J=1,JM),K=1,KM) 
      READ(7)(((TRANS_DYN_VIS(I,J,K),I=1,IM),J=1,JM),K=1,KM)
C
      REWIND(7)
C
      WRITE(6,*) ' RESTART FILE READ IN OK '
      WRITE(6,*)
C**********************************************************************************
C**********************************************************************************
C   SETTING THE SECONDARY VARIABLES FROM THE RESTART FILE OF PRIMARY VARIABLES.
C
      DO 3901 K=1,KM
      DO 3900 J=1,JM
      DO 3900 I=1,IM
      VX(I,J,K)    = ROVX(I,J,K)/RO(I,J,K)
      VR(I,J,K)    = ROVR(I,J,K)/RO(I,J,K)
      VT(I,J,K)    = ROVT(I,J,K)/RO(I,J,K)
      RORVT(I,J,K) = ROVT(I,J,K)*R(J,K)
      EKE = 0.5*(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K) 
     &   +       VT(I,J,K)*VT(I,J,K))
C   THE VALUE OF ROE SENT TO THE PLOT FILE WAS THE TRUE VALUE
C    CHANGE BACK TO THE ARTIFICIAL VALUE "ROSUB*E" .
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6)
     &         ROE(I,J,K) = ROE(I,J,K)*ROSUB(I,J,K)/RO(I,J,K)
C
      IF(IFGAS.EQ.0) THEN
           TSTATIC   = (ROE(I,J,K)/RO(I,J,K) - EKE)/CV
           T_STATIC(I,J,K) = TSTATIC
           HO(I,J,K) = CP*TSTATIC + EKE
      ELSE
           ESTAT     = ROE(I,J,K)/RO(I,J,K)  - EKE
           TSTATIC   = TFROME(ESTAT,TREF,EREF,ET1,ET2,ET3,ET4)
           T_STATIC(I,J,K) = TSTATIC
           HSTAT     = HFROMT(TSTATIC,TREF,HREF,CP1,CP2,CP3)
           HO(I,J,K) = HSTAT + EKE
      END IF
C
      IF(ITIMST.LT.5) THEN
      ROSUB(I,J,K) = RO(I,J,K)
      P(I,J,K)     = RO(I,J,K)*RGAS*TSTATIC
      END IF
C
      IF(ITIMST.GE.5) P(I,J,K) = PO_REF - DP_DRO*(RO_REF - ROSUB(I,J,K))
C
      PEFF(I,J,K)  = P(I,J,K)

 3900 CONTINUE
C
      PDOWN(K) = 0.0
      DO 3902 I=1,IMM1
      PDOWN(K) = PDOWN(K) + FP(I)*0.5*(P(I,JM,K)+P(I+1,JM,K))
 3902 CONTINUE
C
 3901 CONTINUE
C
      IF(IPOUT.GE.1)  GO TO 3700
      IF(IPOUT.EQ.0)  PDIFF =   PD(1)  - PDOWN(1)
      IF(IPOUT.EQ.-1) PDIFF =   PD(KM) - PDOWN(KM) 
C         
      DO 3903 K=1,KM
      PD(K) = PDOWN(K) + PDIFF
 3903 CONTINUE
C
      WRITE(6,*) ' ALL FLOW PROPERTIES SET UP FROM THE RESTART FILE.'
C
C******************************************************************************
C   END OF SETTING UP THE FLOW FROM THE RESTART FILE.
 3700 CONTINUE
C
C******************************************************************************
C******************************************************************************
C        SET THE VALUE OF STEP(I,J,K)= TIMESTEP/VOLUME FOR EACH ELEMENT.
C        STEP(J,K) IS THE MAIN TIME STEP MULTIPLIED BY 1/VOLUME
C        IF ITIMST = 2  FIXED NON UNIFORM TIME STEPS ARE TAKEN WITH DT
C        PROPORTIONAL TO DS. IF ITIMST =3 THE THE TIME STEP IS
C        REGULARLY UPDATED BY SUBROUTINE STEPUP.
C
      WRITE(6,*) ' SETTING THE TIME STEP. '
C
      IF(ITIMST.LT.5) THEN
           IF(IFGAS.EQ.0) THEN
                VSOUND = SQRT(RGAS*GA*TO1(KMID))
           ELSE
                CPNOW = CP1 + CP2*(TO1(KMID)-TREF) 
     &                + CP3*(TO1(KMID)-TREF)*(TO1(KMID)-TREF)
                GAMNOW = CPNOW/(CPNOW-RGAS)
                VSOUND = SQRT(GAMNOW*RGAS*TO1(KMID))
            END IF
      END IF
C
C     SET THE TIME STEP BASED ON THE LOCAL GRID DIMENSION AND A CONSTANT SPEED OF SOUND.
C     THIS IS MODIFIED TO THE TRUE SPEED OF SOUND IN SUBROUTINE "STEPUP" . 
C
      DO 6000 J=2,JM
      DO 6000 K=1,KMM1
      DO 6000 I=1,IMM1
      SMIN = DSMIN
      ATTOT = SQRT(ABT(J,K)*ABT(J,K)+ABX(I,J,K)*ABX(I,J,K)
     &           + ABR(I,J,K)*ABR(I,J,K))
      PERPT = VOL(I,J,K)/ATTOT
      PERPQ = VOL(I,J,K)/AQTOT(I,J,K)
      ASTOT = SQRT(ASX(I,J,K)*ASX(I,J,K)+ ASR(I,J,K)*ASR(I,J,K))
      PERPS = VOL(I,J,K)/ASTOT
      SMIN  = AMIN1(PERPT,PERPQ,PERPS)
      STEP(I,J,K)  = SMIN*CFL/VOL(I,J,K)/VSOUND
      BSTEP(I,J,K) = STEP(I,J,K)
      RSTEP(I,J,K) = 1.0
 6000 CONTINUE
C
C     SET THE TIME STEP TO AN AVERAGE AT THE MIXING PLANE
C
      DO 101 NR=1,NRWSM1
      J = JMIX(NR)
      DO 201 K=1,KMM1
      DO 201 I=1,IMM1
      STEP(I,J+1,K)  = 0.5*(STEP(I,J,K)+STEP(I,J+2,K))
      BSTEP(I,J+1,K) = STEP(I,J+1,K)
      RSTEP(I,J+1,K) = 0.0
  201 CONTINUE
  101 CONTINUE
C
C     CALL STEPUP TO SET THE LOCALLY VARYING TIME STEPS IF ITIMST >= 3 .
C
      IF(ITIMST.GE.3) THEN
           CALL STEPUP(1.0)
      ENDIF
C
C*******************************************************************************
C*******************************************************************************
C
      WRITE(6,*)  ' SETTING UP THE MULTIGRID ARRAYS '
C
C         STARTING TO SET UP ARRAYS FOR USE IN MULTIGRID
C
      IF(JR.EQ.1.AND.IR.EQ.1.AND.KR.EQ.1) GO TO 7100
C
      DO 7000 K=1,KMM1
      KB2(K) =  1 + (K-1)/KRBB
 7000 KB1(K) =  1 + (K-1)/KR
      NKB1 = KB1(KMM1)
      NKB2 = KB2(KMM1)
C
      DO 7002 I=1,IMM1
      IB2(I) = 1 + (I-1)/IRBB
 7002 IB1(I) = 1 + (I-1)/IR
      NIB1 = IB1(IMM1)
      NIB2 = IB2(IMM1)
C
C     Sort the j values into multigrid blocks making sure that
C     No blocks overlap a mixing plane
C
      JSTT = 2
      JB1(1) = 1
      JB1(2) = 1
      DO 7003 J=3,JM
      MARK = 0
      IF(INDMIX(J).EQ.1.OR.J.EQ.JM)       MARK = 1
      IF((MOD((J-JSTT),JR).EQ.0).AND.(MARK.EQ.0)) THEN
           JB1(J) = JB1(J-1) + 1
      ELSE
           JB1(J) = JB1(J-1)
      ENDIF
      IF(INDMIX(J-1).EQ.1) JB1(J) = 0
      IF(INDMIX(J-2).EQ.1) THEN
           JB1(J) = JB1(J-2) + 1
           JSTT   = J
      ENDIF
 7003 CONTINUE
C
      JSTT = 2
      JB2(1) = 1
      JB2(2) = 1
      DO 7004 J=3,JM
      JP1 = J+1
      IF(J.EQ.JM) JP1 = JM
      MARK = 0
      IF(INDMIX(J).EQ.1.OR.J.EQ.JM)       MARK =1
      IF(INDMIX(JP1).EQ.1.OR.(J+1).EQ.JM) MARK =1
      IF((MOD((J-JSTT),JRBB).EQ.0).AND.(MARK.EQ.0)) THEN
           JB2(J) = JB2(J-1) + 1
      ELSE
           JB2(J) = JB2(J-1)
      ENDIF
      IF(INDMIX(J-1).EQ.1) JB2(J) = 0
      IF(INDMIX(J-2).EQ.1) THEN
           JB2(J)=JB2(J-2) + 1
           JSTT = J
      ENDIF
 7004 CONTINUE
C
      NJB1 = JB1(JM)
      NJB2 = JB2(JM)
C
C     CHECK THAT THE MULTIGRID DIMENSIONS ARE NOT TOO LARGE
C
      IF(NIB1.GT.IG1) WRITE(6,*) 'STOPPING BECAUSE NIB1 TOO LARGE.'
      IF(NKB1.GT.KG1) WRITE(6,*) 'STOPPING BECAUSE NKB1 TOO LARGE.'
      IF(NJB1.GT.JG1) WRITE(6,*) 'STOPPING BECAUSE NJB1 TOO LARGE '
      IF(NIB1.GT.IG1.OR.NKB1.GT.KG1.OR.NJB1.GT.JG1) STOP
C
C
      IF(NIB2.GT.IG2) WRITE(6,*) 'STOPPING BECAUSE NIB2 TOO LARGE.'
      IF(NKB2.GT.KG2) WRITE(6,*) 'STOPPING BECAUSE NKB2 TOO LARGE.'
      IF(NJB2.GT.JG2) WRITE(6,*) 'STOPPING BECAUSE NJB2 TOO LARGE '
      IF(NIB2.GT.IG2.OR.NKB2.GT.KG2.OR.NJB2.GT.JG2) STOP

      JSTART1(1)  = 2
      JSTART2(1)  = 2
      DO 7005 J=2,JM
      IF(JB1(J).EQ.0) JB1(J) = NJB1+1
      IF(JB2(J).EQ.0) JB2(J) = NJB2+1
      IF(JB1(J).NE.JB1(J-1)) JSTART1(JB1(J)) = J
      JEND1(JB1(J)) = J
      IF(JB2(J).NE.JB2(J-1)) JSTART2(JB2(J)) = J
      JEND2(JB2(J)) = J
 7005 CONTINUE
C
      JSTART1(NJB1+1) = JM
      JEND1(NJB1+1)   = JM
      JSTART2(NJB2+1) = JM
      JEND2(NJB2+1)   = JM
C
C
 7100 CONTINUE
C
C******************************************************************************
C        SET UP BLOCK SIZES FOR THE SUPERGRID JSBLK(J) IS SUPERBLOCK
C                           INDICATOR
C
      WRITE(6,*)  ' SETTING UP THE SUPER BLOCK ARRAYS',
     &   ' THERE ARE ONLY 4 SUPERBLOCKS PER BLADE ROW.'
      NSB=1
      DO 199 J=1,JM
      JSBLK(J) = NSB
      IF(INDLE(J).EQ.1) NSB=NSB+1
      IF(INDMID(J).EQ.1)NSB=NSB+1
      IF(INDTE(J).EQ.1) NSB=NSB+1
  199 CONTINUE
C      WRITE(6,7102)(JSBLK(J),J=1,JM)
C 7102 FORMAT( ' SUPER BLOCK INDEX= ',20I5)
C
      NSBLK = JSBLK(JM)
C
      IF(NSBLK.GT.JG3) THEN
      WRITE(6,*)    ' STOPPING BECAUSE THE TOTAL NUMBER OF SUPER BLOCKS,
     & NSBLK, IS TOO LARGE.'
      STOP
      ENDIF
C
C******************************************************************************
C     SET THE TIME STEPS FOR THE MULTIGRID BLOCKS
C     FIRST FOR THE LEVEL 1 BLOCKS
C
      DO 7110  I1 = 1,NIB1
      ISTART = (I1-1)*IR + 1
      IEND = ISTART + IR -1
      IF(IEND.GT.IMM1) IEND = IMM1
      DO 7110  J1 = 1,NJB1
      JSTRT = JSTART1(J1)
      JENDD = JEND1(J1)
      DO 7110  K1 = 1,NKB1
      KSTART = (K1-1)*KR + 1
      KEND = KSTART + KR -1
      IF(KEND.GT.KMM1) KEND = KMM1
      PERPT = 0.0
      PERPQ = 0.0
      PERPS = 0.0
      VOLB  = 0.0
      DO 7120 I = ISTART,IEND
      DO 7120 K = KSTART,KEND
      DO 7120 J = JSTRT,JENDD
      ATTOT = SQRT(ABT(J,K)*ABT(J,K)+ABX(I,J,K)*ABX(I,J,K)
     & + ABR(I,J,K)*ABR(I,J,K))
      PERPT = PERPT +  VOL(I,J,K)/ATTOT
      PERPQ = PERPQ + VOL(I,J,K)/AQTOT(I,J,K)
      ASTOT = SQRT(ASX(I,J,K)*ASX(I,J,K)+ ASR(I,J,K)*ASR(I,J,K))
      PERPS = PERPS + VOL(I,J,K)/ASTOT
      VOLB  = VOLB + VOL(I,J,K)
 7120 CONTINUE
      PERPS = PERPS/(IR*JR)
      PERPT = PERPT/(JR*KR)
      PERPQ = PERPQ/(IR*KR)
      PERPMIN = AMIN1(PERPS,PERPT,PERPQ)
      STEP1(I1,J1,K1)     = CFL*FBLK1*PERPMIN/VSOUND/VOLB
      STEP1(I1,NJB1+1,K1) = 0.0
 7110 CONTINUE
C********************************************************************************
C     NOW FOR THE LEVEL  2  BLOCKS
C
      DO 8110  I2 = 1,NIB2
      ISTART = (I2-1)*IRBB + 1
      IEND   = ISTART + IRBB - 1
      IF(IEND.GT.IMM1) IEND = IMM1
      DO 8110  J2 = 1,NJB2
      JSTRT = JSTART2(J2)
      JENDD = JEND2(J2)
      DO 8110  K2 = 1,NKB2
      KSTART = (K2-1)*KRBB + 1
      KEND = KSTART + KRBB - 1
      IF(KEND.GT.KMM1) KEND = KMM1
      PERPT = 0.0
      PERPQ = 0.0
      PERPS = 0.0
      VOLB  = 0.0
      DO 8120 I = ISTART,IEND
      DO 8120 K = KSTART,KEND
      DO 8120 J = JSTRT,JENDD
      ATTOT = SQRT(ABT(J,K)*ABT(J,K)+ABX(I,J,K)*ABX(I,J,K)
     & + ABR(I,J,K)*ABR(I,J,K))
      PERPT = PERPT +  VOL(I,J,K)/ATTOT
      PERPQ = PERPQ + VOL(I,J,K)/AQTOT(I,J,K)
      ASTOT = SQRT(ASX(I,J,K)*ASX(I,J,K)+ ASR(I,J,K)*ASR(I,J,K))
      PERPS = PERPS + VOL(I,J,K)/ASTOT
      VOLB  = VOLB + VOL(I,J,K)
 8120 CONTINUE
      PERPS = PERPS/(IRBB*JRBB)
      PERPT = PERPT/(JRBB*KRBB)
      PERPQ = PERPQ/(IRBB*KRBB)
      PERPMIN = AMIN1(PERPS,PERPT,PERPQ)
      STEP2(I2,J2,K2)     = CFL*FBLK2*PERPMIN/VSOUND/VOLB
      STEP2(I2,NJB2+1,K2) = 0.0
 8110 CONTINUE
C*******************************************************************************
C     NOW FOR THE SUPER BLOCKS
C
      PERPQ = 0.0
      VOLB  = 0.0
      DO 9010 J = 2,JM
      JSB = JSBLK(J)
      DO 9020 I=1,IMM1
      DO 9020 K=1,KMM1
      PERPQ = PERPQ + VOL(I,J,K)/AQTOT(I,J,K)
      VOLB  = VOLB  + VOL(I,J,K)
 9020 CONTINUE
C
      IF(J.EQ.JM) THEN
      PERPAVG      = PERPQ/(IMM1*KMM1)
      STEPSBK(JSB) = CFL*FBLK3*PERPAVG/VSOUND/VOLB
      GO TO 9010
      ENDIF
      IF(JSBLK(J+1).NE.JSB) THEN
      PERPAVG = PERPQ/(IMM1*KMM1)
      STEPSBK(JSB) = CFL*FBLK3*PERPAVG/VSOUND/VOLB
      PERPQ = 0.0
      VOLB  = 0.0
      ENDIF
 9010 CONTINUE
C
C     END OF SETTING THE MULTIGRID TIME STEPS
C
C******************************************************************************
C*****************************************************************************
C     INITIALISE THE SOURCE TERMS TO ZERO
C     NOTE THAT THIS SECTION HAS BEEN MOVED FOR THE SA MODEL BECAUSE IT NEEDS STEP(I,J,K)
C     TO BE SET.
C
      DO 3710 K=1,KM
      DO 3710 J=1,JM
C  TFLOW
      PBLADE(J,K)   = 0.0
      ROVAR_M1(J,K) = 0.0
      BFORCE_X(J,K) = 0.0
      BFORCE_R(J,K) = 0.0
      BFORCE_T(J,K) = 0.0
      BFORCE_Q(J,K) = 0.0
C END TFLOW
      DO 3710 I=1,IM
      XFORCE(I,J,K)   = 0.0
      TFORCE(I,J,K)   = 0.0
      RFORCE(I,J,K)   = 0.0
      QSOURCE(I,J,K)  = 0.0
      SOURCE(I,J,K)   = 0.0
      T_SOURCE(I,J,K) = 0.0
      SGEN(I,J,K)     = 0.0
 3710 CONTINUE
C     
C
      IF(NNEG.GT.0) THEN
           WRITE(6,*)'*******************WARNING***********************'
           WRITE(6,*)   NNEG,' NEGATIVE VOLUMES FOUND ' 
           WRITE(6,*) ' THIS IS VERY LIKELY TO CAUSE FAILURE. '
           WRITE(6,*) ' DO YOU WANT TO CONTINUE DESPITE THIS ?'
           WRITE(6,*) ' ANSWER   Y  or  N '
           READ(1,*)    ANSW
           IF(ANSW.EQ.'N'.OR.ANSW.EQ.'n')  STOP
      END IF
      CLOSE(1)
C*******************************************************************************
C  JDD ADDED 30/9/10. TO SET A LIMIT ON THE VORTICITY
C
      DR1 = R(1,KM)  - R(1,1)
      DX1 = X(1,KM)  - X(1,1)
      DRM = R(JM,KM) - R(JM,1)
      DXM = X(JM,KM) - X(JM,1)
      AVGSPAN  = 0.5*(SQRT(DX1*DX1+DR1*DR1) + SQRT(DXM*DXM+DRM*DRM))
      VORT_MAX = 5000.*VM1(KMID)/AVGSPAN
C
      WRITE(6,*) ' VORTICITY LIMIT, VORT_MAX = ', VORT_MAX
C
C  END OF JDD ADDITION
C******************************************************************************
C
C CALL SET_XLENGTH TO CALCULATE THE WALL DISTANCES AND SET THE MIXING LENGTHS
C
      WRITE(6,*)    ' CALLING  SET_XLENGTH TO SET THE WALL DISTANCES AND 
     & MIXING LENGTHS'
C
       CALL SET_XLENGTH
C
      WRITE(6,*) ' CALLED  SET_XLENGTH '
C
C******************************************************************************
C******************************************************************************
C      CALL LOSS ROUTINES TO INITIALISE THE BODY FORCE TERMS
C
      WRITE(6,*)  ' CALLING THE LOSS ROUTINES FROM SETUP, ILOS = ',ILOS
C
      TEMP_RF_VIS  = RF_VIS
      RF_VIS   = 1.0
      FMIXUP   = 1.0
      NSTEP    = 1 
      IF(ILOS.EQ.10)                  CALL LOSS
      IF(ILOS.GE.100.AND.ILOS.LT.200) CALL NEW_LOSS
      IF(ILOS.GE.200)                 CALL SPAL_LOSS
      RF_VIS   = TEMP_RF_VIS
C
C******************************************************************************
C******************************************************************************
C  INITIALISE THE SHROUD FLOWS, BLEED FLOWS AND COOLING FLOWS TO ZERO.
      DO 9030 J = 1,JM
      SHRDFLOW(J) = 0.0
      SUMBLEED(J) = 0.0
      SUMCWL(J)   = 0.0
 9030 CONTINUE
C
C******************************************************************************
C TFLOW
C******************************************************************************
C
      IF(IM.EQ.2) THEN
C
C  FIND THE CENTRE LINE GRID ANGLE ALPHA_CENT AND SMOOTH ITS STREAMWISE VARIATION.
      DO 9039 K=1,KMM1
      S_DIST(1) = 0.0
      DO 9038  J=2,JM
           ABXAVG     = 0.5*( ABX(1,J,K)  + ABX(2,J,K) )
           ABRAVG     = 0.5*( ABR(1,J,K)  + ABR(2,J,K) )
           VECX       = AQX(1,J,K)/AQTOT(1,J,K)
           VECR       = AQR(1,J,K)/AQTOT(1,J,K)
           ABMERAVG   = ABXAVG*VECX  + ABRAVG*VECR
           CENT_ANGL(J) = ATAN(ABMERAVG/ABT(J,K))
           XDIF  = X(J,K) - X(J-1,K)
           RDIF  = R(J,K) - R(J-1,K)
           S_DIST(J) = S_DIST(J-1) + SQRT(XDIF*XDIF + RDIF*RDIF)
 9038 CONTINUE
C
C     SET THE CENTRE LINE ANGLES UPSTREAM AND DOWNSTREAM OF A BLADE ROW.
      DO 9036 J=1,JM
      NRW   = NROW(J)
      JLEE  = JLE(NRW)
      JTEE  = JTE(NRW) 
C      AVGLE = (CENT_ANGL(JLEE)+CENT_ANGL(JLEE+1)+CENT_ANGL(JLEE+2))/3
C      AVGLE  = 2.0*CENT_ANGL(JLEE+1) - CENT_ANGL(JLEE+2)
C      AVGTE = (CENT_ANGL(JTEE)+CENT_ANGL(JTEE-1)+CENT_ANGL(JTEE-2))/3
C      AVGTE  = 2.0*CENT_ANGL(JTEE-1) - CENT_ANGL(JTEE-2)
C
C   JDD CHANGED THIS TO SETTING THE ANGLES TO EQUAL THE LE AND TE BLADE ANGLES.
C   MAY 2017.
      IF(J.LE.JLEE) CENT_ANGL(J) = CENT_ANGL(J+1)
      IF(J.GE.JTEE) CENT_ANGL(J) = CENT_ANGL(J-1)
 9036 CONTINUE
C
C    SMOOTH THE CENTRE LINE ANGLE VARIATION.
      DO 9035 NRW = 1,NROWS
      JLEE  = JLE(NRW)
      JTEE  = JTE(NRW)      
           CALL SMOOTH(JLEE,JTEE,4,0.25,S_DIST,CENT_ANGL)
 9035 CONTINUE
C
      DO 9037 J=1,JM
      ALPHA_CENT(J,K) = CENT_ANGL(J)
 9037 CONTINUE
C           
 9039 CONTINUE
C
C   INTERPOLATE TO FIND THE DEVIATION ANGLE OR BLADE EXIT FLOW ANGLE DEPENDING
C   ON  ANGL_TYP(NR) .
      DO 9040 NR = 1,NROWS
C
      JTEE = JTE(NR)
      QSPAN(1) = 0.0
      DO 9041 K=2,KM
           XDIF = X(JTEE,K) - X(JTEE,K-1)
           RDIF = R(JTEE,K) - R(JTEE,K-1)
           QSPAN(K) = QSPAN(K-1) + SQRT(XDIF*XDIF + RDIF*RDIF)
 9041 CONTINUE
           QTOT = QSPAN(KM)
      DO 9042 K=1,KM
           QSPAN(K) = QSPAN(K)/QTOT
 9042 CONTINUE
      DO 9043 K=1,NANGLES(NR)
           FSPAN(K)  = FRAC_SPAN(NR,K)
           EXANG(K)  = EXIT_ANGL(NR,K) 
 9043 CONTINUE
      DO 9044 K=1,KM
      CALL INTP(NANGLES(NR),FSPAN,EXANG,QSPAN(K),ANGL_OUT)
      IF(ANGL_TYP(NR).EQ.'A') EXIT_ANGL(NR,K) = ANGL_OUT*DEGRAD
      IF(ANGL_TYP(NR).EQ.'D') DEVN_ANGL(NR,K) = ANGL_OUT*DEGRAD
 9044 CONTINUE
C
      WRITE(6,*)
      WRITE(6,*) ' ROW NUMBER', NR
           WRITE(6,9047) (ALPHA_CENT(JTEE,K)*RADDEG,K=1,KM)
      IF(ANGL_TYP(NR).EQ.'A')
     &     WRITE(6,9045) (EXIT_ANGL(NR,K)*RADDEG,K=1,KM)
      IF(ANGL_TYP(NR).EQ.'D')
     &     WRITE(6,9046) (DEVN_ANGL(NR,K)*RADDEG,K=1,KM)
 9045 FORMAT(' BLADE EXIT FLOW ANGLE IN DEGREES ',/,(10F10.3))
 9046 FORMAT(' EXIT FLOW ANGLE DEVIATION FROM GRID ANGLE, IN DEGREES.',
     & /,(10F10.3))
 9047 FORMAT(' GRID CENTRE LINE ANGLE IN DEGREES AT BLADE EXIT',
     & /,(10F10.3))
C
C   CONVERT THE EXIT ANGLE TO DEVIATION ANGLE IF  "ANGL_TYP" = A .
C
      IF(ANGL_TYP(NR).EQ.'A') THEN 
           DO 9049 K=1,KMM1
           GRID_ANGL       = ALPHA_CENT(JTEE,K)
           DEVN_ANGL(NR,K) = GRID_ANGL - EXIT_ANGL(NR,K)
 9049      CONTINUE
      END IF      
C
 9040 CONTINUE
C
C   END OF  "IF IM = 2" LOOP
      END IF
C******************************************************************************
C  END TFLOW
C******************************************************************************
C
C************IF NOUT(1)=0 WRITE OUT AND PLOT OUT THE INITIAL GUESS OF THE FLOW FIELD.
      IF(NOUT(1).EQ.0) THEN
           WRITE(6,*) ' WRITING A PLOT FILE OF THE INITIAL GUESS '
           CALL OUTPUT
      END IF
C
C
      WRITE(6,*)
      WRITE(6,*)'******************************************************'
      WRITE(6,*) ' LEAVING SUBROUTINE SETUP,  OK SO FAR '
      WRITE(6,*)'******************************************************'
      WRITE(6,*)
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE NEW_INTPOL
C       ====================
C
C       THIS ROUTINE INTERPOLATES VALUES FOR X  ,  R,  RT_UPP  AND
C       RT_THICK FOR THE BLADE SECTIONS ON THE FINAL GRID. 
C       THE GRID IS ADJUSTED TO ALLOW FOR THE TIP GAP AND NUMBER OF POINTRS IN THE GAP.
C
C       ====================
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION   XINT(MAXKI),YINT(MAXKI),RTUP_INT(MAXKI),
     &            RTTK_INT(MAXKI),FRMOD(MAXKI),QODIST(MAXKI)
C
      DOUBLE PRECISION SUMFR,TIPGAP,TIPSPACE,FAC
C
C
      DO 1 K=1,NSECS_IN
      SMERID(1,K) = 0.0
      DO 2 J=2,JM
      XD = XSURF(J,K)  - XSURF(J-1,K)
      RD = RSURF(J,K)  - RSURF(J-1,K)
    2 SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)
    1 CONTINUE
C
C******************************************************************************
C******************************************************************************
C
      DO 1000 J=1,JM
C
      NR = NROW(J)
      J1 = JSTART(NR)
      JL = JLE(NR)
      JT = JTE(NR)
      JE = JMIX(NR)
C
C    SET THE DEFAULT THAT   FRMOD(K)   = FR(K) IF NO TIP GAP .
C
      DO 55 K=1,KMM1
   55 FRMOD(K) = FR(K)
C
C   DO NOT MODIFY THE GRID AND JUMP TO 75 IF THERE IS NO TIP GAP.
C
      IF(KTIPS(NR).LE.0) GO TO 75
C
C******************************************************************************
C   SET FRMOD(K) TO MODIFY THE GRID IF THERE IS A TIP GAP.
C   THE GAP IS VARIED LINEARLY FROM FRACTIP1 AT THE LE TO FRACTIP2  AT THE TE.
C
          KTIP   = NSECS_IN
          IF(KTIPS(NR).EQ.1) KTIP =1
          FAC = 1.0
          IF(J.GE.J1.AND.J.LT.JL) THEN
          FAC = (SMERID(J,KTIP)-SMERID(J1,KTIP))
     &         /(SMERID(JL,KTIP)-SMERID(J1,KTIP)) 
          TIPGAP = FRACTIP1(NR)
          END IF

          IF(J.GE.JL.AND.J.LE.JT) THEN
               FAC  = 1.0
               FRAC = (SMERID(J,KTIP) - SMERID(JL,KTIP))
     &               /(SMERID(JT,KTIP)- SMERID(JL,KTIP))
               TIPGAP = FRACTIP1(NR) + FRAC*(FRACTIP2(NR)-FRACTIP1(NR)) 
          END IF

          IF(J.GT.JT.AND.J.LE.JE) THEN
          FAC = (SMERID(JE,KTIP)-SMERID(J,KTIP))
     &    /(SMERID(JE,KTIP)- SMERID(JT,KTIP))
          TIPGAP = FRACTIP2(NR)
          END IF
C
      SUMFR = 0.0
      NCELL    = KTIPE(NR) - KTIPS(NR)
      TIPSPACE = TIPGAP/NCELL
      DO 50 K=KTIPS(NR),KTIPE(NR)-1
   50 SUMFR = SUMFR + FR(K)
      DO 60 K=KTIPS(NR),KTIPE(NR)-1
   60 FRMOD(K) =  TIPSPACE*FAC + (1.-FAC)*FR(K)
      DO 70 K= KS1(NR),KS2(NR)
   70 FRMOD(K) = FR(K)*((1.-TIPGAP)/(1.-SUMFR)*FAC + (1.-FAC) )
C
C   END OF SETTING FRMOD(K) FOR THE TIP GAP .
C****************************************************************************
C****************************************************************************
C
   75 CONTINUE
C
C     NOW INTERPOLATE IN THE NSECS_IN UNIFORMLY SPACED STREAM SURFACES TO OBTAIN
C     THE COORDINATES ON   KM   STREAM SURFACES FOR THE FINAL GRID.
C
      QODIST(1)   = 0.0
      DO 10 K=1,NSECS_IN
      XINT(K)     = XSURF(J,K)
      YINT(K)     = RSURF(J,K)
      RTUP_INT(K) = RT_UPP(J,K)
      RTTK_INT(K) = RT_THICK(J,K)
      IF (K.GT.1) THEN
         RD = RSURF(J,K)  -  RSURF(J,K-1)
         XD = XSURF(J,K)  -  XSURF(J,K-1)
	 QODIST(K) = QODIST(K-1) + SQRT(XD*XD+RD*RD)
      END IF
   10 CONTINUE
C
      XARG   = 0.0
      QOSPAN = QODIST(NSECS_IN)
C
      DO 30 K=1,KM
      CALL INTP(NSECS_IN,QODIST,XINT,XARG,X(J,K))
      CALL INTP(NSECS_IN,QODIST,YINT,XARG,R(J,K))
      CALL INTP(NSECS_IN,QODIST,RTUP_INT,XARG,RT_UPP(J,K))
      CALL INTP(NSECS_IN,QODIST,RTTK_INT,XARG,RT_THICK(J,K))
      XARG  = XARG + FRMOD(K)*QOSPAN
   30 CONTINUE
C  
 1000 CONTINUE
C
      RETURN
      END
C
C***********************************************************************
C
      SUBROUTINE OLD_INTPOL
C       ====================
C
C       THIS ROUTINE INTERPOLATES VALUES FOR RT_UPP
C       RT_THICK, X  AND R FOR THE REQUIRED CROSS-SECTIONS
C
C       ====================
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION   XINT(MAXKI),YINT(MAXKI),RTUP_INT(MAXKI),
     &            RTTK_INT(MAXKI),FRMOD(MAXKI),QODIST(MAXKI)
C
      DOUBLE PRECISION SUMFR,TIPGAP,TIPSPACE,FAC
C
      QODIST(1) = 0.0
C
      DO 1 K=1,NOSECT
      SMERID(1,K) = 0.0
      DO 2 J=2,JM
      XD = XSURF(J,K)  - XSURF(J-1,K)
      RD = RSURF(J,K)  - RSURF(J-1,K)
    2 SMERID(J,K) = SMERID(J-1,K) + SQRT(XD*XD+RD*RD)
    1 CONTINUE
C
C******************************************************************************
C
      DO 1000 J=1,JM
      NR = NROW(J)
      J1 = JSTART(NR)
      JL = JLE(NR)
      JT = JTE(NR)
      JE = JMIX(NR)
C
C   DO NOT MODIFY THE GRID AND JUMP TO 75 IF NO TIP GAP.
C
      DO 55 K=1,KMM1
   55 FRMOD(K) = FR(K)
C
      IF(KTIPS(NR).LE.0) GO TO 75
C
C******************************************************************************
C   MODIFY THE GRID IF THERE IS A TIP GAP.
C   THE GAP IS VARIED LINEARLY FROM FRACTIP1 AT THE LE TO FRACTIP2  AT THE TE.
C
          KTIP   = NOSECT
          IF(KTIPS(NR).EQ.1) KTIP =1
          FAC = 1.0
          IF(J.GE.J1.AND.J.LT.JL) THEN
          FAC = (SMERID(J,KTIP) -SMERID(J1,KTIP))
     &         /(SMERID(JL,KTIP)-SMERID(J1,KTIP)) 
          TIPGAP = FRACTIP1(NR)
          END IF

          IF(J.GE.JL.AND.J.LE.JT) THEN
               FAC  = 1.0
               FRAC = (SMERID(J,KTIP) - SMERID(JL,KTIP))
     &               /(SMERID(JT,KTIP)- SMERID(JL,KTIP))
               TIPGAP = FRACTIP1(NR) + FRAC*(FRACTIP2(NR)-FRACTIP1(NR)) 
          END IF

          IF(J.GT.JT.AND.J.LE.JE) THEN
          FAC = (SMERID(JE,KTIP)-SMERID(J,KTIP))
     &    /(SMERID(JE,KTIP)- SMERID(JT,KTIP))
          TIPGAP = FRACTIP2(NR)
          END IF
C
C      WRITE(6,*) 'J,J1,JL,JT,JE, FAC, TIPGAP, RATGAP ', J,J1,JL,JT,JE,
C     &            FAC ,TIPGAP, RATGAP
C
C      CHANGE FR(K) TO FRMOD(K) TO SLIGHTLY ADJUST THE BLADE TIP 
C      POSITION SO THAT THE TIP GAP IS FRACTIP.
C      FIRST SET FRMOD(K).
C
      SUMFR = 0.0
      NCELL    = KTIPE(NR) - KTIPS(NR)
      TIPSPACE = TIPGAP/NCELL
      DO 50 K=KTIPS(NR),KTIPE(NR)-1
   50 SUMFR = SUMFR + FR(K)
      DO 60 K=KTIPS(NR),KTIPE(NR)-1
   60 FRMOD(K) =  TIPSPACE*FAC + (1.-FAC)*FR(K)
      DO 70 K= KS1(NR),KS2(NR)
   70 FRMOD(K) = FR(K)*((1.-TIPGAP)/(1.-SUMFR)*FAC + (1.-FAC) )
C
C****************************************************************************
C
   75 CONTINUE
C
      DO 10 K=1,NOSECT
      XINT(K)     = XSURF(J,K)
      YINT(K)     = RSURF(J,K)
      RTUP_INT(K) = RT_UPP(J,K)
      RTTK_INT(K) = RT_THICK(J,K)
      IF (K.GT.1) THEN
         RD = RSURF(J,K)  - RSURF(J,K-1)
         XD = XSURF(J,K)  -  XSURF(J,K-1)
	 QODIST(K) = QODIST(K-1) + SQRT(XD*XD+RD*RD)
      END IF
   10 CONTINUE
C
      DREF = 0.001*QODIST(NOSECT)
      DREF = DREF*DREF
C
C     FIND THE POINT NEAREST TO THE HUB, L1
C
      L1 = 1
      FLAG1 = 1.
      DO 20 L=1,NOSECT-1
      IF ((RSURF(J,L)-R(J,1))*(RSURF(J,L+1)-R(J,1)).LE.DREF.AND.
     *   (XSURF(J,L)-X(J,1))*(XSURF(J,L+1)-X(J,1)).LE.DREF) THEN
         L1=L
         GOTO 21
      END IF
   20 CONTINUE
      FLAG1 = -1.
      IF(INSURF.NE.2) WRITE(6,111) J
  111 FORMAT(//, ' WARNING !!!! THE FIRST BLADE SECTION IS OUTBOARD OF
     &THE HUB',/,' THE EXTRAPOLATION MAY BE VERY INACCURATE ','J=',I3)
C
C     FIND THE POINT NEAREST TO THE CASING, LM
C
   21 CONTINUE
C
      LM = NOSECT
      FLAG2 = 1.
      DO 22 L=NOSECT,2,-1
      IF ((RSURF(J,L)-R(J,KM))*(RSURF(J,L-1)-R(J,KM)).LE.DREF.AND.
     *   (XSURF(J,L)-X(J,KM))*(XSURF(J,L-1)-X(J,KM)).LE.DREF) THEN
         LM=L
         GOTO 25
      END IF
   22 CONTINUE
C
      FLAG2 = -1.
      IF(INSURF.NE.2) WRITE(6,112) J
  112 FORMAT(//,' WARNING !!!! THE LAST BLADE SECTION IS INBOARD OF THE
     &CASING ',/,' THE EXTRAPOLATION MAY BE VERY INACCURATE ','J= ',I3)
C
   25 CONTINUE
C
      RD = RSURF(J,L1) - R(J,1)
      XD = XSURF(J,L1)   - X(J,1)
      QDIST1 = SQRT(RD*RD+XD*XD)*FLAG1
      RD = RSURF(J,LM) - R(J,KM)
      XD = XSURF(J,LM)   - X(J,KM)
      QDISTM = SQRT(RD*RD+XD*XD)*FLAG2
      QOSPAN = QODIST(LM) - QODIST(L1) - QDIST1 - QDISTM
C
      XARG   = QODIST(L1) + QDIST1
C
C      WRITE(6,*)
C      WRITE(6,*) 'K,   X(J,K),   R(J,K),   RT_UPP(J,K),   RT_THICK(J,K)'
      DO 30 K=1,KM
      IF (K.GT.1) XARG = XARG + QOSPAN*FRMOD(K-1)
      IF(INSURF.NE.2.AND.(XARG.LT.QODIST(1).OR.XARG.GT.QODIST(NOSECT)))
     &      THEN
C
      WRITE(6,*) ' WARNING LINEAR EXTRAPOLATION AT J= ',J,'K= ',K
C
      CALL LININT(NOSECT,QODIST,XINT,XARG,X(J,K))
      CALL LININT(NOSECT,QODIST,YINT,XARG,R(J,K))
      CALL LININT(NOSECT,QODIST,RTUP_INT,XARG,RT_UPP(J,K))
      CALL LININT(NOSECT,QODIST,RTTK_INT,XARG,RT_THICK(J,K))
      ELSE
      CALL INTP(NOSECT,QODIST,XINT,XARG,X(J,K))
      CALL INTP(NOSECT,QODIST,YINT,XARG,R(J,K))
      CALL INTP(NOSECT,QODIST,RTUP_INT,XARG,RT_UPP(J,K))
      CALL INTP(NOSECT,QODIST,RTTK_INT,XARG,RT_THICK(J,K))
      ENDIF
C
C      IF(J.EQ.150) WRITE(6,66) K,X(J,K),R(J,K),RT_UPP(J,K),RT_THICK(J,K)
C   66 FORMAT(I5,4F10.5)
C
   30 CONTINUE
C
      DO 40 K=1,KM
      XSURF(J,K)   = X(J,K)
      RSURF(J,K)   = R(J,K)
   40 CONTINUE
 1000 CONTINUE
C
      RETURN
      END
C
C*****************************************************************************
C*****************************************************************************
C
      SUBROUTINE INTP(N,XN,YN,X,Y)
C
C      THIS SUBROUTINE INTERPOLATES IN THE GIVEN TABLE OF YN AS A
C      FUNCTION OF XN TO FIND THE VALUE OF Y AT THE INPUT VALUE
C      OF X.
C
      DIMENSION XN(N),YN(N)
      SPAN=XN(N)-XN(1)
      Y=0.
      L=1
      NM=N
      IF(N.LT.4) GO TO 8
      NM=4
    4 IF(SPAN.GT.0.0.AND.X.LT.XN(L)) GO TO 5
      IF(SPAN.LT.0.0.AND.X.GT.XN(L)) GO TO 5
      IF(L.EQ.N) GO TO 3
      L=L+1
      GO TO 4
    5 IF(L.GT.2) GO TO 6
      L=1
      GO TO 8
    6 IF(L.NE.N) GO TO 7
    3 L=N-3
      GO TO 8
    7 L=L-2
    8 DO 11 L1=1,NM
      CO=1
      DO 10 L2=1,NM
      IF(L1.EQ.L2) GO TO 9
      TEMP=(X-XN(L+L2-1))/(XN(L+L1-1)-XN(L+L2-1))
      GO TO 10
    9 TEMP=1
   10 CO=CO*TEMP
   11 Y=Y+CO*YN(L+L1-1)
      RETURN
      END
C
C*************************************************************************************
C
       SUBROUTINE LININT(NPOINTS,X,Y,XARG,YANS)
C
C      THIS SUBROUTINE INTERPOLATES IN THE GIVEN TABLE OF YN AS A
C      FUNCTION OF X TO FIND THE VALUE OF Y AT THE INPUT VALUE
C      OF X = XARG.
C
C      THIS VERSION USES LINEAR INTERPOLATION TO AVOID ANY POSSIBLE
C      PROBLEMS WITH OVERSHOOTS OR UNDERSHOOTS.
C
      DIMENSION X(NPOINTS),Y(NPOINTS)
C
      IF (X(1).GT.XARG) THEN
      YANS = Y(1) + (XARG-X(1))*(Y(2)-Y(1))/(X(2)-X(1))
      ELSE
      N=2
   10 CONTINUE
      IF(X(N).GT.XARG) GO TO 20
      N=N+1
      IF(N.GT.NPOINTS) GO TO 30
      GO TO 10
   20 YANS = Y(N) + (XARG-X(N))*(Y(N-1)-Y(N))/(X(N-1)-X(N))
      GO TO 40
   30 YANS = Y(NPOINTS) + (XARG-X(NPOINTS))*(Y(NPOINTS)-Y(NPOINTS-1))/
     & (X(NPOINTS)-X(NPOINTS-1))
   40 CONTINUE
C  
      ENDIF
C
      RETURN
      END
C
C******************************************************************************
C 
      SUBROUTINE SUMFLX(BLADE_FLOW,SUMPO,SUMTO,SUMRVT,SUM_ENTPY,
     &                  SUMTSTAT,SUMPSTAT) 
C
C   THIS SUBROUTINE MASS AVERAGES SOME FLOW QUANTITIES AT EVERY  "J"  STATION.
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION BLADE_FLOW(JD),SUM_ENTPY(JD),SUMPO(JD),SUMTO(JD),
     &          SUMRVT(JD),SUMTSTAT(JD),SUMPSTAT(JD)
C
      DO 100 J=1,JM
      BLADE_FLOW(J)  = 0.0
      SUMPO(J)      = 0.0
      SUMTO(J)      = 0.0
      SUMRVT(J)     = 0.0
      SUM_ENTPY(J)  = 0.0
      SUMTSTAT(J)   = 0.0
      SUMPSTAT(J)   = 0.0
      DO 110 K=1,KMM1
      DO 120 I=1,IMM1
      DFLOW        = -FLOWX(I,J,K)*NBLADE(J)
      BLADE_FLOW(J) = BLADE_FLOW(J) + DFLOW
      DFLOW4       = 0.25*DFLOW
      SUMPO(J)     = SUMPO(J) + DFLOW4*(TEMP2(I,J,K)+TEMP2(I+1,J,K)
     &             + TEMP2(I+1,J,K+1)+TEMP2(I,J,K+1))
      SUMTO(J)     = SUMTO(J) + DFLOW4*(TEMP1(I,J,K)+TEMP1(I+1,J,K)
     &             + TEMP1(I+1,J,K+1)+TEMP1(I,J,K+1))
      SUMRVT(J)    = SUMRVT(J)+ DFLOW4*(TEMP3(I,J,K)+TEMP3(I+1,J,K)
     &             + TEMP3(I+1,J,K+1)+TEMP3(I,J,K+1))
      SUM_ENTPY(J) = SUM_ENTPY(J) + DFLOW4*(TEMP4(I,J,K)+TEMP4(I+1,J,K)
     &             + TEMP4(I+1,J,K+1)+TEMP4(I,J,K+1))
      SUMTSTAT(J)  = SUMTSTAT(J) + DFLOW4*(STORE2(I,J,K)+STORE2(I+1,J,K)
     &             + STORE2(I+1,J,K+1)+STORE2(I,J,K+1))
      SUMPSTAT(J)  = SUMPSTAT(J)  + DFLOW4*(P(I,J,K)+P(I+1,J,K)
     &             + P(I+1,J,K+1) + P(I,J,K+1))
  120 CONTINUE
  110 CONTINUE
  100 CONTINUE
      RETURN
      END
C
C******************************************************************************
C
      SUBROUTINE STEPUP(RELAX)
C
C          THIS SUBROUTINE CHANGES THE TIMESTEP IN PROPORTION TO THE LOCAL
C          MACH NUMBER. THE CHANGES ARE RELAXED BY THE FACTOR  "RELAX" .
C
      INCLUDE  'commall-open-18.3'
C
C
      TLIM   = 0.1*TO1(KMID)
      VLIM   = 0.01*VM1(KMID)*VM1(KMID)
      RELAX1 = 1.0 - RELAX
C******************************************************************************
C     CALCULATE THE LOCAL SPEED OF SOUND AND SET THE FACTOR  C/(V+C)
C
      DO 1100 K=1,KM
      DO 1100 J=1,JM
      DO 1100 I=1,IM
      VSQ = VX(I,J,K)*VX(I,J,K)+VT(I,J,K)*VT(I,J,K)+VR(I,J,K)*VR(I,J,K)
      HSTAT  = HO(I,J,K) - 0.5*VSQ
C
      IF(IFGAS.EQ.0) THEN
           TSTATIC = HSTAT*RCP
           GAMNOW  = GA
      ELSE
           TSTATIC = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
           CPNOW   = CP1 + CP2*(TSTATIC-TREF) 
     &             + CP3*(TSTATIC-TREF)*(TSTATIC-TREF)
           GAMNOW  = CPNOW/(CPNOW-RGAS)
      END IF
C
      IF(TSTATIC.LT.TLIM) TSTATIC = TLIM
      WTREL     = WT(I,J,K)
      WSQ       = VSQ - (VT(I,J,K)*VT(I,J,K) - WTREL*WTREL)
      IF(WSQ.LT.VLIM) WSQ = VLIM
C
      IF(ITIMST.GE.5) THEN
           V_SONIC = VSOUND
      ELSE
           V_SONIC = SQRT(GAMNOW*RGAS*TSTATIC)
      END IF
C
      VPLUSC    =  SQRT(WSQ) + V_SONIC
C
 1100 TEMP1(I,J,K) =  VSOUND/VPLUSC
C
C       TEMP1 IS USED AS A TEMPORARY STORE FOR  VSOUND/(W + VSOUND)
C       WHERE VSOUND IS THE STAGNATION SPEED OF SOUND AND W ISTHE RELATIVE VELOCITY.
C
C******************************************************************************
C  USE  BSTEP AND THE LOCAL VELOCITIES TO SET THE TIME STEP/VOLUME   STEP(I,J,K).
C
      DO 10 I=1,IMM1
      IP1=I+1
      DO 10 K=1,KMM1
      KP1=K+1
      AMACHP =  TEMP1(I,1,K)   + TEMP1(IP1,1,K) +
     &          TEMP1(I,1,KP1) + TEMP1(IP1,1,KP1)
      DO 10 J=  2,JM
      AMACHL =  TEMP1(I,J,K)   + TEMP1(I,J,KP1) +
     &          TEMP1(IP1,J,K) + TEMP1(IP1,J,KP1)
      AVMACH =  0.125*(AMACHL  + AMACHP)
      AMACHP =  AMACHL
      STEPNEW      =  BSTEP(I,J,K)*AVMACH
      STEP(I,J,K)  =  RELAX*STEPNEW  +  RELAX1*STEP(I,J,K)
      RSTEP(I,J,K) =  STEP(I,J,K)/BSTEP(I,J,K)
   10 CONTINUE
C
C******************************************************************************
C      SET THE TIME STEP AT THE MIXING PLANE, ie AT J = JMIX + 1..
C
C    JDD MODIFIED THIS  MARCH 2015. SO STEP(JMIX+1) = STEP(JMIX+2).
C    WHICH GIVES BETTER RESULTS
C
      DO 101 NR=1,NRWSM1
      J   = JMIX(NR)
      JP1 = J+1
      JP2 = J+2
      DO 201 K=1,KMM1
      DO 201 I=1,IMM1
      STEP(I,JP1,K)   = STEP(I,JP2,K)
      RSTEP(I,JP1,K)  = 0.0
  201 CONTINUE
  101 CONTINUE
C
C******************************************************************************
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE SETFLO
C
C      THIS SUBROUTINE FORCES THE MASS FLOW TOWARDS AN INPUT VALUE IN IN_FLOW=3
C      OR TOWARDS THE AVERAGE VALUE IF IN_FLOW=2.
C      IT DOES THIS BY MEANS OF A BODY FORCE WHICH GENERATES LOSS IF IN_FLOW=3
C      BUT IN_FLOW =2 SHOULD GIVE IMPROVED CONVERGENCE WITH NO LOSS GENERATION.
C
      INCLUDE  'commall-open-18.3'
C
      DO 100 J=1,JM
      FLOW(J)=0.0
      DO 100 K=1,KMM1
      DO 100 I=1,IMM1
      FLOW(J) = FLOW(J)-FLOWX(I,J,K)*NBLADE(J)/NBLADE(1)
  100 CONTINUE
C
      AVFLOW = FLOWIN/NBLADE(1)
      IF(IN_FLOW.EQ.3) GO TO 300
      AVFLOW=0.0
      DO 200 J=1,JM
      AVFLOW=AVFLOW+FLOW(J)
  200 CONTINUE
      AVFLOW=AVFLOW/JM
  300 CONTINUE
C
      DO 400 J=1,JM
      ROVMER=SQRT(ROVX(IMID,J,KMID)*ROVX(IMID,J,KMID)+ROVR(IMID,J,KMID)
     &  *ROVR(IMID,J,KMID))
      DELTA = (AVFLOW/FLOW(J) -1.0)*ROVMER
      DELTA = DELTA*RFLOW
      DO 400 I=1,IM
      DO 400 K=1,KM
      ROVMER = SQRT(ROVX(I,J,K)*ROVX(I,J,K)+ROVR(I,J,K)*ROVR(I,J,K))
C
C  ADDED   3/9/90 may not be valid for radial flow machines
C
      IF(ROVX(I,J,K).LT.0.0)  ROVMER = -ROVMER
C
      ROVX(I,J,K) = ROVX(I,J,K) + ROVX(I,J,K)/ROVMER*DELTA
      ROVR(I,J,K) = ROVR(I,J,K) + ROVR(I,J,K)/ROVMER*DELTA
      RO_WT       = ROVT(I,J,K) - UBLADE(J,K)*RO(I,J,K)
      ROVT(I,J,K) = ROVT(I,J,K) + RO_WT/ROVMER*DELTA
C
  400 CONTINUE
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE LOSS
C
C            THIS SUBROUTINE COMPUTES A BODY FORCE BASED ON WALL FUNCTIONS FOR THE
C            SURFACE SHEAR STRESS AND A MIXING LENGTH MODEL OF EDDY VISCOSITY.
C
C            THE BODY FORCE IS CALCULATED BY MAKING A THIN SHEAR LAYER APPROXIMATION
C            TO THE N-S EQUATIONS.
C
C            THE METHOD USES WALL FUNCTIONS TO CALCULATE THE SURFACE SHEAR STRESSES
C            IF "YPLUSWALL" < 5 , AND  EVALUATES THE SURFACE SHEAR STRESSES BY ASSUMING
C            THAT THE WALL GRID POINT IS AT "YPLUSWALL" IF YPLUSWALL > 5.
C
C            THE FULL ENERGY EQUATION INCLUDING HEAT CONDUCTION AND VISCOUS WORK IS
C            SOLVED INSTEADY OF ASSUMING THAT THEY CANCEL AS IN EARLIER VERSIONS.
C            ALL SOLID SURFACES ARE ASSUMED TO BE ADIABATIC - NO HEAT FLUX.
C
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION VXAVG(MAXKI),VRAVG(MAXKI),WTAVG(MAXKI),WABS(MAXKI),
     &          ROAVG(MAXKI),XSTRES(MAXKI),RAVG(MAXKI),RSTRES(MAXKI),
     &          TSTRES(MAXKI),AREA(MAXKI),WBOUND(MAXKI),WTB(MAXKI),
     &          VTAVG(MAXKI),TAVG(MAXKI),WVISC(MAXKI),QFLOW(MAXKI),
     &          TEMPP(JD)
C
C
C      CALCULATE THE VISCOSITY OVER THE FIRST QUARTER OF THE STEPS.
C      THEN HOLD IT CONSTANT FOR THE REMAINDER OF THE STEPS.
C
      IF(NSTEP.EQ.1.OR.NSTEP.LT.NSTEPS_MAX/4) THEN
C
      IF(REYNO.GT.100.) THEN
           J1=JLE(1)
           J2=JTE(1)
           IF(JLE(1).GT.JM) J1=1
           IF(JTE(1).GT.JM) J2=JM
           XCHORD = SMERID(J2,KMID)-SMERID(J1,KMID)
           ROW2   = SQRT(ROVX(IMID,J2,KMID)*ROVX(IMID,J2,KMID)
     &            + ROVR(IMID,J2,KMID)*ROVR(IMID,J2,KMID)
     &            + ROWT(IMID,J2,KMID)*ROWT(IMID,J2,KMID))
           VISLAM  = XCHORD*ROW2/REYNO
      END IF
C
      IF(REYNO.GT.0.0.AND.REYNO.LT.99.99) THEN
            VISLAM = REYNO/100000.
      END IF
C
      IF(REYNO.LT.0.0) VISLAM = -REYNO*1.0E-5
C
      TCOND  = CP*VISLAM/PRANDTL
      FTCOND = TCOND/VISLAM
C
C   END OF PART ONLY USED FOR THE FIRST QUARTER OF THE STEPS.
      ENDIF
C******************************************************************************
C    SAVE THE VISCOSITY IT IS ONLY USED TO CALCULATE AND WRITE OUT THE REYNOLDS NUMBER. 
      DO NRW =1,NROWS
      VISCOSY(NRW)  = VISLAM
      END DO
C******************************************************************************
C******************************************************************************
C******************************************************************************
C     DO UP TO STATEMENT 25 ONLY ON THE FIRST CALL TO THE SUBROUTINE
C
      IF(NSTEP.GT.NLOS) GO TO 25
C
C     EVALUATE SOME PARAMETERS NEEDED TO CALCULATE THE VISCOUS STRESSES ON THE
C     STREAMWISE ( K = CONSTANT) SURFACES.
C
      DO 27 NRW = 1,NROWS
C
      SUMF = 0.0
      DO 26 I  = 2,IMM1
      SUMF          = SUMF+FP(I-1)
      XLIM          = XLLIM_I1(NRW)*(1.0-SUMF) + XLLIM_IM(NRW)*SUMF
      FPITCH        = SUMF*(1.-SUMF)
      IF(FPITCH.GT.XLIM) FPITCH = XLIM
      DF            = FP(I)+FP(I-1)
      FILAM(NRW,I)  = FP(I)/DF
      FITURB(NRW,I) = 0.16*FPITCH*FPITCH/(DF*DF)
      FIWAKE(NRW,I) = 0.16*XLLIM_DWN(NRW)*XLLIM_DWN(NRW)/(DF*DF)
      FIUP(NRW,I)   = 0.16*XLLIM_UP(NRW)*XLLIM_UP(NRW)/(DF*DF)
   26 CONTINUE
      DF            = FP(1) + FP(IMM1)
      FILAM(NRW,1)  = FP(1)/DF
      FIWAKE(NRW,1) = 0.16*XLLIM_DWN(NRW)*XLLIM_DWN(NRW)/(DF*DF)
      FIUP(NRW,1)   = 0.16*XLLIM_UP(NRW)*XLLIM_UP(NRW)/(DF*DF)
      FITURB(NRW,1) = FITURB(NRW,2)
      FILAM(NRW,IM) = FP(IMM1)/DF
      FIWAKE(NRW,IM)= 0.16*XLLIM_DWN(NRW)*XLLIM_DWN(NRW)/(DF*DF)
      FIUP(NRW,IM)  = 0.16*XLLIM_UP(NRW)*XLLIM_UP(NRW)/(DF*DF)
      FITURB(NRW,IM)= FITURB(NRW,IMM1)
C
C  Q3D
      IF(KM.EQ.2) GO TO 16
C  END Q3D
C
C     FIND THE MIXING LENGTH LIMITS ON THE HUB AND CASING
C
      AVGSPAN = 0.0
      JMID    = 0.5*(JLE(NRW)+JTE(NRW))
      DO 13 K=2,KM
      XDIF    = X(JMID,K)-X(JMID,K-1)
      RDIF    = R(JMID,K)-R(JMID,K-1)
      AVGSPAN = AVGSPAN + SQRT(XDIF*XDIF+RDIF*RDIF)
   13 CONTINUE
      AVGPIT  = 2*3.1415926*R(JMID,KMID)/NBLADE(JMID)
      XLIMH   = XLLIM_K1(NRW)*AVGPIT/AVGSPAN
      XLIMT   = XLLIM_KM(NRW)*AVGPIT/AVGSPAN
C
C     EVALUATE PARAMETERS NEEDED TO CALCULATE THE VISCOUS STRESSES ON THE
C     BLADEWISE (I = CONSTANT) SURFACES.
C
      SUMF = 0.0
      DO 15 K=2,KMM1
      SUMF     = SUMF + FR(K-1)
      XLIM     = XLIMH*(1.0-SUMF) + XLIMT*SUMF
      FCSPAN   = SUMF*(1.-SUMF)
      IF(FCSPAN.GT.XLIM) FCSPAN = XLIM
      DF       = FR(K) + FR(K-1)
      FKLAM(NRW,K)  =    FR(K)/DF
      FKTURB(NRW,K) =    0.16*FCSPAN*FCSPAN/(DF*DF)
   15 CONTINUE
C
   16 CONTINUE
C
   27 CONTINUE
C
C     END OF THE SETUP USED ONLY ON THE FIRST CALL TO THIS SUBROUTINE.
C
   25 CONTINUE
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C******************************************************************************
C    EVALUATE AND SMOOTH THE PRESSURE GRADIENTS IF YPLUSWALL IS < -10.0.
C
      IF(YPLUSWALL.LT.-10.0) CALL SET_PWALLGRAD
C
C********************************************************************************
C********************************************************************************
C      FIRST WORK OUT THE VISCOUS STRESSES ON THE STREAMWISE (K = CONSTANT) FACES OF THE
C      ELEMENTS IN THE DO 50 J LOOP
C
C  Q3D   
      IF(KM.EQ.2) GO TO 555
C  END Q3D
C
      DO 50 J=2,JM
C
      NRW    = NROW(J)
      J1     = JSTART(NRW)
      JROW   = J - J1 + 1
      JTRHUB = JTRAN_K1(NRW)
      JTRTIP = JTRAN_KM(NRW)
      JLEDGE = JLE(NRW)
      JTEDGE = JTE(NRW)
      WREL   = WRAD(J)
C
C     EVALUATE THE AVERAGE VELOCITIES ETC ON THE STREAMWISE FACES
C     OF THE ELEMENTS.
      DO 40 I=1,IMM1
C
C
      DO 30 K=1,KM
      AREA(K)  = SQRT(ASX(I,J,K)*ASX(I,J,K)+ASR(I,J,K)*ASR(I,J,K))
      VXAVG(K) = 0.25*(VX(I,J,K)+VX(I,J-1,K)+VX(I+1,J,K)+VX(I+1,J-1,K))
      VRAVG(K) = 0.25*(VR(I,J,K)+VR(I,J-1,K)+VR(I+1,J,K)+VR(I+1,J-1,K))
      VTAVG(K) = 0.25*(VT(I,J,K)+VT(I,J-1,K)+VT(I+1,J,K)+VT(I+1,J-1,K))
      ROAVG(K) = 0.25*(RO(I,J,K)+RO(I,J-1,K)+RO(I+1,J,K)+RO(I+1,J-1,K))
      TAVG(K)  = 0.25*(T_STATIC(I,J,K)   + T_STATIC(I,J-1,K)
     &               + T_STATIC(I+1,J,K) + T_STATIC(I+1,J-1,K))
      RAVG(K)  = 0.5*(R(J,K) + R(J-1,K))
      WTAVG(K) = VTAVG(K) - WREL*RAVG(K)
      WABSQ    = VXAVG(K)*VXAVG(K)+VRAVG(K)*VRAVG(K)+WTAVG(K)*WTAVG(K)
      WABS(K)  = SQRT(WABSQ)
      WTB(K)   = VTAVG(K) - WHUB(J)*RAVG(K)
      IF(K.GT.KMID) WTB(K) = VTAVG(K) - WTIP(J)*RAVG(K)
      WBOUND(K)= SQRT(VXAVG(K)*VXAVG(K)+VRAVG(K)*VRAVG(K)+WTB(K)*WTB(K))
C
   30 CONTINUE
C
C     EVALUATE THE VISCOSITY FROM A POWER LAW IF REYNO IS NEGATIVE
C
      IF(REYNO.LT.0.0001) THEN
      VISLAM = (ABS(REYNO)/100000.0) * (TAVG(KMID)/288.0)**0.62
      TCOND  = CP*VISLAM/PRANDTL
      FTCOND = TCOND/VISLAM
      END IF
C
C    SAVE THE VISCOSITY IT IS ONLY USED TO CALCULATE AND WRITE OUT THE REYNOLDS NUMBER.
      IF(I.EQ.IMID.AND.J.EQ.JTEDGE) VISCOSY(NRW) = VISLAM
C
C     CALCULATE THE REYNOLDS NUMBER BASED ON FLOW CONDITIONS AT THE TRAILING EDGE.
C
      IF(I.EQ.IMID.AND.J.EQ.JTEDGE) THEN
           ROVEXIT(NRW)  = ROAVG(KMID)*WABS(KMID)
           VISCOSY(NRW)  = VISLAM
           REYNOLDS(NRW) = CHORD(NRW)*ROVEXIT(NRW)/VISLAM
      END IF
C
C************************************************************************
C      EVALUATE THE VISCOUS STRESSES ON THE HUB AND CASING.
C
C************************************************************************
C      FIRST THE HUB
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           PERPK1    =  VOL(I,J,1)/AREA(1)
           YPLUS_OLD =  YPLUS_K1(I,J)
C
         CALL WALLFUN(I,J,1,1,PERPK1,DPDS_CELL(I,J,1),ROAVG(1),
     &                TWALLK1,YPLUS_OLD,WBOUND(1),YPLUS_NEW)
           YPLUS_K1(I,J) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 345

      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
C
      IF(YPLUSWALL.LT.5.) THEN
C
      PERP   =  VOL(I,J,1)/AREA(1)
      RE     =  PERP*ROAVG(2)*WBOUND(2)/VISLAM
      RELOG  =  1.0/ALOG(RE)
C
C    ALLOW FOR ROUGHNESS
C
      ROUGH = ROUGH_H(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.PERP) ROUGH = PERP
           REK   = RE*ROUGH/PERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C     ALLOW FOR ROUGHNESS IF ROUGH > 1.0E-7 .
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLK1 = 0.5*CF*ROAVG(2)*WBOUND(2)*WBOUND(2)
           VSTAR   = SQRT(TWALLK1/ROAVG(2))
           PLUSK   = ROUGH*VSTAR*ROAVG(2)/VISLAM
           RELOG10 = LOG10(PERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
           TWALLK1    = 0.5*CF*ROAVG(2)*WBOUND(2)*WBOUND(2)
C
      ELSE
C
C   IF YPLUSWALL > 5 USE THE YPLUSWALL VALUE TO CALCULATE THE SKIN FRICTION
C
           TWALLK1 = CFWALL*ROAVG(1)*WBOUND(1)*WBOUND(1)
C
      ENDIF
C
  345 CONTINUE
C
C************************************************************************
C****************************************************************************
C    NEXT ON THE CASING
C
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           PERPKM     =  VOL(I,J,KMM1)/AREA(KM)
           YPLUS_OLD  =  YPLUS_KM(I,J)
C
        CALL WALLFUN(I,J,KMM1,KM,PERPKM,DPDS_CELL(I,J,KMM1),ROAVG(KMM1),
     &              TWALLKM,YPLUS_OLD,WBOUND(KMM1),YPLUS_NEW)
           YPLUS_KM(I,J) = AMIN1(1000.0,YPLUS_NEW)
      GO TO 355
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
C
      IF(YPLUSWALL.LT.5.) THEN
C
      PERP    = VOL(I,J,KMM1)/AREA(KM)
      RE      = PERP*ROAVG(KMM1)*WBOUND(KMM1)/VISLAM
      RELOG   = 1.0/ALOG(RE)
C
C    ALLOW FOR ROUGHNESS
C
      ROUGH = ROUGH_T(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.PERP) ROUGH = PERP
           REK   = RE*ROUGH/PERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF 
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C     ALLOW FOR ROUGHNESS IF ROUGH > 1.0E-7 .
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLKM = 0.5*CF*ROAVG(KMM1)*WBOUND(KMM1)*WBOUND(KMM1)
           VSTAR   = SQRT(TWALLKM/ROAVG(KMM1))
           PLUSK   = ROUGH*VSTAR*ROAVG(KMM1)/VISLAM
           RELOG10 = LOG10(PERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
           TWALLKM  = 0.5*CF*ROAVG(KM M1)*WBOUND(KMM1)*WBOUND(KMM1)
C
      ELSE
C
C   IF YPLUSWALL > 5 USE THE YPLUSWALL VALUE TO CALCULATE THE SKIN FRICTION
C.
           TWALLKM = CFWALL*ROAVG(KM)*WBOUND(KM)*WBOUND(KM)
C
      ENDIF
C
  355 CONTINUE
C
C******************************************************************************
C  JDD ADDITION. JAN/14.  EVALUATE YPLUS ON THE ENDWALLS.
C
      VSTAR     = SQRT(TWALLK1/ROAVG(1))
      DO K=1,KMID
      PERP          = SQRT(DWALLSQ(I,J,K))
      YPLSK1        = VSTAR*PERP*ROAVG(1)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(YPLSK1,1000.0)
      END DO
C
C
      VSTAR     = SQRT(TWALLKM/ROAVG(KM))
      DO K = KMID+1,KMM1
      PERP      = SQRT(DWALLSQ(I,J,K))
      YPLS      = VSTAR*PERP*ROAVG(KM)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(YPLS,1000.0)
      END DO
C
C   END OF JAN/14 ADDITION
C******************************************************************************
C     CALCULATE THE FRICTIONAL STRESSES ON THE HUB
      FMULT      = TWALLK1*AREA(1)/WBOUND(1)
      IF(IBOUND.EQ.1.OR.IBOUND.GT.2) FMULT = 0.0
      XSTRES(1)  = FMULT*VXAVG(1)
      RSTRES(1)  = FMULT*VRAVG(1)
      TSTRES(1)  = FMULT*WTB(1)
C     CALCULATE THE FRICTIONAL STRESSES ON THE CASING
      FMULT      = -TWALLKM*AREA(KM)/WBOUND(KM)
      IF(IBOUND.GE.2) FMULT = 0.0
      XSTRES(KM) = FMULT*VXAVG(KM)
      RSTRES(KM) = FMULT*VRAVG(KM)
      TSTRES(KM) = FMULT*WTB(KM)
C  CALCULATE THE SHEAR WORK AND HEAT FLOW ON THE ENDWALLS.
          WVISC(1)   = TSTRES(1)*WHUB(J)*RAVG(1)
          QFLOW(1)   = 0.0
          WVISC(KM)  = TSTRES(KM)*WTIP(J)*RAVG(KM)
          QFLOW(KM)  = 0.0
C******************************************************************************
C******************************************************************************
C      EVALUATE THE VISCOUS STRESSES ON THE STREAMWISE FACES OF THE
C      ELEMENTS AWAY FROM THE END WALLS.
C
C   JDD ADDED AND CHANGED NEXT SECTION  JAN/14.
C    to reduce the turbulent viscosity in the sublayer and buffer region.
      DO 35 K=2,KMM1
      IF(Y_PLUS(I,J,K).LE.YPLAM) FYPLUS = 0.0
      IF(Y_PLUS(I,J,K).GT.YPLAM)  THEN
               XFAC = (Y_PLUS(I,J,K) - YPLAM)/(YPTURB- YPLAM)
               IF(XFAC.GT.1.0) XFAC = 1.0
               FYPLUS = XFAC*XFAC*(3.0 - 2.0*XFAC)
      END IF
C    
      DPERP   = VOL(I,J,K)/AREA(K)
      VLAM(K) = FKLAM(NRW,K)*VISLAM/DPERP
      VTURB(K)= FKTURB(NRW,K)*ROAVG(K)*ABS(WABS(K+1)-WABS(K-1))*FMIXUP
      VTURB(K)= VTURB(K)*FYPLUS
   35 CONTINUE
C     END JDD JAN/14 ADDITION
C******************************************************************************
C******************************************************************************
C     CHECK FOR TRANSITION
C     FIRST ON THE HUB.
      IF(JROW.LT.JTRHUB) GO TO 34
      RMAX = 0.0
      DO 37 K=2,KMID
      RATVIS = VTURB(K)/VLAM(K)
      IF(RATVIS.LT.RMAX) GO TO 37
      RMAX = RATVIS
   37 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 38
   34 CONTINUE
      DO 39 K = 2,KMID
   39 VTURB(K) = 0.0
   38 CONTINUE
C    NEXT ON THE CASING
      IF(JROW.LT.JTRTIP) GO TO 46
      RMAX = 0.0
      DO 47 K = KMID,KMM1
      RATVIS = VTURB(K)/VLAM(K)
      IF(RATVIS.LT.RMAX) GO TO 47
      RMAX = RATVIS
   47 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 49
   46 CONTINUE
      DO 48 K=KMID,KMM1
   48 VTURB(K) = 0.0
   49 CONTINUE
C******************************************************************************
C     CALCULATE THE STRESSES ON THE STREAMWISE FACES OF THE ELEMENTS.
      DO 41 K=2,KMM1
      FMULT     = AREA(K)*(VLAM(K)+VTURB(K))
      XSTRES(K) = FMULT*(VXAVG(K+1)-VXAVG(K-1))
      RSTRES(K) = FMULT*(VRAVG(K+1)-VRAVG(K-1))
      VISTOT    = VISLAM*(1.0 + VTURB(K)/VLAM(K))
C      VISC_RAT(I,J,K) = VTURB(K)/VLAM(K)
      TSTRES(K) = FMULT*(VTAVG(K+1)-VTAVG(K-1))
     &          - VISTOT*AREA(K)*VTAVG(K)/RAVG(K)
      QFLOW(K)  =  FMULT*FTCOND*(TAVG(K+1)-TAVG(K-1))
      WVISC(K)  =  XSTRES(K)*VXAVG(K) + TSTRES(K)*VTAVG(K)
     &          +  RSTRES(K)*VRAVG(K)
   41 CONTINUE     
C******************************************************************************
C      FORM THE COMPONENT OF VISCOUS FORCE DUE TO THE DIFFERENCE OF THE
C      STRESSES ON ADJACENT STREAMWISE SURFACES.
C
      RF_VIS1 = 1.-RF_VIS
      DO 45 K=1,KMM1
      XFORCE(I,J,K)=RF_VIS1*XFORCE(I,J,K)+RF_VIS*(XSTRES(K)-XSTRES(K+1))
      RFORCE(I,J,K)=RF_VIS1*RFORCE(I,J,K)+RF_VIS*(RSTRES(K)-RSTRES(K+1))
      TFORCE(I,J,K)=RF_VIS1*TFORCE(I,J,K) +
     &              RF_VIS*(TSTRES(K)*RAVG(K)-TSTRES(K+1)*RAVG(K+1))
      QSOURCE(I,J,K) = RF_VIS1*QSOURCE(I,J,K) 
     &     + RF_VIS*(QFLOW(K)  - QFLOW(K+1) + WVISC(K) - WVISC(K+1))
C
   45 CONTINUE
C
C************************************************************************
C     WORK OUT THE ENTROPY GENERATION RATE PER UNIT VOLUME ONLY IF IFEND = 1 OR IFPRINT = 1.
C     FIRST ON THE STREAMWISE SURFACES.
C
      IF(IFPRINT.EQ.1.OR.IFEND.EQ.1) THEN
      DO 57 K = 1,KMM1
      XGEN  = (XSTRES(K)+XSTRES(K+1))*(VXAVG(K+1) -VXAVG(K))
      RGEN  = (RSTRES(K)+RSTRES(K+1))*(VRAVG(K+1) -VRAVG(K))
      TGEN  = (TSTRES(K)+TSTRES(K+1))*(WTAVG(K+1) -WTAVG(K))
      TAVGG = (TAVG(K) + TAVG(K+1))
      QGEN  = 2.0*(QFLOW(K) + QFLOW(K+1))*(TAVG(K+1) - TAVG(K))/TAVGG
      SGEN(I,J,K) = (XGEN + RGEN + TGEN + QGEN)/VOL(I,J,K)/TAVGG
   57 CONTINUE
      END IF
C     END OF ENTROPY GENERATION RATE CALCULATIION ON STREAMWISE SURFACES
C**************************************************************************
C     END OF I LOOP
   40 CONTINUE
C     END OF J LOOP
   50 CONTINUE
C   END OF CALCULATING THE VISCOUS STRESSES ON THE STREAMWISE SURFACES
  555 CONTINUE
C
C********************************************************************************
C********************************************************************************
C********************************************************************************
C      NOW WORK OUT THE VISCOUS STRESSES ON THE BLADES AND BLADEWISE
C      SURFACES IN THE DO 100 J  LOOP.

      DO 100 J=2,JM
C
      NRW    = NROW(J)
      J1     = JSTART(NRW)
      JROW   = J - J1 + 1
      JTRLOW = JTRAN_I1(NRW)
      JTRUP  = JTRAN_IM(NRW)
      JLEDGE = JLE(NRW)
      JTEDGE = JTE(NRW)
      WREL   = WRAD(J)
C
C      FIRST EVALUATE THE AVERAGE VELOCITIES ETC ON THE BLADEWISE SURFACES.
C
      DO 90 K=1,KMM1
      RAVG(K)  = RAVG_CELL(J,K)
      DO 60 I=1,IM
      AREA(I)  = SQRT(ABX(I,J,K)*ABX(I,J,K) + ABR(I,J,K)*ABR(I,J,K)
     &              +ABT(J,K)*ABT(J,K))
      VXAVG(I) = 0.25*(VX(I,J,K)+VX(I,J-1,K)+VX(I,J-1,K+1)+VX(I,J,K+1))
      VRAVG(I) = 0.25*(VR(I,J,K)+VR(I,J-1,K)+VR(I,J-1,K+1)+VR(I,J,K+1))
      VTAVG(I) = 0.25*(VT(I,J,K)+VT(I,J-1,K)+VT(I,J-1,K+1)+VT(I,J,K+1))      
      WTAVG(I) = 0.25*(WT(I,J,K)+WT(I,J-1,K)+WT(I,J-1,K+1)+WT(I,J,K+1))
      ROAVG(I) = 0.25*(RO(I,J,K)+RO(I,J-1,K)+RO(I,J-1,K+1)+RO(I,J,K+1))
      TAVG(I)  = 0.25*(T_STATIC(I,J,K)     + T_STATIC(I,J-1,K)
     &               + T_STATIC(I,J-1,K+1) + T_STATIC(I,J,K+1))
      WABSQ    = VXAVG(I)*VXAVG(I)+VRAVG(I)*VRAVG(I)+WTAVG(I)*WTAVG(I)
      WABS(I)  = SQRT(WABSQ)
      WBOUND(I) = WABS(I)
   60 CONTINUE
C
C**********************************************************************************
C     EVALUATE THE VISCOSITY FROM A POWER LAW IF REYNO IS NEGATIVE
C
      IF(REYNO.LT.0.0001) THEN
      VISLAM = (ABS(REYNO)/100000.0) * (TAVG(IMID)/288.0)**0.62
      TCOND  = CP*VISLAM/PRANDTL
      FTCOND = TCOND/VISLAM
      END IF
C
C**********************************************************************************
C      GO TO 65 TO WORK OUT VISCOUS FORCES BEFORE THE LEADING EDGE AND IN THE WAKES.
C    
      IF(J.LE.JLEDGE.OR.J.GT.JTEDGE) GO TO 65
C
C**********************************************************************************
C**********************************************************************************
C      NOW WORK OUT THE STRESSES ON THE BLADE SURFACES.
C
C      SKIP THE TIP GAP
      IF((KTIPS(NRW).GT.0).AND.(K.GE.KTIPS(NRW)).AND.(K.LT.KTIPE(NRW))) 
     &    GO TO 51
C
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           PERPI1    =  VOL(1,J,K)/AREA(1)
           YPLUS_OLD =  YPLUS_I1(J,K)
C
         CALL WALLFUN(1,J,K,1,PERPI1,DPDS_CELL(1,J,K),ROAVG(1),
     &                TWALLI1,YPLUS_OLD,WBOUND(2),YPLUS_NEW)
          YPLUS_I1(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 365
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
      IF(YPLUSWALL.LT.5.) THEN
C
      PERP   =  VOL(1,J,K)/AREA(1)
      RE     =  ROAVG(2)*WABS(2)*PERP/VISLAM
      RELOG  =  1.0/ALOG(RE)
C
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_L(NRW)
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.PERP) ROUGH = PERP
           REK   = RE*ROUGH/PERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C     ALLOW FOR ROUGHNESS IF ROUGH > 1.0E-7
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLI1    = 0.5*CF*ROAVG(2)*WBOUND(2)*WBOUND(2)
           VSTAR   = SQRT(TWALLI1/ROAVG(2))
           PLUSK   = ROUGH*VSTAR*ROAVG(2)/VISLAM
           RELOG10 = LOG10(PERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF 
C
      TWALLI1    = 0.5*CF*ROAVG(2)*WBOUND(2)*WBOUND(2)
C
      ELSE
C   IF YPLUSWALL > 5 USE THE YPLUSWALL VALUE TO CALCULATE THE SKIN FRICTION
      TWALLI1    = CFWALL*ROAVG(1)*WABS(1)*WABS(1)
C   END OF YPLUSWALL < 5 LOOP
      ENDIF
C
  365 CONTINUE
C
C******************************************************************************
C******************************************************************************
C  NEXT THE I = IM , UPPER, SURFACE
C
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           PERPIM    =  VOL(IMM1,J,K)/AREA(IM)
           YPLUS_OLD =  YPLUS_IM(J,K)
C
           CALL WALLFUN(IMM1,J,K,IM,PERPIM,DPDS_CELL(IMM1,J,K),
     &             ROAVG(IMM1),TWALLIM,YPLUS_OLD,WBOUND(IMM1),YPLUS_NEW)
           YPLUS_IM(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 375
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
C
      IF(YPLUSWALL.LT.5.) THEN
C
      PERP   =  VOL(IMM1,J,K)/AREA(IM)
      RE     =  ROAVG(IMM1)*WABS(IMM1)*PERP/VISLAM
      RELOG  =  1.0/ALOG(RE)
C
C    ALLOW FOR ROUGHNESS
C
       ROUGH = ROUGH_U(NRW)
       IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.PERP) ROUGH = PERP
           REK   = RE*ROUGH/PERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C    ALLOW FOR ROUGHNESS IF  ROUGH > 1.0E-7 .
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLIM = 0.5*CF*ROAVG(IMM1)*WBOUND(IMM1)*WBOUND(IMM1)
           VSTAR   = SQRT(TWALLIM/ROAVG(IMM1))
           PLUSK   = ROUGH*VSTAR*ROAVG(IMM1)/VISLAM
           RELOG10 = LOG10(PERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
           TWALLIM   = 0.5*CF*ROAVG(IMM1)*WBOUND(IMM1)*WBOUND(IMM1)
C
      ELSE
C   IF YPLUSWALL > 5 USE THE YPLUSWALL VALUE TO CALCULATE THE SKIN FRICTION
           TWALLIM   = CFWALL*ROAVG(IM)*WABS(IM)*WABS(IM)
C    END OF YPLUSWALL < 5 LOOP
      ENDIF
C
  375 CONTINUE
C
C******************************************************************************
C     JDD ADDITION JAN/14  . EVALUATE YPLUS ON THE BLADE SURFACES.
C
      VSTAR     = SQRT(TWALLI1/ROAVG(1))
      DO I=1,IMID
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*ROAVG(1)/VISLAM
      YPLSP     = Y_PLUS(I,J,K)
      Y_PLUS(I,J,K) = AMIN1(YPLS,YPLSP)
      END DO
C
C
      VSTAR     = SQRT(TWALLIM/ROAVG(IM))
      DO I = IMID+1,IMM1
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*ROAVG(IM)/VISLAM
      YPLSP     = Y_PLUS(I,J,K)
      Y_PLUS(I,J,K) = AMIN1(YPLS,YPLSP)
      END DO
C   END OF JAN/14 ADDITION
C
      FMULT     = TWALLI1*AREA(1)/WABS(1)
      XSTRES(1) = FMULT*VXAVG(1)
      RSTRES(1) = FMULT*VRAVG(1)
      TSTRES(1) = FMULT*WTAVG(1)
C
      FMULT     = -TWALLIM*AREA(IM)/WABS(IM)
      XSTRES(IM)= FMULT*VXAVG(IM)
      RSTRES(IM)= FMULT*VRAVG(IM)
      TSTRES(IM)= FMULT*WTAVG(IM)
      QFLOW(1)  = 0.0
      WVISC(1)  = TSTRES(1)*WREL*RAVG(K)
      QFLOW(IM) = 0.0
      WVISC(IM) = TSTRES(IM)*WREL*RAVG(K)
C
   51 CONTINUE
C*********************************************************************************
C*********************************************************************************
C     NOW WORK OUT THE STRESSES ON THE BLADEWISE SURFACES OF THE MAIN FLOW.
C     WITHIN THE BLADE PASSAGE
C
C   JDD ADDED AND CHANGED THE NEXT SECTION  JAN/14.
C    reduce the turbulent viscosity in the sublayer and buffer region.
      DO 55 I=2,IMM1     
      IF(Y_PLUS(I,J,K).LE.YPLAM) FYPLUS = 0.0
      IF(Y_PLUS(I,J,K).GT.YPLAM)  THEN
               XFAC = (Y_PLUS(I,J,K) - YPLAM)/(YPTURB- YPLAM)
               IF(XFAC.GT.1.0) XFAC = 1.0
               FYPLUS = XFAC*XFAC*(3.0 - 2.0*XFAC)
      END IF

      DPERP   = VOL(I,J,K)/AREA(I)
      VLAM(I) = FILAM(NRW,I)*VISLAM/DPERP
      VTURB(I)= FITURB(NRW,I)*ROAVG(I)*ABS(WABS(I+1)-WABS(I-1))*FMIXUP
      VTURB(I)= VTURB(I)*FYPLUS
   55 CONTINUE
C   END OF JDD ADDITION JAN/14
C*********************************************************************************
C*********************************************************************************
C     CHECK FOR TRANSITION
C
      IF(JROW.LT.JTRLOW) GO TO 70
      RMAX=0.0
      IMAX=2
      DO 56 I=2,IMID
      RATVIS = VTURB(I)/VLAM(I)
      IF(RATVIS.LT.RMAX) GO TO 56
      RMAX = RATVIS
      IMAX = I
   56 CONTINUE
      IF(RMAX.GT.FTRANS)  GO TO 71
   70 CONTINUE
      DO 72 I = 2,IMID
   72 VTURB(I)= 0.0
   71 CONTINUE
C
      IF(JROW.LT.JTRUP) GO TO 69
      RMAX = 0.0
      IMAX = IMM1
      DO 73 I = IMID,IMM1
      RATVIS  = VTURB(I)/VLAM(I)
      IF(RATVIS.LT.RMAX) GO TO 73
      RMAX = RATVIS
      IMAX = I
   73 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO  77
   69 CONTINUE
      DO 74 I  = IMID,IMM1
   74 VTURB(I) = 0.0
   77 CONTINUE
C
C*********************************************************************************
C     FORM THE VISCOUS STRESSES AWAY FROM THE BLADE SURFACES
      DO 76 I=2,IMM1
      FMULT        = (VLAM(I)+VTURB(I))*AREA(I)
      XSTRES(I)    = FMULT*(VXAVG(I+1)-VXAVG(I-1))
      RSTRES(I)    = FMULT*(VRAVG(I+1)-VRAVG(I-1))
      TSTRES(I)    = FMULT*(WTAVG(I+1)-WTAVG(I-1))
      VISC_RAT(I,J,K) = VTURB(I)/VLAM(I)
      QFLOW(I)     =  FMULT*FTCOND*(TAVG(I+1)-TAVG(I-1))
      WVISC(I)     =  XSTRES(I)*VXAVG(I) + RSTRES(I)*VRAVG(I) 
     &             +  TSTRES(I)*VTAVG(I)
   76 CONTINUE
C*********************************************************************************
C     SET THE STRESSES ON THE PERIODIC BOUNDARIES ABOVE THE TIP
C  Q3D
      IF(KM.EQ.2) GO TO 78
C  END Q3D
C
      IF((KTIPS(NRW).GT.0).AND.(K.GE.KTIPS(NRW)).AND.(K.LT.KTIPE(NRW))) 
     &    THEN
C
      XSTRES(1)  = 0.5*(XSTRES(2)+XSTRES(IMM1))
      RSTRES(1)  = 0.5*(RSTRES(2)+RSTRES(IMM1))
      TSTRES(1)  = 0.5*(TSTRES(2)+TSTRES(IMM1))
      XSTRES(IM) = XSTRES(1)
      TSTRES(IM) = TSTRES(1)
      RSTRES(IM) = RSTRES(1)
      QFLOW(1)   = 0.5*(QFLOW(2) + QFLOW(IMM1))
      WVISC(1)   = 0.5*(WVISC(2) + WVISC(IMM1))
      QFLOW(IM)  = QFLOW(1)
      WVISC(IM)  = WVISC(1)
C
      ENDIF
C
   78 CONTINUE
C
C     END OF THE CALCULATION OF VISCOUS FORCES WITHIN THE BLADE ROW
C*********************************************************************************
C*********************************************************************************
C
      GO TO 75
C
C      FORM THE VISCOUS STRESSES ON THE BLADEWISE SURFACES
C      OUTSIDE A BLADE ROW WHERE THE MIXING LENGTH IS TAKEN TO BE
C      FRACPW or FRACPUP * THE BLADE PITCH
C
   65 CONTINUE
C
C*********************************************************************************
C*********************************************************************************
C     FORM THE VISCOUS STRESSES IN THE WAKE. BLEND THEM TO THE BLADE VALUES USING FACJ
C
      IF(J.GT.JTEDGE) THEN
C
      FACJ = 0.1*(J-JTEDGE)
      IF(FACJ.GT.1.)  FACJ = 1.
C
      DO 66 I=1,IM
      ISUB = I
      IF(I.EQ.IM) ISUB=IMM1
      IM1 = I-1
      IP1 = I+1
      IF(I.EQ.1)  IM1 = IMM1
      IF(I.EQ.IM) IP1 = 2
      FBLEND       = (FACJ*FIWAKE(NRW,I)+(1.-FACJ)*FITURB(NRW,I))*FMIXUP
      DPERP        = VOL(ISUB,J,K)/AREA(I)
      VLAM(I)      = FILAM(NRW,I)*VISLAM/DPERP
      VTURB(I)     = FBLEND*ROAVG(I)*ABS(WABS(IP1)-WABS(IM1))
      VISC_RAT(I,J,K) = VTURB(I)/VLAM(I)
      FMULT        = AREA(I)*(VLAM(I)+VTURB(I))
      XSTRES(I)    = FMULT*(VXAVG(IP1)-VXAVG(IM1))
      RSTRES(I)    = FMULT*(VRAVG(IP1)-VRAVG(IM1))
      TSTRES(I)    = FMULT*(WTAVG(IP1)-WTAVG(IM1))
      QFLOW(I)     = FMULT*FTCOND*(TAVG(IP1) - TAVG(IM1))
      WVISC(I)     = XSTRES(I)*VXAVG(I) + TSTRES(I)*VTAVG(I)
     &             + RSTRES(I)*VRAVG(I)
   66 CONTINUE
C
      ENDIF
C*********************************************************************************
C     FORM THE VISCOUS STRESSES UPSTREAM OF THE LEADING EDGE.
C
      IF(J.LE.JLEDGE) THEN
C
      DO 67 I=1,IM
      ISUB = I
      IF(I.EQ.IM) ISUB=IMM1
      IM1 = I-1
      IP1 = I+1
      IF(I.EQ.1)  IM1 = IMM1
      IF(I.EQ.IM) IP1 = 2
      DPERP       = VOL(ISUB,J,K)/AREA(I)
      VLAM(I)     = FILAM(NRW,I)*VISLAM/DPERP
      VTURB(I)    = FMIXUP*FIUP(NRW,I)*ROAVG(I)*ABS(WABS(IP1)-WABS(IM1))
      VISC_RAT(I,J,K)= VTURB(I)/VLAM(I)
      FMULT       = AREA(I)*(VLAM(I)+VTURB(I))
      XSTRES(I)   = FMULT*(VXAVG(IP1)-VXAVG(IM1))
      RSTRES(I)   = FMULT*(VRAVG(IP1)-VRAVG(IM1))
      TSTRES(I)   = FMULT*(WTAVG(IP1)-WTAVG(IM1))
      QFLOW(I)    = FMULT*FTCOND*(TAVG(IP1) - TAVG(IM1))
      WVISC(I)    = XSTRES(I)*VXAVG(I) + TSTRES(I)*VTAVG(I)
     &            + RSTRES(I)*VRAVG(I)
   67 CONTINUE
C
      ENDIF
C
C   MAKE THE VISCOUS STRESSES PERIODIC UPSTREAM AND DOWNSTREAM OF THE BLADE.
C

      XSTRES(1)  = 0.5*(XSTRES(IM)+XSTRES(1))
      RSTRES(1)  = 0.5*(RSTRES(IM)+RSTRES(1))
      TSTRES(1)  = 0.5*(TSTRES(IM)+TSTRES(1))
      XSTRES(IM) = XSTRES(1)
      RSTRES(IM) = RSTRES(1)
      TSTRES(IM) = TSTRES(1)
      QFLOW(1)   = 0.5*(QFLOW(IM) + QFLOW(1))
      WVISC(1)   = 0.5*(WVISC(IM) + WVISC(1))
      QFLOW(IM)  = QFLOW(1)
      WVISC(IM)  = WVISC(1)
C
C*********************************************************************************
C*********************************************************************************
C   END OF CALCULATION OF VISCOUS FORCES ON THE BLADEWISE (I) SURFACES
C
   75 CONTINUE
C
C      COMPLETE THE VISCOUS FORCE TERMS BY ADDING THE DIFFERENCE OF THE
C      STRESSES ON THE BLADEWISE FACES OF THE ELEMENTS. RELAXING THE CHANGES BY RFVIS .
C
C   Q3D
      IF(KM.NE.2) THEN
C   END Q3D
      DO 80 I=1,IMM1
      XFORCE(I,J,K)  = XFORCE(I,J,K)+RF_VIS*(XSTRES(I)-XSTRES(I+1))
      RFORCE(I,J,K)  = RFORCE(I,J,K)+RF_VIS*(RSTRES(I)-RSTRES(I+1))
      TFORCE(I,J,K)=TFORCE(I,J,K)+RF_VIS*(TSTRES(I)-TSTRES(I+1))*RAVG(K)
      QSOURCE(I,J,K) = QSOURCE(I,J,K) + RF_VIS*(QFLOW(I) - QFLOW(I+1)
     &               + WVISC(I) - WVISC(I+1))    
   80 CONTINUE
C
      ELSE
C   Q3D  .  IF KM = 2
      DO 81 I=1,IMM1
      XFORCE(I,J,K)  = RF_VIS1*XFORCE(I,J,K)
     &               + RF_VIS*(XSTRES(I)-XSTRES(I+1))
      RFORCE(I,J,K)  = RF_VIS1*RFORCE(I,J,K)
     &               + RF_VIS*(RSTRES(I)-RSTRES(I+1))
      TFORCE(I,J,K)  = RF_VIS1*TFORCE(I,J,K)
     &               + RF_VIS*(TSTRES(I)-TSTRES(I+1))*RAVG(K)
      QSOURCE(I,J,K) = RF_VIS1*QSOURCE(I,J,K)
     &   + RF_VIS*(QFLOW(I) - QFLOW(I+1) + WVISC(I) - WVISC(I+1)) 
   81 CONTINUE
C  END Q3D
C
      END IF
C
C
C***************************************************************************
C     COMPLETE THE ENTROPY GENERATION RATE PER UNIT VOLUME IF IFEND = 1 OR IFPRINT = 1.
C     BY ADDING THE GENERATION DUE TO THE BLADEWISE SURFACES.
C     NON DIMENSIONALISE THE ENTROPY GENERATION RATE BY:
C     0.001*(ONE BLADE ROW PRESSURE CHANGE) *(BLADE EXIT VELOCITY)/(BLADE PITCH).
C
      IF(IFPRINT.EQ.1.OR.IFEND.EQ.1) THEN
      DELP_BLADE = ABS((PO1(KMID) - 0.5*(PDOWN_HUB + PDOWN_TIP)))/NROWS
      TAU_REF    = 0.001*DELP_BLADE
      RHO_REF    =  PO1(KMID)/RGAS/TO1(KMID)
      VEL_REF    = SQRT(2*DELP_BLADE/RHO_REF)
      PITCH      = 2*3.14159*R(1,KMID)/NBLADE(1)
      SREF       =  TAU_REF*VEL_REF/PITCH
      DO 91 I = 1,IMM1
      XGEN  = (XSTRES(I)+XSTRES(I+1))*(VXAVG(I+1) -VXAVG(I))
      RGEN  = (RSTRES(I)+RSTRES(I+1))*(VRAVG(I+1) -VRAVG(I))
      TGEN  = (TSTRES(I)+TSTRES(I+1))*(WTAVG(I+1) -WTAVG(I))
      TAVGG = (TAVG(I) + TAVG(I+1))
      QGEN  = 2.0*(QFLOW(I) + QFLOW(I+1))*(TAVG(I+1) - TAVG(I))/TAVGG
      SGEN_NOW     =  (XGEN + RGEN + TGEN + QGEN)/VOL(I,J,K)/TAVGG
      TEMP4(I,J,K) =  (SGEN(I,J,K) + SGEN_NOW)/SREF
   91 CONTINUE
C
      END IF
C
C     END OF ENTROPY GENERATION RATE CALCULATIION ON BLADEWISE SURFACES
C******************************************************************************8
C
C     END OF K LOOP FOR VISCOUS CALCULATIONS ON THE STREAMWISE SURFACES.
   90 CONTINUE
C
C   END OF J LOOP
  100 CONTINUE
C*********************************************************************************
C*********************************************************************************
C    DISTRIBUTE ENTROPY GENERATION RATE TO THE CELL CORNERS.
C
      IF(IFPRINT.EQ.1.OR.IFEND.EQ.1) THEN
C
           DO K=1,KM-1
           DO I=1,IM-1
           TEMP4(I,1,K) = 0.0
           END DO
           END DO
C
           CALL CELL_TO_NODE(TEMP4,SGEN)
C
      END IF
C
C*********************************************************************************
C    END OF SUBROUTINE LOSS 
C*********************************************************************************
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C*****************************************************************************
C******************************************************************************
C
      SUBROUTINE INPINT
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION  SUMFN(MAXKI),SUMF_IN(MAXKI),ANS(MAXKI)
C
C     THIS SUBROUTINE INTERPOLATES IN THE INLET FLOW DATA AT KIN POINTS
C     SPACED BY FR_IN(K) TO PRODUCE NEW DATA AT KM POINTS SPACED VIA FR(K).
C
      WRITE(6,*) '  ENTERED  INPINT '
C
C    MAKE  FR(K)  and FR_IN(K)  BOTH SUM TO 1.0 .
C
      SUMFN(1) = 0.0
      DO 20 K=2,KM
           SUMFN(K) = SUMFN(K-1) + FR(K-1)
   20 CONTINUE
      DO 30 K=2,KM
           SUMFN(K) = SUMFN(K)/SUMFN(KM)
   30 CONTINUE
C
      SUMF_IN(1)=0.0
      DO 40 K=2,KIN
           SUMF_IN(K) = SUMF_IN(K-1) + FR_IN(K-1)
   40 CONTINUE
      DO 50 K=2,KIN
           SUMF_IN(K) = SUMF_IN(K)/SUMF_IN(KIN)
   50 CONTINUE
C
C
      WRITE(6,*) ' SUM FR FOR ACTUAL GRID POINTS '
      WRITE(6,10)(SUMFN(K),K=1,KM)
      WRITE(6,*)
      WRITE(6,*) ' SUM FR  FOR BOUNDARY CONDITION POINTS'
      WRITE(6,10)(SUMF_IN(K),K=1,KIN)
      WRITE(6,*)
   10 FORMAT(8F10.5)
C
C    INTERPOLATE FOR THE STAGNATION PRESSURES AT THE GRID POINTS.
      DO 100 K=1,KM
           CALL INTP(KIN,SUMF_IN,PO1,SUMFN(K),ANS(K))
  100 CONTINUE
      DO 110 K=1,KM
           PO1(K)=ANS(K)
  110 CONTINUE
C
      WRITE(6,*) 'PO1 AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,11)(PO1(K),K=1,KM)
   11 FORMAT(8F10.1)
C
C    INTERPOLATE FOR THE STAGNATION TEMPERATURES AT THE GRID POINTS.
      DO 120 K=1,KM
           CALL INTP(KIN,SUMF_IN,TO1,SUMFN(K),ANS(K))
  120 CONTINUE
      DO 130 K=1,KM
           TO1(K)=ANS(K)
  130 CONTINUE
C
      WRITE(6,*)  'TO1 AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,10)(TO1(K),K=1,KM)
C
C    INTERPOLATE FOR THE TANGENTIAL VELOCITY AT THE GRID POINTS.
      DO 140 K=1,KM
           CALL INTP(KIN,SUMF_IN,VTIN,SUMFN(K),ANS(K))
  140 CONTINUE
      DO 150 K=1,KM
           VTIN(K)=ANS(K)
  150 CONTINUE
C
      WRITE(6,*)  'VTIN AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,10)(VTIN(K),K=1,KM)
C
C    INTERPOLATE FOR THE MERIDIONAL VELOCITY AT THE GRID POINTS.
      DO 160 K=1,KM
           CALL INTP(KIN,SUMF_IN,VM1,SUMFN(K),ANS(K))
  160 CONTINUE
      DO 170 K=1,KM
           VM1(K)=ANS(K)
  170 CONTINUE
C
      WRITE(6,*)  'VM AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,10)(VM1(K),K=1,KM)
C
      IF(IPOUT.EQ.3) THEN
C
C    INTERPOLATE FOR THE STATIC PRESSURES AT THE GRID POINTS.
      DO 180 K=1,KM
           CALL INTP(KIN,SUMF_IN,PD,SUMFN(K),ANS(K))
  180 CONTINUE
      DO 190 K=1,KM
           PD(K)=ANS(K)
  190 CONTINUE
C
      WRITE(6,*)  'P STATIC AT GRID POINTS ON EXIT BOUNDARY.'
      WRITE(6,11)(PD(K),K=1,KM)
C
      END IF
C
C    INTERPOLATE FOR THE MERIDIONAL PITCH ANGLE AT THE GRID POINTS.
      DO 200 K=1,KM
           CALL INTP(KIN,SUMF_IN,BR,SUMFN(K),ANS(K))
  200 CONTINUE
      DO 210 K=1,KM
           BR(K)=ANS(K)
  210 CONTINUE
C
      WRITE(6,*) 'PITCH ANGLE AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,10)(BR(K),K=1,KM)
C
C    INTERPOLATE FOR THE YAW ANGLE AT THE GRID POINTS.
      DO 220 K=1,KM
           CALL INTP(KIN,SUMF_IN,BS,SUMFN(K),ANS(K))
  220 CONTINUE
      DO 230 K=1,KM
           BS(K)=ANS(K)
  230 CONTINUE
C
      WRITE(6,*)  'YAW ANGLE AT GRID POINTS ON INLET BOUNDARY.'
      WRITE(6,10)(BS(K),K=1,KM)
C
C
      RETURN
      END
C
C******************************************************************************C
C******************************************************************************C
C******************************************************************************C
C
      SUBROUTINE GRID_DOWN(K,J1,J2,S1,S2,NGAP,SDIST,XINT,RINT,NEXTRAP,
     &  BETA_DWN1,BETA_DWN2,IFCUSP,ICUSP,LCUSP,LCUSPUP,IFANGLES)
C
C     THIS SUBROUTINE SETS THE GRID DOWNSTREAM OF THE TRAILING EDGE
C     DOWN TO THE MIXING PLANE OR DOWNSTREAM BOUNDARY.
C     IT ALSO FITS A CUSP AT THE TRAILING EDGE AND SETS THE GRID ANGLES
C     DOWNSTREAM OF THE TRAILING EDGE
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION  SDIST(JD),XINT(JD),RINT(JD),SNEW(JD),ANSX(JD),
     &           ANSR(JD),TH_UPP(JD),TH_THICK(JD),TH_MID(JD)
C
C******************************************************************************C
C    CALCULATE THE MERIDIONAL GRID SPACINGS DOWNSTREAM OF THE BLADE. USING
C    A GEOMETRIC PROGRESSION WITH RATIO  "RAT"  WHICH IS SET BY "SOLVE".
C******************************************************************************C
C     SET THE NEW MERIDIONAL GRID USING THE GEOMETRIC PROGRESSION FOUND BY  "SOLVE" .
C
      SNEW(J1) = S1
      SNEW(J2) = S2
      SDIFF    = S2-S1
      NINT = J2-J1
      R1   = SMERID(J1,K)-SMERID(J1-1,K)
      CALL SOLVE(NINT,SDIFF,R1,RAT,K)
      XD=R1
      DO 10 J = J1+1,J2-1
      SNEW(J) = SNEW(J-1) + XD
      XD      = XD*RAT
   10 CONTINUE
C
C******************************************************************************C
C     INTERPOLATE TO FIND XSURF  AND RSURF AT THE NEW GRID POINTS
C
      DO 30 J = J1,J2
      CALL INTP(NGAP,SDIST,XINT,SNEW(J),ANSX(J))
      CALL INTP(NGAP,SDIST,RINT,SNEW(J),ANSR(J))
   30 CONTINUE
      DO 31 J = J1,J2
      XSURF(J,K)   = ANSX(J)
      RSURF(J,K)   = ANSR(J)
   31 CONTINUE
C
C******************************************************************************C
C******************************************************************************C
C     Maintain the cusp (if any) set in the input data if IFCUSP = 0.
C
      IF(IFCUSP.EQ.0) THEN
C

      IF(IFANGLES.EQ.0) THEN
C     RESET THE DOWNSTREAM GRID ANGLES IF THE INPUT VALUE OF IF_ANGLES  IS ZERO.
C
      JSTRT      = J1 - NEXTRAP
      DO 20    J = JSTRT,J2
      XD          = XSURF(J,K) - XSURF(J-1,K)
      RD          = RSURF(J,K) - RSURF(J-1,K)
      SNEW(J)     = SNEW(J-1) + SQRT(XD*XD + RD*RD)
      TH_UPP(J)   = RT_UPP(J,K)/RSURF(J,K) 
      TH_THICK(J) = RT_THICK(J,K)/RSURF(J,K)
      TH_MID(J)   = TH_UPP(J) - 0.5*TH_THICK(J)
   20 CONTINUE
C
      DSMER      = SNEW(J1)  - SNEW(JSTRT)
      SLOPES     = (TH_UPP(J1)   - TH_UPP(JSTRT))/DSMER
      SLOPETK    = (TH_THICK(J1) - TH_THICK(JSTRT))/DSMER
      SLOPEP     = SLOPES - SLOPETK
      SLOPECENT  = 0.5*(SLOPES+SLOPEP)
C   SET BETADWN1  FROM THE SLOPE OF THE BLADE CENTER LINE AT THE TRAILING EDGE.
C   IF IF_ANGLES  WAS ZERO.
      BETADOWN1  = SLOPECENT
      BETADOWN2  = BETADOWN1
C
      ELSE
C
C   USE THE INPUT VALUES OF BETADWN1 AND BETADWN2   IF IF_ANGLES  WAS NOT ZERO.
      BETADOWN1 = TAN(BETA_DWN1*DEGRAD)/RSURF(J1,K)
      BETADOWN2 = TAN(BETA_DWN2*DEGRAD)/RSURF(J2,K)
C
      END IF
C
C******************************************************************************C
C     FIND THE END OF ANY EXISTING CUSP, WHERE THE THICKMESS = 0 .
      JS = J1
      DO 35 J=J1,J2
           IF(RT_THICK(J,K).LT.1.0E-4*SDIFF) THEN
                JS = J
                GO TO 36 
           END IF
   35 CONTINUE
C
   36 CONTINUE
C
C     JS  IS THE END POINT OF ANY EXISTING CUSP.
C     DO NOT CHANGE THE EXISTING CUSP UPSTREAM OF JS .
C     CHANGE FROM R*THETA TO THETA TO GET BETTER EXTRAPOLATION WHEN THE RADIUS CHANGES.
C     SET THE GRID ANGLE, R_THETA, AND THICKNESS DOWNSTREAM OF THE EXISTING CUSP.
C
      DO 40 J = JS+1,J2
      FDOWN     = (SNEW(J) - SNEW(J1))/(SNEW(J2) - SNEW(J1))
      BETADOWN  = BETADOWN1 + FDOWN*(BETADOWN2 - BETADOWN1)
      DSDIST    = SNEW(J) - SNEW(J-1)
      RT_UPP(J,K)   = RT_UPP(J-1,K) + RSURF(J,K)*BETADOWN*DSDIST
      RT_THICK(J,K) = 0.0
   40 CONTINUE
C
      RETURN
C
C   END OF IFCUSP = 0 OPTION. 
      END IF
C
C*******************************************************************************
C*******************************************************************************
C     NOW GENERATE A NEW CUSP  IF IFCUSP IS NOT ZERO
C
C     CHANGE FROM R*THETA TO THETA TO GET BETTER EXTRAPOLATION WHEN THE RADIUS CHANGES.
C     RESET SNEW  FOR THOSE POINTS ON THE BLADE AND DOWNSTREAM WHICH MIGHT BE USED LATER.
      LUP   = 20
      JSTRT = J1 - LUP
      SNEW(JSTRT-1) = 0.0
      DO 100    J = JSTRT,J2
      XD          = XSURF(J,K) - XSURF(J-1,K)
      RD          = RSURF(J,K) - RSURF(J-1,K)
      SNEW(J)     = SNEW(J-1) + SQRT(XD*XD + RD*RD)
      TH_UPP(J)   = RT_UPP(J,K)/RSURF(J,K) 
      TH_THICK(J) = RT_THICK(J,K)/RSURF(J,K)
      TH_MID(J)   = TH_UPP(J) - 0.5*TH_THICK(J)
  100 CONTINUE
C
      JS = JSTRT + 1
      IF(JS.GT.J1)  JS = J1
C
C  SEARCH FOR THE START OF THE TRAILING EDGE THINNING.
      DO 102 J = JS,J1
           RTHIK   = TH_THICK(J)/TH_THICK(J-1)
           IF(RTHIK.LT.0.9) GO TO 103
  102 CONTINUE
C
  103 CONTINUE
C
C   JCS  IS THE START OF THE TRAILING EDGE THINNING
      JCS = J-1 
C
C     EXTAPOLATE THE BLADE SURFACES FROM UPSTREAM, J = JCS, TO THE TRAILING EDGE 
C     SO THAT THE TE THINNING IS REMOVED.
      TKGRAD = (TH_THICK(JCS) - TH_THICK(JCS-2))/(SNEW(JCS)-SNEW(JCS-2))
      YGRAD  = (TH_MID(JCS)   - TH_MID(JCS-2))/(SNEW(JCS) - SNEW(JCS-2))     
      DO 104 J = JCS,J1
           TH_MID(J)   = TH_MID(JCS)   +  YGRAD*(SNEW(J) - SNEW(JCS))
           TH_THICK(J) = TH_THICK(JCS) + TKGRAD*(SNEW(J) - SNEW(JCS))
           IF(TH_THICK(J).LT.0.0) TH_THICK(J) = 0.0
           TH_UPP(J)     = TH_MID(J)  + 0.5*TH_THICK(J)
           RT_UPP(J,K)   = RSURF(J,K)*TH_UPP(J)
           RT_THICK(J,K) = RSURF(J,K)*TH_THICK(J)
  104 CONTINUE
C
C   END OF ADJUSTING THE BLADE UPSTREAM OF THE TRAILING EDGE
C
C**********************************************************************************
C     SET THE SLOPE OF THE BLADE SURFACE AT THE NEW TRAILING EDGE.
C     WHICH IS AT  JSTRT = J1 - LCUSPUP
C
      JSTRT      = J1 - LCUSPUP
      DSMER      = SNEW(JSTRT)  - SNEW(JSTRT-2)
      SLOPES     = (TH_UPP(JSTRT)      - TH_UPP(JSTRT-2))/DSMER
      SLOPETK    = (TH_THICK(JSTRT)    - TH_THICK(JSTRT-2))/DSMER
      SLOPEP     = SLOPES - SLOPETK
      SLOPECENT  = 0.5*(SLOPES+SLOPEP)
C 
      IF(IFANGLES.EQ.0) THEN
C     USE THE EXTRAPOLATED BLADE CENTRE LINE ANGLE IF  IF_ANGLES = ZERO.
	   BETADOWN1 = SLOPECENT
           BETADOWN2 = SLOPECENT
      ELSE
C     USE THE INPUT DOWNSTREAM GRID ANGLES IF IF_ANGLES > ZERO .
	   BETADOWN1 = TAN(BETA_DWN1*DEGRAD)/RSURF(JSTRT,K)
	   BETADOWN2 = TAN(BETA_DWN2*DEGRAD)/RSURF(JSTRT,K)
      ENDIF
C
C**********************************************************************************
C     CALCULATE THE TRAILING EDGE THICKNESS AND CENTRELINE RTHETA VALUE AT THE.
C     NEW TRAILING EDGE, WHICH IS AT  J = J1 - LCUSPUP
C
      YTHICKTE = RT_THICK(JSTRT,K)
      IF(YTHICKTE.LT.0.0)  YTHICKTE = 0.0
      YMIDTE   = RT_UPP(JSTRT,K) - 0.5*YTHICKTE
C
C     Form the cusp of length  LCUSP, Starting LCUSPUP points before the.
C     trailing edge.
C     The cusp is centred on the blade centre line if ICUSP = 0.
C     The cusp makes the I=1 surface continuous on the cusp if  ICUSP =  1.
C     The cusp makes the I=IM surface continuous on the cusp if ICUSP = -1.
C     The cusp extends LCUSPUP upstream onto the solid part of the blade.
C
      DO 500 J = JSTRT+1,J2
C
      FDOWN     = (SNEW(J) - SNEW(JSTRT))/(SNEW(J2) - SNEW(JSTRT))
      BETADOWN  = BETADOWN1 + FDOWN*(BETADOWN2 - BETADOWN1)
      DSDIST    = SNEW(J) - SNEW(J-1)
      YMID      = YMIDTE  + RSURF(J,K)*SLOPECENT*(SNEW(J) - SNEW(JSTRT))
      CLENGTH   = SNEW(JSTRT+LCUSP) - SNEW(JSTRT)
      PLENGTH   = SNEW(J) - SNEW(JSTRT)
      IF(PLENGTH.GT.CLENGTH) PLENGTH = CLENGTH
C
C     SET THE DOWNSTREAM GRID SLOPE ACCORDING TO  "ICUSP" AND "LCUSP".
C     NOTE THAT THE LOCAL RADIUS, RSURF, IS USED TO CHANGE BACK FROM THETA TO RTHETA.
C
      JCUSPEND  = JSTRT+LCUSP
      IF((JCUSPEND-J).GE.0)   THEN
C     IF ON THE CUSP
	   YTNEW         = YTHICKTE*(1. -PLENGTH/CLENGTH)
	   RT_THICK(J,K) = YTNEW
	   IF(ICUSP.EQ.0) RT_UPP(J,K) = YMID + 0.5*YTNEW
	   IF(ICUSP.EQ.1) RT_UPP(J,K) = RT_UPP(J-1,K)
     &                                + RSURF(J,K)*SLOPES*DSDIST
	   IF(ICUSP.EQ.-1) RT_UPP(J,K) = RT_UPP(J-1,K) - RT_THICK(J-1,K) 
     &                   + RSURF(J,K)*SLOPEP*DSDIST  + YTNEW
      ELSE
C     IF PAST END OF CUSP
	   RT_UPP(J,K)   = RT_UPP(J-1,K) + RSURF(J,K)*BETADOWN*DSDIST
	   RT_THICK(J,K) = 0.0
      ENDIF
C
  500 CONTINUE
C
C     END OF NEW CUSP GENERATION. NOVEMBER 2015 .  
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE GRID_UP(K,J1,J2,S1,S2,NGAP,SDIST,XINT,RINT,NEXTRAP,
     &                 BETA_UP,IFANGLES)
C
      INCLUDE  'commall-open-18.3'
C
C     THIS SUBROUTINE SETS THE GRID UPSTREAM OF THE LEADING EDGE,
C     UP TO THE MIXING PLANE OR UPSTREAM BOUNDARY.
C
      DIMENSION SDIST(JD),XINT(JD),RINT(JD),SNEW(JD),ANSX(JD),
     &          ANSR(JD),TH_UPP(JD),TH_THICK(JD)
C
C      FORM THE NEW MERIDIONAL GRID SPACINGS "SNEW" USING A GEOMETRIC PROGRESSION
C      WITH EXPANSION RATIO "RAT" WHICH IS CALCULATED BY SUBROUTINE  "SOLVE" .
C
      SNEW(J1) = S1
      SNEW(J2) = S2
      SDIFF    = S2 - S1
      NINT = J2 - J1
      R1 = SMERID(J2+1,K)-SMERID(J2,K)
      CALL SOLVE(NINT,SDIFF,R1,RAT,K)
      XD = R1
      DO 10 J = J1+1,J2-1
      JSUB    = J2 + J1 - J
      SNEW(JSUB) = SNEW(JSUB+1) - XD
      XD = XD*RAT
   10 CONTINUE
C
C     SET THE GRID POINTS AT THE UPSTREAM MIXING PLANE TO BE VERY CLOSELY,
C     BUT NOT EXACTLY, COINCIDENT, EXCEPT FOR THE FIRST ROW.
C
      IF(J1.NE.1) SNEW(J1) = SNEW(J1) + 0.01*(SNEW(J1+1)-SNEW(J1))
C
C******************************************************************************C
C  INTERPOLATE TO FIND  XSURF   AND  RSURF   AT THE NEW GRID POINTS.
C
      DO 30 J = J1,J2
      CALL INTP(NGAP,SDIST,XINT,SNEW(J),ANSX(J))
      CALL INTP(NGAP,SDIST,RINT,SNEW(J),ANSR(J))
   30 CONTINUE
      DO 35 J = J1,J2
      XSURF(J,K)   = ANSX(J)
      RSURF(J,K)   = ANSR(J)
   35 CONTINUE
C
C******************************************************************************C
C******************************************************************************C
C     CHANGE FROM R*THETA TO THETA TO GET A BETTER EXTRAPOLATION FOR RADIAL BLADES.
C
      DO 100   J = J1,J2 + NEXTRAP
      TH_UPP(J)   = RT_UPP(J,K)/RSURF(J,K)
      TH_THICK(J) = RT_THICK(J,K)/RSURF(J,K)
  100 CONTINUE
C
C     SET THE SLOPE OF THE UPSTREAM GRID AS  DTHETA/DS, OR USE THE INPUT VALUE 
C     OF  BETA_UP  IF  IF_ANGLES  IS GREATER THAN ZERO
C
      IF(IFANGLES.EQ.0) THEN
      SLOPET = (TH_UPP(J2) - 0.5*TH_THICK(J2) - TH_UPP(J2+NEXTRAP)
     & + 0.5*TH_THICK(J2+NEXTRAP))/(SMERID(J2,K) - SMERID(J2+NEXTRAP,K))
      ELSE
      SLOPET = TAN(BETA_UP*DEGRAD)/RSURF(J2,K)
      ENDIF
C
C     USE THE CALCULATED SLOPE TO SET  RTHETA  UPSTREAM OF THE LEADING EDGE.
C     MULTIPLY BY THE LOCAL RADIUS, RSURF, TO CHANGE BACK FROM   THETA  TO  RTHETA.
C
      RTHETA_LE = RT_UPP(J2,K) - 0.5*RT_THICK(J2,K)
      DO 20 J = J1,J2-1
      RT_UPP(J,K) = RTHETA_LE + RSURF(J,K)*SLOPET*(SNEW(J)-SNEW(J2))
   20 CONTINUE
C
C******************************************************************************C
C******************************************************************************C
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************C
C
      SUBROUTINE SOLVE(NPUP,XUP,DX1,RAT,K)
C
C     FIND THE RATIO OF SPACINGS FOR THE UPSTREAM AND DOWNSTREAM GRIDS 
C 
C  
      RUP = XUP/DX1 
      RAT = 1.2
      IF(RUP.LT.NPUP) RAT = 0.9
      M = 1 
  50  RHS = (RAT-1.)*RUP 
      IF(RUP.GE.NPUP) RNEW = (RHS+1.)**(1./NPUP)
      IF(RUP.LT.NPUP) RNEW = (RAT**NPUP -1.)/RUP  + 1.
      DIFF =ABS(RNEW/RAT-1.0)
C      MAY BE MORE STABLE IF USE  "RAT = 0.5*(RNEW+RAT) "
      RAT = RNEW
      M = M+1
      IF(DIFF.LT.0.00001) GO TO 51
      IF(M.GT.250) GO TO 51
      GO TO 50 
   51 CONTINUE
C
      IF(M.GE.250) WRITE(6,*) ' WARNING !!!    GRIDUP/DOWN iteration not
     & converged.'
C
      WRITE(6,1) K, M, RAT
    1 FORMAT(' K= ',I5,'GRID_UP/DOWN ITS= ',I5,'GRID EXPANSION RATIO= ',
     &       F10.5)
C
C
      RETURN
      END
C
C******************************************************************************C
C******************************************************************************
C
      SUBROUTINE NEWGRID(JMROW,J1,J2,JLEROW,JTEROW,JROTHS,JROTHE,
     & JROTTS,JROTTE)
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION UPF(JD),ONF(JD),DOWNF(JD),SOLD(JD),SNEW(JD),
     & XBUF(JD),YBUF(JD),YSNEW(JD), XFRACUP(JD),RELSPUP(JD),
     & XFRACON(JD),RELSPON(JD),XFRACDWN(JD),RELSPDWN(JD)
C
C      THIS SUBOUTINE INTERPOLATES IN THE INPUT DATA TO SET UP A NEW GRID SPACING IN THE STREAMWISE (J)
C      DIRECTION.
C
C      NUP  -   IS NUMBER OF UPSTREAM POINTS INCLUDING THE L.E. POINT
C      		I.E. IT IS THE NUMBER OF UPSTREAM ELEMENTS + 1.
C      NON  -   IS NUMBER OF POINTS ON BLADE INCLUDING THE T.E. POINT BUT NOT
C      		THE L.E. POINT. IT = THE NUMBER OF ELEMENTS ON THE BLADE.
C      NDOWN  - IS THE NUMBER OF POINTS DOWNSTREAM NOT INCLUDING THE T.E. POINT
C     		IT = THE NUMBER OF ELEMENTS DOWNSTREAM OF THE BLADE.
C
C
C      THE RELATIVE SPACINGS ARE ALL ORDERED IN THE DIRECTION OF
C      INCREASING   J  . I.E. UPF(1) WILL TEND TO BE LARGE AND UPF(NUP)
C      WILL TEND TO BE SMALL.
C
C       UPF(NUP) IS NOT USED AS THERE ARE NUP-1 INTERVALS BETWEEN NUP POINTS
C
      READ(5,*)     DUMMY_INPUT
      WRITE(6,*)    DUMMY_INPUT
      WRITE(6,*) 'STARTING TO GENERATE A NEW GRID IN SUBROUTINE NEWGRID'
C
      READ(5,*)     NUP,NON,NDOWN
      WRITE(6,1600) NUP,NON,NDOWN
 1600 FORMAT('  NUP, NON, NDOWN = ',3I10)
C
C******************************************************************************
C     INPUT ONLY A FEW RELATIVE GRID SPACINGS
C     AND GENERATE THE GRID SPACINGS AUTOMATICALLY IF  NUP, NON OR NDOWN = 0.
C
C     READ IN  THE RELATIVE GRID SPACINGS, UPF(J) , ONF(J)  and  DOWNFJ) ,
C     if  NUP, NON or NDOWN > 0 .
C
C******************************************************************************
C******************************************************************************
C
      IF(NUP.EQ.0) THEN
           READ(5,*) NUP
           NINUP = 0
 880       CONTINUE
           NINUP = NINUP+1
           READ(5,*) XFRACUP(NINUP),RELSPUP(NINUP)
           IF(XFRACUP(NINUP).GT.0.9999) GO TO 881
           GO TO 880
 881       CONTINUE
           DO 886 N=1,NUP
                XX = FLOAT(N-1)/(NUP-1)
                CALL INTP(NINUP,XFRACUP,RELSPUP,XX,UPF(N))
 886       CONTINUE
C
      ELSE
           READ(5,*)     (UPF(J),J=1,NUP)
      ENDIF
C
C******************************************************************************
      IF(NON.EQ.0) THEN
           READ(5,*) NON
           NINON = 0
 882       CONTINUE
           NINON = NINON+1
           READ(5,*) XFRACON(NINON),RELSPON(NINON)
           IF(XFRACON(NINON).GT.0.9999) GO TO 883
           GO TO 882
 883       CONTINUE
           DO 887 N=1,NON
                XX = FLOAT(N-1)/(NON-1)
           CALL INTP(NINON,XFRACON,RELSPON,XX,ONF(N))
 887  CONTI   NUE
C
      ELSE
           READ(5,*)     (ONF(J),J=1,NON)
      ENDIF
C
C******************************************************************************
      IF(NDOWN.EQ.0) THEN
           READ(5,*) NDOWN
           NINDWN = 0
 884       CONTINUE
           NINDWN = NINDWN+1
           READ(5,*) XFRACDWN(NINDWN),RELSPDWN(NINDWN)
           IF(XFRACDWN(NINDWN).GT.0.9999) GO TO 885
           GO TO 884
 885       CONTINUE
           DO 888 N=1,NDOWN
                XX = FLOAT(N-1)/(NDOWN-1)
                CALL INTP(NINDWN,XFRACDWN,RELSPDWN,XX,DOWNF(N))
 888       CONTINUE
C
      ELSE
           READ(5,*)     (DOWNF(J),J=1,NDOWN)
      ENDIF
C
C******************************************************************************
C******************************************************************************
C     UPF(J), ONF(J) AND DOWNF(J) ARE NOW SET, PRINT THEM OUT .
C
      WRITE(6,1704) (UPF(J),J=1,NUP)
      WRITE(6,1704) (ONF(J),J=1,NON)
      WRITE(6,1704) (DOWNF(J),J=1,NDOWN)
1704  FORMAT(1H ,12F5.2)
C
C******************************************************************************
C      READ THE CHANGE IN THE UPSTREAM AND DOWNSTREAM EXTENT OF THE GRID
C      THESE MULTIPLY THE ORIGINAL UPSTREAM AND DOWNSTREAM
C      LENGTHS OF THE GRID BY UPEXT AND DWNEXT.
C
      READ(5,*)     UPEXT,DWNEXT
      WRITE(6,1704) UPEXT,DWNEXT
C
C     END OF INPUT DATA FOR GENERATING THE NEW GRID
C******************************************************************************
C
      JMNEW  = NUP + NON + NDOWN
      JLENEW = NUP
      JTENEW = NUP + NON
C
      IF(JROTHS.GT.JMROW) JROTHS = JMROW
      IF(JROTHE.GT.JMROW) JROTHE = JMROW
      IF(JROTTS.GT.JMROW) JROTTS = JMROW
      IF(JROTTE.GT.JMROW) JROTTE = JMROW
C
C******************************************************************************
C******************************************************************************
C   START THE 'DO 5000' LOOP OVER THE NUMBER OF INPUT SECTIONS
C
      NLOOP = NSECS_IN
      IF(NSECS_IN.EQ.1) NLOOP = 2
C
      DO 5000 K=1,NLOOP
C
C     CALCULATE S COORDS. AND PUT IN ARRAY SOLD
C
      SOLD(1) = 0
      DO 601 J=2,JMROW
      JS = J -1 + J1
      SOLD(J) = SOLD(J-1) + SQRT((RSURF(JS,K)-RSURF(JS-1,K))**2
     &                         + (XSURF(JS,K)-XSURF(JS-1,K))**2)
601   CONTINUE
C
C      USE SNEW TO STORE S-COORDS OF NEW GRID LINES
C
C     FIRST CALCULATE UPTOT, ONTOT, AND DOWNTOT:
C
      UPTOT = 0
      DO 605 J=1,NUP-1
      UPTOT = UPTOT + UPF(J)
605   CONTINUE
C
      ONTOT = 0
      DO 607 J=1,NON
      ONTOT = ONTOT + ONF(J)
607   CONTINUE
C
      DOWNTOT = 0
      DO 609 J=1,NDOWN
      DOWNTOT = DOWNTOT + DOWNF(J)
609   CONTINUE
C
C     FIND THE STREAMWISE COORDINATE WHERE THE ENDWALL ROTATIONS START AND END.
C
      IF(K.EQ.1) THEN
      RHSTART= SOLD(JROTHS)*1.0001
      RHEND  = SOLD(JROTHE)*0.9999
      ENDIF
      IF(K.EQ.NSECS_IN) THEN
      RTSTART= SOLD(JROTTS)*1.0001
      RTEND  = SOLD(JROTTE)*0.9999
      ENDIF
C
C     FIND THE STREAMWISE COORDINATES OF THE NEW GRID
C
      SNEW(1) = SOLD(JLEROW)*(1.- UPEXT)
      DO 610 J=2,NUP
      SNEW(J) = SNEW(J-1) + UPF(J-1)/UPTOT*(SOLD(JLEROW) - SNEW(1))
610   CONTINUE
C
      DO 620 J=NUP+1,NUP+NON
      SNEW(J) =S NEW(J-1) + ONF(J-NUP)/ONTOT*(SOLD(JTEROW)-SOLD(JLEROW))
620   CONTINUE
C
      DO 625 J=NUP+NON+1,NUP+NON+NDOWN
      SNEW(J) = SNEW(J-1) + DOWNF(J-NUP-NON)/DOWNTOT*(SOLD(JMROW)
     &        - SOLD(JTEROW))*DWNEXT
625   CONTINUE
C
C     FIND THE COORDINATES WHERE ENDWALL ROTATION STARTS AND ENDS ON THE NEW GRID
C
      DO 626 J=1,JMNEW
      IF(K.EQ.1) THEN
      IF(SNEW(J).LT.RHSTART) JROTHS = J
      IF(SNEW(J).LT.RHEND)   JROTHE = J + 1
      ENDIF
      IF(K.EQ.NSECS_IN) THEN
      IF(SNEW(J).LT.RTSTART) JROTTS = J
      IF(SNEW(J).LT.RTEND)   JROTTE = J + 1
      ENDIF
  626 CONTINUE
C
C******************************************************************************
C******************************************************************************
C     NOW FIND X,R,YSUCT..,THICK.. USING 4 POINT INTERPOLATION
C
C     SET UP BUFFERS XBUF AND YBUF TO CONTAIN THE OLD VALUES OF  S  AND X .
C
      DO 632 JJ=1,JMROW
      JS = JJ-1+J1
      XBUF(JJ) = SOLD(JJ)
      YBUF(JJ) = XSURF(JS,K)
C      WRITE(6,*) 'JROW,JTOT,XBUF,YBUF,SNEW ',
C     &            JJ,JS,XBUF(JJ),YBUF(JJ),SNEW(JJ)
632   CONTINUE
C
C     FIND INTERPOLATED X-VALS, AND PUT THEM IN XIN:
C
      DO 634 J=1,JMNEW
      JS = J-1+J1
      CALL INTP(JMROW,XBUF,YBUF,SNEW(J),XSURF(JS,K))
634   CONTINUE
C
C     FIND NEW R:
C
      DO 636 JJ=1,JMROW
      JS = JJ-1+J1
      YBUF(JJ) = RSURF(JS,K)
636   CONTINUE 
C
      DO 638 J=1,JMNEW
      JS = J-1+J1
      CALL INTP(JMROW,XBUF,YBUF,SNEW(J),RSURF(JS,K))
638   CONTINUE
C
C     NOW  THE BLADE UPPER SURFACE  RT_UPP .
C
C     FIRST UPSTREAM OF THE LEADING EDGE
C
      DO 640 JJ=1,JLEROW
      JS = JJ-1+J1
      YBUF(JJ) = RT_UPP(JS,K)
640   CONTINUE
C
      DO 642 J=1,JLENEW
      JS = J-1+J1
      CALL INTP(JLEROW,XBUF,YBUF,SNEW(J),YSNEW(JS))
642   CONTINUE
C
C    NEXT ON THE BLADE ROW
C
      DO 643 J=JLEROW,JTEROW
      JS =  J-1+J1
      JJ =  J+1-JLEROW
      YBUF(JJ) = RT_UPP(JS,K)
      XBUF(JJ) = SOLD(J)
643   CONTINUE
C    
      NPOINTS= JTEROW - JLEROW + 1
      DO 644 J=JLENEW,JTENEW
      JS = J-1+J1
      CALL INTP(NPOINTS,XBUF,YBUF,SNEW(J),YSNEW(JS))
644   CONTINUE
C
C  NEXT DOWNSTREAM OF THE TRAILING EDGE
C
      DO 645 J=JTEROW,JMROW
      JS = J-1+J1
      JJ =  J+1-JTEROW
      YBUF(JJ) = RT_UPP(JS,K)
      XBUF(JJ) = SOLD(J)
645   CONTINUE
C
      NPOINTS= JMROW - JTEROW + 1
      DO 646 J=JTENEW,JMNEW
      JS = J-1+J1
      CALL INTP(NPOINTS,XBUF,YBUF,SNEW(J),YSNEW(JS))
646   CONTINUE
C
      DO 670 J=1,JMNEW
      JS = J-1+J1
  670 RT_UPP(JS,K) = YSNEW(JS)
C
C     NOW THE BLADE THICKNESS  RT_THICK
C     FIRST ON THE BLADE SURFACE
C
      DO 647 J=JLEROW,JTEROW
      JS =  J-1+J1
      JJ =  J+1-JLEROW
      YBUF(JJ) = RT_THICK(JS,K)
      XBUF(JJ) = SOLD(J)
647   CONTINUE
C
      NPOINTS= JTEROW - JLEROW + 1
      DO 648 J=JLENEW,JTENEW
      JS = J-1+J1
      CALL INTP(NPOINTS,XBUF,YBUF,SNEW(J),RT_THICK(JS,K))
648   CONTINUE
C
C    SET THE THICKNESS TO ZERO UPSTREAM OF THE LEADING EDGE
      DO 649 J=1,JLENEW
      JS = J-1+J1
      RT_THICK(JS,K) = 0.0
649   CONTINUE
C
C  SET THE THICKNESS TO ZERO DOWNSTREAM OF THE TRAILING EDGE
C  LEAVING 2 POINTS FOR A CUSP.
      DO 651 J=JTENEW+2,JMNEW
      JS = J-1+J1
      RT_THICK(JS,K)=0.0
651   CONTINUE
C
C     END OF NEW GRID GENERATION ON ONE STREAM SURFACE.
C
      WRITE(6,995) K
  995 FORMAT(//,' NEW INTERPOLATED GRID POINTS FOR K =',I5,//)
      WRITE(6,997)
  997 FORMAT('         J     X            Y           R          TK ')
      DO 999 J=1,JMNEW
      JS = J - 1 + J1
  999 WRITE(6,998) JS,XSURF(JS,K),RT_UPP(JS,K),RSURF(JS,K),
     &             RT_THICK(JS,K)
  998 FORMAT( I10,4F12.5)
C
C******************************************************************************
C******************************************************************************
C   RETURN TO START A NEW STREAM SURFACE
C
 5000 CONTINUE
C
C******************************************************************************
C******************************************************************************
C    RESET THE J INDICES FOR THE NEW GRID
C
      JMROW  =  JMNEW
      JLEROW = JLENEW
      JTEROW = JTENEW
      IF(JROTHS.GT.JMNEW) JROTHS = JMNEW
      IF(JROTHE.GT.JMNEW) JROTHE = JMNEW
      IF(JROTTS.GT.JMNEW) JROTTS = JMNEW
      IF(JROTTE.GT.JMNEW) JROTTE = JMNEW
C
      WRITE(6,*) ' RETURNING FROM SUBROUTINE NEWGRID '
      WRITE(6,*)
      WRITE(6,*) 'NEW VALUES OF JM,JLE,JTE RELATIVE TO THE START OF THIS
     & BLADE ROW. ' 
      WRITE(6,*)  JMROW,JLEROW,JTEROW
      WRITE(6,*)
      WRITE(6,*) 'NEW VALUES OF JROTHS,JROTHE,JROTTS,JROTTE RELATIVE TO 
     &THE START OF THIS BLADE ROW.'
      WRITE(6,*)  JROTHS,JROTHE,JROTTS,JROTTE
C
C******************************************************************************
C******************************************************************************
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE SHROUDFLOW
C
      INCLUDE  'commall-open-18.3'
C
C
C     READ IN DATA FOR THE SHROUD LEAKAGE MODEL IF KTIPS(NR) < 1 .
C
      RFSHRD = 0.1
C
      IF(NSTEP.EQ.1) THEN
C
      DO 50 NR = 1,NROWS
C
      KSHROUD(NR) = 0
C
      IF(KTIPS(NR).GE.0) GO TO 75
C 
      WRITE(6,*)
      WRITE(6,*)'SHROUD LEAKAGE DATA FOR ROW NUMBER', NR
      READ(5,*)  KSHROUD(NR),JLEAKS(NR),JLEAKE(NR),JLKINS(NR),JLKINE(NR)
      WRITE(6,*)'KSHROUD, JLEAKS, JLEAKE, JLEAKINS, JLEAKINE',
     &           KSHROUD(NR),JLEAKS(NR),JLEAKE(NR),JLKINS(NR),JLKINE(NR) 
      READ(5,*)  SEALGAP(NR),NSEAL(NR),CFSHROUD(NR),CFCASING(NR)
      WRITE(6,*)'SEALGAP, NSEAL, CFSHROUD, CFCASING',
     &           SEALGAP(NR),NSEAL(NR),CFSHROUD(NR),CFCASING(NR)
      READ(5,*)  WCASE(NR),PITCHIN(NR)
      WRITE(6,*)'RPM CASING, INLET PITCH ANGLE', WCASE(NR),PITCHIN(NR)
      WRITE(6,*)
C
      SEALGAP(NR) = SEALGAP(NR)/SQRT(FLOAT(NSEAL(NR)))
      IF(ABS(PITCHIN(NR)).LT.10.0) PITCHIN(NR)= 10.0
      PITCHIN(NR) = PITCHIN(NR)*DEGRAD
      WCASE(NR)   = WCASE(NR)*PI/30
C
   75 CONTINUE
C
      DO 99 J=2,JM
   99 SHRDFLOW(J) = 0.0
C
   50 CONTINUE
C
C     END OF SHROUD LEAKAGE INPUT DATA
C
C     SET RFSHRD = 1 ON THE FIRST ITERATION ONLY
      RFSHRD = 1.
C
C    END OF FIRST ITERATION ONLY PART
C
      ENDIF
C
C
      SWORKTOT = 0.0
C
C     START DO 100 LOOP OVER ALL BLADE ROWS
C
      DO 100 NR =1,NROWS
C
      IF(KTIPS(NR).GE.0) GO TO 100
C
      JS = JSTART(NR)
      K  = KSHROUD(NR)
      J1 = JLEAKS(NR) + JS -1
      J2 = JLEAKE(NR) + JS -1
      J3 = JLKINS(NR) + JS -1
      J4 = JLKINE(NR) + JS -1
      JSHROUD = 0.5*(J2 + J3)
      NBLD    = NBLADE(JSHROUD)
C
C     FIND THE AREA AVERAGE PRESSURE DOWNSTREAM OF THE SEAL
C
      POUTAVG = 0.0
      ASTOT   = 0.0
      DO 200 I = 1,IMM1
      DO 200 J = J3+1,J4
      ASABS    = SQRT(ASX(I,J,K)*ASX(I,J,K) + ASR(I,J,K)*ASR(I,J,K))
      ASTOT    = ASTOT + ASABS
      POUTAVG  = POUTAVG +
     &           ASABS*(P(I,J,K)+P(I+1,J,K)+P(I,J-1,K)+P(I+1,J-1,K))
      AS(I,J)  = ASABS
  200 CONTINUE
      POUTAVG  = 0.25*POUTAVG/ASTOT
      ALEAKIN  = ASTOT
C
C     FIND THE AVERAGE QUANTITIES FOR THE FLOW LEAKING OVER THE SHROUD
C     AND CALCULATE POAX THE PO NEGLECTING SWIRL VELOCITY
C
      DO 250 J=J1,J2
      DO 250 I=1,IM
      VMERSQ    = VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K)
      EKE       = VMERSQ + VT(I,J,K)*VT(I,J,K)
C
      IF(IFGAS.EQ.0) THEN
           TSTAT     = (HO(I,J,K) - 0.5*EKE)/CP
           TOREL     = TSTAT + 0.5*VMERSQ/CP
           POAX(I,J) = P(I,J,K)*(TOREL/TSTAT)**RFGA
      ELSE
           HSTAG     = HO(I,J,K)
           HSTAT     = HSTAG - 0.5*EKE
           HOREL     = HSTAT + 0.5*VMERSQ
           TSTAT     = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
           TOREL     = TFROMH(HOREL,TREF,HREF,HT1,HT2,HT3,HT4)
           TRAT      = TOREL/TSTAT
           PRAT      = PRAT_FROM_TRAT(TSTAT,TRAT,ALPHA,BETA1,BETA2)
           POAX(I,J) = P(I,J,K)*PRAT
      END IF
C
      TOAX(I,J) = TOREL
  250 CONTINUE
C
      POAXAVG = 0.0
      HOAVG   = 0.0
      TOAXAVG = 0.0
      VTAVG   = 0.0
      VXAVG   = 0.0
      VRAVG   = 0.0
      ASTOT   = 0.0
      DO 300 I=1,IMM1
      DO 300 J=J1+1,J2
      ASABS  = SQRT(ASX(I,J,K)*ASX(I,J,K) + ASR(I,J,K)*ASR(I,J,K))
      ASTOT  = ASTOT + ASABS
      POAXAVG = POAXAVG +
     &        ASABS*(POAX(I,J)+POAX(I+1,J)+POAX(I,J-1)+POAX(I+1,J-1))
      TOAXAVG = TOAXAVG +
     &        ASABS*(TOAX(I,J)+TOAX(I+1,J)+TOAX(I,J-1)+TOAX(I+1,J-1))
      HOAVG = HOAVG +
     &        ASABS*(HO(I,J,K)+HO(I+1,J,K)+HO(I,J-1,K)+HO(I+1,J-1,K))
      VTAVG = VTAVG +
     &        ASABS*(VT(I,J,K)+VT(I+1,J,K)+VT(I,J-1,K)+VT(I+1,J-1,K))
      VRAVG = VRAVG +
     &        ASABS*(VR(I,J,K)+VR(I+1,J,K)+VR(I,J-1,K)+VR(I+1,J-1,K))
      VXAVG = VXAVG +
     &        ASABS*(VX(I,J,K)+VX(I+1,J,K)+VX(I,J-1,K)+VX(I+1,J-1,K))
      AS(I,J) = ASABS
 300  CONTINUE
C
      ALEAKOUT  = ASTOT
      POAXAVG   = 0.25*POAXAVG/ASTOT
      HOLEAKOUT = 0.25*HOAVG/ASTOT 
      VTLEAKOUT = 0.25*VTAVG/ASTOT
      VXLEAKOUT = 0.25*VXAVG/ASTOT
      VRLEAKOUT = 0.25*VRAVG/ASTOT
      TOAXAVG   = 0.25*TOAXAVG/ASTOT

C
C     CALCULATE THE LEAKAGE MASS FLOW RATE.
C
      PRAT  = POUTAVG/POAXAVG
      IF(IFGAS.EQ.0) THEN
            TJET       = TOAXAVG*PRAT**FGA
            IF(TJET.GT.TOAXAVG)  TJET = 0.99999*TOAXAVG
            VJET       = SQRT(2*CP*(TOAXAVG - TJET))
      ELSE
            TRAT     = TRAT_FROM_PRAT(TOAXAVG,PRAT,FGAGAS,R_ALPHA,
     &                 BALPHA1,BALPHA2)
            TJET     = TRAT*TOAXAVG
            HJET     = HFROMT(TJET,TREF,HREF,CP1,CP2,CP3)
            HOJET    = HFROMT(TOAXAVG,TREF,HREF,CP1,CP2,CP3)
            IF(HJET.GT.HOJET) HJET = 0.99999*HOJET
            VJET     = SQRT(2.0*(HOJET-HJET))
      END IF
C            
      ROJET      = POUTAVG/(RGAS*TJET)
      ROSHROUD   = POUTAVG/(RGAS*TOAXAVG)
      CONTRACT   = 0.6
      ALEAK      = 2*PI*R(JSHROUD,K)*SEALGAP(NR)*CONTRACT
      SHROUDMASS = ALEAK*VJET*ROJET
      SHROUDFLBL = SHROUDMASS/NBLD
C
C  CALCULATE THE SHROUD SURFACE AREA, ASSUMED FROM LE TO TE
C  THE CASING AREA IS TAKEN AS THIS MULTILPLED BY A FACTOR OF 1.25
C
      RLE   = R(JLE(NR),K)
      XLE   = X(JLE(NR),K)
      RTE   = R(JTE(NR),K)
      XTE   = X(JTE(NR),K)
      DIFR  = RTE - RLE
      DIFX  = XTE - XLE
      DMER  = SQRT(DIFR*DIFR + DIFX*DIFX)
      SAREA = PI*(RLE+RTE)*DMER
      SAREA = SAREA * 1.1
      CAREA = SAREA * 1.25
      CFS   = CFSHROUD(NR)
      CFC   = CFCASING(NR)
      RAVG  = 0.5*(RLE+RTE)
      VTSHROUD    = WRAD(JSHROUD)*RAVG
      VTCASE      = WCASE(NR)*RAVG
      VTFLOW      = VTLEAKOUT
      SHROUDWORK  = 0.0
      CASEWORK    = 0.0
      VLIM        = 2.*AMAX1(ABS(VTFLOW),ABS(VTSHROUD),ABS(VTCASE))
C
      N_INTERVALS = 20
C
      FACS    =    0.5*CFS*SAREA*ROSHROUD/N_INTERVALS
      FACC    =    0.5*CFC*CAREA*ROSHROUD/N_INTERVALS
      DO NS   = 1,N_INTERVALS
      VTRELSHRD   = VTFLOW - VTSHROUD
      VTRELCASE   = VTFLOW - VTCASE
      SHROUDFORCE = FACS*VTRELSHRD*ABS(VTRELSHRD)
      CASEFORCE   = FACC*VTRELCASE*ABS(VTRELCASE)
      VTFLOW      = VTFLOW - (SHROUDFORCE + CASEFORCE)/SHROUDMASS
      IF(ABS(VTFLOW).GT.VLIM) VTFLOW = VLIM*VTFLOW/ABS(VTFLOW)
      SHROUDWORK  = SHROUDWORK + SHROUDFORCE*VTSHROUD
      CASEWORK    = CASEWORK   + CASEFORCE*VTCASE
      END DO
C
      SHROUDWORK  = SHROUDWORK + CASEWORK
      VTLEAKIN    = VTFLOW
      HOLEAKIN    = HOLEAKOUT  - SHROUDWORK/SHROUDMASS
      HLIM = 1.25*HOLEAKOUT
      IF(HOLEAKIN.GT.HLIM) HOLEAKIN = HLIM
      SWORKTOT    = SWORKTOT   + SHROUDWORK
C
C     NOTE "sworktot" IS THE FRICTIONAL WORK DONE on THE SHROUD BY
C     THE LEAKAGE FLOW.
C
C   CALCULATE THE PROPERTIES OF THE RENTRY FLOW
C
      ROVNIN =  SHROUDFLBL/(0.25*ALEAKIN)
      VNIN   =  ROVNIN/ROSHROUD
      IF(VNIN.GT.VJET)  VNIN = VJET
      VSIN   =  VNIN/ABS(TAN(PITCHIN(NR)))
      IF(K.EQ.KM) THEN
      VXIN = VSIN*DIFX/DMER + VNIN*DIFR/DMER
      VRIN = VSIN*DIFR/DMER - VNIN*DIFX/DMER
      ENDIF
      IF(K.EQ.1) THEN
      VXIN = VSIN*DIFX/DMER - VNIN*DIFR/DMER
      VRIN = VSIN*DIFR/DMER + VNIN*DIFX/DMER
      ENDIF
C
C    SET THE WALL FLUXES FOR LEAKAGE AND RE-ENTRY FLOW
C
      FK  = 1.0
      RFSHRD1 = 1.-RFSHRD
      IF(K.EQ.1) FK = -1.0
      FACGR    = FK*SHROUDFLBL/ALEAKOUT
      RIN      = 0.5*(R(J4,K) + R(J3,K))
      ROUT     = 0.5*(R(J1,K) + R(J2,K))
      RVTLKIN  = RIN*VTLEAKIN
      RVTLKOUT = ROUT*VTLEAKOUT
      DO 410 J = J1+1,J2
      SFLOW    = 0.0
      DO 400 I = 1,IMM1
      GROUT         = FACGR*AS(I,J)
      RFGROUT       = RFSHRD*GROUT
      SFLOW         = SFLOW + ABS(GROUT)
      SHROUDGR(I,J) = RFSHRD1*SHROUDGR(I,J)  + RFGROUT
      SHROUDHO(I,J) = RFSHRD1*SHROUDHO(I,J)  + RFGROUT*HOLEAKOUT
      SHROUDVX(I,J) = RFSHRD1*SHROUDVX(I,J)  + RFGROUT*VXLEAKOUT
      SHROUDRVT(I,J)= RFSHRD1*SHROUDRVT(I,J) + RFGROUT*RVTLKOUT
      SHROUDVR(I,J) = RFSHRD1*SHROUDVR(I,J)  + RFGROUT*VRLEAKOUT
  400 CONTINUE
      SHRDFLOW(J)   = SHRDFLOW(J-1) + SFLOW
  410 CONTINUE
C
C     SET THE SHROUD FLOW ABOVE/BELOW THE SHROUDED PART OF THE BLADE.
C
      IF(J3.GT.J2) THEN
      DO 450 J = J2+1,J3
      SHRDFLOW(J) = SHROUDFLBL
  450 CONTINUE
      ENDIF
      IF(J4.LT.J1) THEN
      DO 455 J = J4+1,J1
      SHRDFLOW(J) = -SHROUDFLBL
  455 CONTINUE
      ENDIF
C
      FACGR    = -FK*SHROUDFLBL/ALEAKIN
      DO 510 J = J3+1,J4
      SFLOW    = 0.0
      DO 500 I=1,IMM1
      GROUT         = FACGR*AS(I,J)
      RFGROUT       = RFSHRD*GROUT
      SFLOW         = SFLOW + ABS(GROUT)
      SHROUDGR(I,J) = RFSHRD1*SHROUDGR(I,J)  + RFGROUT
      SHROUDHO(I,J) = RFSHRD1*SHROUDHO(I,J)  + RFGROUT*HOLEAKIN
      SHROUDVX(I,J) = RFSHRD1*SHROUDVX(I,J)  + RFGROUT*VXIN
      SHROUDRVT(I,J)= RFSHRD1*SHROUDRVT(I,J) + RFGROUT*RVTLKIN
      SHROUDVR(I,J) = RFSHRD1*SHROUDVR(I,J)  + RFGROUT*VRIN
  500 CONTINUE
      SHRDFLOW(J)   = SHRDFLOW(J-1) - SFLOW
  510 CONTINUE
C
      SLEAK(NR) = SHROUDMASS
      SWORK(NR) = SHROUDWORK
C
C    END OF LOOP OVER ALL BLADE ROWS
C
  100 CONTINUE
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE SHROUDFLUX(SFLUX)
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION SFLUX(ID,JD)
C
      DO 9555 NR =1,NROWS
      IF(KTIPS(NR).GE.0) GO TO 9555
      JS = JSTART(NR)
      K  = KSHROUD(NR)
      J1 = JLEAKS(NR) + JS -1
      J2 = JLEAKE(NR) + JS -1
      J3 = JLKINS(NR) + JS -1
      J4 = JLKINE(NR) + JS -1
      DO 9554 J=J1,J2
      DO 9554 I=1,IMM1
 9554 RFLUX(I,J,K) = RFLUX(I,J,K) + SFLUX(I,J)
      DO 9553 J=J3,J4
      DO 9553 I=1,IMM1
 9553 RFLUX(I,J,K) = RFLUX(I,J,K) + SFLUX(I,J)
 9555 CONTINUE
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE COOLIN_1
C
C  THIS SUBROUTINE SETS UP FIXED COOLING FLOWS WHICH ARE ONLY SET ONCE AND NOT
C  CHANGED DURING THE ITERATIONS. THE MACH MUMBER OF THE EJECTED FLOW IS SPECIFIED
C  AND THE STAGNATION PRESSURE OF THE COOLANT IS ONLY USED TO CALCULATE THE EFFICIENCY.
C
      INCLUDE  'commall-open-18.3'
C
      WRITE(6,*) ' ENTERED COOLIN_1'
C
      DO N = 1,NSTAGES
           WPUMP(N) = 0.0
      END DO
C
C*************************************************************************
C*************************************************************************
C
C      SET INITIAL COOLING FLOWS TO ZERO AND
C      WORK OUT AREA MAGNITUDES FOR COOLANT FLOWS
C
      DO 400 J=2,JM
      DO 400 K=1,KMM1
      CFLOWI1(J,K) = 0.0
      CFLOWIM(J,K) = 0.0
      ATOTI1(J,K) = SQRT(ABX(1,J,K)*ABX(1,J,K) + ABR(1,J,K)*ABR(1,J,K)
     &                  + ABT(J,K)*ABT(J,K))
      ATOTIM(J,K) = SQRT(ABX(IM,J,K)*ABX(IM,J,K)+ABR(IM,J,K)*ABR(IM,J,K)
     &                  + ABT(J,K)*ABT(J,K))
  400 CONTINUE
C
      DO 500 J=2,JM
      DO 500 I=1,IMM1
      CFLOWK1(I,J) = 0.0
      CFLOWKM(I,J) = 0.0
      ATOTK1(I,J)= SQRT(ASX(I,J,1)*ASX(I,J,1)  + ASR(I,J,1)*ASR(I,J,1))
      ATOTKM(I,J)= SQRT(ASX(I,J,KM)*ASX(I,J,KM)+ASR(I,J,KM)*ASR(I,J,KM))
  500 CONTINUE
C
C***********************************************************************
C***********************************************************************
C
C      NOW EVALUATE THE COOLING MASS FLOWS THROUGH THE BLADE SURFACES.
C      AND THE ASSOCIATED ENERGY AND MOMENTUM FLUXES.
C
C
      DO 600 NCB = 1,NCOOLB
C
      ACOOLB(NCB) = 0.0
      DO 610 J = JCBS(NCB)+1,JCBE(NCB)
      DO 610 K = KCBS(NCB),KCBE(NCB)-1
      IF(IC(NCB).EQ.1)  ACOOLB(NCB) = ACOOLB(NCB) + ATOTI1(J,K)
      IF(IC(NCB).EQ.IM) ACOOLB(NCB) = ACOOLB(NCB) + ATOTIM(J,K)
  610 CONTINUE
C
      DO 620 J = JCBS(NCB)+1,JCBE(NCB)
C
             N_ROW    = NROW(J)
             N_STAGE  = 1 + (N_ROW-1)/2
C
      DO 620 K = KCBS(NCB),KCBE(NCB)-1
C
      IF(IC(NCB).EQ.1)  THEN
           CFLOWIN  = CFLOWB(NCB)*ATOTI1(J,K)/(ACOOLB(NCB)*NBLADE(J))
      ELSE
           CFLOWIN  = CFLOWB(NCB)*ATOTIM(J,K)/(ACOOLB(NCB)*NBLADE(J))
      END IF
C
      TOIN     = TOCOOLB(NCB)
      POIN     = POCOOLB(NCB)
      RAVG     = RAVG_CELL(J,K)
      VBLADE   = WRAD_COOLB(NCB)*RAVG
C  
C     PRE SWIRL INCLUDED 4/12/06  JDD
C
      W_PRESWIRL     = WRAD_COOLB(NCB)*RVT_IN_B(NCB)
      WPUMP(N_STAGE) = WPUMP(N_STAGE) + NBLADE(J)*
     &                 CFLOWIN*(VBLADE*VBLADE - W_PRESWIRL)
C
      IF(IFGAS.EQ.0) THEN
           TOREL       = TOIN  + (0.5*VBLADE*VBLADE - W_PRESWIRL)/CP
           PO_REL      = POIN*(TOREL/TOIN)**RFGA
           TOABS       = TOREL + 0.5*VBLADE*VBLADE/CP
           TIN         = TOREL/(1. + 0.5*(GA-1.)*MACHCOOL(NCB)*
     &                   MACHCOOL(NCB))
           VIN         = SQRT(2.*CP*(TOREL - TIN))
           HOABS       = CP*TOABS
      ELSE
           DELT        = TOIN-TREF
           CPGAS       = CP1 + CP2*DELT + CP3*DELT*DELT
           GAGAS       = CPGAS/(CPGAS - RGAS)
           RFGAGAS     = GAGAS/(GAGAS - 1.)
           HO_IN       = HFROMT(TOIN,TREF,HREF,CP1,CP2,CP3)
           HOREL       = HO_IN  + 0.5*VBLADE*VBLADE - W_PRESWIRL
           HOABS       = HOREL + 0.5*VBLADE*VBLADE
           TOREL       = TFROMH(HOREL,TREF,HREF,HT1,HT2,HT3,HT4)
           PO_REL      = POIN*(TOREL/TOIN)**RFGAGAS
           TIN         = TOREL/(1. + 0.5*(GAGAS-1.)*MACHCOOL(NCB)*
     &                   MACHCOOL(NCB))
           HSTAT       = HFROMT(TIN,TREF,HREF,CP1,CP2,CP3) 
           IF(HSTAT.GT.HOREL) HSTAT = 0.99999*HOREL
           VIN         = SQRT(2.0*(HOREL - HSTAT))
      END IF         
C
      TORELB(NCB) = TOREL
      TSTAT_EJECTB(NCB) = TIN
      PO_EJECTB(NCB)    = PO_REL
      VNIN     = VIN*SIN(SANGLEB(NCB))
      VSIN     = VIN*COS(SANGLEB(NCB))
C
      IF(IC(NCB).EQ.1)  THEN      
           ANORM    = ATOTI1(J,K)
           XNORM    = ABX(1,J,K)/ANORM
           TNORM    = ABT(J,K)/ANORM
           RNORM    = ABR(1,J,K)/ANORM
      ELSE
           ANORM    = ATOTIM(J,K)
           XNORM    = ABX(IM,J,K)/ANORM
           TNORM    = ABT(J,K)/ANORM
           RNORM    = ABR(IM,J,K)/ANORM
      END IF
C
      ANG      =  XANGLEB(NCB)
      XTNORM   =  SQRT(XNORM*XNORM + TNORM*TNORM)
      SX       =  (TNORM*COS(ANG) - RNORM*XNORM*SIN(ANG))/XTNORM
      ST       = -(XNORM*COS(ANG) + RNORM*TNORM*SIN(ANG))/XTNORM
      SR       =  SIN(ANG)*XTNORM
C
      IF(IC(NCB).EQ.1)  THEN
           XVEL          =  VNIN*XNORM + VSIN*SX
           RVEL          =  VNIN*RNORM + VSIN*SR
           TVEL          =  VNIN*TNORM + VSIN*ST
           CFLOWI1(J,K)  = CFLOWIN
           HOCWLI1(J,K)  = CFLOWIN*(HOABS + VBLADE*TVEL)
           VXCWLI1(J,K)  = CFLOWIN*XVEL
           VRCWLI1(J,K)  = CFLOWIN*RVEL
           RVTCWLI1(J,K) = CFLOWIN*(VBLADE  +  TVEL)*RAVG
      ELSE
           XVEL          = -VNIN*XNORM + VSIN*SX
           RVEL          = -VNIN*RNORM + VSIN*SR
           TVEL          = -VNIN*TNORM + VSIN*ST
           CFLOWIM(J,K)  = CFLOWIN
           HOCWLIM(J,K)  = CFLOWIN*(HOABS + VBLADE*TVEL)
           VXCWLIM(J,K)  = CFLOWIN*XVEL
           VRCWLIM(J,K)  = CFLOWIN*RVEL
           RVTCWLIM(J,K) = CFLOWIN*(VBLADE  +  TVEL)*RAVG
      END IF
C
  620 CONTINUE
C
C     END OF BLADE SURFACE COOLING SETUP
C
  600 CONTINUE
C
C**************************************************************************
C**************************************************************************
C
C      NOW EVALUATE THE COOLING MASS FLOWS THROUGH THE HUB AND CASING.
C
      DO 700 NCW = 1,NCOOLW
C
      ACOOLW(NCW) = 0.0
      DO 710 J =  JCWS(NCW)+1,JCWE(NCW)
      DO 710 I =  ICWS(NCW),ICWE(NCW)-1
      IF(KC(NCW).EQ.1)  ACOOLW(NCW) = ACOOLW(NCW) + ATOTK1(I,J)
      IF(KC(NCW).EQ.KM) ACOOLW(NCW) = ACOOLW(NCW) + ATOTKM(I,J)
  710 CONTINUE
C
      DO 720 J = JCWS(NCW)+1,JCWE(NCW)
C
           N_ROW    = NROW(J)
           N_STAGE  = NSTAGE(N_ROW)
C
      DO 720 I = ICWS(NCW),ICWE(NCW)-1
C
      IF(KC(NCW).EQ.1) THEN
           CFLOWIN  = CFLOWW(NCW)*ATOTK1(I,J)/(ACOOLW(NCW)*NBLADE(J))
           RAVG     = 0.5*(R(J,1)+R(J-1,1))
      ELSE
           CFLOWIN  = CFLOWW(NCW)*ATOTKM(I,J)/(ACOOLW(NCW)*NBLADE(J))
           RAVG     = 0.5*(R(J,KM)+R(J-1,KM))
      END IF 
C
      TOIN   = TOCOOLW(NCW)
      POIN   = POCOOLW(NCW)
      VWALL  = WRAD_COOLW(NCW)*RAVG
C
C     PRE SWIRL INCLUDED 4/12/06  JDD
C
      W_PRESWIRL     = WRAD_COOLW(NCW)*RVT_IN_W(NCW)
      WPUMP(N_STAGE) = WPUMP(N_STAGE) + NBLADE(J)*
     &                 CFLOWIN*(VWALL*VWALL - W_PRESWIRL)
C
      IF(IFGAS.EQ.0) THEN
           TOREL  = TOIN  + (0.5*VWALL*VWALL - W_PRESWIRL)/CP
           PO_REL = POIN*(TOREL/TOIN)**RFGA
           TOABS  = TOREL + 0.5*VWALL*VWALL/CP
           TIN    = TOREL/(1. + 0.5*(GA-1.)*MACHCOOL(NCW)*MACHCOOL(NCW))
           VIN    = SQRT(2.*CP*(TOREL - TIN))
           HOABS  = CP*TOABS
      ELSE
           DELT        = TOIN-TREF
           CPGAS       = CP1 + CP2*DELT + CP3*DELT*DELT
           GAGAS       = CPGAS/(CPGAS - RGAS)
           RFGAGAS     = GAGAS/(GAGAS - 1.)
           HO_IN       = HFROMT(TOIN,TREF,HREF,CP1,CP2,CP3)
           HOREL       = HO_IN  + 0.5*VWALL*VWALL - W_PRESWIRL
           HOABS       = HOREL + 0.5*VWALL*VWALL
           TOREL       = TFROMH(HOREL,TREF,HREF,HT1,HT2,HT3,HT4)
           PO_REL      = POIN*(TOREL/TOIN)**RFGAGAS
           TIN         = TOREL/(1. + 0.5*(GAGAS-1.)*MACHCOOL(NCB)*
     &                   MACHCOOL(NCB)) 
           HSTAT       = HFROMT(TIN,TREF,HREF,CP1,CP2,CP3) 
           IF(HSTAT.GT.HOREL) HSTAT = 0.99999*HOREL
           VIN         = SQRT(2.0*(HOREL - HSTAT))
      END IF         
C
      TORELW(NCW) = TOREL
      TSTAT_EJECTW(NCW) = TIN
      PO_EJECTW(NCW)    = PO_REL
      VNIN   = VIN*SIN(SANGLEW(NCW))
      VSIN   = VIN*COS(SANGLEW(NCW))
C
      IF(KC(NCW).EQ.1) THEN
           ANORM  = ATOTK1(I,J)
           XNORM  = ASX(I,J,1)/ANORM
           RNORM  = ASR(I,J,1)/ANORM
      ELSE
           ANORM  = ATOTKM(I,J)
           XNORM  = ASX(I,J,KM)/ANORM
           RNORM  = ASR(I,J,KM)/ANORM
      END IF
C
      ANG    =  TANGLEW(NCW)
      XRNORM =  SQRT(XNORM*XNORM + RNORM*RNORM)
      SX     =  COS(ANG)*RNORM/XRNORM
      ST     =  SIN(ANG)
      SR     = -COS(ANG)*XNORM/XRNORM

C
      IF(KC(NCW).EQ.1) THEN
           XVEL   =  VNIN*XNORM + VSIN*SX
           RVEL   =  VNIN*RNORM + VSIN*SR
           TVEL   =  VSIN*ST
      ELSE
           XVEL   =  -VNIN*XNORM + VSIN*SX
           RVEL   =  -VNIN*RNORM + VSIN*SR
           TVEL   =   VSIN*ST
      END IF
C
      IF(KC(NCW).EQ.1) THEN    
           CFLOWK1(I,J)  = CFLOWIN
           HOCWLK1(I,J)  = CFLOWIN*(HOABS + VWALL*TVEL)
           VXCWLK1(I,J)  = CFLOWIN*XVEL
           VRCWLK1(I,J)  = CFLOWIN*RVEL
           RVTCWLK1(I,J) = CFLOWIN*(VWALL + TVEL)*RAVG
      ELSE
           CFLOWKM(I,J)  = CFLOWIN
           HOCWLKM(I,J)  = CFLOWIN*(HOABS + VWALL*TVEL)
           VXCWLKM(I,J)  = CFLOWIN*XVEL
           VRCWLKM(I,J)  = CFLOWIN*RVEL
           RVTCWLKM(I,J) = CFLOWIN*(VWALL + TVEL)*RAVG
      ENDIF
C
  720 CONTINUE
C
  700 CONTINUE
C
C     FIND THE SUM OF THE COOLING FLOWS ADDED UP TO EACH J STATION
C
      SUMCWL(1) = 0.0
      DO 20 J=2,JM
      COOLADD   = 0.0
      DO 25 K=1,KMM1
   25 COOLADD   = COOLADD  + (CFLOWI1(J,K) + CFLOWIM(J,K))*NBLADE(J)
      DO 30 I=1,IMM1
   30 COOLADD   = COOLADD  + (CFLOWK1(I,J) + CFLOWKM(I,J))*NBLADE(J)
      SUMCWL(J) = SUMCWL(J-1) + COOLADD
   20 CONTINUE
C
C     END OF COOLING FLOW SETUP
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE BLEEDOUT(IBLEED,KBLEED)
C
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION 
     & KBLEEDS(NBLED),KBLEEDE(NBLED),JBLEEDS(NBLED),JBLEEDE(NBLED),
     & IBLEEDS(NBLED),IBLEEDE(NBLED),MBLEED(NBLED)
C
      REAL MASSBLED,MBLEED
C
   99 FORMAT(A72)
      NBLEEDTOT = 0
      IBLEED    = 0
      KBLEED    = 0
C
      DO 100 N_ROW = 1,NROWS
C
      J1 = JSTART(N_ROW)
      READ(5,99)  DUMMY_INPUT
      WRITE(6,99) DUMMY_INPUT
      READ(5,*)   NBLEED
C
C    READ IN THE BLEED FLOW DETAILS. THE J VALUES IN THE INPUT ARE
C    DEFINED FOR THE CURRENT BLADE ROW. IE J = 1 AT THE UPSTREAM MIXING PLANE.
C
      DO 10 NBLD =1,NBLEED
      NBLEEDTOT = NBLEEDTOT + 1
      READ(5,*) IBLDS,IBLDE,JBLDS,JBLDE,KBLDS,KBLDE,MASSBLED
      KBLEEDS(NBLEEDTOT)  = KBLDS
      KBLEEDE(NBLEEDTOT)  = KBLDE
      IBLEEDS(NBLEEDTOT)  = IBLDS
      IBLEEDE(NBLEEDTOT)  = IBLDE
      JBLEEDS(NBLEEDTOT)  = J1 + JBLDS
      JBLEEDE(NBLEEDTOT)  = J1 + JBLDE
      MBLEED(NBLEEDTOT)   = MASSBLED
   10 CONTINUE
C
  100 CONTINUE
C
C      WORK OUT AREA MAGNITUDES FOR BLEED FLOWS
C
      DO 150 J=2,JM
      DO 150 K=1,KMM1
      BLFLOWI1(J,K) = 0.0
      BLFLOWIM(J,K) = 0.0
      ATOTI1(J,K) = SQRT(ABX(1,J,K)*ABX(1,J,K) + ABR(1,J,K)*ABR(1,J,K)
     &                 + ABT(J,K)*ABT(J,K))
      ATOTIM(J,K) = SQRT(ABX(IM,J,K)*ABX(IM,J,K)+ABR(IM,J,K)*ABR(IM,J,K)
     &                 + ABT(J,K)*ABT(J,K))
  150 CONTINUE
C
      DO 500 J=2,JM
      DO 500 I=1,IMM1
      BLFLOWK1(I,J) = 0.0
      BLFLOWKM(I,J) = 0.0
      ATOTK1(I,J) =SQRT(ASX(I,J,1)*ASX(I,J,1) + ASR(I,J,1)*ASR(I,J,1))
      ATOTKM(I,J) =SQRT(ASX(I,J,KM)*ASX(I,J,KM)+ASR(I,J,KM)*ASR(I,J,KM))
  500 CONTINUE
C
C
      DO 200 NBLD =1,NBLEEDTOT
C
      KS = KBLEEDS(NBLD)
      KE = KBLEEDE(NBLD)
      JS = JBLEEDS(NBLD)
      JE = JBLEEDE(NBLD)
      IS = IBLEEDS(NBLD)
      IE = IBLEEDE(NBLD)
C
      IF(KS.EQ.KE.AND. (KS.NE.0) ) THEN
C
C     IF HUB OR CASING BLEED
C
      K = KS
      KBLEED = 1
      ABLEED = 0.0
      DO 210 J = JS+1,JE
      DO 210 I=IS,IE-1
      IF(K.EQ.1)  ABLEED = ABLEED + ATOTK1(I,J)
      IF(K.EQ.KM) ABLEED = ABLEED + ATOTKM(I,J)
  210 CONTINUE
      DO 300 J = JS+1,JE
      NBLAD = NBLADE(J)
      DO 300 I = IS,IE-1
      IF(K.EQ.1)   BLFLOWK1(I,J) = MBLEED(NBLD)*ATOTK1(I,J)/ABLEED/NBLAD
      IF(K.EQ.KM)  BLFLOWKM(I,J) = MBLEED(NBLD)*ATOTKM(I,J)/ABLEED/NBLAD
  300 CONTINUE
C
      ENDIF
C
      IF(IS.EQ.IE.AND. (IS.NE.0) ) THEN
C
C     IF BLADE SURFACE BLEED
C
      I = IS
      IBLEED = 1
      ABLEED = 0.0
      DO 250 J = JS+1,JE
      DO 250 K = KS,KE-1
      IF(I.EQ.1)  ABLEED = ABLEED + ATOTI1(J,K)
      IF(I.EQ.IM) ABLEED = ABLEED + ATOTIM(J,K)
  250 CONTINUE
      DO 270 J = JS+1,JE
      NBLAD = NBLADE(J)
      DO 270 K = KS,KE-1
      IF(I.EQ.1)  BLFLOWI1(J,K) = MBLEED(NBLD)*ATOTI1(J,K)/ABLEED/NBLAD
      IF(I.EQ.IM) BLFLOWIM(J,K) = MBLEED(NBLD)*ATOTIM(J,K)/ABLEED/NBLAD
  270 CONTINUE
C
      ENDIF
C
C     OBTAIN THE TOTAL BLED FLOW FOR USE IN THE CONTINUITY CHECK
C
      SUMBLEED(1) = 0.0
      DO 400 J=2,JM
      NBLAD  = NBLADE(J)
      BLEEDJ = 0.0
      DO 410 I=1,IMM1
  410 BLEEDJ       = BLEEDJ + NBLAD*(BLFLOWK1(I,J) + BLFLOWKM(I,J))
      DO 420 K=1,KMM1
  420 BLEEDJ       = BLEEDJ + NBLAD*(BLFLOWI1(J,K) + BLFLOWIM(J,K))
      SUMBLEED(J)  = SUMBLEED(J-1) + BLEEDJ
  400 CONTINUE
C
C
  200 CONTINUE
C
      RETURN
      END
C
C**********************************************************************
C
      SUBROUTINE EFICOOL(HBLEDTOT)
C
C     THIS SUBROUTINE CALCULATES THE MASS AVERAGED QUANTITIES AND MACHINE
C     EFFICIENCY AND WRITES THEM OUT TO UNIT 6.
C
      INCLUDE  'commall-open-18.3'
C
C
      DIMENSION BLADE_FLOW(JD),SUM_ENTPY(JD),SUMPO(JD),SUMTO(JD),
     &          SUMRVT(JD),SUMTSTAT(JD), ETA_LOSS(JD),SUMPSTAT(JD),
     &          CHECK_FLOW(JD)
C
C      SUM FLUXES OF MASS,STAGNATION PRESURE  ETC EVERY 100 TIME STEPS
C
C      SET TEMP1, TEMP2, TEMP3 ,TEMP4 ,STORE2 TO STAGN TEMP, STAGN PRESS, RVT, ENTROPY,
C      AND STATIC TEMPERATURE. FOR USE IN SUBROUTINE SUMFLX.
C
C    TEMP3  stores RVT
C    TEMP2  stores the stagnation pressure ratio
C    TEMP1  stores the stagnation temperature
C    TEMP4 stores the entropy.
C    STORE2 stores the static temperature
C
      ENTHIN = 0.0
C
      DO 5610 K=1,KM
      DO 5610 J=1,JM
      DO 5610 I=1,IM
      EKE  = 0.5*(VX(I,J,K)*VX(I,J,K) + VT(I,J,K)*VT(I,J,K)
     &          + VR(I,J,K)*VR(I,J,K))
C
      IF(IFGAS.EQ.0) THEN
           TS    = (HO(I,J,K) - EKE)*RCP
           TSTAG = HO(I,J,K)*RCP
      ELSE
           HSTAT  = HO(I,J,K) - EKE
           TS     = TFROMH(HSTAT,TREF,HREF,HT1,HT2,HT3,HT4)
           HSTAG  = HO(I,J,K)
           TSTAG  = TFROMH(HSTAG,TREF,HREF,HT1,HT2,HT3,HT4)
      END IF
C
      IF(TSTAG.LT.0.1*TO1(KMID))  TSTAG = 0.1*TO1(KMID)
C
      IF(TS.LT.0.1*TO1(KMID))      TS   = 0.1*TO1(KMID)
C
      PSTAT = P(I,J,K)
      IF(PSTAT.LT.0.02*PO1(KMID)) PSTAT=0.02*PO1(KMID)
C
      IF(IFGAS.EQ.0) THEN
            PSTAG = PSTAT*(TSTAG/TS)**RFGA
      ELSE
            TRAT  = TSTAG/TS
            PRAT  = PRAT_FROM_TRAT(TS,TRAT,ALPHA,BETA1,BETA2)
            PSTAG = PSTAT*PRAT
      END IF
C
      TRAT    = TS/TO1(KMID)
      PRAT    = PSTAT/PO1(KMID)
      ENTRPY  = 0.0
      IF(PRAT.GT.0.99999.AND.PRAT.LT.1.00001) GO TO 5605
      IF(TRAT.GT.0.99999.AND.TRAT.LT.1.00001) GO TO 5605
C
      IF(IFGAS.EQ.0) THEN
           ENTRPY = CP*LOG(TRAT)- RGAS*LOG(PRAT)
      ELSE
           PRAT   = PRAT_FROM_TRAT(TO1(KMID),TRAT,ALPHA,BETA1,BETA2)
           ENTRPY = PSTAT/(PO1(KMID)*PRAT)
      END IF
C
 5605 CONTINUE
C
      STORE2(I,J,K)= TS
      TEMP1(I,J,K) = TSTAG
      TEMP2(I,J,K) = PSTAG/PO1(KMID)
      TEMP3(I,J,K) = RORVT(I,J,K)/RO(I,J,K)
      TEMP4(I,J,K) = ENTRPY
C
 5610 CONTINUE
C
      CALL SUMFLX(BLADE_FLOW,SUMPO,SUMTO,SUMRVT,SUM_ENTPY,SUMTSTAT,
     &            SUMPSTAT) 
C
C******************************************************************************
C
      IF(IFEND.EQ.1)  THEN
C
C     WRITE OUTPUT FOR GLOBPLOT ONLY AT END OF WHOLE CALCULATION.
C
      WRITE(11)   1,1
      WRITE(11)   (BLADE_FLOW(J), J=1,JM)
      WRITE(11)   (SUMPO(J), J=1,JM)
      WRITE(11)   (SUMTO(J), J=1,JM)
      WRITE(11)   (SUMRVT(J),J=1,JM)
      WRITE(11)   (SUM_ENTPY(J),J=1,JM)

      ENDIF
C********************************************************
C  MASS AVERAGE THE 1D VALUES
C
      DO 5611 J=1,JM
      TFLOW         = BLADE_FLOW(J)
      SUMPO(J)      = SUMPO(J)/TFLOW
      SUMTO(J)      = SUMTO(J)/TFLOW
      SUMRVT(J)     = SUMRVT(J)/TFLOW
      SUM_ENTPY(J)  = SUM_ENTPY(J)/TFLOW
      SUMTSTAT(J)   = SUMTSTAT(J)/TFLOW
      SUMPSTAT(J)   = SUMPSTAT(J)/TFLOW
 5611 CONTINUE
C
      DTREF  = ABS( SUMTO(1) - SUMTO(JM) )
      IFOUT  = 0
      IF(DTREF.LT.0.005*SUMTO(1)) THEN
           DTIN = SUMTO(1)  - SUMTSTAT(1)
           DTEX = SUMTO(JM) - SUMTSTAT(JM)
           DTREF = AMAX1(DTIN,DTEX)
           IFOUT = 1
      END IF
      DHREF = CP*DTREF
      DO J = 1,JM
      ETA_LOSS(J) = SUMTO(JM)*(SUM_ENTPY(J)-SUM_ENTPY(1))/DHREF
      END DO
C
C ********************************************************
C    WRITE A FILE TO PLOT THE ENTROPY LOSS COEFFICIENT OR LOST EFFICIENCY'
C  
      IF(IFEND.EQ.1) THEN   
      OPEN(UNIT=23,FILE='loss-co.plt')
      WRITE(23,*) ' PLOTTING OUTPUT FOR LOST EFFICIENCY '
      WRITE(23,*) ' NUMBER OF OUTPUT POINTS ', JM
      WRITE(23,*) ' MERIDIONAL DISTANCE '
      WRITE(23,5612) (SMERID(J,KMID), J=1,JM)
      WRITE(23,*) ' LOST EFFICIENCY '
      WRITE(23,5612) (ETA_LOSS(J),   J=1,JM)
 5612 FORMAT(10F10.5)
      END IF
C
C******************************************************************************
C
C     WRITE OUT THE MAIN OUTPUT TO UNIT 6
C
      IF(IFOUT.EQ.0) THEN
           WRITE(6,5645)
           WRITE(6,5646)
      ELSE
           WRITE(6,5647)
           WRITE(6,5648)
      END IF
C
 5645 FORMAT( '    J           MASS          LOCAL/INLET         LOST
     &        STAGNATION      STAGNATION       MASS AVG     ')
 5646 FORMAT( '  VALUE      FLOW. KG/S          FLOW          EFFICIENCY
     &        PRESSURE        TEMPERATURE        R*VT       ')
C
 5647 FORMAT( '    J           MASS          LOCAL/INLET       ENTROPY
     &        STAGNATION      STAGNATION       MASS AVG     ')
 5648 FORMAT( '  VALUE      FLOW. KG/S          FLOW         LOSS COEFF.
     &        PRESSURE        TEMPERATURE        R*VT       ')
C
C
      DO 5640  J=1,JM
      CHECK_FLOW(J) =  BLADE_FLOW(J) + SHRDFLOW(J)*NBLADE(J)
     &               +  SUMBLEED(J) - SUMCWL(J)
      RATIO          =  CHECK_FLOW(J)/CHECK_FLOW(1)
      WRITE(6,5630) J,CHECK_FLOW(J),RATIO,ETA_LOSS(J),SUMPO(J),
     &                SUMTO(J),SUMRVT(J)
 5630 FORMAT( I5, 6F17.5)
 5640 CONTINUE
C
C************************************************************************************
C
C     WRITE OUT THE REYNOLDS NUMBER FOR EACH BLADE ROW
C
      WRITE(6,*)
      WRITE(6,*)
      DO NRW = 1,NROWS
           JLEDGE  = JLE(NRW)
           JTEDGE  = JTE(NRW)
           ROWREF = SQRT(ROVX(IMID,JTEDGE,KMID)*ROVX(IMID,JTEDGE,KMID) +
     &                  ROVR(IMID,JTEDGE,KMID)*ROVR(IMID,JTEDGE,KMID) 
     &                + ROWT(IMID,JTEDGE,KMID)*ROWT(IMID,JTEDGE,KMID)  )
           REYNOLDS(NRW) = CHORD(NRW)*ROWREF/VISCOSY(NRW)
      WRITE(6,*) ' BLADE ROW NUMBER ',NRW
      WRITE(6,*) ' CHORD = ',CHORD(NRW),' ROV-EXIT =',ROWREF,
     &           ' VISCOSITY = ',VISCOSY(NRW)
      WRITE(6,*) ' REYNOLDS NUMBER BASED ON BLADE EXIT CONDITIONS = ',
     &             REYNOLDS(NRW)
      END DO
      WRITE(6,*)
C
C***************************************************************************************************
C***************************************************************************************************
C
C     WORK OUT THE EFFICIENCY OF THE WHOLE MACHINE IF THERE ARE COOLING FLOWS PRESENT
C
      IF(IFCOOL.NE.0) THEN
C
C     FIND THE TOTAL WORK DONE IN PUMPING THE COOLANT 
C
      WPUMP_TOT = 0.0
      DO 145 N  = 1,NSTAGES
      WPUMP_TOT = WPUMP_TOT + WPUMP(N)
  145 CONTINUE
      WPUMP(NSTAGES+1) = WPUMP_TOT
C
      NSTGP1  = NSTAGES + 1
      IF(NSTAGES.EQ.1) NSTGP1 = 1
C 
C***********************************************************************************************
C    NOW WRITE OUT THE STAGE PERFORMANCES.
C    NOTE:  THE NSTGP1 STAGE IS THE WHOLE MACHINE
C
      DO 3333  NSTG = 1,NSTGP1
C
      WRITE(6,*)
      WRITE(6,*) '******************************************************
     &******'
      IF(NSTG.EQ.NSTGP1) THEN
           WRITE(6,*) ' THE FOLLOWING OUTPUT IS FOR THE WHOLE MACHINE '
      ELSE
           WRITE(6,*) ' THE FOLLOWING OUTPUT IS FOR STAGE NUMBER ',NSTG
      END IF
      WRITE(6,*) '******************************************************
     &******'
C
      IF(NSTG.LE.NSTAGES) THEN
              J1      = JSTG_START(NSTG)
              J2      = JSTG_END(NSTG) + 1
              J2M1    = J2
C
C     JDD CHANGE  2/7/07  TO INCLUDE THE MIXING LOSS IN THE EFFICIENCY OF THE UPSTREAM BLADE ROW.
C     FOR THE LAST STAGE THERE IS NO MIXING DOWNSTREAM. BASE THE EFFY ON THE LAST BUT ONE POINT.
C
              IF(NSTG.EQ.NSTAGES) J2M1 = JM -1 
C
C     END OF JDD ADDITION
C
C     JMIX_STAGE IS ONLY USED TO FIND THE STAGE REACTION.
              N_ROW      = NROW(J1+1)
              JMIX_STAGE = JMIX(N_ROW) 
      ELSE
C
C   FOR THE WHOLE MACHINE THERE IS NO MIXING DOWNSTREAM. BASE THE EFFY ON THE LAST BUT ONE POINT, J2M1.
C
              J1   = 2
              J2   = JM
              J2M1 = JM - 1
              JMIX_STAGE = (J1+J2)/2
      END IF
C
      WRITE(6,*)
      WRITE(6,*) ' EFFICIENCIES BASED ON MASS AVERAGED VALUES AT J = ',
     &             J1,J2M1
      WRITE(6,*)
C
C     FIND THE STAGE INLET AND OUTLET STAGNATION CONDITIONS
C
      POIN    = SUMPO(J1)*PO1(KMID)
      PO2     = SUMPO(J2M1)*PO1(KMID)
      PORAT   = POIN/PO2
      TO_IN   = SUMTO(J1)
      TO_OUT  = SUMTO(J2M1)
      PSIN    = SUMPSTAT(J1)
      PSOUT   = SUMPSTAT(J2M1)
C
C     IF THERE ARE ANY COOLING FLOWS, SUM THE COOLING MASS FLOWS AND ENTHALPIES
C     AND FIND THE ISENTROPIC WORK AVAILABLE FROM THE COOLING FLOWS.
C
      ENTHIN        = 0.0
      WCOOL_IS      = 0.0
      WCOOL_AT_HOLE = 0.0
      CFLOW_STG     = 0.0
C
      DO 150 NCB = 1,NCOOLB
C
      SUM_PSTAT  = 0.0
      SUM_MAS    = 0.0

      IF(JCBS(NCB).LT.J1.AND.JCBE(NCB).LE.J1) GO TO 150
      IF(JCBS(NCB).GE.J2.AND.JCBE(NCB).GT.J2) GO TO 150

      DO 151 J = JCBS(NCB)+1,JCBE(NCB)

      IF(J.LE.J1.OR.J.GT.J2) GO TO 151

      DO 152 K = KCBS(NCB),KCBE(NCB)-1
      IF(IC(NCB).EQ.1) THEN
           SUM_PSTAT  = SUM_PSTAT + CFLOWI1(J,K)*
     &     (P(1,J-1,K)+P(1,J,K)+P(1,J-1,K+1)+P(1,J,K+1))
           SUM_MAS = SUM_MAS + CFLOWI1(J,K)
      ELSE
           SUM_PSTAT  = SUM_PSTAT + CFLOWIM(J,K)*
     &     (P(IM,J-1,K)+P(IM,J,K)+P(IM,J-1,K+1)+P(IM,J,K+1))
           SUM_MAS = SUM_MAS + CFLOWIM(J,K)
      ENDIF
  152 CONTINUE

  151 CONTINUE
C
C     "PS_EJECT"  IS THE AVERAGE STATIC PRESSURE ON THE SURFACE OF THE COOLING PATCH.
C     "PO_AT_HOLE"  IS THE ABSOLUTE STAGNATION PRESSURE AT WHICH THE COOLING 
C     FLOW WOULD NEED TO BE SUPPLIED TO THE BLADE ROW IF THERE WAS NO LOSS IN
C     THE INTERNAL COOLING PASSAGES.
C
      PS_EJECT  = 0.25*SUM_PSTAT/SUM_MAS
C
      IF(IFGAS.EQ.0) THEN
           GAGAS    = GA
      ELSE
           DELT   = TORELB(NCB) - TREF
           CPGAS  = CP1 + CP2*DELT + CP3*DELT*DELT
           GAGAS  = CPGAS/(CPGAS - RGAS)
      END IF
C
C      TS_EJECT  = TORELB(NCB)
C     &          /(1.+.5*(GAGAS-1)*MACHCOOL(NCB)*MACHCOOL(NCB))
C     JDD CHANGED TO USE STORED VALUE OF TS_EJECT. MAY 2018 .
      TS_EJECT = TSTAT_EJECTB(NCB)
C
      IF(IFGAS.EQ.0) THEN
           PO_AT_HOLE= PS_EJECT*(TOCOOLB(NCB)/TS_EJECT)**RFGA
      ELSE
           TRAT = TOCOOLB(NCB)/TS_EJECT
           PRAT = PRAT_FROM_TRAT(TS_EJECT,TRAT,ALPHA,BETA1,BETA2)
           PO_AT_HOLE= PS_EJECT*PRAT
      END IF
C
      POCOOL     = POCOOLB(NCB)
      TOCOOL     = TOCOOLB(NCB)
      SUM_MAS    = SUM_MAS*NBLADE(JCBS(NCB)+1)
      CFLOW_STG  = CFLOW_STG + SUM_MAS
C
      IF(IFGAS.EQ.0) THEN
           TO2_IS         = TOCOOL*(PO2/POCOOL)**FGA
           WCOOL_IS       = WCOOL_IS + SUM_MAS*CP*(TOCOOL-TO2_IS)
           TO2IS_AT_HOLE  = TOCOOL*(PO2/PO_AT_HOLE)**FGA
           WCOOL_AT_HOLE  = WCOOL_AT_HOLE 
     &                    + SUM_MAS*CP*(TOCOOL - TO2IS_AT_HOLE)
           ENTHIN         = ENTHIN   + SUM_MAS*CP*TOCOOL
      ELSE
           PRAT       = PO2/POCOOL
           TRAT       = TRAT_FROM_PRAT(TOCOOL,PRAT,FGAGAS,R_ALPHA,
     &                  BALPHA1,BALPHA2)
           TO2_IS     = TOCOOL*TRAT
           HOCOOL_IS  = HFROMT(TO2_IS,TREF,HREF,CP1,CP2,CP3)
           HOCOOL     = HFROMT(TOCOOL,TREF,HREF,CP1,CP2,CP3)
           WCOOL_IS   = WCOOL_IS + SUM_MAS*(HOCOOL - HOCOOL_IS)
           PRAT       = PO2/PO_AT_HOLE
           TRAT       = TRAT_FROM_PRAT(TOCOOL,PRAT,FGAGAS,R_ALPHA,
     &                  BALPHA1,BALPHA2)
           TO2IS_AT_HOLE = TOCOOL*TRAT
           HO2IS_AT_HOLE = HFROMT(TO2IS_AT_HOLE,TREF,HREF,CP1,CP2,CP3)
           WCOOL_AT_HOLE = WCOOL_AT_HOLE 
     &                   + SUM_MAS*(HOCOOL - HO2IS_AT_HOLE)
           ENTHIN        = ENTHIN + SUM_MAS*HOCOOL
      END IF
C
      WRITE(6,164) NCB,POCOOL,PO_EJECTB(NCB),PO_AT_HOLE
  164 FORMAT(' BLADE COOLING PATCH NUMBER',I3,' SPECIFIED SUPPLY STAGNAT 
     &ION PRESSURE = ',T67,F10.1,/,' ISENTROPIC STAGNATION PRESSURE WITH
     &IN THE BLADE = ',T67,F10.1,/,' COOLANT STAGNATION PRESSURE AT COOL
     &ING HOLE EXIT= ',T67,F10.1)  
      WRITE(6,*)
C
  150 CONTINUE
C
      DO 160 NCW = 1,NCOOLW
C
      IF(JCWS(NCW).LT.J1.AND.JCWE(NCW).LE.J1) GO TO 160
      IF(JCWS(NCW).GE.J2.AND.JCWE(NCW).GT.J2) GO TO 160
C
      SUM_PSTAT  = 0.0
      SUM_MAS    = 0.0
      DO 161 J = JCWS(NCW)+1,JCWE(NCW)

      IF(J.LE.J1.OR.J.GT.J2) GO TO 161

      DO 162 I=ICWS(NCW),ICWE(NCW)-1    
      IF(KC(NCW).EQ.1) THEN
           SUM_PSTAT  = SUM_PSTAT + CFLOWK1(I,J)*
     &     (P(I,J-1,1)+P(I,J,1)+P(I+1,J-1,1)+P(I+1,J,1))
           SUM_MAS = SUM_MAS + CFLOWK1(I,J)
      ELSE
           SUM_PSTAT  = SUM_PSTAT + CFLOWKM(I,J)*
     &     (P(I,J-1,KM)+P(I,J,KM)+P(I,J-1,KM)+P(IM,J,KM))
           SUM_MAS = SUM_MAS + CFLOWKM(I,J)
      ENDIF
  162 CONTINUE

  161 CONTINUE
C
C     "PS_EJECT"  IS THE AVERAGE STATIC PRESSURE ON THE SURFACE OF THE COOLING PATCH.
C     "PO_AT_HOLE"  IS THE ABSOLUTE STAGNATION PRESSURE AT WHICH THE COOLING 
C     FLOW WOULD NEED TO BE SUPPLIED TO THE BLADE ROW IF THERE WAS NO LOSS IN
C     THE INTERNAL COOLING PASSAGES.
C
      PS_EJECT  = 0.25*SUM_PSTAT/SUM_MAS
C      TS_EJECT  = TORELW(NCW)/(1.+.5*(GA-1)*MACHCOOL(NCW)*MACHCOOL(NCW))
C  JDD CHANGED TO USE STORED VALUE OF TS_EJECT. MAY 2018 .
      TS_EJECT  = TSTAT_EJECTW(NCW)
C
      IF(IFGAS.EQ.0) THEN
           PO_AT_HOLE = PS_EJECT*(TOCOOLW(NCW)/TS_EJECT)**RFGA
      ELSE
           TRAT       = TOCOOLW(NCW)/TS_EJECT
           PRAT       = PRAT_FROM_TRAT(TS_EJECT,TRAT,ALPHA,BETA1,BETA2)
           PO_AT_HOLE = PS_EJECT*PRAT
      END IF
C
      POCOOL     = POCOOLW(NCW)
      TOCOOL     = TOCOOLW(NCW)
      SUM_MAS    = SUM_MAS*NBLADE(JCWS(NCW)+1)
      CFLOW_STG  = CFLOW_STG + SUM_MAS
C
C      IF(IFCOOL.EQ.2) PO_AT_HOLE = POCOOL
C
      IF(IFGAS.EQ.0) THEN
           TO2_IS        = TOCOOL*(PO2/POCOOL)**FGA
           WCOOL_IS      = WCOOL_IS + SUM_MAS*CP*(TOCOOL - TO2_IS)
           TO2IS_AT_HOLE = TOCOOL*(PO2/PO_AT_HOLE)**FGA
           WCOOL_AT_HOLE = WCOOL_AT_HOLE 
     &                   + SUM_MAS*CP*(TOCOOL - TO2IS_AT_HOLE)
           ENTHIN        = ENTHIN   + SUM_MAS*CP*TOCOOL
      ELSE
           PRAT          = PO2/POCOOL
           TRAT          = TRAT_FROM_PRAT(TOCOOL,PRAT,FGAGAS,R_ALPHA,
     &                     BALPHA1,BALPHA2)
           TO2_IS        = TOCOOL*TRAT
           HOCOOL_IS     = HFROMT(TO2_IS,TREF,HREF,CP1,CP2,CP3)
           HOCOOL        = HFROMT(TOCOOL,TREF,HREF,CP1,CP2,CP3)
           WCOOL_IS      = WCOOL_IS + SUM_MAS*(HOCOOL - HOCOOL_IS)
           PRAT          = PO2/PO_AT_HOLE
           TRAT          = TRAT_FROM_PRAT(TOCOOL,PRAT,FGAGAS,R_ALPHA,
     &                     BALPHA1,BALPHA2)
           TO2IS_AT_HOLE = TOCOOL*TRAT
           HO2IS_AT_HOLE = HFROMT(TO2IS_AT_HOLE,TREF,HREF,CP1,CP2,CP3)
           WCOOL_AT_HOLE = WCOOL_AT_HOLE 
     &                  + SUM_MAS*(HOCOOL - HO2IS_AT_HOLE)
           ENTHIN       = ENTHIN + SUM_MAS*HOCOOL
      END IF
C
      WRITE(6,163) NCW,POCOOL,PO_EJECTW(NCW),PO_AT_HOLE
  163 FORMAT(' WALL COOLING PATCH NUMBER',I3,' SPECIFIED SUPPLY STAGNATI 
     &ON PRESSURE  = ',T67,F10.1,/,' ISENTROPIC STAGNATION PRESSURE WITH
     &IN THE BLADE = ',T67,F10.1,/,' COOLANT STAGNATION PRESSURE AT COOL
     &ING HOLE EXIT= ',T67,F10.1)  
      WRITE(6,*)
C
  160 CONTINUE
C
C
C******************************************************************************
C******************************************************************************
C
C     MOD 9/5/2002. TO ALLOW FOR ERRORS IN CONTINUITY IT IS BETTER TO BASE THE ACTUAL
C     WORK ON THE INLET FLOW, ADDED COOLING FLOW, AND THE ACTUAL TEMPERATURE CHANGE.
C     OTHERWISE SMALL CONTINUITY ERRORS CAUSE SIGNIFICANT EFFICIENCY ERRORS.
C
      FLOW_OUT =   BLADE_FLOW(J1) +  CFLOW_STG
C
      IF(IFGAS.EQ.0) THEN
           WNET     = CP*(BLADE_FLOW(J1)*TO_IN - FLOW_OUT*TO_OUT)
     &              + ENTHIN
           TO2_IS   = TO_IN*(PO2/POIN)**FGA
           WMAIN_IS = BLADE_FLOW(J1)*CP*(TO_IN - TO2_IS)
      ELSE
           HO_IN    = HFROMT(TO_IN,TREF,HREF,CP1,CP2,CP3)
           HO_OUT   = HFROMT(TO_OUT,TREF,HREF,CP1,CP2,CP3)
           PRAT     = PO2/POIN
           TRAT     = TRAT_FROM_PRAT(TO_IN,PRAT,FGAGAS,R_ALPHA,
     &                BALPHA1,BALPHA2)
           TO2_IS   = TO_IN*TRAT
           HO2_IS   = HFROMT(TO2_IS,TREF,HREF,CP1,CP2,CP3)
           WMAIN_IS = BLADE_FLOW(J1)*(HO_IN - HO2_IS)
           WNET     = BLADE_FLOW(J1)*HO_IN  + ENTHIN - FLOW_OUT*HO_OUT
      END IF
C
      WTOT_IS  =   WMAIN_IS + WCOOL_IS
      PCCOOL   =   100.*CFLOW_STG/BLADE_FLOW(J1)
C
      IF(WNET.GE.0.0) ETATT      = WNET/WTOT_IS
      IF(WNET.GE.0.0) ETATT_NOP  = (WNET + WPUMP(NSTG))/WTOT_IS
      IF(WNET.LT.0.0) ETATT      = WTOT_IS/WNET
      IF(WNET.LT.0.0) ETATT_NOP  = WTOT_IS/(WNET + WPUMP(NSTG))
C
      WTOT_AT_HOLE = WMAIN_IS + WCOOL_AT_HOLE
      IF(WNET.GE.0.0) ETATT_AT_HOLE =  WNET/WTOT_AT_HOLE
      IF(WNET.GE.0.0) ETA_ATH_NOP   =  (WNET + WPUMP(NSTG))/WTOT_AT_HOLE
      IF(WNET.LT.0.0) ETATT_AT_HOLE =  WTOT_AT_HOLE/WNET
      IF(WNET.LT.0.0) ETA_ATH_NOP   =  WTOT_AT_HOLE/(WNET + WPUMP(NSTG))
C
C     WRITE THE OUPUT FOR THE STAGE OR THE WHOLE MACHINE
C
      WRITE(6,*)
      IF(WNET.GT.0.0) THEN
           WRITE(6,*) ' THIS STAGE IS A TURBINE '
           WRITE(6,*)
           WRITE(6,*) ' PRESSURE RATIO, POinlet/POexit             = ',
     &                  PORAT 
      ELSE
           WRITE(6,*) ' THIS STAGE IS A COMPRESSOR '
           WRITE(6,*)
           WRITE(6,*) ' STAGNATION PRESSURE RATIO, POexit/POinlet  = ',
     &                  1.0/PORAT 
      END IF
C
C     WORK OUT THE GAS PROPERTIES AT STAGE INLET
C
      IF(IFGAS.NE.0) THEN
           DELT  = TO_IN - TREF
           CP = CP1 + CP2*DELT + CP3*DELT*DELT
           GA = CP/(CP - RGAS)
      END IF
C
C
      WRITE(6,*)      ' GAS PROPERTIES AT STAGE INLET, CP = ',CP,
     &                ' GAMMA = ',GA,' R = ',RGAS
      WRITE(6,*)      ' INLET AND EXIT STAGNATION PRESSURES        = ',
     & POIN, PO2
      WRITE(6,*)      ' INLET AND EXIT STATIC PRESSURES            = ',
     &                  PSIN, PSOUT
      WRITE(6,*)      ' INLET AND EXIT STAGNATION TEMPERATURES     = ',
     &                  TO_IN, TO_OUT
      WRITE(6,*)      ' INLET AND OUTLET MASS FLOW RATES           = ',
     &                  BLADE_FLOW(J1), FLOW_OUT, ' Kg/sec .'           
      WRITE(6,*)      ' NET POWER OUTPUT                           = ',
     &                  WNET/1000.,' KW. THIS INCLUDES ANY WORK DONE ON 
     &THE COOLANT BY THE ROTOR- i.e. THE PUMPING WORK.'
      WRITE(6,*)      ' NEGLECTING PUMPING POWER, NET POWER OUTPUT = ', 
     &                  (WNET + WPUMP(NSTG))/1000. ,' KW.'
      WRITE(6,*)      ' COOLANT PUMPING POWER                      = ',
     &                  WPUMP(NSTG)/1000., ' KW.'
      WRITE(6,*)      ' TOTAL COOLANT FLOW ADDED                   = ',
     &                  CFLOW_STG, ' Kg/sec.'
      WRITE(6,*)      ' PERCENTAGE COOLANT FLOW ADDED              = ',
     &                  PCCOOL
C
      IF(NSTG.LE.NSTAGES) THEN
           REACTION = (SUMTSTAT(JMIX_STAGE) - SUMTSTAT(J2M1))
     &               /(SUMTSTAT(J1)         - SUMTSTAT(J2M1))
           IF(WNET.LT.0.0) REACTION = 1.0 - REACTION
           WRITE(6,*) ' MEAN STAGE REACTION                        = ',
     &     REACTION
           WRITE(6,*)
      END IF
C
      WRITE(6,*) '*****************************************************'
      WRITE(6,*)
      WRITE(6,*) ' USING THE VALUES OF COOLANT STAGNATION PRESSURE INPUT
     & AS DATA'
      WRITE(6,*) '- THEN:-  '
      WRITE(6,*)
      WRITE(6,*) ' TOTAL-TO-TOTAL ISENTROPIC EFFICIENCY, ALLOWING FOR TH
     &E POTENTIAL WORK '
      WRITE(6,*) ' OF ALL COOLING FLOWS AND THE PUMPING POWER          =        
     &               ',  ETATT
      WRITE(6,*)
      WRITE(6,*) ' NEGLECTING THE PUMPING POWER, TOTAL-TO-TOTAL ISENTROP
     &IC EFFICIENCY= ', ETATT_NOP
      WRITE(6,*)
      WRITE(6,*) '*****************************************************'
C
      IF(IFCOOL.EQ.1) THEN
      WRITE(6,*) '*****************************************************'
      WRITE(6,*)
      WRITE(6,*) ' ASSUMING NO STAGNATION PRESSURE LOSS OR GAIN IN THE C
     &COOLANT SUPPLY PASSAGES'
      WRITE(6,*) ' AND USING THE STATIC PRESSURE AT EXIT FROM THE EJECTI
     &ON HOLES TO ESTIMATE'
      WRITE(6,*) ' THE COOLANT SUPPLY STAGNATION  PRESSURE'
      WRITE(6,*) ' -THEN:- '
      WRITE(6,*)
      WRITE(6,*) ' TOTAL-TO-TOTAL ISENTROPIC EFFICIENCY, ALLOWING FOR TH
     &E POTENTIAL WORK '
      WRITE(6,*) ' OF ALL COOLING FLOWS AND THE PUMPING POWER          =        
     &               ',  ETATT_AT_HOLE
      WRITE(6,*)
      WRITE(6,*) ' NEGLECTING THE PUMPING POWER, TOTAL-TO-TOTAL ISENTROP
     &IC EFFICIENCY= ', ETA_ATH_NOP
      WRITE(6,*)
      WRITE(6,*) '*****************************************************'
      END IF
C
 3333 CONTINUE
C
      END IF
C
C****************************************************************************
C****************************************************************************
C
C     NOW CALCULATE THE MACHINE PERFORMANCE IF THERE ARE NO COOLING FLOWS
C
      IF(IFCOOL.EQ.0) THEN
C
C******************************************************************************
C******************************************************************************
C
C      IF NO COOLING FLOWS WORK OUT THE  ISENTROPIC EFFICIENCY FOR EACH STAGE
C      THIS IS BASED ON THE MASS AVERAGE QUANTITIES AT INLET AND AT JMM1
C      FIRST FOR EVERY STAGE
C
      NSTGP1    = NSTAGES + 1
      IF(NSTAGES.EQ.1) NSTGP1 = 1
C
C     NSTGP1 IS USED TO CALCULATE THE PERFORMANCE OF THE WHOLE MACHINE
C
C*****************************************************************************************
C
      DO 2222  NSTG = 1,NSTGP1
C
      IF(NSTG.LE.NSTAGES) THEN
C
              J1      = JSTG_START(NSTG)
              J2      = JSTG_END(NSTG) + 1
              J2M1    = J2
C
C     JDD CHANGE  2/7/07  TO INCLUDE THE MIXING LOSS IN THE EFFICIENCY OF THE UPSTREAM BLADE ROW.
C     FOR THE LAST STAGE THERE IS NO MIXING DOWNSTREAM. BASE THE EFFY ON THE LAST BUT ONE POINT.
C
              IF(NSTG.EQ.NSTAGES) J2M1 = JM -1 
C
C     END OF JDD ADDITION
C
C     JMIX_STAGE IS ONLY USED TO FIND THE STAGE REACTION.
              N_ROW      = NROW(J1+1)
              JMIX_STAGE = JMIX(N_ROW) 
      ELSE
C
C   FOR THE WHOLE MACHINE THERE IS NO MIXING DOWNSTREAM. BASE THE EFFY ON THE LAST BUT ONE POINT, J2M1.
C
              J1   = 2
              J2   = JM
              J2M1 = JM - 1
              JMIX_STAGE = (J1+J2)/2
      END IF
C
      IF(NSTGP1.GT.1.AND.NSTG.EQ.NSTGP1) GO TO 123
C
      KMID = KM/2
      PHUB_IN = 0.0
      DO  I = 1,IMM1
      PHUB_IN  = PHUB_IN + 0.5*(P(I,J1,1)+P(I+1,J1,1))*FP(I)
      END DO
C
      PMID_IN = 0.0
      DO  I = 1,IMM1
      PMID_IN  = PMID_IN + 0.5*(P(I,J1,KMID)+P(I+1,J1,KMID))*FP(I)
      END DO
C
C
      PTIP_IN = 0.0
      DO  I = 1,IMM1
      PTIP_IN  = PTIP_IN + 0.5*(P(I,J1,KM)+P(I+1,J1,KM))*FP(I)
      END DO
C
C
      PHUB_MID = 0.0
      DO  I = 1,IMM1
      PHUB_MID  = PHUB_MID 
     &          + 0.5*(P(I,JMIX_STAGE,1)+P(I+1,JMIX_STAGE,1))*FP(I)
      END DO
C
      PMID_MID = 0.0
      DO  I = 1,IMM1
      PMID_MID  = PMID_MID
     &        + 0.5*(P(I,JMIX_STAGE,KMID)+P(I+1,JMIX_STAGE,KMID))*FP(I)
      END DO
C
      PTIP_MID = 0.0
      DO  I = 1,IMM1
      PTIP_MID  = PTIP_MID 
     &          + 0.5*(P(I,JMIX_STAGE,KM)+P(I+1,JMIX_STAGE,KM))*FP(I)
      END DO
C
C
      PHUB_OUT = 0.0
      DO  I = 1,IMM1
      PHUB_OUT  = PHUB_OUT 
     &          + 0.5*(P(I,J2M1,1)+P(I+1,J2M1,1))*FP(I)
      END DO
C
      PMID_OUT = 0.0
      DO  I = 1,IMM1
      PMID_OUT  = PMID_OUT
     &          + 0.5*(P(I,J2M1,KMID)+P(I+1,J2M1,KMID))*FP(I)
      END DO
C
      PTIP_OUT = 0.0
      DO  I = 1,IMM1
      PTIP_OUT  = PTIP_OUT 
     &          + 0.5*(P(I,J2M1,KM)+P(I+1,J2M1,KM))*FP(I)
      END DO
C
C     CALCULATE THE REACTIONS AS FOR A TURBINE
      REAC_HUB  = (PHUB_MID - PHUB_OUT)/(PHUB_IN - PHUB_OUT)
      REAC_MID  = (PMID_MID - PMID_OUT)/(PMID_IN - PMID_OUT)
      REAC_TIP  = (PTIP_MID - PTIP_OUT)/(PTIP_IN - PTIP_OUT)
C
C     CHANGE THE REACTIONS IF A COMPRESSOR
      IF(PMID_OUT.GT.PMID_IN) THEN
           REAC_HUB  = 1.0 - REAC_HUB
           REAC_MID  = 1.0 - REAC_MID
           REAC_TIP  = 1.0 - REAC_TIP
      END IF
C
  123 CONTINUE
C 
      POIN    = SUMPO(J1)*PO1(KMID)
      PO2     = SUMPO(J2M1)*PO1(KMID)
      PSIN    = SUMPSTAT(J1)
      PSOUT   = SUMPSTAT(J2M1)
      TO_IN   = SUMTO(J1)
      TO_OUT  = SUMTO(J2M1)
C
      IF(IFGAS.EQ.0) THEN
           TORAT   = TO_OUT/TO_IN
           DELTO   = TO_OUT - TO_IN
           PORAT   = PO2/POIN
           PSRAT   = PSOUT/POIN
           IF(PORAT.LT.0.01)  PORAT   = 0.01
           TO2IS     = PORAT**FGA *TO_IN
           TSTAT2_IS = PSRAT**FGA *TO_IN
           DELTOIS = TO2IS  - TO_IN
           DELT_TOT_STAT = TSTAT2_IS - TO_IN
           ETAIS   = 1.0
           IF(DELTO.GT.0.0) ETAIS = DELTOIS/DELTO
           IF(DELTO.LT.0.0) ETAIS = DELTO/DELTOIS
           EPOLY   = ALOG(TORAT)/ALOG(PORAT)/FGA
           IF(TORAT.GT.1.0) EPOLY = 1.0/EPOLY
C
           ETATS   = 1.0
           IF(DELTO.GT.0.0) ETATS = DELT_TOT_STAT/DELTO
           IF(DELTO.LT.0.0) ETATS = DELTO/DELT_TOT_STAT
C
      ELSE
           PORAT   = PO2/POIN
           TRAT   = TRAT_FROM_PRAT(TO_IN,PORAT,FGAGAS,R_ALPHA,
     &              BALPHA1,BALPHA2)
           TO_OUT_IS = TO_IN*TRAT
           HO_OUT_IS = HFROMT(TO_OUT_IS,TREF,HREF,CP1,CP2,CP3)
           HO_IN     = HFROMT(TO_IN,TREF,HREF,CP1,CP2,CP3)
           HO_OUT    = HFROMT(TO_OUT,TREF,HREF,CP1,CP2,CP3)
           DELHO     = HO_OUT    - HO_IN
           DELHO_IS  = HO_OUT_IS - HO_IN
           ETAIS     = 1.0
           IF(DELHO.GT.0.0) ETAIS = DELHO_IS/DELHO
           IF(DELHO.LT.0.0) ETAIS = DELHO/DELHO_IS
           EPOLY     = 1.0
           DELT  = TO_IN - TREF
           CP = CP1 + CP2*DELT + CP3*DELT*DELT
           GA = CP/(CP - RGAS)
C
      END IF
C
      WRITE(6,*)'*******************************************************
     &*****************************************************************'
      IF(NSTG.LE.NSTAGES) THEN
           WRITE(6,*) ' RESULTS FOR STAGE NUMBER ',NSTG 
      ELSE
           WRITE(6,*) ' THE FOLLOWING RESULTS ARE FOR THE WHOLE MACHINE'
      END IF
C
      WRITE(6,*)'*******************************************************
     &*****************************************************************'
C
      WRITE(6,*)
      IF(DELTO.LT.0.0) THEN
           WRITE(6,*) ' ********** THIS IS A TURBINE ********** '
           WRITE(6,*)
           WRITE(6,*) ' PRESSURE RATIO, POinlet/POexit            = ',
     &                  1.0/PORAT 
      ELSE
           WRITE(6,*) ' ***********THIS IS A COMPRESSOR ********** '
           WRITE(6,*)
           WRITE(6,*) ' PRESSURE RATIO, POexit/POinlet            = ',
     &                  PORAT 
      END IF
C
      WRITE(6,*)      ' GAS PROPERTIES AT STAGE INLET, CP = ',CP,
     &                ' GAMMA = ',GA,' R = ',RGAS
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,*) ' EFFICIENCIES BASED ON VALUES AT J = ', J1,J2M1
      WRITE(6,*) ' JMIX  = ', JMIX_STAGE
      WRITE(6,*)
      WRITE(6,*)      ' INLET AND EXIT STAGNATION PRESSURES   = ',
     &                  POIN, PO2
      WRITE(6,*)      ' INLET AND EXIT STATIC PRESSURES       = ',
     &                   PSIN,PSOUT
      WRITE(6,*)      ' INLET AND EXIT STAGNATION TEMPERATURES= ',
     &                  TO_IN, TO_OUT
      WRITE(6,*)      ' INLET AND OUTLET MASS FLOW RATES      = ',
     &                  BLADE_FLOW(J1),BLADE_FLOW(J2M1)
C
      WRITE(6,*)
      WRITE(6,*) ' USING THE MASS AVERAGED STAGNATION PRESSURES AND '
      WRITE(6,*) ' TEMPERATURES AT INLET AND JM-1 .'
      WRITE(6,*)
      WRITE(6,*) 
     & ' TOTAL TO TOTAL ISENTROPIC EFFICIENCY       = ', ETAIS 
      IF(IFGAS.EQ.0)  WRITE(6,*) 
     & ' TOTAL TO STATIC ISENTROPIC EFFICIENCY      = ', ETATS 
      IF(IFGAS.EQ.0)  WRITE(6,*)
     & ' TOTAL TO TOTAL POLYTROPIC EFFICIENCY       = ', EPOLY
      WRITE(6,*)
C
      IF(NSTG.LE.NSTAGES) THEN
C
      WRITE(6,*)      ' INLET, MID and EXIT STATIC TEMPERATURES   = ',
     &                 SUMTSTAT(J1),SUMTSTAT(JMIX_STAGE),SUMTSTAT(J2M1)
C
                       REACTION = (SUMTSTAT(JMIX_STAGE)- SUMTSTAT(J2M1))
     &                 /(SUMTSTAT(J1)         - SUMTSTAT(J2M1))
                       IF(PSOUT.GT.PSIN) REACTION = 1.0 - REACTION
C     
      WRITE(6,*)  
      WRITE(6,*) 'AVERAGE STAGE REACTION BASED ON MASS AVERAGED TEMPERAT
     &URES = ',   REACTION
      WRITE(6,*)
C
      WRITE(6,*)
      WRITE(6,*) 'REACTIONS FOR STAGE NUMBER ', NSTG
      WRITE(6,*) 'HUB REACTION BASED ON PRESSURE CHANGES     ',REAC_HUB
      WRITE(6,*) 'MID SPAN REACTION BASED ON PRESSURE CHANGES',REAC_MID
      WRITE(6,*) 'TIP REACTION BASED ON PRESSURE CHANGES     ',REAC_TIP
      WRITE(6,*)
C
      IF(IFGAS.EQ.0) THEN
           WORK = (CP*(BLADE_FLOW(J1)*TO_IN
     &          - BLADE_FLOW(J2M1)*TO_OUT))/1000.
      ELSE
           WORK = (BLADE_FLOW(J1)*HO_IN
     &          -  BLADE_FLOW(J2M1)*HO_OUT)/1000.
      END IF
C
      WRITE(6,*)
      WRITE(6,*)  ' STAGE WORK NEGLECTING ANY BLEED OR COOLING FLOWS = '
     &     , WORK,' KILOWATTS.'
      WRITE(6,*)
C
C
      ELSE
C
C********************************************************************************
C
C     NOW FOR THE WHOLE MACHINE, NSTG = NSTAGES + 1
C
      IF(IFBLEED.EQ.0) HBLEDTOT = 0.0
C
      IF(IFGAS.EQ.0) THEN
           WORK = (CP*(BLADE_FLOW(1)*TO_IN-BLADE_FLOW(JMM1)*TO_OUT)
     &          +  ENTHIN  - HBLEDTOT  )/1000.
      ELSE
           WORK = (BLADE_FLOW(1)*HO_IN  - BLADE_FLOW(JMM1)*HO_OUT 
     &          + ENTHIN  - HBLEDTOT)/1000.
      END IF
C   
      WRITE(6,*)
      WRITE(6,*) ' INLET ENTHALPY FLUX,  KW  = ',
     &             CP*BLADE_FLOW(1)*TO_IN/1000.
      WRITE(6,*) ' OUTLET ENTHALPY FLUX, KW  = ',
     &             CP*BLADE_FLOW(JMM1)*TO_OUT/1000.
      WRITE(6,*) ' ENTHALPY BLED OFF,    KW  = ', HBLEDTOT/1000.
C
      WRITE(6,*)
      WRITE(6,*) ' FOR THE WHOLE MACHINE '
      WRITE(6,*) ' OVERALL POWER OUTPUT             = ',WORK,
     &                ' KILOWATTS.'
      WRITE(6,*)
      WRITE(6,*) ' INLET AND OUTLET MASS FLOW RATES = ',
     &             BLADE_FLOW(1),BLADE_FLOW(JMM1),' Kg/sec.'   
      WRITE(6,*)
C        
      IF(IFBLEED.NE.0) THEN
           WRITE(6,*)
           WRITE(6,*) ' TOTAL MASS FLOW BLED OFF = ',SUMBLEED(JM),
     &                ' Kg/sec'
           WRITE(6,*)
           WRITE(6,*) ' THE POWER OUTPUT ALLOWS FOR THE BLEED FLOWS AND 
     &COOLING FLOWS.'
           WRITE(6,*) ' BUT THE EFFICIENCIES ONLY RELATE THE MASS AVERAG
     & ED STATES OF THE FLUID AT MACHINE ENTRY AND EXIT.'
           WRITE(6,*)
      END IF
C
      END IF
C
 2222 CONTINUE
C
      END IF
C
C
      RETURN
      END
C
C********************************************************************************
C********************************************************************************
C
      SUBROUTINE SMOOTH(N1,N2,NSMOOTH,FSMOOTH,FRAC,VAR)
C
      INCLUDE 'commall-open-18.3'
C 
      DIMENSION FRAC(JD),VAR(JD),TEMP(JD)
C
      DO 10 ITS = 1,NSMOOTH
C
      DO N = N1,N2
      TEMP(N) = VAR(N)
      END DO
C
      DO N = N1+1, N2-1
      FLEFT  = FRAC(N)   - FRAC(N-1)
      FRIGHT = FRAC(N+1) - FRAC(N)
      AVG = ( FLEFT*TEMP(N+1) + FRIGHT*TEMP(N-1) )/(FLEFT + FRIGHT)
      VAR(N) = (1.0-FSMOOTH)*VAR(N) + FSMOOTH*AVG
      END DO
C
   10 CONTINUE
C
      RETURN
      END
C******************************************************************************
C
      SUBROUTINE SET_COEFFS
C
C      THIS SUBROUTINE SETS COEFFICIENTS NEEDED FOR THE NON-PERFECT GAS PROPERTY FUNCTIONS    
C
      INCLUDE 'commall-open-18.3'
C
      HREF = CP1*TREF
      CV1  = CP1 - RGAS
      EREF = CV1*TREF
C
      WRITE(6,*) 'TREF, HREF,CV1, EREF',TREF, HREF,CV1, EREF
C
      HT1 = 1.0/CP1   
      HT2 = -0.5*CP2/CP1**3
      HT3 = (3.0*CP2**2 - 2.0*CP3*CP1)/6/CP1**5
      HT4 = (-15*CP2**3/CP1**7 + 20*CP2*CP3/CP1**6  )/24
      WRITE(6,*) 'HT1,HT2,HT3,HT4',HT1,HT2,HT3,HT4
C
      ET1 = 1.0/CV1
      ET2 = -0.5*CP2/CV1**3
      ET3 = (3.0*CP2**2 - 2.0*CP3*CV1)/6/CV1**5
      ET4 = (-15*CP2**3/CV1**7 + 20*CP2*CP3/CV1**6  )/24
      WRITE(6,*) 'ET1,ET2,ET3,ET4',ET1,ET2,ET3,ET4
C
      ALPHA   = (CP1 - TREF*CP2 + TREF**2*CP3)/RGAS
      R_ALPHA = 1.0/ALPHA
      BETA1   = (CP2 -2.*CP3*TREF)/RGAS
      BETA2   = 0.5*CP3/RGAS
      BALPHA1 = -BETA1/ALPHA
      BALPHA2 = -BETA2/ALPHA 
C
      EPS0    = 1.0/(ALPHA - 1.0)
      EPS1    = -BETA1*EPS0
      EPS2    = -BETA2*EPS0
      WRITE(6,*) ' EPS0,EPS1,EPS2 ',EPS0,EPS1,EPS2
C
      GAGAS   = CP1/CV1
      FGAGAS  = (GAGAS-1.0)/GAGAS
      GA1GAS  = GAGAS - 1.0
      WRITE(6,*) ' GAGAS,FGAGAS,GA1GAS ',GAGAS,FGAGAS,GA1GAS
C
      RETURN
      END
C
C******************************************************************************
C
      FUNCTION HFROMT(TIN,TREF,HREF,CP1,CP2,CP3)
C
      HFROMT = CP1*TIN + CP2*(TIN-TREF)*(TIN-TREF)/2 + 
     &         CP3*(TIN-TREF)*(TIN-TREF)*(TIN-TREF)/3
C
      RETURN
      END  
C
C******************************************************************************
C
      FUNCTION TFROMH(HIN,TREF,HREF,HT1,HT2,HT3,HT4)
C
      DH     = HIN - HREF
      DH2    = DH*DH
      TFROMH = TREF + HT1*DH + HT2*DH2 + HT3*DH*DH2 + HT4*DH2*DH2
C
      RETURN
      END     
C
C******************************************************************************
C
      FUNCTION TFROME(EIN,TREF,EREF,ET1,ET2,ET3,ET4)
C
      DE     = EIN - EREF
      DE2    = DE*DE
      TFROME = TREF + ET1*DE + ET2*DE2 + ET3*DE*DE2 + ET4*DE2*DE2
C
      RETURN
      END  
C
C******************************************************************************
C
      FUNCTION PRAT_FROM_TRAT(TSTART,TRAT,ALPHA,BETA1,BETA2)
C
      T2 = TSTART*TRAT
      BETA  = BETA1*(T2-TSTART) + BETA2*(T2*T2 - TSTART*TSTART)
      PRAT_FROM_TRAT  = TRAT**ALPHA   *   EXP(BETA)
C
      RETURN
      END  
C
C******************************************************************************
C
C
      FUNCTION TRAT_FROM_PRAT(TSTART,PRAT,FGAGAS,R_ALPHA,
     &         BALPHA1,BALPHA2)
C
      PALPHA = PRAT**R_ALPHA
C
      TRAT   = PRAT**FGAGAS
      T2     = TSTART*TRAT
      BETA   = BALPHA1*(T2-TSTART) + BALPHA2*(T2*T2-TSTART*TSTART)
      TRAT   = PALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
C
      BETA   = BALPHA1*(T2-TSTART) + BALPHA2*(T2*T2-TSTART*TSTART)
      TRAT   = PALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
C
      BETA   = BALPHA1*(T2-TSTART) + BALPHA2*(T2*T2-TSTART*TSTART)
      TRAT   = PALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
C
      BETA   = BALPHA1*(T2-TSTART) + BALPHA2*(T2*T2-TSTART*TSTART)
      TRAT   = PALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
C
      TRAT_FROM_PRAT     =  T2/TSTART
      RETURN
      END
C
C******************************************************************************
C
C
      FUNCTION TRAT_FROM_RORAT(TSTART,RORAT,GA1GAS,EPS0,
     &         EPS1,EPS2)
C
      TALPHA = RORAT**EPS0
      TRAT   = RORAT**GA1GAS
      T2     = TSTART*TRAT
      BETA   = EPS1*(T2-TSTART) + EPS2*(T2*T2-TSTART*TSTART)
      TRAT   = TALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
      BETA   = EPS1*(T2-TSTART) + EPS2*(T2*T2-TSTART*TSTART)
      TRAT   = TALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
      BETA   = EPS1*(T2-TSTART) + EPS2*(T2*T2-TSTART*TSTART)
      TRAT   = TALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
      BETA   = EPS1*(T2-TSTART) + EPS2*(T2*T2-TSTART*TSTART)
      TRAT   = TALPHA  *  EXP(BETA)
      T2     = 0.5*(T2 + TSTART*TRAT)
      TRAT_FROM_RORAT     =  T2/TSTART
      RETURN
      END
C
C**********************************************************************
C
C      
      SUBROUTINE MASS_AVG(IM,JM,KM,FCN,FLOWX,SUMFUN,AVGALL,ID,JD,KD)
C
      DIMENSION FCN(ID,JD,KD),FLOWX(ID,JD,KD), SUMFUN(KD)
C
      J = JM
      SUMALL   = 0.0
      FLOWALL  = 0.0
C
      DO 110 K=1,KM-1
      SUMAS    = 0.0
      SUMFUN(K)= 0.0
      DO 120 I=1,IM-1
      DFLOW = -FLOWX(I,J,K)
      IF(DFLOW.LT.0.0) DFLOW = 0.0
      SUMAS     = SUMAS+DFLOW
      SUMFUN(K) = SUMFUN(K) + DFLOW*(FCN(I,J,K)+FCN(I+1,J,K)
     &          + FCN(I+1,J,K+1)+FCN(I,J,K+1))*0.25  
  120 CONTINUE
      SUMALL    = SUMALL  + SUMFUN(K)
      FLOWALL   = FLOWALL + SUMAS
      SUMFUN(K) = SUMFUN(K)/SUMAS
  110 CONTINUE
C
      DO 200 K=2,KM-1
      SUMFUN(K) = 0.5*(SUMFUN(K) + SUMFUN(K-1))
  200 CONTINUE
C
      SUMFUN(1)  = 2.0*SUMFUN(2)    - SUMFUN(3)
      SUMFUN(KM) = 2.0*SUMFUN(KM-1) - SUMFUN(KM-2)
      AVGALL     = SUMALL/FLOWALL
C
      RETURN
      END
C
C************************************************************************************
C************************************************************************************
C
      SUBROUTINE RESTAGG(NSEC,J1,J2,JLEROW,JTEROW,ROTATE,FRACX_ROT)
C
C************************************************************************************
C     THIS SUBROUTINE RESTAGGERS A BLADE SECTION.
C     IT ASSUMES NEGLIGIBLE CHANGE OF RADIUS OF THE STREAM SURFACEAND SO IS NOT USABLE
C     FOR RADIAL FLOW MACHINES.
C
C     IT IS NOT USED IN  VERSION 16.3  AND ABOVE. USE SUBROUTINE "RESTAGGER" INSTEAD .
C************************************************************************************
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  XIN(JD),YUP(JD),YLOW(JD),YMID(JD),
     &           XUP(JD),XMID(JD),XLOW(JD)
C
      ROTATE     = 0.0
      FRACX_ROT  = 0.5
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=1502) ROTATE, FRACX_ROT
 1502 CONTINUE
      WRITE(6,*) ' INPUT FOR RESTAGGERING OPTION'
      WRITE(6,*) ' ROTATION =', ROTATE, 'ABOUT ',FRACX_ROT 
      WRITE(6,*)
C
      RAD_ROT = ROTATE*DEGRAD
C
C     FIND THE Y COORDINATES OF THE UPPER AND LOWER BLADE SURFACES.
C
      NJ  = J2 - J1 + 1
      DO 10 J  = J1,J2
      JLOC  = J - J1 + 1
      XIN(JLOC)   = XSURF(J,NSEC)
      YUP(JLOC)   = RT_UPP(J,NSEC)
      YLOW(JLOC)  = YUP(JLOC) - RT_THICK(J,NSEC)
      YMID(JLOC)  = 0.5*(YUP(JLOC) + YLOW(JLOC))      
   10 CONTINUE
C
C     FIND THE COORDINATES OF THE AXIS OF ROTATION
C
      XROT = XIN(JLEROW) + FRACX_ROT*(XIN(JTEROW) - XIN(JLEROW))
      CALL INTP(NJ,XIN,YMID,XROT,YROT)
C
C     ROTATE ALL THE POINTS CLOCKWISE BY RAD_ROT RADIANS.
C
      DO 20 J = 1,NJ
C
      XREL      = XIN(J) - XROT
      YUPREL    = YUP(J) - YROT
      YMIDREL   = YMID(J)- YROT
      YLOWREL   = YLOW(J)- YROT
      ANGLUP    = ATAN2(YUPREL,XREL) - RAD_ROT
      ANGLMID   = ATAN2(YMIDREL,XREL)- RAD_ROT
      ANGLLOW   = ATAN2(YLOWREL,XREL)- RAD_ROT
      RADUP     = SQRT(XREL*XREL + YUPREL*YUPREL)
      RADMID    = SQRT(XREL*XREL + YMIDREL*YMIDREL)
      RADLOW    = SQRT(XREL*XREL + YLOWREL*YLOWREL)
      XUP(J)    = RADUP*COS(ANGLUP)   + XROT
      XMID(J)   = RADMID*COS(ANGLMID) + XROT
      XLOW(J)   = RADLOW*COS(ANGLLOW) + XROT
      YUP(J)    = RADUP*SIN(ANGLUP)   + YROT
      YMID(J)   = RADMID*SIN(ANGLMID) + YROT
      YLOW(J)   = RADLOW*SIN(ANGLLOW) + YROT
C      
C      WRITE(6,100) J,XMID(J),YUP(J),YMID(J),YLOW(J)
   20 CONTINUE

C
C    INTERPOLATE TO FIND NEW  Y VALUES ON THE BLADE SURFACES AT THE XMID VALUES OF X.
C    SCALE THE MOVEMENT BY FACMOVE SO THAT IT IS REDUCED UPSTREAM AND DOWNSTREAM OF THE
C    BLADE SUCH THAT FIRST AND LAST GRID POINTS DO NOT MOVE UNLESS THEY ARE THE LE AND TE POINTS.
C
      WRITE(6,*)
      WRITE(6,*) ' BLADE RESTAGGERED, NEW VALUES OF X,YUP,YTHICK,FMOVE=' 
      DO 30 J = 1,NJ
      JALL = J1 + J - 1
      FACMOVE = 1.0
      IF(J.LT.JLEROW) FACMOVE=(XMID(J)-XMID(1))/(XMID(JLEROW)-XMID(1))
      IF(J.GT.JTEROW) FACMOVE=(XMID(NJ)-XMID(J))/(XMID(NJ)-XMID(JTEROW)) 
      FM1 = 1.0 - FACMOVE
      XARG = XMID(J)
      CALL INTP(NJ,XUP,YUP,XARG,YUPNEW)
      CALL INTP(NJ,XLOW,YLOW,XARG,YLOWNEW)
      XSURF(JALL,NSEC)    = FM1*XSURF(JALL,NSEC)  + FACMOVE*XARG
      RT_UPP(JALL,NSEC)   = YUPNEW
      RT_THICK(JALL,NSEC) = (YUPNEW - YLOWNEW)
      WRITE(6,100) JALL,XSURF(JALL,NSEC),RT_UPP(JALL,NSEC),
     &             RT_THICK(JALL,NSEC),FACMOVE
   30 CONTINUE
C
  100 FORMAT('J = ',I5,4F15.7)
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE LEAN(K,J1,J2,ANGLEAN)
C
      INCLUDE 'commall-open-18.3'
C******************************************************************************
C******************************************************************************
C     LEAN THE BLADE BY ANGLEAN IF IF_LEAN IS GREATER THAN ZERO.
C     IF "ANGLEAN" IS POSITIVE THE HUB IS HELD FIXED AND THE OTHER SECTIONS ARE 
C     MOVED IN THE POSITIVE THETA DIRECTION .
C
      READ(5,*) DUMMY_INPUT
C
      ANGLEAN = 0.0
      READ(5,*,ERR=20) ANGLEAN
   20 CONTINUE
C
      WRITE(6,*) ' LEANING THE BLADE SECTION BY ANGLEAN = ',
     &             ANGLEAN,' DEGREES RELATIVE TO THE HUB SECTION.'
C
      ANGLEAN = ANGLEAN*DEGRAD
C
      DO 10 J= J1,J2
      JM1 = J-1
      JP1 = J+1
      IF(J.EQ.J1) JM1 = 1
      IF(J.EQ.J2) JP1 = J2
      XDIF = XSURF(JP1,1) - XSURF(JM1,1)
      RDIF = RSURF(JM1,1) - RSURF(JP1,1)
      SDIF = SQRT(XDIF*XDIF + RDIF*RDIF)
      RNORM   =   XDIF/SDIF
      XNORM   =  -RDIF/SDIF
      QDISTR  =   RSURF(J,K) - RSURF(J,1)
      QDISTX  =   XSURF(J,K) - XSURF(J,1)
      QDISTN  =   QDISTR*RNORM + QDISTX*XNORM         
      T_SHIFT =   ANGLEAN*QDISTN
      RT_UPP(J,K)   = T_SHIFT + RT_UPP(J,K)
   10 CONTINUE
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE MIX_BCONDS(IFOUT)
C
C***********************************************************************
C     WRITE THE MASS AVERAGE VALUES AT EXIT TO A FILE 'outbconds'
C***********************************************************************
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION FAVG(KD)
C
      OPEN(UNIT = 12,FILE ='mixbconds')
C
C     CALCULATE AND WRITE OUT THE EXIT STAGNATION PRESSURE.
C
      DO 7777 NR = 1, NROWS + 1
C
      WRITE(12,*)
      WRITE(12,*) ' BLADE ROW NUMBER ', NR, ' JMIX = ',JMIX(NR) 
      WRITE(12,*)
C
      IF(NR.EQ.1) THEN 
C
      J = 1
      IF(IFOUT.EQ.1) THEN
               WRITE(12,*)
               WRITE(12,*) ' INLET BOUNDARY CONDITIONS TO ROW NUMBER 1.'
               WRITE(12,*)
      END IF
C
      ELSE
C
      J = JMIX(NR-1)
      IF(IFOUT.EQ.1) THEN
      WRITE(12,*)  ' ROW NUMBER ', NR-1,' JMIX  =  ', J
      WRITE(12,*) ' EXIT FLOW CONDITIONS FROM THIS BLADE ROW'
      WRITE(12,*) ' THESE CAN BE USED AS THE INLET BOUNDARY CONDITIONS 
     &FOR THE NEXT STAGE.'
      WRITE(12,*)
      END IF
C
      END IF
C
C
      DO 1110 I=1,IM
      DO 1110 K=1,KM
      D=.5*(VX(I,J,K)*VX(I,J,K)+VT(I,J,K)*VT(I,J,K)+VR(I,J,K)*VR(I,J,K))
      TSTAG     = HO(I,J,K)/CP
      TS        = TSTAG - D/CP
      IF(TS.LT.1.) TS = 1.
      STORE(I,J,K) = P(I,J,K)*(TSTAG/TS)**(GA/(GA-1))
 1110 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,PO_OUT_AVG,POAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
             WRITE(12,*) ' STAGNATION PRESSURE '
             WRITE(12,1111) (PO_OUT_AVG(K),K=1,KM)
      END IF
 1111 FORMAT(8F10.1)
C
C     CALCULATE AND WRITE OUT THE EXIT STATIC PRESSURE
C
      CALL MASS_AVG(IM,J,KM,P,FLOWX,FAVG,PSAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' STATIC PRESSURE '
            WRITE(12,1111)(FAVG(K),K=1,KM)
      END IF
C
C
C     CALCULATE AND WRITE OUT THE EXIT STAGNATION TEMP.
C
      DO 1014 I=1,IM
      DO 1014 K=1,KM
      STORE(I,J,K) =  HO(I,J,K)/CP
 1014 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,TO_OUT_AVG,TOAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' STAGNATION TEMPERATURE '
            WRITE(12,1111) (TO_OUT_AVG(K),K=1,KM)
      END IF
 1112 FORMAT(8F10.4)
C
C
C     CALCULATE AND WRITE OUT THE EXIT STATIC TEMPERATURE.
C
      DO 1119 I = 1, IM
      DO 1119 K = 1, KM
      D=.5*(VX(I,J,K)*VX(I,J,K)+VT(I,J,K)*VT(I,J,K)+VR(I,J,K)*VR(I,J,K))
      STORE(I,J,K) = (HO(I,J,K) - D)/CP
 1119 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,FAVG,TSAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' STATIC TEMPERATURE '
            WRITE(12,1111) (FAVG(K),K=1,KM)
      END IF
C
C     WRITE OUT THE EXIT TANGENTIAL VELOCITY.
C
      CALL MASS_AVG(IM,J,KM,VT,FLOWX,VT_OUT_AVG,VTAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' ABSOLUTE TANGENTIAL VELOCITY '
            WRITE(12,1112) (VT_OUT_AVG(K),K=1,KM)
      END IF
C
C     CALCULATE AND WRITE OUT THE EXIT MERIDIONAL VELOCITY
C
      DO 1120 I=1,IM
      DO 1120 K=1,KM
      STORE(I,J,K) = SQRT(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K))
 1120 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,FAVG,VMAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' MERIDIONAL VELOCITY '
            WRITE(12,1112) (FAVG(K),K=1,KM)
      END IF
C
C    CALCULATE AND WRITE OUT THE YAW ANGLE AT EXIT
C
      DO 1130 I=1,IM
      DO 1130 K=1,KM
      VMER = SQRT(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K))
      STORE(I,J,K) = ATAN2(VT(I,J,K),VMER)*RADDEG
 1130 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,YAW_OUT_AVG,YAWAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' YAW ANGLE '
            WRITE(12,1112) (YAW_OUT_AVG(K),K=1,KM)
      END IF 
C
C    CALCULATE AND WRITE OUT THE PITCH ANGLE AT EXIT
C
      DO 1140 I=1,IM
      DO 1140 K=1,KM
      VMER = SQRT(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K))
      STORE(I,J,K) = ATAN2(VR(I,J,K),VMER)*RADDEG
 1140 CONTINUE
C
      CALL MASS_AVG(IM,J,KM,STORE,FLOWX,PITCH_OUT_AVG,
     &              PITCHAVGJM,ID,JD,KD)
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' PITCH ANGLE '
            WRITE(12,1112) (PITCH_OUT_AVG(K),K=1,KM)
      END IF
C
C      WRITE OUT   FR  and   FP  .
C
      IF(IFOUT.EQ.1) THEN
            WRITE(12,*) ' SPANWISE GRID SPACING '
            WRITE(12,1113) (FR(K),K=1,KM-1)
            WRITE(12,*) ' PITCHWISE GRID SPACING '
            WRITE(12,1113) (FP(I),I=1,IM-1)
            FSPAN(1) = 0.0
            DO K=2,KM
            FSPAN(K) = FSPAN(K-1) + FR(K-1)
            END DO
            WRITE(12,*) ' SPANWISE GRID POSITIONS'
            WRITE(12,1113) (FSPAN(K),K=1,KM)            
      END IF
C
 1113 FORMAT(8F10.6)
C
C     END OF ONE BLADE ROW, RETURN FOR NEXT ROW
 7777 CONTINUE
C
C      END OF WRITING FILE TO UNIT 12
C
      CLOSE(12)
C  
C*********************************************************************
C*********************************************************************
      RETURN
      END

C
C******************************************************************************
C
      SUBROUTINE NEWBCONDS 
C
      INCLUDE 'commall-open-18.3'
C
C     THIS SUBROUTINE CHANGES THE INLET BOUNDARY CONDITIONS TO SIMULATE A REPEATING STAGE
C
      RFINBC1  =  1.0 - RFINBC
C
      DO 10 K = 1,KM
      PO1(K)  = RFINBC1*PO1(K) 
     &        + RFINBC*(PO_IN_MID  + PO_OUT_AVG(K)   - PO_OUT_AVG(KMID))
C      TO1(K)  = RFINBC1*TO1(K)
C     &        + RFINBC*(TO_IN_MID  + TO_OUT_AVG(K)   - TO_OUT_AVG(KMID))
C      BS(K)   = RFINBC1*BS(K)
C     &        + RFINBC*(YAW_IN_MID + YAW_OUT_AVG(K) - YAW_OUT_AVG(KMID))
C
       BS(K)   = RFINBC1*BS(K) + RFINBC*YAW_OUT_AVG(K)
C
C      BR(K)   = RFINBC1*BR(K)
C     &      + RFINBC*(PITCH_IN_MID+PITCH_OUT_AVG(K)-PITCH_OUT_AVG(KMID))
      VTIN(K) = RFINBC1*VTIN(K)
     &        + RFINBC*(VT_IN_MID  + VT_OUT_AVG(K)   - VT_OUT_AVG(KMID))
      RPO1(K) = 1./PO1(K)
      BSSIN(K)= SIN(BS(K)*DEGRAD)
      BSCOS(K)= COS(BS(K)*DEGRAD)
      BRSIN(K)= SIN(BR(K)*DEGRAD)
      BRCOS(K)= COS(BR(K)*DEGRAD)
   10 CONTINUE
C
      RETURN
      END
C*****************************************************************************      
C
      SUBROUTINE CELL_TO_NODE(VCELL,VNODE)
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  VCELL(ID,JD,KD), VNODE(ID,JD,KD) 
C    
C     DISTRIBUTE THE CELL CENTRED VARIABLE VCELL 
C     TO THE NODE STORED VARIABLE VNODE
C*****************************************************************************
C
C     FIRST AVERAGE ALONG THE I DIRECTION
C
      IF(IM.GT.2) THEN

      DO 1100 K = 1,KMM1
      DO 1100 J = 1,JMM1
      DO 1101 I = 2,IMM1
      TEMP1(I,J,K)  = (VCELL(I-1,J,K) + VCELL(I,J,K))
 1101 CONTINUE
      TEMP1(1,J,K)  = 3.0*VCELL(1,J,K)    - VCELL(2,J,K)
      TEMP1(IM,J,K) = 3.0*VCELL(IMM1,J,K) - VCELL(IMM2,J,K)
 1100 CONTINUE

      ELSE

      DO 1102 K=1,KMM1
      DO 1102 J=1,JMM1
      TEMP1(1,J,K) = VCELL(1,J,K)
      TEMP1(2,J,K) = VCELL(1,J,K)
 1102 CONTINUE

      END IF
C
C     NEXT AVERAGE ALONG THE K DIRECTION
C  Q3D
C
      IF(KM.GT.2) THEN
C
      DO 1120 J = 1,JMM1
      DO 1120 I = 1,IM
      DO 1121 K = 2,KMM1
      TEMP2(I,J,K) = (TEMP1(I,J,K-1) + TEMP1(I,J,K))
 1121 CONTINUE
      TEMP2(I,J,1)   = 3.0*TEMP1(I,J,1)    - TEMP1(I,J,2)
      TEMP2(I,J,KM)  = 3.0*TEMP1(I,J,KMM1) - TEMP1(I,J,KMM2)
 1120 CONTINUE
C
      ELSE
C
      DO 1122 J = 1,JMM1
      DO 1122 I = 1,IM
      DO 1122 K=1,2
      TEMP2(I,J,K) = TEMP1(I,J,1) + TEMP1(I,J,2)
 1122 CONTINUE
C
      END IF
C
C   END Q3D
C
C      NOW AVERAGE ALONG THE J DIRECTION
C
      DO 1125 K  = 1,KM
      DO 1125 I  = 1,IM
      DO 1126 J  = 2,JMM1
      VNODE(I,J,K)   = 0.125*(TEMP2(I,J-1,K) + TEMP2(I,J,K))
 1126 CONTINUE
      VNODE(I,1,K)   = 0.375*TEMP2(I,1,K)    - 0.125*TEMP2(I,2,K)
      VNODE(I,JM,K)  = 0.375*TEMP2(I,JMM1,K) - 0.125*TEMP2(I,JMM2,K)
 1125 CONTINUE
C
C
      RETURN
      END
C
C**************************************************************************************
C
      SUBROUTINE SET_XLENGTH
C
C********************************************************************************
C********************************************************************************
C    THIS SUBROUTINE CALCULATES THE DISTANCE FROM EVERY GRID POINT TO THE NEAREST
C    SOLID WALL AND USES THIS TO SET THE MIXING LENGTH.
C********************************************************************************
C********************************************************************************
      INCLUDE 'commall-open-18.3'
C
C   THE NEAREST WALL IS SOUGHT OVER THE RANGE  J  +/- JRANGE,  K  +/- KRANGE
C
      WRITE(6,*) ' STARTING SET_XLENGTH. '
C
C     JRANGE AND KRANGE ARE THE RANGE OF POINTS OVER WHICH THE NEAREST 
C     WALL POINT IS SOUGHT
      JRANGE = 20
      KRANGE = 5
C
C
      DO 1000 K=1,KM
      DO 1000 J=1,JM
      NR       = NROW(J)
      PITCH    = 6.283185*R(J,K)/NBLADE(J)
      DIST_REF = 0.5*PITCH
C
C********************************************************************************
C********************************************************************************
      DO 1000 I=1,IM

C     FIND THE DISTANCE TO THE HUB OR CASING
C
      IWALL = IMID
      J1 = J - JRANGE
      IF(J1.LT.2)  J1 = 2
      J2 = J + JRANGE
      IF(J2.GT.JM) J2 = JM
C
      HUBDIST = 1.0E10
      TIPDIST = 1.0E10
      DMIN    = 1.0E10
C
C   FIND THE NEAREST POINT ON THE HUB
      DO 10 JWALL = J1,J2
      XDIF  = X(J,K) - X(JWALL,1)
      RDIF  = R(J,K) - R(JWALL,1)
      DISTSQ = XDIF*XDIF + RDIF*RDIF
      IF(DISTSQ.LT.DMIN) THEN
           DMIN = DISTSQ
           JMIN = JWALL
      END IF
   10 CONTINUE
C
C    FIND THE PERPENDICULAR DISTANCE TO THE NEAREST POINT ON THE HUB
      ATOT  = SQRT(ASX(IWALL,JMIN,1)*ASX(IWALL,JMIN,1)
     &      + ASR(IWALL,JMIN,1)*ASR(IWALL,JMIN,1))
      XNORM = ASX(IWALL,JMIN,1)/ATOT
      RNORM = ASR(IWALL,JMIN,1)/ATOT
      XDIF  = X(J,K) - X(JMIN,1)
      RDIF  = R(J,K) - R(JMIN,1)
      HUBDIST = ABS(XDIF*XNORM + RDIF*RNORM)
      IF(IBOUND.EQ.1.OR.IBOUND.GT.2) HUBDIST = DIST_REF
C
C   FIND THE NEAREST POINT ON THE CASING
      DMIN    = 1.0E10
      DO 15 JWALL = J1,J2
      XDIF  = X(J,K) - X(JWALL,KM)
      RDIF  = R(J,K) - R(JWALL,KM)
      DISTSQ = XDIF*XDIF + RDIF*RDIF
      IF(DISTSQ.LT.DMIN) THEN
           DMIN = DISTSQ
           JMIN = JWALL
      END IF
   15 CONTINUE
C
C    FIND THE PERPENDICULAR DISTANCE TO THE NEAREST POINT ON THE CASING
      ATOT  = SQRT(ASX(IWALL,JMIN,KM)*ASX(IWALL,JMIN,KM)
     &      + ASR(IWALL,JMIN,KM)*ASR(IWALL,JMIN,KM))
      XNORM = ASX(IWALL,JMIN,KM)/ATOT
      RNORM = ASR(IWALL,JMIN,KM)/ATOT
      XDIF  = X(J,K) - X(JMIN,KM)
      RDIF  = R(J,K) - R(JMIN,KM)
      TIPDIST = ABS(XDIF*XNORM + RDIF*RNORM)
      IF(IBOUND.GE.2) TIPDIST = DIST_REF
C
C   STORE THE DISTANCE TO THE NEAREST END WALL
C
       ENDWALL_DIST = HUBDIST*TIPDIST/(HUBDIST + TIPDIST)
C
C*************************************************************************
C*************************************************************************
C     NEXT FIND THE DISTANCE TO THE NEAREST BLADE SURFACE
C
      J1 = J - JRANGE
      IF(J1.LT.2) J1 = 2
      IF(J1.LT.JSTART(NR)) J1 = JSTART(NR)
      IF(J1.GE.JTE(NR))    J1 = JTE(NR) - 1
      J2 = J + JRANGE
      IF(J2.GT.JM) J2 = JM
      IF(J2.GT.JMIX(NR)) J2 = JMIX(NR)
      IF(J2.LE.JLE(NR))  J2 = JLE(NR) + 1
      K1 = K - KRANGE
      IF(K1.LT.1) K1 = 1
      K2 = K + KRANGE
      IF(K2.GT.KM) K2 = KM
C
C    FIND THE NEAREST POINT ON THE  I = 1 BLADE SURFACE
C
      DMIN = 1.0E10
      IF_FOUND = 0
      DO 20 JSURF = J1,J2
      DO 25 KSURF = K1,K2  
      XDIF  = X(J,K)  - X(JSURF,KSURF)
      RDIF  = R(J,K)  - R(JSURF,KSURF)
C
C      TDIF  = RTHETA(I,J,K) - RTHETA(1,JSURF,KSURF)
C      Changed by WS/LX  17 May, 2017  
      TDIFSQ = 4.* R(J,K)*R(JSURF,KSURF)*SIN((THETA(I,J,K)
     &       - THETA(1,JSURF,KSURF))/2.)**2
      DISTSQ = XDIF*XDIF + RDIF*RDIF + TDIFSQ
C
      IF(DISTSQ.LT.DMIN) THEN
           DMIN = DISTSQ
           JMIN = JSURF
           KMIN = KSURF
           IF_FOUND = 1
      END IF
   25 CONTINUE
   20 CONTINUE
C
      IF(IF_FOUND.EQ.1) THEN
C     FIND THE PERPENDICULAR DISTANCE TO THE NEAREST POINT ON THE I=1 BLADE SURFACE.
      KNORM = KMIN
      IF(KMIN.EQ.KM) KNORM = KM-1
      ATOT  = SQRT(ABX(1,JMIN,KNORM)*ABX(1,JMIN,KNORM)  
     &      + ABR(1,JMIN,KNORM)*ABR(1,JMIN,KNORM)
     &      + ABT(JMIN,KNORM)*ABT(JMIN,KNORM))
      XNORM  = ABX(1,JMIN,KNORM)/ATOT
      RNORM  = ABR(1,JMIN,KNORM)/ATOT
      TNORM  = ABT(JMIN,KNORM)/ATOT
      XDIF   = X(J,K)  - X(JMIN,KMIN)
      RDIF   = R(J,K)  - R(JMIN,KMIN)
C
C      TDIF   = RTHETA(I,J,K) - RTHETA(1,JMIN,KMIN)
C      Changed by WS/LX  17 May, 2017  
      TDIF = 2.* SQRT(R(J,K)*R(JMIN,KMIN))*SIN((THETA(I,J,K)
     &       - THETA(1,JMIN,KMIN))/2.) 
C
      SSDIST = ABS(XDIF*XNORM + RDIF*RNORM + TDIF*TNORM)
C   JDD ADDED 11/13 TO IMPROVE TREATMENT OF THIN LEADING AND TRAILING EDGES.
      IF(JMIN.LE.JLE(NR).OR.JMIN.GE.JTE(NR))
     &       SSDIST = SQRT(XDIF*XDIF + RDIF*RDIF+ TDIF*TDIF)
      ELSE
             SSDIST = DIST_REF
      END IF
C
C    FIND THE NEAREST POINT ON THE  I = IM BLADE SURFACE
C
      DMIN = 1.0E10
      IF_FOUND = 0
      DO 30 JSURF = J1,J2
      DO 35 KSURF = K1,K2  
      XDIF  = X(J,K)  - X(JSURF,KSURF)
      RDIF  = R(J,K)  - R(JSURF,KSURF)
C
C      TDIF  = RTHETA(I,J,K) - RTHETA(IM,JSURF,KSURF)  
C      Changed by WS/LX  17 May, 2017  
      TDIFSQ = 4.* R(J,K)*R(JSURF,KSURF)*SIN((THETA(I,J,K)
     &       - THETA(IM,JSURF,KSURF))/2.)**2
      DISTSQ = XDIF*XDIF + RDIF*RDIF + TDIFSQ
C
      IF(DISTSQ.LT.DMIN) THEN
           DMIN = DISTSQ
           JMIN = JSURF
           KMIN = KSURF
           IF_FOUND = 1
      END IF
   35 CONTINUE
   30 CONTINUE
C
C     FIND THE PERPENDICULAR DISTANCE TO THE NEAREST POINT ON THE I=IM BLADE SURFACE.
      KNORM = KMIN
      IF(IF_FOUND.EQ.1) THEN
      IF(KMIN.EQ.KM) KNORM = KM-1
      ATOT = SQRT(ABX(IM,JMIN,KNORM)*ABX(IM,JMIN,KNORM)  
     &     + ABR(IM,JMIN,KNORM)*ABR(IM,JMIN,KNORM)
     &     + ABT(JMIN,KNORM)*ABT(JMIN,KNORM))
      XNORM = ABX(IM,JMIN,KNORM)/ATOT
      RNORM = ABR(IM,JMIN,KNORM)/ATOT
      TNORM = ABT(JMIN,KNORM)/ATOT
      XDIF  = X(J,K)  - X(JMIN,KMIN)
      RDIF  = R(J,K)  - R(JMIN,KMIN)
C
C      TDIF  = RTHETA(I,J,K) - RTHETA(IM,JMIN,KMIN)
C      Changed by WS/LX  17 May, 2017  
      TDIF = 2.* SQRT(R(J,K)*R(JMIN,KMIN))*SIN((THETA(I,J,K)
     &       - THETA(IM,JMIN,KMIN))/2.)
C
      PSDIST = ABS(XDIF*XNORM + RDIF*RNORM + TDIF*TNORM)
C   JDD ADDED 11/13 TO IMPROVE TREATMENT OF THIN LEADING AND TRAILING EDGES.
      IF(JMIN.LE.JLE(NR).OR.JMIN.GE.JTE(NR))
     &      PSDIST = SQRT(XDIF*XDIF + RDIF*RDIF+ TDIF*TDIF)
      ELSE
            PSDIST = DIST_REF
      END IF
C
C     STORE THE DISTANCE TO THE NEAREST BLADE SURFACE.
C
      BLADE_DIST = SSDIST*PSDIST/(SSDIST + PSDIST)
C
C*******************************************************************************  
C*******************************************************************************    
C   SET THE DISTANCE FROM THE NEAREST SOLID SURFACE AS  XDIST*YDIST/SQRT(XDIST**2 + YDIST**2)
C
      IF(K.GT.1.AND.K.LT.KM) THEN
      DISTSQRT  = SQRT(ENDWALL_DIST*ENDWALL_DIST
     &          + BLADE_DIST*BLADE_DIST)
      IF(DISTSQRT.LT.1.0E-10) DISTSQRT = 1.0E-10
      DIST_MIN(I,J,K) = BLADE_DIST*ENDWALL_DIST/DISTSQRT
      ELSE
      DIST_MIN(I,J,K) = 0.0
      END IF
C
C  Q3D
      IF(KM.EQ.2) DIST_MIN(I,J,K) = BLADE_DIST
      IF(IM.EQ.2) DIST_MIN(I,J,K) = ENDWALL_DIST
C  END Q3D
C*******************************************************************************  
C*******************************************************************************
C   SET THE MIXING LENGTH LIMIT, XLIMIT(J) VARYING WITH MERIDIONAL DISTANCE.
C
      IF(I.EQ.IMID.AND.K.EQ.KMID) THEN
            IF(J.LE.JLE(NR))  XLLIM = XLLIM_IN(NR) +
     &      (SMERID(J,K) - SMERID(JSTART(NR),K))
     &     /(SMERID(JLE(NR),K) - SMERID(JSTART(NR),K))
     &     *(XLLIM_LE(NR) - XLLIM_IN(NR))
            IF(J.GT.JLE(NR).AND.J.LE.JTE(NR)) XLLIM = XLLIM_LE(NR) +
     &      (SMERID(J,K) - SMERID(JLE(NR),K))
     &     /(SMERID(JTE(NR),K) - SMERID(JLE(NR),K))
     &     *(XLLIM_TE(NR) - XLLIM_LE(NR))
            IF(J.GT.JTE(NR))  XLLIM = XLLIM_TE(NR) +
     &      (SMERID(J,K) - SMERID(JTE(NR),K))
     &     /(SMERID(JMIX(NR),K) - SMERID(JTE(NR),K))
     &     *(XLLIM_DN(NR) - XLLIM_TE(NR))
C
C     USE A FACTOR OF 2 ON XLIMIT SO THAT IS ROUGHLY THE THROAT WIDTH.
C
            XLIMIT(J) = 2.0*DIST_MIN(I,J,K)*XLLIM
C
      END IF 
C   END OF LOOP OVER ALL  I,J, K,  POINTS .
 1000 CONTINUE
C
C*******************************************************************************  
C*******************************************************************************  
C*******************************************************************************  
C     BLEND THE LINEAR REGION OF THE MIXING LENGTH AND THE CONSTANT REGION AT 
C     RATIO = DBLEND.
C     CHECK IF THE MIXING LENGTH IS GREATER THAN THE MIXING LENGTH LIMIT
C     THEN STORE THE MIXING LENGTH FOR EACH POINT AS STORE(I,J,K)
C
      DBLEND = 0.6666666
      EXPON  = DBLEND/(1.0 - DBLEND)
      FBLEND = (1.0 - DBLEND)*DBLEND**EXPON
C     
      DO 2100 J = 1,JM
      DO 2000 K = 1,KM
      DO 2000 I = 1,IM
      X_LENGTH  = DIST_MIN(I,J,K)
      RATIO   = X_LENGTH/XLIMIT(J)
      IF(RATIO.GT.DBLEND)
     &           X_LENGTH = XLIMIT(J)*(1.0 - FBLEND*RATIO**(-EXPON))         
      STORE(I,J,K) = X_LENGTH
C     TEMP3  IS USED IN THE NEXT SECTION TO DAMP THE FREE STREAM TURBULENCE NEAR A WALL.
      TEMP3(I,J,K) = (X_LENGTH/XLIMIT(J))**4
 2000 CONTINUE
 2100 CONTINUE
C
C*******************************************************************************      
C    NOW AVERAGE THE MIXING LENGTH FOR A CELL AND STORE THE AVERAGE 
C    MULTIPLIED BY THE VON KARMEN CONSTANT  (= 0.41) ALL SQUARED.
C
      DO 2200 K = 1,KMM1
      DO 2200 J = 2,JM
      NR = NROW(J)
      DO 2200 I = 1,IMM1
      AVG_XLENGTH = 0.125*(STORE(I,J,K) + STORE(I,J,K+1)
     &     + STORE(I,J-1,K) + STORE(I,J-1,K+1) + STORE(I+1,J,K)
     &     + STORE(I+1,J,K+1) + STORE(I+1,J-1,K) + STORE(I+1,J-1,K+1))
      XLENGTH(I,J,K) = 0.41*0.41*AVG_XLENGTH*AVG_XLENGTH
C
C     SET THE MIXING LENGTH TO ZERO FOR THE FIRST CELL ON A BLADE ROW.
      IF(J.EQ.JSTART(NR)) XLENGTH(I,J,K) = 0.0 
C   SET THE FUNCTION USED TO DAMP THE FREE STREAM TURBULENCE NEAR A WALL.
      AVG_FSTRAT = 0.125*(TEMP3(I,J,K) + TEMP3(I,J,K+1)
     &     + TEMP3(I,J-1,K) + TEMP3(I,J-1,K+1) + TEMP3(I+1,J,K)
     &     + TEMP3(I+1,J,K+1) + TEMP3(I+1,J-1,K) + TEMP3(I+1,J-1,K+1))
      FST_RAT(I,J,K) = AVG_FSTRAT*FSTURB(NR)
C     SET THE FREE STREAM TURBULENCE TO ZERO FOR THE FIRST CELL ON A BLADE ROW.
      IF(J.EQ.JSTART(NR)) FST_RAT(I,J,K) = 0.0 
C
 2200 CONTINUE
C
C*******************************************************************************  
C     SET THE AVERAGE WALL DISTANCE FOR A CELL , SQUARED .
C
      DO 2300 K=1,KMM1
      DO 2300 J=2,JM
      DO 2300 I=1,IMM1
      AVGDIST = 0.125*(DIST_MIN(I,J,K) + DIST_MIN(I,J,K+1)
     &   + DIST_MIN(I,J-1,K)   + DIST_MIN(I,J-1,K+1) + DIST_MIN(I+1,J,K)
     &   + DIST_MIN(I+1,J,K+1) + DIST_MIN(I+1,J-1,K)
     &   + DIST_MIN(I+1,J-1,K+1))
      DWALLSQ(I,J,K) = AVGDIST*AVGDIST
 2300 CONTINUE
C
C*******************************************************************************  
C
      WRITE(6,*)
      WRITE(6,*) ' LEAVING SUBROUTINE    SET_XLENGTH. '
      WRITE(6,*)
C
      RETURN
      END
C******************************************************************************
C*************************************************************************************
C******************************************************************************
C
      SUBROUTINE NEW_LOSS
C
C******************************************************************************
C************THIS SUBROUTINE COMPUTES A BODY FORCE BASED ON WALL FUNCTIONS FOR THE
C            SURFACE SHEAR STRESS AND A MIXING LENGTH MODEL OF EDDY VISCOSITY.
C
C            NEW MODEL BASED ON TBLOCK MODEL. JANUARY 2010.
C
C
      INCLUDE  'commall-open-18.3'
C
      COMMON/BKSTRESS/  TXX(ID,JD,KD),TXR(ID,JD,KD),
     &                  TXT(ID,JD,KD),TRX(ID,JD,KD),TRR(ID,JD,KD),
     &                  TRT(ID,JD,KD),TTX(ID,JD,KD),TTR(ID,JD,KD),
     &                  TTT(ID,JD,KD),QXX(ID,JD,KD),QRR(ID,JD,KD),
     &                  QTT(ID,JD,KD)
C
      DIMENSION  FORCEX(ID,JD,KD),FORCER(ID,JD,KD),FORCET(ID,JD,KD),
     &           ESOURCE(ID,JD,KD),TEMPP(JD)
C
C      CALCULATE THE VISCOSITY OVER THE FIRST QUARTER OF THE STEPS.
C      THEN HOLD IT CONSTANT FOR THE REMAINDER OF THE STEPS.
C
      IF(NSTEP.EQ.1.OR.NSTEP.LT.NSTEPS_MAX/4) THEN
C
      IF(REYNO.GT.100.) THEN
           J1=JLE(1)
           J2=JTE(1)
           IF(JLE(1).GT.JM) J1=1
           IF(JTE(1).GT.JM) J2=JM
           XCHORD  = SMERID(J2,KMID)-SMERID(J1,KMID)
           ROW2    = SQRT(ROVX(IMID,J2,KMID)*ROVX(IMID,J2,KMID)
     &             + ROVR(IMID,J2,KMID)*ROVR(IMID,J2,KMID) 
     &             + ROWT(IMID,J2,KMID)*ROWT(IMID,J2,KMID) )
           VISLAM  = XCHORD*ROW2/REYNO
      END IF
C
      IF(REYNO.GT.0.0.AND.REYNO.LT.99.99) THEN
            VISLAM = REYNO/100000.
      END IF
C
      IF(REYNO.LT.0.0) VISLAM = -REYNO*1.0E-5
C
      THCOND  = CP*VISLAM/PRANDTL
      FTCOND = THCOND/VISLAM
C
C   END OF PART ONLY CALLED OVER FIRST 1/4 OF THE STEPS.
      ENDIF
C
C********************************************************************************
C     SET THE VISCOUS FORCES AND ENERGY SOURCES TO ZERO
C
      DO 100 K=1,KMM1
      DO 100 J=2,JM
      DO 100 I=1,IMM1
      FORCEX(I,J,K)  = 0.0
      FORCER(I,J,K)  = 0.0
      FORCET(I,J,K)  = 0.0
      ESOURCE(I,J,K) = 0.0
  100 CONTINUE
C
C*******************************************************************************************
C*******************************************************************************************
C     START TO EVALUATE THE TURBULENT VISCOSITY FOR EACH CELL.
C
      DO 150 K=1,KMM1
      KP1 = K+1
      DO 150 J=2,JM
      JM1 = J-1
      NRW = NROW(J)
      JTEDGE = JTE(NRW)
      JLEDGE = JLE(NRW)
      DO 150 I=1,IMM1
      IP1 = I+1
C
C    AVERAGE THE CONDITIONS FOR THE CELLS. Note these are true averages.
C
      WTAVG = 0.125*(WT(I,J,K)+WT(IP1,J,K)+WT(IP1,J,KP1)+WT(I,J,KP1)
     &      + WT(I,JM1,K)+WT(IP1,JM1,K)+WT(IP1,JM1,KP1)+WT(I,JM1,KP1))
      ROAVG = ROAVG_CELL(I,J,K)
      TSAVG = 0.125*(T_STATIC(I,J,K)+T_STATIC(IP1,J,K)
     &        +T_STATIC(IP1,J,KP1)+T_STATIC(I,J,KP1)+T_STATIC(I,JM1,K)
     &        +T_STATIC(IP1,JM1,K)+T_STATIC(IP1,JM1,KP1)
     &        +T_STATIC(I,JM1,KP1))
      RAVG  = RAVG_CELL(J,K)
C
C     CALCULATE THE DERIVATIVES OF THE VELOCITY COMPONENTS AND TEMPERATURE.
C     NOTE THE VORTICITY IS BASED ON THE RELATIVE VELOCITY,  WT.
C     THIS WAS CHANGED FROM THE ABSOLUTE VELOCITY BY JDD AUGUST 2017.
      CALL GRADVEL(I,J,K,VX,DVXDX,DVXDR,DVXDT)
      CALL GRADVEL(I,J,K,VR,DVRDX,DVRDR,DVRDT)
      CALL GRADVEL(I,J,K,WT,DWTDX,DWTDR,DWTDT)
      CALL GRADVEL(I,J,K,T_STATIC,DTEMPDX,DTEMPDR,DTEMPDT)
C
C     CALCULATE THE RATES OF STRAIN FOR CALCULATING THE VISCOUS STRESSES.
C     THESE ARE NOT YET THE STRESSES AS THE FINAL VISCOSITY IS NOT YET KNOWN.
C     
      TXR(I,J,K) = DVXDR + DVRDX
      TXX(I,J,K) = DVXDX
      TXT(I,J,K) = DWTDX + DVXDT
      TRR(I,J,K) = DVRDR
      TRT(I,J,K) = DWTDR + DVRDT - WTAVG/RAVG
      TTT(I,J,K) = DWTDT
      QXX(I,J,K) = DTEMPDX
      QRR(I,J,K) = DTEMPDR
      QTT(I,J,K) = DTEMPDT
C
C        USE THE VORTICITY TO DETERMINE THE TURBULENT VISCOSITY FOR EACH CELL.
C
      VORTX   = DVRDT - DWTDR - WTAVG/RAVG 
      VORTR   = DWTDX - DVXDT
      VORTT   = DVXDR - DVRDX
      VORT    = SQRT(VORTX*VORTX + VORTR*VORTR + VORTT*VORTT)
C  SET A LIMIT TO THE VORTICITY
      IF(VORT.GT.VORT_MAX) VORT = VORT_MAX
C
      ABS_VORT(I,J,K) = VORT
C
C   CALCULATE THE TURBULENT VISCOSITY IN THE BOUNDARY LAYERS.
C
      TURBVIS_WALL = ROAVG*XLENGTH(I,J,K)*VORT*FMIXUP
C.
C    Use FYPLUS to reduce the turbulent viscosity in the buffer region. 
C  
      IF(Y_PLUS(I,J,K).LE.YPLAM) FYPLUS = 0.0
      IF(Y_PLUS(I,J,K).GT.YPLAM)  THEN
               XFAC = (Y_PLUS(I,J,K) - YPLAM)/(YPTURB- YPLAM)
               IF(XFAC.GT.1.0) XFAC = 1.0
               FYPLUS = XFAC*XFAC*(3.0 - 2.0*XFAC)
      END IF
      TURBVIS_WALL = TURBVIS_WALL*FYPLUS
C
C     SET THE LAMINAR VISCOSITY IF IT VARIES WITH TEMPERATURE.
C
      IF(REYNO.LT.0.0001) THEN
            VISLAM = (ABS(REYNO)/100000.0) * (TSAVG/288.0)**0.62
      END IF
C
C    SAVE THE VISCOSITY IT IS ONLY USED TO CALCULATE AND WRITE OUT THE REYNOLDS NUMBER.
      IF(I.EQ.IMID.AND.K.EQ.KMID.AND.J.EQ.JTEDGE) VISCOSY(NRW) = VISLAM
C
C     CALCULATE THE REYNOLDS NUMBER AT MID PASSAGE AT THE TRAILING EDGE.
C
      IF(I.EQ.IMID.AND.K.EQ.KMID.AND.J.EQ.JTEDGE) THEN
           WABS       = SQRT(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K)
     &                + WT(I,J,K)*WT(I,J,K) )
           ROVEXIT(NRW)  = RO(I,J,K)*WABS
           VISCOSY(NRW)  = VISLAM
           REYNOLDS(NRW) = CHORD(NRW)*ROVEXIT(NRW)/VISLAM
      END IF
C
C  SET THE TOTAL VISCOSITY INCLUDING THAT DUE TO FREE STREAM TURBULENCE.
C
      VISTOT  = VISLAM*(1 + FST_RAT(I,J,K) ) + TURBVIS_WALL 
C
C    SET A LIMIT TO THE TURBULENT VISCOSITY
C
      VISLIM = TURBVIS_LIM*VISLAM
      IF(VISTOT.GT.VISLIM) VISTOT = VISLIM
C
C     SAVE THE TURBULENT/LAMINAR VISCOSITY RATIO AS "VISC_RAT"
C
       VISC_RAT(I,J,K)   = VISTOT/VISLAM
       VISC_LAM(I,J,K)   = VISLAM
C
  150 CONTINUE
C
C    END OF SETTING THE TURBULENT VISCOSITY
C****************************************************************************
C****************************************************************************
C     NEXT  CHECK FOR TRANSITION.
C
      DO 160 J = 2,JM
C
      NRW    = NROW(J)
      J1     = JSTART(NRW)
      JROW   = J - J1 + 1
      JTRHUB = JTRAN_K1(NRW)
      JTRTIP = JTRAN_KM(NRW)
      JTRLOW = JTRAN_I1(NRW)
      JTRUP  = JTRAN_IM(NRW)
C
C  Q3D
      IF(KM.EQ.2) GO TO 171
C  END Q3D
C
C   FIRST ON THE HUB
C
      DO 170 I = 1,IMM1
      IF(JROW.LT.JTRHUB) GO TO 34
      RMAX = 0.0
      DO 37 K=1,KMID
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 37
      RMAX = RATVIS
   37 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 38
   34 CONTINUE
      DO 39 K = 1,KMID
   39 VISC_RAT(I,J,K) = 1.0
   38 CONTINUE
C
C  NEXT ON THE CASING
C
      IF(JROW.LT.JTRTIP) GO TO 46
      RMAX = 0.0
      DO 47 K = KMID,KMM1
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 47
      RMAX = RATVIS
   47 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 49
   46 CONTINUE
      DO 48 K=KMID,KMM1
   48 VISC_RAT(I,J,K) = 1.0
   49 CONTINUE
C
  170 CONTINUE
C
  171 CONTINUE
C
C    NEXT ON THE LOWER - I=1, BLADE SURFACE
C
      DO 180  K=1,KMM1
      IF(JROW.LT.JTRLOW) GO TO 70
      RMAX=0.0
      DO 56 I=1,IMID
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 56
      RMAX = RATVIS
   56 CONTINUE
      IF(RMAX.GT.FTRANS)  GO TO 71
   70 CONTINUE
      DO 72 I = 1,IMID
   72 VISC_RAT(I,J,K)= 1.0
   71 CONTINUE
C
C    NEXT ON THE UPPER, I = IM, BLADE SURFACE
C
      IF(JROW.LT.JTRUP) GO TO 69
      RMAX = 0.0
      DO 73 I = IMID,IMM1
      RATVIS  = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 73
      RMAX = RATVIS
   73 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO  77
   69 CONTINUE
      DO 74 I  = IMID,IMM1
   74 VISC_RAT(I,J,K) = 1.0
   77 CONTINUE
C
  180 CONTINUE
C
  160 CONTINUE
C
C********************************************************************************
C********************************************************************************
C
C     NOW THE FINAL VISCOSITY IS FIXED SET THE STRESSES AND HEAT FLOWS IN ALL ELEMENTS
C     JDD WARNING the div V term in the normal stresses is not yet included.
C
      DO 190 K = 1,KMM1
      DO 190 J = 2,JM
      DO 190 I = 1,IMM1
C
      VISLAM  = VISC_LAM(I,J,K)
      VISTOT  = VISC_RAT(I,J,K)*VISLAM
      THCOND  = CP*VISTOT/PRANDTL
      VISTOT2 = 2.0*VISTOT
C     
      TXR(I,J,K) = VISTOT*TXR(I,J,K)
      TXX(I,J,K) = VISTOT2*TXX(I,J,K)
      TXT(I,J,K) = VISTOT*TXT(I,J,K)
      TRX(I,J,K) = TXR(I,J,K)
      TRR(I,J,K) = VISTOT2*TRR(I,J,K)
      TRT(I,J,K) = VISTOT*TRT(I,J,K)
      TTX(I,J,K) = TXT(I,J,K)
      TTT(I,J,K) = VISTOT2*TTT(I,J,K)
      TTR(I,J,K) = TRT(I,J,K)
      QXX(I,J,K) = -THCOND*QXX(I,J,K)
      QRR(I,J,K) = -THCOND*QRR(I,J,K)
      QTT(I,J,K) = -THCOND*QTT(I,J,K)
C
  190 CONTINUE
C
C******************************************************************************
C******************************************************************************
C    EVALUATE AND SMOOTH THE PRESSURE GRADIENTS IF YPLUSWALL IS < -10.0.
C
      IF(YPLUSWALL.LT.-10.0) CALL SET_PWALLGRAD
C
C***************************************************************************************
C****************************************************************************************
C     CALCULATE THE WALL SHEAR STRESSES ON ALL THE SOLID SURFACES.
C
      VISLAM16  = VISLAM*16.0
C
C******************FIRST FOR BLADE SURFACES.****************************
C
      DO 200 K=1,KMM1
      KP1 = K+1
      DO 200 J=2,JM
C
C   SKIP IF NOT ON A BLADE SURFACE
      IF(IND(J).EQ.0)  GO TO 220
C  
      JM1 = J-1
      NRW = NROW(J)
C 
C   ALSO SKIP FOR CELLS ABOVE THE TIP GAP
      IF( (KTIPS(NRW).GT.0).AND. 
     &    (K.GE.KTIPS(NRW)).AND.(K.LT.KTIPE(NRW)) ) 
     &    GO TO 220
C
C*********************FIRST FOR THE I = 1  BLADE SURFACE.*******************************
C                  CALCULATE THE WALL SHEAR STRESS.
C
      VXAVG1 = VX(2,J,K)+VX(2,JM1,K)+VX(2,JM1,KP1)+VX(2,J,KP1)
      VRAVG1 = VR(2,J,K)+VR(2,JM1,K)+VR(2,JM1,KP1)+VR(2,J,KP1)
      WTAVG1 = WT(2,J,K)+WT(2,JM1,K)+WT(2,JM1,KP1)+WT(2,J,KP1)
      ROAVG1 = RO(2,J,K)+RO(2,JM1,K)+RO(2,JM1,KP1)+RO(2,J,KP1)
      WAVG   = SQRT(WTAVG1*WTAVG1 + VXAVG1*VXAVG1 + VRAVG1*VRAVG1)
      AREA1  = SQRT(ABX(1,J,K)*ABX(1,J,K) + ABR(1,J,K)*ABR(1,J,K)
     &            + ABT(J,K)*ABT(J,K)) 
      DPERP  = VOL(1,J,K)/AREA1
C
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN
C
         YPLUS_OLD =  YPLUS_I1(J,K)
         DENSITY   = 0.25*ROAVG1
         WSLIP     = 0.25*WAVG
         CALL WALLFUN(1,J,K,1,DPERP,DPDS_CELL(1,J,K),DENSITY,
     &                TWALLI1,YPLUS_OLD,WSLIP,YPLUS_NEW)
          YPLUS_I1(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 365
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
C
      RE     = DPERP*ROAVG1*WAVG/VISLAM16
      RELOG  = 1./ALOG(RE)
C******************************************************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_L(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C******************************************************************************
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLI1    = CF*ROAVG1*WAVG*WAVG/128
           VSTAR   = SQRT(TWALLI1/(ROAVG1/4))
           PLUSK   = ROUGH*VSTAR*ROAVG1/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C******************************************************************************
C   SET THE SHEAR STRESS ON THE I = 1 BLADE SURFACE.
C
      TWALLI1    = CF*ROAVG1*WAVG*WAVG/128
C
C******************************************************************************
C
  365 CONTINUE
C
C******************************************************************************
C    CALCULATE  YPLUS FOR USE LATER.
      DO I = 1,IMID
      VSTAR     = SQRT(TWALLI1/(0.25*ROAVG1))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVG1)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(1000.0,YPLS)
      END DO  
C    
      FMULT   = -TWALLI1*AREA1/WAVG
C
      XFORCE1 =  FMULT*VXAVG1
      RFORCE1 =  FMULT*VRAVG1
      TFORCE1 =  FMULT*WTAVG1
      WVISC1  =  TFORCE1*WRAD(J)*RAVG_CELL(J,K)
C
      FORCEX(1,J,K)    = FORCEX(1,J,K)    + XFORCE1
      FORCER(1,J,K)    = FORCER(1,J,K)    + RFORCE1
      FORCET(1,J,K)    = FORCET(1,J,K)    + TFORCE1
      ESOURCE(1,J,K)   = ESOURCE(1,J,K)   + WVISC1
C
C*********************NOW FOR THE I = IM BLADE SURFACE************************
C                     CALCULATE THE WALL SHEAR STRESS 
C
      VXAVGIM = VX(IMM1,J,K)+VX(IMM1,JM1,K)
     &        + VX(IMM1,JM1,KP1)+VX(IMM1,J,KP1)
      VRAVGIM = VR(IMM1,J,K)+VR(IMM1,JM1,K)
     &        + VR(IMM1,JM1,KP1)+VR(IMM1,J,KP1)
      WTAVGIM = WT(IMM1,J,K)+WT(IMM1,JM1,K)
     &        + WT(IMM1,JM1,KP1)+WT(IMM1,J,KP1)
      ROAVGIM = RO(IMM1,J,K)+RO(IMM1,JM1,K)
     &        + RO(IMM1,JM1,KP1)+RO(IMM1,J,KP1)
      WAVG    = SQRT(WTAVGIM*WTAVGIM+VXAVGIM*VXAVGIM+VRAVGIM*VRAVGIM)
C 
      AREAIM  = SQRT(ABX(IM,J,K)*ABX(IM,J,K) + ABR(IM,J,K)*ABR(IM,J,K)
     &              +ABT(J,K)*ABT(J,K)) 
      DPERP    = VOL(IMM1,J,K)/AREAIM
C
C******************************************************************************
C******************************************************************************
C  NEXT THE I = IM , UPPER, SURFACE
C
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD =  YPLUS_IM(J,K)
           DENSITY   = 0.25*ROAVGIM
           WSLIP     = 0.25*WAVG
           CALL WALLFUN(IMM1,J,K,IM,DPERP,DPDS_CELL(IMM1,J,K),
     &             DENSITY,TWALLIM,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_IM(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 375
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
      RE       = DPERP*ROAVGIM*WAVG/VISLAM16
      RELOG    = 1./ALOG(RE)
C******************************************************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_U(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C******************************************************************************
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLIM = CF*ROAVGIM*WAVG*WAVG/128
           VSTAR   = SQRT(TWALLIM/(ROAVGIM/4))
           PLUSK   = ROUGH*VSTAR*ROAVGIM/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C******************************************************************************
C   SET THE SHEAR STRESS ON THE I = IM BLADE SURFACE
C
      TWALLIM    = CF*ROAVGIM*WAVG*WAVG/128
C
C******************************************************************************
C******************************************************************************
C 
  375 CONTINUE
C
C*****************************************************************************
C*****************************************************************************   
C     CALCULATE  YPLUS FOR USE LATER.
      DO I = IMID+1,IMM1
      VSTAR     = SQRT(TWALLIM/(0.25*ROAVGIM))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVGIM)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(1000.0,YPLS)
      END DO
C
      FMULT    = -TWALLIM*AREAIM/WAVG
C 
      XFORCEIM =  FMULT*VXAVGIM
      RFORCEIM =  FMULT*VRAVGIM
      TFORCEIM =  FMULT*WTAVGIM
      WVISCIM  =  TFORCEIM*WRAD(J)*RAVG_CELL(J,K)
C
      FORCEX(IMM1,J,K)  = FORCEX(IMM1,J,K)   + XFORCEIM
      FORCER(IMM1,J,K)  = FORCER(IMM1,J,K)   + RFORCEIM
      FORCET(IMM1,J,K)  = FORCET(IMM1,J,K)   + TFORCEIM
      ESOURCE(IMM1,J,K) = ESOURCE(IMM1,J,K)  + WVISCIM
C
  220 CONTINUE
C
C*****************************************************************************
C
C     NOW SET THE FORCES ON THE PERIODIC CELLS, I=1 and I= IM-1,UPSTREAM OF THE
C     LE AND DOWNSTREAM OF THE TE .
C     AND ALSO FOR CELLS IN THE TIP GAP.
C*****************************************************************************
C     DO NOT TREAT CELLS ON THE SOLID BLADE SURFACES AS PERIODIC.
C
       IF( (IND(J).EQ.1).AND.(K.LT.KTIPS(NRW)).OR.(K.GT.KTIPE(NRW)) )
     &      GO TO 225
C
      I   = 1
      IM1 = IM-1
      XFOR = 0.5*((TXR(I,J,K)+TXR(IM1,J,K))*ABR(I,J,K)
     &    +       (TXT(I,J,K)+TXT(IM1,J,K))*ABT(J,K)
     &    +       (TXX(I,J,K)+TXX(IM1,J,K))*ABX(I,J,K))
      FORCEX(IM1,J,K) = FORCEX(IM1,J,K) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR
      AVG = 0.5*(FORCEX(I,J,K) + FORCEX(IM1,J,K))
      FORCEX(I,J,K)   = AVG
      FORCEX(IM1,J,K) = AVG

      RFOR = 0.5*((TRR(I,J,K)+TRR(IM1,J,K))*ABR(I,J,K)
     &    +       (TRT(I,J,K)+TRT(IM1,J,K))*ABT(J,K)
     &    +       (TRX(I,J,K)+TRX(IM1,J,K))*ABX(I,J,K))
      FORCER(IM1,J,K)  = FORCER(IM1,J,K) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR
      AVG = 0.5*(FORCER(I,J,K) + FORCER(IM1,J,K))
      FORCER(I,J,K)   = AVG
      FORCER(IM1,J,K) = AVG

      TFOR = 0.5*((TTR(I,J,K)+TTR(IM1,J,K))*ABR(I,J,K)
     &    +       (TTT(I,J,K)+TTT(IM1,J,K))*ABT(J,K)
     &    +       (TTX(I,J,K)+TTX(IM1,J,K))*ABX(I,J,K))
      FORCET(IM1,J,K)  = FORCET(IM1,J,K) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
      AVG = 0.5*(FORCET(I,J,K) + FORCET(IM1,J,K))
      FORCET(I,J,K)   = AVG
      FORCET(IM1,J,K) = AVG
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(IM1,J,K))*ABX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(IM1,J,K))*ABR(I,J,K)
     &      +       (QTT(I,J,K) + QTT(IM1,J,K))*ABT(J,K))
C
      VXAVG = VX(I,J,K)+VX(I,JM1,K)+VX(I,JM1,KP1)+VX(I,J,KP1)
      VRAVG = VR(I,J,K)+VR(I,JM1,K)+VR(I,JM1,KP1)+VR(I,J,KP1)
      VTAVG = VT(I,J,K)+VT(I,JM1,K)+VT(I,JM1,KP1)+VT(I,J,KP1)
C
      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)
C
      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(IM1,J,K)= ESOURCE(IM1,J,K) - QFLOW - WVISC 
      AVG = 0.5*(ESOURCE(I,J,K) + ESOURCE(IM1,J,K))
      ESOURCE(I,J,K)   = AVG
      ESOURCE(IM1,J,K) = AVG
C
  225 CONTINUE
C
C******************************************************************************
C******************************************************************************
C
C     NOW SUM THE STRESSES ON THE OTHER I = CONSTANT FACES 
C
      DO 250 I=2,IMM1
      IM1  = I-1
      XFOR = 0.5*((TXR(I,J,K)+TXR(IM1,J,K))*ABR(I,J,K)
     &    +       (TXT(I,J,K)+TXT(IM1,J,K))*ABT(J,K)
     &    +       (TXX(I,J,K)+TXX(IM1,J,K))*ABX(I,J,K))
      FORCEX(IM1,J,K) = FORCEX(IM1,J,K) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR

      RFOR = 0.5*((TRR(I,J,K)+TRR(IM1,J,K))*ABR(I,J,K)
     &    +       (TRT(I,J,K)+TRT(IM1,J,K))*ABT(J,K)
     &    +       (TRX(I,J,K)+TRX(IM1,J,K))*ABX(I,J,K))
      FORCER(IM1,J,K)  = FORCER(IM1,J,K) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR

      TFOR = 0.5*((TTR(I,J,K)+TTR(IM1,J,K))*ABR(I,J,K)
     &    +       (TTT(I,J,K)+TTT(IM1,J,K))*ABT(J,K)
     &    +       (TTX(I,J,K)+TTX(IM1,J,K))*ABX(I,J,K))
      FORCET(IM1,J,K)  = FORCET(IM1,J,K) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(IM1,J,K))*ABX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(IM1,J,K))*ABR(I,J,K)
     &      +       (QTT(I,J,K) + QTT(IM1,J,K))*ABT(J,K))
C
      VXAVG = VX(I,J,K)+VX(I,JM1,K)+VX(I,JM1,KP1)+VX(I,J,KP1)
      VRAVG = VR(I,J,K)+VR(I,JM1,K)+VR(I,JM1,KP1)+VR(I,J,KP1)
      VTAVG = VT(I,J,K)+VT(I,JM1,K)+VT(I,JM1,KP1)+VT(I,J,KP1)
C
      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)
C
      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(IM1,J,K)= ESOURCE(IM1,J,K) - QFLOW - WVISC 

  250 CONTINUE

  200 CONTINUE
C
C**************************************************************************************
C***************************************************************************************
C**********CALCULATE THE WALL SHEAR STRESSES ON THE HUB AND CASING SURFACES**********
C
C   Q3D
      IF(KM.EQ.2) GO TO 401
C   END Q3D
C
      DO 400 J=2,JM
      JM1  = J-1
      NRW  = NROW(J)
      DO 400 I=1,IMM1
      IP1 = I+1
C
C     FIRST FOR THE K = 1,  HUB, WALL.
C
C     CALCULATE THE WALL SHEAR STRESS.
C
      VXAVG1 = VX(I,J,2)+VX(IP1,J,2)+VX(IP1,JM1,2)+VX(I,JM1,2)
      VRAVG1 = VR(I,J,2)+VR(IP1,J,2)+VR(IP1,JM1,2)+VR(I,JM1,2)
      VTAVG1 = VT(I,J,2)+VT(IP1,J,2)+VT(IP1,JM1,2)+VT(I,JM1,2)
      ROAVG1 = RO(I,J,2)+RO(IP1,J,2)+RO(IP1,JM1,2)+RO(I,JM1,2)
      RAVG1  = R(J,2)   + R(J,2)    + R(JM1,2)    + R(JM1,2)
      WTAVG1 = VTAVG1 - WHUB(J)*RAVG1
      WAVG   = SQRT(VXAVG1*VXAVG1 + VRAVG1*VRAVG1 + WTAVG1*WTAVG1)
      AREA1  = SQRT(ASX(I,J,1)*ASX(I,J,1) + ASR(I,J,1)*ASR(I,J,1))
      DPERP  = VOL(I,J,1)/AREA1
C
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD =  YPLUS_K1(I,J)
           DENSITY   = 0.25*ROAVG1
           WSLIP     = 0.25*WAVG
         CALL WALLFUN(I,J,1,1,DPERP,DPDS_CELL(I,J,1),DENSITY,
     &                TWALLK1,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_K1(I,J) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 345

      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
      RE     = DPERP*ROAVG1*WAVG/VISLAM16
      RELOG  = 1./ALOG(RE)
C******************************************************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_H(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C******************************************************************************
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLK1     = CF*ROAVG1*WAVG*WAVG/128
           VSTAR   = SQRT(TWALLK1/(ROAVG1/4))
           PLUSK   = ROUGH*VSTAR*ROAVG1/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C******************************************************************************
C  SET THE SHEAR STRESS ON THE K=1, (HUB) , ENDWALL .
C
      TWALLK1    = CF*ROAVG1*WAVG*WAVG/128
C
C******************************************************************************
C
  345 CONTINUE
C
C******************************************************************************
C   SET YPLUS FOR USE LATER
      IF(IBOUND.EQ.0.OR.IBOUND.EQ.2) THEN
      DO K = 1,KMID
      VSTAR     = SQRT(TWALLK1/(0.25*ROAVG1))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVG1)/VISLAM
      YPLSP     = Y_PLUS(I,J,K)
      Y_PLUS(I,J,K) = AMIN1(YPLS,YPLSP)
      END DO
      END IF
C
C******************************************************************************
C
      FMULT   = -TWALLK1*AREA1/WAVG
      IF(IBOUND.EQ.1.OR.IBOUND.GT.2) FMULT = 0.0
C
      XFORCE1 = FMULT*VXAVG1
      RFORCE1 = FMULT*VRAVG1
      TFORCE1 = FMULT*WTAVG1
      WVISC1  = TFORCE1*0.25*WHUB(J)*RAVG1
C
      FORCEX(I,J,1)    = FORCEX(I,J,1)    + XFORCE1
      FORCER(I,J,1)    = FORCER(I,J,1)    + RFORCE1
      FORCET(I,J,1)    = FORCET(I,J,1)    + TFORCE1
      ESOURCE(I,J,1)   = ESOURCE(I,J,1)   + WVISC1
C
C
C**************NOW FOR THE K = KM (CASING) END WALL**********************************
C
C     CALCULATE THE WALL SHEAR STRESS.
C
      VXAVGKM = VX(I,J,KMM1)+VX(IP1,J,KMM1)
     &        + VX(IP1,JM1,KMM1)+VX(I,JM1,KMM1)
      VRAVGKM = VR(I,J,KMM1)+VR(IP1,J,KMM1)
     &        + VR(IP1,JM1,KMM1)+VR(I,JM1,KMM1)
      VTAVGKM = VT(I,J,KMM1)+VT(IP1,J,KMM1)
     &        + VT(IP1,JM1,KMM1)+VT(I,JM1,KMM1)
      ROAVGKM = RO(I,J,KMM1)+RO(IP1,J,KMM1)
     &        + RO(IP1,JM1,KMM1)+RO(I,JM1,KMM1)
      RAVGKM  = R(J,KMM1) + R(J,KMM1) + R(JM1,KMM1) + R(JM1,KMM1)
      WTAVGKM = VTAVGKM   - WTIP(J)*RAVGKM
      WAVG    = SQRT(VXAVGKM*VXAVGKM+VRAVGKM*VRAVGKM+WTAVGKM*WTAVGKM)
      AREAKM  = SQRT(ASX(I,J,KM)*ASX(I,J,KM) + ASR(I,J,KM)*ASR(I,J,KM))
      DPERP   = VOL(I,J,KMM1)/AREAKM
C************************************************************************
C****************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD  =  YPLUS_KM(I,J)
           DENSITY    = 0.25*ROAVGKM
           WSLIP      = 0.25*WAVG
        CALL WALLFUN(I,J,KMM1,KM,DPERP,DPDS_CELL(I,J,KMM1),DENSITY,
     &              TWALLKM,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_KM(I,J) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 355
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************

      RE      = DPERP*ROAVGKM*WAVG/VISLAM16
      RELOG   = 1./ALOG(RE)
C******************************************************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_T(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C******************************************************************************
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLKM = CF*ROAVGKM*WAVG*WAVG/128
           VSTAR   = SQRT(TWALLKM/(ROAVGKM/4))
           PLUSK   = ROUGH*VSTAR*ROAVGKM/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C******************************************************************************
C  SET THE SHEAR STRESS ON THE CASING .
C
      TWALLKM    = CF*ROAVGKM*WAVG*WAVG/128
C
C******************************************************************************
C
  355 CONTINUE
C
C******************************************************************************
C     SET YPLUS FOR USE LATER
      IF(IBOUND.LE.1) THEN
      DO K = KMID+1,KMM1
      VSTAR     = SQRT(TWALLKM/(0.25*ROAVGKM))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVGKM)/VISLAM
      YPLSP     = Y_PLUS(I,J,K)
      Y_PLUS(I,J,K) = AMIN1(YPLS,YPLSP)
      END DO
      END IF
C 
      FMULT   = -TWALLKM*AREAKM/WAVG
      IF(IBOUND.GE.2) FMULT = 0.0 
C
      XFORCEKM = FMULT*VXAVGKM
      RFORCEKM = FMULT*VRAVGKM
      TFORCEKM = FMULT*WTAVGKM
      WVISCKM  =  TFORCEKM*0.25*WTIP(J)*RAVGKM
C
      FORCEX(I,J,KMM1) = FORCEX(I,J,KMM1) + XFORCEKM
      FORCER(I,J,KMM1) = FORCER(I,J,KMM1) + RFORCEKM
      FORCET(I,J,KMM1) = FORCET(I,J,KMM1) + TFORCEKM
      ESOURCE(I,J,KMM1)= ESOURCE(I,J,KMM1)+ WVISCKM
C
C     NOW SUM THE STRESSES ON THE OTHER K = CONSTANT FACES
C
      DO 450 K=2,KMM1
      KM1  = K-1
      XFOR = 0.5*((TXR(I,J,K)+TXR(I,J,KM1))*ASR(I,J,K)
     &    +       (TXX(I,J,K)+TXX(I,J,KM1))*ASX(I,J,K))
      FORCEX(I,J,KM1) = FORCEX(I,J,KM1) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR

      RFOR = 0.5*((TRR(I,J,K)+TRR(I,J,KM1))*ASR(I,J,K)
     &    +       (TRX(I,J,K)+TRX(I,J,KM1))*ASX(I,J,K))
      FORCER(I,J,KM1)  = FORCER(I,J,KM1) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR

      TFOR = 0.5*((TTR(I,J,K)+TTR(I,J,KM1))*ASR(I,J,K)
     &    +       (TTX(I,J,K)+TTX(I,J,KM1))*ASX(I,J,K))
      FORCET(I,J,KM1)  = FORCET(I,J,KM1) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(I,J,KM1))*ASX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(I,J,KM1))*ASR(I,J,K))
C
      VXAVG = VX(I,J,K)+VX(IP1,J,K)+VX(IP1,JM1,K)+VX(I,JM1,K)
      VRAVG = VR(I,J,K)+VR(IP1,J,K)+VR(IP1,JM1,K)+VR(I,JM1,K)
      VTAVG = VT(I,J,K)+VT(IP1,J,K)+VT(IP1,JM1,K)+VT(I,JM1,K)

      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)

      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(I,J,KM1)= ESOURCE(I,J,KM1) - QFLOW - WVISC 
C
  450 CONTINUE
C
  400 CONTINUE
C
  401 CONTINUE
C
C******************************************************************************
C
C     ALL THE VISCOUS FORCES ARE NOW SET. USE THEM TO UPDATE THE
C     GLOBAL BODY FORCE TERMS.
C
C******************************************************************************
C     NOTE THE NEGATIVE SIGN BECAUSE THE SOURCE TERM IS SUBTRACTED IN SUBROUTINE TSTEP.
C 
      RF_VIS1 = 1.0 - RF_VIS
      DO 500 K=1,KMM1
      DO 500 J=2,JM
      DO 500 I=1,IMM1
      XFORCE(I,J,K)   = RF_VIS1*XFORCE(I,J,K)   - RF_VIS*FORCEX(I,J,K)
      TFORCE(I,J,K)   = RF_VIS1*TFORCE(I,J,K)   - RF_VIS*FORCET(I,J,K)
     &                 *RAVG_CELL(J,K)
      RFORCE(I,J,K)   = RF_VIS1*RFORCE(I,J,K)   - RF_VIS*FORCER(I,J,K)
      QSOURCE(I,J,K)  = RF_VIS1*QSOURCE(I,J,K)  - RF_VIS*ESOURCE(I,J,K)
  500 CONTINUE
C
C******************************************************************************
C******************************************************************************
C
      RETURN
      END
C
C******************************************************************************
C
C*************************************************************************************
C
      SUBROUTINE GRADVEL(I,J,K, V, DVDX, DVDR, DVDT)
C
C******************************************************************************************
C     THIS SUBROUTINE CALCULATES THE DERIVATIVES OF THE INPUT FUNCTION  "V " 
C     FOR A CELL WITH CORNER STORAGE USING GAUSS' THEOREM .
C
C     BECAUSE THE AREAS OF THE FACES ARE EVALUATED RELATIVE TO THE RADIAL DIRECTION 
C     AT THE MIDDLE OF THE FACE  AN EXTRA  "V/R"  TERM MUST BE SUBTRACTED FROM THE 
C     RADIAL DERIVATIVE.
C******************************************************************************************
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  V(ID,JD,KD)
C
      IP1 = I+1
      JM1 = J-1
      KP1 = K+1
C  NOTE THAT THE J INDEX OF A CELL IS THE J VALUE ON ITS RIGHT HAND SIDE. 
C  ALL AREAS ARE POSITIVE IN THE DIRECTION OF AN INWARDS NORMAL.    
      RVOL    = 0.25/VOL(I,J,K)
C
      VAVGJ1  = V(I,JM1,K)+ V(IP1,JM1,K)+ V(IP1,JM1,KP1)+ V(I,JM1,KP1)
      VAVGJ2  = V(I,J,K)  + V(IP1,J,K)  + V(IP1,J,KP1)  + V(I,J,KP1)
      VAVGI1  = V(I,J,K)  + V(I,JM1,K)  + V(I,JM1,KP1)  + V(I,J,KP1)
      VAVGI2  = V(IP1,J,K)+ V(IP1,JM1,K)+ V(IP1,JM1,KP1)+ V(IP1,J,KP1)
      VAVGK1  = V(I,J,K)  + V(IP1,J,K)  + V(IP1,JM1,K)  + V(I,JM1,K)
      VAVGK2  = V(I,J,KP1)+ V(IP1,J,KP1)+ V(IP1,JM1,KP1)+ V(I,JM1,KP1)
      VAVGCELL= 0.125*(VAVGJ1 + VAVGJ2)
C
      DVDX = -(VAVGJ2*AQX(I,J,K)   - VAVGJ1*AQX(I,JM1,K)
     &       + VAVGI1*ABX(I,J,K)   - VAVGI2*ABX(IP1,J,K)
     &       + VAVGK1*ASX(I,J,K)   - VAVGK2*ASX(I,J,KP1))*RVOL
C
      DVDR = -(VAVGJ2*AQR(I,J,K)   - VAVGJ1*AQR(I,JM1,K)
     &       + VAVGI1*ABR(I,J,K)   - VAVGI2*ABR(IP1,J,K)
     &       + VAVGK1*ASR(I,J,K)   - VAVGK2*ASR(I,J,KP1))*RVOL
     &       - VAVGCELL/RAVG_CELL(J,K)
C
      DVDT = -(VAVGI1  - VAVGI2)*ABT(J,K)*RVOL
C
      RETURN
      END
C*************************************************************************************
C*************************************************************************************
C
C*************************************************************************************
C
      SUBROUTINE GRADCELL(I,J,K, V, DVDX, DVDR, DVDT)
C
C******************************************************************************************
C     THIS SUBROUTINE CALCULATES THE DERIVATIVES OF THE INPUT FUNCTION  "V" 
C     FOR A CELL WITH CELL CENTRE STORAGE USING GAUSS' THEOREM, WHERE  "V" IS THE AVERAGE
C     VALUE FOR THE CELL .
C
C     BECAUSE THE AREAS OF THE FACES ARE EVALUATED RELATIVE TO THE RADIAL DIRECTION 
C     AT THE MIDDLE OF THE FACE  AN EXTRA  "V/R"  TERM MUST BE SUBTRACTED FROM THE 
C     RADIAL DERIVATIVE.
C******************************************************************************************
C******************************************************************************************
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  V(ID,JD,KD)
C
      IP1 = I+1
      JP1 = J+1
      KP1 = K+1
      IM1 = I-1
      JM1 = J-1
      KM1 = K-1
      IF(I.EQ.1)    IM1 = IMM1
      IF(I.EQ.IMM1) IP1 = 1
      IF(K.EQ.1)    KM1 = 1
      IF(K.EQ.KMM1) KP1 = KMM1
C  NOTE THAT THE J INDEX OF A CELL IS THE J VALUE ON ITS RIGHT HAND SIDE. 
C  ALL AREAS ARE POSITIVE IN THE DIRECTION OF AN INWARDS NORMAL.  
C   1/RVOL IS x 0.5  TO ALLOW FOR AVERAGING.
C  
      RVOL    = 0.5/VOL(I,J,K)
C
      VAVGJ1  = V(I,J,K)  + V(I,JM1,K)
      VAVGJ2  = V(I,J,K)  + V(I,JP1,K)
      VAVGI1  = V(I,J,K)  + V(IM1,J,K)
      VAVGI2  = V(I,J,K)  + V(IP1,J,K)
      VAVGK1  = V(I,J,K)  + V(I,J,KM1)
      VAVGK2  = V(I,J,K)  + V(I,J,KP1)
C   SET ZERO VALUE ON BLADE SURFACES
      IF(IND(J).EQ.1) THEN
           IF(I.EQ.1)    VAVGI1 = 0.0
           IF(I.EQ.IMM1) VAVGI2 = 0.0
      END IF
C   SET ZERO VALUE ON HUB AND CASING
C  JDD CHANGED TO ALLOW FOR NO SHEAR ON HUB OR CASING. JAN/14                
      IF(K.EQ.1.AND.(IBOUND.EQ.0.OR.IBOUND.EQ.2)) VAVGK1 = 0.0
      IF(K.EQ.KMM1.AND.IBOUND.LE.1)               VAVGK2 = 0.0
C
      DVDX = -(VAVGJ2*AQX(I,J,K)   - VAVGJ1*AQX(I,JM1,K)
     &       + VAVGI1*ABX(I,J,K)   - VAVGI2*ABX(IP1,J,K)
     &       + VAVGK1*ASX(I,J,K)   - VAVGK2*ASX(I,J,KP1))*RVOL
C
      DVDR = -(VAVGJ2*AQR(I,J,K)   - VAVGJ1*AQR(I,JM1,K)
     &       + VAVGI1*ABR(I,J,K)   - VAVGI2*ABR(IP1,J,K)
     &       + VAVGK1*ASR(I,J,K)   - VAVGK2*ASR(I,J,KP1))*RVOL
     &       - V(I,J,K)/RAVG_CELL(J,K)
C
      DVDT = -(VAVGI1  - VAVGI2)*ABT(J,K)*RVOL
C
C*************************************************************************************
C*************************************************************************************
C
      RETURN
      END
C
C*************************************************************************************
C*************************************************************************************
C*************************************************************************************
C
      SUBROUTINE SPAL_LOSS
C
C*****************************************************************************************
C     THIS SUBROUTINE COMPUTES A VISCOUS BODY FORCE BASED ON WALL FUNCTIONS FOR THE
C     SURFACE SHEAR STRESS AND THE SPALART - ALLMARAS TURBULENCE MODEL.          
C******************************************************************************************
C
      INCLUDE  'commall-open-18.3'
C
      COMMON/BKSTRESS/  TXX(ID,JD,KD),TXR(ID,JD,KD),
     &                  TXT(ID,JD,KD),TRX(ID,JD,KD),TRR(ID,JD,KD),
     &                  TRT(ID,JD,KD),TTX(ID,JD,KD),TTR(ID,JD,KD),
     &                  TTT(ID,JD,KD),QXX(ID,JD,KD),QRR(ID,JD,KD),
     &                  QTT(ID,JD,KD)
C
      DIMENSION  FORCEX(ID,JD,KD),FORCER(ID,JD,KD),FORCET(ID,JD,KD),
     &           ESOURCE(ID,JD,KD)
C
      RF_VIS1 = 1.0 - RF_VIS
      RF_VIS4 = 0.25*RF_VIS
      NDAMP   = 100000
C
C      CALCULATE THE VISCOSITY OVER THE FIRST QUARTER OF THE STEPS.
C      THEN HOLD IT CONSTANT FOR THE REMAINDER OF THE STEPS.
C
      IF(NSTEP.EQ.1.OR.NSTEP.LT.NSTEPS_MAX/4) THEN
C
      IF(REYNO.GT.100.) THEN
           J1=JLE(1)
           J2=JTE(1)
           IF(JLE(1).GT.JM) J1=1
           IF(JTE(1).GT.JM) J2=JM
           XCHORD = SMERID(J2,KMID)-SMERID(J1,KMID)
           ROW2   = SQRT(ROVX(IMID,J2,KMID)*ROVX(IMID,J2,KMID)
     &            + ROVR(IMID,J2,KMID)*ROVR(IMID,J2,KMID)
     &            + ROWT(IMID,J2,KMID)*ROWT(IMID,J2,KMID) )
           VISLAM = XCHORD*ROW2/REYNO
      END IF
C
      IF(REYNO.GT.0.0.AND.REYNO.LT.99.99) THEN
            VISLAM = REYNO/100000.
      END IF
      IF(REYNO.LT.0.0) VISLAM = -REYNO*1.0E-5
C
C     ENERGY EQN ADDITION
C
      THCOND  = CP*VISLAM/PRANDTL
      FTCOND = THCOND/VISLAM
C
      ENDIF
C
C    END OF SETTING THE VISCOSITY, ETC OVER THE FIRST QUARTER OF THE STEPS.
C
C********************************************************************************
C     INITIALLY SET THE VISCOUS FORCES AND ENERGY SOURCES TO ZERO
C
      DO 100 K=1,KMM1
      DO 100 J=2,JM
      DO 100 I=1,IMM1
      FORCEX(I,J,K)  = 0.0
      FORCER(I,J,K)  = 0.0
      FORCET(I,J,K)  = 0.0
      ESOURCE(I,J,K) = 0.0
  100 CONTINUE
C
C*******************************************************************************************
C*******************************************************************************************
C     EVALUATE THE TURBULENT VISCOSITY FOR EACH CELL in THE DO 150 LOOP.
C     ONLY USE THE PRESSURE GRADIENT TERM IF FAC_PGRAD IS GREATER THAN ZERO.
C
      IF(FAC_PGRAD.GT.0.001) THEN
          CALL SET_PWALLGRAD
      END IF
C
      DO 1000 K=1,KMM1
      KP1 = K+1
      DO 1000 J=2,JM
      JM1 = J-1
      NRW = NROW(J)
      JTEDGE = JTE(NRW)
      JLEDGE = JLE(NRW)
      DO 1000 I=1,IMM1
      IP1 = I+1
C
C    AVERAGE THE CONDITIONS FOR THE CELLS. Note these are true averages.
C
      WTAVG = 0.125*(WT(I,J,K)+WT(IP1,J,K)+WT(IP1,J,KP1)+WT(I,J,KP1)
     &      + WT(I,JM1,K)+WT(IP1,JM1,K)+WT(IP1,JM1,KP1)+WT(I,JM1,KP1))
      ROAVG = ROAVG_CELL(I,J,K)
      TSAVG = 0.125*(T_STATIC(I,J,K)+T_STATIC(IP1,J,K)
     &        +T_STATIC(IP1,J,KP1)+T_STATIC(I,J,KP1)+T_STATIC(I,JM1,K)
     &        +T_STATIC(IP1,JM1,K)+T_STATIC(IP1,JM1,KP1)
     &        +T_STATIC(I,JM1,KP1))
      RAVG  = RAVG_CELL(J,K)
C
C     CALCULATE THE DERIVATIVES OF THE VELOCITY COMPONENTS AND TEMPERATURE.
C     NOTE THE VORTICITY IS BASED ON THE RELATIVE VELOCITY,  WT.
C     THIS WAS CHANGED FROM THE ABSOLUTE VELOCITY BY JDD AUGUST 2017.
      CALL GRADVEL(I,J,K,VX,DVXDX,DVXDR,DVXDT)
      CALL GRADVEL(I,J,K,VR,DVRDX,DVRDR,DVRDT)
      CALL GRADVEL(I,J,K,WT,DWTDX,DWTDR,DWTDT)
C
      CALL GRADVEL(I,J,K,T_STATIC,DTEMPDX,DTEMPDR,DTEMPDT)
C
C     CALCULATE THE RATES OF STRAIN FOR USE WHEN CALCULATING THE VISCOUS STRESSES.
C     THESE ARE NOT YET THE STRESSES AS THE FINAL VISCOSITY IS NOT YET KNOWN.
C     
      TXR(I,J,K) = DVXDR + DVRDX
      TXX(I,J,K) = DVXDX
      TXT(I,J,K) = DWTDX + DVXDT
      TRR(I,J,K) = DVRDR
      TRT(I,J,K) = DWTDR + DVRDT - WTAVG/RAVG
      TTT(I,J,K) = DWTDT
      QXX(I,J,K) = DTEMPDX
      QRR(I,J,K) = DTEMPDR
      QTT(I,J,K) = DTEMPDT
C
C        USE THE VORTICITY TO DETERMINE THE TURBULENT VISCOSITY FOR EACH CELL.
C
      VORTX   = DVRDT - DWTDR - WTAVG/RAVG 
      VORTR   = DWTDX - DVXDT
      VORTT   = DVXDR - DVRDX
      VORT    = SQRT(VORTX*VORTX + VORTR*VORTR + VORTT*VORTT)
      IF(J.EQ.JSTART(NRW)) VORT = 0.0
C
C   SET A LIMIT TO THE VORTICITY
      IF(VORT.GT.VORT_MAX) VORT = VORT_MAX
C
      IF(INDMIX(J).EQ.1) THEN
            VORT = ABS_VORT(I,J-1,K)
      END IF
C
      ABS_VORT(I,J,K) = VORT
C******************************************************************************
C******************************************************************************
C
C     CALCULATE THE FACTOR TO INCREASE THE SOURCE TERM DUE TO STREAMWISE VORTICITY.
C     AS USED BY Lee, Wilson & Vahdati.
C
      IF(FAC_VORT.GT.0.001.OR.FAC_PGRAD.GT.0.001) THEN
C
      WTAVG = (WT(I,J,K)+WT(IP1,J,K)+WT(IP1,J,KP1)+WT(I,J,KP1)
     &      + WT(I,JM1,K)+WT(IP1,JM1,K)+WT(IP1,JM1,KP1)+WT(I,JM1,KP1))
      VXAVG = (VX(I,J,K)+VX(IP1,J,K)+VX(IP1,J,KP1)+VX(I,J,KP1)
     &      + VX(I,JM1,K)+VX(IP1,JM1,K)+VX(IP1,JM1,KP1)+VX(I,JM1,KP1))
      VRAVG =(VR(I,J,K)+VR(IP1,J,K)+VR(IP1,J,KP1)+VR(I,J,KP1)
     &      + VR(I,JM1,K)+VR(IP1,JM1,K)+VR(IP1,JM1,KP1)+VR(I,JM1,KP1)) 
      WREL  = SQRT(VXAVG*VXAVG + VRAVG*VRAVG + WTAVG*WTAVG)
C
      END IF
C
      IF(FAC_VORT.GT.0.001) THEN
C  
      HELICITY  = (VXAVG*VORTX + VRAVG*VORTR + WTAVG*VORTT)
     &            /(WREL*VORT)
      HELFAC    = HELICITY
C
C    THE NEXT COMMENTED OUT LINE IS THE VALUE OF FAC_VORT USED BY Lee, Wislon & Vahdati.
C    VORTICITY_FAC  = FAC_VORT*(1.0 + 0.9191*TANH(3.0*HELFAC*HELFAC))
       VORTICITY_FAC  = (1.0 + FAC_VORT*TANH(3.0*HELFAC*HELFAC))
C    
      END IF
C
C******************************************************************************
C******************************************************************************
C
C     NOW SET THE PRESSURE GRADIENT TERM AS USED BY Lee, Wislon & Vahdati .
C
      IF(FAC_PGRAD.GT.0.001) THEN
C
C  SET A LIMITING RELATIVE VELOCITY = 0.1 x THE LOCAL MID POINT VALUE.
      WLIM = 0.1*SQRT(VX(IMID,J,K)*VX(IMID,J,K)
     &              + VR(IMID,J,K)*VR(IMID,J,K)
     &              + WT(IMID,J,K)*WT(IMID,J,K))
      IF(WREL.LT.WLIM) WREL = WLIM
      WREL  = 0.125*WREL
      PGRAD = DPDS_CELL(I,J,K)
C
C   SET SPAREVAR TO PLOT OUT THE PRESSURE GRADIENT FIELD.
C      SPAREVAR(I,J,K) = PGRAD/1.0E6
C
C    SET THE PRESSURE GRADIENT TERM TO ZERO FOR FAVOURABLE PRESSURE GRADIENTS.
      IF(PGRAD.LT.0.0) PGRAD = 0.0
C
      PTERM = PGRAD*VISLAM/(ROAVG*ROAVG*WREL*WREL*WREL)
C
C   NOTE THE PRESSURE GRADIENT TERM IS SCALED BY REYNOLDS NUMBER SO THAT IT
C   BECOMES INDEPENDENT OF VISCOSITY. THIS IS DIFFERENT TO Lee, Wislon & Vahdati .
      PTERM = PTERM*REYNOLDS(NRW)
C
C  THE NEXT COOMMENTED OUT LINE IS THE VAUE OF FAC_PGRAD USED BY Lee, Wilson & Vahdati.
C      PGRAD_FAC = FAC_PGRAD*(1.0 + 0.6565*TANH(PTERM*PTERM) )
      PGRAD_FAC = (1.0 + FAC_PGRAD*TANH(PTERM*PTERM) )
      SPAREVAR(I,J,K) = PGRAD_FAC
C
      END IF     
       
C******************************************************************************
C     SET A TURBULENT VISCOSITY USING THE MIXING LENGTH MODEL
C
      TURBVIS_WALL = ROAVG*XLENGTH(I,J,K)*VORT*FMIXUP
C
C******************************************************************************
C     SET THE LAMINAR VISCOSITY IF IT IS TO VARY WITH TEMPERATURE
      IF(REYNO.LT.0.0001) THEN
            VISLAM = (ABS(REYNO)/100000.0) * (TSAVG/288.0)**0.62 
      END IF
C
C******************************************************************************
C    SAVE THE VISCOSITY IT IS ONLY USED TO CALCULATE AND WRITE OUT THE REYNOLDS NUMBER.
      IF(I.EQ.IMID.AND.K.EQ.KMID.AND.J.EQ.JTEDGE) VISCOSY(NRW)  = VISLAM
C******************************************************************************
C     CALCULATE THE REYNOLDS NUMBER AT MID PASSAGE AT THE TRAILING EDGE.
C
      IF(I.EQ.IMID.AND.K.EQ.KMID.AND.J.EQ.JTEDGE) THEN
           WABS       = SQRT(VX(I,J,K)*VX(I,J,K) + VR(I,J,K)*VR(I,J,K)
     &                     + WT(I,J,K)*WT(I,J,K) )
           ROVEXIT(NRW)  = RO(I,J,K)*WABS
           VISCOSY(NRW)  = VISLAM
           REYNOLDS(NRW) = CHORD(NRW)*ROVEXIT(NRW)/VISLAM
      END IF
C******************************************************************************
C     ADD THE FREE STREAM TURBULENT VISCOSITY
C
      VISTOT  = VISLAM*(1 + FST_RAT(I,J,K) ) + TURBVIS_WALL 
C
C******************************************************************************
C    SET A LIMIT TO THE TURBULENT VISCOSITY
C
      VISLIM = TURBVIS_LIM*VISLAM
      IF(VISTOT.GT.VISLIM) VISTOT = VISLIM
      VISKLAM  = VISLAM/ROAVG
C
C************************************************************
C     SAVE THE MIXING LENGTH TURBULENT DYNAMIC VISCOSITY. 
C
      DYN_VIS_ML(I,J,K) = (VISTOT - VISLAM) 
C
C******************************************************************************
C******************************************************************************
C     NOW SET THE INITIAL GUESS OF VISCOSITY FOR THE S-A MODEL WHEN NSTEP = 1.
C
      IF(IF_RESTART.EQ.0.AND.NSTEP.EQ.1)THEN
C
           IF(FSTURB(1).LT.1.0) THEN 
           TRANSVISIN  = 4.35*VISLAM*FSTURB(1)**0.25
           ELSE
           TRANSVISIN  = VISLAM*(3.5 + FSTURB(1)*0.85)
           END IF
C
          TRANS_DYN_VIS(I,J,K) = TRANSVISIN
C
      END IF 
C
      IF(NSTEP.EQ.1) TRANS_KVIS(I,J,K) = TRANS_DYN_VIS(I,J,K)/ROAVG
C
C******************************************************************************
C  NOW SET THE FINAL VALUE OF TURBULENT VISCOSITY
C
      TRANSKVIS        = TRANS_KVIS(I,J,K)
      XI               = TRANSKVIS/VISKLAM
C
      IF(XI.LT.1.1) XI = 1.1
      IF(XI.GT.TURBVIS_LIM) XI = TURBVIS_LIM
C
      XI3              = XI*XI*XI
      FV1              = XI3/(XI3 + 357.9)
      TRANSVISLIM      = VISLIM/ROAVG/FV1
      IF(TRANSKVIS.GT.TRANSVISLIM) TRANSKVIS = TRANSVISLIM
C
C*******************************************************************************************
C*******************************************************************************************
C   THIS IS THE MAIN RESULT OF THE SA MODEL. VISTRU IS THE TRUE TURBULENT DYNAMIC VISCOSITY.
C
      VISTRU = TRANSKVIS*FV1*ROAVG
C
C*******************************************************************************************
C
      FV2              = 1.0 - XI/(1.0 + XI*FV1)
C 
      IF(XI.LT.2.0.AND.FV2.GT.0.0) FV2 = 0.0
C
C      TRANS_VORT(I,J,K) = ABS_VORT(I,J,K) 
C     &                  + 5.9488*FV2*TRANSKVIS/DWALLSQ(I,J,K)
C     
C   JDD CHANGED THIS Jan/14.  USE THE TRUE VORTICITY RATHER THAN THE TRANSFORMED ONE.
C   THIS IS DIFFERENT TO THE ORIGINAL SA MODEL WHICH IS COMMENTED OUT ABOVE.
C
      TRANS_VORT(I,J,K) = ABS_VORT(I,J,K)
C      
      VISTOT           = VISTRU + VISLAM
C
      VISLIM = TURBVIS_LIM*VISLAM
      IF(VISTOT.GT.VISLIM) VISTOT = VISLIM
C
C     SAVE THE TURBULENT/LAMINAR VISCOSITY RATIO AS "VISC_RAT" AND THE LAMINAR VISCOSITY AS "VISC_LAM".
C     THESE AND  VISTOT ARE THE DYNAMIC VISCOSITIES. 
C
       VISC_RAT(I,J,K)   = VISTOT/VISLAM
       VISC_LAM(I,J,K)   = VISLAM
       VISC_TURB(I,J,K)  = VISTOT - VISLAM
C
 1000 CONTINUE
C
C******************************************************************************
C******************************************************************************
C   FORM THE GRAD OF TURBULENT VISCOSITY FOR THE SOURCE TERMS.
C
      DO 1500,K=1,KMM1
      DO 1500 J=2,JMM1
      NRW = NROW(J)
      DO 1500,I=1,IMM1
      CALL GRADCELL(I,J,K,TRANS_KVIS,D_TVISDX,D_TVISDR,D_TVISDT)
      GRAD_TVIS(I,J,K) = D_TVISDX*D_TVISDX + D_TVISDR*D_TVISDR
     &                 + D_TVISDT*D_TVISDT
      VISCEFF          = TRANS_KVIS(I,J,K) + VISC_LAM(I,J,K)
      TVISGRADX(I,J,K) = VISCEFF*D_TVISDX
      TVISGRADR(I,J,K) = VISCEFF*D_TVISDR
      TVISGRADT(I,J,K) = VISCEFF*D_TVISDT
      IF(J.EQ.JSTART(NRW).OR.J.EQ.JSTART(NRW)-1) THEN
      GRAD_TVIS(I,J,K) = 0.0
      TVISGRADX(I,J,K) = 0.0
      TVISGRADR(I,J,K) = 0.0
      TVISGRADT(I,J,K) = 0.0
      END IF
 1500 CONTINUE
C******************************************************************************
C
C   
C    CALCULATE THE SOURCE TERMS
C
      CW2    = 0.3
      SIGMA  = 0.6667
      CB2    = 0.622
      CB1    = 0.1355
      CW1    = 0.80607 + 2.433
      FST0   = FAC_ST0*CB1
      FST1   = FAC_ST1/SIGMA
      FST2   = FAC_ST2*CB2/SIGMA
      FST3   = FAC_ST3*CW1
      FSTMIX = FAC_STMIX/NDAMP
C
      DO 1600 K=1,KMM1
      DO 1600 J=2,JM
      DO 1600 I=1,IMM1
C
C   JDD JAN/14 SET FYPLUS TO REDUCE THE SOURCE TERM IN THE SUBLAYER AND BUFFER REGIONS.
C
      IF(Y_PLUS(I,J,K).LE.YPLAM) FYPLUS = 0.0
      IF(Y_PLUS(I,J,K).GT.YPLAM)  THEN
               XFAC = (Y_PLUS(I,J,K) - YPLAM)/(YPTURB- YPLAM)
               IF(XFAC.GT.1.0) XFAC = 1.0
               FYPLUS = XFAC*XFAC*(3.0 - 2.0*XFAC)
      END IF
C
C  END JAN/14 ADDITION.

      CALL GRADCELL(I,J,K,TVISGRADX,TVISDDX,DUMY,DUMMY2) 
      CALL GRADCELL(I,J,K,TVISGRADR,DUMY,TVISDDR,DUMMY2) 
      CALL GRADCELL(I,J,K,TVISGRADT,DUMY,DUMMY2,TVISDDT) 
C
      ROAVG   = ROAVG_CELL(I,J,K)
C
      TVIS    = TRANS_KVIS(I,J,K)

      VISKLAM = VISC_LAM(I,J,K)/ROAVG
C
      ST0     = FST0*TRANS_VORT(I,J,K)*TVIS
      IF(FAC_VORT.GT.0.001)  ST0 = ST0*VORTICITY_FAC
      IF(FAC_PGRAD.GT.0.001) ST0 = ST0*PGRAD_FAC
C
C    JDD ADDITION JAN/14.  REDUCE THE SOURCE TERM IN THE SUBLAYER AND BUFFER LAYER.
      ST0     = ST0*FYPLUS
C  END JAN/14 ADDITION
C
      ST1     = FST1*(TVISDDX + TVISDDR + TVISDDT)
C
      ST2     = FST2*GRAD_TVIS(I,J,K)
C
      TAU     = 5.9488*TVIS/TRANS_VORT(I,J,K)/DWALLSQ(I,J,K)
C
C   SET SOME LIMITS ON TAU
C
      IF(TAU.GT.1.2)   TAU = 1.2
      IF(TAU.LT.0.001) TAU = 0.001
C
      G       = TAU + CW2*(TAU*TAU*TAU*TAU*TAU*TAU - TAU)
      G6CW6   =  G*G*G*G*G*G/64.0
      FW      =  G*((1 + 0.015625)/(1 + G6CW6))**0.166667
C
      ST3     =  FST3*FW*TVIS*TVIS/DWALLSQ(I,J,K)
C
C     SET AN EXTRA SOURCE TERM WHICH DRIVES THE TURBULENT VISCOSITY TOWARDS
C     THE MIXING LENGTH VALUE
C
      ST_MIXL = FSTMIX*(DYN_VIS_ML(I,J,K)-VISC_TURB(I,J,K))/STEP(I,J,K)
C
C******************************************************************************
C  SET THE COMBINED SOURCE TERM , TSOURCE  .
C
      T_SOURCE(I,J,K) = RF_VIS1*T_SOURCE(I,J,K) + 
     & RF_VIS*(ROAVG*(ST0 + ST1 + ST2 - ST3 )*VOL(I,J,K) + ST_MIXL)
C
C******************************************************************************
C  SET THE SOURCE TERM TO ZERO FOR THE DUMMY CELL AT THE MIXING PLANE.
C
      IF(J.GT.2.AND.INDMIX(J-1).EQ.1) T_SOURCE(I,J,K) = 0.0
      IF(J.GT.2) THEN
          IF(INDMIX(J-2).EQ.1)        T_SOURCE(I,J,K) = 0.0
      END IF
C
 1600 CONTINUE
C
C
C    END OF SETTING THE TURBULENT VISCOSITY
C****************************************************************************
C****************************************************************************
C     NEXT  CHECK FOR TRANSITION AND REDUCE THE SOURCE TERM IN LAMINAR REGIONS
C
      IEND1 = IMID/2
      IENDM = IMM1 - IEND1    
      KEND1 = KMID/2
      KENDM = KMM1 - KEND1
C
      DO 160 J = 2,JM
C
      NRW    = NROW(J)
      J1     = JSTART(NRW)
      JROW   = J - J1 + 1
      JTRHUB = JTRAN_K1(NRW)
      JTRTIP = JTRAN_KM(NRW)
      JTRLOW = JTRAN_I1(NRW)
      JTRUP  = JTRAN_IM(NRW)
      JLEROW = JLE(NRW)
C
C   FIRST ON THE HUB
C
C   Q3D
      IF(KM.EQ.2)  GO TO 171
C   END Q3D
      DO 170 I = 1,IMM1
      IF(JTRHUB.LE.1.AND.FTRANS.LT.0.1) GO TO 38
      IF(JROW.LT.JTRHUB) GO TO 34
      RMAX = 0.0
      DO 37 K=1,KEND1
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 37
      RMAX = RATVIS
   37 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 38
   34 CONTINUE
C
C      TVEDGE = VISC_TURB(I,KEND1,J)
C      Changed by WS/LX   17 May, 2017 
      TVEDGE = VISC_TURB(I,J,KEND1)
C
      DO 39 K = 1, KEND1
      FK = FLOAT(KEND1 + 1 -K)/KEND1
      FK  = FK*FK*(3 - 2*FK)
      TVSPEC = FK*TVEDGE
   39 T_SOURCE(I,J,K) = T_SOURCE(I,J,K)
     &    -(VISC_TURB(I,J,K)-TVSPEC)/STEP(I,J,K)/NDAMP
C
   38 CONTINUE

C  NEXT ON THE CASING

      IF(JTRTIP.LE.1.AND.FTRANS.LT.0.1) GO TO 49
      IF(JROW.LT.JTRTIP) GO TO 46
      RMAX = 0.0
      DO 47 K = KENDM, KMM1
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 47
      RMAX = RATVIS
   47 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO 49
   46 CONTINUE
C
C      TVEDGE = VISC_TURB(I,KENDM,J)
C      Changed by WS/LX  17 May, 2017  
      TVEDGE = VISC_TURB(I,J,KENDM)
C
      DO 48 K = KENDM, KMM1
      FK = FLOAT(K + 1 -KENDM)/(KM - KENDM)
      FK  = FK*FK*(3 - 2*FK)
      TVSPEC  = FK*TVEDGE
   48 T_SOURCE(I,J,K) = T_SOURCE(I,J,K)
     &    -(VISC_TURB(I,J,K)-TVSPEC)/STEP(I,J,K)/NDAMP
C
   49 CONTINUE
C
  170 CONTINUE
C
  171 CONTINUE
C
C    NEXT ON THE LOWER - I=1, BLADE SURFACE

      DO 180  K=1,KMM1 
C
      IF(JTRLOW.LE.1.AND.FTRANS.LT.0.1) GO TO 71 
      IF(J.LE.JLEROW) GO TO 71
      IF(JROW.LT.JTRLOW.OR.J.LT.JLEROW) GO TO 70
      RMAX = 0.0
      DO 56 I=1,IEND1
      RATVIS = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 56
      RMAX = RATVIS
   56 CONTINUE
      IF(RMAX.GT.FTRANS)  GO TO 71
   70 CONTINUE
      TVEDGE = VISC_TURB(IEND1,J,K)
      DO 72 I = 1,IEND1
      FI = FLOAT(IEND1 + 1 -I)/IEND1
      FI  = FI*FI*(3 - 2*FI)
      TVSPEC = FI*TVEDGE
   72 T_SOURCE(I,J,K) = T_SOURCE(I,J,K)
     &    -(VISC_TURB(I,J,K)-TVSPEC)/STEP(I,J,K)/NDAMP
C
   71 CONTINUE

C    NEXT ON THE UPPER, I = IM, BLADE SURFACE

      IF(JTRUP.LE.1.AND.FTRANS.LT.0.1) GO TO 77 
      IF(J.LE.JLEROW) GO TO 77
      IF(JROW.LT.JTRUP) GO TO 69
      RMAX = 0.0
      DO 73 I = IENDM,IMM1
      RATVIS  = VISC_RAT(I,J,K)
      IF(RATVIS.LT.RMAX) GO TO 73
      RMAX = RATVIS
   73 CONTINUE
      IF(RMAX.GT.FTRANS) GO TO  77
   69 CONTINUE
      TVEDGE = VISC_TURB(IENDM,J,K)
      DO 74 I  = IENDM,IMM1
      FI = FLOAT(I + 1 -IENDM)/(IM - IENDM)
      FI  = FI*FI*(3 - 2*FI)
      TVSPEC = FI*TVEDGE
   74 T_SOURCE(I,J,K) = T_SOURCE(I,J,K)
     &    -(VISC_TURB(I,J,K)-TVSPEC)/STEP(I,J,K)/NDAMP

   77 CONTINUE
C
  180 CONTINUE
C
  160 CONTINUE
C
C********************************************************************************
C********************************************************************************
C
C     NOW THE FINAL VISCOSITY IS FIXED SET THE STRESSES AND HEAT FLOWS IN THE ELEMENT
C     JDD WARNING the div V term in the normal stresses is not yet included.
C
      DO 190 K = 1,KMM1
      DO 190 J = 2,JM
      DO 190 I = 1,IMM1
C
      VISLAM  = VISC_LAM(I,J,K)
      VISTOT  = VISC_RAT(I,J,K)*VISLAM
      THCOND  = CP*VISTOT/PRANDTL
      VISTOT2 = 2.0*VISTOT
C     
      TXR(I,J,K) = VISTOT*TXR(I,J,K)
      TXX(I,J,K) = VISTOT2*TXX(I,J,K)
      TXT(I,J,K) = VISTOT*TXT(I,J,K)
      TRX(I,J,K) = TXR(I,J,K)
      TRR(I,J,K) = VISTOT2*TRR(I,J,K)
      TRT(I,J,K) = VISTOT*TRT(I,J,K)
      TTX(I,J,K) = TXT(I,J,K)
      TTT(I,J,K) = VISTOT2*TTT(I,J,K)
      TTR(I,J,K) = TRT(I,J,K)
      QXX(I,J,K) = -THCOND*QXX(I,J,K)
      QRR(I,J,K) = -THCOND*QRR(I,J,K)
      QTT(I,J,K) = -THCOND*QTT(I,J,K)
C
  190 CONTINUE
C
C
C***************************************************************************************
C****************************************************************************************
C     CALCULATE THE WALL SHEAR STRESSES ON ALL THE SOLID SURFACES.
C
      VISLAM16  = VISLAM*16.0
C
C******************FIRST FOR BLADE SURFACES.****************************
C******************************************************************************
C******************************************************************************
C    EVALUATE AND SMOOTH THE PRESSURE GRADIENTS IF YPLUSWALL IS < -10.0
C
      IF(YPLUSWALL.LT.-10.0) CALL SET_PWALLGRAD
C
C***************************************************************************************
C****************************************************************************************
      DO 200 K=1,KMM1
      KP1 = K+1
      DO 200 J=2,JM
C
C   SKIP IF NOT ON A BLADE SURFACE
      IF(IND(J).EQ.0)  GO TO 220
C   
      JM1 = J-1
      NRW = NROW(J)
C   
C   ALSO SKIP FOR CELLS ABOVE THE TIP GAP
      IF( (KTIPS(NRW).GT.0).AND. 
     &    (K.GE.KTIPS(NRW)).AND.(K.LT.KTIPE(NRW)) ) 
     &    GO TO 220
C
C******************FIRST FOR THE I = 1  BLADE SURFACE.*******************************
C                  CALCULATE THE WALL SHEAR STRESS.
C
      VXAVG1 = VX(2,J,K)+VX(2,JM1,K)+VX(2,JM1,KP1)+VX(2,J,KP1)
      VRAVG1 = VR(2,J,K)+VR(2,JM1,K)+VR(2,JM1,KP1)+VR(2,J,KP1)
      WTAVG1 = WT(2,J,K)+WT(2,JM1,K)+WT(2,JM1,KP1)+WT(2,J,KP1)
      ROAVG1 = RO(2,J,K)+RO(2,JM1,K)+RO(2,JM1,KP1)+RO(2,J,KP1)
      WAVG   = SQRT(WTAVG1*WTAVG1 + VXAVG1*VXAVG1 + VRAVG1*VRAVG1)
      AREA1  = SQRT(ABX(1,J,K)*ABX(1,J,K) + ABR(1,J,K)*ABR(1,J,K)
     &              +ABT(J,K)*ABT(J,K)) 
      DPERP  = VOL(1,J,K)/AREA1
C
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN
C
         YPLUS_OLD =  YPLUS_I1(J,K)
         DENSITY   = 0.25*ROAVG1
         WSLIP     = 0.25*WAVG
         CALL WALLFUN(1,J,K,1,DPERP,DPDS_CELL(1,J,K),DENSITY,
     &                TWALLI1,YPLUS_OLD,WSLIP,YPLUS_NEW)
          YPLUS_I1(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 365
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************

      RE     = DPERP*ROAVG1*WAVG/VISLAM16
      RELOG  = 1./ALOG(RE)
C
C*******************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_L(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLI1    = CF*ROAVG1*WAVG*WAVG/128
C
C          VSTAR   = SQRT(TWALLI1/ROAVG1/4)  
C      Changed by WS/LX  17 May, 2017 
           VSTAR   = SQRT(TWALLI1/(0.25*ROAVG1))
C
           PLUSK   = ROUGH*VSTAR*ROAVG1/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C*********************************************************************************
C     SET THE SHEAR STRESS ON THE  I = 1 BLADE SURFACE
C
      TWALLI1 = CF*ROAVG1*WAVG*WAVG/128
C
C*********************************************************************************
C
  365 CONTINUE
C
C*********************************************************************************
C    CALCULATE YPLUS FOR USE IN SETTING  FYPLUS.
C
      VSTAR   = SQRT(TWALLI1/(0.25*ROAVG1))
      DO I = 1,IMID
C
C      YPLS = RO(I,J,K)*VSTAR*SQRT(DWALLSQ(I,J,K))/VISLAM 
C      Changed by WS/LX  17 May, 2017 
      YPLS = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVG1)/VISLAM
C
      Y_PLUS(I,J,K) = AMIN1(1000.0,YPLS)
      END DO
C        
      FMULT   = -TWALLI1*AREA1/WAVG
C
      XFORCE1 =  FMULT*VXAVG1
      RFORCE1 =  FMULT*VRAVG1
      TFORCE1 =  FMULT*WTAVG1
      WVISC1  =  TFORCE1*WRAD(J)*RAVG_CELL(J,K)
C
      FORCEX(1,J,K)    = FORCEX(1,J,K)    + XFORCE1
      FORCER(1,J,K)    = FORCER(1,J,K)    + RFORCE1
      FORCET(1,J,K)    = FORCET(1,J,K)    + TFORCE1
      ESOURCE(1,J,K)   = ESOURCE(1,J,K)   + WVISC1
C
C***********************************************************************************
C***************************NOW FOR THE I = IM BLADE SURFACE************************
C                     CALCULATE THE WALL SHEAR STRESS.
C
      VXAVGIM = VX(IMM1,J,K)+VX(IMM1,JM1,K)
     &        + VX(IMM1,JM1,KP1)+VX(IMM1,J,KP1)
      VRAVGIM = VR(IMM1,J,K)+VR(IMM1,JM1,K)
     &        + VR(IMM1,JM1,KP1)+VR(IMM1,J,KP1)
      WTAVGIM = WT(IMM1,J,K)+WT(IMM1,JM1,K)
     &        + WT(IMM1,JM1,KP1)+WT(IMM1,J,KP1)
      ROAVGIM = RO(IMM1,J,K)+RO(IMM1,JM1,K)
     &        + RO(IMM1,JM1,KP1)+RO(IMM1,J,KP1)
      WAVG    = SQRT(WTAVGIM*WTAVGIM+VXAVGIM*VXAVGIM+VRAVGIM*VRAVGIM)
      AREAIM  = SQRT(ABX(IM,J,K)*ABX(IM,J,K) + ABR(IM,J,K)*ABR(IM,J,K)
     &              +ABT(J,K)*ABT(J,K)) 
      DPERP   = VOL(IMM1,J,K)/AREAIM
C******************************************************************************
C******************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD =  YPLUS_IM(J,K)
           DENSITY   = 0.25*ROAVGIM
           WSLIP     = 0.25*WAVG
           CALL WALLFUN(IMM1,J,K,IM,DPERP,DPDS_CELL(IMM1,J,K),
     &             DENSITY,TWALLIM,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_IM(J,K) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 375
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
      RE      = DPERP*ROAVGIM*WAVG/VISLAM16
      RELOG   = 1./ALOG(RE)
C
C*******************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_U(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLIM     = CF*ROAVGIM*WAVG*WAVG/128
C
C          VSTAR   = SQRT(TWALLIM/ROAVGIM/4)  
C      Changed by WS/LX  17 May, 2017  
           VSTAR   = SQRT(TWALLIM/(0.25*ROAVGIM))
C
           PLUSK   = ROUGH*VSTAR*ROAVGIM/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C*********************************************************************************
C    SET THE SHEAR STRESS ON THE  I = IM BLADE SURFACE
C
      TWALLIM  = CF*ROAVGIM*WAVG*WAVG/128
C
C*********************************************************************************
C
  375 CONTINUE
C
C*********************************************************************************
C     CALCULATE YPLUS FOR USE IN SETTING  FYPLUS.
C
      VSTAR   = SQRT(TWALLIM/(0.25*ROAVGIM))
      ISTRT = IMID + 1
      DO I = ISTRT,IM
C
C      YPLS = RO(I,J,K)*VSTAR*SQRT(DWALLSQ(I,J,K))/VISLAM
C      Changed by WS/LX  17 May, 2017  
      YPLS = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVGIM)/VISLAM
C
      Y_PLUS(I,J,K) = AMIN1(1000.,YPLS)
      END DO
C
      FMULT    = -TWALLIM*AREAIM/WAVG
C 
      XFORCEIM =  FMULT*VXAVGIM
      RFORCEIM =  FMULT*VRAVGIM
      TFORCEIM =  FMULT*WTAVGIM
      WVISCIM  =  TFORCEIM*WRAD(J)*RAVG_CELL(J,K) 
C
      FORCEX(IMM1,J,K)  = FORCEX(IMM1,J,K)   + XFORCEIM
      FORCER(IMM1,J,K)  = FORCER(IMM1,J,K)   + RFORCEIM
      FORCET(IMM1,J,K)  = FORCET(IMM1,J,K)   + TFORCEIM
      ESOURCE(IMM1,J,K) = ESOURCE(IMM1,J,K)  + WVISCIM
C
  220 CONTINUE
C
C*****************************************************************************
C
C     NOW SET THE FORCES ON THE PERIODIC CELLS, I=1 and I= IM-1,UPSTREAM OF THE 
C     LE AND DOWNSTREAM OF THE TE AND ALSO FOR CELLS IN THE TIP GAP.
C     DO NOT TREAT CELLS ON THE SOLID BLADE SURFACES AS PERIODIC.
C*********************************************************************************
C
       IF( (IND(J).EQ.1).AND.(K.LT.KTIPS(NRW)).OR.(K.GT.KTIPE(NRW)) )
     &      GO TO 225
C
      I   = 1
      IM1 = IMM1
      XFOR = 0.5*((TXR(I,J,K)+TXR(IM1,J,K))*ABR(I,J,K)
     &    +       (TXT(I,J,K)+TXT(IM1,J,K))*ABT(J,K)
     &    +       (TXX(I,J,K)+TXX(IM1,J,K))*ABX(I,J,K))
      FORCEX(IM1,J,K) = FORCEX(IM1,J,K) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR
      AVG = 0.5*(FORCEX(I,J,K) + FORCEX(IM1,J,K))
      FORCEX(I,J,K)   = AVG
      FORCEX(IM1,J,K) = AVG

      RFOR = 0.5*((TRR(I,J,K)+TRR(IM1,J,K))*ABR(I,J,K)
     &    +       (TRT(I,J,K)+TRT(IM1,J,K))*ABT(J,K)
     &    +       (TRX(I,J,K)+TRX(IM1,J,K))*ABX(I,J,K))
      FORCER(IM1,J,K)  = FORCER(IM1,J,K) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR
      AVG = 0.5*(FORCER(I,J,K) + FORCER(IM1,J,K))
      FORCER(I,J,K)   = AVG
      FORCER(IM1,J,K) = AVG

      TFOR = 0.5*((TTR(I,J,K)+TTR(IM1,J,K))*ABR(I,J,K)
     &    +       (TTT(I,J,K)+TTT(IM1,J,K))*ABT(J,K)
     &    +       (TTX(I,J,K)+TTX(IM1,J,K))*ABX(I,J,K))
      FORCET(IM1,J,K)  = FORCET(IM1,J,K) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
      AVG = 0.5*(FORCET(I,J,K) + FORCET(IM1,J,K))
      FORCET(I,J,K)   = AVG
      FORCET(IM1,J,K) = AVG
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(IM1,J,K))*ABX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(IM1,J,K))*ABR(I,J,K)
     &      +       (QTT(I,J,K) + QTT(IM1,J,K))*ABT(J,K))
C
      VXAVG = VX(I,J,K)+VX(I,JM1,K)+VX(I,JM1,KP1)+VX(I,J,KP1)
      VRAVG = VR(I,J,K)+VR(I,JM1,K)+VR(I,JM1,KP1)+VR(I,J,KP1)
      VTAVG = VT(I,J,K)+VT(I,JM1,K)+VT(I,JM1,KP1)+VT(I,J,KP1)
C
      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)
C
      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(IM1,J,K)= ESOURCE(IM1,J,K) - QFLOW - WVISC 
      AVG = 0.5*(ESOURCE(I,J,K) + ESOURCE(IM1,J,K))
      ESOURCE(I,J,K)   = AVG
      ESOURCE(IM1,J,K) = AVG
C
  225 CONTINUE
C
C*********************************************************************************
C*********************************************************************************
C     NOW SUM THE STRESSES ON THE OTHER I = CONSTANT FACES 
C
      DO 250 I=2,IMM1
      IM1  = I-1
      XFOR = 0.5*((TXR(I,J,K)+TXR(IM1,J,K))*ABR(I,J,K)
     &    +       (TXT(I,J,K)+TXT(IM1,J,K))*ABT(J,K)
     &    +       (TXX(I,J,K)+TXX(IM1,J,K))*ABX(I,J,K))
      FORCEX(IM1,J,K) = FORCEX(IM1,J,K) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR

      RFOR = 0.5*((TRR(I,J,K)+TRR(IM1,J,K))*ABR(I,J,K)
     &    +       (TRT(I,J,K)+TRT(IM1,J,K))*ABT(J,K)
     &    +       (TRX(I,J,K)+TRX(IM1,J,K))*ABX(I,J,K))
      FORCER(IM1,J,K)  = FORCER(IM1,J,K) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR

      TFOR = 0.5*((TTR(I,J,K)+TTR(IM1,J,K))*ABR(I,J,K)
     &    +       (TTT(I,J,K)+TTT(IM1,J,K))*ABT(J,K)
     &    +       (TTX(I,J,K)+TTX(IM1,J,K))*ABX(I,J,K))
      FORCET(IM1,J,K)  = FORCET(IM1,J,K) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(IM1,J,K))*ABX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(IM1,J,K))*ABR(I,J,K)
     &      +       (QTT(I,J,K) + QTT(IM1,J,K))*ABT(J,K))
C
      VXAVG = VX(I,J,K)+VX(I,JM1,K)+VX(I,JM1,KP1)+VX(I,J,KP1)
      VRAVG = VR(I,J,K)+VR(I,JM1,K)+VR(I,JM1,KP1)+VR(I,J,KP1)
      VTAVG = VT(I,J,K)+VT(I,JM1,K)+VT(I,JM1,KP1)+VT(I,J,KP1)
C
      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)
C
      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(IM1,J,K)= ESOURCE(IM1,J,K) - QFLOW - WVISC 

  250 CONTINUE

  200 CONTINUE
C
C**************************************************************************************
C***************************************************************************************
C**********CALCULATE THE WALL SHEAR STRESSES ON THE HUB AND CASING SURFACES**********
C
C   Q3D
      IF(KM.EQ.2) GO TO 401
C   END Q3D
      DO 400 J=2,JM
      JM1  = J-1
      DO 400 I=1,IMM1
      IP1 = I+1
C
C     FIRST CALCULATE THE WALL SHEAR STRESS.FOR THE K = 1,  HUB, WALL.   
C
      VXAVG1 = VX(I,J,2)+VX(IP1,J,2)+VX(IP1,JM1,2)+VX(I,JM1,2)
      VRAVG1 = VR(I,J,2)+VR(IP1,J,2)+VR(IP1,JM1,2)+VR(I,JM1,2)
      VTAVG1 = VT(I,J,2)+VT(IP1,J,2)+VT(IP1,JM1,2)+VT(I,JM1,2)
      ROAVG1 = RO(I,J,2)+RO(IP1,J,2)+RO(IP1,JM1,2)+RO(I,JM1,2)
      RAVG1  = R(J,2)   + R(J,2)    + R(JM1,2)    + R(JM1,2)
      WTAVG1 = VTAVG1   - WHUB(J)*RAVG1
      WAVG   = SQRT(VXAVG1*VXAVG1 + VRAVG1*VRAVG1 + WTAVG1*WTAVG1)
      AREA1  = SQRT(ASX(I,J,1)*ASX(I,J,1) + ASR(I,J,1)*ASR(I,J,1))
      DPERP  = VOL(I,J,1)/AREA1
C
C**********************************************************************************
C**********************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD =  YPLUS_K1(I,J)
           DENSITY   = 0.25*ROAVG1
           WSLIP     = 0.25*WAVG
         CALL WALLFUN(I,J,1,1,DPERP,DPDS_CELL(I,J,1),DENSITY,
     &                TWALLK1,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_K1(I,J) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 345

      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
      RE     = DPERP*ROAVG1*WAVG/VISLAM16
      RELOG  = 1./ALOG(RE)
C
C**************************************************************************************
C
      ROUGH  = ROUGH_H(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLK1     = CF*ROAVG1*WAVG*WAVG/128
C
C		   VSTAR   = SQRT(TWALLK1/ROAVG1/4)  
C      Changed by WS/LX  17 May, 2017  
           VSTAR   = SQRT(TWALLK1/(0.25*ROAVG1))
C
           PLUSK   = ROUGH*VSTAR*ROAVG1/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
C
      END IF
C**************************************************************************************
C    SET THE SHEAR STRESS ON THE HUB SURFACE, K=1.
C
      TWALLK1    = CF*ROAVG1*WAVG*WAVG/128
C
C**************************************************************************************
C
  345 CONTINUE
C
C**************************************************************************************
C     CALCULATE YPLUS FOR USE IN SETTING  FYPLUS
      IF(IBOUND.EQ.0.OR.IBOUND.EQ.2) THEN
      DO K = 1,KMID
      VSTAR     = SQRT(TWALLK1/(0.25*ROAVG1))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVG1)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(YPLS,Y_PLUS(I,J,K))
      END DO
      END IF
C
      FMULT   = -TWALLK1*AREA1/WAVG
      IF(IBOUND.EQ.1.OR.IBOUND.GT.2) FMULT = 0.0
C
      XFORCE1 = FMULT*VXAVG1
      RFORCE1 = FMULT*VRAVG1
      TFORCE1 = FMULT*WTAVG1
      WVISC1  = TFORCE1*0.25*WHUB(J)*RAVG1
C
      FORCEX(I,J,1)    = FORCEX(I,J,1)    + XFORCE1
      FORCER(I,J,1)    = FORCER(I,J,1)    + RFORCE1
      FORCET(I,J,1)    = FORCET(I,J,1)    + TFORCE1
      ESOURCE(I,J,1)   = ESOURCE(I,J,1)   + WVISC1
C
C**************************************************************************************
C**************NOW SET THE SHEAR STRESS ON  THE CASING, K = KM ************************
C
      VXAVGKM = VX(I,J,KMM1)+VX(IP1,J,KMM1)
     &        + VX(IP1,JM1,KMM1)+VX(I,JM1,KMM1)
      VRAVGKM = VR(I,J,KMM1)+VR(IP1,J,KMM1)
     &        + VR(IP1,JM1,KMM1)+VR(I,JM1,KMM1)
      VTAVGKM = VT(I,J,KMM1)+VT(IP1,J,KMM1)
     &        + VT(IP1,JM1,KMM1)+VT(I,JM1,KMM1)
      ROAVGKM = RO(I,J,KMM1)+RO(IP1,J,KMM1)
     &        + RO(IP1,JM1,KMM1)+RO(I,JM1,KMM1)
      RAVGKM  = R(J,KMM1) + R(J,KMM1) + R(JM1,KMM1) + R(JM1,KMM1)
      WTAVGKM = VTAVGKM   - WTIP(J)*RAVGKM
      WAVG    = SQRT(VXAVGKM*VXAVGKM+VRAVGKM*VRAVGKM+WTAVGKM*WTAVGKM)
      AREAKM  = SQRT(ASX(I,J,KM)*ASX(I,J,KM) + ASR(I,J,KM)*ASR(I,J,KM))
      DPERP   = VOL(I,J,KMM1)/AREAKM
C
C************************************************************************
C****************************************************************************
C    USE THE Shih et al  WALL FUNCTIONS IF YPLUSWALL IS NEGATIVE.
      IF(YPLUSWALL.LT.-0.001) THEN

           YPLUS_OLD  =  YPLUS_KM(I,J)
           DENSITY    = 0.25*ROAVGKM
           WSLIP      = 0.25*WAVG
        CALL WALLFUN(I,J,KMM1,KM,DPERP,DPDS_CELL(I,J,KMM1),DENSITY,
     &              TWALLKM,YPLUS_OLD,WSLIP,YPLUS_NEW)
           YPLUS_KM(I,J) = AMIN1(1000.0,YPLUS_NEW)
C
      GO TO 355
C
      END IF
C    END OF Shih et al WALLFUNCTIONS
C**********************************************************************************
C**********************************************************************************
C
      RE      = DPERP*ROAVGKM*WAVG/VISLAM16
      RELOG   = 1./ALOG(RE)
C
C**************************************************************************************
C    ALLOW FOR ROUGHNESS
C
      ROUGH  = ROUGH_T(NRW)
C
      IF(ROUGH.GT.1.0E-7) THEN
           IF(ROUGH.GT.DPERP) ROUGH = DPERP
           REK   = RE*ROUGH/DPERP
           REK   = REK - 80.0
           IF(REK.LT.0.0)  REK = 0.0
           REKSQ = REK*REK
           A1 = -.00178493 + .0000814923*REK + .000000150445*REKSQ
           A2 = .029072 - .001584*REK - .00000225194*REKSQ
           A3=.270313+.0091409*REK+.00000451537*REKSQ +
     &     .00000000464767*REKSQ*REK
           CF    = A1 + A2*RELOG + A3*RELOG*RELOG
      ELSE
C
C    END ROUGHNESS, NEXT EQN FOR SMOOTH SURFACES.
C
           CF = -0.00178493 + 0.029072*RELOG + 0.270313*RELOG*RELOG
C
      END IF
C
C   TAKE CF AS THE MAX OF THE LAMINAR AND TURBULENT VALUES.
C
      CFLAM = 2.0/RE
      CF    = AMAX1(CF,CFLAM)
      IF(RE.LT.125.) CF = CFLAM
C
C  CHECK FOR FULLY ROUGH SURFACES
C
      IF(ROUGH.GT.1.0E-7) THEN
           TWALLKM    = CF*ROAVGKM*WAVG*WAVG/128
C
C		   VSTAR   = SQRT(TWALLKM/ROAVGKM/4)  
C      Changed by WS/LX  17 May, 2017  
           VSTAR   = SQRT(TWALLKM/(0.25*ROAVGKM))
C
           PLUSK   = ROUGH*VSTAR*ROAVGKM/VISLAM/4
           RELOG10 = LOG10(DPERP/ROUGH)
           FUNCN   = 5.75*RELOG10 + 8.5
           CF_FULL_ROUGH = 2.0/(FUNCN*FUNCN)
           IF(PLUSK.GT.45.AND.CF.GT.CF_FULL_ROUGH) CF = CF_FULL_ROUGH
      END IF
C
C**************************************************************************************
C   SET THE SHEAR STRESS ON THE CASING , K = KM, SURFACE
C
      TWALLKM    = CF*ROAVGKM*WAVG*WAVG/128
C
C**************************************************************************************
C
  355 CONTINUE
C
C**************************************************************************************
C
C     CALCULATE YPLUS FOR USE IN SETTING FYPLUS.
      IF(IBOUND.EQ.0.OR.IBOUND.EQ.1) THEN
      DO K = KMID+1,KMM1
      VSTAR     = SQRT(TWALLKM/(0.25*ROAVGKM))
      YPLS      = VSTAR*SQRT(DWALLSQ(I,J,K))*(0.25*ROAVGKM)/VISLAM
      Y_PLUS(I,J,K) = AMIN1(YPLS,Y_PLUS(I,J,K))
      END DO
      END IF
C
      FMULT   = -TWALLKM*AREAKM/WAVG
      IF(IBOUND.GE.2) FMULT = 0.0 
C
      XFORCEKM = FMULT*VXAVGKM
      RFORCEKM = FMULT*VRAVGKM
      TFORCEKM = FMULT*WTAVGKM
      WVISCKM  =  TFORCEKM*0.25*WTIP(J)*RAVGKM
C
      FORCEX(I,J,KMM1) = FORCEX(I,J,KMM1) + XFORCEKM
      FORCER(I,J,KMM1) = FORCER(I,J,KMM1) + RFORCEKM
      FORCET(I,J,KMM1) = FORCET(I,J,KMM1) + TFORCEKM
      ESOURCE(I,J,KMM1)= ESOURCE(I,J,KMM1)+ WVISCKM
C
C     NOW SUM THE STRESSES ON THE OTHER K = CONSTANT FACES
C
      DO 450 K=2,KMM1
      KM1  = K-1
      XFOR = 0.5*((TXR(I,J,K)+TXR(I,J,KM1))*ASR(I,J,K)
     &    +       (TXX(I,J,K)+TXX(I,J,KM1))*ASX(I,J,K))
      FORCEX(I,J,KM1) = FORCEX(I,J,KM1) + XFOR
      FORCEX(I,J,K)   = FORCEX(I,J,K)   - XFOR

      RFOR = 0.5*((TRR(I,J,K)+TRR(I,J,KM1))*ASR(I,J,K)
     &    +       (TRX(I,J,K)+TRX(I,J,KM1))*ASX(I,J,K))
      FORCER(I,J,KM1)  = FORCER(I,J,KM1) + RFOR
      FORCER(I,J,K)    = FORCER(I,J,K)   - RFOR

      TFOR = 0.5*((TTR(I,J,K)+TTR(I,J,KM1))*ASR(I,J,K)
     &    +       (TTX(I,J,K)+TTX(I,J,KM1))*ASX(I,J,K))
      FORCET(I,J,KM1)  = FORCET(I,J,KM1) + TFOR
      FORCET(I,J,K)    = FORCET(I,J,K)   - TFOR
C
C     SET THE HEAT FLOW AND VISCOUS WORK
C
      QFLOW =  0.5*((QXX(I,J,K) + QXX(I,J,KM1))*ASX(I,J,K)
     &      +       (QRR(I,J,K) + QRR(I,J,KM1))*ASR(I,J,K))
C
      VXAVG = VX(I,J,K)+VX(IP1,J,K)+VX(IP1,JM1,K)+VX(I,JM1,K)
      VRAVG = VR(I,J,K)+VR(IP1,J,K)+VR(IP1,JM1,K)+VR(I,JM1,K)
      VTAVG = VT(I,J,K)+VT(IP1,J,K)+VT(IP1,JM1,K)+VT(I,JM1,K)

      WVISC = -0.25*(XFOR*VXAVG + RFOR*VRAVG + TFOR*VTAVG)

      ESOURCE(I,J,K)  = ESOURCE(I,J,K)   + QFLOW + WVISC
      ESOURCE(I,J,KM1)= ESOURCE(I,J,KM1) - QFLOW - WVISC 
C
  450 CONTINUE
C
  400 CONTINUE
C  RE ENTER HERE IF Q3D
  401 CONTINUE
C   END Q3D 
C**************************************************************************************
C**************************************************************************************
C     ALL THE VISCOUS FORCES ARE NOW SET. USE THEM TO UPDATE THE GLOBAL BODY FORCE TERMS.
C     NOTE THE NEGATIVE SIGN BECAUSE THE SOURCE TERM IS SUBTRACTED IN SUBROUTINE TSTEP.
C 
      DO 500 K=1,KMM1
      DO 500 J=2,JM
      DO 500 I=1,IMM1
      XFORCE(I,J,K)   = RF_VIS1*XFORCE(I,J,K)   - RF_VIS*FORCEX(I,J,K)
      TFORCE(I,J,K)   = RF_VIS1*TFORCE(I,J,K)   - RF_VIS*FORCET(I,J,K)
     &                                           *RAVG_CELL(J,K)
      RFORCE(I,J,K)   = RF_VIS1*RFORCE(I,J,K)   - RF_VIS*FORCER(I,J,K)
      QSOURCE(I,J,K)  = RF_VIS1*QSOURCE(I,J,K)  - RF_VIS*ESOURCE(I,J,K)
  500 CONTINUE
C
C**************************************************************************************
C**************************************************************************************
C
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C    
      SUBROUTINE SMOOTH_RESID(VAR,SMFAC,NSMTH)
C
C******************************************************************************************
C     THIS SUBROUTINE APPLIES A SIMPLE LINEAR SMOOTHING TO THE
C     VARIABLE VAR(I,J,K). A ZERO GRADIENT CONDITION IS APPLIED AT THE BOUNDARIES.
C     WITH NO SMOOTHING ACROSS THE MIXING PLANE.
C******************************************************************************************
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION TEMPVAR(JD),VAR(ID,JD,KD)
C
      SFAC1 = 1.0 -SMFAC
      SFACH = 0.5*SMFAC
C 
      DO 200, NSM=1,NSMTH
C
C     SMOOTH IN THE J DIRECTION
C
      DO 100 NR = 1,NROWS
      J1   = JSTART(NR)
      J1P1 = J1+1
      J1P2 = J1+2
      J2   = JMIX(NR)
      J2M1 = J2-1
C
      DO 40 K=1,KMM1
      DO 40 I=1,IMM1
      DO 50 J=J1P2,J2M1
      TEMPVAR(J) = SFAC1*VAR(I,J,K)
     &           + SFACH*(VAR(I,J+1,K)+VAR(I,J-1,K))
   50 CONTINUE
      TEMPVAR(J1P1) = VAR(I,J1P1,K)+ SFACH*(VAR(I,J1P2,K)-VAR(I,J1P1,K))
      TEMPVAR(J2)   = VAR(I,J2,K)  + SFACH*(VAR(I,J2M1,K)-VAR(I,J2,K))
      DO 60 J = J1P1,J2
      VAR(I,J,K) = TEMPVAR(J)
   60 CONTINUE
   40 CONTINUE
C
C     SMOOTH IN THE I DIRECTION
C     TFLOW ADDITION
      IF(IM.EQ.2) GO TO 11

      DO 10 K=1,KMM1
      DO 10 J=J1P1,J2
      DO 20 I=2,IMM2
      TEMPVAR(I) = SFAC1*VAR(I,J,K)
     &           + SFACH*(VAR(I+1,J,K)+VAR(I-1,J,K))
   20 CONTINUE
      TEMPVAR(1)   = VAR(1,J,K)    + SFACH*(VAR(2,J,K)   -VAR(1,J,K))
      TEMPVAR(IMM1)= VAR(IMM1,J,K) + SFACH*(VAR(IMM2,J,K)-VAR(IMM1,J,K))
      DO 30 I=1,IMM1
      VAR(I,J,K) = TEMPVAR(I)
   30 CONTINUE
   10 CONTINUE

   11 CONTINUE
C
C  Q3D
      IF(KM.EQ.2) GO TO 100
C  END Q3D
C
C     SMOOTH IN THE K DIRECTION
C
      DO 70 I=1,IMM1
      DO 70 J=J1P1,J2
      DO 80 K=2,KMM2
      TEMPVAR(K) = SFAC1*VAR(I,J,K) 
     &           + SFACH*(VAR(I,J,K+1)+VAR(I,J,K-1))
   80 CONTINUE
      TEMPVAR(1)   = VAR(I,J,1)    + SFACH*(VAR(I,J,2)   -VAR(I,J,1))
      TEMPVAR(KMM1)= VAR(I,J,KMM1) + SFACH*(VAR(I,J,KMM2)-VAR(I,J,KMM1))
      DO 90 K=1,KMM1
      VAR(I,J,K) = TEMPVAR(K)
   90 CONTINUE
   70 CONTINUE
C
  100 CONTINUE
C
  200 CONTINUE
C
C******************************************************************************************
C******************************************************************************************
      RETURN
      END
C
C*****************************************************************************
C*****************************************************************************
C
      SUBROUTINE NEW_MIXPLAN
C
C******************************************************************************************
C     THIS SUBROUTINE RELAXES THE PROPERTIES DOWNSTREAM OF THE MIXING PLANE TOWARDS PITCHWISE
C     UNIFORM ISENTROPIC VALUES.AND ALSO ALLOWS THE FLOW ANGLES TO BE EXTRAPOLATED FROM DOWNSTREAM.
C     THE FLOW CAN CROSS THE MIXING PLANE IN EITHER DIRECTION.
C******************************************************************************************
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION  ALPHAA(ID),GAMA(ID),WREF_JP2(ID),ROREF_JP2(ID),
     &           HSTAT(ID),PISENT(ID),ROISENT(ID),
     &           WMACH_JP2(ID),TISENT(ID),WIS_JP2(ID),RO_SUB_JP2(ID)
C
      RFMIX1   = 1.0 - RFMIX
C
      DO 500 NR = 1,NRWSM1  
C
      DO 400  K = 1,KM
C
      KSUB = K
      IF(K.EQ.KM) KSUB = KM-1
C
C   CHECK THE FLOW DIRECTION AT THE MIXING PLANE.
C
           J     = JMIX(NR)
      IF(FLOWX(IMID,J,KSUB).LT.0.0) THEN 
           JP1   = J+1
           JP2   = J+2  
           JP3   = J+3 
           JP4   = J+4
           F_ANGLE  = FANGLE
      ELSE
           J   = JMIX(NR) + 1
           JP1 = J-1
           JP2 = J-2
           JP3 = J-3
           JP4 = J-4
           F_ANGLE  = 0.0
      END IF
C
           F_ANGLE1 = 1.0 - F_ANGLE
C
C*************************************************************************************************
C   CALCULATE THE MIXED OUT FLOW CONDITION AT THE MIXING PLANE. THE FLOW THERE IS PITCHWISE UNIFORM
C   SO WE ONLY NEED TO USE THE MID-PITCH VALUES.
C
      RJP1     = R(JP1,K)
      UJP1     = UBLADE(JP1,K)
      ROMIX    = RO(IMID,JP1,K)
      VXMIX    = ROVX(IMID,JP1,K)/ROMIX
      VRMIX    = ROVR(IMID,JP1,K)/ROMIX
      VTMIX    = RORVT(IMID,JP1,K)/ROMIX/RJP1
      ROEMIX   = ROE(IMID,JP1,K)
      HOMIX    = HO(IMID,JP1,K)   
      VMSQ     = VXMIX*VXMIX + VRMIX*VRMIX
      VMMIX    = SQRT(VMSQ)
      VABSQ    = VMSQ + VTMIX*VTMIX
      WTMIX    = VTMIX - UJP1
      WMIXSQ   = VMSQ  + WTMIX*WTMIX
      WMIX     = SQRT(WMIXSQ)
      ALPHAMIX = ASIN(WTMIX/WMIX)
      GAMAMIX  = ATAN2(VRMIX,VXMIX)

C
      IF(IFGAS.EQ.0) THEN
           TMIX      = (ROEMIX/ROMIX - 0.5*VABSQ)/CV
C
           IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) 
     &      TMIX = (ROEMIX/ROSUB(IMID,JP1,K) - 0.5*VABSQ)/CV
C
           ROTHALPY  = CP*TMIX + 0.5*(WMIXSQ - UJP1*UJP1)
           GAGAS     = GA
      ELSE
            EMIX     = ROEMIX/ROMIX - 0.5*VABSQ
            TMIX     = TFROME(EMIX,TREF,EREF,ET1,ET2,ET3,ET4)
            HSTATIC  = HFROMT(TMIX,TREF,HREF,CP1,CP2,CP3)
            ROTHALPY = HSTATIC + 0.5*(WMIXSQ - UJP1*UJP1)
            DELT     =  TMIX - TREF
            CPGAS    = CP1 + CP2*DELT + CP3*DELT*DELT
            GAGAS    = CPGAS/(CPGAS-RGAS)
      END IF
C
      PMIX    = ROMIX*RGAS*TMIX
C
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) 
     &   PMIX = PO_REF + DP_DRO*(ROSUB(IMID,JP1,K) - RO_REF)
C
      IF(ITIMST.EQ.6) PORELMIX = PMIX + 0.5*DENSTY*WMIXSQ
C
C
C*********************************************************************************************
C
      NMACH = 0
      DO 350 I = 1,IM
C
C     CALCULATE THE VELOCITIES AT JMIX + 2
C
      RJP2       = R(JP2,K)
      UNOWJP2    = UBLADE(JP2,K)
      ROJP2      = RO(I,JP2,K)
      XVELJP2    = ROVX(I,JP2,K)/ROJP2
      RVELJP2    = ROVR(I,JP2,K)/ROJP2
      TVELJP2    = RORVT(I,JP2,K)/RJP2/ROJP2
      WTJP2      = TVELJP2 - UNOWJP2
      VMSQJP2    = XVELJP2*XVELJP2 + RVELJP2*RVELJP2
      VELSQJP2   = VMSQJP2 + TVELJP2*TVELJP2
      WACT_SQ_JP2= VMSQJP2 + WTJP2*WTJP2
      WACT_JP2   = SQRT(WACT_SQ_JP2)
      RO_ACT_JP2 = ROJP2
      T_ACT_JP2  = (ROE(I,JP2,K)/ROJP2 - 0.5*VELSQJP2)/CV
C
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) 
     &       T_ACT_JP2 = (ROE(I,JP2,K)/ROSUB(I,JP2,K) - 0.5*VELSQJP2)/CV
C
C
C**********************************************************************************************
C    OBTAIN VELOCITIES AND ANGLES AT JMIX + 3. ONLY THE ANGLES ARE USED.
C
      RJP3       = R(JP3,K)
      UNOWJP3    = UBLADE(JP3,K)
      ROJP3      = RO(I,JP3,K)
      ROVXJP3    = ROVX(I,JP3,K)
      ROVRJP3    = ROVR(I,JP3,K)
      ROVTJP3    = RORVT(I,JP3,K)/RJP3
      ROWTJP3    = ROVTJP3 - UNOWJP3*ROJP3
      ROVMSQJP3  = ROVXJP3*ROVXJP3 + ROVRJP3*ROVRJP3
      ROWRELSQJP3= ROVMSQJP3 + ROWTJP3*ROWTJP3
      ROWRELJP3  = SQRT(ROWRELSQJP3)
      ALPHAJP3   = ASIN(ROWTJP3/ROWRELJP3)
      GAMAJP3    = ATAN2(ROVRJP3,ROVXJP3)
C
C*********************************************************************************************
C    OBTAIN VELOCITIES AND ANGLES AT JMIX + 4. ONLY THE ANGLES ARE USED.
C
      RJP4       = R(JP4,K)
      UNOWJP4    = UBLADE(JP4,K)
      ROJP4      = RO(I,JP4,K)
      ROVXJP4    = ROVX(I,JP4,K)
      ROVRJP4    = ROVR(I,JP4,K)
      ROVTJP4    = RORVT(I,JP4,K)/RJP4
      ROWTJP4    = ROVTJP4 - UNOWJP4*ROJP4
      ROVMSQJP4    = ROVXJP4*ROVXJP4 + ROVRJP4*ROVRJP4
      ROWRELSQJP4  = ROVMSQJP4 + ROWTJP4*ROWTJP4
      ROWRELJP4    = SQRT(ROWRELSQJP4)
      ALPHAJP4     = ASIN(ROWTJP4/ROWRELJP4)
      GAMAJP4      = ATAN2(ROVRJP4,ROVXJP4)
C
C*********************************************************************************************
C
C   NEXT THE ISENTROPIC VALUES AT JMIX+2 ASSUMING AN ISENTROPIC EXPANSION FROM THE MIXING PLANE CONDITIONS
C   TO  PSTAT  AT JMIX+2 .
C
      IF(IFGAS.EQ.0) THEN
           PSTAT      = ROJP2*RGAS*T_ACT_JP2

           IF(ITIMST.EQ.5.OR.ITIMST.EQ.6)
     &     PSTAT   = PO_REF  + DP_DRO*(ROSUB(I,JP2,K) - RO_REF)

           TISENT_JP2 = TMIX*(PSTAT/PMIX)**FGA
           WIS_SQ_JP2 = 2*(ROTHALPY- CP*TISENT_JP2+ 0.5*UNOWJP2*UNOWJP2)
C
           IF(ITIMST.EQ.6) WIS_SQ_JP2 = 2.0*(PORELMIX - PSTAT)/DENSTY
C
      ELSE
           ESTAT  = ROE(I,JP2,K)/ROJP2 - 0.5*VELSQJP2
           TSTAT  = TFROME(ESTAT,TREF,EREF,ET1,ET2,ET3,ET4)
           T_ACT_JP2 = TSTAT
           PSTAT  = ROJP2*RGAS*TSTAT
           PRAT   = PSTAT/PMIX
           TRAT   = TRAT_FROM_PRAT
     &     (TMIX,PRAT,FGAGAS,R_ALPHA,BALPHA1,BALPHA2)
           TISENT_JP2 = TMIX*TRAT
           HSTAT(I)  = HFROMT(TISENT_JP2,TREF,HREF,CP1,CP2,CP3)
           WIS_SQ_JP2  = 2*(ROTHALPY - HSTAT(I) + 0.5*UNOWJP2*UNOWJP2)
      END IF
C
      ROISENT_JP2  = PSTAT/RGAS/TISENT_JP2
C
      IF(ITIMST.EQ.6) ROISENT_JP2 = DENSTY
C
      IF(ITIMST.EQ.5.OR.ITIMST.EQ.6) 
     &       RO_SUB_JP2(I)= RO_REF + (PSTAT - PO_REF)/DP_DRO

      ROISENT(I)  = ROISENT_JP2
      PISENT(I)   = PSTAT
      TISENT(I)   = TISENT_JP2
C
      IF(WIS_SQ_JP2.LT.1.0) WIS_SQ_JP2 = 1.0
      WIS_JP2(I)     = SQRT(WIS_SQ_JP2)
C
C**********************************************************************************************
C     SET THE CURRENT VALUES AT JMIX+2 FOR USE IN FINDING THE MACH NUMBER AND ANGLE VARIATION.
C
      WREF_JP2(I)    =  WACT_JP2 
      ROREF_JP2(I)   =  RO_ACT_JP2
      WMACH_JP2(I)   =  WACT_JP2/SQRT(GAGAS*RGAS*T_ACT_JP2)
      IF(WMACH_JP2(I).GT.1.0) NMACH = NMACH + 1
C
C***********************************************************************************************
C    SET THE FLOW ANGLES AT JMIX+2 WEIGHTING THE FLOW DIRECTION ACCORDING TO F_ANGLE
C
C      ALPHAA(I) = F_ANGLE1*ALPHAMIX  + F_ANGLE*(2*ALPHAJP3 - ALPHAJP4)
C      GAMA(I)   = F_ANGLE1*GAMAMIX   + F_ANGLE*(2*GAMAJP3  - GAMAJP4)
      ALPHAA(I) = F_ANGLE1*ALPHAMIX  + F_ANGLE*ALPHAJP3
      GAMA(I)   = F_ANGLE1*GAMAMIX   + F_ANGLE*GAMAJP3
C      ALPHAA(I) = F_ANGLE1*ALPHAMIX  + F_ANGLE*(ALPHAJP3 +ALPHAJP4)/2
C      GAMA(I)   = GAMAMIX
C
  350 CONTINUE
C
C*********************************************************************************
C  
      IF(ITIMST.GE.5) GO TO 250
      IF(NMACH.EQ.0)  GO TO 250
C
C********************************************************************************
C     APPLY CHARACTERISTICS RELATIONSHIP BETWEEN THE YAW ANGLE AND MACH No. IF
C     THE RELATIVE FLOW IS SUPERSONIC
C
C**********************************************************************************************
C    SET FROTN TO ALLOW FOR LEFT OR RIGHT RUNNING MACH WAVES
C      
      FROTN  =  1.0
      IF(WTMIX.LT.0.0) FROTN = -1.0
C
C   IF THE RELATIVE FLOW IS SUPERSONIC,CALCULATE THE PITCHWISE VARIATION IN ANGLE NEEDED
C   TO SATISFY THE PRANDTL-MEYER RELATIONSHIP.
C
      DO 211 I=IMID+1,IM 
      IM1 = I-1  
      IF(WMACH_JP2(I).GT.1.001.AND.WMACH_JP2(IM1).GT.1.001) THEN      
      DELP     = PISENT(I) - PISENT(IM1)   
      DELTAI   = SQRT(WMACH_JP2(I)*WMACH_JP2(I) - 1 )/
     &           (ROREF_JP2(I)*WREF_JP2(I)*WREF_JP2(I))     
      DELTAIM1 = SQRT(WMACH_JP2(IM1)*WMACH_JP2(IM1) -1 )/
     &           (ROREF_JP2(IM1)*WREF_JP2(IM1)*WREF_JP2(IM1))
      ALPHAA(I) = ALPHAA(IM1) + FROTN*DELP*0.5*(DELTAI + DELTAIM1)
      END IF       
  211 CONTINUE  
C
      DO 212 I=IMID-1,1,-1   
      IP1  = I + 1
      IF(WMACH_JP2(I).GT.1.001.AND.WMACH_JP2(IP1).GT.1.001) THEN      
      DELP     = PISENT(I) - PISENT(IP1)   
      DELTAI   = SQRT(WMACH_JP2(I)*WMACH_JP2(I) - 1 )/
     &           (ROREF_JP2(I)*WREF_JP2(I)*WREF_JP2(I))     
      DELTAIP1 = SQRT(WMACH_JP2(IP1)*WMACH_JP2(IP1) -1 )/
     &           (ROREF_JP2(IP1)*WREF_JP2(IP1)*WREF_JP2(IP1))
      ALPHAA(I) = ALPHAA(IP1) + FROTN*DELP*0.5*(DELTAI + DELTAIP1)
      END IF       
  212 CONTINUE        
C
C    ENFORCE PERIODICITY OF ANGLE
C  
      ROTN  = ALPHAA(IM) - ALPHAA(1)
      SUMFP = 0.0
      DO 230 I=1,IM
      ALPHAA(I) = ALPHAA(I) - ROTN*(SUMFP - 0.5)
      SUMFP = SUMFP + FP(I)
  230 CONTINUE
C
C******************************************************************************
  250 CONTINUE
C******************************************************************************
C 
C    STORE THE VALUE OF PITCHWISE FLOW ANGLE ON THE FIRST TIME STEP
C
C      IF(NSTEP.EQ.1) THEN
C      DO 240 I=1,IM
C      ALPHA_S(I,K,NR) = ALPHAA(I)
C  240 CONTINUE
C      END IF
C****************************************************************************************
C****************************************************************************************
C
      DO 300 I = 1,IM
C
C     RELAX CHANGES TO THE VALUE OF ALPHAA BY RFMIX
C
C      ALPHAA(I)   = RFMIX1*ALPHA_S(I,K,NR) + RFMIX*ALPHAA(I)
C      ALPHA_S(I,K,NR) = ALPHAA(I)
C
C     SET THE ISENTROPIC VELOCITIES AT JMIX+2. BASED ON THE CALCULATED RELATIVE VELOCITY MAGNITUDE
C     AND THE RELATIVE FLOW DIRECTION.
C
      RO_JP2    = ROISENT(I)      
      WT_JP2    = WIS_JP2(I)*SIN(ALPHAA(I))
      VM_JP2    = WIS_JP2(I)*COS(ALPHAA(I))
      VX_JP2    = VM_JP2*COS(GAMA(I))
      VR_JP2    = VM_JP2*SIN(GAMA(I))
      VM_JP2SQ  = VX_JP2*VX_JP2 + VR_JP2*VR_JP2
      VT_JP2    = WT_JP2 + UNOWJP2
      V_JP2_SQ  = VM_JP2*VM_JP2 + VT_JP2*VT_JP2
C
      IF(IFGAS.EQ.0) THEN
              ROE_JP2   = ROISENT(I)*(CV*TISENT(I)  + 0.5*V_JP2_SQ)
C
              IF(ITIMST.EQ.5.OR.ITIMST.EQ.6)
     &        ROE_JP2   = RO_SUB_JP2(I)*(CV*TISENT(I)  + 0.5*V_JP2_SQ)
C
      ELSE
              EISENT    = HSTAT(I) - RGAS*TISENT(I) + 0.5*V_JP2_SQ
              ROE_JP2   = EISENT*ROISENT(I)
      END IF
C
C*****************************************************************************************
C*****************************************************************************************  
C     RELAX THE PRIMARY VARIABLES AT J = JMIX + 2 TOWARDS THE "ISENTROPIC"  VALUES 
C
      RO(I,JP2,K)    = RFMIX1*RO(I,JP2,K)   + RFMIX*RO_JP2
      ROVX(I,JP2,K)  = RFMIX1*ROVX(I,JP2,K) + RFMIX*RO_JP2*VX_JP2
      ROVR(I,JP2,K)  = RFMIX1*ROVR(I,JP2,K) + RFMIX*RO_JP2*VR_JP2
      RORVT(I,JP2,K) = RFMIX1*RORVT(I,JP2,K)+ RFMIX*RO_JP2*VT_JP2*RJP2
      ROE(I,JP2,K)   = RFMIX1*ROE(I,JP2,K)  + RFMIX*ROE_JP2
C
  300 CONTINUE      
C
C******************************************************************************************
C
  400 CONTINUE
C
C******************************************************************************************
C
  500 CONTINUE
C
C******************************************************************************************
C
      RETURN
      END
C
C******************************************************************************************
C******************************************************************************************
C******************************************************************************************
C
      SUBROUTINE INSECT(JD,MAXKI,IIN,XI,YI,JIN,XJ,YJ,XINSECT,YINSECT)
C
C   THIS SUBROUTINE FINDS THE INTERSECTION POINT OF TWO LINES DEFINED BY THEIR COORDINATES
C   XI,YI   and  XJ,YJ .
C
      DIMENSION   XI(JD),YI(JD),XJ(MAXKI),YJ(MAXKI),DPERP(500)
C
C
C   FIRST FIND THE TWO CLOSEST POINTS ON THE TWO CURVES.
C
      DMINOV = 1.0E12
      DO I=1,IIN
      DMIN_I= 1.0E12
C
      DO J=1,JIN
      XDIF = XJ(J) - XI(I)
      YDIF = YJ(J) - YI(I)
      DIST = XDIF*XDIF + YDIF*YDIF
      IF(DIST.LT.DMIN_I) THEN
           JMIN   = J
           DMIN_I = DIST
      END IF
      END DO
C
      IF(DMIN_I.LT.DMINOV) THEN
           IMINOV = I
           JMINOV = JMIN
           DMINOV = DMIN_I
      END IF 
      END DO
C
C      WRITE(6,*) ' THE CLOSEST DISTANCE BETWEEN POINTS =', SQRT(DMINOV)
C      WRITE(6,*) ' AT  I = ',IMINOV, ' J= ', JMINOV
C
C    FORM A SIMPLE ESTIMATE OF THE SLOPE OF THE CURVE AT THE CLOSEST POINT.
C
      I = IMINOV
      IM1 = I-1
      IF(I.EQ.1)   IM1 = 1
      IP1 = I+1
      IF(I.EQ.IIN) IP1 = IIN
      DX = XI(IP1) - XI(IM1)
      DY = YI(IP1) - YI(IM1)
      DYDX = DY/DX
      DS   = SQRT(DX*DX + DY*DY)
      VECX =  DX/DS
      VECY =  DY/DS
C
C   DECIDE WHICH IS THE NEXT TO CLOSEST POINT
C   AND FORM A BETTER ESTIMATE OF SLOPE OF THE CURVE AT THE CLOSEST POINT.
C
C
      PROJ = (XJ(JMINOV)-XI(I))*VECX + (YJ(JMINOV)-YI(I))*VECY
C
      IF(PROJ.GT.0.0) THEN
               DX = XI(IP1) - XI(I)
               DY = YI(IP1) - YI(I)
      ELSE
               DX = XI(I) - XI(IM1)
               DY = YI(I) - YI(IM1)
      END IF
C
C    NOW FIND A UNIT NORMAL VECTOR TO THE "I"  LINE AT THE CLOSEST POINT.
C
      DYDX = DY/DX
      DS   = SQRT(DX*DX + DY*DY)
      VECX =  DX/DS
      VECY =  DY/DS
      VNORMX = -VECY
      VNORMY =  VECX
C
C    NOW FIND THE PERPENDICULAR DISTANCE BETWEEN THE CLOSEST POINT ON THE "I" LINE
C    AND ALL POINTS ON THE "J" LINE
C
      DO 20 J=1,JIN
      XDIF = XJ(J) - XI(I)
      YDIF = YJ(J) - YI(I)
      DPERP(J) = XDIF*VNORMX + YDIF*VNORMY
   20 CONTINUE
C
C     INTERPOLATE TO FIND THE POINT WHERE THE PERPENDICULAR DISTANCE = ZERO.
C
      CALL INTP(JIN,DPERP,XJ,0.0,XINSECT)
      CALL INTP(JIN,DPERP,YJ,0.0,YINSECT)
C
C      WRITE(6,*) ' XINSECT, YINSECT = ', XINSECT,YINSECT
C
      RETURN
C
      END
C
C******************************************************************************************
C******************************************************************************************
C******************************************************************************************
C  
      SUBROUTINE SMOOTH_VAR(D)   
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION D(ID,JD,KD),AVG(JD),CURVE(JD),SCURVE(JD)
      NUM3  = 3
C
C******************************************************************************************
C
C     APPLY COMBINED 2nd and 4th ORDER STREAMWISE SMOOTHING TO THE VARIABLE D
C
C******************************************************************************************
C
      IF(SFX.LT.0.0001) GO TO 4001
C
C      STREAMWISE SMOOTHING WITH CONSTANT COEFFICIENT, SFX
C
      DO 3990 NR = 1,NROWS
           J1 = JSTART(NR)
           J2 = JMIX(NR)
      DO 3900 J = J1+3,J2-3
      DO 3900 K=1,KM
      DO 3900 I=1,IM
           AVERG = 0.5*(D(I,J+1,K)+D(I,J-1,K))
           CURV  = AVERG - 0.5*D(I,J,K) - 0.25*(D(I,J-2,K)+D(I,J+2,K))
           STORE(I,J,K) = SFX1*D(I,J,K) + SFX*(AVERG + FAC_4TH*CURV)
 3900 CONTINUE
C******************************************************************************************
C******************************************************************************************
C      EXTRA SMOOTHING AT THE UPSTREAM, DOWNSTREAM AND AND MIXING PLANE BOUNDARIES.
C      THIS MODIFIED BY JDD  MARCH 2015. NEW VERSION GIVES BETTER RESULTS FOR SUPERSONIC
C      FLOWS AT THE MIXING PLANE
C
      DO 3905 K=1,KM
      DO 3905 I=1,IM
      STORE(I,J2,K)   = D(I,J2,K)
      STORE(I,J1,K)   = D(I,J1,K)
 3905 CONTINUE
C
      DO 3910 K=1,KM
      DO 3910 I=1,IM
C
      STORE(I,J1+2,K) =   SFXB1*D(I,J1+2,K) 
     &                +  SFXBH*(D(I,J1+3,K) + D(I,J1+1,K)  )
C
      STORE(I,J1+1,K) = SFXB1*D(I,J1+1,K) + SFXB*(1.25*D(I,J1+2,K)
     &                 + 0.5*D(I,J1+3,K) -  0.75*D(I,J1+4,K) )
C 
      STORE(I,J2-2,K) =  SFXB1*D(I,J2-2,K)
     &                + SFXBH*(D(I,J2-1,K)  + D(I,J2-3,K)  )
C
      STORE(I,J2-1,K) = SFXB1*D(I,J2-1,K) + SFXB*(1.25*D(I,J2-2,K)
     &                  + 0.5*D(I,J2-3,K) - 0.75*D(I,J2-4,K) )
C
 3910 CONTINUE
C
 3990 CONTINUE
C
C     END OF STREAMWISE SMOOTHING
C
C*******************************************************************************
C*******************************************************************************
C
C      NOW APPLY SPECIAL SMOOTHING AROUND THE LEADING EDGE
C
      DO 4006 NR = 1,NROWS
      JS = JLED(NR)-2
      JE = JS + 1
      DO 4004 J = JS,JE
      DO 4004 K=1,KM
      STORE(1,J,K)  = SFXHM1*D(1,J,K) + SFX14*
     &                         (D(2,J,K) + D(IMM1,J,K))
      STORE(IM,J,K) = SFXHM1*D(IM,J,K) + SFX14*
     &                         (D(2,J,K)+D(IMM1,J,K))
 4004 CONTINUE
      J = JLED(NR)
      DO 4007 K=1,KM
      STORE(1,J,K)  = SFXHM1*D(1,J,K)+SFX14*
     &                      (D(1,J+1,K)+D(IM,J,K))
      STORE(IM,J,K) = SFXHM1*D(IM,J,K)+SFX14*
     &                       (D(IM,J+1,K)+D(1,J,K))
 4007 CONTINUE
 4006 CONTINUE
C*******************************************************************************
C     RE-SET THE VARIABLE   "D"  TO THE NEW SMOOTHED VALUE..
C
      DO 4015 NR = 1,NROWS
      J1 = JSTART(NR)
      J2 = JMIX(NR)
      DO 4015 J  = J1,J2
      DO 4015 K  = 1,KM
      DO 4015 I  = 1,IM
      D(I,J,K)   = STORE(I,J,K)
 4015 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C     RE ENTER HERE IF NO STREAMWISE SMOOTHING.
C
 4001 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C
      IF(SFT.LT.0.0001) GO TO 9500
C      
C          NOW PERFORM THE PITCHWISE AND SPANWISE SMOOTHING.
C
      DO 9000 J=2,JM
C  TFLOW ADDITION
      IF(IM.EQ.2) GO TO 9101
C
      DO 9100 K=1,KM
      DO 9110 I=2,IMM1
      AVG(I)    = FD(I)*D(I+1,J,K)+FU(I)*D(I-1,J,K)
      CURVE(I)  = D(I,J,K)-AVG(I)
 9110 CONTINUE
C
      IF(IND(J).EQ.1.OR.INDLE(J).EQ.1) GO TO 9112
C    SET VALUES ON THE PERIODIC BOUNDARIES
      AVG(1)    = 0.5*(D(2,J,K)+D(IMM1,J,K))
      AVG(IM)   = AVG(1)
      CURVE(1)  = D(1,J,K)-AVG(1)
      CURVE(IM) = CURVE(1)
      SCURVE(1) = 0.5*(CURVE(2)+CURVE(IMM1))
      SCURVE(IM)= SCURVE(1)  
      GO TO 9113
C
 9112 CONTINUE
C     SET VALUES ON THE BLADE SURFACES
      AVG(1)    = D(2,J,K)    + FU(1)*(D(2,J,K)-D(NUM3,J,K))
      AVG(IM)   = D(IMM1,J,K) + FU(IM)*(D(IMM1,J,K)-D(IMM2,J,K))
      CURVE(1)  = 0.0
      CURVE(IM) = 0.0
      SCURVE(1) = 0.0
      SCURVE(IM)= 0.0
C
 9113 CONTINUE
C   SMOOTH THE SECOND DERIVATIVE, CURVE, TO FORM SCURVE .
      DO 9120 I=2,IMM1
      SCURVE(I) = FU(I)*CURVE(I-1)+FD(I)*CURVE(I+1)
 9120 CONTINUE
C
C     APPLY THE FINAL PITCHWISE SMOOTHING.
      DO 9130 I=1,IM
      D(I,J,K) = SFT1*D(I,J,K) +SFT*(AVG(I) + FAC_4TH*SCURVE(I))
 9130 CONTINUE
C
 9100 CONTINUE
C  TFLOW ADDITION
 9101 CONTINUE
C*******************************************************************************
C     START TO APPLY THE SPANWISE SMOOTHING
C  Q3D  
      IF(KM.EQ.2) GO TO 9301
C  END Q3D
C
      DO 9300 I=1,IM
      DO 9200 K=2,KMM1
      AVG(K)   = FKD(K)*D(I,J,K+1) + FKU(K)*D(I,J,K-1)
      CURVE(K) = D(I,J,K)-AVG(K)
 9200 CONTINUE
C     FORM VALUES ON THE END WALLS
      AVG(1)   = D(I,J,2) +0.5*FKU(1)*(D(I,J,2)-D(I,J,3))
      AVG(KM)  = D(I,J,KMM1)+0.5* FKU(KM)*(D(I,J,KMM1)-D(I,J,KMM2))
      CURVE(1) = 0.0
      CURVE(KM)= 0.0
      SCURVE(1)= 0.0
      SCURVE(KM)=0.0
C     SMOOTH THE SECOND DERIVATIVE
      DO 9210 K=2,KMM1
      SCURVE(K) = FKU(K)*CURVE(K-1) + FKD(K)*CURVE(K+1)
 9210 CONTINUE
C
C    APPLY THE FINAL SPANWISE SMOOTHING
      DO 9230 K=1,KM
      D(I,J,K)  = SFT1*D(I,J,K) + SFT*(AVG(K)+FAC_4TH*SCURVE(K))
 9230 CONTINUE
C
 9300 CONTINUE
C
 9301 CONTINUE
C
C*******************************************************************************
C   SMOOTHING OF THE CORNER POINTS USING SFT .
C
      D(1,J,1)  = SFT1*D(1,J,1)  + SFTH*(D(2,J,1)+D(1,J,2))
      D(IM,J,1) = SFT1*D(IM,J,1) + SFTH*(D(IMM1,J,1)+D(IM,J,2))
      D(1,J,KM) = SFT1*D(1,J,KM) + SFTH*(D(1,J,KMM1)+D(2,J,KM))
      D(IM,J,KM)= SFT1*D(IM,J,KM)+ SFTH*(D(IMM1,J,KM)+D(IM,J,KMM1))
 9000 CONTINUE
C
C*******************************************************************************
C
C      RE ENTER HERE IF NO PITCHWISE OR SPANWISE SMOOTHING
 9500 CONTINUE
C
C*******************************************************************************
C*******************************************************************************
C      SMOOTH THE EXIT FLOW IF REQUIRED, NSFEXIT POINTS ARE SMOOTHED BY SFEXIT.
C
      IF(SFEXIT.LT.0.0001) GO TO 8700
C      
      JSSTART = JM + 1 - NSFEXIT
      DO 8502 J= JSSTART,JM
C   SMOOTH IN THE "I" DIRECTION
C
      IF(IM.EQ.2) GO TO 8505
C
      DO 8500 K=1,KM
      DO 8501 I=2,IMM1
      D(I,J,K) = SFEX1*D(I,J,K)+SFEXH*(D(I+1,J,K)+D(I-1,J,K))
 8501 CONTINUE
      D(1,J,K) = D(1,J,K) +SFEXIT*(2.*D(2,J,K)-D(NUM3,J,K)-D(1,J,K))
      D(IM,J,K)= D(IM,J,K)+SFEXIT*(2.*D(IMM1,J,K)-D(IMM2,J,K)-D(IM,J,K))
 8500 CONTINUE
C
 8505 CONTINUE
C
C  Q3D
      IF(KM.EQ.2) GO TO 8602
C  END Q3D    
C     SMOOTH IN THE "K" DIRECTION
      DO 8600 I=1,IM
      DO 8601 K=2,KMM1
      D(I,J,K) = SFEX1*D(I,J,K)+SFEXH*(D(I,J,K-1) +D(I,J,K+1))
 8601 CONTINUE
      D(I,J,1) = D(I,J,1) +SFEXIT*(2.*D(I,J,2)-D(I,J,3)-D(I,J,1))
      D(I,J,KM)= D(I,J,KM)+SFEXIT*(2.*D(I,J,KMM1)-D(I,J,KMM2)-D(I,J,KM))
 8600 CONTINUE
C
 8602 CONTINUE
C
 8502 CONTINUE
C
C
 8700 CONTINUE
C
C******************************************************************************
C     END OF EXIT FLOW SMOOTHING
C******************************************************************************
      RETURN
      END
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
      SUBROUTINE RE_DESIGN(NR,K,J1,J2,JLEROW,JTEROW) 
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  SMER(JD),FRACNEW(JD),BETANEW(JD),
     &           SLOPE(JD),THICKUP(JD),THICKLOW(JD),TKUP(JD),TKLOW(JD),
     &		 XSS(JD),RSS(JD),S_SS(JD),RELSPACE(JD),S_REL(JD),
     &           SPACE(JD),FRAC_CHORD(JD),XNEW(JD),RNEW(JD),
     &           SSDIST(JD)
C
C*******************************************************************************
C     SUBROUTINE TO DESIGN A NEW BLADE SECTION ON THE STREAM SURFACE INPUT BELOW.
C*******************************************************************************
C
      WRITE(6,*)
      WRITE(6,*) ' USING THE BLADE REDESIGN OPTION ON SECTION NUMBER ',K
      WRITE(6,*) ' J1, J2,  JLEROW, JTEROW = ', J1,J2,JLEROW,JTEROW
C
      JLEE   = J1+JLEROW-1
      JTEE   = J1+JTEROW-1
      JMROW  = J2 - J1 + 1
      NXTRAP = 4
C
      WRITE(6,*) ' ABSOLUTE VALUES OF JLEE and JTEE= ',JLEE,JTEE
      WRITE(6,*) ' JM FOR THE ROW BEING DESIGNED IS  ', JMROW
C
C      INPUT THE NEW STREAM SURFACE COORDINATES AND RELATIVE GRID SPACING.
C      MARK THE POSITION OF THE BLADE LEADING EDGE AND THE TRAILING EDGE.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*) N_SS, N_LE, N_TE
C
      READ(5,*)   DUMMY_INPUT
      DO N = 1,N_SS
      READ(5,*) XSS(N),RSS(N),RELSPACE(N)
      END DO
C  
C
C   INPUT THE NUMBER OF POINTS AT WHICH A NEW BLADE CAMBER LINE AND THICKNESS WILL DE INPUT.
C     ALSO THE NUMBER OF TIMES THEY WILL BE SMOOTHED.
C
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    NNEW, NSMOOTH
      WRITE(6,*)
C
C   INPUT THE NEW CAMBER LINE SLOPE AND UPPER AND LOWER TANGENTIAL THICKNESS.
C   AS FRACTIONS OF THE MERIDIONAL CHORD. 
C
      READ(5,*)    DUMMY_INPUT
      DO  NN = 1,NNEW
      READ(5,*)     FRACNEW(NN),BETANEW(NN),THICKUP(NN),THICKLOW(NN)
      WRITE(6,*)  ' FRACNEW(NN),BETANEW(NN),THICKUP(NN),THICKLOW(NN) ',
     &              FRACNEW(NN),BETANEW(NN),THICKUP(NN),THICKLOW(NN)
      END DO
      WRITE(6,*)
C
      READ(5,*)   DUMMY_INPUT
      READ(5,*)   FRAC_CHORD_UP, FRAC_CHORD_DWN, RTHETA_MID
      WRITE(6,*) 'FRAC_CHORD_UP, FRAC_CHORD_DWN, RTHETA_MID',
     &            FRAC_CHORD_UP, FRAC_CHORD_DWN, RTHETA_MID
      WRITE(6,*)
C
C    SMOOTH THE INPUT VALUES OF BLADE CAMBER LINE ANGLE AND THICKNESS.      
      FSMOOTH = 0.25
      CALL SMOOTH(1,NNEW,NSMOOTH,FSMOOTH,FRACNEW,BETANEW)
      CALL SMOOTH(1,NNEW,NSMOOTH,FSMOOTH,FRACNEW,THICKUP)
      CALL SMOOTH(1,NNEW,NSMOOTH,FSMOOTH,FRACNEW,THICKLOW)
C
C    FIND THE MERIDIONAL DISTANCE OF THE POINTS ON THE STREAM SURFACE.
      S_SS(1) = 0.0
      DO N = 2,N_SS
      XDIF = XSS(N) - XSS(N-1)
      RDIF = RSS(N) - RSS(N-1)
      S_SS(N) = S_SS(N-1) + SQRT(XDIF*XDIF + RDIF*RDIF)
      END DO
C
C   SET THE MERIDIONAL CHORD
      S_MERID  = S_SS(N_TE) - S_SS(N_LE)
C
C     SCALE THE THICKNESS BY THE MERIDIONAL CHORD
      DO NN = 1,NNEW 
           TKUP(NN)  = THICKUP(NN)*S_MERID
           TKLOW(NN) = THICKLOW(NN)*S_MERID
      END DO
C
C    NON-DIMENSIONALISE THE MERIDIONAL DISTANCE BY THE MERIDIONAL CHORD
C    S_REL(N_LE) = 0, S_REL(N_TE) = 1.
      DO N = 1,N_SS
      S_REL(N) = ( S_SS(N) - S_SS(N_LE) )/S_MERID
      END DO
C
C    INTERPOLATE TO FIND THE RELATIVE GRID SPACINGS AT THE FINAL GRID POINTS.
C    THIS IS FOR ALL GRID POINTS.
      DO J = 1,JMROW
      FRACJ = FLOAT(J-JLEROW)/FLOAT(JTEROW - JLEROW)
      CALL INTP(N_SS,S_REL,RELSPACE,FRACJ,SPACE(J))
      END DO
C
C    FIND THE FRACTION OF MERIDIONAL CHORD AT THE STREAM SURFACE POINTS INPUT
C    ON THE BLADE SURFACE.
      FRAC_CHORD(JLEROW) = 0.0
      DO J = JLEROW+1,JTEROW
      FRAC_CHORD(J) = FRAC_CHORD(J-1) + 0.5*(SPACE(J)+SPACE(J-1))
      END DO
C     MAKE SURE THAT FRAC_CHORD = 0 AT THE  LE  AND  = 1 AT THE TE.
      DO J = JLEROW,JTEROW
      FRAC_CHORD(J) = FRAC_CHORD(J)/FRAC_CHORD(JTEROW)
      END DO
C
C     SET THE FINAL MERIDIONAL POSITIONS OF THE GRID POINTS ON THE BLADE SURFACE
      DO J = JLEROW,JTEROW
      SSDIST(J) = S_SS(N_LE) + S_MERID*FRAC_CHORD(J)
      END DO
C
C    FIND THE FRACTION OF MERIDIONAL CHORD AT THE STREAM SURFACE POINTS INPUT
C    UPSTREAM OF THE LE.
      FRAC_CHORD(1) = 0.0
      DO J = 2,JLEROW
      FRAC_CHORD(J) = FRAC_CHORD(J-1) + 0.5*(SPACE(J)+SPACE(J-1))
      END DO

C     MAKE SURE THAT FRAC_CHORD = 0 AT THE  LE  AND  = 1 AT THE LEADING EDGE.
      DO J = 1,JLEROW
      FRAC_CHORD(J) = FRAC_CHORD(J)/FRAC_CHORD(JLEROW)
      END DO
C
C     SET THE FINAL MERIDIONAL POSITIONS OF THE GRID POINTS UPSTREAM OF THE LE.
      DO J = 1,JLEROW-1
      SSDIST(J) = S_SS(N_LE) + FRAC_CHORD_UP*S_MERID*(FRAC_CHORD(J)-1.0)
      END DO
C
C    FIND THE FRACTION OF MERIDIONAL CHORD AT THE STREAM SURFACE POINTS INPUT
C    DOWNSTREAM OF THE TRAILING EDGE.
      FRAC_CHORD(JTEROW) = 0.0
      DO J = JTEROW+1,JMROW
      FRAC_CHORD(J) = FRAC_CHORD(J-1) + 0.5*(SPACE(J)+SPACE(J-1))
      END DO

C     MAKE SURE THAT FRAC_CHORD = 0 AT THE  TE  AND  = 1 AT THE EXIT.
      SUM_FRAC = FRAC_CHORD(JMROW) - FRAC_CHORD(JTEROW)
      DO J = JTEROW,JMROW
      FRAC_CHORD(J) = (FRAC_CHORD(J) - FRAC_CHORD(JTEROW))/SUM_FRAC
      END DO
C
C     SET THE FINAL MERIDIONAL POSITIONS OF THE GRID POINTS DOWNSTREAM OF THE TE.
      DO J = JTEROW+1,JMROW
      SSDIST(J) = SSDIST(JTEROW) + FRAC_CHORD_DWN*S_MERID*FRAC_CHORD(J)
      END DO
C
C   INTERPOLATE TO FIND THE X  AND  R  COORDINATES OF THE FINAL GRID.
      DO J = 1,JMROW
      JALL  = J1 + J - 1
      CALL INTP(N_SS,S_SS,XSS,SSDIST(J),XSURF(JALL,K))
      CALL INTP(N_SS,S_SS,RSS,SSDIST(J),RSURF(JALL,K))
      END DO 
C
C   INTERPOLATE IN THE NEW BLADE TO OBTAIN THE COORDINATES AT THE GRID POINTS.
C 
      DO J = JLEROW,JTEROW
      JALL = J + J1 - 1
      FRACS = (SSDIST(J)-SSDIST(JLEROW))/(SSDIST(JTEROW)-SSDIST(JLEROW))
      CALL INTP(NNEW,FRACNEW, BETANEW,FRACS, SLOPE(JALL))
      CALL INTP(NNEW,FRACNEW, TKUP, FRACS, THICKUP(JALL))
      CALL INTP(NNEW,FRACNEW, TKLOW,FRACS, THICKLOW(JALL))
      END DO
C
C     FORM THE GRID UPSTREAM OF THE LEADING EDGE
      DTHDX  = TAN(SLOPE(JLEE+NXTRAP)*DEGRAD)/RSURF(JLEE+NXTRAP,K)
      DO J = J1,JLEE
      DISTUP        = SSDIST(J) - SSDIST(JLEROW)
      RT_UPP(J,K)   = DISTUP*DTHDX*RSURF(J,K)
      RT_THICK(J,K) = 0.0
      END DO
C
C    FORM THE GRID ON THE BLADE
      DTHDXP = TAN(SLOPE(JLEE)*DEGRAD)/RSURF(JLEE,K)
      TH_MID_P = RT_UPP(JLEE,K)/RSURF(JLEE,K)
      DO J   = JLEE + 1,JTEE
      DTHDX  = TAN(SLOPE(J)*DEGRAD)/RSURF(J,K)
      TH_MID = TH_MID_P + 0.5*(DTHDX+DTHDXP)*(SSDIST(J) - SSDIST(J-1) )
      RT_UPP(J,K)   = TH_MID*RSURF(J,K)  + THICKUP(J)
      RT_THICK(J,K) = THICKUP(J) + THICKLOW(J)
      DTHDXP        = DTHDX
      TH_MID_P      = TH_MID
      END DO
C
C     FORM THE NEW GRID DOWNSTREAM OF THE TRAILING EDGE
      DTHDX  = TAN(SLOPE(JTEE-NXTRAP)*DEGRAD)/RSURF(JTEE-NXTRAP,K)
      DO J = JTEE + 1, J2
      DISTDWN       = SSDIST(J) - SSDIST(JTEROW)
      RT_UPP(J,K)   = RT_UPP(JTEE,K)  + DISTDWN*DTHDX*RSURF(J,K)
      RT_THICK(J,K) = 0.0
      END DO                 
C
C     SET RTHETA = RTHETA_MID  AT THE MID J POINT OF THIS SECTION
C
      JMDD  = (JLEE + JTEE)/2
      THSHIFT = (RTHETA_MID
     &       - (RT_UPP(JMDD,K)- 0.5*RT_THICK(JMDD,K)))/RSURF(JMDD,K)
      DO J = J1,J2
      RT_UPP(J,K) = RT_UPP(J,K) + THSHIFT*RSURF(J,K)
      END DO
C
C    WRITE OUT THE NEW BLADE COORDINATES.
C
      WRITE(6,*)
      WRITE(6,*)'NEW BLADE COORDINATES, J   X , RT_UPP, RT_THICK, RSURF'
      DO J=J1,J2
      WRITE(6,666) J, XSURF(J,K),RT_UPP(J,K),RT_THICK(J,K),RSURF(J,K)
      END DO
  666 FORMAT(I5,4F20.5)
C
C*********************************************************************************** 
C     END OF BLADE SECTION REDESIGN
C***********************************************************************************     
C 
      RETURN
      END
C
C*****************************************************************************
C
C*****************************************************************************
      SUBROUTINE SET_SSTHICK(J1,J2)
C
      INCLUDE 'commall-open-18.3'
C
      DIMENSION FRACSS(20),TKSS(20),TK_SS(JD),SDIST(JD)
      SAVE      SRANGE1
C
      WRITE(6,*) '  INPUTTING THE DATA FOR A Q3D CALCULATION '
      WRITE(6,*)
      READ(5,*)    DUMMY_INPUT
      READ(5,*)    Q3DFORCE
      WRITE(6,*) ' Q3DFORCE = ', Q3DFORCE
      READ(5,*)    NSS
      WRITE(6,*) ' STREAM SURFACE DATA '
      READ(5,*)   (FRACSS(N), N=1,NSS)
      WRITE(6,*)  'FRACSS',( FRACSS(N), N=1,NSS)
      READ(5,*)   (TKSS(N),N=1,NSS)
      WRITE(6,*)  '  TKSS',( TKSS(N),N=1,NSS)
      WRITE(6,*)
C
      SDIST(J1) = 0.0
      DO 10 J = J1+1, J2
           XD = XSURF(J,1)  - XSURF(J-1,1)
           RD = RSURF(J,1)  - RSURF(J-1,1)
           SDIST(J) = SDIST(J-1) + SQRT(XD*XD+RD*RD)
   10 CONTINUE
      SRANGE   = SDIST(J2) - SDIST(J1)
      IF(J1.EQ.1) SRANGE1 = SRANGE
C
      DO 20  J = J1,J2
      SDIST(J) = (SDIST(J)-SDIST(J1))/SRANGE
   20 CONTINUE
C
      TKSS1 = TKSS(1)
      DO 25 N=1,NSS
      TKSS(N) = 0.05*TKSS(N)/TKSS1*SRANGE1
   25 CONTINUE
C   
      DO 30 J = J1,J2
      CALL INTP(NSS,FRACSS,TKSS,SDIST(J),TK_SS(J+1-J1))
   30 CONTINUE
C
C  SMOOTH THE STREAM SURFACE THICKNESS
      NVAL = J2+1-J1
      CALL SMOOTH(1,NVAL,5,0.25,SDIST,TK_SS)
C
      TK_SS(1)     = TK_SS(2)
      TK_SS(NVAL)  = TK_SS(NVAL-1) 
C
      DO 40 J=J1,J2
      JP1 = J+1
      IF(J.EQ.J2) JP1 = J2
      JM1 = J-1
      IF(J.EQ.J1) JM1 = J1
      DRDS = (RSURF(JP1,1)
     &     -  RSURF(JM1,1))/( SRANGE*(SDIST(JP1) - SDIST(JM1)) )
      DXDS = (XSURF(JP1,1)
     &     -  XSURF(JM1,1))/( SRANGE*(SDIST(JP1) - SDIST(JM1)) )
      RSURF(J,2)    = RSURF(J,1) + TK_SS(J+1-J1)*DXDS
      XSURF(J,2)    = XSURF(J,1) - TK_SS(J+1-J1)*DRDS
      RT_UPP(J,2)   = RT_UPP(J,1)*RSURF(J,2)/RSURF(J,1)
      RT_THICK(J,2) = RT_THICK(J,1)
   40 CONTINUE  
C      
      RETURN
      END
C************************************************************************************
C************************************************************************************
C************************************************************************************
C
      SUBROUTINE RESTAGGER(K,J1,J2,JLEROW,JTEROW,ROTATE,FRACX_ROT)
C
C************************************************************************************
C     THIS SUBROUTINE RESTAGGERS A BLADE SECTION.
C     THE ROTATION IS APPLIED ON A STREAM SURFACE SO THAT IT IS APPLICABLE
C     TO BOTH RADIAL FLOW AND AXIAL FLOW MACHINES.
C************************************************************************************
      INCLUDE 'commall-open-18.3'
C
      DIMENSION  SDIS_T(JD),THET_A(JD),SLOPE(JD)
C
      ROTATE     = 0.0
      FRACX_ROT  = 0.5
      READ(5,*) DUMMY_INPUT
      READ(5,*,ERR=1502) ROTATE, FRACX_ROT
 1502 CONTINUE
      WRITE(6,*) ' INPUT FOR RESTAGGERING OPTION'
      WRITE(6,*) ' ROTATION =', ROTATE, 'ABOUT ',FRACX_ROT 
      WRITE(6,*)
C
      RAD_ROT = ROTATE*DEGRAD 
C
C  CALCULATE THE MERIDIONAL DISTANCES AND ANGLES.
C
      SDIS_T(J1) = 0.0  
      THET_A(J1) =  (RT_UPP(J1,K) - 0.5*RT_THICK(J1,K))/RSURF(J1,K)   
      DO 10 J = J1+1,J2
      XDIF = XSURF(J,K) - XSURF(J-1,K)
      RDIF = RSURF(J,K) - RSURF(J-1,K)
      SDIS_T(J) = SDIS_T(J-1) + SQRT(XDIF*XDIF + RDIF*RDIF)
      THET_A(J) = (RT_UPP(J,K)- 0.5*RT_THICK(J,K))/RSURF(J,K) 
   10 CONTINUE
C
C     FIND THE J VALUE AT THE CENTRE OF ROTATION, JROT .
C
      JROT = 1
      JLEE = J1+JLEROW-1
      JTEE = J1+JTEROW-1
      S_MERID = SDIS_T(J2)
      DO 20 J = J1,J2
      FRACRD = (SDIS_T(J) - SDIS_T(JLEE))/ (SDIS_T(JTEE) - SDIS_T(JLEE))
      IF(FRACRD.GT.FRACX_ROT)  GO TO 30
   20 CONTINUE
   30 CONTINUE
      JROT = J
      RTROT = RT_UPP(JROT,K)- 0.5*RT_THICK(JROT,K)
C
C     CALCULATE THE CENTRE LINE SLOPES ON THE STREAM SURFACE
C     DEFINED AS AS  " R D(theta)/DS  "
C       
      DO 40 J = J1+1,J2-1
      R_DT_DS = RSURF(J,K)*(THET_A(J+1) - THET_A(J-1) )
     &         / (SDIS_T(J+1) - SDIS_T(J-1) )
      SLOPE(J) = ATAN(R_DT_DS)
   40 CONTINUE
      SLOPE(J1) = SLOPE(J1+1)
      SLOPE(J2) = SLOPE(J2-1)
C
C     ROTATE THE BLADE ON THE STREAM SURFACE SO THAT ALL CENTRE LINE
C     SLOPES CHANGE BY "ROTATE" DEGREES IN A CLOCKWISE SENSE.
C     MAINTAIN THE SAME PERPENDICULAR THICKNESS.
C
      THMID   = THET_A(J1)
      DO 50 J = J1,J2
      JM1 = J-1
      IF(J.EQ.1) JM1 = 1
      AVGSLOPE = 0.5*(SLOPE(J) + SLOPE(JM1)) - RAD_ROT
      ANGNEW   = SLOPE(J) - RAD_ROT
      RT_THICK(J,K) = RT_THICK(J,K)*COS(SLOPE(J))/COS(ANGNEW)
      THMID = THMID 
     &      + (SDIS_T(J) - SDIS_T(JM1))*TAN(AVGSLOPE)/RSURF(JM1,K)
      RT_UPP(J,K) = RSURF(J,K)*THMID + 0.5*RT_THICK(J,K)
   50 CONTINUE
C
C      ROTATE AS A SOLID BODY SO THAT THE CENTRE OF ROTATION, JROT, DOES NOT MOVE.
C 
      THSHIFT = (RTROT
     &        - (RT_UPP(JROT,K)- 0.5*RT_THICK(JROT,K)))/RSURF(JROT,K)
      DO 60 J = J1,J2
      RT_UPP(J,K) = RT_UPP(J,K) + THSHIFT*RSURF(J,K)
   60 CONTINUE
C
      RETURN
      END
C
C*********************************************************************************
C*********************************************************************************
C*********************************************************************************
      SUBROUTINE WALLFUN(I,J,K,NWALL,PERP,DPDS,DENSITY,TWALL,
     &                   YPLUS_OLD,WREL_WALL,YPLUS_NEW)
C
C   THIS SUBROUTINE APPLIES THE WALL FUNCTIONS PROPOSED BY Shih et al. NASA/TM-1999-209398 
C
      INCLUDE  'commall-open-18.3'
C
C    DO NOT USE THE PRESSURE GRADIENT TERM UNLESS YPLUSWALL IS LESS THAN -10.
      IF(YPLUSWALL.GT.-10.0) GO TO 100
C
      IF(ABS(DPDS).LT.1.0) DPDS = 1.0
      SGN_DPDS = ABS(DPDS)/DPDS
C
      VPRES    = (VISLAM/DENSITY/DENSITY* ABS(DPDS))**0.33333
      VSTAR    = YPLUS_OLD*VISLAM/(DENSITY*PERP)
C
C     SET LIMITS TO THE PRESSURE GRADIENT VELOCITY TERM.
      IF(SGN_DPDS.LT.0.0.AND.VPRES.GT.0.33*VSTAR) VPRES = 0.33*VSTAR
      IF(SGN_DPDS.GT.0.0.AND.VPRES.GT.2.0*VSTAR)  VPRES = 2.0*VSTAR
C 
      YPLUS_PRES = VPRES*PERP*DENSITY/VISLAM
C  
C     SET THE PRESSURE GRADIENT TERM FROM YPLUS_PRES.
      YP1     = YPLUS_PRES
      YP2     = YP1*YP1
      YP3     = YP2*YP1
      YP4     = YP2*YP2
      IF(YP1.LE.4.0)                  FYPLUSP = 0.5*YP2 - 7.31E-3*YP3
      IF(YP1.GT.4.0.AND.YP1.LE.15.0)  FYPLUSP = -15.138 + 8.4688*YP1
     &                 - 0.81976*YP2 + 3.7292E-2*YP3 - 6.3866E-4*YP4
      IF(YP1.GT.15.0.AND.YP1.LE.30.0) FYPLUSP = 11.925 + 0.934*YP1
     &             - 2.7805E-2*YP2 + 4.6262E-4*YP3 - 3.1442E-6*YP4
      IF(YP1.GT.30.0)                 FYPLUSP = 5*ALOG(YP1) + 8.0
C
  100 CONTINUE
C
C    RE-ENTER HERE IF YPLUSWALL IS BETWEEN -5 and -10.
C    SET THE WALL DISTANCE-VELOCITY  TERM FROM YPLUS_OLD
      YP1      = YPLUS_OLD
      YP2      = YP1*YP1
      YP3      = YP2*YP1
      YP4      = YP2*YP2 
      IF(YP1.LE.5.0) FYPLUS = YP1 + 1.0E-2*YP2 -2.9E-03*YP3
      IF(YP1.GT.5.0.AND.YP1.LE.30.0) FYPLUS = -0.872 + 1.465*YP1
     &                 -7.02E-2*YP2 + 1.66E-03*YP3 - 1.495E-5*YP4
      IF(YP1.GT.30.AND.YP1.LE.140.0)  FYPLUS = 8.6 + 0.1864*YP1
     &             - 2.006E-3*YP2 + 1.144E-5*YP3 -2.551E-8*YP4
      IF(YP1.GT.140.0)  FYPLUS   = 2.439*ALOG(YP1) + 5.0
C
C
      IF(YPLUSWALL.GT.-10.0) THEN
           FRACP = 0.0
      ELSE
           FRACP = SGN_DPDS*(VPRES/VSTAR) * FYPLUSP/FYPLUS
      END IF
C
      VPLUS     = FYPLUS*(1.0 + FRACP)
C
      VSTAR     = WREL_WALL/VPLUS
C      
      TWALL     = DENSITY*VSTAR*VSTAR
C
C  SET THE NEW VALUE OF YPLUS. THIS IS USED ON THE NEXT CALL OF THIS SUBROUTINE
C  AND THE ITERATION  SHOULD CONVERGE TO THE FINAL VALUE.
      YPLUS_NEW = DENSITY*VSTAR*PERP/VISLAM
C
      RETURN
      END
C
C********************************************************************************* 
C*********************************************************************************
C      
      SUBROUTINE SET_PWALLGRAD
C
C    THIS SUBROUTINE EVALUATRES THE PRESSURE GRADIENTS ON THW WALLS FOR USE WITH
C    THE Shih et al WALL FUNCTIONS
C
      INCLUDE  'commall-open-18.3'
C
      DIMENSION  TEMPP(JD)
C
C   THE LIMIT ON PRESSURE GRADIENT IS IF AVERAGE PRESSURE CHANGE PER BLADE ROW
C   OCCURS OVER 10% OF THE CHORD.  
      DPDS_LIM = 20.0*(PO1(KMID) - 0.5*(PDOWN_HUB+PDOWN_TIP)/NROWS/
     &                (CHORD(1)  + CHORD(NROWS)) )
      DPDS_LIM =  ABS(DPDS_LIM)
C
C     CALCULATE THE PRESSURE GRADIENT ALONG A STREAMLINE
      DO 226 J = 2,JM
      DO 226 K = 1,KMM1
      DO 226 I = 1,IMM1
      CALL GRADVEL(I,J,K,PEFF,DPDX,DPDR,DPDT)
      VXCELL = VX(I,J,K)+VX(I,J-1,K)+VX(I+1,J,K)+VX(I+1,J-1,K) 
     &       + VX(I,J,K+1)+VX(I,J-1,K+1)+VX(I+1,J,K+1)+VX(I+1,J-1,K+1)
      VRCELL = VR(I,J,K)+VR(I,J-1,K)+VR(I+1,J,K)+VR(I+1,J-1,K)
     &       + VR(I,J,K+1)+VR(I,J-1,K+1)+VR(I+1,J,K+1)+VR(I+1,J-1,K+1)
      WTCELL = WT(I,J,K)+WT(I,J-1,K)+WT(I+1,J,K)+WT(I+1,J-1,K)
     &       + WT(I,J,K+1)+WT(I,J-1,K+1)+WT(I+1,J,K+1)+WT(I+1,J-1,K+1)
      WREL_CELL = SQRT(VXCELL*VXCELL + VRCELL*VRCELL + WTCELL*WTCELL)
      DPDS      = (DPDX*VXCELL+DPDR*VRCELL+DPDT*WTCELL)/WREL_CELL
      IF(ABS(DPDS).GT.DPDS_LIM) DPDS = DPDS_LIM*DPDS/ABS(DPDS)
      DPDS_CELL(I,J,K) = DPDS      
  226 CONTINUE
C 
C   SMOOTH THE PRESSURE GRADIENT ALONG A GRID LINE.
      DO 250 NR = 1,NROWS
C
      J1 = JSTART(NR)
      J2 = JMIX(NR)
      DO 252 K=1,KMM1
      DO 252 I=1,IMM1
C
      DPDS_CELL(I,J2+1,K) = 0.0
      DPDS_CELL(I,J2+2,K) = 2*DPDS_CELL(I,J2+3,K) - DPDS_CELL(I,J2+4,K)
      DPDS_CELL(I,J2,K)   = 2*DPDS_CELL(I,J2-1,k) - DPDS_CELL(I,J2-2,K)
C
      DO 251 J = J1+2,J2-1
      TEMPP(J) = 0.5*DPDS_CELL(I,J,K)
     &         + 0.25*(DPDS_CELL(I,J-1,K) + DPDS_CELL(I,J+1,K))  
  251 CONTINUE 
      TEMPP(J2+1) = 0.0
      DO 253 J = J1+3,J2-2
      DPDS_CELL(I,J,K)= 0.5*TEMPP(J) + 0.25*(TEMPP(J-1)+TEMPP(J+1))
  253 CONTINUE
C
      DPDS_CELL(I,J2+1,K) = 0.0
      DPDS_CELL(I,J2+2,K) = 2*DPDS_CELL(I,J2+3,K) - DPDS_CELL(I,J2+4,K)
      DPDS_CELL(I,J2,K)   = 2*DPDS_CELL(I,J2-1,k) - DPDS_CELL(I,J2-2,K)
C
  252 CONTINUE
C 
  250 CONTINUE 
C
      RETURN
      END
C******************************************************************************
C****************************************************************************** 
C
      SUBROUTINE COOL_INPUT
C
C  THIS SUBROUTINE READS IN THE DATA FORANY COOLING FLOWS.
C
      INCLUDE  'commall-open-18.3'
C
  99  FORMAT(A72)
C
      NCOOLB     = 0
      NCOOLW     = 0
      NCWLBLADEP = 0
      NCWLWALLP  = 0
      PI_180     = 3.1415926/180
      PI_30      = 3.1415926/30
C******************************************************************************
C******************************************************************************
      DO 100 N_ROW = 1,NROWS
C
      J1      =  JSTART(N_ROW)
C
      WRITE(6,*) ' READING THE COOLING FLOW DATA FOR ROW No.', N_ROW
      READ(5,99)  DUMMY_INPUT
      WRITE(6,99) DUMMY_INPUT
C
      READ(5,*)  NCWLBLADE, NCWLWALL
      WRITE(6,*) NCWLBLADE, NCWLWALL
C
      IF(NCWLBLADE.EQ.0.AND.NCWLWALL.EQ.0) GO TO 100
C
      NCOOLB  = NCWLBLADEP + NCWLBLADE
      NCOOLW  = NCWLWALLP  + NCWLWALL
C
C******************************************************************************
C******************************************************************************
C     READ IN THE COOLING FLOW DETAILS FOR THE BLADE SURFACES.
C
      IF(NCWLBLADE.NE.0) THEN
C
      WRITE(6,*) ' READING THE BLADE SURFACE COOLING DATA FOR ROW No.',
     &             N_ROW
C
      READ(5,99)  DUMMY_INPUT
      WRITE(6,99) DUMMY_INPUT
C
      DO 200 NBSURF = NCWLBLADEP+1,NCOOLB
      READ(5,*)  IC(NBSURF),JCBS(NBSURF),JCBE(NBSURF),KCBS(NBSURF),
     &           KCBE(NBSURF)
      READ(5,*)  CFLOWB(NBSURF),TOCOOLB(NBSURF),POCOOLB(NBSURF),
     &           MACHCOOL(NBSURF),SANGLEB(NBSURF),XANGLEB(NBSURF),
     &           RVT_IN_B(NBSURF),RPM_COOL
      WRITE(6,*) IC(NBSURF),JCBS(NBSURF),JCBE(NBSURF),KCBS(NBSURF),
     &           KCBE(NBSURF)
      WRITE(6,*) CFLOWB(NBSURF),TOCOOLB(NBSURF),POCOOLB(NBSURF),
     &           MACHCOOL(NBSURF),SANGLEB(NBSURF),XANGLEB(NBSURF),
     &           RVT_IN_B(NBSURF),RPM_COOL
C
      SANGLEB(NBSURF)    = SANGLEB(NBSURF)*PI_180
      XANGLEB(NBSURF)    = XANGLEB(NBSURF)*PI_180
      WRAD_COOLB(NBSURF) = RPM_COOL*PI_30
      JCBS(NBSURF)       = JCBS(NBSURF) + J1 -1
      JCBE(NBSURF)       = JCBE(NBSURF) + J1 -1
C
      IF(TOCOOLB(NBSURF).LT.TCOOL_MIN) TCOOL_MIN = TOCOOLB(NBSURF)
C
  200 CONTINUE
C
      END IF
C
C      JCBS & JCBE are the J values where cooling starts and ends for
C      this patch, the J indices are defined relative to the start
C      of the current blade row. ie the upstream mixing plane is J = 1.
C
C      CFLOWB  is the coolant mass flow through this patch in Kg/s.
C
C      TOCOOLB is the ABSOLUTE stagnation temperature of the coolant when
C      it is first supplied to the blade row. (NOTE ABSOLUTE).
C
C      MACHCOOL is the RELATIVE Mach number of the coolant as it leaves the
C      blade and enters the mainstream flow. (NOTE RELATIVE).
C
C      SANGLEB is the angle between the coolant jet and the blade surface.
C
C      XANGLEB is the angle between the projection of the cooling jet onto
C      the blade surface and the intersection of the blade surface with
C      a surface of constant radius (cylindrical surface).
C
C      RVT_IN_B  is the angular momentum with which the cooling flow is supplied
C      disk chamber before entering the blades.
C
C      RPM_COOL is the rotational speed of the surface through which the
C      coolant is ejected. The coolant velocity and directions are relative
c      to this rotational speed.
C******************************************************************************
C******************************************************************************
C     READ IN THE COOLING FLOW DETAILS FOR THE ENDWALL SURFACES.
C******************************************************************************
C******************************************************************************
C
      IF(NCWLWALL.NE.0) THEN
C
      WRITE(6,*) ' READING THE ENDWALL COOLING DATA FOR ROW No.', N_ROW
C
      READ(5,99)  DUMMY_INPUT
      WRITE(6,99) DUMMY_INPUT
C
      DO 300 NWSURF = NCWLWALLP+1,NCOOLW
      READ(5,*)  KC(NWSURF),JCWS(NWSURF),JCWE(NWSURF),ICWS(NWSURF),
     &           ICWE(NWSURF)
      READ(5,*)  CFLOWW(NWSURF),TOCOOLW(NWSURF),POCOOLW(NWSURF),
     &           MACHCOOL(NWSURF),SANGLEW(NWSURF),TANGLEW(NWSURF),
     &           RVT_IN_W(NWSURF),RPM_COOL
      WRITE(6,*) KC(NWSURF),JCWS(NWSURF),JCWE(NWSURF),ICWS(NWSURF),
     &           ICWE(NWSURF)
      WRITE(6,*) CFLOWW(NWSURF),TOCOOLW(NWSURF),POCOOLW(NWSURF),
     &           MACHCOOL(NWSURF),SANGLEW(NWSURF),TANGLEW(NWSURF),
     &           RVT_IN_W(NWSURF),RPM_COOL
C
      SANGLEW(NWSURF)    = SANGLEW(NWSURF)*PI_180
      TANGLEW(NWSURF)    = TANGLEW(NWSURF)*PI_180
      WRAD_COOLW(NWSURF) = RPM_COOL*PI_30
      JCWS(NWSURF)       = JCWS(NWSURF) + J1 -1
      JCWE(NWSURF)       = JCWE(NWSURF) + J1 -1
C
      IF(TOCOOLW(NWSURF).LT.TCOOL_MIN) TCOOL_MIN = TOCOOLW(NWSURF)
C
  300 CONTINUE
C
      END IF
C
C     ASSUME THE MAXI<UM COOLANT EJECTION MACH NUMBER = 1.0 .
      TCOOL_MIN = TCOOL_MIN/(1.+0.5*(GA-1.0))
C   END OF ENDWALL COOLING FLOW DATA INPUT.

C
C      JCWS & JCWE are the J values where cooling starts and ends for
C      this patch, the J indices are defined relative to the start
C      of the current blade row. ie the upstream mixing plane is J = 1.
C
C      CFLOWW is the coolant mass flow through this patch in Kg/s.
C
C      TOCOOLW is the ABSOLUTE stagnation temperature of the coolant when
C      it is first supplied to the blade row. (NOTE ABSOLUTE).
C
C      MACHCOOL is the RELATIVE Mach number of the coolant as it leaves the
C      endwall and enters the mainstream flow. (NOTE RELATIVE).
C
C      SANGLEW is the angle between the coolant jet and the endwall surface.
C
C      TANGLEW is the angle between the projection of the cooling jet onto
C      the endwall surface and the intersection of the endwall surface with
C      the axial- radial plane ( i.e with the plane theta = constant).
C
C      RVT_IN_W  is the angular momentum with which the cooling flow is supplied
C      disk chamber before entering the blades.
C
C      RPM_COOL is the rotational speed of the surface through which the
C      coolant is ejected. The coolant velocity and directions are relative
c      to this rotational speed.
C
C

C
C******************************************************************************
C******************************************************************************
C
      NCWLBLADEP = NCOOLB
      NCWLWALLP  = NCOOLW
C
      WRITE(6,*)'BLADE ROW NUMBER ', N_ROW
      WRITE(6,*)'NUMBER OF COOLING PATCHES THIS ROW,ON BLADE,ON WALL =',
     &           NCWLBLADE, NCWLWALL
      WRITE(6,*)'TOTAL COOLING PATCHES UP TO NOW, NCOOLB,NCOOLW = ',
     &           NCOOLB, NCOOLW
C
  100 CONTINUE
C
C    END OF COOLANT DATA INPUT LOOP.
C*************************************************************************
C*************************************************************************
      WRITE(6,*) ' COOLING FLOW DATA INPUT COMPLETE. '
      WRITE(6,*) ' TOTAL NUMBER OF BLADE COOLING PATCHES   = ',NCOOLB
      WRITE(6,*) ' TOTAL NUMBER OF ENDWALL COOLING PATCHES = ',NCOOLW
C
C*************************************************************************
C*************************************************************************
      RETURN
      END
C*************************************************************************
C*************************************************************************
C*************************************************************************
C
      SUBROUTINE COOLIN_2
C
C  THIS SUBROUTINE SETS THE COOLING FLOWS BASED ON THE INPUT COOLANT STAGNATION TEMPERATURE
C  AND PRESSURE AND THE LOCAL STATIC PRESSURE AT THE POINT OF EJECTION. THE COOLING FLOWS
C  THEREFORE CHANGE DURING CONVERGENCE AND MUST BE RESET EVERY FEW ITERATIONS. THE COOLANT
C  EJECTION MACH NUMBER IS INPUT BUT NOT USED.

      INCLUDE  'commall-open-18.3'
C
      DIMENSION HO_ABS(JD,MAXKI),VNIN(JD,MAXKI),VSIN(JD,MAXKI),
     &          CELL_FLOW(JD,MAXKI)
C
      DO N = 1,NSTAGES
           WPUMP(N) = 0.0
      END DO
C
C*************************************************************************
C*************************************************************************
C
C      SET THE INITIAL COOLING FLOWS TO ZERO AND
C      WORK OUT AREA MAGNITUDES FOR COOLANT FLOWS
C
      IF(NCOOLB.EQ.0) GO TO 401
      DO 400 J=2,JM
      DO 400 K=1,KMM1
      CFLOWI1(J,K) = 0.0
      CFLOWIM(J,K) = 0.0
      ATOTI1(J,K) = SQRT(ABX(1,J,K)*ABX(1,J,K) + ABR(1,J,K)*ABR(1,J,K)
     &                  + ABT(J,K)*ABT(J,K))
      ATOTIM(J,K) = SQRT(ABX(IM,J,K)*ABX(IM,J,K)+ABR(IM,J,K)*ABR(IM,J,K)
     &                  + ABT(J,K)*ABT(J,K))
  400 CONTINUE
  401 CONTINUE
C
C
      IF(NCOOLW.EQ.0) GO TO 501
      DO 500 J=2,JM
      DO 500 I=1,IMM1
      CFLOWK1(I,J) = 0.0
      CFLOWKM(I,J) = 0.0
      ATOTK1(I,J)= SQRT(ASX(I,J,1)*ASX(I,J,1)  + ASR(I,J,1)*ASR(I,J,1))
      ATOTKM(I,J)= SQRT(ASX(I,J,KM)*ASX(I,J,KM)+ASR(I,J,KM)*ASR(I,J,KM))
  500 CONTINUE
  501 CONTINUE
C
C******************************************************************************
C******************************************************************************
C
C      NOW EVALUATE THE COOLING MASS FLOWS THROUGH THE BLADE SURFACES.
C      AND THE ASSOCIATED ENERGY AND MOMENTUM FLUXES.
C
C******************************************************************************
C******************************************************************************
      DO 600 NCB = 1,NCOOLB
C
      IF(IFGAS.EQ.0) THEN
           CPGAS       = CP
           GAGAS       = GA
           FGAGAS      = FGA
      ELSE
           DELT        = TOIN - TREF
           CPGAS       = CP1  + CP2*DELT + CP3*DELT*DELT
           GAGAS       = CPGAS/(CPGAS - RGAS)
           FGAGAS      = (GAGAS - 1.0)/GAGAS
      END IF 
C
C
      TOIN        = TOCOOLB(NCB)
      POIN        = POCOOLB(NCB)
      W_PRESWIRL  = WRAD_COOLB(NCB)*RVT_IN_B(NCB)
      SSANGLEB    = SIN(SANGLEB(NCB))
      CSANGLEB    = COS(SANGLEB(NCB))
      PATCH_FLOW  = 0.0
C
C***********************************************************************
C***********************************************************************
C
      DO 620 J = JCBS(NCB)+1,JCBE(NCB)
             N_ROW    = NROW(J)
             N_STAGE  = NSTAGE(N_ROW)
      DO 620 K = KCBS(NCB),KCBE(NCB)-1
C 
      IF(IC(NCB).EQ.1)
     &     PAVG_IN  = 0.25*(P(1,J,K)+P(1,J-1,K)+P(1,J,K+1)+P(1,J-1,K+1))
      IF(IC(NCB).EQ.IM)
     &     PAVG_IN  = 0.25*(P(IM,J,K)+P(IM,J-1,K)+P(IM,J,K+1)
     &                +P(IM,J-1,K+1))
C
      RAVG        = RAVG_CELL(J,K)
      VBLADE      = WRAD_COOLB(NCB)*RAVG
      TOREL       = TOIN  + (0.5*VBLADE*VBLADE - W_PRESWIRL)/CPGAS
      TOABS       = TOREL + 0.5*VBLADE*VBLADE/CPGAS
      PO_REL      = POIN*(TOREL/TOIN)**(1./FGAGAS)
      PO_EJECTB(NCB) = PO_REL
      PCOOL_RAT    = PAVG_IN/PO_REL
      IF(PCOOL_RAT.LT.0.5)    PCOOL_RAT = 0.5
      IF(PCOOL_RAT.GT.0.9999) PCOOL_RAT = 0.9999
      PAVG_IN     = PCOOL_RAT*PO_REL
      TCOOL_RAT   = PCOOL_RAT**FGAGAS
      TIN         = TOREL*TCOOL_RAT
      ROIN        = PAVG_IN/RGAS/TIN           
      VIN         = SQRT(2.*CPGAS*(TOREL - TIN))
      HO_ABS(J,K) = CPGAS*TOABS
      TORELB(NCB) = TOREL
      TSTAT_EJECTB(NCB) = TIN
C
C   EVALUATE THE NORMAL AND TANGENTIAL COOLANT RELATIVE VELOCITIES
      VNIN(J,K)     = VIN*SSANGLEB
      VSIN(J,K)     = VIN*CSANGLEB
C
C      EVALUATE THE  MASS FLOW THROUGH THE PATCH IF IT FILLS THE WHOLE CELL SURFACE AREA.
C
      IF(IC(NCB).EQ.1)  CELL_FLOW(J,K) = VNIN(J,K)*ROIN*ATOTI1(J,K)
      IF(IC(NCB).EQ.IM) CELL_FLOW(J,K) = VNIN(J,K)*ROIN*ATOTIM(J,K)
C   CALCULATE THE TOTAL FLOW THROUGH THE PATCH IF THE FULL AREA IS USED.
      PATCH_FLOW = PATCH_FLOW + CELL_FLOW(J,K)*NBLADE(J) 
C
C   END OF  J  , K  LOOPS
  620 CONTINUE
C
C    EVALUATE THE FRACTION OF CELL AREA USED BY THE COOLING FLOW IN ORDER TO OBTAIN THE SPECIFIED FLOW.
       FRAC_FLOW = CFLOWB(NCB)/PATCH_FLOW
C******************************************************************************
C******************************************************************************
C
C   NOW RECALCULATE THE COOLANT FLOWS SO THAT THEY SUM TO THE SPECIFIED TOTAL FLOW
      DO 650 J = JCBS(NCB)+1,JCBE(NCB)
      DO 650 K = KCBS(NCB),KCBE(NCB)-1
C
C   SET THE COOLANT FLOW THROUGH THE CELL FACES SO THAT THEY SUM TO THE SPECIFIED FLOW.
      CFLOWIN        = CELL_FLOW(J,K) * FRAC_FLOW
C
      RAVG           = RAVG_CELL(J,K)
      VBLADE         = WRAD_COOLB(NCB)*RAVG
      WPUMP(N_STAGE) = WPUMP(N_STAGE) + NBLADE(J)*
     &                 CFLOWIN*(VBLADE*VBLADE - W_PRESWIRL)
C
      IF(IC(NCB).EQ.1)  THEN      
           ANORM    = ATOTI1(J,K)
           XNORM    = ABX(1,J,K)/ANORM
           TNORM    = ABT(J,K)/ANORM
           RNORM    = ABR(1,J,K)/ANORM
      ENDIF
      IF(IC(NCB).EQ.IM) THEN
           ANORM    = ATOTIM(J,K)
           XNORM    = ABX(IM,J,K)/ANORM
           TNORM    = ABT(J,K)/ANORM
           RNORM    = ABR(IM,J,K)/ANORM
      END IF
C
      ANG      =  XANGLEB(NCB)
      XTNORM   =  SQRT(XNORM*XNORM + TNORM*TNORM)
      SX       =  (TNORM*COS(ANG) - RNORM*XNORM*SIN(ANG))/XTNORM
      ST       = -(XNORM*COS(ANG) + RNORM*TNORM*SIN(ANG))/XTNORM
      SR       =  SIN(ANG)*XTNORM
C
      IF(IC(NCB).EQ.1)  THEN
           XVEL          =  VNIN(J,K)*XNORM + VSIN(J,K)*SX
           RVEL          =  VNIN(J,K)*RNORM + VSIN(J,K)*SR
           TVEL          =  VNIN(J,K)*TNORM + VSIN(J,K)*ST
           CFLOWI1(J,K)  = CFLOWIN
           HOCWLI1(J,K)  = CFLOWIN*(HO_ABS(J,K) + VBLADE*TVEL)
           VXCWLI1(J,K)  = CFLOWIN*XVEL
           VRCWLI1(J,K)  = CFLOWIN*RVEL
           RVTCWLI1(J,K) = CFLOWIN*(VBLADE  +  TVEL)*RAVG
      ELSE
           XVEL          = -VNIN(J,K)*XNORM + VSIN(J,K)*SX
           RVEL          = -VNIN(J,K)*RNORM + VSIN(J,K)*SR
           TVEL          = -VNIN(J,K)*TNORM + VSIN(J,K)*ST
           CFLOWIM(J,K)  = CFLOWIN
           HOCWLIM(J,K)  = CFLOWIN*(HO_ABS(J,K) + VBLADE*TVEL)
           VXCWLIM(J,K)  = CFLOWIN*XVEL
           VRCWLIM(J,K)  = CFLOWIN*RVEL
           RVTCWLIM(J,K) = CFLOWIN*(VBLADE  +  TVEL)*RAVG
      END IF
C   END OF THE  J,K  LOOP
  650 CONTINUE
C
C******************************************************************************
C     END OF SETTING THE BLADE SURFACE COOLING FLUXES.
  600 CONTINUE
C
C**************************************************************************
C**************************************************************************
C**************************************************************************
C**************************************************************************
C**************************************************************************
C
C      NOW EVALUATE THE COOLING MASS FLOWS THROUGH THE ENDWALL SURFACES.
C      AND THE ASSOCIATED ENERGY AND MOMENTUM FLUXES.
C
      DO 700 NCW = 1,NCOOLW
C
      IF(IFGAS.EQ.0) THEN
           CPGAS       = CP
           GAGAS       = GA
           FGAGAS      = FGA
      ELSE
           DELT        = TOIN-TREF
           CPGAS       = CP1 + CP2*DELT + CP3*DELT*DELT
           GAGAS       = CPGAS/(CPGAS - RGAS)
           FGAGAS      = (GAGAS - 1.0)/GAGAS
      END IF 
C
C
      TOIN        = TOCOOLW(NCW)
      POIN        = POCOOLW(NCW)
      W_PRESWIRL  = WRAD_COOLW(NCW)*RVT_IN_W(NCW)
      SSANGLEW    = SIN(SANGLEW(NCW))
      CSANGLEW    = COS(SANGLEW(NCW))
      PATCH_FLOW  = 0.0
C
C***********************************************************************
C***********************************************************************
C
      DO 720 J = JCWS(NCW)+1,JCWE(NCW)
             N_ROW    = NROW(J)
             N_STAGE  = NSTAGE(N_ROW)
      DO 720 I = ICWS(NCW),ICWE(NCW)-1
C 
      IF(KC(NCW).EQ.1) THEN
           PAVG_IN  = 0.25*(P(I,J,1)+P(I,J-1,1)+P(I+1,J,1)+P(I+1,J-1,1))
           RAVG     = 0.5*(R(J,1)+R(J-1,1))
      END IF
      IF(KC(NCW).EQ.KM) THEN
           PAVG_IN  = 0.25*(P(I,J,KM)+P(I,J-1,KM)+P(I+1,J,KM)
     &                    +P(I+1,J-1,KM))
           RAVG     = 0.5*(R(J,KM) + R(J-1,KM))
      END IF

      VWALL       = WRAD_COOLW(NCW)*RAVG
      TOREL       = TOIN  + (0.5*VWALL*VWALL - W_PRESWIRL)/CPGAS
      TOABS       = TOREL + 0.5*VWALL*VWALL/CPGAS
      PO_REL      = POIN*(TOREL/TOIN)**(1.0/FGAGAS)
      PO_EJECTW(NCW) = PO_REL
      PCOOL_RAT   = PAVG_IN/PO_REL
      IF(PCOOL_RAT.LT.0.5)    PCOOL_RAT = 0.5
      IF(PCOOL_RAT.GT.0.9999) PCOOL_RAT = 0.9999
      PAVG_IN     = PCOOL_RAT*PO_REL
      TCOOL_RAT   = PCOOL_RAT**FGAGAS
      TIN         = TOREL*TCOOL_RAT
      ROIN        = PAVG_IN/RGAS/TIN           
      VIN         = SQRT(2.*CPGAS*(TOREL - TIN))
      HO_ABS(J,I) = CPGAS*TOABS
C
      TORELW(NCW)       = TOREL
      TSTAT_EJECTW(NCW) = TIN
C   EVALUATE THE NORMAL AND TANGENTIAL COOLANT RELATIVE VELOCITIES
      VNIN(J,I)     = VIN*SSANGLEW
      VSIN(J,I)     = VIN*CSANGLEW
C
C      EVALUATE THE  MASS FLOW THROUGH THE PATCH IF IT FILLS THE WHOLE CELL SURFACE AREA.
C      THIS IS NOT THE FINAL COOLANT FLOW.
      IF(KC(NCW).EQ.1)  CELL_FLOW(J,I) = VNIN(J,I)*ROIN*ATOTK1(I,J)
      IF(KC(NCW).EQ.KM) CELL_FLOW(J,I) = VNIN(J,I)*ROIN*ATOTKM(I,J)
C   CALCULATE THE TOTAL FLOW THROUGH THE PATCH IF THE FULL AREA IS USED.
      PATCH_FLOW = PATCH_FLOW + CELL_FLOW(J,I)*NBLADE(J) 
C
C   END OF  I, J   LOOPS
  720 CONTINUE
C
C    EVALUATE THE FRACTION OF CELL AREA USED BY THE COOLING FLOW IN ORDER TO OBTAIN
C    THE SPECIFIED COOLANT FLOWRATE.
       FRAC_FLOW = CFLOWW(NCW)/PATCH_FLOW
C******************************************************************************
C******************************************************************************
C
C   NOW RECALCULATE THE COOLANT FLOWS SO THAT THEY SUM TO THE SPECIFIED TOTAL FLOW
      DO 750 J = JCWS(NCW)+1,JCWE(NCW)
      DO 750 I = ICWS(NCW),ICWE(NCW)-1
C
C   SET THE COOLANT FLOW THROUGH THE CELL FACES SO THAT THEY SUM TO THE SPECIFIED FLOW.
      CFLOWIN        = CELL_FLOW(J,I) * FRAC_FLOW
C
      IF(KC(NCW).EQ.1)  RAVG = 0.5*(R(J,1)  + R(J-1,1))
      IF(KC(NCW).EQ.KM) RAVG = 0.5*(R(J,KM) + R(J-1,KM))
      VWALL          = WRAD_COOLW(NCW)*RAVG
      WPUMP(N_STAGE) = WPUMP(N_STAGE) + NBLADE(J)*
     &                 CFLOWIN*(VWALL*VWALL - W_PRESWIRL)
C
      IF(KC(NCW).EQ.1) THEN
           ANORM  = ATOTK1(I,J)
           XNORM  = ASX(I,J,1)/ANORM
           RNORM  = ASR(I,J,1)/ANORM
      ELSE
           ANORM  = ATOTKM(I,J)
           XNORM  = ASX(I,J,KM)/ANORM
           RNORM  = ASR(I,J,KM)/ANORM
      END IF
C
      ANG    =  TANGLEW(NCW)
      XRNORM =  SQRT(XNORM*XNORM + RNORM*RNORM)
      SX     =  COS(ANG)*RNORM/XRNORM
      ST     =  SIN(ANG)
      SR     = -COS(ANG)*XNORM/XRNORM
C
      IF(KC(NCW).EQ.1)  THEN
           XVEL          =  VNIN(J,I)*XNORM + VSIN(J,I)*SX
           RVEL          =  VNIN(J,I)*RNORM + VSIN(J,I)*SR
           TVEL          =  VSIN(J,I)*ST
           CFLOWK1(I,J)  = CFLOWIN
           HOCWLK1(I,J)  = CFLOWIN*(HO_ABS(J,I) + VWALL*TVEL)
           VXCWLK1(I,J)  = CFLOWIN*XVEL
           VRCWLK1(I,J)  = CFLOWIN*RVEL
           RVTCWLK1(I,J) = CFLOWIN*(VWALL  +  TVEL)*RAVG
      ELSE
           XVEL          =  -VNIN(J,I)*XNORM + VSIN(J,I)*SX
           RVEL          =  -VNIN(J,I)*RNORM + VSIN(J,I)*SR
           TVEL          =   VSIN(J,I)*ST
           CFLOWKM(I,J)  = CFLOWIN
           HOCWLKM(I,J)  = CFLOWIN*(HO_ABS(J,I) + VWALL*TVEL)
           VXCWLKM(I,J)  = CFLOWIN*XVEL
           VRCWLKM(I,J)  = CFLOWIN*RVEL
           RVTCWLKM(I,J) = CFLOWIN*(VWALL  +  TVEL)*RAVG
      END IF
C   END OF THE  I, J  LOOP
  750 CONTINUE
C
C******************************************************************************
C     END OF SETTING THE BLADE SURFACE COOLING FLUXES.
  700 CONTINUE
C
C******************************************************************************
C******************************************************************************
C******************************************************************************
C******************************************************************************
C
C     FIND THE SUM OF THE COOLING FLOWS ADDED UP TO EACH J STATION
C
      SUMCWL(1) = 0.0
      DO 20 J=2,JM
      COOLADD   = 0.0
      DO 25 K=1,KMM1
   25 COOLADD   = COOLADD  + (CFLOWI1(J,K) + CFLOWIM(J,K))*NBLADE(J)
      DO 30 I=1,IMM1
   30 COOLADD   = COOLADD  + (CFLOWK1(I,J) + CFLOWKM(I,J))*NBLADE(J)
      SUMCWL(J) = SUMCWL(J-1) + COOLADD
   20 CONTINUE
C
C     END OF COOLING FLOW CALCULATION.
C******************************************************************************
C******************************************************************************
C
      RETURN
      END
C
C**********************************************************************
