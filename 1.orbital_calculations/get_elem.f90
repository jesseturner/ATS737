PROGRAM GET_ELEM
! 
! --------------------------------------------------------------------------------------------------
! Gets satellite orbital elements from various sources and
! creates an output file (ORB_ELEM.DAT) with the elements in standard format
! for use by succeeding programs
! 
! Translated from Fortran 77 on 11 May 2020 by SQK
! Annotated with equation numbers from Kidder & Vonder Haar 11 Nov 2022
! Note: K&VH = Kidder, S. Q., and T. H. Vonder Haar (1995). Satellite Meteorology: An Introduction. 
!              Academic Press, San Diego, 466 pp.
! 
! The output elements:
!    
!   sat_name    =  name of the satellite
!   edate       =  utc epoch date (yyyymmdd)
!   etime       =  utc epoch time (fraction of a day)
!   a           =  semimajor axis (km)
!   inc         =  inclination angle (deg)
!   e           =  eccentricity (unitless)
!   ra          =  right ascension of ascending node (deg)
!   radot       =  rate of change of ra (deg/day)
!   ap          =  argument of perigee (deg)
!   apdot       =  rate of change of the arg. of perigee (deg/day)
!   ma          =  mean anomaly at epoch time (deg)
!   madot       =  rate of change of the mean anomaly (deg/day)
!   gra         =  Greenwich right ascension at epoch (deg)
!   gradot      =  rate of change of gra (deg/day)
! --------------------------------------------------------------------------------------------------

! Modules

   USE Time_Module
   
! Variables

   IMPLICIT NONE
   
   ! Parameters
   INTEGER, PARAMETER :: GYR  =  1990
   INTEGER, PARAMETER :: GMO  =  1
   INTEGER, PARAMETER :: GDY  =  1
!   DOUBLE PRECISION, PARAMETER :: gradot  =  360.9856507D0
   DOUBLE PRECISION, PARAMETER :: GMe  =  3.986005D14
   DOUBLE PRECISION, PARAMETER :: GRA0  =  100.38641D0
   DOUBLE PRECISION, PARAMETER :: TWOPI  =  8.D0*ATAN(1.D0)
   INTEGER, PARAMETER :: NS  =  19
   CHARACTER(LEN = 7), PARAMETER :: SAT(NS)    =  (/ &
                                  'TIROS N',   'NOAA 6 ',   'NOAA 7 ',   'NOAA 8 ',   &
                                  'NOAA 9 ',   'NOAA 10',   'NOAA 11',   'NOAA 12',   &
                                  'NOAA 13',   'NOAA 14',   'NOAA 15',   'NOAA 16',   &
                                  'NOAA 17',   'NOAA 18',   'NOAA 19',   'SNPP   ',   &
                                  'NOAA 20',   'NOAA 21',   'UNKNOWN'/)
   CHARACTER(LEN = 9), PARAMETER :: SATID(NS)  =  (/ &
                                  '1978 096A', '1979 057A', '1981 085A', '1983 022A', &
                                  '1984 123A', '1986 073A', '1988 089A', '1991 032A', &
                                  '1993 050A', '1994 089A', '1998 030A', '2000 055A', &
                                  '2002 032A', '2005 018A', '2009 005A', '2011 061A', &
                                  '2017 073A', '2022 150A', '         ' /)
   
   ! Namelist
   CHARACTER(LEN = 50) :: sat_name
   INTEGER :: edate
   DOUBLE PRECISION :: etime, a, inc, e, ra, radot, ap, apdot, ma, madot, gra
   NAMELIST /ORB_ELEM/ sat_name, edate, etime, a, inc, e, ra, radot, &
                       ap, apdot, ma, madot, gra, gradot
                       
   ! Local Variables
   DOUBLE PRECISION :: gradot = 360.9856507D0
   CHARACTER(LEN = 255) :: fn
   CHARACTER(LEN = 80)  :: line1, line2, line3, line4, line5
   CHARACTER(LEN = 1)   :: ca, cr, cm
   INTEGER :: j, k, kk
   INTEGER :: iyr, imo, idy
   INTEGER :: jdy
   INTEGER :: ihr, imin
   DOUBLE PRECISION :: sec
   INTEGER :: Jg, Je
   DOUBLE PRECISION :: dt
   INTEGER :: nsats
   INTEGER :: isat, jsat
   DOUBLE PRECISION :: revs
   DOUBLE PRECISION :: ehr, emin, esec
   INTEGER :: ios
!
! Select input method
!
   10 PRINT 1
    1 FORMAT(' Select the input method:'          &
             '   1. Manual input'                 &
             '   2. 2 - line elements from file'  &
             '   3. TBUS Bulletin from file'      &
             '                               >')
      READ *, k
      IF(k == 1) GOTO 100
      IF(k == 2) GOTO 200
      IF(k == 3) GOTO 300
      GOTO 10
!
! Manual input
!
  100 PRINT 102
  102 FORMAT(' Satellite name >')
      READ (*,'(A50)') sat_name
      PRINT 103
  103 FORMAT(' UTC Epoch date (YYYY MM DD) >')
      READ *, iyr, imo, idy
      edate  =  iyr*10000 + imo*100 + idy
      PRINT 104
  104 FORMAT(' UTC Epoch time (HH MM SS.SSS) >')
      READ *, ihr, imin, sec
      ETIME = IHR/24. + IMIN/1440. + SEC/86400.
      PRINT 105
  105 FORMAT(' Semimajor axis (km) >')
      READ *, a
      PRINT 155
  155 FORMAT(' Eccentricity >')
      READ *, e
      PRINT 106
  106 FORMAT(' Inclination angle (deg) >')
      READ *, inc
      PRINT 107
  107 FORMAT(' Right ascension of ascending node (deg)>')
      READ *, ra
      ra = MOD(ra, 360.D0)
      IF(RA < 0.) ra = ra + 360.D0
      PRINT 108
  108 FORMAT(' Argument of perigee (deg) >')
      READ *, ap
      ap = MOD(ap, 360.D0)
      IF(AP < 0.) ap = ap + 360.D0
      PRINT 109
  109 FORMAT(' Mean anomaly (deg) >')
      READ *, ma
      ma = MOD(ma, 360.D0)
      IF(ma < 0.) ma = ma + 360.D0

      ! Calculate rates of change
      CALL DOT_CALC(a, e, inc, madot, radot, apdot)

      ! Calculate Greenwich right ascension
      Jg = Julian_Date(GYR, GMO, GDY)
      Je = Julian_Date(iyr, imo, idy)
!      CALL JDAY(JG,GYR,GMO,GDY,IW)
!      CALL JDAY(JE,IYR,IMO,IDY,IW)
      dt = Je + etime - Jg
      gra = MOD(GRA0 + GRADOT*dt, 360.D0)
      GOTO 1000
!
! 2-line elements (actualy 3-line elements)
!
  200 PRINT 201 
  201 FORMAT(' 2-Line element file name >')
      READ *, fn
      OPEN(UNIT = 1, FILE = fn, STATUS = 'OLD', ERR = 200)

      ! Read file and present list of satellites
      nsats = 0
  210 READ(1, 202, END = 230) line1
      READ(1, 202, END = 230) line2
  220 READ(1, 202, END = 230) line3
  202 FORMAT(A80)
      IF(line2(1:2) == '1 '.AND. line3(1:2) == '2 ') THEN
         nsats = nsats + 1
         PRINT 203, FLOAT(nsats), line1(1:22)
  203    FORMAT(' ',F4.0,2X,A22)
         GOTO 210
      ELSE
         line1 = line2
         line2 = line3
         GOTO 220
      END IF
  230 IF(nsats < 1) THEN
         WRITE(*,*) 'NO VALID SATELLITES FOUND'
         STOP
      END IF
  240 WRITE(*,204)
  204 FORMAT(' Enter number of desired satellite (0 to exit) >')
      READ *, isat
      IF(isat == 0)  STOP
      IF(isat < 1 .OR. isat > nsats) GOTO 240
      
      ! Read the data
      REWIND 1
      jsat = 0
  250 READ(1,202) line1
      READ(1,202) line2
  260 READ(1,202) line3
      IF(line2(1:2) == '1 '.AND. line3(1:2) == '2 ') THEN
         jsat = jsat + 1
         IF(jsat < isat) GOTO 250
      ELSE
         line1 = line2
         line2 = line3
         GOTO 260
      END IF
      sat_name = line1(1:22)
      READ(line2, 205) iyr, jdy, etime
  205 FORMAT(18X,I2,I3,F9.8)
      READ(line3, 206) inc,ra,e,ap,ma,revs
  206 FORMAT(8X,F8.0,1X,F8.0,1X,F7.7,1X,F8.0,1X,F8.0,1X,F11.0)
      iyr = iyr + 1900
      IF(iyr < 1957) iyr = iyr + 100
      CALL MODAY(jdy, iyr, imo, idy)
      edate = iyr*10000 + imo*100 + idy
      
      ! Find semimajor axis--start with Keplerian approximation
      a = (GMe*(86400.D0/revs/TWOPI)**2)**(1./3.)  ! K&VH Eq. 2.9
      a = a/1000.
!      PRINT *, "Incrementing a:"
      k = 0
!      PRINT *, k, a
      DO k = 1, 4
         CALL DOT_CALC(a, e, inc, madot, radot, apdot)
         a = a*(madot/360./revs)**(2./3.)  ! K&VH Eq. 2.9
!         PRINT *, k, a
      END DO

      ! Calculate final values
      CALL DOT_CALC(a, e, inc, madot, radot, apdot)

      ! Calculate Greenwich right ascension
      Jg = Julian_Date(GYR, GMO, GDY)
      Je = Julian_Date(iyr, imo, idy)
      dt = Je + etime - Jg
      gra = MOD(GRA0 + GRADOT*dt, 360.D0)
!      PRINT *, 'A  = ',A
      CLOSE(1)
      GOTO 1000
!
! TBUS Bulletin input
!
      ! Open TBUS file
  300 WRITE(*, 301) 
  301 FORMAT(' Name of the file which contains the TBUS bulletin >'$)
      READ(*,'(A50)') fn
      OPEN(UNIT = 1, FILE = fn, STATUS = 'OLD', ERR = 300)

      ! Find the sections labeled "PART IV"
      nsats = 0
  310 READ(1, '(A80)', END = 320) line1
      k = INDEX(line1, "PART IV")
      IF(k == 0) GOTO 310
      nsats = nsats + 1
      READ(1,'(A80)') line1
      j = NS
      DO kk = 1, NS - 1
         k = INDEX(line1, SATID(kk))
         IF(k > 0) j = kk
      END DO
      PRINT 302, FLOAT(nsats), SAT(j)
  302 FORMAT(' ',F4.0,2X,A7)
      GOTO 310
  320 IF(nsats < 1) THEN
         WRITE(*,*) 'NO VALID SATELLITES FOUND'
         STOP
      END IF
  330 WRITE(*, 303)
  303 FORMAT(' Enter number of desired satellite (0 to exit) >')
      READ *, isat
      IF(isat == 0)  STOP
      IF(isat < 0 .OR. isat > nsats) GOTO 330
      REWIND 1
      jsat = 0
  340 READ(1, '(A80)') line1
      k = INDEX(line1, "PART IV")
      IF(k == 0) GOTO 340
      jsat = jsat + 1
      IF(jsat < isat) GOTO 340
      READ(1,'(A80)') line1
      READ(1,'(A80)') line2
      READ(1,'(A80)') line3
      READ(1,'(A80)') line4
      READ(1,'(A80)') line5
      j = NS
      DO kk = 1, NS - 1
         k = INDEX(line1, SATID(kk))
         IF(k > 0) j = kk
      END DO
      sat_name = SAT(j)
      line1 = ADJUSTL(line1)
      READ(line1, 304) iyr, imo, idy, ehr, emin, esec, gra
  304 FORMAT(29X,3I2,2F2.0,F5.3,F8.4)
      iyr = iyr + 1900
      if(iyr < 1957) iyr = iyr + 100
      edate = iyr*10000 + imo*100 + idy
      etime = ehr/24. + emin/1440. + esec/86400.
      line2 = ADJUSTL(line2)
      READ(line2, 305) e, ap, ra, inc
  305 FORMAT(18X,F8.8,3(1X,F8.5))
      line3 = ADJUSTL(line3)
      READ(line3, 306) ma, a
  306 FORMAT(F8.5,1X,F8.3)
      line5 = ADJUSTL(line5)
      READ(line5,307) ca, apdot, cr, radot, cm, madot
  307 FORMAT(11X,A1,F8.5,1X,A1,F8.5,1X,A1,F8.2)
      IF(ca == 'M') apdot =  -apdot
      IF(cr == 'M') radot =  -radot
      IF(cm == 'M') madot =  -madot
!
! Write elements to file and exit
!
1000  OPEN(UNIT = 2, FILE = 'ORB_ELEM.DAT', DELIM = 'QUOTE', STATUS = 'REPLACE')
      CALL NO_SPACES(SAT_NAME)
      WRITE(2, NML=ORB_ELEM, IOSTAT = ios)
      CLOSE(2)
      IF(ios == 0) THEN
         PRINT *, 'Elements written to file ORB_ELEM.DAT'
      ELSE
         PRINT *, "Error writing to ORB_ELEM.DAT"
         PRINT *, "IOSTAT = ", ios
      END IF
      WRITE(*, NML=ORB_ELEM, IOSTAT = ios)
      IF(ios /= 0) PRINT *, "IO Error: IOSTAT = ", ios
      
CONTAINS      

   SUBROUTINE DOT_CALC(A, E, INC, MADOT, RADOT, APDOT)
!
! CALCULATES TIME RATES OF CHANGE OF THE MEAN ANOMALY,
! RIGHT ASCENSION, AND ARG. OF PERIGEE
!
      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(IN)  :: a, e, inc
      DOUBLE PRECISION, INTENT(OUT) :: madot, radot, apdot
      
      ! Parameters
      DOUBLE PRECISION, PARAMETER :: GMe = 3.986005D14
      DOUBLE PRECISION, PARAMETER :: J2 = 1.08263D-3
      DOUBLE PRECISION, PARAMETER :: Ree = 6.378137D6
      DOUBLE PRECISION, PARAMETER :: RCON = 45.D0/ATAN(1.D0)
      
      ! Local variables
      DOUBLE PRECISION :: mmc, fac, sigma, ammc
      
! Executable statements      

      mmc = SQRT(GMe/(a*1000.D0)**3)  ! Mean motion constant, K&VH Eq. 2.9
      fac = 1.D0 - e**2
      sigma = 1.5D0*J2*(Ree/(a*1000.D0)/fac)**2
      ammc = (1.D0 + sigma*SQRT(fac)*(1.D0 - 1.5D0*SIN(inc/RCON)**2))*mmc ! Anomalistic MMC, K&VH Eq. 2.12
      madot = ammc*rcon*86400.D0                                          ! Converting AMMC to deg/day from rad/sec.
      radot = -sigma*COS(inc/RCON)*ammc*RCON*86400.D0                     ! K&VH Eq. 2.13, deg/day
      apdot = sigma*(2.D0 - 2.5D0*SIN(inc/RCON)**2)*ammc*RCON*86400.D0    ! K&VH Eq. 2.14, deg/day

   END SUBROUTINE DOT_CALC

   SUBROUTINE NO_SPACES(C)
   
      IMPLICIT NONE
      
! Argument
      CHARACTER(LEN=*), INTENT(INOUT) :: C
      
! Local variables
      INTEGER :: i, n
      
! Executable statements
      
      n = LEN_TRIM(C)
      DO i = 1, n
         IF(C(i:i) == ' ') C(i:i) = '_'
      END DO

   END SUBROUTINE NO_SPACES

END PROGRAM GET_ELEM      
