MODULE Orbits_Module
!
!---------------------------------------------------------------------------------
!
! Written by Stan Kidder, CIRA/CSU, 27 May 2020
!
! Subroutines:
!   SUBROUTINE DECODEL(line1, line2)
!   SUBROUTINE POSIT(jyr, jmo, jdy, time, lat, lon, r, x, y, z)
!   SUBROUTINE KEPLER(ma, epsilon, theta)
!   SUBROUTINE ROT(x, y, z, angle, axis)
!   SUBROUTINE DOT_CALC(A, E, INC, MADOT, RADOT, APDOT)
!
! Functions:
!   DOUBLE PRECISION FUNCTION TIMEDIF(iyr, imo, idy, time1, jyr, jmo, jdy, time2)
!   DOUBLE PRECISION FUNCTION GREENWICH_RA(iyr, imo, idy, time)
! 
!
!---------------------------------------------------------------------------------
!

   USE Time_Module

   IMPLICIT NONE
   
   PRIVATE
   PUBLIC :: DECODEL, TIMEDIF, GREENWICH_RA, POSIT, KEPLER, ROT

! Global Parameters
   DOUBLE PRECISION, PARAMETER, PUBLIC :: GMe = 3.986004418D14 ! WGS 84 Orbital constant (m**3 s**-2)
   DOUBLE PRECISION, PARAMETER, PUBLIC :: Ree = 6.378137D3     ! WGS 84 Equatorial radius of Earth (km)
!   DOUBLE PRECISION, PARAMETER, PUBLIC :: OMEGA = 7.292115D-5  ! WGS 84 Rotation rate of Earth (radians/s)
   DOUBLE PRECISION, PARAMETER, PUBLIC :: J2 = 1.08263D-3      ! Quadrupole coefficient of Earth (unitless)
   DOUBLE PRECISION, PARAMETER, PUBLIC :: FLAT = 1.D0/298.257222101D0 ! WGS 84 flattening of Earth (unitless)
   DOUBLE PRECISION, PARAMETER, PUBLIC :: PI = 4.D0*ATAN(1.D0) 
   DOUBLE PRECISION, PARAMETER, PUBLIC :: TWOPI = 2.D0*PI
   DOUBLE PRECISION, PARAMETER, PUBLIC :: RCON = 180.D0/PI
   DOUBLE PRECISION, PARAMETER, PUBLIC :: SEC_PER_DAY = 86400.D0
      
! Global variables
   CHARACTER(LEN=50), PUBLIC :: sat_id        ! satellite number (YYNNNAA)
   DOUBLE PRECISION, PUBLIC  :: a             ! semimajor axis (km)
   DOUBLE PRECISION, PUBLIC  :: inc           ! inclination angle (degrees)
   DOUBLE PRECISION, PUBLIC  :: eccen         ! eccentricity (unitless)
   DOUBLE PRECISION, PUBLIC  :: ma0           ! mean anomaly at epoch (degrees)
   DOUBLE PRECISION, PUBLIC  :: ra0           ! right ascension of ascending node at epoch (degrees)
   DOUBLE PRECISION, PUBLIC  :: ap0           ! argument of perigee at epoch (degrees)
   DOUBLE PRECISION, PUBLIC  :: madot         ! rate of change of mean anomaly (degrees/day)
   DOUBLE PRECISION, PUBLIC  :: radot         ! rate of change of right ascension (degrees/day)
   DOUBLE PRECISION, PUBLIC  :: apdot         ! rate of change of argument of perigee (degrees/day)
   DOUBLE PRECISION, PUBLIC  :: ttilda        ! nodal (or synodic) period 
   INTEGER, PUBLIC           :: eyr           ! epoch time year
   INTEGER, PUBLIC           :: emo           ! epoch time month
   INTEGER, PUBLIC           :: edy           ! epoch time day-of-month
   DOUBLE PRECISION, PUBLIC  :: etime         ! epoch time (fraction of a day)
   INTEGER, PUBLIC           :: numrev        ! orbit number
   DOUBLE PRECISION, PUBLIC  :: revs          ! revolutions per day (

CONTAINS

   SUBROUTINE DECODEL(line1, line2)
!      
! Decodes two-line elements
!
! Variables

      IMPLICIT NONE  
      
      ! Arguments      
      CHARACTER(LEN = 80), INTENT(IN) :: line1, line2
      
      ! Local variables
      INTEGER :: k
      INTEGER :: numsat
      INTEGER :: doy     ! day of year
      DOUBLE PRECISION :: dt, dt2, drag
      CHARACTER(LEN=1) :: etype
      INTEGER :: nelem
      DOUBLE PRECISION :: i
      DOUBLE PRECISION :: fac
      
! Executable Statements

      ! Read the two-line elements
      READ(line1, 1) numsat, sat_id, eyr, doy, etime, dt, dt2, drag, etype, nelem
    1 FORMAT(2X,I5,2X,A8,1X,I2,I3,F9.8,1X,F10.0,1X,F8.0,1X,F8.0,1X,A1,1X,I4)
      READ(LINE2, 2) numsat, inc, ra0, eccen, ap0, ma0, revs, numrev
    2 FORMAT(2X,I5,1X,F8.0,1X,F8.0,1X,F7.7,1X,F8.0,1X,F8.0,1X,F11.0,I5)
!
! Find semimajor axis--start with Keplerian approximation
!
      i = inc/RCON
      IF(eyr < 57) THEN
         eyr = eyr + 2000
      ELSE
         eyr = eyr + 1900
      END IF
      CALL MODAY(doy, eyr, emo, edy)
      
      ! Find semimajor axis--start with Keplerian approximation
      a = (GMe*(SEC_PER_DAY/revs/TWOPI)**2)**(1./3.)
      a = a/1000.
      ! Increment a. Five interations is arbitrary, but enough.
      DO k = 1, 5
         CALL DOT_CALC(a, eccen, inc, madot, radot, apdot)
         a = a*(madot/360./revs)**(2./3.)
      END DO

      ! Calculate final values
      CALL DOT_CALC(a, eccen, inc, madot, radot, apdot)
      ttilda = 360.D0/(madot + apdot)

!        tbar = 360.D0/madot
!        relrot = omega*60.D0*rcon - radot/1440.D0
!        dlon = relrot*ttilda
!        orbday = 1440.D0/ttilda
!        ddlon = 360.D0 - dlon*ANINT(orbday)
!        retper = abs(dlon/ddlon)
        
   END SUBROUTINE DECODEL
   
   SUBROUTINE DOT_CALC(a, eccen, inc, madot, radot, apdot)
!
! Calculates time rates of change of the mean anomaly,
! right ascension, and arg. of perigee

! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(IN)  :: a     ! semimajor axis (km)
      DOUBLE PRECISION, INTENT(IN)  :: eccen ! orbit eccentriciy (unitless)
      DOUBLE PRECISION, INTENT(IN)  :: inc   ! inclination angle (degrees)
      DOUBLE PRECISION, INTENT(OUT) :: madot ! time rate of change of mean anomaly (degrees/day)
      DOUBLE PRECISION, INTENT(OUT) :: radot ! time rate of change of right ascension of ascending node (degrees/day)
      DOUBLE PRECISION, INTENT(OUT) :: apdot ! rate of change of argument of perigee (degrees/day)
      
      ! Local variables
      DOUBLE PRECISION :: mmc   ! mean motion constant (n, radians/second)
      DOUBLE PRECISION :: ammc  ! anomalistic mean motion constant (n bar, radians/second)
      DOUBLE PRECISION :: fac   ! useful factor (unitless)
      DOUBLE PRECISION :: sigma ! anoher useful factor (unitless)
      
! Executable statements      

      mmc = SQRT(GMe/(a*1000.D0)**3)
      fac = 1.D0 - eccen**2
      sigma = 1.5D0*J2*((Ree/a)**2)/fac**2
      ammc = mmc*(1.D0 + sigma*SQRT(fac)*(1.D0 - 1.5D0*SIN(inc/RCON)**2))
      madot =  ammc*RCON*SEC_PER_DAY
      radot = -ammc*sigma*COS(inc/RCON)*RCON*SEC_PER_DAY
      apdot =  ammc*sigma*(2.D0 - 2.5D0*SIN(inc/RCON)**2)*RCON*SEC_PER_DAY

   END SUBROUTINE DOT_CALC
!
   DOUBLE PRECISION FUNCTION TIMEDIF(yr1, mo1, dy1, t1, yr2, mo2, dy2, t2)
!
!       Calculates the time difference (time2 - time1) in days
!
! Variables

      IMPLICIT NONE
      
      ! Arguments
      INTEGER, INTENT(IN)  :: yr1, mo1, dy1  ! year/month/day of first time
      DOUBLE PRECISION, INTENT(IN) :: t1     ! first time (fraction of a day)
      INTEGER, INTENT(IN)  :: yr2, mo2, dy2  ! year/month/day of second time 
      DOUBLE PRECISION, INTENT(IN) :: t2     ! second time (fraction of a day)
      
      ! Local variables
      INTEGER :: jd1, jd2
      
! Executable statements      
      jd1 = Julian_Date(yr1, mo1, dy1)
      jd2 = Julian_Date(yr2, mo2, dy2)
      TIMEDIF = DBLE(jd2) + t2 - DBLE(jd1) - t1

   END FUNCTION TIMEDIF

   DOUBLE PRECISION FUNCTION GREENWICH_RA(iyr, imo, idy, time)
!
! Calculates the right ascension of Greenwich (degrees)
!
! Variables

      IMPLICIT NONE
      
      ! Arguments
      INTEGER, INTENT(IN)  :: iyr, imo, idy
      DOUBLE PRECISION, INTENT(IN) :: time
      
      ! Local parameters
      INTEGER, PARAMETER :: GYR   = 2020
!      INTEGER, PARAMETER :: GYR   = 1990
      INTEGER, PARAMETER :: GMO   = 1
      INTEGER, PARAMETER :: GDY   = 1
      DOUBLE PRECISION, PARAMETER :: GTIME = 0.D0
      DOUBLE PRECISION, PARAMETER :: GRA0 = 100.16113D0
!      DOUBLE PRECISION, PARAMETER :: GRA0 = 100.38641D0  ! for 1/1/1990
      DOUBLE PRECISION, PARAMETER :: GRADOT = 360.9856507D0
      
      ! Local Variables
      DOUBLE PRECISION :: dt
      
! Executable statements 
     
        dt = TIMEDIF(GYR, GMO, GDY, GTIME, iyr, imo, idy, time)
        GREENWICH_RA = MOD(GRA0 + dt*GRADOT, 360.D0)

   END FUNCTION GREENWICH_RA
!
   SUBROUTINE POSIT(yr, mo, dy, time, lat, lon, r, x, y, z)
   
! Calculates the position (latitude, longitude, r) and (x, y, z) of a satellite  
! on date (year, month, day) and time (fraction of a day).
   
! Variables

      IMPLICIT NONE
      
      ! Arguments
      INTEGER, INTENT(IN)  :: yr, mo, dy
      DOUBLE PRECISION, INTENT(IN) :: time
      DOUBLE PRECISION, INTENT(OUT) :: lat, lon, r
      DOUBLE PRECISION, INTENT(OUT) :: x, y, z
            
      ! Local variables
      DOUBLE PRECISION :: dt, ma, ra, ap
      DOUBLE PRECISION :: theta ! true anomaly
      DOUBLE PRECISION :: gra   ! Greenwich right ascension
      
! Executable statements    
  
      dt = TIMEDIF(eyr, emo, edy, etime, yr, mo, dy, time)
      !  Calculate secular orbital elements
      ma = mod(ma0 + madot*dt, 360.d0)
      ra = mod(ra0 + radot*dt, 360.d0)
      ap = mod(ap0 + apdot*dt, 360.d0)

      ! Calculate true anomaly
      CALL KEPLER(ma/RCON, eccen, theta)

      ! Calculate position of satellite in ellipse
      r = a*(1.D0 - eccen**2)/(1.D0 + eccen*COS(theta))
      x = r*COS(theta)
      y = r*SIN(theta)
      z = 0

      ! Rotate ellipse to proper orientation in Right Ascension - Declination
      ! Coordinate System.
      CALL ROT(x, y, z, ap/RCON, 'Z')
      CALL ROT(x, y, z, inc/RCON, 'X')
      CALL ROT(x, y, z, ra/RCON, 'Z')
      
      ! Calculate longitude and geocentric latitude (declination) of satellite
      gra = GREENWICH_RA(yr, mo, dy, time)
      lon = MOD(ATAN2(y, x)*RCON - gra, 360.D0)
      if(lon >  180.D0) lon = lon - 360.D0
      if(lon < -180.D0) lon = lon + 360.D0
      lat = ASIN(z/r)*RCON
      
   END SUBROUTINE POSIT

   DOUBLE PRECISION FUNCTION GEODETIC_LATITUDE(declination)
!
! Calculates the geodetic latitude given the geocentric latitude (aka declination)
!  
! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(IN) :: declination
      
! Executable statements 
     
      Geodetic_Latitude = ATAN(TAN(declination/RCON)/(1 - FLAT)**2)*RCON

   END FUNCTION GEODETIC_LATITUDE
   
   SUBROUTINE KEPLER(ma, epsilon, theta)
!
!     Returns the true anomaly (theta) given the mean anomaly (ma) and the
!     eccentricity (epsilon).  All angles are in radians.
!
! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(IN)  :: ma        ! mean anomaly (radians)
      DOUBLE PRECISION, INTENT(IN)  :: epsilon   ! eccentricity (unitless)
      DOUBLE PRECISION, INTENT(OUT) :: theta     ! true anomaly (radians)
      
      ! Local variables
      DOUBLE PRECISION :: mm
      DOUBLE PRECISION :: mguess
      DOUBLE PRECISION :: e       ! eccentric anomaly
      INTEGER :: iter
      DOUBLE PRECISION :: de
      
! Executable statements

      ! Iterate e using Newton-Raphson approach.  Stop when e increment is less
      ! than 1.E-6*e or iter > 10.

      mm = MOD(ma, TWOPI) ! mean anomaly modulo two pi
      e = mm              ! first guess for e
      iter = 0
      DO
         mguess = e - epsilon*SIN(e)
         de = (mguess - mm)/(1.D0 - epsilon*COS(e))
         e = e - de
         iter = iter + 1
         IF(ABS(de) < 1.D-6 .OR. iter > 10) EXIT
      END DO
!
!     Calculate true anomaly (THETA) from eccentric anomaly (E).
!
      theta = ACOS((COS(e) - epsilon)/(1.D0 - epsilon*COS(e)))
      IF(mm > PI) theta = TWOPI - theta
      
   END SUBROUTINE KEPLER

   SUBROUTINE ROT(x, y, z, angle, axis)
 
! Rotates (x, y, z) by "angle" about "axis".
! angle is in radians  
 
! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(INOUT) :: x, y, z
      DOUBLE PRECISION, INTENT(IN)    :: angle      ! radians
      CHARACTER(LEN=1), INTENT(IN)    :: axis 
      
      ! Local variables 
      INTEGER :: i, j 
      DOUBLE PRECISION :: cosa, sina
      
      ! Local arrays
      DOUBLE PRECISION :: vector(3)
      DOUBLE PRECISION :: R(3, 3)
      
! Executable statements
      
      ! Set R equal to the identity matrix.
      R = 0.
      R(1, 1) = 1.D0
      R(2, 2) = 1.D0
      R(3, 3) = 1.D0

      ! Set up rotation matrix R depending on AXIS.
      cosa = COS(angle)
      sina = SIN(angle)
      IF(AXIS == 'X' .OR. AXIS =='x') THEN
         i = 2
         j = 3
      ELSE IF(AXIS == 'Y' .OR. AXIS == 'y') THEN
         i = 1
         j = 3
      ELSE ! Default to rotation about the z axis.
         i = 1
         j = 2
      END IF
      R(i, i) =  cosa
      R(j, j) =  cosa
      R(i, j) = -sina
      R(j, i) =  sina

      ! Perform rotation
      vector(1) = x
      vector(2) = y
      vector(3) = z
      vector = MATMUL(R, vector)
      x = vector(1)
      y = vector(2)
      z = vector(3)
      
   END SUBROUTINE ROT

END MODULE Orbits_Module