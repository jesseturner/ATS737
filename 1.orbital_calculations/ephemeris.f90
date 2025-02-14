PROGRAM ephemeris

! --------------------------------------------------------------------------------------------------
! This program uses orbital elements from file ORB_ELEM.DAT
! to create an ephemeris file named EPHEMERIS.DAT
! 
! EPHEMERIS.FOR written 6/24/1997 by SQK
! Translated to Fortran 90 (ephemeris.f90) 5/13/2020 by SQK
! Annotated with equation numbers from Kidder & Vonder Haar 11 Nov 2022
! Note: K&VH = Kidder, S. Q., and T. H. Vonder Haar (1995). Satellite Meteorology: An Introduction. 
!              Academic Press, San Diego, 466 pp.

! The orbital elements:
!    
!   sat_name   = name of the satellite
!   edate      = utc epoch date (yyyymmdd)
!   etime      = utc epoch time (fraction of a day)
!   a          = semimajor axis (km)
!   inc        = inclination angle (deg)
!   e          = eccentricity (unitless)
!   ra         = right ascension of ascending node (deg)
!   radot      = rate of change of ra (deg/day)
!   ap         = argument of perigee (deg)
!   apdot      = rate of change of the arg. of perigee (deg/day)
!   ma         = mean anomaly at epoch time (deg)
!   madot      = rate of change of the mean anomaly (deg/day)
!   gra        = greenwich right ascension at epoch (deg)
!   gradot     = rate of change of gra (deg/day)
! --------------------------------------------------------------------------------------------------

! Modules

      USE Time_Module
   
! Variables

      IMPLICIT NONE
   
      ! Parameters
      DOUBLE PRECISION, PARAMETER :: PI    =  4.D0*ATAN(1.D0)
      DOUBLE PRECISION, PARAMETER :: RCON  =  180.D0/PI
      
      ! Namelist
      CHARACTER(LEN = 50) :: sat_name
      INTEGER :: edate
      DOUBLE PRECISION :: etime, a, inc, e, ra, radot, ap, apdot, ma, madot, gra, gradot
      NAMELIST /ORB_ELEM/ sat_name, edate, etime, a, inc, e, ra, radot, &
                          ap, apdot, ma, madot, gra, gradot

      DOUBLE PRECISION :: satvec(3)
      DOUBLE PRECISION :: x, y, z
      EQUIVALENCE (satvec(1),x),(satvec(2),y),(satvec(3),z)
      DOUBLE PRECISION :: r
      DOUBLE PRECISION :: lat, lon
      DOUBLE PRECISION :: ma0, ra0, ap0, gra0
      INTEGER :: eyr, emo, edy
      INTEGER :: iyr, imo, idy
      INTEGER :: jyr, jmo, jdy
      INTEGER :: j1, j2, je
      INTEGER :: ihr, imin
      INTEGER :: jhr, jmin
      INTEGER :: khr, kmin
      DOUBLE PRECISION :: seci, secj, seck
      DOUBLE PRECISION :: t1, t2
      DOUBLE PRECISION :: delta
      INTEGER :: nt
      DOUBLE PRECISION :: stime, ftime
      INTEGER :: istime, iftime
      INTEGER :: it
      DOUBLE PRECISION :: time, dt
      DOUBLE PRECISION :: theta
!
! Get orbital elements
!
      OPEN(UNIT = 1,  FILE = 'ORB_ELEM.DAT',  STATUS = 'OLD')
      READ(1, NML = ORB_ELEM)
      CLOSE(1)
      WRITE(*, NML = ORB_ELEM)
!
! Get start and end times for ephemeris
!
      PRINT 1
    1 FORMAT(' START DATE (YYYY MM DD)   >')
      READ *,  iyr, imo, idy
      j1 = Julian_Date(iyr, imo, idy)
      PRINT 2
    2 FORMAT(' START TIME (HH MM S.SSS)  >')
      READ *,  ihr, imin, seci
      t1 = ihr/24.D0 + imin/1440.D0 + seci/86400.D0
      PRINT 3
    3 FORMAT(' END DATE (YYYY MM DD)     >')
      READ *,  jyr, jmo, jdy
      j2 = Julian_Date(jyr, jmo, jdy)
      PRINT 4
    4 FORMAT(' END TIME (HH MM S.SSS)    >')
      READ *,  jhr, jmin, secj
      t2 = jhr/24.D0 + jmin/1440.D0 + secj/86400.D0
      PRINT 5
    5 FORMAT(' TIME INCREMENT (HH MM S.SSS) >')
      READ *,  khr, kmin, seck
      eyr = edate/10000
      emo = (edate - eyr*10000)/100
      edy = edate - eyr*10000 - emo*100
      je = Julian_Date(eyr, emo, edy)
      t1 = (j1 - je) + t1
      t2 = (j2 - je) + t2
      delta = khr/24.D0 + kmin/1440.D0 + seck/86400.D0
      nt = NINT((t2 - t1)/delta)
      ma0 = ma
      ra0 = ra
      ap0 = ap
      gra0 = gra
!
! Open ephemeris file
!
      IF(nt < 1) THEN
         PRINT *,  'TIMES INCORRECT -  - NO LINES WRITTEN'
         STOP
      END IF
      OPEN(UNIT = 2, FILE = 'EPHEMERIS.DAT', STATUS = 'REPLACE') 
      WRITE(2, '(A50)') sat_name
      stime = ihr*10000 + imin*100 + seci
      istime = int(stime)
      ftime = jhr*10000 + jmin*100 + secj
      iftime = int(ftime)
      WRITE(2, 6) edate, iyr*10000 + imo*100 + idy, istime, &
                 stime - istime, jyr*10000 + jmo*100 + jdy, &
                 iftime, ftime - iftime, khr, kmin, seck
    6 FORMAT(I8, 2(I9, I7.6, F4.3), 2I4, F10.3)
      WRITE(2, 7) nt + 1
    7 FORMAT(I5, ' LINES')
      WRITE(2, 8)
    8 FORMAT("     TIME (days)      Geocentric Latitude       Longitude      ", &
             "       R (km)               X (km)               Y (km)               Z (km)")
!
! Loop through all times
!
      DO it = 0, nt
!
!        Update orbital elements and Greenwich right ascension
!
         time = t1 + it*delta    ! TIME is measured in days from 0000 UTC on the epoch date
         dt = time - etime
         ma =  MOD(ma0  + madot *dt, 360.D0)
         IF(ma < 0.) ma = ma + 360.D0
         ra =  MOD(ra0  + radot *dt, 360.D0)
         IF(ra < 0.) ra = ra + 360.D0
         ap =  MOD(ap0  + apdot *dt, 360.D0)
         IF(ap < 0.) ap = ap + 360.D0
         gra = MOD(gra0 + gradot*dt, 360.D0)
         IF(gra < 0.) gra = gra + 360.D0
!
!        Calculate true anomaly
!
         CALL Kepler(ma/RCON, e, theta)  
!
!        Calculate position of satellite in ellipse
!
         r = a*(1 - e**2)/(1 + e*cos(theta))  ! K&VH Eq. 2.7
         x = r*COS(theta)
         y = r*SIN(theta)
         z = 0
!
!        Rotate ellipse to proper orientation in Right Ascension - Declination
!        Coordinate System.
!
         CALL ROT(satvec, ap/RCON,  'z')
         CALL ROT(satvec, inc/RCON, 'x')
         CALL ROT(satvec, ra/RCON,  'z')
!
!        Calculate latitude and longitude of satellite (assuming spherical earth)
!
         lat = ASIN(z/r)*RCON  ! K&VH Eq. 2.26b
         lon = MOD(ATAN2(y, x)*RCON - gra, 360.D0)  ! K&VH Eq. 2.26c and Fig. 2.12
         if(lon >  180.D0) lon = lon - 360.D0
         if(lon < -180.D0) lon = lon + 360.D0
!
!        Write results
!
         WRITE(2, 9) time, lat, lon, r, x, y, z
    9    FORMAT(7(1PE21.13))
    
      END DO
      CLOSE(2)
      PRINT *,  'Output written to file EPHEMERIS.DAT'
      
CONTAINS

   SUBROUTINE Kepler(ma, eccen, theta)

! --------------------------------------------------------------------------------------------------
!  Solves Kepler's Equation (K&VH Eq. 2.8)
! 
!  o Returns the true anomaly (theta) given the mean anomaly (ma) and the
!    eccentricity (eccen).  
!  o All angles are in radians.
! --------------------------------------------------------------------------------------------------

! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(IN)  :: ma        ! mean anomaly
      DOUBLE PRECISION, INTENT(IN)  :: eccen     ! eccentricity
      DOUBLE PRECISION, INTENT(OUT) :: theta     ! true anomaly
      
      ! Parameters
      DOUBLE PRECISION, PARAMETER :: PI = 4.D0*ATAN(1.D0)
      DOUBLE PRECISION, PARAMETER :: TWOPI = 2.D0*PI
      
      ! Local variables
      DOUBLE PRECISION :: mm
      DOUBLE PRECISION :: mguess
      DOUBLE PRECISION :: e       ! eccentric anomaly (K&VH Eqs. 2.10a & 2.10b, Fig. 2.3)
      INTEGER :: iter
      DOUBLE PRECISION :: de
      
! Executable statements

      mm = MOD(ma, TWOPI) ! mean anomaly modulo two pi
      e = mm              ! first guess for e

      ! Calcule mguess from current estimate of e.
      iter = 0
   10 mguess = e - eccen*SIN(e)
!
!     Iterate e using Newton-Raphson approach.  Stop when e increment is less
!     than 1.E-6*e or 10 iterations.
!
      de = (mguess - mm)/(1.D0 - eccen*COS(e))
      e = e - de
      IF(ABS(de) < 1.D-6) GOTO 20
      iter = iter + 1
      IF(iter >= 10) GOTO 20
      GOTO 10
!
!     Calculate true anomaly (theta) from eccentric anomaly (e).
!
   20 theta = ACOS((COS(e) - eccen)/(1.D0 - eccen*COS(e)))  ! K&VH Eq. 2.10a
      IF(mm > PI) theta = TWOPI - theta
      
   END SUBROUTINE Kepler

   SUBROUTINE ROT(vector, angle, axis)
 
! --------------------------------------------------------------------------------------------------
! Rotates "vector" by "angle" about "axis". (Generalized from K&VH Eqs. 2.23 & 2.24)
!
!  o angle is in radians  
!  o axis = 'X', 'Y', or 'Z'
! --------------------------------------------------------------------------------------------------
 
! Variables

      IMPLICIT NONE
      
      ! Arguments
      DOUBLE PRECISION, INTENT(INOUT) :: vector(3)
      DOUBLE PRECISION, INTENT(IN)    :: angle      ! radians
      CHARACTER(LEN=1), INTENT(IN)    :: axis 
      
      ! Local variables 
      INTEGER :: i, j 
      DOUBLE PRECISION :: cosa, sina
      
      ! Local arrays
      DOUBLE PRECISION :: R(3, 3), X(3)
      
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
      ELSE IF(AXIS == 'Z' .OR. AXIS == 'z') THEN
         i = 1
         j = 2
      ELSE
      	 PRINT *, "Subroutine ROT Error: Must rotate about x, y, or z axis."
      	 CALL EXIT
      END IF
      	! Default to rotation about the z axis.
      R(i, i) =  cosa
      R(j, j) =  cosa
      R(i, j) = -sina
      R(j, i) =  sina

      ! Perform rotation
      vector = MATMUL(R, vector)
      
   END SUBROUTINE ROT
      
END PROGRAM ephemeris
