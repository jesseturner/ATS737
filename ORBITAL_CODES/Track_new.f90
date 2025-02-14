PROGRAM Track

! This program reads an ephemeris file (EPHEMERIS.DAT), asks for a latitude and longitude
! of the observing location, and calculates the zenith, azimuth and elevation angles to 
! the satellite. (From K&VH Section 2.5.2)
! 
! NOTE: Slight errors are caused by the assumptions of (1) a spherical Earth and (2) that 
!       the elevation of the observing location is zero (i.e., it is at sea level).

! Modules

   USE Time_Module
   
! Variables
   
   IMPLICIT NONE
   
   ! Parameters
   DOUBLE PRECISION, PARAMETER :: Earth_Radius = 6371.D0 ! km
   DOUBLE PRECISION, PARAMETER :: PI = 4.D0*ATAN(1.D0)
   DOUBLE PRECISION, PARAMETER :: RCON = 180.D0/PI
   
   ! Local Variables
   CHARACTER(LEN = 50) :: sat_name
   CHARACTER(LEN = 255) :: input_line
   INTEGER :: ios
   INTEGER :: i
   INTEGER :: n
   REAL :: year, month, day, hour, minute, second
   REAL :: e_year, e_month, e_day
   DOUBLE PRECISION :: t, t0
   DOUBLE PRECISION :: lat0, lon0
   DOUBLE PRECISION :: zenith, azimuth, elevation
   DOUBLE PRECISION :: LEN_rD, LEN_re, LEN_rN, LEN_rH
   
   ! Local Arrays
   DOUBLE PRECISION, ALLOCATABLE :: time(:), lat(:), lon(:), r(:), x(:), y(:), z(:)
   DOUBLE PRECISION :: rs(3), re(3), rD(3), rN(3), rH(3)
   DOUBLE PRECISION :: rs_hat(3), re_hat(3), rD_hat(3), rN_hat(3), rH_hat(3)

! Executable statements

!   PRINT *, "RCON = ", RCON
!   PRINT *, "PI   = ", PI
   PRINT *

   ! Get observing point latitude and longitude
   PRINT '("Enter observing point latitude (-90 to 90) and longitude (-180 to 180) >"\)'
   READ *, lat0, lon0
!   PRINT *, "Latitude  = ", lat0
!   PRINT *, "Longitude = ", lon0
   PRINT *

      
   ! Read ephemeris file (always named EPHEMERIS.DAT)
   OPEN(UNIT = 1, FILE = "EPHEMERIS.DAT", STATUS = "OLD", IOSTAT = ios)
   IF(ios /= 0) THEN
      PRINT *, "EPHEMERIS.DAT not found"
      CALL EXIT
   END IF
   OPEN(UNIT = 2, FILE = "Track.txt", STATUS = "REPLACE")
   READ(1, '(A)') sat_name
   WRITE(2, '(A)') TRIM(sat_name)
   PRINT '(A)', TRIM(sat_name)
   READ(1, '(F4.0,2F2.0)') e_year, e_month, e_day
!   PRINT '(I4,2I2.2)', INT(e_year), INT(e_month), INT(e_day)
   READ(1, '(I5)') n
!   PRINT '(I5," LINES")', n
   READ(1, '(A)') input_line
!   PRINT '(A)', TRIM(input_line)
   ALLOCATE(time(n), lat(n), lon(n), r(n), x(n), y(n), z(n))
!   e_date = time_in_seconds(e_year, e_month, e_day, 0., 0., 0.)
   t0 = time_in_seconds(e_year, e_month, e_day, 0., 0., 0.)
   DO i = 1, n
      READ(1,'(7(1PE21.13))') time(i), lat(i), lon(i), r(i), x(i), y(i), z(i)
   END DO
   CLOSE(1)
   
   ! Location vector
   re(1) = Earth_Radius*COS(lat0/RCON)*COS(lon0/RCON)  ! K&VH Eq. 2.32
   re(2) = Earth_Radius*COS(lat0/RCON)*SIN(lon0/RCON)  ! K&VH Eq. 2.32
   re(3) = Earth_Radius*SIN(lat0/RCON)                 ! K&VH Eq. 2.32
   LEN_re = SQRT(DOT_PRODUCT(re, re))
   re_hat = re/LEN_re                                  ! K&VH Eq. 2.35a
   
   ! North vector
   rN_hat(1) = -SIN(lat0/RCON)*COS(lon0/RCON)          ! K&VH Eq. 2.34
   rN_hat(2) = -SIN(lat0/RCON)*SIN(lon0/RCON)          ! K&VH Eq. 2.34
   rN_hat(3) =  COS(lat0/RCON)                         ! K&VH Eq. 2.34
   
   ! Calculate zenith, azimuth, and elevation angles
   ! Print if zenith angle < 90 degrees
   WRITE(2,'(A)') " DATE & TIME (UTC)     ZENITH   AZIMUTH ELEVATION"
   PRINT '(A)', " DATE & TIME (UTC)     ZENITH   AZIMUTH ELEVATION"
   DO i = 1, n
   
      ! Satellite vector
      rs(1) = r(i)*COS(lat(i)/RCON)*COS(lon(i)/RCON)   ! K&VH Eq. 2.31
      rs(2) = r(i)*COS(lat(i)/RCON)*SIN(lon(i)/RCON)   ! K&VH Eq. 2.31
      rs(3) = r(i)*SIN(lat(i)/RCON)                    ! K&VH Eq. 2.31
      
      ! Difference vector
      rD = rs - re                                     ! bottom of K&VH page 34
      LEN_rD = SQRT(DOT_PRODUCT(rD, rD))
      rD_hat = rD/LEN_rD                               ! K&VH Eq. 2.35b
      
      ! Zenith angle
      zenith = RCON*ACOS(DOT_PRODUCT(re_hat, rD_hat))  ! K&VH Eq. 2.33
      elevation = 90. - zenith
      IF(zenith >= 90.) CYCLE                          ! Satellite below horizon
      
      ! Azimuth angle
      rH = rD - DOT_PRODUCT(re_hat, rD)*re_hat         ! K&VH Eq. 2.36
      LEN_rH = SQRT(DOT_PRODUCT(rH, rH))
      rH_hat = rH/LEN_rH
      azimuth = RCON*ACOS(DOT_PRODUCT(rN_hat, rH_hat)) ! K&VH Eq. 2.37
      IF(lon(i) < lon0) azimuth = 360.D0 - azimuth     ! See caution below K&VH Eq. 2.37
      
      ! Convert time and print result
      t = t0 + time(i)*86400.D0
      CALL traditional_time(t, year, month, day, hour, minute, second)
!      PRINT '(I4,4I3.2,F7.3)', INT(year), INT(month), INT(day), INT(hour), INT(minute), second
      PRINT 2, INT(month), INT(day), INT(year), INT(hour), INT(minute), INT(second), &
               zenith, azimuth, elevation
      WRITE(2,2) INT(month), INT(day), INT(year), INT(hour), INT(minute), INT(second), &
                 zenith, azimuth, elevation
    2 FORMAT(2(I2.2,"/"),I4,1X,2(I2.2,":"),I2.2, 3F10.3)
    
   END DO
   CLOSE(2)
   PRINT *
   PRINT *, "Output in file Track.txt"

END PROGRAM Track

