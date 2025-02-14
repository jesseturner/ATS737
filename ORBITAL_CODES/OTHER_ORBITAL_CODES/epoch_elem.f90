PROGRAM EPOCH_ELEM
!
! Program to decode two-line elements.  
! Output is a TAB - delimited file suitable for import into Excel.  
!
! Modules

   USE Orbits_Module
   USE Time_Module
      
! Variables

   IMPLICIT NONE
   
   ! Parameters
   CHARACTER(LEN = 1), PARAMETER :: TAB = CHAR(9)
   
   ! Local Variables
   CHARACTER(LEN = 255) :: TLE_file
   CHARACTER(LEN = 255) :: output_fn
   INTEGER :: kyr, kmo, kdy, khr, kmin
   DOUBLE PRECISION :: ksec
   DOUBLE PRECISION :: ktime
   INTEGER :: iostat
   INTEGER :: nsats
   CHARACTER(LEN = 80) :: line0, line1, line2
   CHARACTER(LEN = 80)  :: sat_name
	CHARACTER(LEN = 7) :: ID
	INTEGER :: ehr, emin
	DOUBLE PRECISION :: esec
	DOUBLE PRECISION :: dt
	DOUBLE PRECISION :: ma, ra, ap
	DOUBLE PRECISION :: theta          ! true anomaly (radians)
	DOUBLE PRECISION :: true_anom      ! true anomaly (degrees)
	DOUBLE PRECISION :: arg_lat        ! argument of latitude (degrees)
	DOUBLE PRECISION :: lat, lon, r, x, y, z
	INTEGER :: k
	DOUBLE PRECISION :: gra
   
! Executable Statements

   ! Open TLE file.
   WRITE(*, '(" Name of the two-line element file >"\)')
   READ(*, '(A)') TLE_file
   OPEN(UNIT = 1, FILE = TLE_file, STATUS = 'OLD', IOSTAT = iostat)
   IF(iostat /= 0) THEN
      PRINT *, "File not found."
      STOP
   END IF
   
   ! Form output file name
!   k = INDEX(TLE_file, '\', BACK=.TRUE.)  ! Strip direcory
   k = INDEX(TLE_file, '"\"', BACK=.TRUE.)  ! Strip direcory
   output_fn = TLE_file(k+1:)
   k = INDEX(output_fn, ".", BACK=.TRUE.) ! replace file type (.xxx) with  ".tsv"
   IF(k == 0) THEN
      output_fn = TRIM(output_fn) // "_epoch.tsv"
   ELSE
      output_fn = output_fn(1:k-1) // "_epoch.tsv"
   END IF

   ! Open output file and write header
   OPEN(UNIT = 2, FILE = output_fn, STATUS = 'REPLACE')
   
   WRITE(2, 3) 'SATELLITE'//TAB, 'ID'//TAB, 'EPOCH_DATE/TIME'//TAB, TAB, TAB, TAB, TAB, TAB,        &
               'SEMI-MAJOR_AXIS'//TAB, 'INCLINATION'//TAB, 'ECCENTRICITY'//TAB,                     &
               'RIGHT_ASCENSION'//TAB, 'MOTION'//TAB, 'ARG._OF_PERIGEE'//TAB, 'MOTION'//TAB,        &            
               'MEAN_ANOMALY'//TAB, 'MOTION'//TAB, 'NODAL_PERIOD'//TAB, 'ARG._OF_LATITUDE'//TAB, &
               'GEOCENTRIC_LATITUDE'//TAB, 'LONGITUDE'//TAB, 'R'//TAB, 'X'//TAB, 'Y'//TAB, 'Z'
               
                   
   WRITE(2, 3) TAB, TAB, 'year'//TAB, 'month'//TAB, 'day'//TAB, 'hour'//TAB, 'minute'//TAB, &
               'second'//TAB, 'km'//TAB, 'deg'//TAB, 'unitless'//TAB,   'deg'//TAB,         &
               'deg/day'//TAB, 'deg'//TAB,  'deg/day'//TAB, 'deg'//TAB, 'deg/day'//TAB,     &
               'min'//TAB, 'deg'//TAB, 'deg'//TAB, 'deg'//TAB, 'km'//TAB, 'km'//TAB, 'km'//TAB, 'km'
               
 3 FORMAT(35A)
 
   ! Process the file
   nsats = 0
   DO
   
      ! Read next TLE
      READ(UNIT=1, FMT='(A80)', END=50) line0
      READ(UNIT=1, FMT='(A80)', END=50) line1
      DO
         READ(UNIT=1, FMT='(A80)', END=50) line2
         IF(line1(1:2) == "1 " .AND. line2(1:2) == "2 ") EXIT
         line0 = line1
         line1 = line2
      END DO
      
      ! Decode the TLE
      nsats = nsats + 1
      CALL DECODEL(line1, line2)
	   ID = line1(10:16)
      ehr = INT(24*etime)
      emin = INT((24*etime - ehr)*60)
      esec = ((24*etime - ehr)*60 - emin)*60.
      
      ! Calculate orbital elements at epoch time 
      ma = ma0
      ra = ra0
      ap = ap0
   
      ! Calculate argument of latitude
      CALL KEPLER(ma/RCON, eccen, theta)
      true_anom = theta*RCON
      arg_lat = MOD(true_anom + ap, 360.D0)

      ! Calculate satellite position at epoch time
      CALL POSIT(eyr, emo, edy, etime, lat, lon, r, x, y, z)
      
      ! Extract satellite name
      sat_name = line0
      IF(sat_name(1:2) == "0 ") sat_name =  sat_name(3:)
      k = INDEX(sat_name, "[")
      IF(k > 1) sat_name = sat_name(:k-1)

	   WRITE(2, 5) TRIM(sat_name), TAB, TRIM(ID), TAB, &
	               eyr, TAB, emo, TAB, edy, TAB, ehr, TAB, emin, TAB, esec, TAB, &
                  a, TAB, inc, TAB, eccen, TAB, ra0, TAB, radot, TAB, ap0, TAB,   &
                  apdot, TAB, ma0, TAB, madot, TAB, ttilda*1440.D0, TAB, arg_lat, TAB, &
                  lat, TAB, lon, TAB, r, TAB, x, TAB, y, TAB, z
               
    5 FORMAT(2(A, A1), I4, A1, 4(I2, A1), F6.3, A1, F9.3, A1,  &
             F7.3, A1, F9.7, A1, 3(0PF8.3, A1, 1PE14.6, A1), 4(0PF10.3, A1), 4(0PF11.3, A1))

!      ! Write record
!	   WRITE(2, 5) TRIM(sat_name), TAB, TRIM(ID), TAB, &
!	               eyr, TAB, emo, TAB, edy, TAB, ehr, TAB, emin, TAB, esec, TAB,       &
!	               kyr, TAB, kmo, TAB, kdy, TAB, khr, TAB, kmin, TAB, ksec, TAB,       &
!                  a, TAB, inc, TAB, eccen, TAB, ra, TAB, radot, TAB, ap, TAB,         &
!                  apdot, TAB, ma, TAB, madot, TAB, ttilda*1440.D0, TAB, arg_lat, TAB, &
!                  lat, TAB, lon, TAB, r, TAB, x, TAB, y, TAB, z
!               
!    5 FORMAT(2(A, A1), 2(I4, A1, 4(I2, A1), F6.3, A1), F9.3, A1, F7.3, A1, F9.7, A1,  &
!             3(F8.3, A1, F10.3, A1), 4(F10.3, A1), 4(F11.3, A1))
          
   END DO ! next TLE
 
!   ! Process the file
!   nsats = 0
!20 READ(1, 4, END = 1000) line0
!   READ(1, 4, END = 1000) line1
!30 READ(1, 4, END = 1000) line2
! 4 FORMAT(A80)
!   IF(line1(1:2) /= "1 " .OR. line2(1:2) /= "2 ") THEN
!      line0 = line1
!      line1 = line2
!      GOTO 30
!   END IF
!   nsats = nsats + 1
!   sat_name = line0
!   IF(sat_name(1:2) == "0 ") sat_name = sat_name(3:)
!   k = INDEX(sat_name, "[")
!   IF(k > 1) sat_name = sat_name(1:k-1)
!   CALL DECODEL(line1, line2)
!	ID = line1(10:16)
!   ehr = INT(24*etime)
!   emin = INT((24*etime - ehr)*60)
!   esec = ((24*etime - ehr)*60 - emin)*60.
!!   dt = TIMEDIF(eyr, emo, edy, etime, kyr, kmo, kdy, ktime)
!!   ma = mod(ma0 + madot*dt, 360.d0)
!!   IF(ma < 0.) ma = ma + 360.D0
!!   ra = mod(ra0 + radot*dt, 360.d0)
!!   IF(ra < 0.) ra = ra + 360.D0
!!   ap = mod(ap0 + apdot*dt, 360.d0)
!!   IF(ap < 0.) ap = ap + 360.D0
!!   PRINT *, ttilda, revs, ttilda*revs
!   CALL POSIT(eyr, emo, edy, etime, lat, lon, r, x, y, z)
!
!	WRITE(2, 5) TRIM(sat_name), TAB, TRIM(ID), TAB, &
!	            eyr, TAB, emo, TAB, edy, TAB, ehr, TAB, emin, TAB, esec, TAB, &
!               a, TAB, inc, TAB, eccen, TAB, ra0, TAB, radot, TAB, ap0, TAB,   &
!               apdot, TAB, ma0, TAB, madot, TAB, ttilda*1440.D0, TAB, &
!               lat, TAB, lon, TAB, r, TAB, x, TAB, y, TAB, z
!               
! 5 FORMAT(2(A, A1), I4, A1, 4(I2, A1), F6.3, A1, F9.3, A1,  &
!          F7.3, A1, F9.7, A1, 3(0PF8.3, A1, 1PE14.6, A1), 3(0PF8.3, A1), 4(0PF10.3, A1))
!          
!        GOTO 20

50 CLOSE(1)
   CLOSE(2)
   PRINT *,  nsats, " satellites processed"
   PRINT *,  "Data written to file " // TRIM(output_fn)
        
END PROGRAM EPOCH_ELEM

