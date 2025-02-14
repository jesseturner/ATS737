MODULE Time_Module
!
!---------------------------------------------------------------------------------
!
! Written by Stan Kidder, CIRA/CSU
!
! Subroutines:
!
!    SUBROUTINE Calendar_Date(julday, y, m, d)
!    SUBROUTINE MODAY(doy, y, m, d)
!    SUBROUTINE traditional_time(time_in_seconds, year, month, day, hour, minute, second)
!    SUBROUTINE UNPACKTIMEGMT(timedate, iyr, imo, idy, ihr, imin, isec)
!
! Functions:
!
!    INTEGER FUNCTION Julian_Date(y, m, d)
!    INTEGER FUNCTION day_of_week(y, m, d)
!    INTEGER FUNCTION day_of_year(y, m, d)
!    REAL*8 FUNCTION time_in_seconds(year, month, day, hour, minute, second)
!    REAL*8 FUNCTION current_time()
!    REAL*8 FUNCTION UTC_difference()
!    REAL*8 FUNCTION UTC_correction(iyr, imo, idy, ihr)
!    REAL*8 FUNCTION convert_date(idate, itime)
!
!---------------------------------------------------------------------------------
!

   IMPLICIT NONE

! Global variables

   INTEGER, ALLOCATABLE :: jd_table(:)
   INTEGER :: nday(12) =(/31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
   INTEGER, PARAMETER :: year_min = -45 ! = 46 BCE
   INTEGER, PARAMETER :: year_max = 4000


CONTAINS

   INTEGER FUNCTION Julian_Date(y, m, d)
!
!-------------------------------------------------------------------------------
!
!     Function Julian_Date calculates the Julian Date (number of days since
!     1 Jan 4713 BCE) given the month, day, and year.  
!
!     y = year (YYYY). Negative years indicate BCE. There was no year zero, so 
!         don't try to use it. (It returns a big, bad number.)
!     m = month.
!     d = day.
!
!     NOTE: This routine works from 1 Jan 46 BCE (the beginning of the Julian 
!           Calendar) to 31 Dec 4000, when the Gregorian Calendar will get into 
!           trouble. The switch to the Gregorian Calendar is 14 Sep 1752 (the 
!           previous day is 2 Sep 1752), which is when the switch occurred in 
!           England and its colonies. Pope Gregory decreed that the calendar
!           switch on 15 Oct 1582 (the previous day was 5 Oct 1582), but it
!           took the rest of the world a while to catch up. The last places 
!           converted around 1925.
!
!-------------------------------------------------------------------------------
!
      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN)  :: y, m, d


! Variables

      INTEGER :: yy
      INTEGER :: i

! Executable statements

! Set up table

      IF(.NOT.ALLOCATED(jd_table)) CALL jd_table_setup() 

! Check date range

      Julian_Date = HUGE(Julian_Date) 
      IF(y == 0) RETURN  ! No year zero
      yy = y
      IF(yy < 0) yy = yy + 1  ! shift BCE years
      IF(yy < year_min .OR. yy > year_max) RETURN
      IF(yy == 1752) THEN
         IF(m == 9 .AND. d > 2 .AND. d < 14) RETURN  ! Missing 11 days, Julian to Gregorian
      END IF
      IF(m < 1 .OR. m > 12) RETURN
      nday(2) = 28
      IF(leap_year(y)) nday(2) = 29
      IF(d < 1 .OR. d > nday(m)) RETURN

! Calculate the Julian Date

      Julian_Date = jd_table(yy) + d - 1
      DO i = 1, m - 1
         Julian_Date = Julian_Date + nday(i)
      END DO
      IF(yy == 1752) THEN
         IF(m > 9 .OR. (m == 9 .AND. d > 13)) Julian_Date = Julian_Date - 11
      END IF

   END FUNCTION Julian_Date

   SUBROUTINE jd_table_setup()
!
!-------------------------------------------------------------------------------
!
!      This routine sets up a table (jd_table) which contains the Julian Date 
!      for the first day of each year in the table.
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Parameters

      INTEGER, PARAMETER :: N0 = 2448989 ! Julian Date for 1 Jan 1993
      INTEGER, PARAMETER :: Y0 = 1993

! Variables

      INTEGER :: yy

! Executable statements

      ALLOCATE(jd_table(year_min:year_max))
      jd_table(Y0) = N0   ! Julian Date for 1 Jan 1993
      DO yy = Y0 + 1, year_max
         jd_table(yy) = jd_table(yy - 1) + 365
         IF(yy < 2) THEN
            IF(leap_year(yy - 2)) jd_table(yy) = jd_table(yy) + 1   
         ELSE
            IF(leap_year(yy - 1)) jd_table(yy) = jd_table(yy) + 1   
         END IF                      
      END DO
      DO yy = Y0 - 1, year_min, -1
         IF(yy == 1752) THEN
            jd_table(yy) = jd_table(yy + 1) - 355
         ELSE
            jd_table(yy) = jd_table(yy + 1) - 365
            IF(yy < 1) THEN
               IF(leap_year(yy - 1)) jd_table(yy) = jd_table(yy) - 1
            ELSE
               IF(leap_year(yy)) jd_table(yy) = jd_table(yy) - 1
            END IF
         END IF      
      END DO        

   END SUBROUTINE jd_table_setup

   LOGICAL FUNCTION leap_year(y)
!
!-------------------------------------------------------------------------------
!
!     Returns .TRUE. if it is a leap year, .FALSE. if it isn't
!     NOTE: Assumes that BCE years have been shifted by 1 (i.e., year -1 is 
!           shifted to zero.
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Argument

      INTEGER, INTENT(IN) :: y

! Variables

      INTEGER :: yy

! Executable statements

      leap_year = .FALSE.
      yy = y
      IF(yy < 0) yy = yy + 1                 ! no year zero
      IF(MOD(yy, 4) == 0) leap_year = .TRUE.
      IF(yy <= 1752) RETURN
      IF(MOD(yy, 100) == 0 .AND. MOD(yy, 400) /= 0) leap_year = .FALSE. ! Gregorian correction

   END FUNCTION leap_year


   SUBROUTINE Calendar_Date(julday, y, m, d)
!
!-------------------------------------------------------------------------------
!
!     Subroutine Calendar_Date calculates the month, day, and year given the 
!     Julian Date.
!
!     julday = Julian Date
!     y      = year (YYYY).
!     m      = month.
!     d      = day.
!
!     NOTE: This routine works only from 1 March 1900 to 28 February 2100
!
!-------------------------------------------------------------------------------
!
      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN)  :: julday
      INTEGER, INTENT(OUT) :: y, m, d

! Parameters

      INTEGER, PARAMETER :: jd_min = 1704622 ! 1 Jan 46 BCE
      INTEGER, PARAMETER :: jd_max = 3182395 ! 31 Dec 4000 CE

! Variables

      INTEGER :: iy, n

! Executable statements

! Set up table

      IF(.NOT.ALLOCATED(jd_table)) CALL jd_table_setup() 
      
! Check date range

      IF(julday < jd_min .OR. julday > jd_max)  THEN
         y = HUGE(y)
         m = HUGE(m)
         d = HUGE(d)
         RETURN
      END IF

! Calculate year

      iy = INT(FLOAT(year_max + 1 - year_min)*(julday - jd_min)/(jd_max + 1 - jd_min)) + year_min
      DO
         IF(julday >= jd_table(iy)) EXIT
         iy = iy - 1
      END DO
      IF(iy < year_max) THEN
         DO
            IF(julday < jd_table(iy + 1)) EXIT
            iy = iy + 1
         END DO
      END IF

      y = iy
      IF(y < 1) y = y - 1  ! No year zero

! Set up date calculation

      nday(2) = 28
      IF(leap_year(iy)) nday(2) = 29

! Calculate date

      n = julday - jd_table(iy) + 1 ! day of year
      IF(iy == 1752 .AND. julday > 2361221) n = n + 11  ! Gregorian conversion
      DO m = 1, 12
         IF(n <= nday(m)) EXIT
         n = n - nday(m)
      END DO
      d = n

   END SUBROUTINE Calendar_Date

   INTEGER FUNCTION day_of_week(y, m, d)
!
!-------------------------------------------------------------------------------
!
! day_of week returns the day of the week: Sunday = 1, ..., Saturday = 7
!
!     y  = year (YYYY).
!     m  = month.
!     d  = day.
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN) :: y, m, d

! Variables

      INTEGER :: julday

! Executable Statements

      day_of_week = HUGE(day_of_week)
      julday = Julian_Date(y, m, d)
      IF(julday == HUGE(julday)) RETURN 
      day_of_week = MOD(julday + 1, 7) + 1

   END FUNCTION day_of_week

   INTEGER FUNCTION day_of_year(y, m, d)
!
!-------------------------------------------------------------------------------
!
!     Returns the day of year (1-366) given the month, day, and year.
!
!     y = year (YYYY)
!     m = month
!     d = day
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN)  :: y, m, d

! Variables

      INTEGER :: julday

! Executable statements

      day_of_year = HUGE(day_of_year)
      IF(y == HUGE(y)) RETURN
      julday = Julian_Date(y, m, d)
      IF(julday == HUGE(julday)) RETURN
      day_of_year = julday - Julian_Date(y, 1, 1) + 1

   END FUNCTION day_of_year


   SUBROUTINE MODAY(doy, y, m, d)
!
!-------------------------------------------------------------------------------
!
!     MODAY calculates the month and day given the year and day of year.
!
!     doy = day of year
!     y = year (YYYY)
!     m = month
!     m = day
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN)  :: doy, y
      INTEGER, INTENT(OUT) :: m, d

! Variables

      INTEGER :: julday, iy

! Executable statements

   julday = Julian_Date(y, 1, 1) + doy - 1
   CALL Calendar_Date(julday, iy, m, d)
   
   END SUBROUTINE MODAY


   REAL*8 FUNCTION time_in_seconds(year, month, day, hour, minute, second)
!
!-------------------------------------------------------------------------------
!
! Returns the number of seconds since midnight, January 1, 1993
! Leap seconds are not included.
!
! All input parameters are real*4 values. 
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      REAL, INTENT(IN) :: year, month, day, hour, minute, second

! Parameters

      INTEGER, PARAMETER :: jd_1_jan_1993 = 2448989

! Variables

      INTEGER :: jd

! Executable statements

      time_in_seconds = HUGE(time_in_seconds)

! Get Julian Day

      jd = Julian_Date(INT(year), INT(month), INT(day))
      IF(jd == HUGE(jd)) RETURN
      jd = jd - jd_1_jan_1993

      time_in_seconds = 86400.d0*jd + 3600.d0*MOD(INT(hour),24) &
                        + 60.d0*MOD(INT(minute),60) + MOD(second,60.)

   END FUNCTION time_in_seconds

   SUBROUTINE traditional_time(time_in_seconds, year, month, day, hour, minute, second)
!
!-------------------------------------------------------------------------------
!
! Returns year, month, day, hour, minute, and second given time_in_seconds
! Leap seconds are not included.
!
! time_in_seconds = the time in seconds since midnight, January 1, 1993 (real*8, input)
! year, month, day, hour, minute, second = real*4, output
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      REAL*8, INTENT(IN)  :: time_in_seconds
      REAL*4, INTENT(OUT) :: year, month, day, hour, minute, second

! Parameters

      INTEGER, PARAMETER :: jd_1_jan_1993 = 2448989

! Variables

      INTEGER :: jd, iyr, imo, idy
      REAL*8  :: part_of_day

! Executable statements

      IF(time_in_seconds == HUGE(time_in_seconds)) THEN   ! Check value
         year = HUGE(year)  
         month = HUGE(month)
         day = HUGE(day)
         hour = HUGE(hour)
         minute = HUGE(minute)
         second = HUGE(second)
         RETURN
      END IF

      jd = INT(time_in_seconds/86400.d0)           ! Calculate Julian Date
      part_of_day = time_in_seconds - jd*86400.d0
      IF(part_of_day < 0) THEN
         part_of_day = part_of_day + 86400.d0
         jd = jd - 1
      END IF       
      jd = jd + jd_1_jan_1993

      CALL Calendar_Date(jd, iyr, imo, idy)

      year = iyr
      month = imo
      day = idy
      hour = INT(part_of_day/3600.d0)
      minute = INT((part_of_day - 3600.d0*hour)/60.d0)
      second = SNGL(part_of_day - 3600.d0*hour - 60.d0*minute)
      IF(second >= 60.) THEN
         second = second - 60.
         minute = minute + 1.
         IF(minute >= 60.) THEN
            minute = minute - 60.
            hour = hour + 1.
            IF(hour >= 24.) THEN
               hour = hour - 24.
               jd = Julian_Date(iyr, imo, idy) + 1
               CALL Calendar_Date(jd, iyr, imo, idy)
               year = iyr
               month = imo
               day = idy
            END IF
         END IF
      END IF

   END SUBROUTINE traditional_time

   REAL*8 FUNCTION current_time()
!
!-------------------------------------------------------------------------------
!
! Returns the current time (UTC) as the number of seconds since midnight, 1 January 1993
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Variables

      INTEGER :: values(8)
      REAL    :: year, month, day, hour, minute, second

! Executable statements

      CALL DATE_AND_TIME(VALUES=values)
      year = values(1)
      month = values(2)
      day = values(3)
      hour = values(5)
      minute = values(6)
      second = values(7) + values(8)/1000.
      current_time = time_in_seconds(year, month, day, hour, minute, second) &
                     - values(4)*60.d0
!     print *, values

   END FUNCTION current_time

   REAL*8 FUNCTION UTC_difference()
!
!-------------------------------------------------------------------------------
!
! Returns local time minus UTC (in seconds) 
!
! Subtract this number from local time to get UTC
! Add this number to UTC to get local time 
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Variables

      INTEGER :: values(8)

! Executable statements

      CALL DATE_AND_TIME(VALUES=values)
      UTC_difference = values(4)*60.d0
!     print *, values

   END FUNCTION UTC_difference
   
   LOGICAL FUNCTION daylight_saving_time(iyr, imo, idy, ihr)

!
!-------------------------------------------------------------------------------
!
! Determines whether daylight saving time is in effect. Assumes that daylight  
! time starts at 2:00 AM on the first Sunday in April and ends at 2:00 am on the  
! last Sunday in October.
!
! NOTE: On August 8, 2005, President George W. Bush signed the Energy Policy 
!       Act of 2005. This Act changed the time change dates for Daylight Saving
!       Time in the U.S. Beginning in 2007, DST will begin on the second Sunday 
!       of March and end the first Sunday of November. The Secretary of Energy 
!       will report the impact of this change to Congress. Congress retains 
!       the right to revert the Daylight Saving Time back to the 2005 time 
!       schedule once the Department of Energy study is complete. THIS ROUTINE
!       WAS MODIFIED ON 27 FEB 2007 TO TAKE THE CHANGE INTO ACCOUNT.
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN) :: iyr, imo, idy, ihr

! Variables

      INTEGER :: jdy, iw, jsun

! Executable statements

      IF(iyr < 2007) THEN
         IF(imo < 4 .OR. imo > 10) THEN
            daylight_saving_time = .FALSE.
         ELSE IF(imo > 4 .AND. imo < 10) THEN
            daylight_saving_time = .TRUE.
         ELSE IF(imo == 4) THEN
            jdy = 1 ! find first Sunday in April
            DO 
               iw = day_of_week(iyr, imo, jdy)
               IF(iw == 1) EXIT
               jdy = jdy + 1
            END DO
            if(idy < jdy .OR. (idy == jdy .AND. ihr < 2)) THEN
               daylight_saving_time = .FALSE.
            ELSE
               daylight_saving_time = .TRUE.
            END IF
         ELSE
            jdy = 31 ! find last Sunday in October
            DO 
               iw = day_of_week(iyr, imo, jdy)
               IF(iw == 1) EXIT
               jdy = jdy - 1
            END DO
            IF(idy < jdy .OR. (idy == jdy .AND. ihr < 2)) THEN
               daylight_saving_time = .TRUE.
            ELSE
               daylight_saving_time = .FALSE.
            END IF
         END IF
      ELSE
         IF(imo < 3 .OR. imo > 11) THEN
            daylight_saving_time = .FALSE.
         ELSE IF(imo > 3 .AND. imo < 11) THEN
            daylight_saving_time = .TRUE.
         ELSE IF(imo == 3) THEN
            jsun = 0 ! find second Sunday in March
            jdy = 1
            DO 
               iw = day_of_week(iyr, imo, jdy)
               IF(iw == 1) THEN
                  jsun = jsun + 1
                  IF(jsun == 2) EXIT
               END IF
               jdy = jdy + 1
            END DO
            IF(idy < jdy .OR. (idy == jdy .AND. ihr < 2)) THEN
               daylight_saving_time = .FALSE.
            ELSE
               daylight_saving_time = .TRUE.
            END IF
         ELSE
            jdy = 1 ! find first Sunday in November
            DO 
               iw = day_of_week(iyr, imo, jdy)
               IF(iw == 1) EXIT
               jdy = jdy + 1
            END DO
            IF(idy < jdy .OR. (idy == jdy .AND. ihr < 2)) THEN
               daylight_saving_time = .TRUE.
            ELSE
               daylight_saving_time = .FALSE.
            END IF
         END IF
      END IF
   
   END FUNCTION daylight_saving_time   

   REAL*8 FUNCTION convert_date(idate, itime)

!
!-------------------------------------------------------------------------------
!
! Converts the date in CCYYDDD, HHMMSS format to time in seconds (since 1 Jan 1993)
! If CC = 0, the year is assumed to be between 1998 and 2097
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN) :: idate  ! YYDDD
      INTEGER, INTENT(IN) :: itime  ! HHMMSS

! Variables

!     REAL :: hour, minute, second
      INTEGER :: doy, iyr, imo, idy, ihr, imin, isec

! Executable statements

      iyr = idate/1000          ! get the year
      IF(iyr >= 0 .AND. iyr < 98) THEN
         iyr = iyr + 2000
      ELSE IF(iyr == 98 .OR. iyr == 99) THEN
         iyr = iyr + 1900
      END IF
      doy = MOD(idate,1000)     ! get the day of year
      CALL MODAY(doy, iyr, imo, idy)  ! convert to calendar date
      ihr = itime/10000               ! get the time
      imin = MOD(itime,10000)/100
      isec = MOD(itime,100)

      ! Convert to time in seconds since 1 Jan 1993
      !
      convert_date = time_in_seconds(FLOAT(iyr), FLOAT(imo), FLOAT(idy), &
                                     FLOAT(ihr), FLOAT(imin), FLOAT(isec))

   END FUNCTION convert_date

   REAL*8 FUNCTION UTC_correction(iyr, imo, idy, ihr)

!
!-------------------------------------------------------------------------------
!
! Returns the correction (seconds) to be added to the input (local) time to 
! compute UTC.
!
!-------------------------------------------------------------------------------
!

      IMPLICIT NONE

! Arguments

      INTEGER, INTENT(IN) :: iyr, imo, idy, ihr

! Variables

      INTEGER :: values(8)
      LOGICAL :: d_then, d_now

! Executable statements

      CALL DATE_AND_TIME(VALUES=values)
      UTC_correction = -values(4)*60.d0

      d_now  = daylight_saving_time(values(1), values(2), values(3), values(5))
      d_then = daylight_saving_time(iyr, imo, idy, ihr)

      IF(d_now .EQV. d_then) RETURN

      IF(d_now) THEN
         UTC_correction = UTC_correction + 3600.d0 
      ELSE
         UTC_correction = UTC_correction - 3600.d0
      END IF

   END FUNCTION UTC_correction

!   SUBROUTINE UNPACKTIMEGMT(timedate, iyr2, imo2, idy2, ihr2, imin2, isec2)
!!
!!-------------------------------------------------------------------------------
!!
!! Subroutine UNPACKTIMEGMT is a substitute for Portability Subroutine 
!! UNPACKTIMEQQ, which doesn't work correctly. It relies on the fact that
!! timedate is the number of seconds from 00:00:00 GMT 1 January 1970. 
!!
!! NOTE: An important difference between UNPACKTIMEQQ and UNPACKTIMEGMT is that 
!! the output time in UNPACKTIMEGMT is GMT, whereas the output time in 
!! UNPACKTIMEQQ is in the time zone of the computer (I think).
!!
!!-------------------------------------------------------------------------------
!!
!
!      USE IFPORT
!      IMPLICIT NONE
!
!! Arguments
!
!      INTEGER(KIND=4), INTENT(IN)  :: timedate
!      INTEGER(KIND=2), INTENT(OUT) :: iyr2, imo2, idy2, ihr2, imin2, isec2
!      
!! Local Variables
!
!      INTEGER(KIND=4) :: iyr4, imo4, idy4, ihr4, imin4, isec4
!      INTEGER(KIND=4) :: jd1970, jd1993
!      REAL(KIND=8) :: t, t_offset, t_daylight
!      REAL(KIND=4) :: year, month, day, hour, minute, second
!      INTEGER :: values(8)
!      LOGICAL :: d_now
!
!! Executable statements
!
!      ! timedate is supposed to be in number of seconds since 00:00:00 GMT 1 January 1970.
!      ! So, the first step is to convert timedate into HDFEOS time, which is the number of
!      ! seconds since 00:00:00 GMT 1 January 1993. 
!
!      jd1970 = Julian_Date(1970, 1, 1)
!      jd1993 = Julian_Date(1993, 1, 1)
!      t_offset = 86400.D0*(jd1993 - jd1970)
!      t = timedate
!      t = t - t_offset
!      
!      ! However, GETFILEINFOQQ returns timedate + 3600 seconds if daylight saving time is
!      ! *CURRENTLY* in effect. The following corrects this.
!      CALL DATE_AND_TIME(VALUES=values)
!      d_now  = daylight_saving_time(values(1), values(2), values(3), values(5))
!      t_daylight = 0.
!      IF(d_now) t_daylight = 3600.D0
!      t = t - t_daylight
!      
!      CALL Traditional_Time(t, year, month, day, hour, minute, second)
!      iyr2 = NINT(year)
!      imo2 = NINT(month)
!      idy2 = NINT(day)
!      ihr2 = NINT(hour)
!      imin2 = NINT(minute)
!      isec2 = NINT(second)     
!   
!   END SUBROUTINE UNPACKTIMEGMT

END MODULE Time_Module