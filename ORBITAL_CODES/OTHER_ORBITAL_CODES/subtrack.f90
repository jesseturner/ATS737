PROGRAM SUBTRACK
!
! This program plots a satellite subtrack. Input data
! are from file EPHEMERIS.DAT.
!
   USE Graphics_Module
   
! Variables

   IMPLICIT NONE
   
   ! Parameters
   INTEGER, PARAMETER :: NX = 1100
   INTEGER, PARAMETER :: NY = INT(0.55*NX)

   ! Local Variables
   CHARACTER(LEN=50) :: sat_name
   CHARACTER(LEN=80) :: label
   CHARACTER(LEN=255) :: input_line
   LOGICAL :: penup
   REAL :: u, v
   REAL :: uu, vv
   REAL :: umin, umax, vmin, vmax
   REAL :: ulast, vlast
   REAL :: u1, v1, u2, v2
   REAL :: time
   INTEGER :: nseg, idls, idrs, npts
   INTEGER :: i, m
   INTEGER :: nvec
   INTEGER :: n, nn
   INTEGER :: nlines
      
   ! Local Arrays
   REAL :: pts(200)
      
! Executable statements

   !  Open ephemeris file
   OPEN(UNIT = 1, FILE = "EPHEMERIS.DAT", STATUS = "OLD")
   READ(1, "(A50)") sat_name
   READ(1, "(A80)") label
   READ(1, "(I5)") nlines
   READ(1, "(A)") input_line

   !  Set up plot file
   CALL SETUP(NX, NY)
   umin = -180.
   umax = 180.
   vmin = -90.
   vmax = 90.
   CALL SET(.05, .95, .10, .90, umin, umax, vmin, vmax)
   CALL SET_BACKGROUND(255)

   ! Plot world map
   nvec = 0
   CALL SET_COLOR(127)
   CALL SET_WIDTH(1)
   OPEN(UNIT = 27, FILE = "Continental_Outlines.DAT", STATUS = "OLD")
   DO
      IF(EOF(27)) EXIT
      READ(27, '(4I6/(8F10.4))') nseg, idls, idrs, npts, (pts(m), m = 1, 2*npts)
      penup = .TRUE.
      DO i = 1, npts
         v = pts(2*i-1)
         u = pts(2*i)
         IF(u >= umin .and. u <= umax .and. v >= vmin .and. v <= vmax) THEN
            IF(i == 1) THEN
               CALL FRSTPT(u, v)
               penup = .FALSE.
            ELSE IF(penup) THEN
               CALL WINDO(ulast, vlast, u, v, umin, vmin, umax, vmax, uu, vv)
               CALL LINE(uu, vv, u, v)
               penup = .FALSE.
            ELSE
               CALL VECTOR(u, v)
               nvec = nvec + 1
            END IF
         ELSE IF(i > 1) THEN
            IF(.NOT.penup) THEN
               CALL WINDO(u, v, ulast, vlast, umin, vmin, umax, vmax, uu, vv)
               CALL VECTOR(uu, vv)
               nvec = nvec + 1
               penup = .TRUE.
            END IF
         END IF
         ulast = u
         vlast = v
      END DO
   END DO
   CLOSE(27)
 
!   ! Plot important latitudes
!   CALL DASHD(#CCCC)
!   CALL LINED(-180., 66.55, 180., 66.55)    ! Arctic Circle
!   CALL LINED(-180., 23.45, 180., 23.45)    ! Tropic of Cancer
!   CALL LINED(-180., 0., 180., 0.)          ! Equator
!   CALL LINED(-180., -23.45, 180., -23.45)  ! Tropic of Capricorn
!   CALL LINED(-180., -66.55, 180., -66.55)  ! Antarctic Circle

   ! Plot lat/lon
   CALL DASHD(#CCCC)
   DO i = -150, 150, 30
      CALL LINED(REAL(i), -90., REAL(i), 90.)
   END DO
   DO i = -60, 60, 30
      CALL LINED(-180., REAL(i), 180., REAL(i))
   END DO
   CALL SET_COLOR(248)
   CALL PERIM(1, 1, 1, 1)
   
   ! Plot subtrack
   ulast = 9999.
   CALL SET_WIDTH(3)
   DO nn = 1, nlines
      READ(1, '(3F21.0)') time, v, u
      IF(ABS(u - ulast) > 180.) THEN
         IF(nn == 1) THEN
            CALL FRSTPT(u, v)
         ELSE
            CALL DATE_LINE(u, v, ulast, vlast, u1, v1, u2, v2)
            CALL VECTOR(u1, v1)
            CALL FRSTPT(u2, v2)
            CALL VECTOR(u, v)
         END IF
      ELSE
         CALL VECTOR(u, v)
      END IF
      IF(nn == 3) THEN
         CALL ARROW_HEAD(u, v, ulast, vlast, u1, v1, u2, v2)
         CALL VECTOR(u1, v1)
         CALL FRSTPT(u2, v2)
         CALL VECTOR(u, v)
      END IF
      ulast = u
      vlast = v
   END DO
   CALL SET_WIDTH(1)
   
   ! Label the plot
   CALL SET(0., 1., 0., 1., 0., 800., 0., 440.)
   n = LEN_TRIM(sat_name)
   DO i = 1, n
      IF(sat_name(i:i) == "_") sat_name(i:i) = " "
   END DO
   CALL PWRIT(400., 420., sat_name, n, 2, 0, 0)
   CALL PWRIT(40., 20., label(10:28), 19, 1, 0, -1)
   CALL PWRIT(400., 20., label(30:48), 19, 1, 0, 0)
   CALL PWRIT(760., 20., label(50:67), 18, 1, 0, 1)

   ! Finish
   CLOSE(1)
   CALL FRAME("SUBTRACK.TIF")
   PRINT *,  nlines, " vectors plotted."
   PRINT *,  "Output written to file SUBTRACK.TIF"
      
CONTAINS

   SUBROUTINE DATE_LINE(U, V, ULAST, VLAST, U1, V1, U2, V2)
   
! Handle line segments that cross the Date Line

! Variables

      IMPLICIT NONE
      
      ! Arguments
      REAL, INTENT(IN)  :: u, v, ulast, vlast
      REAL, INTENT(OUT) :: u1, v1, u2, v2
      
      ! Local Variables
      REAL :: d, f
      
! Executable Statements
         
      IF(u > ulast) THEN
         d = 360. - (u - ulast)
         f = (180. - u)/d
         v1 = (1 - f)*v + f*vlast
         v2 = v1
         u1 = -180.
         u2 = 180.
      ELSE
         d = 360. - (ulast - u)
         f = (180. - ulast)/d
         v1 = (1 - f)*vlast + f*v
         v2 = v1
         u1 = 180.
         u2 = -180.
      END IF
      
   END SUBROUTINE DATE_LINE

   SUBROUTINE ARROW_HEAD(u, v, ulast, vlast, u1, v1, u2, v2)
   
! Calculates the arrowhead to draw near the beginning of the ephemeris

! Variables

      IMPLICIT NONE
      
      ! Arguments
      REAL, INTENT(IN)  :: u, v, ulast, vlast
      REAL, INTENT(OUT) :: u1, v1, u2, v2
      
      ! Parameters
      REAL, PARAMETER :: RCON = 45./ATAN(1.)
      REAL, PARAMETER :: ANGLE = 150.
      INTEGER, PARAMETER :: L = 5
      
      ! Local Variables
      REAL :: s, c, x, y, d
      
! Executable Statements

      s = SIN(ANGLE/RCON)
      c = COS(ANGLE/RCON)
      x = u - ulast
      y = v - vlast
      d = SQRT(x*x + y*y)
      x = L*x/d
      y = L*y/d
      u1 = u + c*x - s*y
      v1 = v + s*x + c*y
      u2 = u + c*x + s*y
      v2 = v - s*x + c*y
      
   END SUBROUTINE ARROW_HEAD

END PROGRAM SUBTRACK

