MODULE auxiliary

	! Date		Programmer		Changes
	! ====		==========		=======
	! Mar/13 	E. Martinez		Original
	
	! Description:
	!
	! Auxiliary subroutines for the DMST main program:
	!	1.- Simpson's rule for torque integration.
	!	2.- A parsing routine to fill drag and lift coefficient tables.
	!	3.- A routine that collects a csv line's values into an array.
	!	4.- A 2D linear interpolator routine to obtain CD and CL.
	!   5.- A parser for the dat table.

	IMPLICIT NONE


	CONTAINS


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 1.- Simpson's rule for torque integration.
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	SUBROUTINE integrate_torque (Q, NZ, NT, dz, dt, cq)
		IMPLICIT NONE

		! Data dictionary:
		INTEGER, INTENT(IN) :: NZ 					! Number of vertical levels
		INTEGER, INTENT(IN) :: NT					! Number of sectors
		REAL, INTENT(IN) :: dz 						! Height delta [-]
		REAL, INTENT(IN) :: dt 						! Sector delta [-]
		REAL, INTENT(OUT) :: cq 					! Torque coefficient output

		! The Q array's bounds have been set according to the C notation
		! which is from 0 to N (non-inclusive):
		REAL, INTENT(IN), DIMENSION (0:NZ-1, 0:NT-1) :: Q

		REAL, DIMENSION (0:NT-1) :: QTHETA			! Array of sector torque
		REAL :: acc 								! Temporary accumulator
		INTEGER :: i, j								! 2D array indices
		INTEGER :: l, m, n 							! Auxiliary indices

		! Get each angular sector torque vertically:
		DO j = 0, NT - 1

			acc = 0.	! Initialize accumulator to zero

			DO  i = 1, (NZ - 1) / 2

				! Dummy indices:
				l = (2 * i) - 2		! f(a)
				m = (2 * i) - 1		! f((a + b) / 2)
				n = (2 * i) 		! f(b)

				! Increase the accumulator value:
				acc = acc + Q(l, j) + (4.0 * Q(m, j)) + Q(n, j)

			END DO

			QTHETA(j) = (dz / 3.0) * acc

		END DO

		! Integrate total angular torque:
		acc = 0.

		DO  j = 1, (NT - 1) / 2

			! Dummy indices:
			l = (2 * j) - 2		! f(a)
			m = (2 * j) - 1		! f((a + b) / 2)
			n = (2 * j) 		! f(b)

			! Increase the accumulator value:
			acc = acc + QTHETA(l) + (4.0 * QTHETA(m)) + QTHETA(n)

		END DO

		! The total angular torque is multiplied by dt / 3.0
		cq = acc * (dt / 3.0)

	END SUBROUTINE integrate_torque


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 2.- A parsing routine to fill drag and lift coefficient tables.
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	SUBROUTINE parse_csv (filename, ARRAY, NR, NC)
		IMPLICIT NONE

		! Data dictionary (subroutine arguments):
		CHARACTER(LEN=33), INTENT(IN) :: filename 		! File name with .csv
		INTEGER, INTENT(IN) :: NR, NC					! Rows and columns
		REAL, INTENT(OUT), DIMENSION(NR, NC) :: ARRAY 	! 2D array to be filled

		! Data dictionary (subroutine temporary variables):
		INTEGER :: stat 								! I/0 status
		CHARACTER(LEN=333) :: line 						! Current file line
		CHARACTER(LEN=66) :: error            			! Error message
		INTEGER :: row, col 							! Array parser indices
		CHARACTER(LEN=33), DIMENSION(NC) :: tokens	    ! Array of numbers

		! Open the file (do not use UNIT=6):
		OPEN(UNIT=10, FILE=filename, STATUS='OLD', ACTION='READ')

		! Begin reading the columns row by row:
		row = 1

		DO

			READ(10, '(A)', IOSTAT=stat) line			! Read the line
			IF (stat .NE. 0) EXIT 						! Exit if error occurs

			CALL tokenize(line, tokens, NC)				! Fill the tokens

			DO col = 1, NC

				! Insert tokens into the array columns:
				READ(tokens(col), '(F9.6)', IOSTAT=stat, IOMSG=error) ARRAY(row, col)

				IF (stat .NE. 0) WRITE (*, *) 'Error: ', error

			END DO

			row = row + 1								! Move to the next row

			IF (row .GT. NR) EXIT

		END DO

		CLOSE(UNIT=10)

	END SUBROUTINE parse_csv


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 3.- A routine that collects a csv line's values into an array.
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

    SUBROUTINE tokenize (line, tokens, NCOLS)
		IMPLICIT NONE

		! Data dictionary (subroutine arguments):
        CHARACTER(LEN=333), INTENT(IN) :: line
        CHARACTER(LEN=33), DIMENSION(NCOLS), INTENT(OUT) :: tokens
        INTEGER, INTENT(IN) :: NCOLS
        INTEGER, PARAMETER :: N = 33

        ! Data dictionary (procedure):
        INTEGER :: chcount          	! Counter index for temporary string ch
        INTEGER :: chpos 				! Each character's position in the line
        INTEGER :: tokenidx 			! Token index, each token is a value
        INTEGER :: length				! Length of the line of characters
        INTEGER :: tsidx				! Temporary string index
        CHARACTER(LEN=1) :: ch 			! Current character of the line
        CHARACTER(LEN=N) :: tempstring  ! Temporary string

        length = LEN_TRIM(line) 		! Get the length of the line w/o spaces
        chcount = 1						! Initialize character counter
        tokenidx = 1 					! Initialize tokens array index
        chpos = 1 						! Initialize character position

		! Purge the temporary string first:
		DO tsidx = 1, N
			tempstring(tsidx:tsidx) = ""
		END DO

		! Loop while the current position number is less than the length
        ! of the line of characters.

        DO WHILE ( (chpos .LE. length) )

            READ(line(chpos:chpos), '(A1)') ch

			! If a comma is found increase the token index:
            IF (ch .EQ. ',') THEN

				! Pass the character string to the current token,
				! increase the tokens array index and reset
				! temporary string character count:
				tokens(tokenidx) = tempstring
                tokenidx = tokenidx + 1
                chcount = 1

                ! Purge the temporary string:
                DO tsidx = 1, N
					tempstring(tsidx:tsidx) = ""
                END DO

            ! If it's the last character, catch it in the temporaty string
            ! and then pass the temporary string to the current token:
            ELSE IF (chpos .EQ. length) THEN

                tempstring(chcount:chcount) = ch
                tokens(tokenidx) = tempstring

			! If not a comma, store the character in the temporary string
			! array, and increase the temporary string character count:
            ELSE

                tempstring(chcount:chcount) = ch
                chcount = chcount + 1

            END IF

			! Move on to the next character in the line:
            chpos = chpos + 1

        END DO

    END SUBROUTINE tokenize


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 4.- A 2D linear interpolator routine to obtain CD and CL.
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	SUBROUTINE interpolate_table (aa_in, re_in, C, A, R, NR, NC, val)
		IMPLICIT NONE

		! Data dictionary (subroutine arguments):
		REAL, INTENT(IN) :: aa_in 						! Angle of attack input
		REAL, INTENT(IN) :: re_in 						! Local Re input
		INTEGER, INTENT(IN) :: NR 						! Number of rows
		INTEGER, INTENT(IN) :: NC 						! Number of columns
		REAL, INTENT(IN), DIMENSION(0:NR-1) :: A 		! Table of angles
		REAL, INTENT(IN), DIMENSION(0:NC-1) :: R 		! Table of Re numbers

		! Table of coefficients:
		REAL, INTENT(IN), DIMENSION(0:NR-1, 0:NC-1) :: C

		REAL, INTENT(OUT) :: val 						! Output value

		! Data dictionary ( dummy variables):
		REAL :: todeg = 180. / (4. * ATAN(1.)) 			! To degrees factor
		REAL :: aa 										! Angle in degrees
		REAL :: re 										! Reynolds in 1E-6

		INTEGER :: i0 = 0 								! Row index 0
		INTEGER :: i1 = 0 								! Row index 1
		INTEGER :: j0 = 0 								! Column index 0
		INTEGER :: j1 = 0								! Column index 1

		REAL :: aax = 0. 								! aa location [0:1]
		REAL :: rex = 0.								! re location [0:1]

		REAL :: rej0 = 0.								! re @ aa, < re
		REAL :: rej1 = 0.								! re @ aa, > re

		INTEGER :: k 									! Array index

		! Convert aa to degrees:
		aa = aa_in * todeg

		! Normalize Re:
		re = re_in / 1E6

		! Find i0 and i1
		DO k = 0, NR - 1

			IF ( (aa .GE. A(k)) .AND. (aa .LE. A(k + 1)) ) THEN

				i0 = k
				i1 = k + 1
				aax = ( aa - A(k) ) / ( A(k + 1) - A(k) )

				EXIT

			END IF

		END DO

		! Find j0 and j1
		DO k = 0, NC - 1

			IF ( (re .GE. R(k)) .AND. (re .LE. R(k + 1)) ) THEN

				j0 = k
				j1 = k + 1
				rex = ( re - R(k) ) / ( R(k + 1) - R(k) )

				EXIT

			END IF

		END DO

		! Coefficients at aa:
		rej0 = C(i0, j0) + aax * ( C(i1, j0) - C(i0, j0) )
		rej1 = C(i0, j1) + aax * ( C(i1, j1) - C(i0, j1) )

		! Coefficient at re:
		val = rej0 + rex * (rej1 - rej0)

	END SUBROUTINE interpolate_table


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 5.- A parser for the angle of attack table
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	SUBROUTINE parse_dat (f, ARR, L)
		IMPLICIT NONE

		! Data dictionary (subroutine arguments):
		INTEGER, INTENT(IN) :: L 					! Length of the array
		REAL, INTENT(OUT), DIMENSION(L) :: ARR		! Array of reals
		CHARACTER(LEN=33), INTENT(IN) :: f 			! File name

		! Data dictionary (procedure):
		INTEGER :: i 								! Array index
		INTEGER :: stat								! I/O status
		CHARACTER(LEN=66) :: msg 					! Error message

		! Open the file, check for errors and open it
		OPEN (UNIT=11, FILE=f, STATUS='OLD', ACTION='READ', IOSTAT=stat, IOMSG=msg)

		! Check if there are no errors after opening the file:
		IF (stat .EQ. 0) THEN

			! OPEN was ok. Proceed to read values.

			! Array index begins at 1:
			i = 1
			DO
				READ (11, *, IOSTAT=stat) ARR(i) 	! Get current line
				i = i + 1 							! Move to the next index

				IF ( stat .NE. 0) EXIT				! Exit do loop

			END DO

			! An error ocurred, we are here because the do loop was interrupted.

			! stat > 0, an error occurred:
			IF (stat .GT. 0) THEN
				WRITE (*, *) "An error occurred: "
			END IF

		ELSE

			WRITE (*, 100) stat
			100 FORMAT ('Error opening file: IOSTAT = ', I6)
			WRITE (*, 200) TRIM(msg)
			200 FORMAT (A)

		END IF

		CLOSE (UNIT=11)

	END SUBROUTINE parse_dat
	

END MODULE auxiliary 
