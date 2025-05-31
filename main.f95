PROGRAM main
	! Date		Programmer		Changes
	! ====		==========		=======
	! Mar/30 	E. Martinez		Original
	
	! Description:
	!
	! Main Double Multiple Streamtube Model.

	USE auxiliary

	IMPLICIT NONE


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 1.- Data dictionary
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	! Pi (constant):
	REAL, PARAMETER :: PI = 4.0 * ATAN(1.0)

	! Convergence residual (constant):
	REAL, PARAMETER :: CR = 0.01

	! Air's kinematic viscosity:
	REAL, PARAMETER :: NU = 1.5E-5

	! Arrays' dimensions:
	INTEGER, PARAMETER :: ROWS = 117			! Number of rows
	INTEGER, PARAMETER :: COLS = 11 			! Number of columns
	INTEGER, PARAMETER :: NZ = 21 				! Number of vertical levels
	INTEGER, PARAMETER :: NT = 35 				! Number of tubes

	! Arrays (from files):
	REAL, DIMENSION(ROWS) :: AA_TABLE 			! Angle of attack data array
 	REAL, DIMENSION(COLS) :: RE_TABLE 			! Re data array
	REAL, DIMENSION(ROWS, COLS) :: CL_TABLE 	! CL data array
	REAL, DIMENSION(ROWS, COLS) :: CD_TABLE		! CD data array

	! Arrays (procedure):
	REAL, DIMENSION(NZ, NT) :: UU 		! Upwind:	Induction factor complement
	REAL, DIMENSION(NZ, NT) :: UD 		! Downwind: Induction factor complement
	REAL, DIMENSION(NZ, NT) :: WU 		! Upwind: 	Dimensionless relative speed
	REAL, DIMENSION(NZ, NT) :: WD 		! Downwind: Dimensionless relative speed
	REAL, DIMENSION(NZ, NT) :: AAU 		! Upwind: 	Angle of attack
	REAL, DIMENSION(NZ, NT) :: AAD 		! Downwind: Angle of attack
	REAL, DIMENSION(NZ, NT) :: REU 		! Upwind: 	Local Reynolds number
	REAL, DIMENSION(NZ, NT) :: RED 		! Downwind: Local Reynolds number
	REAL, DIMENSION(NZ, NT) :: CLU 		! Upwind: 	Lift coefficient
	REAL, DIMENSION(NZ, NT) :: CLD 		! Downwind: Lift coefficient
	REAL, DIMENSION(NZ, NT) :: CDU 		! Upwind: 	Drag coefficient
	REAL, DIMENSION(NZ, NT) :: CDD 		! Downwind: Drag coefficient
	REAL, DIMENSION(NZ, NT) :: CNU 		! Upwind: 	Normal force coefficient
	REAL, DIMENSION(NZ, NT) :: CND 		! Downwind: Normal force coefficient
	REAL, DIMENSION(NZ, NT) :: CTU 		! Upwind: 	Tangential force coefficient
	REAL, DIMENSION(NZ, NT) :: CTD 		! Downwind: Tangential force coefficient
	REAL, DIMENSION(NZ, NT) :: QU 		! Upwind: 	Local torque
	REAL, DIMENSION(NZ, NT) :: QD 		! Downwind: Local torque

	! Coefficients file name:
	CHARACTER(LEN=33) :: filename 				! File name with .csv suffix

	! User's input: [Geometrical]
	REAL :: r 								    ! Turbine's equatorial radius
	REAL :: h 									! Turbine's half height
	REAL :: z0									! Turbine's lowest point
	REAL :: chord_eq 							! Turbine's equatorial chord
	REAL :: chord_tip 							! Turbine's tip chord
	INTEGER :: nb								! Turbine's number of blades
	INTEGER :: naca 							! Turbine's blade section
	INTEGER :: geometry 						! Turbine's shape

	! User's input: [Operational]
	INTEGER :: tsr_i 							! Initial tip speed ratio
	INTEGER :: tsr_f 							! Final tip speed ratio
	REAL :: rpm 								! Revolutions per minute
	REAL :: alpha 								! Wind shear exponent

	! Computed based on user's input:
	REAL :: re_global 							! Global Reynolds number
	REAL :: aspect_ratio 						! Turbine's aspect ratio
	REAL :: swept_area 							! Turbine's swept area
	REAL :: z_eq								! Turbine's equatorial height
	REAL :: dz 									! Height step size
	REAL :: dt 									! Sector step size

	! Inside procedure:
	REAL :: z 									! Local height
	REAL :: zeta 								! Dimensionless height
	REAL :: eta 								! Dimensionless radius
	REAL :: chord 								! Local chord
	REAL :: delta 								! Blade's angle
	REAL :: theta 								! Local sector angle
	REAL :: a0									! Induction factor
	REAL :: a1 									! New induction factor
	REAL :: tsr 								! Local tip speed ratio
	REAL :: free_stream_ratio 					! Free stream ratio
	INTEGER :: iter								! Iterations

	! Indices:
	INTEGER :: i 								! Height index
	INTEGER :: j 								! Sector index
	INTEGER :: k 								! Tip speed ratio index

	! Torque and power coefficients:
	REAL :: cq1 								! Upwind torque coefficient
	REAL :: cq2 								! Downwind torque coefficient
	REAL :: cp1 								! Upwind power coefficient
	REAL :: cp2 								! Downwind power coefficient
	REAL :: cp 									! Total power coefficient


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 2.- Read user's input
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	WRITE (*, *)
	WRITE (*, *)
	WRITE (*, *) "* * * * Double Multiple Stream-Tube Model * * * *"
	WRITE (*, *) "* * * *                                   * * * *"
	WRITE (*, *) "* * * * * * * * * * * * * * * * * * * * * * * * *"
	WRITE (*, *)

	! Prompt for geometrical parameters:
	WRITE (*, '(A)', ADVANCE="NO") "* * * Shape: (Parabola[0], Straight[1])......... "
	READ (*, '(I2)') geometry

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's radius in [m]................... "
	READ (*, '(F5.2)') r

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's half length in [m].............. "
	READ (*, '(F5.2)') h

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's lowest point in [m]............. "
	READ (*, '(F5.2)') z0

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's equatorial chord in [m]......... "
	READ (*, '(F5.2)') chord_eq

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's tip chord in [m]................ "
	READ (*, '(F5.2)') chord_tip

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's number of blades................ "
	READ (*, '(I2)') nb

	WRITE (*, '(A)', ADVANCE="NO") "* * * NACA 00 (12, 15, 18)...................... "
	READ (*, '(I2)') naca

	WRITE (*, *)

	! Prompt for operational parameters:
	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's initial tip-speed ratio......... "
	READ (*, '(I5)') tsr_i

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's final tip-speed ratio........... "
	READ (*, '(I5)') tsr_f

	WRITE (*, '(A)', ADVANCE="NO") "* * * Turbine's rpm............................. "
	READ (*, '(F5.2)') rpm

	WRITE (*, '(A)', ADVANCE="NO") "* * * Wind shear exponent alpha................. "
	READ (*, '(F5.2)') alpha


	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 3.- Initialization of variables
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	! Global Reynolds number:
	re_global = (2.0 * PI / 60.0) * ( (REAL(rpm) * REAL(r)) * REAL(chord_eq) ) / NU

	! Turbine's aspect ratio:
	aspect_ratio = REAL(r) / REAL(h)

	! Equatorial height:
	z_eq = REAL(z0) + REAL(h)

	! Height step size:
	dz = 2.0 / NZ

	! Sector step size:
	dt = PI / NT

	! Swept area:
	IF (geometry .EQ. 0) THEN

		! Parabola:
		swept_area = (2.0 / 3.0) * (4.0 * REAL(h) * REAL(r))

	ELSE IF (geometry .EQ. 1) THEN

		! Straight:
		swept_area = 4.0 * REAL(h) * REAL(r)

	END IF

	! Angle of attack and Reynolds tables:
	filename = "aa.dat"
	CALL parse_dat (filename, AA_TABLE, ROWS)

	filename = "re.dat"
	CALL parse_dat (filename, RE_TABLE, COLS)

	! Lift and drag coefficients:
	SELECT CASE (naca)
	CASE (12)
		filename = "naca0012cl.csv"
		CALL parse_csv (filename, CL_TABLE, ROWS, COLS)
		filename = "naca0012cd.csv"
		CALL parse_csv (filename, CD_TABLE, ROWS, COLS)
	CASE (15)
		filename = "naca0015cl.csv"
		CALL parse_csv (filename, CL_TABLE, ROWS, COLS)
		filename = "naca0015cd.csv"
		CALL parse_csv (filename, CD_TABLE, ROWS, COLS)
	CASE (18)
		filename = "naca0018cl.csv"
		CALL parse_csv (filename, CL_TABLE, ROWS, COLS)
		filename = "naca0018cd.csv"
		CALL parse_csv (filename, CD_TABLE, ROWS, COLS)
	CASE DEFAULT
		WRITE (*, *) "Wrong choice!"
	END SELECT

	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	!
	! 4.- Tip speed ratio loop
	!
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

	! Print the power coefficients headers:
	WRITE (*, *)
	WRITE (*, *)

	WRITE (*, 100) "TSR", "CP1", "CP2", "CPT" 
	100 FORMAT (T6, A6, 6X, A6, 6X, A6, 7X, A6, 3X)

	WRITE (*, 110) "***", "***", "***", "***" 
	110 FORMAT (T6, A6, 6X, A6, 6X, A6, 7X, A6, 3X)
	WRITE (*, *)

	! Tip speed ratio index is "k":
	power: DO k = tsr_i, tsr_f

		! Height sector index is "i":
		vertical: DO i = 1, NZ

			! Dimensionless height zeta [-0.9, 0.9]:
			zeta = -0.9 + (1.8 * (i - 1) / (NZ - 1))

			! Local height in [m]:
			z = (z0 + 0.1 * h) + (1.8 * (i - 1) * h / (NZ - 1))

			! Free stream ratio:
			free_stream_ratio = (z / z_eq) ** (alpha)

			! eta and delta:
			IF (geometry .EQ. 0) THEN

				! Parabola:
				eta = 1.0 - zeta ** 2
				delta = ATAN(2.0 * aspect_ratio * zeta)

			ELSE IF (geometry .EQ. 1) THEN

				! Straight:
				eta = 1.0
				delta = 0.0

			END IF

			! Chord of the airfoil:
			IF (zeta .GT. 0.0) THEN

				chord = zeta * (chord_tip - chord_eq) + chord_eq

			ELSE IF (zeta .LE. 0.0) THEN

				chord = zeta * (chord_eq - chord_tip) + chord_eq

			END IF

			! Angular sector index is "j":

			! * * *
			!	Upstream:
			upstream: DO j = 1, NT

				! theta is the current angular position.
				!	Starts at -PI/2 + DT/2.
				!	Ends at PI/2 - DT/2.
				theta = ((-PI / 2.0) + (dt / 2.0)) + ((j - 1) * (PI - dt) / (NT - 1))

				! Initialize induction factor:
				a0 = 0.0

				! Initialize number of iterations:
				iter = 0

				! Iterative loop:
				iterative1: DO

					! Complement of the induction factor:
					UU(i, j) = 1.0 - a0

					! Local tip speed ratio:
					tsr = (k * eta) / (free_stream_ratio * UU(i, j))

					! Local relative velocity:
					WU(i, j) = SQRT( ((tsr - SIN(theta)) ** 2) + ((COS(theta) ** 2) * (COS(delta) ** 2)) )

					! Local angle of attack:
					AAU(i, j) = ASIN( COS(theta) * COS(delta) / WU(i, j) )

					! Local Reynolds number:
					REU(i, j) = re_global * eta * WU(i, j) / tsr

					! Lift and drag coefficients:
					CALL interpolate_table (AAU(i, j), REU(i, j), CL_TABLE, AA_TABLE, RE_TABLE, ROWS, COLS, CLU(i, j))
					CALL interpolate_table (AAU(i, j), REU(i, j), CD_TABLE, AA_TABLE, RE_TABLE, ROWS, COLS, CDU(i, j))

					! Normal and tangential coefficients:
					CNU(i, j) = CLU(i, j) * COS(AAU(i, j)) + CDU(i, j) * SIN(AAU(i, j))
					CTU(i, j) = CLU(i, j) * SIN(AAU(i, j)) - CDU(i, j) * COS(AAU(i, j))

					! New induction factor:
					a1 = 0.0
					a1 = ( nb * chord / (8.0 * PI * r * eta) )
					a1 = a1 * ( ( CNU(i, j) * COS(theta) / ABS(COS(theta)) ) + ( CTU(i, j) * SIN(theta) / (ABS(COS(theta)) * COS(delta)) ) )
					a1 = a1 * ( WU(i, j) ** 2 )
					a1 = a1 * ( UU(i, j) ** 2 )

					! Glauert factor correction:
					IF (a0 .LE. 0.33) THEN

						a1 = a1 + (a0 ** 2)

					ELSE IF (a0 .GT. 0.33) THEN

						a1 = a1 + ( (1.0 / 4.0) * (5.0 - 3.0 * a0) * (a0 ** 2) )

					END IF

					! Increase iteration counter:
					iter = iter + 1

					! Break out of loop if number of iterations is large:
					IF (iter .GE. 50) THEN

						! Shift the sign:
						UU(i, j) = -1.0 * UU(i, j)
						a1 = 1.0 - UU(i, j)
						EXIT

					END IF

					! Break out of loop if induction factor is above 1.0:
					IF (a1 .GE. 1.0) THEN

						UU(i, j) = 0.0
						a1 = 1.0 - UU(i, j)
						EXIT

					END IF

					! Break out of loop if convergence is achieved:
					IF ( (ABS(a1 - a0) / a0) .LE. CR ) THEN
						EXIT
					END IF

					! If the loop did not break, update the induction factor:
					a0 = a1

				END DO iterative1

				! Collect the torque:
				QU(i, j) = (nb * chord * h) / (2.0 * PI * swept_area)
				QU(i, j) = QU(i, j) * CTU(i, j)
				QU(i, j) = QU(i, j) * ( WU(i, j) * UU(i, j) * free_stream_ratio ) ** 2
				QU(i, j) = QU(i, j) * (eta / COS(delta))

			END DO upstream

			! * * *
			!	Downstream:
			downstream: DO j = 1, NT

				! theta is the current angular position.
				!	Starts at PI/2 + DT/2.
				!	Ends at 3 PI/2 - DT/2.
				theta = ((PI / 2.0) + (dt / 2.0)) + ((j - 1) * (PI - dt) / (NT - 1))

				! Initialize induction factor based on the upstream tube:
				a0 = 1.0 - UU(i, NT - j + 1)

				! Initialize number of iterations:
				iter = 0

				! Iterative loop:
				iterative2: DO

					! Complement of the induction factor:
					UD(i, j) = 1.0 - a0

					! Local tip speed ratio:
					tsr = (k * eta) / ( free_stream_ratio * UD(i, j) * (2.0 * UU(i, NT - j + 1) - 1.0) )

					! Local relative velocity:
					WD(i, j) = SQRT( ((tsr - SIN(theta)) ** 2) + ((COS(theta) ** 2) * (COS(delta) ** 2)) )

					! Local angle of attack:
					AAD(i, j) = ASIN( COS(theta) * COS(delta) / WD(i, j) )

					! Local Reynolds number:
					RED(i, j) = re_global * eta * WD(i, j) / tsr

					! Lift and drag coefficients:
					CALL interpolate_table (AAD(i, j), RED(i, j), CL_TABLE, AA_TABLE, RE_TABLE, ROWS, COLS, CLD(i, j))
					CALL interpolate_table (AAD(i, j), RED(i, j), CD_TABLE, AA_TABLE, RE_TABLE, ROWS, COLS, CDD(i, j))

					! Normal and tangential coefficients:
					CND(i, j) = CLD(i, j) * COS(AAD(i, j)) + CDD(i, j) * SIN(AAD(i, j))
					CTD(i, j) = CLD(i, j) * SIN(AAD(i, j)) - CDD(i, j) * COS(AAD(i, j))

					! New induction factor:
					a1 = 0.0
					a1 = ( nb * chord / (8.0 * PI * r * eta) )
					a1 = a1 * ( ( CND(i, j) * COS(theta) / ABS(COS(theta)) ) + ( CTD(i, j) * SIN(theta) / (ABS(COS(theta)) * COS(delta)) ) )
					a1 = a1 * ( WD(i, j) ** 2 )
					a1 = a1 * ( UD(i, j) ** 2 )

					! Glauert factor correction:
					IF (a0 .LE. 0.33) THEN

						a1 = a1 + (a0 ** 2)

					ELSE IF (a0 .GT. 0.33) THEN

						a1 = a1 + ( (1.0 / 4.0) * (5.0 - 3.0 * a0) * (a0 ** 2) )

					END IF

					! Increase iteration counter:
					iter = iter + 1

					! Break out of loop if number of iterations is large:
					IF (iter .GE. 50) THEN

						! Shift the sign:
						UD(i, j) = -1.0 * UD(i, j)
						a1 = 1.0 - UD(i, j)
						EXIT

					END IF

					! Break out of loop if induction factor is above 1.0:
					IF (a1 .GE. 1.0) THEN

						UD(i, j) = 0.0
						a1 = 1.0 - UD(i, j)
						EXIT

					END IF

					! Break out of loop if convergence is achieved:
					IF ( (ABS(a1 - a0) / a0) .LE. CR ) THEN
						EXIT
					END IF

					! If the loop did not break, update the induction factor:
					a0 = a1

				END DO iterative2

				! Collect the torque:
				QD(i, j) = (nb * chord * h) / (2.0 * PI * swept_area)
				QD(i, j) = QD(i, j) * CTD(i, j)
				QD(i, j) = QD(i, j) * ( WD(i, j) * UD(i, j) * (2.0 * UU(i, NT - j + 1) - 1.0) * free_stream_ratio ) ** 2
				QD(i, j) = QD(i, j) * (eta / COS(delta))

			END DO downstream

		END DO vertical

		! Compute torque and power coefficients at the current tip speed ratio:
		CALL integrate_torque (QU, NZ, NT, dz, dt, cq1)
		CALL integrate_torque (QD, NZ, NT, dz, dt, cq2)

		cp1 = k * cq1
		cp2 = k * cq2

		! The total power coefficient is the sum of the upwind and downwind power coefficients:
		cp = cp1 + cp2

		! Print the power coefficients headers:
		WRITE (*, 200) REAL(k), cp1, cp2, cp 
		200 FORMAT (T6, F6.2, 6X, F6.2, 6X, F6.2, 7X, F6.2, 3X)

	END DO power

	WRITE (*, *)
	WRITE (*, *)


END PROGRAM main
