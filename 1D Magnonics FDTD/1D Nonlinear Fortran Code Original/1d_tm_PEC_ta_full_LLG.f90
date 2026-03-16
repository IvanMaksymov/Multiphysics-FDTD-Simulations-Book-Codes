	PROGRAM FDTD !big
    
    IMPLICIT NONE
    
    INTEGER MAX_ITER, NY
    INTEGER STIMULUS_SIZE
    !PARAMETER (MAX_ITER =  2D7)
    PARAMETER (NY = 76)
    REAL*8 PI
    REAL*8 LIGHT_SPEED
    REAL*8 MU_0
    REAL*8 EPSILON_0
    INTEGER J, JB, J1, J2, K, I, IH0, J_SOURCE
    INTEGER ITERATION
    REAL*8 STIMULUS
    REAL*8 CURRENTSIMULATEDTIME
    REAL*8 TOTALSIMULATEDTIME
    REAL*8 OMEGA
    REAL*8 LAMBDA
    REAL*8 DY
    REAL*8 DT
    REAL*8 TW, T0
    REAL*8 DTDY
    
    REAL*8 EPSILON(1:NY+1), EPS_INF, SGM(1:NY+1)
    REAL*8 EZ(1:NY+1), JZ(1:NY+1)
	REAL*8 HX(1:NY), HY(1:NY-1) !step n+1/2
	REAL*8 HXN(1:NY), HYN(1:NY-1) !step n
	REAL*8 oldHXN(1:NY) !step n-1/2		
	REAL*8 BX(1:NY) !step n+1/2
	REAL*8 oldBX(1:NY) !step n-1/2	
	REAL*8 BXN(1:NY) !step n
	
	REAL*8 MX(1:NY+1), MY(1:NY+1), MZ(1:NY+1) !magnetization step n+1/2
	REAL*8 MXN(1:NY+1), MYN(1:NY+1), MZN(1:NY+1) !magnetization step n
	REAL*8 HEX_X(1:NY+1), HEX_Y(1:NY+1), HEX_Z(1:NY+1) !excahnge field
	REAL*8 oldMXN(1:NY+1), oldMYN(1:NY+1), oldMZN(1:NY+1) !magnetization step n, previous iteration

	REAL*8 DFT_EZ1(1:NY+1), DFT_JZ1(1:NY+1)
	REAL*8 DFT_HX1(1:NY)
	REAL*8 DFT_MX1(1:NY+1), DFT_MY1(1:NY+1) 
	REAL*8 DFT_HEX_X1(1:NY+1), DFT_HEX_Y1(1:NY+1)

	REAL*8 DFT_EZ2(1:NY+1), DFT_JZ2(1:NY+1)
	REAL*8 DFT_HX2(1:NY)
	REAL*8 DFT_MX2(1:NY+1), DFT_MY2(1:NY+1) 
	REAL*8 DFT_HEX_X2(1:NY+1), DFT_HEX_Y2(1:NY+1)
	
	REAL*8 SUM_MX(1:NY+1), SUM_MY(1:NY+1), SUM_HX(1:NY)
		    
    INTEGER PLOT_MODULUS, DUMMY, SAMPLING
    REAL*8 WORK_FREQ, FT_FREQ(3)
    REAL*8 R, W_CO, W_PT, D, MS1, MS2, A1, A2, A12, H0, ALPHA1, ALPHA2, GAMMA1, GAMMA2, SIGMA1, SIGMA2
    REAL*8 MS, A, ALPHA, GAMMA
    REAL*8 OE, GAUSS, ERG_over_CM, TESLA, ERG_over_CM2
    REAL*8 AZ, X0_X, X0_Y, ACC
	
    INTEGER CENTER_Y
    
    REAL*8 C1,C2,TMP1,TMP2,TMP,WINDOW,A_VEC(3),X_VEC(3),A_TIMES_X_VEC(3),M0_VEC(3),H0_VEC(3),HEFF(3),SUM_VEC(3),TMP_VEC(3)
    REAL*8 M_VEC(3), A_M_VEC(3), A_M_DOT
    REAL*8 VEC_ABS, CDOT, SCALING, U1, U2, U3, U4, U5, U6, U7
    INTEGER FLAG, FLAG_2_FILMS, NH0, ITMP
    REAL*8 T1, T2, TA_LAYER
    
! CPML PARAMETERS
	REAL*8 M
! POLINOMIAL GRADING
	REAL*8 D2Y
! PML THICKNESS IN METERS
	REAL*8 SIGMA_MAX
! MAX VALUE OF SIGMA PARAMETER
	REAL*8 A_MAX
! MAX VALUE OF A PARAMETER
	REAL*8 KMAX
! MAX VALUE OF K PARAMETER
    
    OPEN (UNIT=2,FILE='CPML.OUT',STATUS='UNKNOWN')
    OPEN (UNIT=53,FILE='DETECTOR4.OUT',STATUS='UNKNOWN')
    OPEN (UNIT=54,FILE='DETECTOR5.OUT',STATUS='UNKNOWN')
    
    PLOT_MODULUS = 1D7
    SAMPLING = 1	
    EPS_INF = 1.0D0
    CENTER_Y = NY/2
    
    ACC = 1.0D-6
    
    SCALING = 1.0D0
    FT_FREQ(1) = 18D9
    FT_FREQ(2) = FT_FREQ(1)
    FT_FREQ(3) = FT_FREQ(1)
        
    PI = 3.1415926D0
    LIGHT_SPEED = 299792458.0D0
    MU_0 = 1.256637D-6
    EPSILON_0 =  8.8541878D-12
    
    OE = 1.0D3/4.0D0/PI
	GAUSS = 1.0D3
	ERG_over_CM = 1.0D-7/1.0D-2
	ERG_over_CM2 = 1.0D-7/1.0D-4
	TESLA = 795.744D3

    W_CO = 10.0D-9
    W_PT = 90.0D-9
    
    FLAG_2_FILMS = 1
    
    NH0 = 100
    
    DO IH0 = 1, NH0+1
    
    call cpu_time(t1)
    
    H0 = 4000.0D0*OE*(1.0D0 - DBLE(IH0-1)/DBLE(NH0))
    !H0 = 450.0D0*OE

    WRITE(*,*) 'H0: (A/m, Oe) = ', H0, H0/OE      
	    
    DO J = 1,NY+1
    
        EPSILON(J) = EPS_INF
        HEX_X(J) = 0.0d0
        HEX_Y(J) = 0.0d0        
        SGM(J) = 0.0d0
        EZ(J) = 0.0d0
        JZ(J) = 0.0d0		
		MX(J) = 0.0d0
		MY(J) = 0.0d0
		MZ(J) = 0.0d0
		MXN(J) = 0.0d0
		MYN(J) = 0.0d0
		MZN(J) = 0.0d0
		HEX_X(J) = 0.0d0
		HEX_Y(J) = 0.0d0
		HEX_Z(J) = 0.0d0
		oldMXN(J) = 0.0d0
		oldMYN(J) = 0.0d0
		oldMZN(J) = 0.0d0
		SUM_MX(J) = 0.0d0
		SUM_MY(J) = 0.0d0
    
    ENDDO

    DO J = 1,NY    
    
		HX(J) = 0.0d0
		HXN(J) = 0.0d0
		BX(J) = 0.0d0
		oldBX(J) = 0.0d0
		oldHXN(J) = 0.0d0
		BXN(J) = 0.0d0
				
		DFT_HX1(J) = 0.0D0
		DFT_HX2(J) = 0.0D0
					
		DFT_MX1(J) = 0.0D0
		DFT_MX2(J) = 0.0D0

		DFT_MY1(J) = 0.0D0
		DFT_MY2(J) = 0.0D0

		DFT_EZ1(J) = 0.0D0
		DFT_EZ2(J) = 0.0D0

		DFT_JZ1(J) = 0.0D0
		DFT_JZ2(J) = 0.0D0
			
		DFT_HEX_X1(J) = 0.0D0
		DFT_HEX_X2(J) = 0.0D0

		DFT_HEX_Y1(J) = 0.0D0
		DFT_HEX_Y2(J) = 0.0D0
					
	ENDDO
	
	DO J = 1,NY - 1 
	
		HYN(J) = 0.0D0
		HY(J) = 0.0D0
	
	ENDDO

    DY = 1.6667D-9
    DT = 1.0D0*DY*DSQRT(EPS_INF)/LIGHT_SPEED
    PRINT 5, DY, DT
    5 FORMAT (1X,' DY=',E10.4,' DT=',E10.4)
    
    MAX_ITER = INT(50.0D0/FT_FREQ(3)/DT)
        
    IF(FLAG_2_FILMS.EQ.1) THEN
    
 		A12 = 2D-6*ERG_over_CM2/DY
		D = 0.0D0 !pinning parameter
		MS1 = 15080.0D0*GAUSS/4.0D0/PI
		A1 = 1.0D-6*ERG_over_CM                     
		MS2 = 8042.0D0*GAUSS/4.0D0/PI
		A2 = 0.55D-6*ERG_over_CM
		ALPHA1 = 0.016
		ALPHA2 = 0.008
		SIGMA1 = 1.8D7
		SIGMA2 = 4.5D6
		GAMMA1 = 2.0D0*PI*2.92D6/OE
		GAMMA2 = GAMMA1

	ELSE
		
		D = 0.0D0 !pinning parameter
		MS1 = 8042.0D0*GAUSS/4.0D0/PI
		A1 = 0.55D-6*ERG_over_CM
		MS2 = MS1
		A2 = A1
		A12 = A1/DY
		ALPHA1 = 0.008
		ALPHA2 = ALPHA1
		SIGMA1 = 4.5D6
		SIGMA2 = SIGMA1
		GAMMA1 = 2.0D0*PI*2.92D6/OE 
		GAMMA2 = GAMMA1

	ENDIF
	
	WRITE(*,*) 'Ms1: (A/m, Gauss) = ', MS1, MS1*4.0D0*PI/GAUSS
	WRITE(*,*) 'Ms2: (A/m, Gauss) = ', MS2, MS2*4.0D0*PI/GAUSS
	WRITE(*,*) 'A1: (J/m, erg/cm) = ', A1, A1/ERG_over_CM	
	WRITE(*,*) 'A2: (J/m, erg/cm) = ', A2, A2/ERG_over_CM
	WRITE(*,*) 'A12: (J/m^2, erg/cm^2) = ', A12, A12/ERG_over_CM2
	WRITE(*,*) 'Sigma1: (S/m) = ', SIGMA1	
	WRITE(*,*) 'Sigma2: (S/m) = ', SIGMA2
	WRITE(*,*) 'Alpha1:  = ', ALPHA1	
	WRITE(*,*) 'Alpha2:  = ', ALPHA2
	WRITE(*,*) 'Gamma1: (Hz m/A, MHz/Oe/2pi) = ', GAMMA1, GAMMA1*OE/1D6/2.0D0/PI
	WRITE(*,*) 'Gamma2: (Hz m/A, MHz/Oe/2pi) = ', GAMMA2, GAMMA2*OE/1D6/2.0D0/PI
	    		
	CALL CONFIG(NY,EPSILON,SGM,0,DY,W_CO,W_PT,JB,J1,J2,SIGMA1,SIGMA2,EPS_INF,TA_LAYER)
	
	WRITE(*,*) 'Ta LAYER ', TA_LAYER
	WRITE(*,*) 'Co LAYER', W_CO
	WRITE(*,*) 'Py LAYER', W_PT

	!INITIAL CONDITION FOR MAGNETIZATION
	!MX = 0, MY = 0, MZ = MS
	DO J=J1, JB
		
		MX(J) = 0.0D0
		MY(J) = 0.0D0
		MZ(J) = MS1				
				
	ENDDO
			
	DO J=JB+1, J2

		MX(J) = 0.0D0
		MY(J) = 0.0D0			
		MZ(J) = MS2
					
	ENDDO	
    
    DO ITERATION = 0, MAX_ITER
    
		CURRENTSIMULATEDTIME = DT*DBLE(ITERATION)
    
		IF (MOD(ITERATION, PLOT_MODULUS).EQ.0) THEN
    
			PRINT*, 'STEP NUMBER ', ITERATION, ' FROM', MAX_ITER, ' TIME, SECONDS ', CURRENTSIMULATEDTIME		
			
			call cpu_time(t2)
			print *,'CPU time in seconds ',t2-t1
			
		ENDIF
    
		! COMPUTE THE STIMULUS:
    
		! ANGULAR FREQUENCY IN RADIANS/SECOND:
    	WORK_FREQ = FT_FREQ(3)
    	OMEGA = 2.0D0*PI*WORK_FREQ
 	
		STIMULUS = DCOS(OMEGA*CURRENTSIMULATEDTIME)

    
! UPDATE THE BX VALUES -- STEP 1

		DO J=1, NY-1

			DTDY=DT/DY
			oldBX(J) = BX(J)! step n-1/2
			BX(J) = BX(J) + ( - DTDY*(EZ(J+1) - EZ(J))) !step n+1/2
			BXN(J) = 0.5D0*(BX(J) + oldBX(J)) !step n
			oldHXN(J) = HXN(J)
			
		ENDDO

		! STEP 2. No spatial interpolation of BXN for 1D Model because bx(n-1/2) is already at the same point as m

		! STEP 3
		K = 0
		DO !iterations
			
			DO J=J1, J2
			
				IF(J.LE.JB) THEN
				
					MS = MS1
					A = A1
					GAMMA = GAMMA1
					ALPHA = ALPHA1
					
				ELSE
				
					MS = MS2
					A = A2
					GAMMA = GAMMA2
					ALPHA = ALPHA2					

				ENDIF									
				
				HEFF(1) = HXN(J) + HEX_X(J)
				HEFF(2) = HYN(J) + HEX_Y(J)
				HEFF(3) = H0 + HEX_Z(J)
				
				A_VEC(1) = -(0.5D0*DABS(GAMMA)*DT*HEFF(1) + ALPHA*MX(J)/MS)
				A_VEC(2) = -(0.5D0*DABS(GAMMA)*DT*HEFF(2) + ALPHA*MY(J)/MS)
				A_VEC(3) = -(0.5D0*DABS(GAMMA)*DT*HEFF(3) + ALPHA*MZ(J)/MS)
												
				oldMXN(J) = MXN(J) !previous iteration
				oldMYN(J) = MYN(J) !previous iteration
				oldMZN(J) = MZN(J) !previous iteration
				
				M_VEC(1) = MX(J)
				M_VEC(2) = MY(J)
				M_VEC(3) = MZ(J)
				
				CALL CTIMES(A_VEC, M_VEC, A_M_VEC)
				
				MXN(J) = (M_VEC(1) + CDOT(A_VEC, M_VEC)*A_VEC(1) - A_M_VEC(1))/(1.0D0 + VEC_ABS(A_VEC)**2)
				MYN(J) = (M_VEC(2) + CDOT(A_VEC, M_VEC)*A_VEC(2) - A_M_VEC(2))/(1.0D0 + VEC_ABS(A_VEC)**2)
				MZN(J) = (M_VEC(3) + CDOT(A_VEC, M_VEC)*A_VEC(3) - A_M_VEC(3))/(1.0D0 + VEC_ABS(A_VEC)**2)
							
			ENDDO
			
			!BOUNDARY CONDITION FOR MXN and MYN
			
			!Outer surface of film 1
			MYN(J1) = MYN(J1+1)/(1.0D0 - D*DY)
			MXN(J1) = MXN(J1+1)
			MZN(J1) = MZN(J1+1)/(1.0D0 - D*DY) !?

			!Outer surface of film 2
			MYN(J2) = MYN(J2-1)/(1.0D0 - D*DY)
			MXN(J2) = MXN(J2-1)
			MZN(J2) = MZN(J2-1)/(1.0D0 - D*DY) !?
			
			!Inner interface
			MXN(JB) = (MXN(JB-1) + DY*(A12*MS1/(A1*MS2))*MXN(JB+1))/(1.0D0 + DY*A12/A1)
			MYN(JB) = (MYN(JB-1) + DY*(A12*MS1/(A1*MS2))*MYN(JB+1))/(1.0D0 + DY*A12/A1)
			MZN(JB) = (MZN(JB-1) + DY*(A12*MS1/(A1*MS2))*MZN(JB+1))/(1.0D0 + DY*A12/A1)
			 
			MXN(JB+1) = (MXN(JB+2) + DY*(A12*MS2/(A2*MS1))*MXN(JB))/(1.0D0 + DY*A12/A2)
			MYN(JB+1) = (MYN(JB+2) + DY*(A12*MS2/(A2*MS1))*MYN(JB))/(1.0D0 + DY*A12/A2)
			MZN(JB+1) = (MZN(JB+2) + DY*(A12*MS2/(A2*MS1))*MZN(JB))/(1.0D0 + DY*A12/A2)
						
			!EXCHANGE FIELD FILM 1
			DO J=J1+1, JB
			
				MS = MS1
				A = A1
									
				HEX_X(J) = (2.0D0*A/(MU_0*MS**2))*(MXN(J+1)-2.0D0*MXN(J)+MXN(J-1))/DY**2
				HEX_Y(J) = (2.0D0*A/(MU_0*MS**2))*(MYN(J+1)-2.0D0*MYN(J)+MYN(J-1))/DY**2
				HEX_Z(J) = (2.0D0*A/(MU_0*MS**2))*(MZN(J+1)-2.0D0*MZN(J)+MZN(J-1))/DY**2   				
				
			ENDDO
			
			! EXTRAPOLATION TO FIND HEX IN J1
			HEX_X(J1) = 2.0D0*HEX_X(J1+1)-HEX_X(J1+2)
			HEX_Y(J1) = 2.0D0*HEX_Y(J1+1)-HEX_Y(J1+2)
			HEX_Z(J1) = 2.0D0*HEX_Z(J1+1)-HEX_Z(J1+2)
							
			!EXCHANGE FIELD FILM 2
			DO J=JB+1, J2-1
			
				MS = MS2
				A = A2
				
				HEX_X(J) = (2.0D0*A/(MU_0*MS**2))*(MXN(J+1)-2.0D0*MXN(J)+MXN(J-1))/DY**2
				HEX_Y(J) = (2.0D0*A/(MU_0*MS**2))*(MYN(J+1)-2.0D0*MYN(J)+MYN(J-1))/DY**2
				HEX_Z(J) = (2.0D0*A/(MU_0*MS**2))*(MZN(J+1)-2.0D0*MZN(J)+MZN(J-1))/DY**2 								
				
			ENDDO
			
			! EXTRAPOLATION TO FIND HEX IN J2
			HEX_X(J2) = 2.0D0*HEX_X(J2-1)-HEX_X(J2-2)
			HEX_Y(J2) = 2.0D0*HEX_Y(J2-1)-HEX_Y(J2-2)
			HEX_Z(J2) = 2.0D0*HEX_Z(J2-1)-HEX_Z(J2-2)			
							
			!MAGNETIC FIELD X step n
			DO J=1, NY			
				
				HXN(J) = BXN(J)/MU_0 - MXN(J) !step n, current iteration
				
			ENDDO
			
			!MAGNETIC FIELD Y step n
 
			DO J=1, NY-1

				!!!-------->HYN(J) = - MYN(J) WAS WRONG!!!
				HYN(J) =  - 0.5D0*(MYN(J+1)+MYN(J))

			ENDDO
						
			TMP1 = 0.0D0
			TMP2 = 0.0D0
		
			DO J=J1, J2
			
				TMP1 = TMP1 + DSQRT(MXN(J)**2+MYN(J)**2+MZN(J)**2)
				TMP2 = TMP2 + DSQRT(oldMXN(J)**2+oldMYN(J)**2+oldMZN(J)**2)	
				
			ENDDO						
						
			IF(TMP1/TMP2-1.0D0.GT.ACC) THEN

				FLAG = 1
					
			ELSE
				
				FLAG = 0
										
			ENDIF			
			
			K = K + 1
			
			IF(FLAG.EQ.0) EXIT
						
		ENDDO!WHILE
		
		IF(K.GT.100) THEN
			WRITE(*,*) '#ITERS DONE = ', K
			WRITE(*,*) 'TOO MANY ITERATIONS! SOMETHING IS WRONG! STOP'
			STOP 				
		ENDIF					

		! Step 5
33		DO J=J1, J2

			MX(J) = 2.0D0*MXN(J) - MX(J) !step n+1/2
			MY(J) = 2.0D0*MYN(J) - MY(J) !step n+1/2
			MZ(J) = 2.0D0*MZN(J) - MZ(J) !step n+1/2
			
		ENDDO	
				
		!Step 6: interpolation is not required for the 1D model

		DO J=1, NY-1

! UPDATE THE HX VALUES
			
			HX(J) = BX(J)/MU_0 - MX(J)			
			
!	UPDATE THE HY VALUES

			HY(J) = -0.5D0*(MY(J+1)+MY(J))

		ENDDO

		! EXCITATION
		J_SOURCE = J1-CEILING(TA_LAYER/DY)-1
		JZ(J_SOURCE) = STIMULUS

    		
! UPDATE THE EZ VALUES

		DO J=2, NY

         	C1=2.0D0*EPSILON(J)*EPSILON_0-SGM(J)*DT
    		C2=2.0D0*EPSILON(J)*EPSILON_0+SGM(J)*DT

			EZ(J) = (C1/C2)*EZ(J)+ (2.0D0*DT/C2)*(- (HX(J) - HX(J-1))/DY - JZ(J))		
	
			!ON_FLY FT	
			K = 1
			DFT_HX1(J) = DFT_HX1(J) + HX(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_HX2(J) = DFT_HX2(J) + HX(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
		
			DFT_MX1(J) = DFT_MX1(J) + MX(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_MX2(J) = DFT_MX2(J) + MX(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING

			DFT_MY1(J) = DFT_MY1(J) + MY(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_MY2(J) = DFT_MY2(J) + MY(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING

			DFT_EZ1(J) = DFT_EZ1(J) + EZ(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_EZ2(J) = DFT_EZ2(J) + EZ(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING

			DFT_JZ1(J) = DFT_JZ1(J) + JZ(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_JZ2(J) = DFT_JZ2(J) + JZ(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			
			DFT_HEX_X1(J) = DFT_HEX_X1(J) + HEX_X(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_HEX_X2(J) = DFT_HEX_X2(J) + HEX_X(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING		

			DFT_HEX_Y1(J) = DFT_HEX_Y1(J) + HEX_Y(J)*DCOS(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING
			DFT_HEX_Y2(J) = DFT_HEX_Y2(J) + HEX_Y(J)*DSIN(2.0D0*PI*FT_FREQ(K)*CURRENTSIMULATEDTIME)*DT*SAMPLING		

			SUM_MX(J) = SUM_MX(J) + MX(J)
			SUM_MY(J) = SUM_MY(J) + MY(J)
			SUM_HX(J) = SUM_HX(J) + HX(J)

		ENDDO
    
    ENDDO!FDTD

	IF(IH0.EQ.0) THEN
  	  	
		CALL PLOT_DFT(NY,DFT_HX1,DFT_MX1,DFT_MY1,DFT_EZ1,DFT_JZ1,DFT_HEX_X1,DFT_HEX_Y1,ITERATION,1)
		CALL PLOT_DFT(NY,DFT_HX2,DFT_MX2,DFT_MY2,DFT_EZ2,DFT_JZ2,DFT_HEX_X2,DFT_HEX_Y2,ITERATION,2)
		CALL PLOT(NY,SUM_HX,SUM_MX,SUM_MY,EZ,JZ,HEX_X,HEX_Y,0)
	
	ENDIF
  	
  	U1 = 0.0D0
  	U2 = 0.0D0
  	U3 = 0.0D0
  	U4 = 0.0D0
  	U5 = 0.0D0
  	U6 = 0.0D0
  	U7 = 0.0D0

	DO J = J1, J2
		U1 = U1 + DFT_MX1(J)
		U2 = U2 + DFT_MX2(J)
		U3 = U3 + DFT_MY1(J)
		U4 = U4 + DFT_MY2(J)
		U5 = U5 + SUM_MX(J)
		U6 = U6 + SUM_MY(J)
	ENDDO
	
	U7 = U7 + EZ(J_SOURCE)
	
	WRITE(54,'(10(1X,E15.8))') H0,U1,U2,U3,U4,U5,U6,DFT_EZ1(J_SOURCE),DFT_EZ2(J_SOURCE),U7

	ENDDO!H0
    
    CLOSE(2)
    CLOSE(53)
    
END
    
    
    

