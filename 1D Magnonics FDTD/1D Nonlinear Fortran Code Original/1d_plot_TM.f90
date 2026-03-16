   	SUBROUTINE PLOT(NY,HX,MX,MY,EZ,JZ,HEX_X,HEX_Y,ITERATION)
	
	IMPLICIT NONE
	
	INTEGER NY
	REAL*8 EZ(1:NY+1), JZ(1:NY+1)
	REAL*8 HX(1:NY) !STEP N+1/2
	REAL*8 MX(1:NY+1), MY(1:NY+1) !MAGNETIZATION STEP N+1/2
	REAL*8 HEX_X(1:NY+1), HEX_Y(1:NY+1) !excahnge field
	CHARACTER*32 FILENAME
	INTEGER I, J, K, CENTER_Y, ITERATION
				
	CENTER_Y = NY/2
	

    WRITE(FILENAME,'(A,I0,A)')'HEX_X_',ITERATION,'.OUT'

	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  HEX_X(J)
        ENDDO
        CLOSE(34)	

    WRITE(FILENAME,'(A,I0,A)')'HEX_Y_',ITERATION,'.OUT'

	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  HEX_Y(J)
        ENDDO
        CLOSE(34)

	WRITE(FILENAME,'(A,I0,A)')'MX_',ITERATION,'.OUT'
		
	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  MX(J)
        ENDDO
        CLOSE(34)


	WRITE(FILENAME,'(A,I0,A)')'MY_',ITERATION,'.OUT'
		
	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  MY(J)
        ENDDO
        CLOSE(34)

		WRITE(FILENAME,'(A,I0,A)')'HX_',ITERATION,'.OUT'
		
	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY
                WRITE(34,"(1X, E18.8E3,$)")  HX(J)
        ENDDO
        CLOSE(34)

		WRITE(FILENAME,'(A,I0,A)')'EZ_',ITERATION,'.OUT'
		
	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  EZ(J)
        ENDDO
        CLOSE(34)

		WRITE(FILENAME,'(A,I0,A)')'JZ_',ITERATION,'.OUT'
		
	OPEN (UNIT=34,FILE=filename,STATUS='UNKNOWN')
        DO J=1, NY+1
                WRITE(34,"(1X, E18.8E3,$)")  JZ(J)
        ENDDO
        CLOSE(34)
        
	END
