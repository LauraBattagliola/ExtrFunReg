      SUBROUTINE metrl2(A,B,M,N,D,METR)
c     computes a fast approximation of the L2 metric between A and B
c     supporting functions on regular grids only
c     used in the h-mode depth 
        
      INTEGER N,M,D
      double precision A(M*D),B(N*D),METR(M*N)
      INTEGER I,J,K

      DO 15 I=1,M
            DO 10 J=1,N
            METR((J-1)*M+I) = 0.0
                  DO 5 K=1,D
                  METR((J-1)*M+I) = METR((J-1)*M+I) + (A((K-1)*M+I)-
     +B((K-1)*N+J))**2
5                 CONTINUE
            METR((J-1)*M+I) = sqrt(METR((J-1)*M+I) - 
     +((A((0)*M+I)-B((0)*N+J))**2+(A((D-1)*M+I)-B((D-1)*N+J))**2)
     +/(2.0))
10          CONTINUE
15    CONTINUE
      RETURN
      END 
      
      SUBROUTINE metrl2B(A,B,M,N,D,METR)
c     computes a fast approximation of the L2 metric between A and B
c     supporting functions on regular grids only
c     non-corrected for boundary effects version
        
      INTEGER N,M,D
      double precision A(M*D),B(N*D),METR(M*N)
      INTEGER I,J,K

      DO 15 I=1,M
            DO 10 J=1,N
            METR((J-1)*M+I) = 0.0
                  DO 5 K=1,D
                  METR((J-1)*M+I) = METR((J-1)*M+I) + (A((K-1)*M+I)-
     +B((K-1)*N+J))**2
5                 CONTINUE
            METR((J-1)*M+I) = sqrt(METR((J-1)*M+I))
10          CONTINUE
15    CONTINUE
      RETURN
      END 