C PROGRAM wire
C Pulse basis and point matching solution for a thin wire
C Copyright (c) 2015 by Robert Paknys
C Requires: constants_m.f moler_m.f
      MODULE simpdata_mod !k,a shared data for the integrand
      IMPLICIT NONE
      REAL:: k,a
      END MODULE simpdata_mod
C---
      MODULE simp_mod
      CONTAINS
C--- Function to be integrated
      FUNCTION g(z,zprime)
      USE constants_mod
      USE simpdata_mod
      IMPLICIT NONE
      COMPLEX:: g
      REAL, INTENT(IN):: z,zprime
      REAL:: RR
      RR=sqrt((z-zprime)**2+a**2)
      g= -cj*eta0/k*EXP(-cj*k*RR)/(4*pi*RR**5)*
     ? ((1+cj*k*RR)*(2*RR**2-3*a**2)+(k*a*RR)**2 )
      !print *, "z = ", z, "   ,z_prime = ", zprime, "  , g=",g
      RETURN
      END FUNCTION g
C--- Apply complex Simpson's rule to g(z,z') dz prime
      FUNCTION csimp(a,b,n_slices,z)
      IMPLICIT NONE
      REAL, INTENT(IN):: a,b,z
      INTEGER, INTENT(IN):: n_slices
      COMPLEX:: csimp
      REAL:: h
      COMPLEX:: area
      INTEGER:: n, i
      IF(n_slices .LE. 0)THEN
        WRITE(0,*)' Stopped: csimp, n_slices .LE. 0'
        STOP
      ENDIF
      n=n_slices
      IF(n/2*2 .NE. n)n=n+1
      h=(b-a)/n
      !print *, "hello before the two gs"
      area=g(z,a)+g(z,b)
      !print *, "hello after the two g and befre the odd loop"
      DO i=1, n-1, 2
      area=area+4*g(z,a+i*h)
      !print *,"i= ", n
      ENDDO
      !print *, "hello after the odd loop and before the even loop"
      DO i=2, n-2, 2
      area=area+2*g(z,a+i*h)
      !print *,"n = ", n
      ENDDO
      !print *, "hello after the even loop"

      area=area*h/3
      csimp=area
      RETURN
      END FUNCTION csimp
C--- Write g(z,z') to a file, for plotting
      SUBROUTINE cplot_integrand(a,b,n_slices,outfile,z)
      IMPLICIT NONE
      REAL, INTENT(IN):: a,b,z
      INTEGER, INTENT(IN):: n_slices
      REAL:: h,zp
      COMPLEX:: fx
      INTEGER:: i
      CHARACTER(LEN=*), INTENT(IN):: outfile
      h=(b-a)/n_slices
      OPEN(UNIT=51, FILE=outfile, STATUS='UNKNOWN')
       DO i=0,n_slices
       zp=a+i*h
       WRITE(51,*)zp, REAL(g(z,zp)), IMAG(g(z,zp))
       ENDDO
      CLOSE(51)
      RETURN
      END SUBROUTINE cplot_integrand
C---
      SUBROUTINE getzmn(k0,a0,n_seg,delta,length,zmn)
      USE constants_mod
      USE simpdata_mod !k,a for the integrand
      IMPLICIT NONE
      REAL, INTENT(IN):: k0,a0,delta,length
      INTEGER, INTENT(IN):: n_seg
      COMPLEX, DIMENSION(:,:), INTENT(OUT):: zmn
c
      COMPLEX, DIMENSION(SIZE(zmn(1,:))):: zmn1 !automatic array
      REAL:: zn1,zn2,zn,thetan1,thetan2,zm
      INTEGER:: nsr,nsi,n_slices,m,n,i
      IF(k0*a0 .LT. tpi*0.0005 .OR. delta .GT. 0.25)THEN
      WRITE(0,*)'Stopped: SUBROUTINE getzmn not tested for these values'
      STOP
      ENDIF
      k=k0
      a=a0
      zm=-length/2+delta/2 !match point position
      n_slices= 100*(.005/a)*(delta*80.) !ok if delta=1/80, a=0.005
      !n_slices = 6
      print *,"n_slice = ", n_slices
      loop1: DO n=1,n_seg
      zn=-length/2+delta/2+(n-1)*delta
      zn1=zn-delta/2
      zn2=zn+delta/2
      zmn1(n) = -csimp(zn1,zn2,n_slices,zm) !First row of zmn
      !print *, "from the zmn determination not CSIMP"
      !print *, "n=,",n, "    ,zn1=",zn1,"  ,   zn2=",zn2,"    , zm=",zm
      !print *, "zmn1(n) = ", zmn1(n)
c      IF(n .EQ. 1) CALL cplot_integrand(zn1,zn2,n_slices,'int.dat',zm)
      ENDDO loop1
c
      DO m=1,n_seg
      DO n=1,n_seg
      zmn(m,n)=zmn1(ABS(m-n)+1) !Use first row to fill Toeplitz zmn
      ENDDO
      ENDDO
c
      RETURN
      END SUBROUTINE getzmn
C---
      END MODULE simp_mod
C--- Main Program -----------------
C Thin wire pulse basis point matching solution.
      PROGRAM wire
      USE constants_mod
      USE simp_mod
      USE moler_mod
      IMPLICIT NONE
      INTEGER, PARAMETER:: n_dim=251
      COMPLEX, DIMENSION(n_dim,n_dim)::zmn
      COMPLEX, DIMENSION(n_dim)::current, voltage
      COMPLEX:: zin
      REAL:: lambda,k,a,delta,length,zp
      INTEGER::m,n,n_seg,n_gen, i,j, unitNumber
      lambda=1.
      k=tpi/lambda
      a=0.005*lambda
      length= 0.47*lambda
      n_seg=41 !number of wire segments
      IF(n_seg .EQ. n_seg/2*2)n_seg=n_seg+1
      delta=length/n_seg !wire segment length
c
      CALL getzmn(k,a,n_seg,delta,length,zmn) !Find zmn(m,n)
      DO i = 1, n_seg   
        DO j = 1, n_seg
            !PRINT *, "zmn(", i, ",", j, ") = ", zmn(i,j)
        ENDDO
      ENDDO
      !PRINT *,"zmn(n_seg,n_seg)=", zmn(n_seg,n_seg)
      n_gen=n_seg/2+1
      voltage=(0.,0.)
      voltage(n_gen)= 1/delta
c
      CALL molerLU(n_seg,zmn)
      CALL molerSolve(n_seg,zmn,voltage,current) !current = a_n
c
      Do i = 1, n_seg
           PRINT *, "current(i",i,")=",ABS(current(i))
      ENDDO

      OPEN(UNIT=102, FILE='currentsMoM.txt', STATUS='UNKNOWN')
      
      DO i = 1, n_seg
            
            zp=-length/2.0+length/(2.0*n_seg)+(i-1)*length/n_seg
            WRITE(102, *) zp, ABS(current(i))
        ENDDO
        
        ! Close the file
        CLOSE(102)

      print *, "current(n_gen)=",current(n_gen)
      zin=1./current(n_gen)
      PRINT *, "zin = ", zin
      WRITE(6,1200)zin
      STOP
 1200 FORMAT(' ','R=',G11.3,3X,'X=',G11.3)
      STOP
      END PROGRAM wire
c
