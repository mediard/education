C MODULE constants_mod
C Various constants
C Copyright (c) 2015 by Robert Paknys
C Requires: none
      MODULE constants_mod
      COMPLEX, PARAMETER:: cj=(0.,1.)
      REAL, PARAMETER :: pi=3.141592653589793238462643, tpi=2*pi,
     ? sqrtpi=1.772453850905516027298167,
     ? dpr=180/pi,
     ? mu0=4E-7*pi, cspeed=2.99792458E8, epsilon0=1/(mu0*cspeed**2),
     ? eta0=1/(cspeed*epsilon0)
      END MODULE constants_mod
