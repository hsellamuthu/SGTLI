C     -------------------------------------------------------------
C     Orbit generation using Numerical integration for lunisolar gra
C     -vity, j2 to j4, drag and SRP  effect using KS elements (Earth).
C     -------------------------------------------------------------
C     If perigee altitude is less than 500 km drag is on
C     -------------------------------------------------------------
c     Code purpose:
C     Sun gravity-assist to trans-lunar injection orbits
c     July 31, 2018
C     -------------------------------------------------------------
C     Authors: Harishkumar Sellamuthu, Dr. Ram Krishan Sharma
C     Department of Aerospace Engineering, Karunya Institute of Technology and Sciences,
C     Coimbatore, India.
C     -------------------------------------------------------------
c     Updates:
c     date is entered in a new line (2nd) in the dinp.dat,
c     printing intermediate output data,
c     Modified based on E0-Ef with respect to Ek0-EKf,
c     Writing only KS and BR,
c     DeltaV - 4D map,
c     Direct Transfer to TLI,
c     fixing only perigee, apogee, eccentricity and mean anomaly,
c     WRITING ONLY MAX. ECCENTRICITY VALUES
C     -------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      COMMON /UW/ EC,E,DE
      COMMON /XY/ PI,CD,CR,RE
      COMMON /SUN/ ANG,DIP,GAM
      COMMON /TINT/ AJ2,AJ3,AJ4,AJ5,AK2,AKS2,AKM2
      common /rarp/ rpin,rain      
      DIMENSION U(4),US(4),Y(10),Q(4),Z(10),SP(3),X(3),P(3),XD(3)
      DIMENSION QQ(4),UU(4),PS(4),PL(4),SL(3),PQ(4) !
      DIMENSION ALT(300),SCH(300),DEN(300),HHT(22),RRA(22),O(3)
      DIMENSION DS(4),SRP(4),XR(3),PR(3)
C     -------------------------------------------------------------      
      DIMENSION DKEP(1000,14),DPV(1000,6),DSO(1000,5),DLU(1000,5)
      DIMENSION DEC(1000),DVEL(2)
C     -------------------------------------------------------------
      
      OPEN(UNIT = 1, FILE='ATM.DAT')
      OPEN(UNIT = 2, FILE='CONST.DAT')
      OPEN(UNIT = 3, FILE='DINP.DAT')
      OPEN(UNIT = 5, FILE='PV.OUT')
      OPEN(UNIT = 8, FILE='KEPNEW.OUT')
      OPEN(UNIT = 10, FILE = 'LEPH.OUT')
      OPEN(UNIT = 11, FILE = 'SANG.OUT')
      
C     *** CONSTANTS  ***
      N=10
      NHT=22
      NDIM=291
      N=10
      IO=1
      LP=0
      kl=0

!     ATMOSPHERIC PARAMETERS
      ALT(1)=60.D0
      DO 102 I=2,141
         I1=I-1
         ALT(I)=ALT(I1)+1.D0
 102  CONTINUE

      DO 103 I=142,NDIM
         I1=I-1
         ALT(I)=ALT(I1)+2.D0
 103  CONTINUE

      READ(1,100)(SCH(I),I=1,61)
      READ(1,100)(SCH(I),I=62,NDIM)
      READ(1,101)(DEN(I),I=1,61)
      READ(1,101)(DEN(I),I=62,NDIM)

      READ(2,*) SK,WE,EPS,ACD

!     GRAVITATIONAL PARAMETERS
      READ(2,*) AK2,AKS2,AKM2,RE
      READ(2,*) AJ2,AJ3,AJ4,AJ5

!     SOLAR PARAMETERS
      READ(2,*) PSR

      DO 105 I=1,NDIM
         ALT(I)=ALT(I)*1000.D0
         DEN(I)=DEN(I)*1.0D10
         SCH(I)=SCH(I)*1000.0D0
 105  CONTINUE

      AK=SQRT(AK2)
      PI=4.D0*DATAN(1.D0)
      REQ=RE*1.D3
      AMUE= AK2*1.0D9         !3.986004415D14 
      
      RE2=RE*RE
      RE3=RE*RE2
      RE4=RE2*RE2
      RE5=RE3*RE2
      CR=PI/180.0D0 
      CD=180.0D0/PI 

      
 301  IF(IO.EQ.1) GO TO 5
      READ (3,*,END=99) X1,X2,X3,X1DOT,X2DOT,X3DOT,JYR,JMON,
     * JDAY,JH,JMIN,JSEC,MM

      CALL CPVOE(SA,EC,AI,OMGA,SOMGA,EA,AM)

      GO TO 832
      
 5    READ(3,*,END=99) HSQ,EC,AI,OMGA,SOMGA,AM
      READ(3,*) JYR,JMON,JDAY,JH,JMIN,JSEC
      CALL JDATE(JYR,JMON,JDAY,JH,JMIN,JSEC,DDJ)

      READ(3,*)BN,ISTEP,NREV
      READ(3,*)FR,IDRAG
      READ(3,*)CRR,AMR,IPSR,ISHAD

      READ(3,*) ISUN,IMOON,IJ2,IJ3,IJ4,IJ5
      
      SA=(HSQ+RE)/(1.D0-EC)
      IF(IDRAG.EQ.0) GO TO 832
 833  SA=(HSQ+RE*(1.D0-EPS*(DSIN(AI*CR)**2*DSIN(SOMGA*CR)**2)))
     *     /(1.D0-EC)

 832  CONTINUE


C     SWITCH ON/OFF J2..J4, LUNISOLAR,SRP
      IF(ISUN.EQ.0) AKS2=0.D0
      IF(IMOON.EQ.0) AKM2=0.D0
      IF(IJ2.EQ.0) AJ2=0.D0
      IF(IJ3.EQ.0) AJ3=0.D0
      IF(IJ4.EQ.0) AJ4=0.D0
      IF(IJ5.EQ.0) AJ5=0.D0
      IF(IPSR.EQ.0) PSR=0.D0

      ddjin=ddj
      SAIN=SA
      ECIN=EC
      omgin=omga
      aiin=ai
      somgin=somga
      amin=am
      ha=(hsq+re)*(1.d0+ec)/(1.d0-ec)-re
      rpin = sain*(1.d0-ecin)
      rain = sain*(1.d0+ecin)
      aiin0=aiin
      omgin0=omgin
      somgin0=somgin
c      WRITE(8,39)T,IEDG,sAin,ECin,AIin,OMGin,SOMGin,amin,hsq,HA,LP      
 999  CALL COEPV(SAin,ECin,AIin,OMGin,SOMGin,AMin)
c      print *,omgin
      SA0=SA*1.D3
      XX0=SA0*EC
      
 64   E=0.0D0
      AT=0.0
      DE=2.0D0*PI/FLOAT(ISTEP)
      E=0.D0
      T=AT
      TT=0.0D0
      R=VMN(X1,X2,X3)
      R2=R*R
      RI=1.0D0/R
      RI2=RI*RI
      RI3=RI/R2
      RI4=RI2*RI2
      RI5=RI3*RI2
      RI6=RI3*RI3
      R4=R2*R2
      RDOT=VMN(X1DOT,X2DOT,X3DOT)
      RDOT2=RDOT*RDOT
      X32=X3*X3
      X3R=X3/R

      X32R2=X3R*X3R
      XR3=X3R*X32R2
      XR4=X32R2*X32R2
      XR5=XR3*X32R2

c      CALL CPVOE(SAIN,ECIN,AIIN,OMGAIN,SOMAIN,EAIN,AMIN)
      X(1)=X1
      X(2)=X2
      X(3)=X3

C     *** SUN & MOON  POSITION VECTORS ***

      CALL SOLARNPV(DDJ)
      CALL LUNARPV(DDJ)

C     *** SUN AND MOON GRAVITY POTENTIAL
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS3=AS2*AS
      AS4=AS3*AS
      AS41=1.0D0/AS4
      AS5=AS3*AS2
      AS51=1.0D0/AS5
      X1XS=X1*XS
      X2YS=X2*YS
      X3ZS=X3*ZS
      VS1=X1XS+X2YS+X3ZS
      VS2=VS1*VS1

      VAM=VMN(XM,YM,ZM)
      VAM2=VAM*VAM
      VAM3=VAM2*VAM
      VAM4=VAM3*VAM
      VAM41=1.0D0/VAM4
      VAM5=VAM3*VAM2
      VAM51=1.0D0/VAM5
      X1XM=X1*XM
      X2YM=X2*YM
      X3ZM=X3*ZM
      VM1=X1XM+X2YM+X3ZM
      VM2=VM1*VM1

C     *** PERTURBING POTENTIAL ***
      V=0.5d0*(-AKS2*AS51*(3.0D0*VS2-R2*AS2)
     *     -(AKM2*VAM51*(3.0D0*VM2-R2*VAM2))
     *     +AK2*RE2*AJ2*RI3*(3.D0*X32R2-1.D0)
     *     +AK2*RE3*AJ3*RI4*(5.D0*XR3-3.D0*X3R)
     *     +AK2*RE4*AJ4*0.25D0*RI5*(35.D0*XR4-30.D0*X32R2-3.D0)
     *     +AK2*RE5*AJ5*0.25D0*RI6*(63.D0*XR5-70.D0*XR3+15.D0*X3R))

      W=AK2/R-RDOT2/2.0D0-V 
      W2=W/2.0D0 
      W=DSQRT(W2)
!     Modification
      wk=dsqrt((AK2/R-RDOT2/2.0D0)/2.d0)
      W3=W2*W
      gam=w/wk
      IF (X1) 12,11,11
 11   U(4)=0.0D0
      RPX1=R+X1 
      U(1)=.5D0*RPX1-U(4)*U(4)
      U(1)=DSQRT(U(1)) 
      U(2)=(X2*U(1)+X3*U(4))/RPX1 
      U(3)=(X3*U(1)-X2*U(4))/RPX1 
      GO TO 13
 12   U(3)=0.0D0
      RMX1=R-X1 
      U(2)=DSQRT(.5D0*RMX1)
      U(1)=(X2*U(2)+X3*U(3))/RMX1 
      U(4)=(X3*U(2)-X2*U(3))/RMX1 
 13   CNST=1.0D0/(4.0D0*W)
      US(1)=(U(1)*X1DOT+U(2)*X2DOT+U(3)*X3DOT)*CNST 
      US(2)=(-U(2)*X1DOT+U(1)*X2DOT+U(4)*X3DOT)*CNST
      US(3)=(-U(3)*X1DOT-U(4)*X2DOT+U(1)*X3DOT)*CNST
      US(4)=(U(4)*X1DOT-U(3)*X2DOT+U(2)*X3DOT)*CNST
      
      TOW=T+(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W 

      NT=0

      Y(1)=TOW
      Y(10)=W
      DO 21 I=1,4
         Y(I+1)=U(I)
         Y(I+5)=2.0*US(I) 
 21   CONTINUE


C     *** 1/2W^2 *** CHANGED FROM 1/4W2 TO 1/2W2 (100)
      CH1=1.0D0/(2.0D0*W2)
      CH3=2.D0/W

C     *** 1/8W^3 ***
      CH2=1.0D0/(8.0D0*W3)
      EEE=EAIN*CR+E
      SIE=DSIN(E/2.0D0)
      CSE=DCOS(E/2.0D0)

C     CHECKS DURING INTEGRATION
      CHECK=Y(5)*Y(6)-Y(4)*Y(7)+Y(3)*Y(8)-Y(2)*Y(9) 
      DDJ1=DDJ
      X(1)=X1
      X(2)=X2
      X(3)=X3
C     ---------------------------- 
      CALL CPVOE(A,EC,AI,OMGA,SOMGA,EA,AM)
      RP=A*(1.D0-EC)
      RA=A*(1.D0+EC)
      HP=RP-RE
      HA=RA-RE
      EDG=0
      IEDG=(EDG)

      de=de*gam

      DO 107 LP=1,(NREV)

      CALL CPVOE(A,EC,AI,OMGA,SOMGA,EA,AM)
      OPEC=1.D0+EC
      E2=EC*EC
      EA=EA*CR
      VIPP=4.D0
      IF(EA.LT.PI) VIPP=2.D0
      DE=(VIPP*PI-EA)/FLOAT(ISTEP)
      A=A*1000
      AI=AI*CR
      CSI=DCOS(AI)
      SI=DSIN(AI)
      SOMGA=SOMGA*CR
      SOMGA2=2.0D0*SOMGA
      SN2W=DSIN(SOMGA2)
      CS2W=DCOS(SOMGA2)
      SIOMGA=DSIN(SOMGA)
      REQ=RE*1.D3
      RREQ=REQ*(1.0D0-SI*SI*SIOMGA*SIOMGA*EPS)
      SQ=A*(1.0D0-EC)
      HSQ=SQ-RREQ
C     PRINT *,HSQ

      IF(HSQ.GT.500.D3) C=0.D0
      IF(HSQ.GT.500.D3) GO TO 120
      CALL INTPOL(ALT,DEN,NDIM,HSQ,RHOPI)
      RHOP=1.D0
      RHOPI=RHOPI*RHOP
      CALL INTPOL(ALT,SCH,NDIM,HSQ,H)
      BETA=1.D0/H
      RPO=SQ
      VPO=AMUE*(2.0/RPO-1.0/A)
      VPO=DSQRT(VPO)
      F=1.0D0-RPO*WE*FR*CSI/VPO
      F=F*F
      XX=A*EC
      XXX=BETA*XX 
      C=EPS*BETA*RPO*SI*SI/2.0D0
 120  CONTINUE

      DO 122 KK=1,ISTEP
         CALL CPVOE(SAO,ECO,AIO,OMGAO,SOMGAO,EAO,AMO)

 111     USTAR2=0.0D0
         
         DO 31 I=1,4
            USTAR2=USTAR2+US(I)*US(I)
 31      CONTINUE
         

C     ENERGY EQUATION 
      EEQ=(R+4.0D0*USTAR2)/(AK2-R*V)-CH1 
      AH0=AK2/R-.5D0*RDOT2-V

 30   FORMAT(10X,"EEQ =",E20.10,10X,"BR =",E20.10)

      IF(ISHAD.EQ.0) THEN
         ROHP=0.D0
         GO TO 536
      ENDIF
      IF(ISHAD.EQ.2) GO TO 534
      ROHP=SHADFNCYL(X1,X2,X3,XS,YS,ZS)
      GO TO 536
 534  ROHP=SHADFNCONE(X1,X2,X3,XS,YS,ZS)
 536  CONTINUE


C     -----------------------------
C     R-K-G 4(5) INTEGRATOR LOOP 
C     -----------------------------
      
 60   CONTINUE

      XD(1)= X1DOT
      XD(2)= X2DOT
      XD(3)= X3DOT

C     DRAG WITH OBLATE ATMOSPHERIC MODEL TERMS

      SEA=DSIN(EA)
      CEA=DCOS(EA)
      TER=1.D0-EC*CEA
      E2=EC*EC
      CSEF=(CEA-EC)/TER
      SIEF=SEA*DSQRT(1.0D0-E2)/TER
      FA=ATAN3(SIEF,CSEF)
      TER2=DCOS(SOMGA)**2-DCOS(SOMGA+FA)**2
      TER=DEXP(-2.0D0*C*TER2-XXX*(1.0D0-CEA))
      RR=A*(1.0D0-EC*CEA)
      TER1=(RPO-RR)/H
      DENS=RHOPI*TER
      DENS=-F*DENS*RDOT/(BN*2.D0)
      DENS=DENS*1.D-7
      RQ=R-RE

!     DRAG ON/OFF CONDITION
      IF(IDRAG.EQ.0) DENS=0.D0

C     SOLAR RADIATION PRESSURE TERMS

      SRPC=CRR*AMR
      XR(1)=XS-X1
      XR(2)=YS-X2
      XR(3)=ZS-X3
      XRR=VMN(XR(1),XR(2),XR(3))

      DO I=1,3
         PR(I)=-ROHP*PSR*SRPC*XR(I)/XRR*1.D-3
      END DO

C     SOLAR GRAVITY TERMS
      SP(1)=XS
      SP(2)=YS
      SP(3)=ZS
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS3=AS2*AS
      AS4=AS3*AS
      AS41=1.0D0/AS4
      AS5=AS3*AS2
      AS51=1.0D0/AS5
      VS=X1*XS+X2*YS+X3*ZS
      VS2=VS*VS
      
      CNTS1=3.D0*VS2-3.d0*R2*AS2
      CNTS2=6.D0*R*VS
      CNTS3=CNTS1+2.D0*R2*AS2

C     LUNAR GRAVITY TERMS
      SL(1)=XM
      SL(2)=YM
      SL(3)=ZM
      VAM=VMN(XM,YM,ZM)
      VAM2=VAM*VAM
      VAM3=VAM2*VAM
      VAM4=VAM3*VAM
      VAM41=1.0D0/VAM4
      VAM5=VAM3*VAM2
      VAM51=1.0D0/VAM5
      VM=X1*XM+X2*YM+X3*ZM
      VM2=VM*VM
      
      CNTM1=3.D0*VM2-3.d0*R2*VAM2
      CNTM2=6.D0*R*VM
      CNTM3=CNTM1+2.D0*R2*VAM2

C     SOLAR GRAVITY POTENTIAL, PS      
      PS(1)=U(1)*SP(1)+U(2)*SP(2)+U(3)*SP(3)
      PS(2)=-U(2)*SP(1)+U(1)*SP(2)+U(4)*SP(3) 
      PS(3)=-U(3)*SP(1)-U(4)*SP(2)+U(1)*SP(3) 
      PS(4)=U(4)*SP(1)-U(3)*SP(2)+U(2)*SP(3)

C     LUNAR GRAVITY POTENTIAL, PL
      PL(1)=U(1)*SL(1)+U(2)*SL(2)+U(3)*SL(3)
      PL(2)=-U(2)*SL(1)+U(1)*SL(2)+U(4)*SL(3) 
      PL(3)=-U(3)*SL(1)-U(4)*SL(2)+U(1)*SL(3) 
      PL(4)=U(4)*SL(1)-U(3)*SL(2)+U(2)*SL(3)

C     SOLAR RADIATION PRESSURE, SRP
      SRP(1)=U(1)*PR(1)+U(2)*PR(2)+U(3)*PR(3)
      SRP(2)=-U(2)*PR(1)+U(1)*PR(2)+U(4)*PR(3) 
      SRP(3)=-U(3)*PR(1)-U(4)*PR(2)+U(1)*PR(3) 
      SRP(4)=U(4)*PR(1)-U(3)*PR(2)+U(2)*PR(3)

      DO I=1,3
         P(I)=DENS*XD(I)
      END DO

      QQ(1)=U(1)*P(1)+U(2)*P(2)+U(3)*P(3)
      QQ(2)=-U(2)*P(1)+U(1)*P(2)+U(4)*P(3)
      QQ(3)=-U(3)*P(1)-U(4)*P(2)+U(1)*P(3)
      QQ(4)=U(4)*P(1)-U(3)*P(2)+U(2)*P(3)

C     J2, J3, J4 POTENTIAL TERMS
c$$$      CNTJ1=6.D0*X3R
c$$$      CNTJ2=-12.D0*X32R2
c$$$      CNTK1=30.D0*X32R2-6.D0
c$$$      CNTK2=12.D0*X3R-30.D0*XR3
c$$$      CNTL1=-1.5D0*X3R+35.D0*XR3
c$$$      CNTL2=-70.D0*XR4+45.D0*X32R2-3.D0

      
!     ORIGINAL
      CNTJ1=3.D0*X3R            !
      CNTJ2=-6.D0*X32R2         !
      CNTK1=7.5D0*X32R2-1.5D0   !
      CNTK2=6.D0*X3R-15.D0*XR3  !
      CNTL1=(17.5D0*XR3-7.5D0*X3R) !
      CNTL2=-35.D0*XR4+22.5D0*X32R2-1.5D0 !
      CNTN1=39.375D0*XR4-26.25D0*X32R2+1.875D0 !
      CNTN2=-78.75D0*XR5+70.D0*XR3-11.25D0*X3R
!     CNTL1=(-15.D0*X3R+35.D0*XR3)*0.5d0!
!     CNTL2=-35.D0*XR4+22.5D0*X32R2-1.5D0!
!     CNTL1=(-15.D0*X3R+35.D0*XR3)*0.5d0!
!     CNTL2=-35.D0*XR4+22.5D0*X32R2-1.5D0!

      
      CNTJ=AK2*RE2*AJ2*RI3      !
      CNTK=AK2*RE3*AJ3*RI4      !
      CNTL=AK2*RE4*AJ4*RI5      !
      CNTN=AK2*RE5*AJ5*RI6      !
      
      UU(1)=U(3)
      UU(2)=U(4)
      UU(3)=U(1)
      UU(4)=U(2)
      
      RD2=R/2.D0
      
      C1=1.0D0/(2.0D0*W2) 
      C2=2.0D0/W

C     ----------------------------------      
C     DALFA/DE,DW/DE,DTAU/DE,DW/DE
C     ----------------------------------

      Z(10)=US(1)*QQ(1)+US(2)*QQ(2)+US(3)*QQ(3)+US(4)*QQ(4)
      Z(10)=Z(10)+(US(1)*SRP(1)+US(2)*SRP(2)+US(3)*SRP(3)+US(4)*SRP(4))
      Z(10)=-.5D0*Z(10)/W

      DO 52 I=1,4
         O(I)=2.D0*Z(10)*US(I)/W
 52   CONTINUE

      Z(1)=(AK2+2*AKS2*AS51*R*CNTS3+2*AKM2*VAM51*R*CNTM3
     *     +R*(CNTJ*0.5D0*(3.D0*X32R2-1.D0)+
     *     CNTK*0.5D0*(5.D0*XR3-3.D0*X3R)+
     *     CNTL*3.D0*0.125D0*(35.D0*XR4-30.D0*X32R2+3.D0)+
     *     CNTN*4.D0*0.125D0*(63.D0*XR5-70.D0*XR3+15.D0*X3R))+
     *     R*(U(1)*QQ(1)+U(2)*QQ(2)+U(3)*QQ(3)+U(4)*QQ(4))+
     *     R*(U(1)*SRP(1)+U(2)*SRP(2)+U(3)*SRP(3)+U(4)*SRP(4))) 
     *     /(8.D0*W*W2)-(Z(10)*2.D0/(W2))*(U(1)*US(1)+U(2)*US(2)
     *     +U(3)*US(3)+U(4)*US(4))
      

!     Original
      DO 51 I=1,4 

      Z(I+1)=((-AKS2*AS51*(CNTS1*U(I)+CNTS2*PS(I)))*0.5D0
     *        -(AKM2*VAM51*(CNTM1*U(I)+CNTM2*PL(I)))*0.5D0+
     *        (CNTJ*((1.D0+CNTJ2)*U(I)+CNTJ1*UU(I)))+
     *        CNTK*(U(I)*CNTK2+UU(I)*CNTK1)+
     *        CNTL*(U(I)*CNTL2+UU(I)*CNTL1)+
     *        CNTN*(U(I)*CNTN2+UU(I)*CNTN1))-R*QQ(I)-R*SRP(I)
      
      Z(I+5)=((-AKS2*AS51*(CNTS1*U(I)+CNTS2*PS(I)))*0.5D0
     *     -(AKM2*VAM51*(CNTM1*U(I)+CNTM2*PL(I)))*0.5D0+
     *     (CNTJ*((1.D0+CNTJ2)*U(I)+CNTJ1*UU(I)))+
     *     CNTK*(U(I)*CNTK2+UU(I)*CNTK1)+
     *     CNTL*(U(I)*CNTL2+UU(I)*CNTL1)+
     *     CNTN*(U(I)*CNTN2+UU(I)*CNTN1))-R*QQ(I)-R*SRP(I)
         
 51   CONTINUE

      DO I=1,4
         Z(I+1)=(Z(I+1)*0.25/W2+O(I))*SIE
         Z(I+5)=-(Z(I+5)*0.25/W2+O(I))*CSE
      END DO

C     ----------------------------------

      S=RKGIL(N,Y,Z,E,DE,NT)

      CSE=DCOS(E/2.D0)
      SIE=DSIN(E/2.D0)

      DO 71 I=1,4 
         U(I)=Y(I+1)*CSE+Y(I+5)*SIE
         US(I)=-.5D0*Y(I+1)*SIE+.5D0*Y(I+5)*CSE
 71   CONTINUE
      
      X1=U(1)*U(1)-U(2)*U(2)-U(3)*U(3)+U(4)*U(4)
      X2=2.0D0*(U(1)*U(2)-U(3)*U(4))
      X3=2.0D0*(U(1)*U(3)+U(2)*U(4))
      R=U(1)*U(1)+U(2)*U(2)+U(3)*U(3)+U(4)*U(4)

      R2=R*R 
      R4=R2*R2
      CST1=4.0D0*W/R
      X32=X3*X3
      X3R=X3/R
      X32R2=X3R*X3R
      XR3=X3R*X32R2
      XR4=X32R2*X32R2
      XR5=XR3*X32R2

      X1DOT=CST1*(U(1)*US(1)-U(2)*US(2)-U(3)*US(3)+U(4)*US(4))
      X2DOT=CST1*(U(2)*US(1)+U(1)*US(2)-U(4)*US(3)-U(3)*US(4))
      X3DOT=CST1*(U(3)*US(1)+U(4)*US(2)+U(1)*US(3)+U(2)*US(4))
      RDOT=VMN(X1DOT,X2DOT,X3DOT)
      RDOT2=RDOT*RDOT
      X(1)=X1
      X(2)=X2
      X(3)=X3

      RI=1.0D0/R
      RI2=RI*RI
      RI3=RI*RI2
      RI4=RI2*RI2
      RI5=RI3*RI2
      RI6=RI3*RI3
      W=Y(10)
      W2=W*W
      T=Y(1)-(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W
      CALL CTHMS(T,DH,DM,SECONDS)
      CALL CSMHJD(DH,DM,SECONDS,DJT)
      DDJ1=DJT+DDJ

      CALL SOLARNPV(DDJ1)
      CALL LUNARPV(DDJ1)

C     *** SUN AND MOON POSITION VECTORS ***  
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS5=AS**5
      AS51=1.0D0/AS5
      X1XS=X1*XS
      X2YS=X2*YS
      X3ZS=X3*ZS
      VS1=X1XS+X2YS+X3ZS
      VS2=VS1**2
      RRS2=(R*AS)**2

      VAM=VMN(XM,YM,ZM)
      VAM2=VAM*VAM
      VAM5=VAM**5
      VAM51=1.0D0/VAM5
      X1XM=X1*XM
      X2YM=X2*YM
      X3ZM=X3*ZM
      VM1=X1XM+X2YM+X3ZM
      VM2=VM1**2
      RRM2=(R*VAM)**2

      V=0.5d0*(-AKS2*AS51*(3.0D0*VS2-R2*AS2)
     *     -(AKM2*VAM51*(3.0D0*VM2-R2*VAM2))
     *     +AK2*RE2*AJ2*RI3*(3.D0*X32R2-1.D0)
     *     +AK2*RE3*AJ3*RI4*(5.D0*XR3-3.D0*X3R)
     *     +AK2*RE4*AJ4*0.25D0*RI5*(35.D0*XR4-30.D0*X32R2-3.D0))
     *     +AK2*RE5*AJ5*0.25D0*RI6*(63.D0*XR5-70.D0*XR3+15.D0*X3R)

      
      IF (S-1.) 60,60,70
  70  CONTINUE

C     ---------------------------------------
C     *** TIME ***
      T=Y(1)-(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W

      CALL CTHMS(T,DH,DM,SECONDS)
      CALL CSMHJD(DH,DM,SECONDS,DJT)
      DDJ1=DJT+DDJ

      CALL SOLARNPV(DDJ1)
      CALL LUNARPV(DDJ1)
      EAR=EA
      CALL CPVOE(A,EC,AI,OMGAD,SOMGAD,EA,AM)

      AH=AK2/R-.5D0*RDOT2-V 
      CHECK=Y(5)*Y(6)-Y(4)*Y(7)+Y(3)*Y(8)-Y(2)*Y(9)
      EDG=E*CD

      IEDG=ANINT(EDG)
c      EA=EAR+DE

      RP=A*(1.D0-EC)
      RA=A*(1.D0+EC)
      HP=RP-RE
      HA=RA-RE

      IF((R-RE).LE.80) GO TO 99

 122  CONTINUE
      adj=ddj-ddjin 

      DKEP(LP,1)=T
      DKEP(LP,2)=ADJ
      DKEP(LP,3)=AIIN
      DKEP(LP,4)=OMGIN
      DKEP(LP,5)=SOMGIN
      DKEP(LP,6)=A
      DKEP(LP,7)=EC
      DKEP(LP,8)=AI
      DKEP(LP,9)=OMGAD
      DKEP(LP,10)=SOMGAD
      DKEP(LP,11)=AM
      DKEP(LP,12)=HP
      DKEP(LP,13)=HA
      DKEP(LP,14)=LP
      CALL SOLARNPV(DDJ)     
      DSO(LP,1)=GAM
      DSO(LP,2)=DDJ
      DSO(LP,3)=DDJ1
      DPV(LP,1)=X1
      DPV(LP,2)=X2
      DPV(LP,3)=X3
      DPV(LP,4)=X1DOT
      DPV(LP,5)=X2DOT
      DPV(LP,6)=X3DOT
      DEC(LP)=EC-ECIN

 107  CONTINUE


      INDE=1
      DECM=DEC(1)
      DO LPS=2,NREV
         IF(DEC(LPS).GT.DECM) THEN
            DECM=DEC(LPS)
            INDE=LPS
         ENDIF
      END DO

      CALL LUNEPR(DPV(INDE,:),DDJ1,SCM,DVEL,KIO)

      IF(KIO.EQ.1) THEN
         WRITE(8,39) DKEP(INDE,1),DKEP(INDE,2)*12.D0,
     *        ((INT(DKEP(INDE,INJ))),INJ=3,5),
     *        (DKEP(INDE,INK),INK=6,13),INT(DKEP(INDE,14))
         WRITE(11,1261) (DSO(INDE,INL),INL=1,3),SCM,DVEL(1),DVEL(2)
      ENDIF

!     Looping
C     ddj=ddj+1.d0
      ddj=ddj+0.041666667D0   !every hour
      if(adj.gt.365.d0) go to 1001
      go to 999
 1001 ddj=ddjin
      aiin=aiin+1.d0
      if(aiin.gt.28.d0) go to 1002
      go to 999
 1002 aiin=aiin0
      omgin=omgin+1.d0
      if(omgin.gt.359.d0) go to 1003
      go to 999
 1003 omgin=omgin0
      print *,somgin
      somgin=somgin+1.d0
      if(somgin.gt.359.d0) go to 301

      go to 999

!     FORMAT SPECIFIERS
 10   FORMAT(F18.10,1x,F18.10,1x,F18.10,1x,F20.8)
 1261 FORMAT(F18.10,1x,F18.10,1x,F18.10,1x,F18.10,1x,F18.10,1x,F18.10)      
 15   FORMAT(2X,4F20.5) 
 20   FORMAT(2X,3F18.10,4F15.8)                              
 39   FORMAT(F20.8,1X,F20.8,3I10.1,1X,F18.10,F14.10,F16.10,3F20.8
     & ,2F18.10,I10.1)
 40   FORMAT(F18.10,1X,F18.10,1X,F18.10,1X,F18.10,1X,F18.10,1X,
     &     F18.10,1X)
 84   FORMAT(2X,E20.6,5F15.4,E20.5)
 100  FORMAT(8F10.3)
 101  FORMAT(8E10.3)
 117  FORMAT(2X,6F18.9)
 127  FORMAT(1X,F8.2,F16.10,F18.10,F15.10,F14.10,4F16.10)
 212  FORMAT(F20.9,1X,F20.9,1X,F20.9,1X,F20.9)
 213  FORMAT(6F20.9)
 239  FORMAT(F20.8,1X,I10.1,1X,4F18.10)
 333  FORMAT(10X,F10.4,2I10,2F9.2)
 768  FORMAT(1X,7F18.10)
 946  FORMAT(F20.8,1X,9F20.9,E20.5)
 1010 FORMAT(A,6I6.2)
 1011 FORMAT(A,I10.2)
      
 99   STOP

      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(5)
      CLOSE(8)
      CLOSE(10)
      CLOSE(11)

      END

C     -------------------------------------------------------------
      FUNCTION RKGIL(N,Y,F,X,H,NT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(10),F(10),Q(10) 
      NT=NT+1 
      GO TO (1,2,3,4),NT
    1 DO 11 J=1,N 
   11 Q(J)=0.0D0
         A=0.5D0
      X=X+H/2.0D0 
      GO TO 5 
    2 A=0.29289321881D0 
      GO TO 5 
    3 A=1.7071067812D0
      X=X+H/2.0D0 
      GO TO 5 
    4 DO 41 I=1,N 
  41  Y(I)=Y(I)+H*F(I)/6.0D0-Q(I)/3.0D0 
      NT=0
      RKGIL=2.0D0 
      GO TO 6 
    5 DO 51 L=1,N 
      Y(L)=Y(L)+A*(H*F(L)-Q(L)) 
  51  Q(L)=2.0D0*A*H*F(L)+(1.0D0-3.0D0*A)*Q(L)
      RKGIL=1.0D0 
    6 CONTINUE
      RETURN
      END 
C     -------------------------------------------------------------
      SUBROUTINE INTPOL(X,Y,M1,X1,Y1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION X(M1),Y(M1)
      IF(X1.LT.X(1).OR.X1.GT.X(M1)) GO TO 50
      DO 11 I=1,M1
         IF(X1-X(I))11,13,10
 10      IF(X1-X(I+1))15,11,11
 15      Y1=Y(I)+(Y(I+1)-Y(I))*(X1-X(I))/(X(I+1)-X(I))
         GO TO 12
 13      Y1=Y(I)
 11   CONTINUE
 12   RETURN
 50   CONTINUE!WRITE(7,600)X
 600  FORMAT(/3X,"X1 IS OUT OF RANGE",E15.6)
      STOP
      END
C     -------------------------------------------------------------
      SUBROUTINE CPVOE(SA,E,AI,OMGAD,SOMGAD,EA,AM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3),XD(3),A(3),B(3)
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      COMMON /XY/ PI,CD,CR,RE
      X(1)= X1
      X(2)= X2
      X(3)= X3
      XD(1)= X1DOT
      XD(2)= X2DOT
      XD(3)= X3DOT
      AK2=3.986004415D05
      AK=DSQRT(AK2)
      R=VMN(X(1),X(2),X(3))
      R2=R*R
      V=VMN(XD(1),XD(2),XD(3))
      V2=V*V
      
      XXD=X(1)*XD(1)+X(2)*XD(2)+X(3)*XD(3)
      SA=2.0D0/R-V2/AK2 
      SA=1.0D0/SA 
      RTA=DSQRT(SA)
      RA1=1.0D0-R/SA
      E=RA1*RA1+XXD*XXD/(AK2*SA)
      E=DSQRT(E) 
      CSE=RA1/E 
      SIE=XXD/(AK*RTA*E)
      EA=ATAN3(SIE,CSE) 
      AM=EA-E*SIE 
      EA=EA*CD 
      AM=AM*CD
      CF1=CSE/R 
      CF2=RTA*SIE/AK
      E21=DSQRT(1.0D0-E*E) 
      CF3=SIE/(R*E21) 
      CF4=(CSE-E)*RTA/(AK*E21)
      DO 22 I=1,3 
      A(I)=CF1*X(I)-CF2*XD(I) 
      B(I)=CF3*X(I)+CF4*XD(I) 
   22 CONTINUE
      A1B2=A(1)+B(2)
      A2B1=A(2)+B(1)
      A2MB1=A(2)-B(1) 
      A1MB2=A(1)-B(2) 
      ALMD=A1B2*A1B2+A2MB1*A2MB1
      ALMD=DSQRT(ALMD) 
      AMUE=A1MB2*A1MB2+A2B1*A2B1
      AMUE=DSQRT(AMUE) 
      CSI=(ALMD-AMUE)/2.0D0 
      SI=DSQRT(A(3)*A(3)+B(3)*B(3))
       AI=ATAN3(SI,CSI) 
      AI=AI*CD
      SSOMGA=A(3)/SI
      CSOMGA=B(3)/SI
      SOMGA=ATAN3(SSOMGA,CSOMGA)
      SOMGAD=SOMGA*CD 
      OMGAS=(A(1)*A(3)+B(1)*B(3))/((A(2)*B(1)-A(1)*B(2))*SI)
      OMGAC=(A(2)*A(3)+B(2)*B(3))/((A(1)*B(2)-A(2)*B(1))*SI)
      OMGA=ATAN3(OMGAS,OMGAC) 
      OMGAD=OMGA*CD 
      RETURN
      END 
C     -------------------------------------------------------------
      FUNCTION ATAN3(A,B) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PI=4.0D0*DATAN(1.0D0)
      ATAN3=DATAN(DABS(A/B))
      IF((B.LT.0.).AND.(A.GT.0.))ATAN3=PI-ATAN3 
      IF((B.LT.0.).AND.(A.LT.0.))ATAN3=PI+ATAN3 
      IF((B.GT.0.).AND.(A.LT.0.))ATAN3=2.0D0*PI-ATAN3 
      RETURN
      END 
C     -------------------------------------------------------------
      SUBROUTINE COEPV(SA,E,AI,OMGA,SOMGA,AM) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     CONVERSION OF ORBITAL ELEMENTS TO POSITION AND VELOCITY COMPONENTS
      DIMENSION X(3),XD(3) 
      DIMENSION A(3),B(3)
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      COMMON /XY/ PI,CD,CR,RE
      CR=PI/180. 
      AK=3.986004415D05
      AK=DSQRT(AK) 
      AIR=AI*CR 
      OMGAR=OMGA*CR 
      SOMGAR=SOMGA*CR 
      CSOMGA=DCOS(OMGAR) 
      SIOMGA=DSIN(OMGAR) 
      SOMGAC=DCOS(SOMGAR)
      SOMGAS=DSIN(SOMGAR)
      CSI=DCOS(AIR)
      SI=DSIN(AIR) 
      A(1)=CSOMGA*SOMGAC-SIOMGA*SOMGAS*CSI
      A(2)=SIOMGA*SOMGAC+CSOMGA*SOMGAS*CSI
      A(3)=SOMGAS*SI
      B(1)=-CSOMGA*SOMGAS-SIOMGA*SOMGAC*CSI 
      B(2)=-SIOMGA*SOMGAS+CSOMGA*SOMGAC*CSI 
      B(3)=SOMGAC*SI
      AMR=AM*CR
      CALL ECCENT(AMR,E,EA) 
      EAD=EA*CD 
C     PRINT 33,EAD
   33 FORMAT(10X,"EA =",F20.5)
      CSE=DCOS(EA) 
      SIE=DSIN(EA) 
      R=SA*(1.0D0-E*CSE)
      RTE=DSQRT(1.0D0-E*E) 
      DIN=AK*DSQRT(SA)/R 
      AA=SA*(CSE-E) 
      BB=SA*RTE*SIE 
      DO 11 I=1,3 
      X(I)=AA*A(I)+BB*B(I)
 11   CONTINUE 
      DO 22 I=1,3 
         XD(I)  =DIN*(RTE*CSE*B(I)-SIE*A(I)) 
 22   CONTINUE
C     PRINT 40,(X(I),I=1,3) 
C     PRINT 50,(XD(I),I=1,3)
 40   FORMAT(5X,"POSITION",3F20.10)
 50   FORMAT(5X,"VELOCITY  ",3F20.10) 
      X1=X(1)
      X2=X(2) 
      X3=X(3) 
      X1DOT=XD(1) 
      X2DOT=XD(2) 
      X3DOT=XD(3) 
      RETURN
      END
C     -------------------------------------------------------------
      SUBROUTINE ECCENT(A,B,C)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EL=B
      XINOT=0.0 
      K=1 
      ANOMO=ABS(DSIN(A)) 
    1 CONTINUE
      XIONE=B*DSIN(A+XINOT)
      EL=EL*B 
      REMIND=(EL/(1.0-B))*ANOMO 
      IF(REMIND .LE. 0.00000000000001) GO TO 2
      K=K+1 
      XINOT=XIONE 
      GO TO 1 
    2 CONTINUE
      C=A+XIONE 
      RETURN
      END
C     -------------------------------------------------------------
C     CONVERSION FROM CALENDAR DATE AND TIME TO JULIAN DATE
      SUBROUTINE JDATE(IYR,IMON,IDAY,IH,IMIN,ISEC,DDJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTMON1=(IMON+9.0D0)/12.0D0
      INTMON2=(275.0D0*IMON)/9.0D0
      JYR=(7.0D0*(IYR+INTMON1))/4.0D0
      TYR1=367.0D0*IYR
      FT1=(ISEC/60.0D0)
      FT2=FT1+IMIN
      FT3=FT2/60.0D0
      FT4=FT3+IH
      FTIME=FT4/24.0D0
      CCDAY1=IDAY+1721013.5D0
      DDJ=TYR1-JYR+INTMON2+CCDAY1+FTIME
      RETURN
      END
C     -------------------------------------------------------------
C     CONVERSION FROM  JULIAN DATE TO CALENDAR DAY
c     Meeus - Astronomical Algorithm
      SUBROUTINE JUL2D(DDJ,IYR,IMON,IDAY,IH,IMIN,ISEC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      ddj5=ddj+0.5d0
      fint = int(ddj5)
      fdec=ddj5-fint
      if(fint.lt.2299161) then
         acc=fint
      else
         alfa=int((fint-1867216.25)/36524.25)
         acc=fint+1+alfa-int(alfa/4)
      end if
      bcc=acc+1524
      iccc=int((bcc-122.1)/365.25)
      idcc=int(365.25*iccc)
      iecc=int((bcc-idcc)/30.6001)
      day=bcc-idcc-int(30.6001*iecc)+fdec
      if(iecc.lt.14) then
         imon=iecc-1
      else
         imon=iecc-13
      end if
      if(imon.gt.2) then
         iyr=iccc-4716
      else
         iyr=iccc-4715
      end if
      iday=int(day)
      call day2hms(day,ih,imin,isec)
      RETURN
      END
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C     CONVERSION FROM  JULIAN DATE TO CALENDAR DAY
      SUBROUTINE JUL3D(DDJ,IYR,IMON,IDAY,IH,IMIN,ISEC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      fint=int(ddj)
      fdec=ddj-fint
      if(fint.lt.2299161) then
         acc=fint
      else
         alfa=int((fint-1867216.25)/36524.25)
         acc=fint+1+alfa-int(alfa/4)
      end if
      bcc=acc+1524
      iccc=int((bcc-122.1)/365.25)
      idcc=int(365.25*iccc)
      iecc=int((bcc-idcc)/30.6001)
      day=bcc-idcc-int(30.6001*iecc)+fdec
      if(iecc.lt.14) then
         imon=iecc-1
      else
         imon=iecc-13
      end if
      if(imon.gt.2) then
         iyr=iccc-4716
      else
         iyr=iccc-4715
      end if
      iday=int(day)
      call day2hms(day,ih,imin,isec)
      
      RETURN
      END
C     -------------------------------------------------------------
C     DAYS TO HOURS, MINUTES AND SECONDS
      SUBROUTINE DAY2HMS(D,IH,IMN,ISC)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      iday=int(d)
      frday=d-iday
      ih=int(frday*24)
      amn=frday*24-ih
      imn=int(amn*60)
      asc=amn*60-imn
      isc=int(asc*60)
      if(isc.eq.60) then
         isc=0
         imn=imn+1
      end if
      if(imn.eq.60) then
         imn=0
         ih=ih+1
      end if
      RETURN
      END
C     -------------------------------------------------------------      
C     SECONDS TO HOURS, MINUTES AND SECONDS
      SUBROUTINE CTHMS(TIME,DHOUR,DMIN,SECS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DHOUR=0.0D0
C      IHOUR=INT(TIME/3600)
C      DHOUR=IHOUR
      TIME=TIME-60.0D0*DHOUR
C      WRITE(9,*) DHOUR,TIME
      DMIN=DINT(TIME/60.0D0)
C      PRINT *,DMIN
      SECS=TIME-60.0D0*DMIN
C      PRINT *,SECS
C      IF(IHOUR.EQ.1) STOP
      RETURN
      END
C     -------------------------------------------------------------
C     SECONDS,MINUTES,HOURS TO JULIAN DATE
      SUBROUTINE CSMHJD(DHOUR,DMINUTE,SECONDS,DJT)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DHOUR=0.0D0
      S1=SECONDS/60.0D0
      S2=(S1+DMINUTE)/60.0D0
      DJT=(S2+DHOUR)/24.0D0
      RETURN
      END
C     -------------------------------------------------------------
      FUNCTION ATAN2(A,B) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     PI=3.1415926535897932384D0 
      PI=4.0D0*DATAN(1.0D0)
      PIB4=PI/4.0D0
      PI3B4=PIB4*3.0D0
      ABSA=DABS(A)
      IF(B.GE.0.0D0) ATAN2=PIB4-PIB4*((B-ABSA)/(B+ABSA))
      IF(B.LT.0.0D0) ATAN2=PI3B4-PIB4*((B+ABSA)/(ABSA-B))
      IF(A.GT.0.0D0) ATAN2=-ATAN2
      RETURN
      END 
C     -----------------------------------------------------------
C     LUNAR VECTOR COMPUTATION
C     1999 GSFC ISSFD PAPER
      SUBROUTINE LUNARPV(DJ)  
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      COMMON /XY/ PI,CD,CR,RE
   
C     TIME IN JULIAN CENTURIES
      TLUN=(DJ-2451545.0D0)/36525.0D0

C     X COORDINATE
      A11=383.0D0*DSIN(8399.685D0*TLUN+5.381D0)
      A12=31.5D0*DSIN(70.99D0*TLUN+6.619D0)
      A13=10.6D0*DSIN(16728.377D0*TLUN+1.453D0)
      A14=6.2D0*DSIN(1185.622D0*TLUN+0.481D0)
      A15=3.2D0*DSIN(7143.070D0*TLUN+5.017D0)
      A16=2.3D0*DSIN(15613.745D0*TLUN+0.857D0)
      A17=0.8D0*DSIN(8467.263D0*TLUN+1.01D0)

      XM=(A11+A12+A13+A14+A15+A16+A17)*1000.0D0

C     Y COORDINATE
      A21=351.0D0*DSIN(8399.687D0*TLUN+3.811D0)
      A22=28.9D0*DSIN(70.997D0*TLUN+4.596D0)
      A23=13.7D0*DSIN(8433.466D0*TLUN+4.766D0)
      A24=9.7D0*DSIN(16728.380D0*TLUN+6.165D0)
      A25=5.7D0*DSIN(1185.667D0*TLUN+5.164D0)
      A26=2.9D0*DSIN(7143.058D0*TLUN+0.3D0)
      A27=2.1D0*DSIN(15613.755D0*TLUN+5.565D0)

      YM=(A21+A22+A23+A24+A25+A26+A27)*1000.0D0

C     Z COORDINATE
      A31=153.2D0*DSIN(8399.672D0*TLUN+3.807D0)
      A32=31.5D0*DSIN(8433.464D0*TLUN+1.629D0)
      A33=12.5D0*DSIN(70.996D0*TLUN+4.595D0)
      A34=4.2D0*DSIN(16728.364D0*TLUN+6.162D0)
      A35=2.5D0*DSIN(1185.645D0*TLUN+5.167D0)
      A36=3.0D0*DSIN(104.881D0*TLUN+2.555D0)
      A37=1.8D0*DSIN(8399.116D0*TLUN+6.248D0)

      ZM=(A31+A32+A33+A34+A35+A36+A37)*1000.0D0

      RETURN
      END
C     ----------------------------------------
C     SOLAR VECTOR COMPUTATION
C     COMPUTING SOLAR VECTOR PAPER (PSA ALGORITHM)
      SUBROUTINE SOLARNPV(DJ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z,M)
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      COMMON /XY/ PI,CD,CR,RE
      COMMON /SUN/ ANG,DIP,GAM
      
C     JULIAN DAY
      DJN=(DJ-2451545.0D0)

C     OMEGA
      OMEGA=2.1429D0-0.0010394594D0*DJN

C     MEAN LONGITUDE OF THE SUN
      TLAM=(4.8950630D0+0.017202791698D0*DJN)

C     MEAN ANOMALY OF THE SUN
      TMMM=(6.24006D0+0.0172019699D0*DJN)

C     ECLIPTIC LONGITUDE

      TECL=TLAM+0.03341607D0*DSIN(TMMM)   
      TECL=TECL+0.00034894D0*DSIN(2*TMMM)
      TECL=TECL-0.0001134D0
      TECL=TECL-0.0000203D0*DSIN(OMEGA)

C     OBLIQUITY OF THE ECLIPTIC
      TOBL=0.4090928D0
      TOBL=TOBL-(6.214D-09)*DJN
      TOBL=TOBL+0.0000396D0*DCOS(OMEGA)

C     RIGHT ASCENSION AND DECLINATION
      DEO=DSIN(TECL)*DCOS(TOBL)
      DEC=DCOS(TECL)
      RALFA = ATAN2(DEO,DEC)
      IF(RALFA.LT.0.0D0) RALFA=RALFA+2.0D0*PI
      RDECL=DASIN(DSIN(TOBL)*DSIN(TECL))

      SRALFA=DSIN(RALFA)
      CRALFA=DCOS(RALFA)
      SRDECL=DSIN(RDECL)
      CRDECL=DCOS(RDECL)
      TRS=1.495978707D08
      XS=TRS*CRALFA*CRDECL
      YS=TRS*SRALFA*CRDECL
      ZS=TRS*SRDECL

      CALL CPVOE(A,EC,AI,OMGAD,SOMGAD,EA,AM)
      AI=AI*CR
      OMG=OMGAD*CR
      SOMG=SOMGAD*CR
      
      CSW=DCOS(SOMG)
      SSW=DSIN(SOMG)
      CBW=DCOS(OMG)
      SBW=DSIN(OMG)
      SI=DSIN(AI)
      CI=DCOS(AI)

      AL1=CSW*CBW-CI*SSW*SBW
      AM1=CSW*SBW+CI*SSW*CBW
      AN1=SSW*SI
      
      AL2=-SSW*CBW-CI*CSW*SBW
      AM2=CI*CSW*CBW-SSW*SBW
      AN2=CSW*SI

      AL3=SI*SBW
      AM3=-SI*CBW
      AN3=CI

      ALS=CRALFA*CRDECL
      AMS=SRALFA*CRDECL
      ANS=SRDECL

      AA=AL1*ALS+AM1*AMS+AN1*ANS
      BB=AL2*ALS+AM2*AMS+AN2*ANS
      CC=AL3*ALS+AM3*AMS+AN3*ANS

      ANG=DACOS(AA)*CD
      DIP=DASIN(CC)*CD
      GAM=DATAN2(bb,aa)
      IF(GAM.LT.0.D0) GAM=GAM+2.D0*PI
      GAM=GAM*CD
      

      RETURN
      END
C     ----------------------------------------
      FUNCTION SHADFNCYL(X1,X2,X3,XS,YS,ZS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z,M)
      COMMON /XY/ PI,CD,CR,RE
C     ----------------------------------------      
C     CYLINDRICAL SHADOWN FUNCTION OBTAINED FROM
C     HUBAUX ET AL., 2012
C     ----------------------------------------      
      X=X1
      Y=X2
      Z=X3
      GAM=1.D9
      RS=VMN(XS,YS,ZS)
      R=VMN(X,Y,Z)
      COSPHI=(X*XS+Y*YS+Z*ZS)/RS
      SCR=COSPHI+DSQRT(R*R-RE*RE)
      SHADFNCYL=0.5D0*(1.D0+DTANH(GAM*SCR))
      RETURN
      END
C     ----------------------------------------
      FUNCTION SHADFNCONE(X1,X2,X3,XS,YS,ZS)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z,M)
      COMMON /XY/ PI,CD,CR,RE
C     ----------------------------------------
C     CONICAL SHADOWN FUNCTION OBTAINED FROM
C     DG COOK, 2001 MASTERS THESIS
C     ----------------------------------------
      RDS=696342.D0
      X=X1
      Y=X2
      Z=X3
      XP=XS-X
      YP=YS-Y
      ZP=ZS-Z
      RS=VMN(XS,YS,ZS)
      R=VMN(X,Y,Z)
      P=VMN(XP,YP,ZP)
      ET=DACOS(-(X*XP+Y*YP+Z*ZP)/(R*P))
      ET2=ET*ET
      RHOE=DASIN(RE/R)
      RHOS=DASIN(RDS/P)
      RHOE2=RHOE*RHOE
      RHOS2=RHOS*RHOS
      RH1=DABS(RHOS2-RHOE2)
      RH2=DABS(RHOS2+RHOE2)

      IF(ET.GE.(RHOS+RHOE)) GO TO 1
      IF(ET.LE.(RHOE-RHOS)) GO TO 2
      IF(ET.GT.(RHOE-RHOS).AND.ET.LT.(RHOE+RHOS)) GO TO 3

      GO TO 4
 1    SHADFNCONE=1.D0
      GO TO 4
 2    SHADFNCONE=0.D0
      GO TO 4
 3    AL=0.5D0*(RHOS+RHOE+ET)
      Q=2.D0*DSQRT(AL*(AL-ET)*(AL-RHOS)*(AL-RHOE))/ET
      IF(RH1.LE.ET2) T=RHOS2*DASIN(Q/RHOS)+RHOE2*DASIN(Q/RHOE)-Q*ET
      IF(RH1.GT.ET2) T=RHOS2*(PI-DASIN(Q/RHOS))+RHOE2*DASIN(Q/RHOE)-Q*ET
      SHADFNCONE=1.D0-T/(PI*RHOS2)

 4    CONTINUE
      RETURN
      END      
C     ----------------------------------------
      FUNCTION VMN(X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      VMN=DSQRT(X*X+Y*Y+Z*Z)
      RETURN
      END
C     ----------------------------------------
      SUBROUTINE QUAD(A,B,C,R1,R2)
c     Quadratic Equation Roots
      IMPLICIT DOUBLE PRECISION(A-H,O-Z,M)
      D=B*B-4.D0*A*C
      IF(D.GT.0.D0) THEN
         R1=-B+DSQRT(D)/(2.D0*A)
         R2=-B-DSQRT(D)/(2.D0*A)
      ELSEIF(D.EQ.0.D0) THEN
         R1=-B/(2.D0*A)
         R2=R1
      ENDIF         
      RETURN
      END
C     ----------------------------------------      

      SUBROUTINE LUNEPR(DPV,DDJ,SCM,DVEL,IO)
c     propagator for half a revolution
      IMPLICIT DOUBLE PRECISION(A-H,O-Z,M)
!      PARAMETER(REM = 384400)
      COMMON /XY/ PI,CD,CR,RE
      COMMON /TINT/ AJ2,AJ3,AJ4,AJ5,AK2,AKS2,AKM2
      COMMON /UV/ X1,X2,X3,X1DOT,X2DOT,X3DOT,XS,YS,ZS,XM,YM,ZM
      common /rarp/ rpin,rain
      DIMENSION DPV(1,6),DVEL(2)
      DIMENSION U(4),US(4),Y(10),Q(4),Z(10),SP(3),X(3),P(3),XD(3)
      DIMENSION UU(4),PS(4),PL(4),SL(3),PQ(4),O(3) !

      ISTEP = 180
      DE = 1.D0*cr
      E=0.0D0
      AT=0.0
      T=AT
      TT=0.0D0
            
      X1 = DPV(1,1)
      X2 = DPV(1,2)
      X3 = DPV(1,3)
      X1DOT = DPV(1,4)
      X2DOT = DPV(1,5)
      X3DOT = DPV(1,6)
      CALL LUNARPV(DDJ)
      VAM=VMN(XM,YM,ZM)      
      call cpvoe(a,ec,ai,omg,somg,ea,am)

      ra = a*(1.d0+ec)
      rp = a*(1.d0-ec)      
      REM = VAM
      VPERI1=VPERI(rpin,rain)
      VPERI2=VPERI(rpin,REM)
      DVEL(1)=VPERI2-VPERI1

C     DELTAV AFTER SOLAR PERTURBATION
      VPERI1=VPERI(rp,ra)
      VPERI2=VPERI(rp,REM)
      DVEL(2)=VPERI2-VPERI1
      
      a = (rem+rp)*0.5d0
      ec = (rem-rp)/(rem+rp)
      
      call coepv(a,ec,ai,omg,somg,am)
      
      RE2=RE*RE
      RE3=RE*RE2
      RE4=RE2*RE2
      RE5=RE3*RE2

      
      R=VMN(X1,X2,X3)
      R2=R*R
      RI=1.0D0/R
      RI2=RI*RI
      RI3=RI/R2
      RI4=RI2*RI2
      RI5=RI3*RI2
      RI6=RI3*RI3
      R4=R2*R2
      RDOT=VMN(X1DOT,X2DOT,X3DOT)
      RDOT2=RDOT*RDOT
      X32=X3*X3
      X3R=X3/R

      X32R2=X3R*X3R
      XR3=X3R*X32R2
      XR4=X32R2*X32R2
      XR5=XR3*X32R2

      X(1)=X1
      X(2)=X2
      X(3)=X3

C     *** SUN & MOON  POSITION VECTORS ***

      CALL SOLARNPV(DDJ)
      CALL LUNARPV(DDJ)

C     *** SUN AND MOON GRAVITY POTENTIAL
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS3=AS2*AS
      AS4=AS3*AS
      AS41=1.0D0/AS4
      AS5=AS3*AS2
      AS51=1.0D0/AS5
      X1XS=X1*XS
      X2YS=X2*YS
      X3ZS=X3*ZS
      VS1=X1XS+X2YS+X3ZS
      VS2=VS1*VS1
      
      VAM=VMN(XM,YM,ZM)
      VAM2=VAM*VAM
      VAM3=VAM2*VAM
      VAM4=VAM3*VAM
      VAM41=1.0D0/VAM4
      VAM5=VAM3*VAM2
      VAM51=1.0D0/VAM5
      X1XM=X1*XM
      X2YM=X2*YM
      X3ZM=X3*ZM
      VM1=X1XM+X2YM+X3ZM
      VM2=VM1*VM1

C     *** PERTURBING POTENTIAL ***
      V=0.5d0*(-AKS2*AS51*(3.0D0*VS2-R2*AS2)
     *     -(AKM2*VAM51*(3.0D0*VM2-R2*VAM2))
     *     +AK2*RE2*AJ2*RI3*(3.D0*X32R2-1.D0)
     *     +AK2*RE3*AJ3*RI4*(5.D0*XR3-3.D0*X3R)
     *     +AK2*RE4*AJ4*0.25D0*RI5*(35.D0*XR4-30.D0*X32R2-3.D0)
     *     +AK2*RE5*AJ5*0.25D0*RI6*(63.D0*XR5-70.D0*XR3+15.D0*X3R))

      W=AK2/R-RDOT2/2.0D0-V 
      W2=W/2.0D0 
      W=DSQRT(W2)

!     Modification
      wk=dsqrt((AK2/R-RDOT2/2.0D0)/2.d0)
      W3=W2*W
      gam=w/wk

      IF (X1) 12,11,11
 11   U(4)=0.0D0
      RPX1=R+X1 
      U(1)=.5D0*RPX1-U(4)*U(4)
      U(1)=DSQRT(U(1)) 
      U(2)=(X2*U(1)+X3*U(4))/RPX1 
      U(3)=(X3*U(1)-X2*U(4))/RPX1 
      GO TO 13
 12   U(3)=0.0D0
      RMX1=R-X1 
      U(2)=DSQRT(.5D0*RMX1)
      U(1)=(X2*U(2)+X3*U(3))/RMX1 
      U(4)=(X3*U(2)-X2*U(3))/RMX1 
 13   CNST=1.0D0/(4.0D0*W)
      US(1)=(U(1)*X1DOT+U(2)*X2DOT+U(3)*X3DOT)*CNST 
      US(2)=(-U(2)*X1DOT+U(1)*X2DOT+U(4)*X3DOT)*CNST
      US(3)=(-U(3)*X1DOT-U(4)*X2DOT+U(1)*X3DOT)*CNST
      US(4)=(U(4)*X1DOT-U(3)*X2DOT+U(2)*X3DOT)*CNST
      
      TOW=T+(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W 

      NT=0

      Y(1)=TOW
      Y(10)=W
      DO 21 I=1,4
         Y(I+1)=U(I)
         Y(I+5)=2.0*US(I) 
 21   CONTINUE

C     *** 1/2W^2 *** CHANGED FROM 1/4W2 TO 1/2W2 (100)
      CH1=1.0D0/(2.0D0*W2)
      CH3=2.D0/W

C     *** 1/8W^3 ***
      CH2=1.0D0/(8.0D0*W3)
      EEE=EAIN*CR+E
      SIE=DSIN(E/2.0D0)
      CSE=DCOS(E/2.0D0)

C     CHECKS DURING INTEGRATION
      CHECK=Y(5)*Y(6)-Y(4)*Y(7)+Y(3)*Y(8)-Y(2)*Y(9)
      DDJ1=DDJ
      X(1)=X1
      X(2)=X2
      X(3)=X3
C     ---------------------------- 
      CALL CPVOE(A,EC,AI,OMGA,SOMGA,EA,AM)

      de=de*gam
      
      DO 122 KK = 1,ISTEP
C     -----------------------------
C     R-K-G 4(5) INTEGRATOR LOOP 
C     -----------------------------
      
 60   CONTINUE

      XD(1)= X1DOT
      XD(2)= X2DOT
      XD(3)= X3DOT

C     SOLAR GRAVITY TERMS
      SP(1)=XS
      SP(2)=YS
      SP(3)=ZS
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS3=AS2*AS
      AS4=AS3*AS
      AS41=1.0D0/AS4
      AS5=AS3*AS2
      AS51=1.0D0/AS5
      VS=X1*XS+X2*YS+X3*ZS
      VS2=VS*VS
      
      CNTS1=3.D0*VS2-3.d0*R2*AS2
      CNTS2=6.D0*R*VS
      CNTS3=CNTS1+2.D0*R2*AS2

C     LUNAR GRAVITY TERMS
      SL(1)=XM
      SL(2)=YM
      SL(3)=ZM
      AM=VMN(XM,YM,ZM)

      AM2=AM*AM
      AM3=AM2*AM
      AM4=AM3*AM
      AM41=1.0D0/AM4
      AM5=AM3*AM2
      AM51=1.0D0/AM5
      VM=X1*XM+X2*YM+X3*ZM
      VM2=VM*VM
      
      CNTM1=3.D0*VM2-3.d0*R2*AM2
      CNTM2=6.D0*R*VM
      CNTM3=CNTM1+2.D0*R2*AM2

C     SOLAR GRAVITY POTENTIAL, PS      
      PS(1)=U(1)*SP(1)+U(2)*SP(2)+U(3)*SP(3)
      PS(2)=-U(2)*SP(1)+U(1)*SP(2)+U(4)*SP(3) 
      PS(3)=-U(3)*SP(1)-U(4)*SP(2)+U(1)*SP(3) 
      PS(4)=U(4)*SP(1)-U(3)*SP(2)+U(2)*SP(3)

C     LUNAR GRAVITY POTENTIAL, PL
      PL(1)=U(1)*SL(1)+U(2)*SL(2)+U(3)*SL(3)
      PL(2)=-U(2)*SL(1)+U(1)*SL(2)+U(4)*SL(3) 
      PL(3)=-U(3)*SL(1)-U(4)*SL(2)+U(1)*SL(3) 
      PL(4)=U(4)*SL(1)-U(3)*SL(2)+U(2)*SL(3)

C     J2, J3, J4 POTENTIAL TERMS      
!     ORIGINAL
      CNTJ1=3.D0*X3R            !
      CNTJ2=-6.D0*X32R2         !
      CNTK1=7.5D0*X32R2-1.5D0   !
      CNTK2=6.D0*X3R-15.D0*XR3  !
      CNTL1=(17.5D0*XR3-7.5D0*X3R) !
      CNTL2=-35.D0*XR4+22.5D0*X32R2-1.5D0 !
      CNTN1=39.375D0*XR4-26.25D0*X32R2+1.875D0 !
      CNTN2=-78.75D0*XR5+70.D0*XR3-11.25D0*X3R
      
      CNTJ=AK2*RE2*AJ2*RI3      !
      CNTK=AK2*RE3*AJ3*RI4      !
      CNTL=AK2*RE4*AJ4*RI5      !
      CNTN=AK2*RE5*AJ5*RI6      !
      
      UU(1)=U(3)
      UU(2)=U(4)
      UU(3)=U(1)
      UU(4)=U(2)
      
      RD2=R/2.D0
      
      C1=1.0D0/(2.0D0*W2) 
      C2=2.0D0/W

C     ----------------------------------      
C     DALFA/DE,DW/DE,DTAU/DE,DW/DE
C     ----------------------------------

      Z(10) = 0.D0
      DO 52 I=1,4
         O(I)=2.D0*Z(10)*US(I)/W
 52   CONTINUE

      Z(1)=(AK2+2*AKS2*AS51*R*CNTS3+2*AKM2*AM51*R*CNTM3
     *     +R*(CNTJ*0.5D0*(3.D0*X32R2-1.D0)+
     *     CNTK*0.5D0*(5.D0*XR3-3.D0*X3R)+
     *     CNTL*3.D0*0.125D0*(35.D0*XR4-30.D0*X32R2+3.D0)+
     *     CNTN*4.D0*0.125D0*(63.D0*XR5-70.D0*XR3+15.D0*X3R)))
     *     /(8.D0*W*W2)-(Z(10)*2.D0/(W2))*(U(1)*US(1)+U(2)*US(2)
     *     +U(3)*US(3)+U(4)*US(4))
      

!     Original
      DO 51 I=1,4 

      Z(I+1)=((-AKS2*AS51*(CNTS1*U(I)+CNTS2*PS(I)))*0.5D0
     *        -(AKM2*AM51*(CNTM1*U(I)+CNTM2*PL(I)))*0.5D0+
     *        (CNTJ*((1.D0+CNTJ2)*U(I)+CNTJ1*UU(I)))+
     *        CNTK*(U(I)*CNTK2+UU(I)*CNTK1)+
     *        CNTL*(U(I)*CNTL2+UU(I)*CNTL1)+
     *        CNTN*(U(I)*CNTN2+UU(I)*CNTN1))
      
      Z(I+5)=((-AKS2*AS51*(CNTS1*U(I)+CNTS2*PS(I)))*0.5D0
     *     -(AKM2*AM51*(CNTM1*U(I)+CNTM2*PL(I)))*0.5D0+
     *     (CNTJ*((1.D0+CNTJ2)*U(I)+CNTJ1*UU(I)))+
     *     CNTK*(U(I)*CNTK2+UU(I)*CNTK1)+
     *     CNTL*(U(I)*CNTL2+UU(I)*CNTL1)+
     *     CNTN*(U(I)*CNTN2+UU(I)*CNTN1))
         
 51   CONTINUE

      DO I=1,4
         Z(I+1)=(Z(I+1)*0.25/W2+O(I))*SIE
         Z(I+5)=-(Z(I+5)*0.25/W2+O(I))*CSE
      END DO

C     ----------------------------------

      S=RKGIL(N,Y,Z,E,DE,NT)

      CSE=DCOS(E/2.D0)
      SIE=DSIN(E/2.D0)

      DO 71 I=1,4 
         U(I)=Y(I+1)*CSE+Y(I+5)*SIE
         US(I)=-.5D0*Y(I+1)*SIE+.5D0*Y(I+5)*CSE
 71   CONTINUE
      
      X1=U(1)*U(1)-U(2)*U(2)-U(3)*U(3)+U(4)*U(4)
      X2=2.0D0*(U(1)*U(2)-U(3)*U(4))
      X3=2.0D0*(U(1)*U(3)+U(2)*U(4))
      R=U(1)*U(1)+U(2)*U(2)+U(3)*U(3)+U(4)*U(4)

      R2=R*R 
      R4=R2*R2
      CST1=4.0D0*W/R
      X32=X3*X3
      X3R=X3/R
      X32R2=X3R*X3R
      XR3=X3R*X32R2
      XR4=X32R2*X32R2
      XR5=XR3*X32R2

      X1DOT=CST1*(U(1)*US(1)-U(2)*US(2)-U(3)*US(3)+U(4)*US(4))
      X2DOT=CST1*(U(2)*US(1)+U(1)*US(2)-U(4)*US(3)-U(3)*US(4))
      X3DOT=CST1*(U(3)*US(1)+U(4)*US(2)+U(1)*US(3)+U(2)*US(4))
      RDOT=VMN(X1DOT,X2DOT,X3DOT)
      RDOT2=RDOT*RDOT
      X(1)=X1
      X(2)=X2
      X(3)=X3

      RI=1.0D0/R
      RI2=RI*RI
      RI3=RI*RI2
      RI4=RI2*RI2
      RI5=RI3*RI2
      RI6=RI3*RI3
      W=Y(10)
      W2=W*W
      T=Y(1)-(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W
      CALL CTHMS(T,DH,DM,SECONDS)
      CALL CSMHJD(DH,DM,SECONDS,DJT)
      DDJ1=DJT+DDJ

      CALL SOLARNPV(DDJ1)
      CALL LUNARPV(DDJ1)

C     *** SUN AND MOON POSITION VECTORS ***  
      AS=VMN(XS,YS,ZS)
      AS2=AS*AS
      AS5=AS**5
      AS51=1.0D0/AS5
      X1XS=X1*XS
      X2YS=X2*YS
      X3ZS=X3*ZS
      VS1=X1XS+X2YS+X3ZS
      VS2=VS1**2
      RRS2=(R*AS)**2

      AM=VMN(XM,YM,ZM)
      AM2=AM*AM
      AM5=AM**5
      AM51=1.0D0/AM5
      X1XM=X1*XM
      X2YM=X2*YM
      X3ZM=X3*ZM
      VM1=X1XM+X2YM+X3ZM
      VM2=VM1**2
      RRM2=(R*AM)**2

      V=0.5d0*(-AKS2*AS51*(3.0D0*VS2-R2*AS2)
     *     -(AKM2*AM51*(3.0D0*VM2-R2*AM2))
     *     +AK2*RE2*AJ2*RI3*(3.D0*X32R2-1.D0)
     *     +AK2*RE3*AJ3*RI4*(5.D0*XR3-3.D0*X3R)
     *     +AK2*RE4*AJ4*0.25D0*RI5*(35.D0*XR4-30.D0*X32R2-3.D0))
     *     +AK2*RE5*AJ5*0.25D0*RI6*(63.D0*XR5-70.D0*XR3+15.D0*X3R)

      
      IF (S-1.) 60,60,70
  70  CONTINUE

C     ---------------------------------------
C     *** TIME ***
      T=Y(1)-(U(1)*US(1)+U(2)*US(2)+U(3)*US(3)+U(4)*US(4))/W

      CALL CTHMS(T,DH,DM,SECONDS)
      CALL CSMHJD(DH,DM,SECONDS,DJT)
      DDJ1=DJT+DDJ
      CALL SOLARNPV(DDJ1)
      CALL LUNARPV(DDJ1)
      EAR=EA
      CALL CPVOE(A,EC,AI,OMGAD,SOMGAD,EA,AM)      
      
C     STOPPING CONDITIONS
      DXM=(X1-XM)**2
      DYM=(X2-YM)**2
      DZM=(X3-ZM)**2

      SCM=DSQRT(DXM+DYM+DZM)
      IF(SCM.LT.50000.D0) THEN
         IO = 1
         GO TO 422
      ENDIF
      IO = 0
 122  CONTINUE
 422  CONTINUE
      RETURN
      END

      FUNCTION VPERI(RP,RA)
C     COMPUTE VELOCITY AT PERIGEE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      AK2=3.986004415D05
      VP2=2.D0*AK2*(1.D0/RP-1.D0/(RA+RP))
      VPERI=DSQRT(VP2)
      RETURN
      END
      
