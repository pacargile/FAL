      PROGRAM ROTATE
c     revised 4nov14  constants given D exponents
c     revised 18jan05
c     default radius is 100 pixels instead of 50
c     differential rotation put in using solar expression
c     VROT(LAT2)=(462-75*SIN(LAT)**2)-50*SIN(LAT)**4)*2*PI*RSUN/1.E9/1.E5 KM/S
c     VROT(0)= 2.020 km/s   equator
c     VROT(1)= 1.474 km/s   pole
c     VROT(LAT)/VROT(0) = (1-75./462.*SIN(LAT)**2)-50./462.*SIN(LAT)**4)
c     from Libbrecht, K.G. and Morrow, C.A. The solar rotation. pp. 479-500 in
c     The Solar Interior and Atmosphere, eds. A.N. Cox, W.C. Livingston, and
c     M. Matthews, Tucson: University of Arizona Press, 1991.
c     All input rotation velocities are equitorial.  Differential rotation
c     velocities are specified by making the velocity negative.
c     Thus 2 produces the approximate solar rotation, -2 produces the
c     approximate solar differential rotation, and -2.020
c     matches the solar differential rotation expression above.
c
      parameter (npiece=2000,npiece2=npiece*2,npiece3=npiece*3)
      PARAMETER (kw=99)
      COMMON /HROT/H(500),HROT(npiece3)
      COMMON /WT/MUNWT(10000),IVNWT(10000),WTNWT(10000)
      DIMENSION CONT(npiece2)
      DIMENSION WTMU(100)
      DIMENSION XMU100(100),INT100(102)
      EQUIVALENCE (INT100(101),FLUX),(INT100(102),CONTIN)
      DIMENSION R(25),INTEN(26),XX(26)
      REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,XMU(20),WLEDGE(377)
      REAL*8 QMU(40),Q2(2)
C      REAL*8 LINDAT8(14)
C      REAL*4 LINDAT(28)
CCCCCC
      COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),
     1        WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,
     2      NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,
     3        DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
     4 ,ALINEC(kw)
      REAL*8 LINDAT8(14)
      REAL*4 LINDAT(28)
      EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
      REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
      REAL*8 LABEL,LABELP,OTHER1,OTHER2
C     CHANGE THE FOLLOWING TO DOUBLES?
      REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW
      REAL*4 REF,X1,X2,ELO,GF,GS,GR,GW
      REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
      REAL*4 ALINEC
CCCCCC
      REAL INT100
      REAL INTEN
C      REAL*8 WEND
      REAL*8 WEND,RATIO
      DIMENSION VROT(25)
      EQUIVALENCE (DUMMY,IDUMMY)
C      CHARACTER*5 ROTNAME(25)
      CHARACTER*4 ROTNAME(9)
      CHARACTER*25 TESTX
      DIMENSION APLOT(101)
      DATA APLOT/101*1H /
C      DATA ROTNAME/'ROT1','ROT2','ROT3','ROT4','ROT5','ROT6','ROT7',
C     1     'ROT8','ROT9','ROT10','ROT11','ROT12','ROT13','ROT14',
C     $     'ROT15',
C     2     'ROT16','ROT17','ROT18','ROT19','ROT20','ROT21','ROT22',
C     3     'ROT23','ROT24','ROT25'/
      DATA ROTNAME/'ROT1','ROT2','ROT3','ROT4','ROT5','ROT6','ROT7',
     1     'ROT8','ROT9'/
C      CHARACTER*9 RUNID
C      CALL GETARG(1,RUNID)
C     CALL BEGTIME
      WRITE(6,*)"STARTING ROTATE"
      DO 7 I=1,100
    7 XMU100(I)=FLOAT(I)*.01-.005
C       LINOUT=30000
      LINOUT=100000000
C       LINOUT=300
      WRITE(6,*)"OPENING FORT.19"
      OPEN(UNIT=19,FILE='fort.19',STATUS='OLD',
     $     FORM='UNFORMATTED', POSITION='REWIND')
      READ(5,*)TESTX
      WRITE(6,*)TESTX
      READ(5,1001)NROT,NRADIUS,(VROT(IROT),IROT=1,NROT)
 1001 FORMAT(I5,I5/(8F10.1))
      WRITE(6,*)NRADIUS,NROT
C     IF(NRADIUS.EQ.0)NRADIUS=50
      IF(NRADIUS.EQ.0)NRADIUS=100
      WRITE(6,1002)NROT,NRADIUS,(VROT(IROT),IROT=1,NROT)
 1002 FORMAT(18H1ROTATION          ,I3,I5,10F6.1/(10F6.1))
      OPEN(UNIT=1,FILE='fort.1',FORM='UNFORMATTED',STATUS='OLD',
     $     POSITION='REWIND')
C      REWIND 1
      READ(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,
     $     NEDGE,WLEDGE
C     IFSURF=3 FOR ROTATED SPECTRUM
      IFSURF=3
      WRITE(6,1010)TEFF,GLOG,TITLE
 1010 FORMAT(  5H TEFF,F7.0,7H   GRAV,F7.3/7H TITLE ,74A1)
      WRITE(6,1007)NMU,(XMU(IMU),IMU=1,NMU)
 1007 FORMAT(18H SURFACE INTENSITY,I3,10F6.3/10F6.3)
      RATIO=1.+1./RESOLU
      WEND=WBEGIN*RATIO**(NWL-1)
      VSTEP=299792.458D0/RESOLU
      WRITE(6,1005)WBEGIN,WEND,RESOLU,VSTEP
 1005 FORMAT(2F12.5,F12.1,F12.5)
      NMU2=NMU+NMU
C
      XX(1)=0.
      NM1=NMU+1
      DO 11 MU=1,NMU
      NN=NMU-MU+2
   11 XX(NN)=XMU(MU)
      CALL WTROT(0.,0.,0,NWT,WTMU,NRADIUS)
      WRITE(6,777)WTMU
  777 FORMAT(1P10E12.3)
      INTEN(1)=0.
      DO 19 IWL=1,NWL
      READ(1)(QMU(I),I=1,NMU2)
      FLUX=0.
      CONTIN=0.
      DO 13 MU=1,NMU
      NN=NMU-MU+2
   13 INTEN(NN)=QMU(MU+NMU)
      IDUMMY=MAP1(XX,INTEN,NM1,XMU100,INT100,100)
      DO 14 I=1,100
   14 CONTIN=CONTIN+INT100(I)*WTMU(I)
      DO 15 MU=1,NMU
      NN=NMU-MU+2
   15 INTEN(NN)=QMU(MU)
      IDUMMY=MAP1(XX,INTEN,NM1,XMU100,INT100,100)
      DO 16 I=1,100
   16 FLUX=FLUX+INT100(I)*WTMU(I)
      WRITE(19)INT100
   19 CONTINUE
      NMU=1
C
      DO 500 IROT=1,NROT
         OPEN(UNIT=9,FILE=ROTNAME(IROT),
     $        FORM='UNFORMATTED',STATUS='OLD',POSITION='REWIND')
         REWIND 19
         VEL=ABS(VROT(IROT))
         NV=VEL/VSTEP+1.5
         NAV=NV/5+1
         NAVWT=NAV
         ENDWT=0.
         IF(MOD(NAV,2).EQ.0)ENDWT=.5
         IF(MOD(NAV,2).EQ.0)NAV=NAV+1
         NAV100=500-NAV/2
         NAVNAV=NAV100+NAV-1
         WRITE(6,1011)VEL,NV
 1011    FORMAT(5H1VROT,F10.1,I5)
         WRITE(6,7171)
 7171    FORMAT(3HNAV,1H,6HNAV100,1H,6HNAVNAV)
         write(6,778)NAV,NAV100,NAVNAV
 778     FORMAT(10I10)
C        CHECKS FOR VEL = 0
         IF(VEL.EQ.0.)GO TO 50
C        IF WE HAVE A VROT, DO THE BROADENING
         CALL WTROT(VEL,VSTEP,NV,NWT,WTMU,NRADIUS)
         WRITE(6,1013)NWT
 1013    FORMAT(4H NWT,I6)
         DO 29 IWL=1,npiece3
 29         HROT(IWL)=0.
C     
         WRITE(6,1117)
 1117    FORMAT(1H1)
         WRITE(9)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,
     1        NEDGE,WLEDGE,VEL,NV
         REWIND 1
         DO 40 IWL=npiece+1,NWL+npiece,npiece
            MAX=MIN0(npiece2,NWL+npiece2-IWL+1)
            DO 30 J=npiece+1,MAX
               KWL=IWL+J-npiece2-1
               READ(19)INT100
               CONT(J)=CONTIN
               DO 25 I=1,NWT
                  MU=MUNWT(I)
                  IV=IVNWT(I)
                  W=WTNWT(I)*INT100(MU)
                  HROT(J-IV)=HROT(J-IV)+W
 25               HROT(J+IV)=HROT(J+IV)+W
 30         CONTINUE
            IF(IWL.EQ.npiece+1)GO TO 37
            DO 33 J=1,npiece
               QH=-(H(J+NAV100)+H(J+NAVNAV))*ENDWT
               DO 330 K=NAV100,NAVNAV
 330              QH=QH+H(J+K)
               Q2(1)=QH/FLOAT(NAVWT)
               Q2(2)=CONT(J)
               WRITE(9)Q2
               JWL=IWL+J-npiece2-1
               IF(JWL.GT.LINOUT)GO TO 33
               WAVE=WBEGIN*RATIO**(JWL-1)
               RESID=Q2(1)/Q2(2)
               IRESID=RESID*1000.+.5
               IPLOT=RESID*100.+1.5
               IPLOT=MAX0(1,MIN0(101,IPLOT))
               APLOT(IPLOT)=1HX
C                WRITE(6,2300)JWL,WAVE,IRESID,APLOT
 2300          FORMAT(1H ,I5,F11.4,I7,101A1)
               APLOT(IPLOT)=(1H )
 33         CONTINUE
 37         DO 34 J=1,npiece
 34            CONT(J)=CONT(J+npiece)
            DO 350 J=1,500
 350           H(J)=HROT(J+npiece-500)            
            DO 35 J=1,npiece2
 35            HROT(J)=HROT(J+npiece)
            DO 36 J=npiece2+1,npiece3
 36            HROT(J)=0.
            IF(KWL.LT.NWL)GO TO 40
            MAX=MIN0(npiece,NWL+npiece-IWL+1)
            DO 38 J=1,MAX
               QH=-(H(J+NAV100)+H(J+NAVNAV))*ENDWT
C     H SEEMS TO BE THE PROBLEM ON ODY, EDIT BELOW
C     WORKS ON STAMPEDE!!!!
               DO 380 K=NAV100,NAVNAV
 380              QH=QH+H(J+K)
C               DO 380 K=NAV100,NAVNAV
C 380              QH=QH+HROT(J+K-500)
               Q2(1)=QH/FLOAT(NAVWT)
               Q2(2)=CONT(J)
               WRITE(9)Q2
               JWL=IWL+J-npiece-1
               IF(JWL.GT.LINOUT)GO TO 38
               WAVE=WBEGIN*RATIO**(JWL-1)
               RESID=Q2(1)/Q2(2)
               IRESID=RESID*1000.+.5
               IPLOT=RESID*100.+1.5
               IPLOT=MAX0(1,MIN0(101,IPLOT))
               APLOT(IPLOT)=1HX
C               WRITE(6,2300)JWL,WAVE,IRESID,APLOT
               APLOT(IPLOT)=(1H )
 38         CONTINUE
 40      CONTINUE
         GO TO 400
 50      WRITE(9)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,
     1        NEDGE,WLEDGE,VEL,NV
         WRITE(6,1117)
         DO 55 IWL=1,NWL
            READ(19)INT100
            Q2(1)=FLUX
            Q2(2)=CONTIN
            WRITE(9)Q2
            IF(IWL.GT.LINOUT)GO TO 55
            WAVE=WBEGIN*RATIO**(IWL-1)
            RESID=FLUX/CONTIN
            IRESID=RESID*1000.+.5
            IPLOT=RESID*100.+1.5
            IPLOT=MAX0(1,MIN0(101,IPLOT))
            APLOT(IPLOT)=1HX
C             WRITE(6,2300)IWL,WAVE,IRESID,APLOT
            APLOT(IPLOT)=(1H )
 55      CONTINUE
 400     REWIND 1
         READ(1)
         DO 42 I=1,NWL
 42         READ(1)
         READ(1)NLINES
         WRITE(9)NLINES
         WRITE(6,*)"NLINES FROM FORT.1:"
         WRITE(6,*)NLINES
         WRITE(6,*)"LINE INFO:"
         DO 41 I=1,NLINES
            READ(1)LINDAT8,LINDAT
            WRITE(9)LINDAT8,LINDAT
            WRITE(6,*)CODE,WL,DWL,GFLOG,DGFLOG,GW,DGAMMAW
C            IF(WL.GT.1561.11.AND.WL.LT.1561.12) THEN
C            WRITE(6,*)CODE,WL,DWL,GFLOG,DGFLOG,GW,DGAMMAW
C            END IF
 41      CONTINUE
 500     CLOSE(UNIT=9)
      CLOSE(UNIT=2,DISPOSE='DELETE')
      CALL EXIT
      END
      SUBROUTINE WTROT(VEL,VSTEP,NV,NWT,WTMU,NRAD)
      COMMON /WT/MUNWT(10000),IVNWT(10000),WTNWT(10000)
      DIMENSION WTMU(100)
      REAL*4 LAT
      DO 1 MU=1,100
 1       WTMU(MU)=0.
C     SYMMETRY ABOUT THE EQUATOR AND AXIS
C      NRAD=100
C      NRAD=50
      RADIUS=NRAD
C     CHOSEN TO YIELD HNU
      W=4./4./3.14159/RADIUS**2
C     CENTER
      CX=.5
      CY=.5
      N3=0
      DO 100 IX=1,NRAD
         DO 100 IY=1,NRAD
C     R IS THE PROJECTED RADIUS
            R=SQRT((IX-CX)**2+(IY-CY)**2)
            IF(R.GT.RADIUS)GO TO 100
            XMU=SQRT(RADIUS**2-R**2)/RADIUS
            MU=XMU*100.+.9999999
            IF(MU.EQ.0)GO TO 100
            WTMU(MU)=WTMU(MU)+W
            IF(VEL.EQ.0.)GO TO 100
C     RX IS THE RADIUS OF THE LATITUDE CIRCLE
            RX=SQRT(RADIUS**2-(IY-CY)**2)
C     VLAT IS THE VELOCITY AT THE LATITUDE
            VLAT=RX/RADIUS*ABS(VEL)
            IF(VEL.LT.0.)THEN
               LAT=ACOS(RX/RADIUS)
               VLAT=VLAT*(1.-75./462.*SIN(LAT)**2-50./462.*SIN(LAT)**4)
            ENDIF
C     VX IS THE PROJECTED VELOCITY
            VX=(IX-CX)/RX*VLAT
            IV=VX/VSTEP+.5
            IVMU=IV*1000+MU
            N3=N3+1
            MUNWT(N3)=IVMU
  100 CONTINUE
      IF(VEL.EQ.0.)RETURN
      CALL INTSORT(MUNWT,N3)
      ISAVE=-1
      NWT=0
C     POSITIVE AND NEGATIVE DOPPLER SHIFTS
      W=W*.5
      DO 300 I=1,N3
         IVMU=MUNWT(I)
         IF(IVMU.EQ.ISAVE)GO TO 310
         ISAVE=IVMU
         IV=IVMU/1000
         MU=IVMU-IV*1000
         NWT=NWT+1
         IF(NWT.GT.10000)STOP 'MORE THAN 10000 POINTS'
         MUNWT(NWT)=MU
         IVNWT(NWT)=IV
         WTNWT(NWT)=W
         GO TO 300
 310     WTNWT(NWT)=WTNWT(NWT)+W
  300 CONTINUE
      RETURN
      END
      FUNCTION MAP1(XOLD,FOLD,NOLD,XNEW,FNEW,NNEW)
      DIMENSION XOLD(1),FOLD(1),XNEW(1),FNEW(1)
      L=2
      LL=0
      DO 50 K=1,NNEW
   10 IF(XNEW(K).LT.XOLD(L))GO TO 20
      L=L+1
      IF(L.GT.NOLD)GO TO 30
      GO TO 10
   20 IF(L.EQ.LL)GO TO 50
      IF(L.EQ.2)GO TO 30
      IF(L.EQ.3)GO TO 30
      L1=L-1
      IF(L.GT.LL+1.OR.L.EQ.3)GO TO 21
      IF(L.GT.LL+1.OR.L.EQ.4)GO TO 21
      CBAC=CFOR
      BBAC=BFOR
      ABAC=AFOR
      IF(L.EQ.NOLD)GO TO 22
      GO TO 25
   21 L2=L-2
      D=(FOLD(L1)-FOLD(L2))/(XOLD(L1)-XOLD(L2))
      CBAC=FOLD(L)/((XOLD(L)-XOLD(L1))*(XOLD(L)-XOLD(L2)))+
     1(FOLD(L2)/(XOLD(L)-XOLD(L2))-FOLD(L1)/(XOLD(L)-XOLD(L1)))/
     2(XOLD(L1)-XOLD(L2))
      BBAC=D-(XOLD(L1)+XOLD(L2))*CBAC
      ABAC=FOLD(L2)-XOLD(L2)*D+XOLD(L1)*XOLD(L2)*CBAC
      IF(L.LT.NOLD)GO TO 25
   22 C=CBAC
      B=BBAC
      A=ABAC
      LL=L
      GO TO 50
   25 D=(FOLD(L)-FOLD(L1))/(XOLD(L)-XOLD(L1))
      CFOR=FOLD(L+1)/((XOLD(L+1)-XOLD(L))*(XOLD(L+1)-XOLD(L1)))+
     1(FOLD(L1)/(XOLD(L+1)-XOLD(L1))-FOLD(L)/(XOLD(L+1)-XOLD(L)))/
     2(XOLD(L)-XOLD(L1))
      BFOR=D-(XOLD(L)+XOLD(L1))*CFOR
      AFOR=FOLD(L1)-XOLD(L1)*D+XOLD(L)*XOLD(L1)*CFOR
      WT=0.
      IF(ABS(CFOR).NE.0.)WT=ABS(CFOR)/(ABS(CFOR)+ABS(CBAC))
      A=AFOR+WT*(ABAC-AFOR)
      B=BFOR+WT*(BBAC-BFOR)
      C=CFOR+WT*(CBAC-CFOR)
      LL=L
      GO TO 50
   30 IF(L.EQ.LL)GO TO 50
      L=AMIN0(NOLD,L)
      C=0.
      B=(FOLD(L)-FOLD(L-1))/(XOLD(L)-XOLD(L-1))
      A=FOLD(L)-XOLD(L)*B
      LL=L
   50 FNEW(K)=A+(B+C*XNEW(K))*XNEW(K)
      MAP1=LL-1
      RETURN
      END
      SUBROUTINE INTSORT(DATA,N)
      INTEGER X,Z,DATA(1)
      NTRY=0
      N1=2
   15 DO 1 J=N1,N
      Z=DATA(J)
      IF(J-2)1,2,3
    2 IF(Z-DATA(1))4,1,1
    4 DATA(2)=DATA(1)
      DATA(1)=Z
      GO TO 1
    3 K7=J-1
      IF(Z-DATA(K7))5,1,1
    5 LFST=1
      LAST=K7
    6 MID=(LFST+LAST)/2
      IF(Z-DATA(MID))7,8,9
    7 IF(MID-LAST)10,8,8
   10 LAST=MID
      GO TO 6
    8 NSTART=MID
      GO TO 11
    9 IF(LFST-MID)12,13,13
   12 LFST=MID
      GO TO 6
   13 NSTART=MID+1
   11 DO 14 I=NSTART,K7
      K9=J+NSTART-I
   14 DATA(K9)=DATA(K9-1)
      DATA(NSTART)=Z
    1 CONTINUE
      NTRY=NTRY+1
      DO 16 I=2,N
      IF(DATA(I)-DATA(I-1))17,16,16
   17 N1=I
      IF(NTRY-5)15,15,18
   16 CONTINUE
   18 RETURN
      END
