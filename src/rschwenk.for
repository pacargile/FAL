      PROGRAM RTIOSCHWENKE
C     READS PACKED BINARY VERSION OF SCHWENKE'S TIO LINELIST
C     THROWS AWAY LEVEL INFORMATION IF LINOUT < 0
      PARAMETER (kw=99)
      REAL*4 XJTIO(269300)
      REAL*8 ETIO(269300,5),STATETIO(269300,5)
      COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),
     1        WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,
     2      NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,
     3        DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
      REAL*8 LINDAT8(14)
      REAL*4 LINDAT4(28)
      EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
      REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
      REAL*8 LABEL,LABELP,OTHER1,OTHER2,LABELISO(5)
      CHARACTER*10 COTHER1,COTHER2
      EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
      INTEGER TYPE
      EQUIVALENCE (GF,G,CGF),(TYPE,NLAST)
      REAL*8 RESOLU,RATIO,RATIOLG,WLBEG,WLEND,RATIOLOG
      REAL*4 DECKJ(7,kw),XISO(5),X2ISO(5)
      REAL*4 TABLOG(32768)
      REAL*8 AIRSHIFT(60000)
      INTEGER*2 IELION,IELO,IGFLOG,IGR,IGS,IGW
      COMMON /IIIIIII/IWL,IELION,IELO,IGFLOG,IGR,IGS,IGW
      INTEGER*4 IIIIIII(4)
      EQUIVALENCE (IIIIIII(1),IWL)
      BYTE IIBYTE(16),ONEBYTE
      EQUIVALENCE (IIIIIII(1),IIBYTE(1))
C               46TiO 47TiO 48TiO 49TiO 50TiO
      DATA XISO/.0793,.0728,.7394,.0551,.0534/                                                              
      DATA X2ISO/-1.101,-1.138,-0.131,-1.259,-1.272/
      DATA LABELISO/2H46,2H47,2H48,2H49,2H50/
C
      DO 1 I=1,32768
    1 TABLOG(I)=10.**((I-16384)*.001)
      IF(IFPRED.NE.1)CALL TABVACAIR(AIRSHIFT)
      RATIOLOG=LOG(1.D0+1.D0/2000000.D0)

      OPEN(UNIT=11,file='fort.11',STATUS='OLD',READONLY,
     1FORM='UNFORMATTED',RECORDTYPE='FIXED',BLOCKSIZE=8000,
     1RECORDSIZE=4,ACCESS='DIRECT')
      OPEN(UNIT=12,FILE='fort.12',STATUS='OLD',FORM='UNFORMATTED'
     $     ,ACCESS='APPEND')
      OPEN(UNIT=14,FILE='fort.14',STATUS='OLD',FORM='UNFORMATTED'
     $     ,ACCESS='APPEND')
      OPEN(UNIT=93,FILE="fort.93",STATUS='OLD',
     $     FORM='UNFORMATTED',POSITION='REWIND')
      READ(93)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,
     1     WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT
      CLOSE(UNIT=93)

      IF(DEXP(IXWLBEG*RATIOLG).LT.WLBEG)IXWLBEG=IXWLBEG+1
C
C      OPEN(UNIT=48,FILE='fort.48',STATUS='OLD',READONLY,
C     1FORM='UNFORMATTED',RECORDTYPE='FIXED',BLOCKSIZE=269300,
C     2RECORDSIZE=5)
      READ(48)ETIO,XJTIO,STATETIO
      CLOSE(UNIT=48)

      N14=0
      NBLO=0                                                                    
      NBUP=0                                                                    
      OTHER1(1)=(8H        )                                                    
      OTHER1(2)=(2H  )                                                          
      OTHER2(1)=(8H        )                                                    
      OTHER2(2)=(2H  )                                                          
      LABEL(2)=(2H  )
      LABELP(2)=(2H  )
      REF=(4HSCHW)                                                                  
      ISO1=16
      X1=0.
      CODE=822.
      DWL=0.
      DGFLOG=0.
      DGAMMAR=0.
      DGAMMAS=0.
      DGAMMAW=0.
      DWLISO=0.
      ISOSHIFT=0
C
      ISTART=DLOG(WLBEG-1.)/RATIOLOG+.5
      ISTOP=DLOG(WLEND+1.)/RATIOLOG+.5
      N=0
      READ(11,REC=1)IWL1
      IF(IWL1.GT.ISTOP)GO TO 21
C     FIND NUMBER OF LINES
      LIMITBLUE=1
      LIMITRED=50000000
    8 NEWLIMIT=(LIMITRED+LIMITBLUE)/2
      READ(11,REC=NEWLIMIT,ERR=9)IWL
      LIMITBLUE=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.EQ.1)GO TO 11
      GO TO 8
    9 LIMITRED=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.EQ.1)GO TO 11
      GO TO 8
   11 LENGTHFILE=LIMITBLUE
      WLVAC=EXP(IWL1*RATIOLOG)
      PRINT 3334,IWL1,WLVAC
 3334 FORMAT(' FIRST LINE IS        1','  IWL',I10,'   WL',F12.4)
      READ(11,REC=LENGTHFILE)IWL
      WLVAC=EXP(IWL*RATIOLOG)
      PRINT 3335,LENGTHFILE,IWL,WLVAC
 3335 FORMAT(' LAST LINE IS ',I9,'  IWL',I10,'   WL',F12.4)
      IF(IWL.LT.ISTART)GO TO 21
C     FIND THE FIRST LINE AFTER ISTART
      LIMITBLUE=1
      LIMITRED=LENGTHFILE
   12 NEWLIMIT=(LIMITRED+LIMITBLUE)/2
      PRINT 3333,LIMITBLUE,NEWLIMIT,LIMITRED
 3333 FORMAT(3I10)
      READ(11,REC=NEWLIMIT)IWL
C     IF COMPUTER REQUIRES BYTE ROTATION
C      DO 17 I=1,4,2
C      ONEBYTE=IIBYTE(I)
C      IIBYTE(I)=IIBYTE(I+1)
C   17 IIBYTE(I+1)=ONEBYTE
C      ONEBYTE=IIBYTE(1)
C      IIBYTE(1)=IIBYTE(3)
C      IIBYTE(3)=ONEBYTE
C      ONEBYTE=IIBYTE(2)
C      IIBYTE(2)=IIBYTE(4)
C      IIBYTE(4)=ONEBYTE
      IF(IWL.LT.ISTART)GO TO 13
      LIMITRED=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.LE.1)GO TO 14
      GO TO 12
   13 LIMITBLUE=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.LE.1)GO TO 14
      GO TO 12
   14 ISTART=NEWLIMIT
      PRINT 3333,LIMITBLUE,LIMITRED,NEWLIMIT
      WRITE(6,6)ISTART
    6 FORMAT(I10,14H IS FIRST LINE)
      DO 20 ILINE=ISTART,LENGTHFILE
      READ(11,REC=ILINE)IIIIIII
C     IF COMPUTER REQUIRES BYTE ROTATION
C      DO 18 I=1,16,2
C      ONEBYTE=IIBYTE(I)
C      IIBYTE(I)=IIBYTE(I+1)
C   18 IIBYTE(I+1)=ONEBYTE
C      ONEBYTE=IIBYTE(1)
C      IIBYTE(1)=IIBYTE(3)
C      IIBYTE(3)=ONEBYTE
C      ONEBYTE=IIBYTE(2)
C      IIBYTE(2)=IIBYTE(4)
C      IIBYTE(4)=ONEBYTE
      IF(IWL.GT.ISTOP)GO TO 21
      ISO=ABS(IELION)-8949
      ISO2=ISO+45      
      X2=X2ISO(ISO)
C     NELIONNEW=ABS(IELION)/10
C     NELION=NELIONOLD(NELIONNEW)
C     IF(NELION.EQ.0)GO TO 20
      NELION=366
      WLVAC=EXP(IWL*RATIOLOG)
      KWL=WLVAC*10.D0+.5D0
      WL=WLVAC+AIRSHIFT(KWL)
      IF(IFVAC.NE.1)WLVAC=WL
      ELO=TABLOG(IELO)
      IXWL=DLOG(WLVAC)/RATIOLG+.5D0
      NBUFF=IXWL-IXWLBEG+1
      FREQ=2.99792458D17/WLVAC
      CONGF=.01502D0*TABLOG(IGFLOG)/FREQ*XISO(ISO)
c     Reduce loggf by 3.0
c      CONGF=CONGF/3.0
c
      GFLOG=((IGFLOG-16384D0)*.001D0)
c      print 2,WLVAC,ELO,GFLOG
c 2    FORMAT(F12.6,F12.3,F10.4)
      FRQ4PI=FREQ*12.5664D0
      KGW=IGW
      KGS=IGS
      LEVELLO=KGS*10+MOD(ABS(KGW),10)
      LEVELUP=KGW/10+LEVELLO
C     GAMMAS=0
C     LOG GAMMAW=-7
      IGS=1
      IGW=9384
      GAMRF=TABLOG(IGR)/FRQ4PI
      GAMSF=TABLOG(IGS)/FRQ4PI
      GAMWF=TABLOG(IGW)/FRQ4PI
      WRITE(6,*)'NBUFF',NBUFF
      WRITE(12)NBUFF,CONGF,NELION,ELO,GAMRF,GAMSF,GAMWF
      NLINES=NLINES+1
      IF(N.EQ.0)WRITE(6,19)WLVAC
   19 FORMAT(F12.4)
      N=N+1
      IF(LINOUT.LT.0)GO TO 20
      E=ETIO(LEVELLO,ISO)
      EP=ETIO(LEVELUP,ISO)
      XJ=XJTIO(LEVELLO)
      XJP=XJTIO(LEVELUP)
      LABEL(1)=STATETIO(LEVELLO,ISO)
      LABELP(1)=STATETIO(LEVELUP,ISO)
      LABELP(2)=LABELISO(ISO)
      GFLOG=(IGFLOG-16384D0)*.001D0
      GF=TABLOG(IGFLOG)
      GR=(IGR-16384D0)*.001D0
C     GS=-16.383
      GS=-9.99
      GW=-7.
      WRITE(14)LINDAT8,LINDAT4
      N14=N14+1
   20 CONTINUE
   21 N=N-1
      print *,n14
      WRITE(6,22)N
   22 FORMAT(I10,13H IS LAST LINE)
      WRITE(6,19)WLVAC
   25 WRITE(6,26)NLINES
   26 FORMAT(I10,25H LINES WRITTEN ON TAPE 12)
C
      OPEN(UNIT=93,FILE="fort.93",STATUS='OLD',
     $     FORM='UNFORMATTED',POSITION='REWIND')
      WRITE(93)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,
     1     WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT
      CLOSE(UNIT=93)
C
      CALL EXIT
      END
      SUBROUTINE TABVACAIR(AIRSHIFT)
      REAL*8 AIRSHIFT(60000)
      REAL*8 WLVAC,VACAIR
      DO 1 IWL=1,1999
    1 AIRSHIFT(IWL)=0.
      DO 2 IWL=2000,60000
      WLVAC=IWL*.1D0
    2 AIRSHIFT(IWL)=VACAIR(WLVAC)-WLVAC
      RETURN
      END
      FUNCTION VACAIR(W)
      IMPLICIT REAL*8 (A-H,O-Z)
C     W IS VACUUM WAVELENGTH IN NM
      WAVEN=1.D7/W
      VACAIR=W/(1.0000834213D0+
     1 2406030.D0/(1.30D10-WAVEN**2)+15997.D0/(3.89D9-WAVEN**2))
C    1(1.000064328+2949810./(1.46E10-WAVEN**2)+25540./(4.1E9-WAVEN**2))         
      RETURN
      END
