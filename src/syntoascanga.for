      PROGRAM SYNTOASCang
C     ******************************************************************
C     
C     Linux port by L. Sbordone, P. Bonifacio and F. Castelli
C     
C     -------------------------------------------------------
C     
C     - March 2004: Initial Linux port by L.S. and P.B.
C     
C     -------------------------------------------------------
C     
C     Please aknowledge the use of this code by citing:
C     
C     * Kurucz, R. 1993, ATLAS9 Stellar Atmosphere Programs and 2 km/s
C     grid. Kurucz CD-ROM No. 13. Cambridge, Mass.: 
C     Smithsonian Astrophysical
C     Observatory, 1993., 13
C     
C     * Sbordone, L., Bonifacio, P., Castelli, F., & Kurucz, R. L. 
C     2004a, Memorie
C     della Societa Astronomica Italiana Supplement, 5, 93
C     
C     --------------------------------------------------------
C     
C     For updates, documentation, utilities and needed files please 
C     refer to:
C     www.******.it
C     
C     ******************************************************************
c     implicit     REAL*8 (A-H,O-Z)
C     
C     
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     
C     DERIVED FROM SYNTOASC
C     PURPOSE: WRITE A SYNTHETIC SPECTRUM AND LINE DATA INTO TWO
C     SEPARATE FILES  
C     
C     P. BONIFACIO
C     
C     
C     
C     LIKE SYNTOASC BUT OUTPUT WAVELENGTHS ARE IN ANGSTROEMS
C     
C     
C     MARCH 1993
C     
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     
C     TAPE1=SPECTRUM INPUT
C     TAPE2=SPECTRUM OUTPUT IN ASCII 
C     TAPE3=LINE DATA
C     TAPE4=HEADER FILE FOR TEFF,GLOG ETC.
C     TAPE6=OUTPU
C********
C     revision for IFC 8.0, 23012004 LS
c     use ifport
C     *******
      PARAMETER (kw=99)
      INTEGER NEDGE
C      DIMENSION LINDAT(24)
C      EQUIVALENCE (LINDAT(1),WL)
      DIMENSION XMU(20),QMU(40),WLEDGE(200),TITLE(74)
      REAL*8 TEFF,GLOG,TITLE,WBEGIN,RESOLU,XMU,WLEDGE,RATIO
      REAL*8 QMU
      double precision WAVE,wend,wcen,vstep,resid
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
      DIMENSION APLOT(101)
      DATA APLOT/101*1H  /
c     interface
c     subroutine exit(status)
c     integer(4),optional,intent(in)::status
c     end subroutine
c     subroutine abort(string)
c     !ms$attributes alias:'abort_'::abort
c     character(len=*),optional,intent(in)::string
c     end subroutine
c     end interface
C     K=1.38054E-16

c     
      linout=1000000
c     OPEN(UNIT=2,STATUS='NEW')
c     OPEN(UNIT=3,STATUS='NEW')
c     OPEN(UNIT=4,STATUS='NEW')
      OPEN(UNIT=1,FILE='fort.1',STATUS='OLD',FORM='UNFORMATTED'
     $     , POSITION='REWIND')
      OPEN(UNIT=2,FILE='specfile.dat',STATUS='NEW',
     $     FORM='FORMATTED')
      OPEN(UNIT=3,FILE='lineinfo.dat',STATUS='NEW',
     $     FORM='FORMATTED')
      OPEN(UNIT=4,FILE='headinfo.dat',STATUS='NEW',
     $     FORM='FORMATTED')
      READ(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,
     1     WLEDGE
      WRITE(4,2233)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU
     &     ,NEDGE,WLEDGE
 2233 FORMAT(F10.1,F10.3/6HTITLE ,74A1/F10.3,F10.1,I10,I5,I5/
     1     10F8.4/10F8.4/I10/(5F16.5))

      RATIO=1.+1./RESOLU
      WEND=WBEGIN*RATIO**(NWL-1)
      WCEN=(WBEGIN+WEND)*.5
      VSTEP=2.99792458D5/RESOLU
      WRITE(4,*) "WBEGIN,WCEN,WEND,RESOLU,VSTEP,NWL"
      WRITE(4,*)WBEGIN,WCEN,WEND,RESOLU,VSTEP,NWL

      WRITE(6,1010)TEFF,GLOG,TITLE
 1010 FORMAT(  5H TEFF,F7.0,7H   GRAV,F7.3/7H TITLE ,74A1)
C      WRITE(6,1007)NMU,(XMU(IMU),IMU=1,NMU)
 1007 FORMAT(I4,20F6.3)
C     FOR FLUX SPECTRA NMU IS 1
      IF(IFSURF.EQ.3) NMU=1
      NMU1=NMU+1
      NMU2=NMU+NMU
      DO 70 IWL=1,NWL
         READ(1)(QMU(IMU),IMU=1,NMU2)
         IWLNMU=(IWL+9999)*NMU
c     IF(IWL.GT.LINOUT)GO TO 63                                         
         WAVE=WBEGIN*RATIO**(IWL-1)
         RESID=QMU(1)/QMU(NMU1)                                         
         IRESID=RESID*1000.+.5                                          
C         WRITE(6,2300)IWL,WAVE,IRESID,APLOT                           
 2300    FORMAT(1H ,I5,F11.4,I7,101A1)
c     
c     convert wavelengths to angstroms
c     
         FREQTOWAVE=2.99792458D17/WAVE**2
C         WAVEA=WAVE*10.0
         WRITE(2,2301)WAVE,QMU(1)*FREQTOWAVE,QMU(NMU1)*FREQTOWAVE
 2301    FORMAT(F17.8,2E13.6)
c     write(2,2301)wave,QMU(1),QMU(NMU1),resid
c     2301	format(f11.4,2E20.8,F13.6)
c     APLOT(IPLOT)=(1H )
 63      CONTINUE                                                       
c     68 CONTINUE                                                       
 70   CONTINUE                                                          
      READ(1)NLINES
      WRITE(4,2244)NLINES
 2244 FORMAT(1X,'NLINES= ',I10)
      IF(NLINES.EQ.0)GO TO 99
      DO 9 I=1,NLINES
         READ(1)LINDAT8,LINDAT
         resid=center/concen
         WRITE(3,140)WL,DWL,GFLOG,DGFLOG,CODE,E,XJ,LABEL,
     1        EP,XJP,LABELP,GR,DGAMMAR,GS,DGAMMAS,GW,DGAMMAW,WAVENO,
     2        REF,NBLO,NBUP,ISO1,X1,ISO2,X2,OTHER1,OTHER2,ISOSHIFT,
     $        NELION,resid
 140     FORMAT(F11.4,F7.4,2F7.3,F8.2,F12.3,F5.1,1X,A8,A2,
     1        F12.3,F5.1,1X,A8,A2,6F6.2,F11.3,
     2        1X,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,I6,I4,2X,f10.6)
C         WRITE(3,8)WL,GFLOG,XJ,E,XJP,EP,CODE,REF,resid
C 8       FORMAT(F10.4,F7.3,F5.1,F12.3,F5.1,F12.3,F9.2,2x,A4,1x,f8.4)
    9 CONTINUE
 99   CONTINUE
      CLOSE(UNIT=1)
      CLOSE(UNIT=2)
      CLOSE(UNIT=3)
      CLOSE(UNIT=4)
      CALL EXIT
      END
