      PROGRAM RPREDICT
c     Fior-1jun04
C     eguale a RLOWLINES.FOR salvo che skippa le rigke con ielion.gt.0
C     si usa con le files gf100
      PARAMETER (kw=99)
      COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),
     1        WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,
     2      NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,
     3        DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
      REAL*8 LINDAT8(14)
      REAL*4 LINDAT4(28)
      EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)	
      REAL*8 RESOLU,RATIO,RATIOLG,WLBEG,WLEND,RATIOLOG
      REAL*8 WL,E,EP,WLVAC,START,STOP
      REAL*8 CENTER,CONCEN,LABEL,LABELP,OTHER1,OTHER2
      REAL*4 DECKJ(7,kw)
      REAL*4 TABLOG(32768)
      real*8 vacair
      INTEGER*2 IELION,IELO,IGFLOG,IGR,IGS,IGW
      COMMON /IIIIIII/IWL,IELION,IELO,IGFLOG,IGR,IGS,IGW
      INTEGER*4 IIIIIII(4)
      EQUIVALENCE (IIIIIII(1),IWL)
      DIMENSION NELIONOLD(1005)
      DIMENSION NELIONOLDA(209),NELIONOLDB(286),NELIONOLDC(95)
      DIMENSION NELIONOLDD( 95),NELIONOLDE( 95),NELIONOLDF(60)
      DIMENSION NELIONOLDG(165)
      EQUIVALENCE (NELIONOLD(  1),NELIONOLDA(1))
      EQUIVALENCE (NELIONOLD(210),NELIONOLDB(1))
      EQUIVALENCE (NELIONOLD(496),NELIONOLDC(1))
      EQUIVALENCE (NELIONOLD(591),NELIONOLDD(1))
      EQUIVALENCE (NELIONOLD(686),NELIONOLDE(1))
      EQUIVALENCE (NELIONOLD(781),NELIONOLDF(1))
      EQUIVALENCE (NELIONOLD(841),NELIONOLDG(1))
      DATA NELIONOLDA/
     1   1,  2,
     2   7,  8,  9,
     3  13, 14, 15, 16,
     4  19, 20, 21, 22, 23,
     5  25, 26, 27, 28, 29, 30,
     6  31, 32, 33, 34, 35, 36,0,
     7  37, 38, 39, 40, 41, 42,0,0,
     8  43, 44, 45, 46, 47, 48,0,0,0,
     9  49, 50, 51, 52, 53, 54,0,0,0,0,
     A  55, 56, 57, 58, 59, 60,0,0,0,0,0,
     1  61, 62, 63, 64, 65, 66,0,0,0,0,0,0,
     2  67, 68, 69, 70, 71, 72,0,0,0,0,0,0,0,
     3  73, 74, 75, 76, 77, 78,0,0,0,0,0,0,0,0,
     4  79, 80, 81, 82, 83, 84,0,0,0,0,0,0,0,0,0,
     5  85, 86, 87, 88, 89, 90,0,0,0,0,0,0,0,0,0,0,
     6  91, 92, 93, 94, 95, 96,0,0,0,0,0,0,0,0,0,0,0,
     7  97, 98, 99,100,101,  0,0,0,0,0,0,0,0,0,0,0,0,0,
     8 103,104,105,106,107,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     9 109,110,111,112,113,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/
      DATA NELIONOLDB/
     A 115,116,117,118,119,120,299,359,419,479, 11*0,
     1 121,122,123,124,125,126,305,365,425,485, 12*0,
     2 127,128,129,130,131,132,311,371,431,491, 13*0,
     3 133,134,135,136,137,138,317,377,437,497, 14*0,
     4 139,140,141,142,143,144,323,383,443,503, 15*0,
     5 145,146,147,148,149,150,329,389,449,509, 16*0,
     6 151,152,153,154,155,156,335,395,455,515, 17*0,
     7 157,158,159,160,161,162,341,401,461,521, 18*0,
     8 163,164,165,166,167,168,347,407,467,527, 19*0,
     9 169,170,171,             27*0,
     A 175,176,177,             28*0/
      DATA NELIONOLDC/
     1 181,182,183,  0,  0,
     2 187,188,189,  0,  0,
     3 193,194,195,  0,  0,
     4 199,200,201,  0,  0,
     5 205,206,207,  0,  0,
     6 211,212,213,  0,  0,
     7 217,218,219,  0,  0,
     8 223,224,225,  0,  0,
     9 229,230,231,  0,  0,
     A 235,236,237,  0,  0,
     1 241,242,243,  0,  0,
     2 247,248,249,  0,  0,
     3 253,254,255,  0,  0,
     4 259,260,261,  0,  0,
     5 265,266,267,  0,  0,
     6 271,272,273,  0,  0,
     7 277,278,279,  0,  0,
     8 283,284,285,  0,  0,
     9 289,290,291,  0,  0/
      DATA NELIONOLDD/
     A 295,296,297,  0,  0,
     1 301,302,303,  0,  0,
     2 307,308,309,  0,  0,
     3 313,314,315,  0,  0,
     4 319,320,321,  0,  0,
     5 325,326,327,  0,  0,
     6 331,332,333,  0,  0,
     7 337,338,339,  0,  0,
     8 343,344,345,  0,  0,
     9 349,350,351,  0,  0,
     A 355,356,357,  0,  0,
     1 361,362,363,  0,  0,
     2 367,368,369,  0,  0,
     3 373,374,375,  0,  0,
     4 379,380,381,  0,  0,
     5 385,386,387,  0,  0,
     6 391,392,393,  0,  0,
     7 397,398,399,  0,  0,
     8 403,404,405,  0,  0/
      DATA NELIONOLDE/
     9 409,410,411,  0,  0,
     A 415,416,417,  0,  0,
     1 421,422,423,  0,  0,
     2 427,428,429,  0,  0,
     3 433,434,435,  0,  0,
     4 439,440,441,  0,  0,
     5 445,446,447,  0,  0,
     6 451,452,453,  0,  0,
     7 457,458,459,  0,  0,
     8 463,464,465,  0,  0,
     9 469,470,471,  0,  0,
     A 475,476,477,  0,  0,
     1 481,482,483,  0,  0,
     2 487,488,489,  0,  0,
     3 493,494,495,  0,  0,
     4 499,500,501,  0,  0,
     5 505,506,507,  0,  0,
     6 511,512,513,  0,  0,
     7 517,518,519,  0,  0/
      DATA NELIONOLDF/
     8 523,524,525,  0,  0,
     9 529,530,531,  0,  0,
     A 535,536,537,  0,  0,
     1 541,542,543,  0,  0,
     2 547,548,549,  0,  0,
     3 553,554,555,  0,  0,
     4 559,560,561,  0,  0,
     5 565,566,567,  0,  0,
     6 571,572,573,  0,  0,
     7 577,578,579,  0,  0,
     8 583,584,585,  0,  0,
     9 589,590,591,  0,  0/
c     H2     841   240
c     HeH    842
c     LiH    843   378
c     BeH    844   384
c     BH     845   390
c     CH     846   246
c     NH     847   252
c     OH     848   258
c     HF     849   396
c     NaH    850
c     MgH    851   300
c     AlH    852   306
c     SiH    853   312
c     PH     854   402
c     HS     855   336
c     HCl    856   408
c     KH     857
c     CaH    858   342
c     ScH    859   414
c     TiH    860   420
c     VH     861   426
c     CrH    862   432
c     MnH    863   438
c     FeH    864   444
c     CoH    865   558
c     NiH    866   564
c     CuH    867   570
c     C2     868   264
c     CN     869   270
c     CO     870   276
c     CF     871
c     SiC    872
c     CP     873
c     CS     874
c     N2     875   282
c     NO     876   288
c     NF     877
c     SiN    878
c     PN     879
c     NS     880
c     LiO    881
c     BeO    882   492
c     BO     883   498
c     O2     884   294
c     FO     885
c     NaO    886
c     MgO    887   318
c     AlO    888   324
c     SiO    889   330
c     PO     890   504
c     SO     891   348
c     ClO    892   510
c     CaO    893   354
c     ScO    894   360
c     TiO    895   366
c     VO     896   372
c     CrO    897   516
c     MnO    898   522
c     FeO    899   528
c     CoO    900   576
c     NiO    901   582
c     CuO    902   588
c     GeO    903
c     SrO    904
c     YO     905
c     ZrO    906
c     NbO    907
c     Si2    908
c     SiS    909
c     S2     910
c     TiS    911
c     ZrS    912
c     H2+    913
c     HeH+   914
c     LiH+   915
c     CH+    916   450
c     NH+    917   456
c     OH+    918   462
c     HF+    919
c     NeH+   920
c     MgH+   921   468
c     AlH+   922   474
c     SiH+   923   480
c     PH+    924
c     SH+    925
c     HCl+   926
c     CaH+   927   486
c     He2+   928
c     C2+    929
c     CN+    930
c     CO+    931   594
c     N2+    932
c     NO+    933
c     NS+    934
c     O2+    935
c     SiO+   936
c     PO+    937
c     SO+    938
c     S2+    939
c     H2O    940   534
c     CO2    941   540
c     CH2    942   546
c     C2H    943
c     C2N    944
c     C3     945   552
c     O3     946
c     NO2    947
c     N2O    948
c     NH2    949
c     HCO    950
c     HCN    951
c     HNO    952
c     SiC2   953
c     NaOH   954
c     MgOH   955
c     AlOH   956
c     KOH    957
c     CaOH   958
c     AlOF   959
c     AlOCl  960
c     Al2O   961
c     SH2    962
c     CaF2   963
c     CaCl2  964
c     COS    965
c     SiO2   966
c     SO2    967
c     TiO2   968
c     VO2    969
c     NH3    970
c     CH3    971
c     C2H2   972
c     C3H    973
c     C2N2   974
c     CH4    975
c     H-     976
c     Li-    977
c     C-     978
c     O-     979
c     F-     980
c     Na-    981
c     Al-    982
c     Si-    983
c     P-     984
c     S-     985
c     Cl-    986
c     K-     987
c     Sc-    988
c     Ti-    989
c     V-     990
c     Cr-    991
c     Fe-    992
c     Co-    993
c     Ni-    994
c     Cu-    995
c     C2-    996
c     CH-    997
c     CN-    998
c     CO-    999
c     N2-   1000
c     NO-   1001
c     OH-   1002
c     O2-   1003
c     S2-   1004
c     SH-   1005
C
      DATA NELIONOLDG/
     1 240,  0,378,384,390,246,252,258,396,  0,
     2 300,306,312,402,336,408,  0,342,414,420,
     3 426,432,438,444,558,564,570,264,270,276,
     4   0,  0,  0,  0,282,288,  0,  0,  0,  0,
     5   0,492,498,294,  0,  0,318,324,330,504,
     6 348,510,354,360,366,372,516,522,528,576,
     7 582,588,  0,  0,  0,  0,  0,  0,  0,  0,
     8   0,  0,  0,  0,  0,450,456,462,  0,  0,
     9 468,474,480,  0,  0,  0,486,  0,  0,  0,
     A 594,  0,  0,  0,  0,  0,  0,  0,  0,534,
     1 540,546,  0,  0,552,  0,  0,  0,  0,  0,
     2   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     3   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     4   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     5   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     6   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     7   0,  0,  0,  0,  0/
      DO 1 I=1,32768
    1 TABLOG(I)=10.**((I-16384)*.001)
      RATIOLOG=LOG(1.D0+1.D0/2000000.D0)

      OPEN(UNIT=11,file='fort.11',STATUS='OLD',READONLY,
     1FORM='UNFORMATTED',RECORDTYPE='FIXED',BLOCKSIZE=8000,
     1RECORDSIZE=4,ACCESS='DIRECT')

C     CHANGED THE FOLLOWING OPEN STATEMENTS
      OPEN(UNIT=12,FILE='fort.12',STATUS='OLD',FORM='UNFORMATTED'
     $     ,ACCESS='APPEND')
      OPEN(UNIT=14,FILE='fort.14',STATUS='OLD',FORM='UNFORMATTED'
     $     ,ACCESS='APPEND')
      OPEN(UNIT=93,FILE="fort.93",STATUS='OLD',
     $     FORM='UNFORMATTED',POSITION='REWIND')
      READ(93)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,
     1     WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT
      CLOSE(UNIT=93)

      IXWLBEG=DLOG(WLBEG)/RATIOLG
      IF(DEXP(IXWLBEG*RATIOLG).LT.WLBEG)IXWLBEG=IXWLBEG+1
C
      REF=2HK 
      NBLO=0
      NBUP=0
      ISO1=0
      X1=0.
      ISO2=0
      X2=0.
      XJ=0.
      XJP=0.
      OTHER1(1)=(8H        )
      OTHER1(2)=(2H  ) 
      OTHER2(1)=(8H        )
      OTHER2(2)=(2H  ) 
      LABEL(1)=(8H        )
      LABEL(2)=(2H  ) 
      LABELP(1)=(8H        )
      LABELP(2)=(2H  ) 
      ISTART=DLOG(WLBEG-1.)/RATIOLOG+.5
      ISTOP=DLOG(WLEND+1.)/RATIOLOG+.5
      N=0
      READ(11,REC=1)IWL1
      IF(IWL1.GT.ISTOP)GO TO 21
C     FIND NUMBER OF LINES
      LIMITBLUE=1
      LIMITRED=500000000
    8 NEWLIMIT=(LIMITRED+LIMITBLUE)/2
      READ(11,REC=NEWLIMIT,ERR=9)IWL
      LIMITBLUE=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.EQ.1)GO TO 11
      GO TO 8
    9 LIMITRED=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.EQ.1)GO TO 11
      GO TO 8
   11 LENGTHFILE=LIMITBLUE
      PRINT 3334,IWL1
 3334 FORMAT(' FIRST LINE IS        1','  IWL',I10)
      READ(11,REC=LENGTHFILE)IWL
      PRINT 3335,LENGTHFILE,IWL
 3335 FORMAT(' LAST LINE IS ',I9,'  IWL',I10)
      IF(IWL.LT.ISTART)GO TO 21
C     FIND THE FIRST LINE AFTER ISTART
      LIMITBLUE=1
      LIMITRED=LENGTHFILE
   12 NEWLIMIT=(LIMITRED+LIMITBLUE)/2
      PRINT 3333,LIMITBLUE,NEWLIMIT,LIMITRED
 3333 FORMAT(3I10)
      READ(11,REC=NEWLIMIT)IWL
      IF(IWL.LT.ISTART)GO TO 13
      LIMITRED=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.LE.1)GO TO 14
      GO TO 12
   13 LIMITBLUE=NEWLIMIT
      IF(LIMITRED-LIMITBLUE.LE.1)GO TO 14
      GO TO 12
   14 NSTART=NEWLIMIT
      PRINT 3333,LIMITBLUE,LIMITRED,NEWLIMIT
      WRITE(6,6)NSTART
    6 FORMAT(I10,14H IS FIRST LINE)

      WRITE(6,*)NSTART,LENGTHFILE,NLINES

      DO 20 N=NSTART,LENGTHFILE
      READ(11,rec=n)IIIIIII
      IF(IWL.GT.ISTOP)GO TO 21
C     REMOVE POSITIVE LINES WHICH ARE ALSO IN OTHER FILES
C     THE NEW PRED LINE FILE DOES NOT HAVE THIS CONVENTION
C      IF(IELION.GT.0)GO TO 20
      NELIONNEW=ABS(IELION)/10
      NELION=NELIONOLD(NELIONNEW)
      IF(NELION.EQ.0)GO TO 20
      NELEM=INT(NELION/6)+1
      IZ=NELION-(NELEM-1)*6
      ICHARGE=IZ-1
      CHARGE=FLOAT(ICHARGE)
      CODE=FLOAT(NELEM)+CHARGE*0.01
      WLVAC=EXP(IWL*RATIOLOG)
      WL=WLVAC
      IF(IFVAC.EQ.0)WL=VACAIR(WL)
      GF=TABLOG(IGFLOG)
      GFLOG=LOG10(GF)
      FREQ=2.99792458D17/WLVAC
C     CONGF=.026538/1.77245*GF/FREQ
      CONGF=.01502D0*GF/FREQ
      ELO=TABLOG(IELO)
      E=ELO
      EP=E+1.D7/WLVAC
      GAMMAR=TABLOG(IGR)
      IF(LOG10(GAMMAR).EQ.0.0)WRITE(6,*)GAMMAR
      GAMMAS=TABLOG(IGS)
      GAMMAW=TABLOG(IGW)
      GR = LOG10(real(GAMMAR))
      GS = LOG10(real(GAMMAS))
      GW = LOG10(real(GAMMAW))
      IF(GR.EQ.0.0)WRITE(6,*)GAMMAR
      IF(GR.EQ.0.)THEN
      GAMMAR=2.223D13/WLVAC**2
      GR=LOG10(GAMMAR)
      ENDIF
      FRQ4PI=FREQ*12.5664D0
      GAMRF=TABLOG(IGR)/FRQ4PI
      GAMSF=TABLOG(IGS)/FRQ4PI
      GAMWF=TABLOG(IGW)/FRQ4PI
      IF(GAMRF.EQ.0.0)WRITE(6,*)GAMMAR
      alpha=0.
      IXWL=DLOG(WLVAC)/RATIOLG+.5D0
      if(ifvac.eq.0)ixwl=dlog(wl)/ratiolg+.5D0
      NBUFF=IXWL-IXWLBEG+1
      WRITE(12)NBUFF,CONGF,NELION,ELO,GAMRF,GAMSF,GAMWF,alpha
      WRITE(14)LINDAT8,LINDAT4
      NLINES=NLINES+1
   20 CONTINUE

   21 N=N-1
      WRITE(6,22)N
   22 FORMAT(I10,13H IS LAST LINE)
   25 WRITE(6,26)NLINES
   26 FORMAT(I10,25H LINES WRITTEN ON TAPE 12)
      OPEN(UNIT=93,FILE="fort.93",STATUS='OLD',
     $     FORM='UNFORMATTED',POSITION='REWIND')
      WRITE(93)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,
     1     WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT
      CLOSE(UNIT=93)
      CALL EXIT
      END
      FUNCTION VACAIR(W)
c      IMPLICIT REAL*8 (A-H,O-Z)
C     W IS VACUUM WAVELENGTH IN NM
      real*8 vacair,w,waven
      WAVEN=1.D7/W
      VACAIR=W/(1.0000834213D0+
     1 2406030.D0/(1.30D10-WAVEN**2)+15997.D0/(3.89D9-WAVEN**2))
      RETURN
      END
