module f_wrapper

use iso_c_binding, only: c_double, c_int, c_char, c_null_char

implicit none

contains

function c_to_f_string(s) result(str)
  use iso_c_binding
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  integer i, nchars
  i = 1
  do
     if (s(i) == c_null_char) exit
     i = i + 1
  end do
  nchars = i - 1  ! Exclude null character from Fortran string
  allocate(character(len=nchars) :: str)
  str = transfer(s(1:nchars), str)
end function c_to_f_string

subroutine readoutspecbin(&
  s, NWLi, NLINESi,&
  wli, qmu1i, qmu2i, &
  WLin,DWLin,GFLOGin,DGFLOGin,CODEin,Ein,XJin,&
  LABELin) bind(c, name='readoutspecbin')
  use iso_c_binding, only: c_double, c_int, c_char, c_null_char
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str

  integer(c_int), intent(in), value :: NWLi
  integer(c_int), intent(in), value :: NLINESi

  real(c_double), intent(out) :: wli(NWLi)
  real(c_double), intent(out) :: qmu1i(NWLi)
  real(c_double), intent(out) :: qmu2i(NWLi)

  real(c_double), intent(out) :: WLin(NLINESi)
  real(c_double), intent(out) :: DWLin(NLINESi)
  real(c_double), intent(out) :: GFLOGin(NLINESi)
  real(c_double), intent(out) :: DGFLOGin(NLINESi)
  real(c_double), intent(out) :: CODEin(NLINESi)
  real(c_double), intent(out) :: Ein(NLINESi)
  real(c_double), intent(out) :: XJin(NLINESi)
  character(kind=c_char,len=1),  intent(out) :: LABELin(NLINESi,8)
  ! real(c_double), intent(out) :: EPin(NLINESi)
  ! real(c_double), intent(out) :: XJPin(NLINESi)
  ! character(kind=c_char,len=1),   intent(out) :: LABELPin(NLINESi)
  ! real(c_double), intent(out) :: GRin(NLINESi)
  ! real(c_double), intent(out) :: DGAMMARin(NLINESi)
  ! real(c_double), intent(out) :: GSin(NLINESi)
  ! real(c_double), intent(out) :: DGAMMASin(NLINESi)
  ! real(c_double), intent(out) :: GWin(NLINESi)
  ! real(c_double), intent(out) :: DGAMMAWin(NLINESi)
  ! real(c_double), intent(out) :: WAVENOin(NLINESi)
  ! character(kind=c_char,len=1),   intent(out) :: REFin(NLINESi)
  ! integer(c_int),    intent(out) :: NBLOin(NLINESi)
  ! integer(c_int),    intent(out) :: NBUPin(NLINESi)
  ! integer(c_int),    intent(out) :: ISO1in(NLINESi)
  ! real(c_double), intent(out) :: X1in(NLINESi)
  ! integer(c_int),    intent(out) :: ISO2in(NLINESi)
  ! real(c_double), intent(out) :: X2in(NLINESi)
  ! character(kind=c_char,len=1),   intent(out) :: OTHER1in(NLINESi)
  ! character(kind=c_char,len=1),   intent(out) :: OTHER2in(NLINESi)
  ! integer(c_int),    intent(out) :: ISOSHIFTin(NLINESi)
  ! integer(c_int),    intent(out) :: NELIONin(NLINESi)
  ! real(c_double), intent(out) :: RESIDin(NLINESi)

  COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
        WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
        NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
        DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3,ALINEC


  REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,WLEDGE,RATIO,SWL
  REAL*8 QMU(40),XMU(20),NEDGE,IFSURF,NMU
  REAL*8 NAV,NMU2,FREQTOWAVE

  REAL*8 LINDAT8(14)
  REAL*4 LINDAT(28)
  EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
  REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN,RESID
  REAL*8 REF,LABEL,LABELP,OTHER1,OTHER2
  REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW,WAVENO
  REAL*4 X1,X2,ELO,GF,GS,GR,GW
  REAL*4 NBLO, NBUP, ISO1, ISO2, ISOSHIFT
  REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3,NELION
  REAL*4 ALINEC

  character(len=8) :: SLABEL
  ! character(len=1) :: REF,OTHER1(2),OTHER2(2)

  INTEGER IWL, NWL, I, NLINESO

  open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
  read(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,WLEDGE
  
  NAV=IFSURF/10
  IFSURF=IFSURF-NAV*10
  IF(NAV.EQ.0.)NAV=1
  IF(IFSURF.EQ.3)NMU=1
  NMU2=NMU+NMU
  RATIO=1.D0+1.D0/RESOLU
  DO 6 IWL=1,NWL
     SWL=WBEGIN*RATIO**(IWL-1)
     FREQTOWAVE=2.99792458D17/SWL**2
     wli(IWL) = SWL
     READ(1)(QMU(NLINESi),I=1,NMU2)
     qmu1i(IWL) = QMU(1)*FREQTOWAVE
     qmu2i(IWL) = QMU(2)*FREQTOWAVE
6 CONTINUE
  READ(1)NLINESO
140     FORMAT(F11.4,F7.4,2F7.3,F8.2,F12.3,F5.1,1X,A8,A2,&
  F12.3,F5.1,1X,A8,A2,6F6.2,F11.3,&
  1X,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,I6,I4,2X,f8.4)

  DO 9 I=1,NLINESO
     READ(1)LINDAT8,LINDAT
     IF(I.EQ.1)WRITE(6,140)WL,DWL,GFLOG,DGFLOG,CODE,E,XJ,LABEL,&
     EP,XJP,LABELP,GR,DGAMMAR,GS,DGAMMAS,GW,DGAMMAW,WAVENO,&
     REF,NBLO,NBUP,ISO1,X1,ISO2,X2,OTHER1,OTHER2,ISOSHIFT,&
     NELION,resid

     resid=center/concen
     WLin(I) = WL
     DWLin(I) = DWL
     GFLOGin(I) = GFLOG
     DGFLOGin(I) = DGFLOG
     CODEin(I) = CODE
     Ein(I) = E
     XJin(I) = XJ
     WRITE(SLABEL,'(A8)') LABEL(1)
     print *, SLABEL
     LABELin(I,1) = SLABEL
     ! EPin(I) = EP
     ! XJPin(I) = XJP
     ! LABELPin(I) = LABELP(1)
     ! GRin(I) = GR
     ! DGAMMARin(I) = DGAMMAR
     ! GSin(I) = GS
     ! DGAMMASin(I) = DGAMMAS
     ! GWin(I) = GW
     ! DGAMMAWin(I) = DGAMMAW
     ! WAVENOin(I) = WAVENO
     ! REFin(I) = REF
     ! NBLOin(I) = NBLO
     ! NBUPin(I) = NBUP
     ! ISO1in(I) = ISO1
     ! X1in(I) = X1
     ! ISO2in(I) = ISO2
     ! X2in(I) = X2
     ! OTHER1in(I) = OTHER1(1)
     ! OTHER2in(I) = OTHER2(1)
     ! ISOSHIFTin(I) = ISOSHIFT
     ! NELIONin(I) = NELION
     ! RESIDin(I) = RESID
9 CONTINUE

  CLOSE(UNIT=1)

end subroutine readoutspecbin

! wrap more functions here
! ...

end module
