module f_wrapper

use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_PTR, C_LOC, c_float, c_long

! implicit none

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

!** Convert a Fortran string to a C string
function f_to_c_string(f_string) result(c_string)
    character(len=*), intent(in) :: f_string
    character(len=1, kind=c_char) :: c_string(len_trim(f_string)+1)
    type(C_PTR) :: str
    integer :: N, i

    N = len_trim(f_string)
    do i = 1, N
        c_string(i) = f_string(i:i)
    end do
    c_string(n + 1) = c_null_char
    str = transfer(c_string,str)
end function f_to_c_string

subroutine readoutspecbin(&
  s, NWLi, NLINESi,&
  wli, qmu1i, qmu2i, &
  WLin,DWLin,GFLOGin,DGFLOGin,CODEin,Ein,XJin,&
  LABELin,EPin,XJPin,LABELPin,&
  GRin,DGAMMARin,GSin,DGAMMASin,GWin,DGAMMAWin,&
  WAVENOin,REFin,NBLOin,NBUPin,ISO1in,X1in,ISO2in,X2in,&
  OTHER1in,OTHER2in,ISOSHIFTin,NELIONin,RESIDin) bind(c, name='readoutspecbin')
  use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str
  character(len=25) :: SLABEL

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
  character(kind=c_char,len=1),  intent(out) :: LABELin(10,NLINESi)
  real(c_double), intent(out) :: EPin(NLINESi)
  real(c_double), intent(out) :: XJPin(NLINESi)
  character(kind=c_char,len=1),   intent(out) :: LABELPin(10,NLINESi)
  real(c_double), intent(out) :: GRin(NLINESi)
  real(c_double), intent(out) :: DGAMMARin(NLINESi)
  real(c_double), intent(out) :: GSin(NLINESi)
  real(c_double), intent(out) :: DGAMMASin(NLINESi)
  real(c_double), intent(out) :: GWin(NLINESi)
  real(c_double), intent(out) :: DGAMMAWin(NLINESi)
  real(c_double), intent(out) :: WAVENOin(NLINESi)
  character(kind=c_char,len=1),   intent(out) :: REFin(5,NLINESi)
  integer(c_long), intent(out) :: NBLOin(NLINESi)
  integer(c_long), intent(out) :: NBUPin(NLINESi)
  integer(c_long), intent(out) :: ISO1in(NLINESi)
  real(c_double), intent(out) :: X1in(NLINESi)
  integer(c_long), intent(out) :: ISO2in(NLINESi)
  real(c_double), intent(out) :: X2in(NLINESi)
  character(kind=c_char,len=1),   intent(out) :: OTHER1in(10,NLINESi)
  character(kind=c_char,len=1),   intent(out) :: OTHER2in(10,NLINESi)
  integer(c_long), intent(out) :: ISOSHIFTin(NLINESi)
  integer(c_long), intent(out) :: NELIONin(NLINESi)
  real(c_double), intent(out) :: RESIDin(NLINESi)

  REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,WLEDGE,RATIO,SWL
  REAL*8 QMU(40),XMU(20),NEDGE,IFSURF,NMU
  REAL*8 NAV,NMU2,FREQTOWAVE

  COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
          WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
        NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
          DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3&
   ,ALINEC(99)
  REAL*8 LINDAT8(14)
  REAL*4 LINDAT(28)
  EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
  REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
  REAL*8 LABEL,LABELP,OTHER1,OTHER2
  REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW
  REAL*4 REF,X1,X2,ELO,GF,GS,GR,GW
  REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
  REAL*4 ALINEC

  ! REAL*8 LINDAT8(14)
  ! REAL*4 LINDAT(28)
  ! EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
  ! REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN,RESID
  ! REAL*8 REF,LABEL,LABELP,OTHER1,OTHER2
  ! REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW,WAVENO
  ! REAL*4 X1,X2,ELO,GF,GS,GR,GW
  ! REAL*4 NBLO, NBUP, ISO1, ISO2, ISOSHIFT
  ! REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3,NELION
  ! REAL*4 ALINEC

  INTEGER IWL, NWL, I, INMU, NLINESO, J

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
     READ(1)(QMU(INMU),INMU=1,NMU2)
     qmu1i(IWL) = QMU(1)*FREQTOWAVE
     qmu2i(IWL) = QMU(2)*FREQTOWAVE
6 CONTINUE
  READ(1)NLINESO
140 FORMAT(F11.4,F7.4,2F7.3,F8.2,F12.3,F5.1,1X,A8,A2,&
  F12.3,F5.1,1X,A8,A2,6F6.2,F11.3,&
  1X,A4,I2,I2,I3,F6.3,I3,F6.3,A8,A2,A8,A2,I6,I4,2X,f8.4)

141 FORMAT(F11.4,F8.2,1X,A8,A2,1X,A8,A2)

  DO I=1,NLINESO
     READ(1)LINDAT8,LINDAT

     ! IF(I.EQ.6)WRITE(6,141)WL,CODE,LABEL,LABELP

     ! IF(I.EQ.6)WRITE(6,140)WL,DWL,GFLOG,DGFLOG,CODE,E,XJ,LABEL,&
     ! EP,XJP,LABELP,GR,DGAMMAR,GS,DGAMMAS,GW,DGAMMAW,WAVENO,&
     ! REF,NBLO,NBUP,ISO1,X1,ISO2,X2,OTHER1,OTHER2,ISOSHIFT,&
     ! NELION,resid

     resid=center/concen
     WLin(I) = WL
     DWLin(I) = DWL
     GFLOGin(I) = GFLOG
     DGFLOGin(I) = DGFLOG
     CODEin(I) = CODE
     ! WRITE(SLABEL,'(A8)') CODE
     ! SLABEL = SLABEL//c_null_char
     ! DO J=1,8
     !   CODEin(J,I) = SLABEL(J:J)
     ! END DO
     Ein(I) = E
     XJin(I) = XJ
     WRITE(SLABEL,'(X1,A8,A2)') LABEL
     SLABEL = SLABEL//c_null_char
     DO J=1,10
       LABELin(J,I) = SLABEL(J:J)
     END DO
     EPin(I) = EP
     XJPin(I) = XJP
     WRITE(SLABEL,'(X1,A8,A2)') LABELP
     SLABEL = SLABEL//c_null_char
     DO J=1,10
       LABELPin(J,I) = SLABEL(J:J)
     END DO
     GRin(I) = GR
     DGAMMARin(I) = DGAMMAR
     GSin(I) = GS
     DGAMMASin(I) = DGAMMAS
     GWin(I) = GW
     DGAMMAWin(I) = DGAMMAW
     WAVENOin(I) = WAVENO
     WRITE(SLABEL,'(X1,A4)') REF
     SLABEL = SLABEL//c_null_char
     DO J=1,5
       REFin(J,I) = SLABEL(J:J)
     END DO
     NBLOin(I) = INT(NBLO)
     NBUPin(I) = INT(NBUP)
     ISO1in(I) = INT(ISO1)
     X1in(I) = X1
     ISO2in(I) = INT(ISO2)
     X2in(I) = X2
     WRITE(SLABEL,'(A8,A2)') OTHER1
     SLABEL = SLABEL//c_null_char
     DO J=1,10
       OTHER1in(J,I) = SLABEL(J:J)
     END DO
     WRITE(SLABEL,'(A8,A2)') OTHER2
     SLABEL = SLABEL//c_null_char
     DO J=1,10
       OTHER2in(J,I) = SLABEL(J:J)
     END DO
     ISOSHIFTin(I) = INT(ISOSHIFT)
     NELIONin(I) = INT(NELION)
     RESIDin(I) = RESID
  END DO

  CLOSE(UNIT=1)

end subroutine readoutspecbin

! wrap more functions here
! ...

end module
