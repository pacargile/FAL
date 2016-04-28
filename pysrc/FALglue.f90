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

subroutine readoutspecbin(s, N, wli, qmu1i, qmu2i) bind(c, name='readoutspecbin')
  use iso_c_binding, only: c_double, c_int, c_char, c_null_char
  integer(c_int), intent(in), value :: N
  real(c_double), intent(out) :: wli(N)
  real(c_double), intent(out) :: qmu1i(N)
  real(c_double), intent(out) :: qmu2i(N)
  character(kind=c_char,len=1), intent(in) :: s(*)
  character(len=:), allocatable :: str

  REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,WLEDGE,RATIO,WL
  REAL*8 QMU(40),XMU(20),NEDGE,IFSURF,NMU
  REAL*8 NAV,NMU2,FREQTOWAVE
  INTEGER IWL, NWL, I

  open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
  read(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,WLEDGE
  
  NAV=IFSURF/10
  IFSURF=IFSURF-NAV*10
  IF(NAV.EQ.0.)NAV=1
  IF(IFSURF.EQ.3)NMU=1
  NMU2=NMU+NMU
  RATIO=1.D0+1.D0/RESOLU
  DO 6 IWL=1,NWL
     WL=WBEGIN*RATIO**(IWL-1)
     FREQTOWAVE=2.99792458D17/WL**2
     wli(IWL) = WL
     READ(1)(QMU(I),I=1,NMU2)
     qmu1i(IWL) = QMU(1)*FREQTOWAVE
     qmu2i(IWL) = QMU(2)*FREQTOWAVE
6 CONTINUE
  CLOSE(UNIT=1)

end subroutine readoutspecbin

! wrap more functions here
! ...

end module
