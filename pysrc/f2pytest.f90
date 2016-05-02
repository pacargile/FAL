subroutine chararr_test(strlen,nlines,chararr)
  implicit None
  !f2py intent(in) :: strlen
  !f2py intent(in) :: nlines
  !f2py intent(out) :: chararr
  integer :: strlen
  integer :: nlines
  character(len=strlen) :: chararr(nlines)
  integer :: i
  do i = 1, nlines
    chararr(i) = 'AAAA'
    end do
end subroutine chararr_test
