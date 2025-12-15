program test_odepack_sub1_dewset
  implicit none
  integer, parameter :: n = 4
  double precision :: rtol(n), atol(n), ycur(n), ewt(n)
  integer :: unit, itol

  rtol = (/ 1.0d-3, 2.0d-3, 3.0d-3, 4.0d-3 /)
  atol = (/ 1.0d-6, 2.0d-6, 3.0d-6, 4.0d-6 /)
  ycur = (/ 1.0d0, -2.0d0, 0.5d0, -1.5d0 /)

  open(newunit=unit, file="dewset.bin", access="stream", form="unformatted", status="replace")
  do itol = 1, 4
     call dewset(n, itol, rtol, atol, ycur, ewt)
     write(unit) ewt
  end do
  close(unit)
end program test_odepack_sub1_dewset
