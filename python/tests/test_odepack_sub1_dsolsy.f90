program test_odepack_sub1_dsolsy
  use odepack_common
  implicit none
  integer, parameter :: n = 3, lda = 3
  double precision :: a(lda,n), b(n)
  integer :: ipiv(n), info, unit
  type(odepack_common_data) :: common

  a = reshape((/ 4d0, 1d0, 0d0, 2d0, 5d0, 1d0, 0d0, 2d0, 3d0 /), (/ lda, n /))
  b = (/ 1d0, 2d0, 3d0 /)
  call dgetrf(n, n, a, lda, ipiv, info)
  call dgetrs('N', n, 1, a, lda, ipiv, b, n, info)

  open(newunit=unit, file="dsolsy_dense.bin", access="stream", form="unformatted", status="replace")
  write(unit) info
  write(unit) ipiv
  write(unit) a
  write(unit) b
  close(unit)
end program test_odepack_sub1_dsolsy
