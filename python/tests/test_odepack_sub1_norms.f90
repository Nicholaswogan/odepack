program test_odepack_sub1_norms
  implicit none
  integer, parameter :: n = 4, ml = 1, mu = 1, nra = ml + mu + 1
  double precision :: v(n), w(n), a(n,n), ab(nra,n)
  double precision :: r_dm, r_df, r_db
  double precision dmnorm, dfnorm, dbnorm
  external dmnorm, dfnorm, dbnorm
  integer :: i, j, unit

  v = (/ 1.0d0, -2.0d0, 0.5d0, -1.5d0 /)
  w = (/ 1.0d0, 0.5d0, 2.0d0, 1.5d0 /)

  a = reshape((/ &
       2.0d0, -1.0d0,  0.0d0,  0.0d0, &
      -1.0d0,  2.0d0, -1.0d0,  0.0d0, &
       0.0d0, -1.0d0,  2.0d0, -1.0d0, &
       0.0d0,  0.0d0, -1.0d0,  2.0d0 /), (/ n, n /))

  ab = 0.0d0
  do j = 1, n
     do i = max(1, j-mu), min(n, j+ml)
        ab(mu+1 + i - j, j) = a(i, j)
     end do
  end do

  r_dm = dmnorm(n, v, w)
  r_df = dfnorm(n, a, w)
  r_db = dbnorm(n, ab, nra, ml, mu, w)

  open(newunit=unit, file="norms.bin", access="stream", form="unformatted", status="replace")
  write(unit) r_dm
  write(unit) r_df
  write(unit) r_db
  close(unit)
end program test_odepack_sub1_norms
