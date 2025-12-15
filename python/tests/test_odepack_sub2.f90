program odepack_sub2_driver
  implicit none
  external xerrwd, ixsav
  integer :: unit

  unit = 77
  call ixsav(1, unit, .true.)
  call ixsav(2, 1, .true.)
  open(unit=unit, file="odepack_sub2.txt", status="replace")

  call xerrwd('TEST MSG', 8, 0, 0, 0, 0, 0, 0, 0d0, 0d0)
  call xerrwd('WITH I', 6, 0, 0, 1, 123, 0, 0, 0d0, 0d0)
  call xerrwd('WITH I12', 7, 0, 0, 2, 12, 34, 0, 0d0, 0d0)
  call xerrwd('WITH R', 6, 0, 0, 0, 0, 0, 1, 1.23456789d0, 0d0)
  call xerrwd('WITH R12', 7, 0, 0, 0, 0, 0, 2, 1.0d0, -2.5d0)

  close(unit)
end program odepack_sub2_driver
