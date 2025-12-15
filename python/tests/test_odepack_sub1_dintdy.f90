program odepack_sub1_dintdy_driver
  use iso_c_binding, only: c_int
  use odepack_common
  implicit none
  integer, parameter :: n = 2, l = 3, nyh = 2
  double precision :: yh(nyh,l), dky(n)
  integer :: iflag, k, iunit
  type(odepack_common_data) :: common_data

  yh(:,1) = (/1.0d0, 2.0d0/)
  yh(:,2) = (/0.1d0, 0.2d0/)
  yh(:,3) = (/0.01d0, 0.02d0/)

  common_data%DLS001%ints(19) = l      ! L
  common_data%DLS001%ints(32) = n      ! N
  common_data%DLS001%ints(33) = 2      ! NQ
  common_data%DLS001%reals(212) = 0.5d0  ! H
  common_data%DLS001%reals(215) = 0.5d0  ! HU
  common_data%DLS001%reals(217) = 1.0d0  ! TN
  common_data%DLS001%reals(218) = 1.0d-16 ! UROUND

  open(newunit=iunit, file="dintdy_ref.bin", access="stream", form="unformatted", status="replace")
  do k = 0,2
     call dintdy(1.0d0, k, yh, nyh, dky, iflag, common_data)
     write(iunit) k
     write(iunit) iflag
     write(iunit) dky
  end do
  close(iunit)
end program odepack_sub1_dintdy_driver
