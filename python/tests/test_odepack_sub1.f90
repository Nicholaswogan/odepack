program odepack_sub1_driver
  implicit none
  double precision dumach
  external dumach, dcfode
  integer :: i, j, unit
  double precision :: elco1(13,12), tesco1(3,12), elco2(13,12), tesco2(3,12)

  do i = 1, 13
     do j = 1, 12
        elco1(i,j) = 7.7d0
        elco2(i,j) = 7.7d0
     end do
  end do
  do i = 1, 3
     do j = 1, 12
        tesco1(i,j) = 8.8d0
        tesco2(i,j) = 8.8d0
     end do
  end do

  call dcfode(1, elco1, tesco1)
  call dcfode(2, elco2, tesco2)

  open(newunit=unit, file="odepack_sub1.bin", access="stream", form="unformatted", status="replace")
  write(unit) dumach()
  write(unit) elco1
  write(unit) tesco1
  write(unit) elco2
  write(unit) tesco2
  close(unit)
end program odepack_sub1_driver
