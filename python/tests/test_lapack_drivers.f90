program lapack_drivers
  implicit none
  call run_dense()
  call run_band()
contains
  subroutine run_dense()
    integer, parameter :: n = 3, lda = 3
    double precision :: a(lda,n), b(n)
    integer :: ipiv(n), info, unit
    a = reshape((/ 4d0, 1d0, 0d0, 2d0, 5d0, 1d0, 0d0, 2d0, 3d0 /), (/ lda, n /))
    b = (/ 1d0, 2d0, 3d0 /)
    call dgetrf(n, n, a, lda, ipiv, info)
    call dgetrs('N', n, 1, a, lda, ipiv, b, n, info)
    open(newunit=unit, file="dgetrs_ref.bin", access="stream", form="unformatted", status="replace")
    write(unit) info
    write(unit) ipiv
    write(unit) b
    close(unit)
  end subroutine run_dense

  subroutine run_band()
    integer, parameter :: n = 4, kl = 1, ku = 1, ldab = 2*kl+ku+1
    double precision :: ab(ldab,n), b(n), adense(n,n)
    integer :: ipiv(n), info, unit, i
    adense = 0.0d0
    do i = 1, n
       adense(i,i) = 2.0d0
    end do
    adense(1,2) = -1.0d0
    adense(2,1) = -1.0d0
    adense(2,3) = -1.0d0
    adense(3,2) = -1.0d0
    adense(3,4) = -1.0d0
    adense(4,3) = -1.0d0
    b = (/ 1d0, 0d0, 0d0, 1d0 /)
    call dgesv(n, 1, adense, n, ipiv, b, n, info)
    open(newunit=unit, file="dgbtrs_ref.bin", access="stream", form="unformatted", status="replace")
    write(unit) info
    write(unit) b
    close(unit)
  end subroutine run_band
end program lapack_drivers
