module utils_solver
  implicit none
  integer, parameter :: dp = kind(1.0d0)
contains
  ! Método de Thomas para sistema tridiagonal
  subroutine thomas_solver(n, a, b, c, d, x)
    integer, intent(in) :: n
    real(dp), intent(inout) :: a(:), b(:), c(:), d(:)
    real(dp), intent(out) :: x(:)
    integer :: i
    real(dp), allocatable :: bb(:), dd(:)
    allocate(bb(n), dd(n))
    bb = b; dd = d
    do i = 2, n
      bb(i) = bb(i) - a(i-1)*c(i-1)/bb(i-1)
      dd(i) = dd(i) - a(i-1)*dd(i-1)/bb(i-1)
    end do
    x(n) = dd(n)/bb(n)
    do i = n-1, 1, -1
      x(i) = (dd(i) - c(i)*x(i+1)) / bb(i)
    end do
    deallocate(bb, dd)
  end subroutine thomas_solver
end module utils_solver
