program experiment_A
  use utils_solver
  implicit none
  integer :: i, n, m, nt
  real(dp) :: L, alpha, T, h, k, lambda
  real(dp), allocatable :: x(:), u(:), u_new(:), a(:), b(:), c(:), rhs(:)
  real(dp) :: start, finish
  integer, dimension(4) :: m_list = [10,20,40,80]

  L = 1.0_dp; alpha = 1.0_dp; T = 0.1_dp; k = 0.01_dp
  print *, "===== Experimento A (k fijo = 0.01) ====="

  do n = 1, size(m_list)
    m = m_list(n)
    h = L / m
    lambda = alpha**2 * k / (h**2)
    nt = int(T / k)
    allocate(x(0:m), u(0:m), u_new(0:m), a(m-2), b(m-1), c(m-2), rhs(m-1))
    
    do i = 0, m
      x(i) = i*h
      u(i) = sin(acos(-1.0_dp)*x(i))
    end do
    
    a = -lambda; b = 1.0_dp + 2.0_dp*lambda; c = -lambda
    
    call cpu_time(start)
    do i = 1, nt
      rhs = u(1:m-1)
      call thomas_solver(m-1, a, b, c, rhs, u_new(1:m-1))
      u_new(0) = 0.0_dp; u_new(m) = 0.0_dp
      u = u_new
    end do
    call cpu_time(finish)

    print '(A,I4,2F10.4, A, F8.6)', "m=", m, h, lambda, "  tiempo=", finish - start

    deallocate(x,u,u_new,a,b,c,rhs)
  end do
end program experiment_A
