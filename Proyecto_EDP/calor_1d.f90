! ============================================================================
! MÉTODO EXPLÍCITO PARA LA ECUACIÓN DEL CALOR 1D
! Esquema: u_i^(n+1) = λ·u_(i-1)^n + (1-2λ)·u_i^n + λ·u_(i+1)^n
! Análisis de estabilidad: λ ≤ 0.5
! ============================================================================

module heat_equation_explicit
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  
contains
  
  ! --------------------------------------------------------------------------
  ! Solución exacta: u(x,t) = sin(πx/L) · exp(-α²π²t/L²)
  ! --------------------------------------------------------------------------
  function exact_solution(x, t, L, alpha) result(u)
    real(dp), intent(in) :: x, t, L, alpha
    real(dp) :: u
    u = sin(pi * x / L) * exp(-alpha**2 * pi**2 * t / L**2)
  end function exact_solution
  
  ! --------------------------------------------------------------------------
  ! Condición inicial: u(x,0) = sin(πx/L)
  ! --------------------------------------------------------------------------
  subroutine initialize_condition(u, x, m, L)
    integer, intent(in) :: m
    real(dp), intent(in) :: L
    real(dp), intent(inout) :: u(0:m), x(0:m)
    integer :: i
    
    do i = 0, m
      x(i) = real(i, dp) * L / real(m, dp)
      u(i) = sin(pi * x(i) / L)
    end do
    
    ! Condiciones de frontera
    u(0) = 0.0_dp
    u(m) = 0.0_dp
  end subroutine initialize_condition
  
  ! --------------------------------------------------------------------------
  ! Método explícito: Un paso temporal
  ! --------------------------------------------------------------------------
  subroutine explicit_step(u_old, u_new, m, lambda)
    integer, intent(in) :: m
    real(dp), intent(in) :: lambda
    real(dp), intent(in) :: u_old(0:m)
    real(dp), intent(out) :: u_new(0:m)
    integer :: i
    
    ! Puntos interiores
    do i = 1, m-1
      u_new(i) = lambda * u_old(i-1) + (1.0_dp - 2.0_dp*lambda) * u_old(i) + lambda * u_old(i+1)
    end do
    
    ! Condiciones de frontera
    u_new(0) = 0.0_dp
    u_new(m) = 0.0_dp
  end subroutine explicit_step
  
  ! --------------------------------------------------------------------------
  ! Resolver la ecuación completa
  ! --------------------------------------------------------------------------
  subroutine solve_heat_equation(u, x, m, nt, L, alpha, T, lambda, elapsed_time)
    integer, intent(in) :: m, nt
    real(dp), intent(in) :: L, alpha, T, lambda
    real(dp), intent(inout) :: u(0:m), x(0:m)
    real(dp), intent(out) :: elapsed_time
    real(dp), allocatable :: u_temp(:)
    integer :: n
    real(dp) :: start_time, end_time
    
    allocate(u_temp(0:m))
    
    ! Inicializar condición inicial
    call initialize_condition(u, x, m, L)
    
    ! Medir tiempo
    call cpu_time(start_time)
    
    ! Iterar en el tiempo
    do n = 1, nt
      call explicit_step(u, u_temp, m, lambda)
      u = u_temp
    end do
    
    call cpu_time(end_time)
    elapsed_time = end_time - start_time
    
    deallocate(u_temp)
  end subroutine solve_heat_equation
  
  ! --------------------------------------------------------------------------
  ! Calcular errores respecto a la solución exacta
  ! --------------------------------------------------------------------------
  subroutine compute_errors(u, x, m, t, L, alpha, error_max, error_l2)
    integer, intent(in) :: m
    real(dp), intent(in) :: u(0:m), x(0:m), t, L, alpha
    real(dp), intent(out) :: error_max, error_l2
    real(dp) :: u_exact, diff, h
    integer :: i
    
    h = L / real(m, dp)
    error_max = 0.0_dp
    error_l2 = 0.0_dp
    
    do i = 0, m
      u_exact = exact_solution(x(i), t, L, alpha)
      diff = abs(u(i) - u_exact)
      error_max = max(error_max, diff)
      error_l2 = error_l2 + diff**2
    end do
    
    error_l2 = sqrt(error_l2 * h)
  end subroutine compute_errors
  
  ! --------------------------------------------------------------------------
  ! Verificar condición de estabilidad
  ! --------------------------------------------------------------------------
  function check_stability(lambda) result(is_stable)
    real(dp), intent(in) :: lambda
    logical :: is_stable
    is_stable = (lambda <= 0.5_dp)
  end function check_stability
  
end module heat_equation_explicit

! ============================================================================
! PROGRAMA PRINCIPAL: EXPERIMENTOS
! ============================================================================

program main
  use heat_equation_explicit
  implicit none
  
  ! Parámetros físicos
  real(dp) :: L, alpha, T
  
  ! Parámetros numéricos
  integer :: m, nt
  real(dp) :: h, k, lambda
  
  ! Arrays
  real(dp), allocatable :: u(:), x(:)
  
  ! Variables de análisis
  real(dp) :: error_max, error_l2, elapsed_time
  logical :: is_stable
  
  ! Contadores
  integer :: i, experiment
  
  print *, "========================================================================"
  print *, "   MÉTODO DE DIFERENCIAS FINITAS EXPLÍCITO - ECUACIÓN DEL CALOR 1D"
  print *, "========================================================================"
  print *, ""
  
  ! Parámetros del problema
  L = 1.0_dp
  alpha = 1.0_dp
  T = 0.1_dp
  
  ! =========================================================================
  ! EXPERIMENTO 1: ANÁLISIS DE ESTABILIDAD (λ variable, m fijo)
  ! =========================================================================
  
  print *, "EXPERIMENTO 1: ANÁLISIS DE ESTABILIDAD"
  print *, "----------------------------------------------------------------------"
  print *, "Condición de estabilidad: λ ≤ 0.5"
  print *, ""
  
  m = 40
  h = L / real(m, dp)
  
  print '(A)', "  Lambda    Estable?    h         k           nt       Error Max      Error L2     Tiempo(s)"
  print '(A)', "  --------  ----------  --------  ----------  -------  -------------  -----------  ---------"
  
  ! Probar diferentes valores de lambda
  do experiment = 1, 6
    select case(experiment)
      case(1); lambda = 0.3_dp
      case(2); lambda = 0.45_dp
      case(3); lambda = 0.5_dp
      case(4); lambda = 0.55_dp
      case(5); lambda = 0.75_dp
      case(6); lambda = 1.0_dp
    end select
    
    k = lambda * h**2 / alpha**2
    nt = int(T / k)
    
    allocate(u(0:m), x(0:m))
    
    call solve_heat_equation(u, x, m, nt, L, alpha, T, lambda, elapsed_time)
    call compute_errors(u, x, m, T, L, alpha, error_max, error_l2)
    is_stable = check_stability(lambda)
    
    if (is_stable) then
      print '(2X, F6.3, 4X, A, 3X, F8.6, 2X, F10.8, 2X, I7, 2X, ES13.6, 2X, ES11.4, 2X, F9.6)', &
            lambda, "SI", h, k, nt, error_max, error_l2, elapsed_time
    else
      print '(2X, F6.3, 4X, A, 3X, F8.6, 2X, F10.8, 2X, I7, 2X, ES13.6, 2X, ES11.4, 2X, F9.6)', &
            lambda, "NO", h, k, nt, error_max, error_l2, elapsed_time
    end if
    
    deallocate(u, x)
  end do
  
  print *, ""
  print *, "OBSERVACIONES:"
  print *, "  - Para λ ≤ 0.5: El método es estable y convergente"
  print *, "  - Para λ > 0.5: El método es INESTABLE (errores crecen)"
  print *, ""
  
  ! =========================================================================
  ! EXPERIMENTO 2: CONVERGENCIA (m variable, λ fijo)
  ! =========================================================================
  
  print *, "========================================================================"
  print *, "EXPERIMENTO 2: ESTUDIO DE CONVERGENCIA"
  print *, "----------------------------------------------------------------------"
  print *, "λ fijo = 0.4 (estable)"
  print *, ""
  
  lambda = 0.4_dp
  
  print '(A)', "    m      h         k           nt       Error Max      Error L2     Tiempo(s)   Nodos·Pasos"
  print '(A)', "  -----  --------  ----------  -------  -------------  -----------  ---------  -------------"
  
  do experiment = 1, 6
    select case(experiment)
      case(1); m = 10
      case(2); m = 20
      case(3); m = 40
      case(4); m = 80
      case(5); m = 160
      case(6); m = 320
    end select
    
    h = L / real(m, dp)
    k = lambda * h**2 / alpha**2
    nt = int(T / k)
    
    allocate(u(0:m), x(0:m))
    
    call solve_heat_equation(u, x, m, nt, L, alpha, T, lambda, elapsed_time)
    call compute_errors(u, x, m, T, L, alpha, error_max, error_l2)
    
    print '(2X, I5, 2X, F8.6, 2X, F10.8, 2X, I7, 2X, ES13.6, 2X, ES11.4, 2X, F9.6, 2X, I13)', &
          m, h, k, nt, error_max, error_l2, elapsed_time, (m+1)*nt
    
    deallocate(u, x)
  end do
  
  print *, ""
  print *, "OBSERVACIONES:"
  print *, "  - Al reducir h, el error disminuye como O(h²)"
  print *, "  - El número de pasos temporales crece como nt ~ 1/h²"
  print *, "  - Costo computacional total ~ O(1/h³)"
  print *, ""
  
  ! =========================================================================
  ! EXPERIMENTO 3: LÍMITE DE ENMALLADO
  ! =========================================================================
  
  print *, "========================================================================"
  print *, "EXPERIMENTO 3: ANÁLISIS DE LÍMITE DE ENMALLADO"
  print *, "----------------------------------------------------------------------"
  print *, "Probando mallas muy refinadas (λ = 0.45)"
  print *, ""
  
  lambda = 0.45_dp
  
  print '(A)', "    m       Nodos      nt        k          Ops·10⁶   Tiempo(s)   Memoria(KB)"
  print '(A)', "  ------  ---------  --------  ----------  ---------  ---------  ------------"
  
  do experiment = 1, 5
    select case(experiment)
      case(1); m = 100
      case(2); m = 200
      case(3); m = 500
      case(4); m = 1000
      case(5); m = 2000
    end select
    
    h = L / real(m, dp)
    k = lambda * h**2 / alpha**2
    nt = int(T / k)
    
    allocate(u(0:m), x(0:m))
    
    call solve_heat_equation(u, x, m, nt, L, alpha, T, lambda, elapsed_time)
    
    print '(2X, I6, 2X, I9, 2X, I8, 2X, F10.8, 2X, F9.3, 2X, F9.6, 2X, F12.2)', &
          m, m+1, nt, k, real((m-1)*nt*5)/1.0e6_dp, elapsed_time, real((m+1)*8)/1024.0_dp
    
    deallocate(u, x)
  end do
  
  print *, ""
  print *, "LÍMITE PRÁCTICO:"
  print *, "  - Para m ~ 2000: nt ~ 1.6·10⁶ pasos temporales"
  print *, "  - Costo: ~ 16 millones de operaciones"
  print *, "  - Limitación principal: Restricción de estabilidad k ~ h²"
  print *, ""
  
  ! =========================================================================
  ! EXPERIMENTO 4: COMPARACIÓN DE ESQUEMAS
  ! =========================================================================
  
  print *, "========================================================================"
  print *, "EXPERIMENTO 4: MÉTODO EXPLÍCITO vs CRANK-NICOLSON"
  print *, "----------------------------------------------------------------------"
  print *, ""
  print *, "Características del método EXPLÍCITO:"
  print *, "  ✓ Implementación muy simple"
  print *, "  ✓ No requiere resolver sistemas de ecuaciones"
  print *, "  ✓ Bajo costo por paso temporal"
  print *, "  ✗ CONDICIONALMENTE ESTABLE (λ ≤ 0.5)"
  print *, "  ✗ Requiere muchos pasos temporales para mallas finas"
  print *, "  ✗ Error temporal O(k)"
  print *, ""
  print *, "Características del método CRANK-NICOLSON:"
  print *, "  ✓ Incondicionalmente estable (cualquier λ)"
  print *, "  ✓ Mayor precisión temporal O(k²)"
  print *, "  ✓ Permite pasos temporales más grandes"
  print *, "  ✗ Requiere resolver sistemas tridiagonales"
  print *, "  ✗ Implementación más compleja"
  print *, ""
  
  print *, "========================================================================"
  print *, "RESUMEN DE CONCLUSIONES"
  print *, "========================================================================"
  print *, ""
  print *, "1. ESTABILIDAD:"
  print *, "   - El método explícito es estable SOLO si λ ≤ 0.5"
  print *, "   - Para λ > 0.5, las soluciones divergen (inestabilidad numérica)"
  print *, ""
  print *, "2. CONVERGENCIA:"
  print *, "   - Error espacial: O(h²)"
  print *, "   - Error temporal: O(k)"
  print *, "   - Para estabilidad: k ≤ h²/(2α²)"
  print *, ""
  print *, "3. COSTO COMPUTACIONAL:"
  print *, "   - Operaciones por paso: O(m)"
  print *, "   - Pasos temporales: O(1/h²)"
  print *, "   - Costo total: O(m/h²) = O(1/h³)"
  print *, ""
  print *, "4. LÍMITE PRÁCTICO:"
  print *, "   - Mallado fino requiere MUCHOS pasos temporales"
  print *, "   - Para m=2000: nt ≈ 1.6 millones de pasos"
  print *, "   - Alternativa: Usar métodos implícitos (Crank-Nicolson)"
  print *, ""
  print *, "========================================================================"
  
end program main