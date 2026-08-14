program navier_stokes_python_triangulos
  
  ! SOLUCIONADOR DE NAVIER-STOKES EN MALLA NO ESTRUCTURADA DE TRIÁNGULOS esta raro xd
  ! Implementa el método de proyección de Chorin (PDF, Sección 2) sobre una
  ! malla de triángulos rectángulos, usando elementos finitos lineales (P1).
  !definicion del cubo 201
  implicit none

  !Parámetros de la malla 
  integer, parameter :: Nx = 200, Ny = 100 !para el paso 360000 mas èqueño
  integer, parameter :: n_nodos     = Nx * Ny
  integer, parameter :: n_triangulos = (Nx-1)*(Ny-1)*2
  
  !Coordenadas y conectividad de la malla
  real :: x(n_nodos), y(n_nodos)
  integer :: tri_nodos(3, n_triangulos)
  real :: area_tri(n_triangulos)
  
  ! Vecindad y matrices esparcidas 
  integer, parameter :: max_vecinos = 8
  integer :: n_vecinos(n_nodos)
  integer :: vecinos(n_nodos, max_vecinos)
  real    :: S_val(n_nodos, max_vecinos)
  real    :: S_diag(n_nodos)
  real    :: M_lump(n_nodos)
  
  ! Campos físicos 
  real :: u(n_nodos), v(n_nodos), p(n_nodos)
  real :: u_star(n_nodos), v_star(n_nodos), p_new(n_nodos)
  real :: grad_u_x(n_nodos), grad_u_y(n_nodos)
  real :: grad_v_x(n_nodos), grad_v_y(n_nodos)
  real :: grad_p_x(n_nodos), grad_p_y(n_nodos)
  real :: div_star(n_nodos)
  real :: conv_u(n_nodos), conv_v(n_nodos)
  
  ! Parámetros físicos 
  real, parameter :: Lx = 14.85, Ly = 4.9
  real, parameter :: dx = Lx / (Nx-1), dy = Ly / (Ny-1)
  real, parameter :: dt = 0.0002 !para la nueva malla 200-100 0.0005
  real, parameter :: rho = 1.0
  real, parameter :: nu = 0.01
  real, parameter :: U_inf = 1.0 !tiene que ser mas pequeña 
  integer, parameter :: MAX_ITER_P = 100 !cuando aumentamos el paso del tiempo 
  
  ! VARIABLES AUXILIARES 
  integer :: i, j, k, tri, iter, iter_p, frame, v_idx
  real :: speed
  character(len=32) :: nombre_archivo

  integer :: n_sw, n_se, n_ne, n_nw
  integer :: n1, n2, n3, nodos(3)
  real :: x1, y1, x2, y2, x3, y3
  real :: b(3), c(3), area, contrib
  integer :: r, s_idx, gi, gj
  logical :: encontrado
  integer :: pos
  real :: dif_u, dif_v
  real :: suma, b_i, suma_off
  
  !!!!! CUBO Y FUEGO !!!!!

  ! Original: cubo_centro_y = Ly/1.0  (centro en y=4.9, flotando)
  ! Nuevo: cubo_centro_y = 1.0 (centro en y=1, base toca el suelo en y=0)
  
  real, parameter :: cubo_centro_x = Lx/2.0
  !real, parameter :: cubo_centro_y = Ly/1.0  ! ORIGINAL (comenta)
  real, parameter :: cubo_centro_y = 1.0      ! NUEVO: cubo en el suelo
  real, parameter :: cubo_mitad_ancho = 2.0, cubo_mitad_alto = 1.0 !anchura del cubo
  
 
  ! MODIFICACIÓN 3: VELOCIDAD DE INYECCIÓN DE COMBUSTIBLE AJUSTABLE
  ! Parámetro nuevo para controlar la velocidad de salida del combustible
  ! Valores sugeridos: 0.1 (lento), 0.5 (moderado), 1.0 (rápido), 2.0 (muy rápido)
  
  real, parameter :: V_combustible = 0.5  ! velocidad de inyección (m/s)
  
  logical :: es_obstaculo(n_nodos)
  real :: T_cubo
  real, parameter :: T_ambiente = 300.0 !temperatura ambiente (K)
  real, parameter :: rho_cubo = 8000.0, cp_cubo = 500.0 !acero 
  real, parameter :: V_cubo = (2.0*cubo_mitad_ancho)*(2.0*cubo_mitad_alto)*1.0
  real, parameter :: A_cubo = 2.0*(2.0*cubo_mitad_ancho + 2.0*cubo_mitad_alto)*1.0
  real, parameter :: h_conv = 50.0
  real :: T_fluido_promedio

  real :: T(n_nodos), T_new(n_nodos)
  real :: YF(n_nodos), YF_new(n_nodos)
  real :: YO(n_nodos), YO_new(n_nodos)
  real, parameter :: D_thermal = 0.05    ! Difusividad termica aumentada para estabilidad
  real, parameter :: D_species = 0.05    ! Difusividad masicaaumentada para estabilidad
  real, parameter :: s = 8.0 !relacion esquetiometrica para combustion de metano con el oxigeno  8 kg de O₂ por 1 kg de combustible
  real, parameter :: Q_comb =2.0e7 !calor liverado 
  real, parameter :: cp_aire = 1004.0 !Calor de combusti´on
  real, parameter :: tau_comb = 0.005  ! Tiempo caracterıstico de la reaccion  ligeramente mayor para suavizar
  real :: omega

  real, parameter :: beta = 1.0 / T_ambiente !Coeficiente de expansion termica
  real, parameter :: g_acel = -9.81 ! Aceleracion de la gravedad (-?)
  real :: fuente_flotacion

  real :: grad_T_x(n_nodos), grad_T_y(n_nodos)
  real :: grad_YF_x(n_nodos), grad_YF_y(n_nodos)
  real :: grad_YO_x(n_nodos), grad_YO_y(n_nodos)
  real :: lap_T(n_nodos), lap_YF(n_nodos), lap_YO(n_nodos)
  real :: adv_T(n_nodos), adv_YF(n_nodos), adv_YO(n_nodos)
  integer :: idx_vecino, k_count
  ! FIN CUBO Y FUEGO !!

  ! CORRECCIÓN ESTABILIDAD- límites para la temperatura !!!!!
  real, parameter :: DT_MAX = 50.0   ! incremento máximo de T por paso (K)
  ! FIN CORRECCIÓN !!!!!

  ! 1. GENERACIÓN DE LA MALLA 
  do j = 1, Ny
    do i = 1, Nx
      k = (j-1)*Nx + i
      x(k) = (i-1) * dx
      y(k) = (j-1) * dy
    end do
  end do
  
  tri = 0
  do j = 1, Ny-1
    do i = 1, Nx-1
      n_sw = (j-1)*Nx + i
      n_se = n_sw + 1
      n_ne = j*Nx + i + 1
      n_nw = j*Nx + i
      
      tri = tri + 1
      tri_nodos(1, tri) = n_sw
      tri_nodos(2, tri) = n_se
      tri_nodos(3, tri) = n_ne
      
      tri = tri + 1
      tri_nodos(1, tri) = n_sw
      tri_nodos(2, tri) = n_ne
      tri_nodos(3, tri) = n_nw
    end do
  end do
  
  do tri = 1, n_triangulos
    n1 = tri_nodos(1, tri); n2 = tri_nodos(2, tri); n3 = tri_nodos(3, tri)
    x1 = x(n1); y1 = y(n1)
    x2 = x(n2); y2 = y(n2)
    x3 = x(n3); y3 = y(n3)
    area_tri(tri) = 0.5 * abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) )
  end do
  
  ! 2. CONSTRUCCIÓN DE VECINDADES Y MATRIZ DE RIGIDEZ
  n_vecinos(:) = 0
  vecinos(:,:) = 0
  S_val(:,:) = 0.0
  S_diag(:) = 0.0
  M_lump(:) = 0.0
  
  do tri = 1, n_triangulos
    nodos(1) = tri_nodos(1, tri)
    nodos(2) = tri_nodos(2, tri)
    nodos(3) = tri_nodos(3, tri)
    x1 = x(nodos(1)); y1 = y(nodos(1))
    x2 = x(nodos(2)); y2 = y(nodos(2))
    x3 = x(nodos(3)); y3 = y(nodos(3))
    area = area_tri(tri)
    
    b(1) = y2 - y3;  c(1) = x3 - x2
    b(2) = y3 - y1;  c(2) = x1 - x3
    b(3) = y1 - y2;  c(3) = x2 - x1
    
    do r = 1, 3
      gi = nodos(r)
      M_lump(gi) = M_lump(gi) + area / 3.0
    end do
    
    do r = 1, 3
      gi = nodos(r)
      do s_idx = 1, 3
        if (r == s_idx) then
          S_diag(gi) = S_diag(gi) + (b(r)*b(r) + c(r)*c(r)) / (4.0 * area)
        else
          gj = nodos(s_idx)
          contrib = (b(r)*b(s_idx) + c(r)*c(s_idx)) / (4.0 * area)
          
          encontrado = .false.
          do pos = 1, n_vecinos(gi)
            if (vecinos(gi, pos) == gj) then
              encontrado = .true.
              exit
            end if
          end do
          if (.not. encontrado) then
            n_vecinos(gi) = n_vecinos(gi) + 1
            if (n_vecinos(gi) > max_vecinos) then
              print *, "Demasiados vecinos en nodo", gi
              stop
            end if
            vecinos(gi, n_vecinos(gi)) = gj
            S_val(gi, n_vecinos(gi)) = contrib
          else
            S_val(gi, pos) = S_val(gi, pos) + contrib
          end if
        end if
      end do
    end do
  end do
  
  do i = 1, n_nodos
    suma_off = 0.0
    do j = 1, n_vecinos(i)
      suma_off = suma_off + S_val(i, j)
    end do
    S_diag(i) = -suma_off
  end do
  
  ! Definición del cubo  
  es_obstaculo = .false.
  do i = 1, n_nodos
      if (x(i) >= cubo_centro_x - cubo_mitad_ancho .and. x(i) <= cubo_centro_x + cubo_mitad_ancho .and. &
          y(i) >= cubo_centro_y - cubo_mitad_alto  .and. y(i) <= cubo_centro_y + cubo_mitad_alto) then
          es_obstaculo(i) = .true.
      end if
  end do

  ! 3. INICIALIZACIÓN 
  u(:) = U_inf; v(:) = 0.0; p(:) = 0.0
  u_star(:) = U_inf; v_star(:) = 0.0; p_new(:) = 0.0
  frame = 0

  T(:) = T_ambiente; YF(:) = 0.0; YO(:) = 0.23 !este queda al eliminar el if del 552
  T_cubo = T_ambiente

  do i = 1, n_nodos
      if (es_obstaculo(i)) then
          u(i) = 0.0; v(i) = 0.0 !cubo quieto 
      end if
  end do

  print *, "Malla generada:", n_nodos, "nodos,", n_triangulos, "triángulos"
  print *, "Nodos dentro del cubo:", count(es_obstaculo)

  open(unit=12, file='malla_nodos.dat', status='replace')
  do i = 1, n_nodos
    write(12, '(2(ES14.6,1X))') x(i), y(i)
  end do
  close(12)
  
  open(unit=13, file='malla_triangulos.dat', status='replace')
  do tri = 1, n_triangulos
    write(13, '(3(I8,1X))') tri_nodos(1,tri)-1, tri_nodos(2,tri)-1, tri_nodos(3,tri)-1
  end do
  close(13)

  ! 5. BUCLE PRINCIPAL 
  do iter = 1, 90000

    ! - TRANSPORTE DE ESCALARES 
    call calcular_gradiente(T,  grad_T_x,  grad_T_y,  n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YF, grad_YF_x, grad_YF_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YO, grad_YO_x, grad_YO_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)

    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            adv_T(i)  = 0.0; lap_T(i)  = 0.0 !no hay transporte de calor en el cubo
            adv_YF(i) = 0.0; lap_YF(i) = 0.0
            adv_YO(i) = 0.0; lap_YO(i) = 0.0
            cycle
        end if

        adv_T(i)  = u(i)*grad_T_x(i)  + v(i)*grad_T_y(i)
        adv_YF(i) = u(i)*grad_YF_x(i) + v(i)*grad_YF_y(i)
        adv_YO(i) = u(i)*grad_YO_x(i) + v(i)*grad_YO_y(i)

        lap_T(i)  = 0.0; lap_YF(i) = 0.0; lap_YO(i) = 0.0
        do idx_vecino = 1, n_vecinos(i)
            j = vecinos(i, idx_vecino)
            lap_T(i)  = lap_T(i)  + S_val(i, idx_vecino) * (T(j)  - T(i))
            lap_YF(i) = lap_YF(i) + S_val(i, idx_vecino) * (YF(j) - YF(i))
            lap_YO(i) = lap_YO(i) + S_val(i, idx_vecino) * (YO(j) - YO(i))
        end do
        lap_T(i)  = -lap_T(i)  / M_lump(i)
        lap_YF(i) = -lap_YF(i) / M_lump(i)
        lap_YO(i) = -lap_YO(i) / M_lump(i)
    end do

    
    ! MODIFICACIÓN 2: INYECCIÓN DE COMBUSTIBLE POR LA PARTE SUPERIOR
    ! Original: inyección por debajo del cubo (y <= cubo_centro_y - cubo_mitad_alto + 0.2)
    ! Nuevo: inyección por encima del cubo (abs(y - (cubo_centro_y + cubo_mitad_alto)) < 0.3)
    ! También se asigna velocidad vertical v = V_combustible
   
        ! INICIO BLOQUE ORIGINAL COMEN-
        !if ( .not. es_obstaculo(i) .and. y(i) <= cubo_centro_y - cubo_mitad_alto + 0.2 .and. &
        !     x(i) >= cubo_centro_x - 0.5 .and. x(i) <= cubo_centro_x + 0.5 ) then
        !    YF(i) = 0.01   ! antes 0.05
        !    YO(i) = 0.0
        !end if
        !- FIN BLOQUE ORIGINAL COMENTADO 
        
        !- INICIO NUEVO BLOQUE (MODIFICACIÓN 2 y 3) 
        if ( .not. es_obstaculo(i) .and. &
             abs(y(i) - (cubo_centro_y + cubo_mitad_alto)) < 0.3 .and. &
             x(i) >= cubo_centro_x - 0.5 .and. x(i) <= cubo_centro_x + 0.5 ) then
            YF(i) = 0.01
            YO(i) = 0.0
            v(i) = V_combustible  ! MODIFICACIÓN 3: velocidad de inyección hacia arriba
        end if
        !---- FIN NUEVO BLOQUE ----
    end do

    ! Tasa de reacción con limitador
    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            omega = 0.0
        else
            omega = min(YF(i), YO(i)/s) / tau_comb
            ! Limitar omega para no consumir más de una fracción máxima por paso
            omega = min(omega, 0.1 / dt)   ! no consumir más de 0.1 fracción másica por dt
        end if

        T_new(i)  = T(i)  + dt * ( -adv_T(i)  + D_thermal * lap_T(i)  + Q_comb * omega / cp_aire )
        YF_new(i) = YF(i) + dt * ( -adv_YF(i) + D_species * lap_YF(i) - omega )
        YO_new(i) = YO(i) + dt * ( -adv_YO(i) + D_species * lap_YO(i) - s * omega )
    end do

    ! Actualizar escalares con límite de incremento de T
    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            T(i)  = T_cubo  ! El cubo impone su temperatura al fluido en la frontera
            YF(i) = 0.0 ! No hay combustible dentro del sólido
            YO(i) = 0.0   ! No hay oxidante dentro del sólido
        else
            ! Limitar incremento de temperatura
            if (T_new(i) - T(i) > DT_MAX) T_new(i) = T(i) + DT_MAX
            if (T(i) - T_new(i) > DT_MAX) T_new(i) = T(i) - DT_MAX

            T(i)  = max(T_ambiente, T_new(i))
            YF(i) = max(0.0, min(1.0, YF_new(i)))
            YO(i) = max(0.0, min(1.0, YO_new(i)))
        end if
    end do

    ! Calentamiento del cubo (sin cambios)
    T_fluido_promedio = 0.0
    k_count = 0
    do i = 1, n_nodos
        if (.not. es_obstaculo(i)) then
            do idx_vecino = 1, n_vecinos(i)
                j = vecinos(i, idx_vecino)
                if (es_obstaculo(j)) then
                    T_fluido_promedio = T_fluido_promedio + T(i)
                    k_count = k_count + 1
                    exit
                end if
            end do
        end if
    end do
    if (k_count > 0) T_fluido_promedio = T_fluido_promedio / real(k_count)
    T_cubo = T_cubo + dt * (h_conv * A_cubo) / (rho_cubo * cp_cubo * V_cubo) * (T_fluido_promedio - T_cubo)
    do i = 1, n_nodos
        if (es_obstaculo(i)) T(i) = T_cubo
    end do

    !CORRECCIÓN ESTABILIDAD: condiciones de contorno para escalares !
    call condiciones_contorno_escalares(T, YF, YO, x, y, n_nodos, Lx, Ly, T_ambiente)
    !FIN CORRECCI !!!!!

    ! FIN ESCALARES ----


    !  Paso 5.1: Gradientes de velocidad
    call calcular_gradiente(u, grad_u_x, grad_u_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(v, grad_v_x, grad_v_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    
    do i = 1, n_nodos
      conv_u(i) = u(i) * grad_u_x(i) + v(i) * grad_u_y(i)
      conv_v(i) = u(i) * grad_v_x(i) + v(i) * grad_v_y(i)
    end do
    
    !  Paso 5.2: Predicción con flotación
    do i = 1, n_nodos
      if (es_obstaculo(i)) then
          u_star(i) = 0.0; v_star(i) = 0.0 ! Velocidad nula en el sólido
          cycle
      end if

      dif_u = 0.0; dif_v = 0.0
      do v_idx = 1, n_vecinos(i)
        j = vecinos(i, v_idx)
        dif_u = dif_u + S_val(i, v_idx) * (u(j) - u(i))
        dif_v = dif_v + S_val(i, v_idx) * (v(j) - v(i))
      end do
      dif_u = -nu * dif_u / M_lump(i)
      dif_v = -nu * dif_v / M_lump(i)
      
      fuente_flotacion = -g_acel * beta * (T(i) - T_ambiente)
      
      u_star(i) = u(i) + dt * ( -conv_u(i) + dif_u )
      v_star(i) = v(i) + dt * ( -conv_v(i) + dif_v + fuente_flotacion )
    end do
    
    call condiciones_contorno_velocidad(u_star, v_star, x, y, n_nodos, Lx, Ly, U_inf)
    
    ! Paso 5.3: Divergencia
    call calcular_gradiente(u_star, grad_u_x, grad_u_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(v_star, grad_v_x, grad_v_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    do i = 1, n_nodos
      div_star(i) = grad_u_x(i) + grad_v_y(i)
    end do
    
    ! Paso 5.4: Poisson
    p_new = p
    do iter_p = 1, MAX_ITER_P
      do i = 1, n_nodos
        if (abs(x(i) - Lx) < 1.0e-6) then
          p_new(i) = 0.0
          cycle
        end if
        if (es_obstaculo(i)) then
            p_new(i) = p(i)
            cycle
        end if
        
        b_i = -(rho / dt) * div_star(i) * M_lump(i)
        suma = b_i
        do v_idx = 1, n_vecinos(i)
          j = vecinos(i, v_idx)
          suma = suma - S_val(i, v_idx) * p(j)
        end do
        
        if (abs(S_diag(i)) > 1.0e-12) then
          p_new(i) = suma / S_diag(i)
        else
          p_new(i) = 0.0
        end if
      end do
      p = p_new
    end do
    
    ! Paso 5.5: Corrección
    call calcular_gradiente(p, grad_p_x, grad_p_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    do i = 1, n_nodos
      u(i) = u_star(i) - (dt / rho) * grad_p_x(i)
      v(i) = v_star(i) - (dt / rho) * grad_p_y(i)
    end do
    
    call condiciones_contorno_velocidad(u, v, x, y, n_nodos, Lx, Ly, U_inf)

    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            u(i) = 0.0; v(i) = 0.0
        end if
    end do

    ! Detección de NaN jeje
    do i = 1, n_nodos
        if (u(i) /= u(i) .or. v(i) /= v(i) .or. p(i) /= p(i) .or. T(i) /= T(i)) then
            print *, "NaN detectado en el paso", iter
            stop
        end if
    end do
    
    ! Escritura cada 500 pasos
    if (mod(iter, 500) == 0) then
      frame = frame + 1
      write(nombre_archivo, '(A,I4.4,A)') 'campo_', frame, '.dat'
      open(unit=10, file=nombre_archivo, status='replace')
      do i = 1, n_nodos
        speed = sqrt(u(i)**2 + v(i)**2)
        write(10, '(4(ES14.6,1X))') u(i), v(i), p(i), speed
      end do
      close(10)

      write(nombre_archivo, '(A,I4.4,A)') 'temperatura_', frame, '.dat'
      open(unit=11, file=nombre_archivo, status='replace')
      do i = 1, n_nodos
        write(11, '(3(ES14.6,1X))') x(i), y(i), T(i)
      end do
      close(11)

      print *, "Paso ", iter, " T_cubo =", T_cubo
    end if
    
  end do
  
  print *, "Simulación finalizada."

contains

  subroutine calcular_gradiente(f, df_dx, df_dy, n_nodos, n_tri, tri_nodos, x, y, area_tri, M_lump)
    integer, intent(in) :: n_nodos, n_tri
    real, intent(in)    :: f(n_nodos), x(n_nodos), y(n_nodos)
    integer, intent(in) :: tri_nodos(3, n_tri)
    real, intent(in)    :: area_tri(n_tri), M_lump(n_nodos)
    real, intent(out)   :: df_dx(n_nodos), df_dy(n_nodos)
    
    real :: suma_dx(n_nodos), suma_dy(n_nodos), suma_area(n_nodos)
    integer :: tri, r, nodo, n1, n2, n3, i
    real :: x1, y1, x2, y2, x3, y3, area, b(3), c(3)
    real :: f1, f2, f3, gx, gy
    
    suma_dx(:) = 0.0; suma_dy(:) = 0.0; suma_area(:) = 0.0
    
    do tri = 1, n_tri
      n1 = tri_nodos(1,tri); n2 = tri_nodos(2,tri); n3 = tri_nodos(3,tri)
      x1 = x(n1); y1 = y(n1)
      x2 = x(n2); y2 = y(n2)
      x3 = x(n3); y3 = y(n3)
      area = area_tri(tri)
      
      b(1) = y2 - y3; c(1) = x3 - x2
      b(2) = y3 - y1; c(2) = x1 - x3
      b(3) = y1 - y2; c(3) = x2 - x1
      
      f1 = f(n1); f2 = f(n2); f3 = f(n3)
      gx = (f1*b(1) + f2*b(2) + f3*b(3)) / (2.0 * area)
      gy = (f1*c(1) + f2*c(2) + f3*c(3)) / (2.0 * area)
      
      do r = 1, 3
        nodo = tri_nodos(r, tri)
        suma_dx(nodo) = suma_dx(nodo) + area * gx
        suma_dy(nodo) = suma_dy(nodo) + area * gy
        suma_area(nodo) = suma_area(nodo) + area
      end do
    end do
    
    do i = 1, n_nodos
      if (suma_area(i) > 1e-12) then
        df_dx(i) = suma_dx(i) / suma_area(i)
        df_dy(i) = suma_dy(i) / suma_area(i)
      else
        df_dx(i) = 0.0; df_dy(i) = 0.0
      end if
    end do
  end subroutine calcular_gradiente

  subroutine condiciones_contorno_velocidad(u, v, x, y, n_nodos, Lx, Ly, U_inf)
    integer, intent(in) :: n_nodos
    real, intent(inout) :: u(n_nodos), v(n_nodos)
    real, intent(in)    :: x(n_nodos), y(n_nodos)
    real, intent(in)    :: Lx, Ly, U_inf
    
    integer :: i, j, v_idx
    logical :: found
    
    do i = 1, n_nodos
      if (abs(x(i)) < 1.0e-6) then
        u(i) = U_inf; v(i) = 0.0
      end if
      if (abs(y(i)) < 1.0e-6 .or. abs(y(i) - Ly) < 1.0e-6) then
        u(i) = 0.0; v(i) = 0.0
      end if
      if (abs(x(i) - Lx) < 1.0e-6) then
        found = .false. 
        do v_idx = 1, n_vecinos(i)
          j = vecinos(i, v_idx)
          if (x(j) < Lx - 1.0e-6) then
            u(i) = u(j); v(i) = v(j)
            found = .true.
            exit
          end if
        end do
        if (.not. found) then
          u(i) = U_inf * 0.5; v(i) = 0.0
        end if
      end if
    end do
  end subroutine condiciones_contorno_velocidad

  !NUEVA SUBRUTINA: condiciones de contorno para escalares !!!!!
  subroutine condiciones_contorno_escalares(T, YF, YO, x, y, n_nodos, Lx, Ly, T_amb)
    integer, intent(in) :: n_nodos
    real, intent(inout) :: T(n_nodos), YF(n_nodos), YO(n_nodos)
    real, intent(in)    :: x(n_nodos), y(n_nodos)
    real, intent(in)    :: Lx, Ly, T_amb
    integer :: i

    do i = 1, n_nodos
      ! Entrada: valores fijos de ambiente
      if (abs(x(i)) < 1.0e-6) then
          T(i) = T_amb
          YF(i) = 0.0
          YO(i) = 0.23 !repocicion del oxigeno comentar para quitar 
      end if
      ! Paredes superior e inferior: gradiente nulo (aproximación: copiar del interior si no es obstáculo)
      if (abs(y(i)) < 1.0e-6 .or. abs(y(i) - Ly) < 1.0e-6) then
          ! Se impone adiabático: no hacer nada, el gradiente nulo se debería resolver vía la matriz.
          ! Para forzarlo explícitamente, se podría igualar al primer nodo interior, pero no l.
          continue
      end if
      ! Salida: gradiente nulo (ya se maneja en el cálculo de laplaciano, pero forzamos por seguri
      if (abs(x(i) - Lx) < 1.0e-6) then
          ! no se modifica, se mantiene el valor calculado
          continue
      end if
    end do
  end subroutine condiciones_contorno_escalares

end program navier_stokes_python_triangulos