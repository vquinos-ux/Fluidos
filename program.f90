program navier_stokes_python_triangulos

  ! SOLUCIONADOR DE NAVIER-STOKES EN MALLA NO ESTRUCTURADA DE TRIÁNGULOS
  ! Implementa el método de proyección de Chorin (PDF, Sección 2) sobre una
  ! malla de triángulos rectángulos, usando elementos finitos lineales (P1).
  ! Versión mejorada con termodinámica real, densidad variable, especies y sublimación real, con viscosidad contante.
 
  implicit none

  ! Parámetros de la malla
  integer, parameter :: Nx = 200, Ny = 100
  integer, parameter :: n_nodos     = Nx * Ny
  integer, parameter :: n_triangulos = (Nx-1)*(Ny-1)*2

  ! Coordenadas y conectividad de la malla
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
  real, parameter :: Lx = 14.85, Ly = 4.9 !TAMAÑO DEL DOMINIO
  real, parameter :: dx = Lx / (Nx-1), dy = Ly / (Ny-1)
  real, parameter :: dt = 0.0003
  real, parameter :: rho = 1.0            ! densidad de referencia (usada en Poisson)
  real, parameter :: nu = 0.01            !Viscosidad cinemática (m²/s) deberia variar con la T
  real, parameter :: U_inf = 2.0          ! se rompe en 5 o mayor a 2 velocidad alta del oxidante
  integer, parameter :: MAX_ITER_P = 200

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
  ! Parámetros geométricos del cubo
  real, parameter :: x_izq_cubo = 5.425        ! borde izquierdo original
  !real, parameter :: x_inicio = 5.0 !no funciona agregar parámetro para definir el inicio del cubo en x
  !real, parameter :: cubo_centro_x = Lx/2.0 funciona para el primero
  real, parameter :: cubo_centro_x = (x_izq_cubo + Lx) / 2.0
   real, parameter :: cubo_mitad_ancho = (Lx - x_izq_cubo) / 2.0
   real, parameter :: cubo_centro_y = 1.0      ! base en y=0
  !real, parameter :: cubo_centro_y = 1.0      ! base en y=0
  !real, parameter :: cubo_centro_y =(x_inicio + Lx) / 2.0 !no funciona
  !real, parameter :: cubo_mitad_ancho = 2.0
 ! real, parameter :: cubo_mitad_ancho = (Lx - x_inicio) / 2.0
  real :: cubo_mitad_alto = 1.0               ! se reducirá con la regresión
  real, parameter :: cubo_alto_inicial = 1.0
  ! Parámetros del metano sólido
  !el rho_solido es igual al rho_cubo, que es el del metano sólido, y cp_solido es igual al cp_cubo
  real, parameter :: rho_solido = 100.0       ! kg/m³
  real, parameter :: cp_solido  = 2000.0      ! J/kg·K
  real, parameter :: T_sub      = 150.0       ! K (idealizado)
  real, parameter :: L_sub      = 50000.0    ! J/kg !antes 500000 calor latenete de sublimacion 
  real :: r_reg                               ! velocidad de regresión (m/s) - ahora variable

  logical :: es_obstaculo(n_nodos)
  logical :: es_obstaculo_anterior(n_nodos)  ! para detectar nodos que pasan de sólido a fluido
  real :: T_cubo
  real, parameter :: T_ambiente = 300.0
  real, parameter :: rho_cubo = 100.0, cp_cubo = 2000.0   ! metano sólido
  real :: A_cubo                                 ! área expuesta (variable)
  real :: masa_cubo
  real :: flujo_calor_total                       ! W total hacia el cubo
  real :: tasa_sublimacion                        ! kg/s total desde el cubo
  real :: v_inyeccion                             ! velocidad de inyección calculada (m/s)

  ! Variables de especies y termodinámica
  real :: T(n_nodos), T_new(n_nodos)
  real :: YF(n_nodos), YF_new(n_nodos)      ! CH4
  real :: YO(n_nodos), YO_new(n_nodos)      ! O2
  real :: YCO2(n_nodos), YCO2_new(n_nodos)  ! CO2
  real :: YH2O(n_nodos), YH2O_new(n_nodos)  ! H2O
  real :: H(n_nodos), H_new(n_nodos)        ! entalpía total de la mezcla (J/kg)
  real, parameter :: D_thermal = 0.05 !difusividad térmica efectiva
  real, parameter :: D_species = 0.05 !difusividad másica efectiva:
  real, parameter :: tau_comb = 0.005 !tiempo característico de reacción:
  real :: omega                       !tasa de reacción másica

  real, parameter :: beta = 1.0 / T_ambiente !coeficiente de expansión térmica (1/K) para Boussinesq
  real, parameter :: g_acel = -9.81
  real :: fuente_flotacion                   !érmino fuente de flotación en la ecuación de momento vertical:

  real :: grad_T_x(n_nodos), grad_T_y(n_nodos)
  real :: grad_YF_x(n_nodos), grad_YF_y(n_nodos)
  real :: grad_YO_x(n_nodos), grad_YO_y(n_nodos)
  real :: grad_YCO2_x(n_nodos), grad_YCO2_y(n_nodos)
  real :: grad_YH2O_x(n_nodos), grad_YH2O_y(n_nodos)
  real :: grad_H_x(n_nodos), grad_H_y(n_nodos)
  real :: lap_T(n_nodos), lap_YF(n_nodos), lap_YO(n_nodos)
  real :: lap_YCO2(n_nodos), lap_YH2O(n_nodos)
  real :: adv_T(n_nodos), adv_YF(n_nodos), adv_YO(n_nodos)
  real :: adv_YCO2(n_nodos), adv_YH2O(n_nodos)
  real :: adv_H(n_nodos), lap_H(n_nodos)
  integer :: idx_vecino, k_count

  real, parameter :: DT_MAX = 50.0 !Temperatura maxima en los pasos de tiempo.

  ! ============ VARIABLES PARA TERMODINÁMICA ============
  integer, parameter :: n_tabla = 9601
  real :: T_tabla(n_tabla)
  real :: cp_O2_tabla(n_tabla), h_O2_tabla(n_tabla), s_O2_tabla(n_tabla)
  real :: cp_CH4_tabla(n_tabla), h_CH4_tabla(n_tabla), s_CH4_tabla(n_tabla)
  real :: cp_CO2_tabla(n_tabla), h_CO2_tabla(n_tabla), s_CO2_tabla(n_tabla)
  real :: cp_H2O_tabla(n_tabla), h_H2O_tabla(n_tabla), s_H2O_tabla(n_tabla)

  ! Entalpías de formación (J/kg) para calor de reacción================peso molecular de cada compuesto
  real, parameter :: hf_CH4 = -74600.0 / 16.04246 * 1000.0
  real, parameter :: hf_O2  = 0.0
  real, parameter :: hf_CO2 = -393510.0 / 44.00950 * 1000.0
  real, parameter :: hf_H2O = -241826.0 / 18.01528 * 1000.0
  !pesos moleculares
  real, parameter :: W_CH4 = 16.04246
  real, parameter :: W_O2  = 31.99880
  real, parameter :: W_CO2 = 44.00950
  real, parameter :: W_H2O = 18.01528
  real, parameter :: Ru = 8.31446261815324
  real, parameter :: p_atm = 101325.0   ! presión de referencia para densidad

  ! Estequiometría másica corregida no se xd
  real, parameter :: s     = 2.0 * W_O2  / W_CH4   ! ≈ 3.989
  real, parameter :: r_CO2 = W_CO2 / W_CH4         ! ≈ 2.743
  real, parameter :: r_H2O = 2.0 * W_H2O / W_CH4   ! ≈ 2.246

  ! Variables para densidad local
  real :: rho_local(n_nodos)
  real :: lambda_ef = 1.0   ! con 5 funciona bien conductividad térmica efectiva (W/m·K) antes 0.1

  ! Variables escalares (movidas fuera del bucle)
  real :: suma_Yk_Wk
  real :: cp_CH4, cp_O2, cp_CO2, cp_H2O, cp_mezcla
  real :: calor_reaccion
  ! ========================================================

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

  ! Definición inicial del cubo (metano sólido)===========condiciones 
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

  T(:) = T_ambiente !Parece que T(:) deberia ser diferente 
  YF(:) = 0.0
  YO(:) = 1.0
  YCO2(:) = 0.0
  YH2O(:) = 0.0
  H(:) = 0.0
  T_cubo = 300.0!T_ambiente
  masa_cubo = rho_solido * (2.0*cubo_mitad_ancho) * (2.0*cubo_mitad_alto)

  do i = 1, n_nodos
      if (es_obstaculo(i)) then
          u(i) = 0.0; v(i) = 0.0
      end if
  end do

  print *, "Malla generada:", n_nodos, "nodos,", n_triangulos, "triángulos"
  print *, "Nodos dentro del cubo:", count(es_obstaculo)

  ! ============ LECTURA DE TABLAS TERMODINÁMICAS ============
  call leer_tabla('O2.txt',  T_tabla, cp_O2_tabla,  h_O2_tabla,  s_O2_tabla,  n_tabla)
  call leer_tabla('CH4.txt', T_tabla, cp_CH4_tabla, h_CH4_tabla, s_CH4_tabla, n_tabla)
  call leer_tabla('CO2.txt', T_tabla, cp_CO2_tabla, h_CO2_tabla, s_CO2_tabla, n_tabla)
  call leer_tabla('H2O.txt', T_tabla, cp_H2O_tabla, h_H2O_tabla, s_H2O_tabla, n_tabla)
  print *, "Tablas termodinámicas cargadas."

  ! Calcular entalpía inicial (ahora usando h de las tablas,cuidado aca)
  call calcular_entalpia_mezcla(T, YF, YO, YCO2, YH2O, n_nodos, H, &
                                h_CH4_tabla, h_O2_tabla, h_CO2_tabla, h_H2O_tabla, &
                                T_tabla, n_tabla)

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
  do iter = 1, 600000 !llega hasta 90000 despues de 3 horas. 

    ! ============ SUBLIMACIÓN REAL ============
    flujo_calor_total = 0.0
    do i = 1, n_nodos
      if (es_obstaculo(i)) then
        do idx_vecino = 1, n_vecinos(i)
          j = vecinos(i, idx_vecino)
          if (.not. es_obstaculo(j)) then
            flujo_calor_total = flujo_calor_total + lambda_ef * (T(j) - T(i)) / sqrt(dx**2 + dy**2)
          end if
        end do
      end if
    end do

    tasa_sublimacion = flujo_calor_total / L_sub
    if (tasa_sublimacion < 0.0) tasa_sublimacion = 0.0

    A_cubo = 2.0 * cubo_mitad_ancho
    r_reg = tasa_sublimacion / (rho_solido * A_cubo)

    cubo_mitad_alto = cubo_mitad_alto - r_reg * dt
    if (cubo_mitad_alto < 0.01) cubo_mitad_alto = 0.01

    masa_cubo = rho_solido * (2.0*cubo_mitad_ancho) * (2.0*cubo_mitad_alto)
    if (masa_cubo < 1.0e-6) masa_cubo = 0.0

    es_obstaculo_anterior = es_obstaculo

    es_obstaculo = .false.
    do i = 1, n_nodos
        if (x(i) >= cubo_centro_x - cubo_mitad_ancho .and. x(i) <= cubo_centro_x + cubo_mitad_ancho .and. &
            y(i) >= cubo_centro_y - cubo_mitad_alto  .and. y(i) <= cubo_centro_y + cubo_mitad_alto) then
            es_obstaculo(i) = .true.
        end if
    end do

    do i = 1, n_nodos
      if (es_obstaculo_anterior(i) .and. .not. es_obstaculo(i)) then
        T(i) = T_ambiente
        YF(i) = 0.0; YO(i) = 1.0; YCO2(i) = 0.0; YH2O(i) = 0.0
        call interpolar_tabla(T_ambiente, T_tabla, h_O2_tabla, n_tabla, H(i))
      end if
    end do

    if (tasa_sublimacion > 0.0 .and. masa_cubo > 0.0) then
      v_inyeccion = tasa_sublimacion / (rho * A_cubo)
    else
      v_inyeccion = 0.0
    end if
    v_inyeccion = min(v_inyeccion, 0.5 * min(dx, dy) / dt)

    ! - TRANSPORTE DE ESCALARES Y ENERGÍA
    call calcular_gradiente(T,  grad_T_x,  grad_T_y,  n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YF, grad_YF_x, grad_YF_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YO, grad_YO_x, grad_YO_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YCO2, grad_YCO2_x, grad_YCO2_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(YH2O, grad_YH2O_x, grad_YH2O_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(H,  grad_H_x,  grad_H_y,  n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)

    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            adv_T(i)  = 0.0; lap_T(i)  = 0.0
            adv_YF(i) = 0.0; lap_YF(i) = 0.0
            adv_YO(i) = 0.0; lap_YO(i) = 0.0
            adv_YCO2(i)= 0.0; lap_YCO2(i)= 0.0
            adv_YH2O(i)= 0.0; lap_YH2O(i)= 0.0
            adv_H(i)  = 0.0; lap_H(i)  = 0.0
            cycle
        end if

        adv_T(i)  = u(i)*grad_T_x(i)  + v(i)*grad_T_y(i)
        adv_YF(i) = u(i)*grad_YF_x(i) + v(i)*grad_YF_y(i)
        adv_YO(i) = u(i)*grad_YO_x(i) + v(i)*grad_YO_y(i)
        adv_YCO2(i)= u(i)*grad_YCO2_x(i) + v(i)*grad_YCO2_y(i)
        adv_YH2O(i)= u(i)*grad_YH2O_x(i) + v(i)*grad_YH2O_y(i)
        adv_H(i)  = u(i)*grad_H_x(i)  + v(i)*grad_H_y(i)

        lap_T(i)  = 0.0; lap_YF(i) = 0.0; lap_YO(i) = 0.0
        lap_YCO2(i) = 0.0; lap_YH2O(i) = 0.0; lap_H(i) = 0.0
        do idx_vecino = 1, n_vecinos(i)
            j = vecinos(i, idx_vecino)
            lap_T(i)  = lap_T(i)  + S_val(i, idx_vecino) * (T(j)  - T(i))
            lap_YF(i) = lap_YF(i) + S_val(i, idx_vecino) * (YF(j) - YF(i))
            lap_YO(i) = lap_YO(i) + S_val(i, idx_vecino) * (YO(j) - YO(i))
            lap_YCO2(i)= lap_YCO2(i) + S_val(i, idx_vecino) * (YCO2(j) - YCO2(i))
            lap_YH2O(i)= lap_YH2O(i) + S_val(i, idx_vecino) * (YH2O(j) - YH2O(i))
            lap_H(i)  = lap_H(i)  + S_val(i, idx_vecino) * (H(j)  - H(i))
        end do
        lap_T(i)  = -lap_T(i)  / M_lump(i)
        lap_YF(i) = -lap_YF(i) / M_lump(i)
        lap_YO(i) = -lap_YO(i) / M_lump(i)
        lap_YCO2(i)= -lap_YCO2(i) / M_lump(i)
        lap_YH2O(i)= -lap_YH2O(i) / M_lump(i)
        lap_H(i)  = -lap_H(i)  / M_lump(i)
    end do

    ! Inyección de metano gaseoso desde los nodos del cubo
    do i = 1, n_nodos
      if (es_obstaculo(i)) then
        do idx_vecino = 1, n_vecinos(i)
          j = vecinos(i, idx_vecino)
          if (.not. es_obstaculo(j)) then
            if (masa_cubo > 0.0 .and. tasa_sublimacion > 0.0) then
              YF(j) = 0.1
              YO(j) = 0.0
              YCO2(j) = 0.0
              YH2O(j) = 0.0
              v(j) = v_inyeccion
              u(j) = 0.0
            end if
          end if
        end do
      end if
    end do

    ! Tasa de reacción
    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            omega = 0.0
        else
            omega = min(YF(i), YO(i)/s) / tau_comb
            omega = min(omega, 0.1 / dt)
        end if

        suma_Yk_Wk = YF(i)/W_CH4 + YO(i)/W_O2 + YCO2(i)/W_CO2 + YH2O(i)/W_H2O
        if (suma_Yk_Wk < 1.0e-12) suma_Yk_Wk = 1.0/W_O2
        rho_local(i) = p_atm / (Ru * T(i) * suma_Yk_Wk)

        call interpolar_tabla(T(i), T_tabla, cp_CH4_tabla, n_tabla, cp_CH4)
        call interpolar_tabla(T(i), T_tabla, cp_O2_tabla,  n_tabla, cp_O2)
        call interpolar_tabla(T(i), T_tabla, cp_CO2_tabla, n_tabla, cp_CO2)
        call interpolar_tabla(T(i), T_tabla, cp_H2O_tabla, n_tabla, cp_H2O)
        cp_mezcla = YF(i)*cp_CH4 + YO(i)*cp_O2 + YCO2(i)*cp_CO2 + YH2O(i)*cp_H2O
        if (cp_mezcla < 1.0e-10) cp_mezcla = 1004.0

        calor_reaccion = (hf_CH4 + 2.0*hf_O2 - hf_CO2 - 2.0*hf_H2O) * omega

        YF_new(i)  = YF(i)  + dt * ( -adv_YF(i)  + D_species * lap_YF(i)  - omega )
        YO_new(i)  = YO(i)  + dt * ( -adv_YO(i)  + D_species * lap_YO(i)  - s * omega )
        YCO2_new(i)= YCO2(i)+ dt * ( -adv_YCO2(i)+ D_species * lap_YCO2(i)+ r_CO2 * omega )
        YH2O_new(i)= YH2O(i)+ dt * ( -adv_YH2O(i)+ D_species * lap_YH2O(i)+ r_H2O * omega )

        H_new(i) = H(i) + dt * ( -adv_H(i) + D_thermal * lap_H(i) + calor_reaccion )
    end do
 !aqui cambio cuidado aca
    if (mod(iter, 500) == 0) then
    print *, "max YF =", maxval(YF), " max YO =", maxval(YO)
end if
!aqui termina 
    ! Aplicar límites y actualizar
    do i = 1, n_nodos
        if (es_obstaculo(i)) then
            T(i)  = T_cubo
            YF(i) = 0.0
            YO(i) = 0.0
            YCO2(i) = 0.0
            YH2O(i) = 0.0
            H(i)  = 0.0
        else
            YF(i)  = max(0.0, min(1.0, YF_new(i)))
            YO(i)  = max(0.0, min(1.0, YO_new(i)))
            YCO2(i)= max(0.0, min(1.0, YCO2_new(i)))
            YH2O(i)= max(0.0, min(1.0, YH2O_new(i)))
            H(i)   = H_new(i)
        end if
    end do

    ! Convertir entalpía H a temperatura T usando tablas
    call entalpia_a_temperatura(H, YF, YO, YCO2, YH2O, n_nodos, T, T_tabla, &
                                cp_CH4_tabla, cp_O2_tabla, cp_CO2_tabla, cp_H2O_tabla, &
                                h_CH4_tabla, h_O2_tabla, h_CO2_tabla, h_H2O_tabla, &
                                n_tabla, T_ambiente)

    ! Calentamiento del cubo
    T_cubo = T_cubo + dt * (flujo_calor_total / (masa_cubo * cp_cubo))
    do i = 1, n_nodos
        if (es_obstaculo(i)) T(i) = T_cubo
    end do

    call condiciones_contorno_escalares(T, YF, YO, YCO2, YH2O, H, x, y, n_nodos, Lx, Ly, T_ambiente)

    ! Paso 5.1: Gradientes de velocidad
    call calcular_gradiente(u, grad_u_x, grad_u_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)
    call calcular_gradiente(v, grad_v_x, grad_v_y, n_nodos, n_triangulos, tri_nodos, x, y, area_tri, M_lump)

    do i = 1, n_nodos
      conv_u(i) = u(i) * grad_u_x(i) + v(i) * grad_u_y(i)
      conv_v(i) = u(i) * grad_v_x(i) + v(i) * grad_v_y(i)
    end do

    ! Paso 5.2: Predicción con flotación
    do i = 1, n_nodos
      if (es_obstaculo(i)) then
          u_star(i) = 0.0; v_star(i) = 0.0
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
        if (abs(x(i) - Lx) < 0.5*dx) then
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

    ! Detección de NaN Dice donde ocurrio el fallo 
    do i = 1, n_nodos
        if (u(i) /= u(i) .or. abs(u(i)) > 1.0e15) then
            print *, "Problema en u: nodo", i, "x,y=", x(i), y(i), "valor=", u(i), "paso", iter
            stop
        end if
        if (v(i) /= v(i) .or. abs(v(i)) > 1.0e15) then
            print *, "Problema en v: nodo", i, "x,y=", x(i), y(i), "valor=", v(i), "paso", iter
            stop
        end if
        if (p(i) /= p(i) .or. abs(p(i)) > 1.0e15) then
            print *, "Problema en p: nodo", i, "x,y=", x(i), y(i), "valor=", p(i), "paso", iter
            stop
        end if
        if (T(i) /= T(i) .or. abs(T(i)) > 1.0e8) then
            print *, "Problema en T: nodo", i, "x,y=", x(i), y(i), "T=", T(i), "H=", H(i), &
                     "YF,YO,YCO2,YH2O=", YF(i), YO(i), YCO2(i), YH2O(i), "paso", iter
            stop
        end if
        if (H(i) /= H(i) .or. abs(H(i)) > 1.0e15) then
            print *, "Problema en H: nodo", i, "x,y=", x(i), y(i), "H=", H(i), &
                     "es_obstaculo=", es_obstaculo(i), "paso", iter
            stop
        end if
        if (YF(i) /= YF(i) .or. YO(i) /= YO(i) .or. YCO2(i) /= YCO2(i) .or. YH2O(i) /= YH2O(i)) then
            print *, "Problema en especies: nodo", i, "x,y=", x(i), y(i), &
                     "YF,YO,YCO2,YH2O=", YF(i), YO(i), YCO2(i), YH2O(i), "paso", iter
            stop
        end if
    end do

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

      write(nombre_archivo, '(A,I4.4,A)') 'especies_', frame, '.dat'
      open(unit=14, file=nombre_archivo, status='replace')
      do i = 1, n_nodos
        write(14, '(5(ES14.6,1X))') x(i), y(i), YF(i), YO(i), YCO2(i)
      end do
      close(14)

      print *, "Paso ", iter, " T_cubo =", T_cubo, " masa_cubo =", masa_cubo, &
               " altura_cubo =", cubo_mitad_alto, " r_reg =", r_reg, " v_iny =", v_inyeccion
      print *, "  [diag] max|u|=", maxval(abs(u)), " max|v|=", maxval(abs(v)), &
               " max T=", maxval(T), " min T=", minval(T), &
               " max|H|=", maxval(abs(H)), " flujo_calor_total=", flujo_calor_total
    end if

  end do

  print *, "Simulación finalizada."

contains

  ! ============ SUBRUTINA: leer_tabla ============
  subroutine leer_tabla(nombre, T_tabla, cp_tabla, h_tabla, s_tabla, n)
    character(len=*), intent(in) :: nombre
    real, intent(out) :: T_tabla(n), cp_tabla(n), h_tabla(n), s_tabla(n)
    integer, intent(in) :: n
    integer :: unidad, i, io
    character(len=200) :: linea
    real :: dummy1, dummy2, dummy3

    open(newunit=unidad, file=nombre, status='old', action='read', iostat=io)
    if (io /= 0) then
      print *, "Error abriendo tabla ", nombre
      stop
    end if

    do
      read(unidad, '(A)', iostat=io) linea
      if (io /= 0) exit
      if (len_trim(linea) > 0) then
        if (linea(1:1) /= '#') exit
      end if
    end do

    backspace(unidad)
    do i = 1, n
      read(unidad, *, iostat=io) T_tabla(i), cp_tabla(i), h_tabla(i), s_tabla(i), &
                                  dummy1, dummy2, dummy3
      if (io /= 0) exit
    end do
    close(unidad)
  end subroutine leer_tabla

  ! ============ SUBRUTINA: interpolar_tabla ============
  subroutine interpolar_tabla(T, T_tabla, prop_tabla, n, prop)
    real, intent(in) :: T, T_tabla(n), prop_tabla(n)
    integer, intent(in) :: n
    real, intent(out) :: prop
    integer :: lo, hi, mid
    real :: frac

    if (T <= T_tabla(1)) then
      prop = prop_tabla(1)
      return
    end if
    if (T >= T_tabla(n)) then
      prop = prop_tabla(n)
      return
    end if

    lo = 1
    hi = n
    do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (T_tabla(mid) <= T) then
        lo = mid
      else
        hi = mid
      end if
    end do

    if (abs(T_tabla(hi) - T_tabla(lo)) > 1.0e-8) then
      frac = (T - T_tabla(lo)) / (T_tabla(hi) - T_tabla(lo))
      prop = prop_tabla(lo) + frac * (prop_tabla(hi) - prop_tabla(lo))
    else
      prop = prop_tabla(lo)
    end if
  end subroutine interpolar_tabla

  ! ============ SUBRUTINA: calcular_entalpia_mezcla ============
  ! Ahora interpola directamente la entalpía de las tablas.
  subroutine calcular_entalpia_mezcla(T, YF, YO, YCO2, YH2O, n_nodos, H, &
                                      h_CH4_tabla, h_O2_tabla, h_CO2_tabla, h_H2O_tabla, &
                                      T_tabla, n_tabla)
    integer, intent(in) :: n_nodos, n_tabla
    real, intent(in) :: T(n_nodos), YF(n_nodos), YO(n_nodos), YCO2(n_nodos), YH2O(n_nodos)
    real, intent(out) :: H(n_nodos)
    real, intent(in) :: h_CH4_tabla(n_tabla), h_O2_tabla(n_tabla), &
                        h_CO2_tabla(n_tabla), h_H2O_tabla(n_tabla)
    real, intent(in) :: T_tabla(n_tabla)
    integer :: i
    real :: h_CH4, h_O2, h_CO2, h_H2O

    do i = 1, n_nodos
      call interpolar_tabla(T(i), T_tabla, h_CH4_tabla, n_tabla, h_CH4)
      call interpolar_tabla(T(i), T_tabla, h_O2_tabla,  n_tabla, h_O2)
      call interpolar_tabla(T(i), T_tabla, h_CO2_tabla, n_tabla, h_CO2)
      call interpolar_tabla(T(i), T_tabla, h_H2O_tabla, n_tabla, h_H2O)
      H(i) = YF(i)*h_CH4 + YO(i)*h_O2 + YCO2(i)*h_CO2 + YH2O(i)*h_H2O
    end do
  end subroutine calcular_entalpia_mezcla

  ! ============ SUBRUTINA: entalpia_a_temperatura ============
  subroutine entalpia_a_temperatura(H, YF, YO, YCO2, YH2O, n_nodos, T, T_tabla, &
                                    cp_CH4_tabla, cp_O2_tabla, cp_CO2_tabla, cp_H2O_tabla, &
                                    h_CH4_tabla, h_O2_tabla, h_CO2_tabla, h_H2O_tabla, &
                                    n_tabla, T_ambiente)
    integer, intent(in) :: n_nodos, n_tabla
    real, intent(inout) :: H(n_nodos)
    real, intent(in) :: YF(n_nodos), YO(n_nodos), YCO2(n_nodos), YH2O(n_nodos)
    real, intent(inout) :: T(n_nodos)
    real, intent(in) :: T_tabla(n_tabla)
    real, intent(in) :: cp_CH4_tabla(n_tabla), cp_O2_tabla(n_tabla), cp_CO2_tabla(n_tabla), cp_H2O_tabla(n_tabla)
    real, intent(in) :: h_CH4_tabla(n_tabla), h_O2_tabla(n_tabla), h_CO2_tabla(n_tabla), h_H2O_tabla(n_tabla)
    real, intent(in) :: T_ambiente
    integer :: i, lo, hi, mid
    real :: h_lo, h_hi, h_mid, frac, h_target

    do i = 1, n_nodos
      h_target = H(i)
      if (h_target < -1.0e10) then
        T(i) = T_ambiente
        cycle
      end if

      h_lo = YF(i)*h_CH4_tabla(1) + YO(i)*h_O2_tabla(1) + YCO2(i)*h_CO2_tabla(1) + YH2O(i)*h_H2O_tabla(1)
      h_hi = YF(i)*h_CH4_tabla(n_tabla) + YO(i)*h_O2_tabla(n_tabla) + YCO2(i)*h_CO2_tabla(n_tabla) + YH2O(i)*h_H2O_tabla(n_tabla)

      if (h_target <= h_lo) then
        T(i) = T_tabla(1)
        H(i) = h_lo
        cycle
      end if
      if (h_target >= h_hi) then
        T(i) = T_tabla(n_tabla)
        H(i) = h_hi
        cycle
      end if

      lo = 1
      hi = n_tabla
      do while (hi - lo > 1)
        mid = (lo + hi) / 2
        h_mid = YF(i)*h_CH4_tabla(mid) + YO(i)*h_O2_tabla(mid) + YCO2(i)*h_CO2_tabla(mid) + YH2O(i)*h_H2O_tabla(mid)
        if (h_mid <= h_target) then
          lo = mid
        else
          hi = mid
        end if
      end do

      h_lo = YF(i)*h_CH4_tabla(lo) + YO(i)*h_O2_tabla(lo) + YCO2(i)*h_CO2_tabla(lo) + YH2O(i)*h_H2O_tabla(lo)
      h_hi = YF(i)*h_CH4_tabla(hi) + YO(i)*h_O2_tabla(hi) + YCO2(i)*h_CO2_tabla(hi) + YH2O(i)*h_H2O_tabla(hi)
      if (abs(h_hi - h_lo) > 1.0e-8) then
        frac = (h_target - h_lo) / (h_hi - h_lo)
      else
        frac = 0.0
      end if
      T(i) = T_tabla(lo) + frac * (T_tabla(hi) - T_tabla(lo))
    end do
  end subroutine entalpia_a_temperatura

  ! ============ SUBRUTINAS ORIGINALES ============
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
      if (abs(x(i)) < 0.5*dx) then
        u(i) = U_inf; v(i) = 0.0
      end if
      if (abs(y(i)) < 0.5*dy .or. abs(y(i) - Ly) < 0.5*dy) then
        u(i) = 0.0; v(i) = 0.0
      end if
      if (abs(x(i) - Lx) < 0.5*dx) then
        found = .false. 
        do v_idx = 1, n_vecinos(i)
          j = vecinos(i, v_idx)
          if (x(j) < Lx - 0.5*dx) then
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

  subroutine condiciones_contorno_escalares(T, YF, YO, YCO2, YH2O, H, x, y, n_nodos, Lx, Ly, T_amb)
    integer, intent(in) :: n_nodos
    real, intent(inout) :: T(n_nodos), YF(n_nodos), YO(n_nodos), YCO2(n_nodos), YH2O(n_nodos), H(n_nodos)
    real, intent(in)    :: x(n_nodos), y(n_nodos)
    real, intent(in)    :: Lx, Ly, T_amb
    integer :: i
    real :: h_O2_in

    do i = 1, n_nodos
      if (abs(x(i)) < 0.5*dx) then
          T(i) = T_amb
          YF(i) = 0.0
          YO(i) = 1.0
          YCO2(i) = 0.0
          YH2O(i) = 0.0
!==========================================================================
          ! Usar directamente la entalpía de la tabla para O2 puro
          call interpolar_tabla(T_amb, T_tabla, h_O2_tabla, n_tabla, h_O2_in)
          H(i) = YO(i) * h_O2_in
      end if
      if (abs(y(i)) < 0.5*dy .or. abs(y(i) - Ly) < 0.5*dy) then
          continue
      end if
      if (abs(x(i) - Lx) < 0.5*dx) then
          continue
      end if
    end do
  end subroutine condiciones_contorno_escalares

end program navier_stokes_python_triangulos
