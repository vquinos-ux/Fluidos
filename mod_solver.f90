module mod_solver
  use mod_precision
  use mod_termo
  use mod_cinetica
  use mod_malla
  implicit none

  integer  :: diag_k_poisson  = 0                                     ! barridos Poisson [Ec.48]
  logical  :: diag_poisson_ok = .true.                                ! convergio Poisson
  real(dp) :: diag_res_p      = 0.0_dp                                ! residuo [Ec.48]
  real(dp) :: diag_cont_ant   = 0.0_dp                                ! cont antes [Ec.48]
  real(dp) :: diag_cont_des   = 0.0_dp                                ! cont despues [Ec.48]
  integer  :: diag_i_cont = 0, diag_j_cont = 0                        ! pos max cont
  real(dp) :: diag_umax = 0.0_dp                                      ! max |U|
  integer  :: diag_i_umax = 0, diag_j_umax = 0                        ! pos max |U|
  integer  :: un_diag = 0                                             ! unidad log
  logical  :: diag_log_abierto = .false.                              ! log abierto
  logical  :: diag_armado = .true.                                    ! rearma avisos
  integer  :: diag_n_avisos = 0                                       ! contador avisos
  real(dp), parameter :: U_AVISO = 1000.0_dp                          ! umbral |U| [m/s] [Ec.30]
  integer,  parameter :: N_DIAG  = 10                                 ! iter entre lineas

  logical, parameter :: MOMENTO_CONSISTENTE = .true.                  ! rho^n+1 u^n [Ec.32]
  real(dp), parameter :: EPS4 = 1.0_dp/64.0_dp                        ! filtro par-impar (sin ec)
  logical, parameter :: SIN_REFLUJO = .true.                          ! anti-reflujo salida

contains

  function masa_combustible() result(m)                               ! masa grano [g]
    implicit none
    real(dp) :: m, dxi
    integer  :: k
    m = 0.0_dp
    do k = 1, NX
       if (.not. es_grano(k)) cycle
       if (k == 1) then
          dxi = 0.5_dp*(xg(2) - xg(1))                                ! dx borde izq
       else if (k == NX) then
          dxi = 0.5_dp*(xg(NX) - xg(NX-1))                            ! dx borde der
       else
          dxi = 0.5_dp*(xg(k+1) - xg(k-1))                            ! dx interior
       end if
       m = m + max(esp(k) - e_pared, 0.0_dp)*dxi                      ! integral e-e_pared
    end do
    m = m*rho_s*ancho_z*1000.0_dp                                     ! kg -> g
  end function masa_combustible

  subroutine mezcla_cp_h(Y, T, cpm, hm)                               ! [Ec.6]
    implicit none
    real(dp), intent(in)  :: Y(N_ESP), T
    real(dp), intent(out) :: cpm, hm
    cpm = cp_mezcla(T, Y)                                             ! [Ec.6] cp_mix
    hm  = entalpia_mezcla(T, Y)                                       ! [Ec.6] h_mix
  end subroutine mezcla_cp_h

  subroutine construir_malla()                                        ! [Ec.37, 38]
    implicit none
    integer :: i, j, jg
    do i = 1, NX
       do j = 1, NYS
          yg(i,j) = esp(i)*fsol(j)                                    ! [Ec.37] solido
       end do
       do jg = 2, NYG
          j = NYS + jg - 1
          yg(i,j) = esp(i) + (Ly - esp(i))*fgas(jg)                   ! [Ec.37] gas
       end do
    end do
    do j = 1, NY
       dYdx(1,j)  = (yg(2,j)-yg(1,j))/(xg(2)-xg(1))                   ! [Ec.38] metrica borde
       dYdx(NX,j) = (yg(NX,j)-yg(NX-1,j))/(xg(NX)-xg(NX-1))           ! [Ec.38]
       do i = 2, NX-1
          dYdx(i,j) = (yg(i+1,j)-yg(i-1,j))/(xg(i+1)-xg(i-1))         ! [Ec.38] M = dY/dx
       end do
    end do
  end subroutine construir_malla

  subroutine actualizar_entrada()                                     ! U_ent, G_ox segun modo
    implicit none
    real(dp) :: A_in, h_abierta
    cos_ent = cos(ang_entrada*PI_M/180.0_dp)                          ! cos ang
    sin_ent = sin(ang_entrada*PI_M/180.0_dp)                          ! sin ang
    if (usar_ranura) then
       h_abierta = max(y_ranura - esp(1), 1.0e-4_dp)                  ! ranura
    else
       h_abierta = max(Ly - esp(1), 1.0e-4_dp)                        ! hueco completo
    end if
    A_in   = h_abierta*ancho_z                                        ! area paso
    rho_in = P_OPER*W_ESP(I_O2)/(RU_KMOL*T_in)                        ! [Ec.9] rho O2
    select case (modo_entrada)
    case (1)
       G_ox  = mdot_ox/A_in                                           ! mdot fijo
       U_ent = G_ox/rho_in
    case (3)
       G_ox  = G_ox_fijo                                              ! G_ox fijo
       U_ent = G_ox/rho_in
    case default
       U_ent = U_in_fijo                                              ! U fija
       G_ox  = rho_in*U_ent
    end select
  end subroutine actualizar_entrada

  pure logical function entrada_abierta(j)                            ! ranura: fila abierta?
    implicit none
    integer, intent(in) :: j
    if (usar_ranura) then
       entrada_abierta = (yg(1,j) <= y_ranura)                        ! bajo techo ranura
    else
       entrada_abierta = .true.                                       ! todo abierto
    end if
  end function entrada_abierta

  subroutine inicializar()                                            ! malla + estado inicial
    implicit none
    integer  :: i, j, jg
    real(dp) :: t, cpm, hm
    do i = 1, NX
       t = real(i-1,dp)/real(NX-1,dp)
       xg(i) = Lx*t**pot_x                                            ! x concentrado
    end do
    do j = 1, NYS
       t = real(j-1,dp)/real(NYS-1,dp)
       fsol(j) = 1.0_dp - (exp(b_sol*(1.0_dp-t))-1.0_dp)/(exp(b_sol)-1.0_dp)  ! perfil solido
    end do
    fsol(1) = 0.0_dp;  fsol(NYS) = 1.0_dp
    do jg = 1, NYG
       t = real(jg-1,dp)/real(NYG-1,dp)
       fgas(jg) = (exp(b_gas*t)-1.0_dp)/(exp(b_gas)-1.0_dp)           ! perfil gas
    end do
    fgas(1) = 0.0_dp;  fgas(NYG) = 1.0_dp
    do i = 1, NX
       esp(i)      = perfil_espesor(xg(i))                            ! [Ec.36]
       es_grano(i) = (xg(i) >= x_ini_gr .and. xg(i) <= x_fin_gr)      ! mascara grano
    end do
    agotado = .false.
    call construir_malla()                                            ! [Ec.37,38]
    call actualizar_entrada()

    Ysp = 0.0_dp
    Ysp(:,:,I_O2) = 1.0_dp                                            ! O2 puro
    Tg = T_in                                                         ! T entrada
    u  = 0.0_dp;  v = 0.0_dp;  pr = 0.0_dp
    do j = NYS+1, NY
       u(:,j) = U_ent*cos_ent                                         ! semilla axial
    end do
    do j = 1, NY
       do i = 1, NX
          call mezcla_cp_h(Ysp(i,j,:), Tg(i,j), cpm, hm)              ! [Ec.6]
          hmix(i,j) = hm
       end do
    end do
    call propiedades()                                                ! [Ec.9-14]
    ru_ = rho*u;  rv_ = rho*v                                         ! [Ec.32]
    rho_o = rho

    Ts = T_amb                                                        ! solido T_amb
    Tw = T_amb
    rdot = 0.0_dp;  mdot = 0.0_dp;  qconv = 0.0_dp;  qrad = 0.0_dp
    qvol = 0.0_dp
  end subroutine inicializar

  subroutine propiedades()                                            ! [Ec.9-14] por celda gas
    implicit none
    integer  :: i, j, s
    real(dp) :: Yl(N_ESP), D_k(N_ESP)
    real(dp) :: sumYW, Tl, mu_mix, lam_mix, cp_mix
    !$omp parallel do default(shared) &
    !$omp   private(i,j,s,Yl,D_k,sumYW,Tl,mu_mix,lam_mix,cp_mix)
    do j = NYS, NY
       do i = 1, NX
          Tl = min(max(Tg(i,j), T_TAB_MIN), T_TAB_MAX)                ! clip T
          Yl = Ysp(i,j,:)
          do s = 1, N_ESP
             if (Yl(s) < 0.0_dp) Yl(s) = 0.0_dp
          end do
          sumYW = 0.0_dp
          do s = 1, N_ESP
             sumYW = sumYW + Yl(s)/W_ESP(s)                           ! sum Yk/Wk [Ec.8]
          end do
          if (sumYW < 1.0e-12_dp) sumYW = 1.0_dp/W_ESP(I_O2)
          rho(i,j) = P_OPER/(RU_KMOL*Tl*sumYW)                        ! [Ec.9] rho
          call propiedades_mezcla(Tl, Yl, mu_mix, lam_mix, cp_mix, D_k)  ! [Ec.10-14,18]
          if (cp_mix < 200.0_dp)  cp_mix  = 200.0_dp                  ! proteccion
          if (mu_mix < 1.0e-7_dp) mu_mix  = 1.0e-7_dp
          if (lam_mix< 1.0e-4_dp) lam_mix = 1.0e-4_dp
          cp_(i,j)  = cp_mix                                          ! [Ec.6]
          mu_(i,j)  = mu_mix                                          ! [Ec.13]
          lam_(i,j) = lam_mix                                         ! [Ec.14]
          Dm_(i,j) = lam_(i,j)/cp_(i,j)                               ! rho D = lam/cp (Lewis)
       end do
    end do
    !$omp end parallel do
  end subroutine propiedades

  subroutine paso_tiempo()                                            ! [Ec.50-52]
    implicit none
    integer  :: i, j
    real(dp) :: dtc, dtd, dxm, dym, vrel, gmx
    dtc = dt_max;  dtd = dt_max
    do j = NYS+1, NY
       do i = 2, NX-1
          dxm = 0.5_dp*(xg(i+1)-xg(i-1))                              ! dx centrado
          if (j == NY) then
             dym = yg(i,NY)-yg(i,NY-1)                                ! dy borde sup
          else
             dym = 0.5_dp*(yg(i,j+1)-yg(i,j-1))                       ! dy centrado
          end if
          vrel = abs(v(i,j) - u(i,j)*dYdx(i,j))                       ! [Ec.39] v_r
          if (abs(u(i,j)) > 1.0e-9_dp) dtc = min(dtc, CFL*dxm/abs(u(i,j)))  ! [Ec.50] dx
          if (vrel > 1.0e-9_dp)        dtc = min(dtc, CFL*dym/vrel)   ! [Ec.50] dy
          gmx = max(mu_(i,j), lam_(i,j)/cp_(i,j), Dm_(i,j))/rho(i,j)  ! gamma [Ec.51]
          dtd = min(dtd, FOU*dxm*dxm/gmx)                             ! [Ec.51]
       end do
    end do
    dt = min(dtc, dtd, dt_max)                                        ! [Ec.52]
  end subroutine paso_tiempo

  pure function d_dx_g(f,i,j) result(r)                               ! d/dx a j=cte
    implicit none
    real(dp), intent(in) :: f(NX,NY)
    integer,  intent(in) :: i, j
    real(dp) :: r
    if (i == 1) then
       r = (f(2,j)-f(1,j))/(xg(2)-xg(1))                              ! borde izq
    else if (i == NX) then
       r = (f(NX,j)-f(NX-1,j))/(xg(NX)-xg(NX-1))                      ! borde der
    else
       r = (f(i+1,j)-f(i-1,j))/(xg(i+1)-xg(i-1))                      ! centrada
    end if
  end function d_dx_g

  pure function d_dy_g(f,i,j) result(r)                               ! d/dy fisico
    implicit none
    real(dp), intent(in) :: f(NX,NY)
    integer,  intent(in) :: i, j
    real(dp) :: r
    if (j == NYS) then
       r = (f(i,NYS+1)-f(i,NYS))/(yg(i,NYS+1)-yg(i,NYS))              ! borde pared
    else if (j == NY) then
       r = (f(i,NY)-f(i,NY-1))/(yg(i,NY)-yg(i,NY-1))                  ! borde techo
    else
       r = (f(i,j+1)-f(i,j-1))/(yg(i,j+1)-yg(i,j-1))                  ! centrada
    end if
  end function d_dy_g

  pure function d_dx_fis(f,i,j) result(r)                             ! d/dx fisico a y=cte
    implicit none
    real(dp), intent(in) :: f(NX,NY)
    integer,  intent(in) :: i, j
    real(dp) :: r
    r = d_dx_g(f,i,j) - dYdx(i,j)*d_dy_g(f,i,j)                       ! [Ec.38] metrica
  end function d_dx_fis

  pure function adv_up(f,i,j,uu,vr) result(r)                         ! upwind 1er orden
    implicit none
    real(dp), intent(in) :: f(NX,NY), uu, vr
    integer,  intent(in) :: i, j
    real(dp) :: r, dfx, dfy
    if (uu >= 0.0_dp) then
       if (i > 1) then;  dfx = (f(i,j)-f(i-1,j))/(xg(i)-xg(i-1))
       else;             dfx = (f(2,j)-f(1,j))/(xg(2)-xg(1));  end if
    else
       if (i < NX) then; dfx = (f(i+1,j)-f(i,j))/(xg(i+1)-xg(i))
       else;             dfx = (f(NX,j)-f(NX-1,j))/(xg(NX)-xg(NX-1));  end if
    end if
    if (vr >= 0.0_dp) then                                            ! [Ec.39] v_r
       if (j > NYS) then; dfy = (f(i,j)-f(i,j-1))/(yg(i,j)-yg(i,j-1))
       else;              dfy = (f(i,NYS+1)-f(i,NYS))/(yg(i,NYS+1)-yg(i,NYS));  end if
    else
       if (j < NY) then;  dfy = (f(i,j+1)-f(i,j))/(yg(i,j+1)-yg(i,j))
       else;              dfy = (f(i,NY)-f(i,NY-1))/(yg(i,NY)-yg(i,NY-1));  end if
    end if
    r = uu*dfx + vr*dfy
  end function adv_up

  pure function difx(f,G,i,j) result(r)                               ! d/dx(G df/dx)
    implicit none
    real(dp), intent(in) :: f(NX,NY), G(NX,NY)
    integer,  intent(in) :: i, j
    real(dp) :: r, dxp, dxm, dxc, Gp, Gm
    if (i <= 1 .or. i >= NX) then
       r = 0.0_dp;  return                                            ! borde: 0
    end if
    dxp = xg(i+1)-xg(i);  dxm = xg(i)-xg(i-1);  dxc = 0.5_dp*(dxp+dxm)  ! pasos
    Gp = 0.5_dp*(G(i+1,j)+G(i,j));  Gm = 0.5_dp*(G(i,j)+G(i-1,j))     ! G cara
    r = (Gp*(f(i+1,j)-f(i,j))/dxp - Gm*(f(i,j)-f(i-1,j))/dxm)/dxc     ! flujo neto
  end function difx

  subroutine thomas(a,b,c,d,n,x)                                      ! [Ec.46,47]
    implicit none
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a(n), b(n), c(n), d(n)
    real(dp), intent(out) :: x(n)
    integer  :: k
    real(dp) :: cp2(n), dp2(n), m
    cp2(1) = c(1)/b(1);  dp2(1) = d(1)/b(1)                           ! pivote inicial
    do k = 2, n
       m = b(k) - a(k)*cp2(k-1)                                       ! [Ec.47] denominador
       if (abs(m) < 1.0e-300_dp) m = 1.0e-300_dp
       cp2(k) = c(k)/m
       dp2(k) = (d(k) - a(k)*dp2(k-1))/m
    end do
    x(n) = dp2(n)
    do k = n-1, 1, -1
       x(k) = dp2(k) - cp2(k)*x(k+1)                                  ! sustitucion atras
    end do
  end subroutine thomas

  subroutine dif_y_impl(i, f, G, f_pared, neumann_arriba, f_arriba, dtl)  ! [Ec.46,47]
    implicit none
    integer,  intent(in)    :: i
    real(dp), intent(inout) :: f(NX,NY)
    real(dp), intent(in)    :: G(NX,NY), f_pared, f_arriba, dtl
    logical,  intent(in)    :: neumann_arriba
    integer  :: j, k, n
    real(dp) :: aa(NYG), bb(NYG), cc(NYG), dd(NYG), xx(NYG)
    real(dp) :: dyp, dym, dyc, Gp, Gm
    n = NYG
    aa(1) = 0.0_dp;  bb(1) = 1.0_dp;  cc(1) = 0.0_dp;  dd(1) = f_pared  ! Dirichlet pared
    do k = 2, n-1
       j   = NYS + k - 1
       dyp = yg(i,j+1)-yg(i,j);  dym = yg(i,j)-yg(i,j-1)              ! pasos y
       dyc = 0.5_dp*(dyp+dym)
       Gp  = 0.5_dp*(G(i,j+1)+G(i,j));  Gm = 0.5_dp*(G(i,j)+G(i,j-1)) ! G cara
       aa(k) = -dtl*Gm/(dym*dyc)/rho(i,j)                             ! [Ec.47] a_k
       cc(k) = -dtl*Gp/(dyp*dyc)/rho(i,j)                             ! [Ec.47] c_k
       bb(k) = 1.0_dp - aa(k) - cc(k)                                 ! [Ec.47] b_k
       dd(k) = f(i,j)                                                 ! [Ec.46] f*
    end do
    if (neumann_arriba) then
       aa(n) = -1.0_dp;  bb(n) = 1.0_dp;  cc(n) = 0.0_dp;  dd(n) = 0.0_dp  ! Neumann arriba
    else
       aa(n) = 0.0_dp;   bb(n) = 1.0_dp;  cc(n) = 0.0_dp;  dd(n) = f_arriba  ! Dirichlet
    end if
    call thomas(aa,bb,cc,dd,n,xx)                                     ! [Ec.46] resolver
    do k = 1, n
       f(i, NYS+k-1) = xx(k)
    end do
  end subroutine dif_y_impl

  subroutine transporte_escalares()                                   ! [Ec.34,35]
    implicit none
    integer  :: i, j, s
    real(dp) :: uu, vr, adv, dif, h_pared
    Ynew = Ysp;  hnew = hmix
    !$omp parallel do default(shared) private(i,j,s,uu,vr,adv,dif)
    do j = NYS+1, NY
       do i = 2, NX
          uu = u(i,j)
          vr = v(i,j) - u(i,j)*dYdx(i,j)                              ! [Ec.39] v_r
          do s = 1, N_ESP
             adv = adv_up(Ysp(:,:,s), i, j, uu, vr)                   ! adveccion Yk
             dif = difx(Ysp(:,:,s), Dm_, i, j)                        ! difusion x Yk
             Ynew(i,j,s) = Ysp(i,j,s) + dt*(-adv + dif/rho(i,j))      ! [Ec.34] transporte Yk
          end do
          adv = adv_up(hmix, i, j, uu, vr)                            ! adveccion h
          dif = difx(hmix, Dm_, i, j)                                 ! difusion x h
          hnew(i,j) = hmix(i,j) + dt*(-adv + dif/rho(i,j) + qvol(i,j)/rho(i,j))  ! [Ec.35]
       end do
    end do
    !$omp end parallel do
    Ysp = Ynew;  hmix = hnew
    !$omp parallel do default(shared) private(i,s,h_pared)
    do i = 2, NX
       do s = 1, N_ESP
          call dif_y_impl_sp(i, s, dt)                                ! [Ec.34] difusion y
       end do
       h_pared = hmix(i,NYS)                                          ! copia local
       call dif_y_impl(i, hmix, Dm_, h_pared, .true., 0.0_dp, dt)     ! [Ec.35] difusion y
    end do
    !$omp end parallel do
    do j = NYS, NY
       do i = 1, NX
          do s = 1, N_ESP
             if (Ysp(i,j,s) < 0.0_dp) Ysp(i,j,s) = 0.0_dp             ! clip Y
             if (Ysp(i,j,s) > 1.0_dp) Ysp(i,j,s) = 1.0_dp
          end do
          Ysp(i,j,:) = Ysp(i,j,:)/max(sum(Ysp(i,j,:)), 1.0e-12_dp)    ! normalizacion
       end do
    end do
  end subroutine transporte_escalares

  subroutine dif_y_impl_sp(i, s, dtl)                                 ! [Ec.46,47] por especie
    implicit none
    integer,  intent(in) :: i, s
    real(dp), intent(in) :: dtl
    integer  :: j, k, n
    real(dp) :: aa(NYG), bb(NYG), cc(NYG), dd(NYG), xx(NYG)
    real(dp) :: dyp, dym, dyc, Gp, Gm
    n = NYG
    aa(1) = 0.0_dp;  bb(1) = 1.0_dp;  cc(1) = 0.0_dp;  dd(1) = Ysp(i,NYS,s)  ! Dirichlet pared
    do k = 2, n-1
       j   = NYS + k - 1
       dyp = yg(i,j+1)-yg(i,j);  dym = yg(i,j)-yg(i,j-1)
       dyc = 0.5_dp*(dyp+dym)
       Gp  = 0.5_dp*(Dm_(i,j+1)+Dm_(i,j));  Gm = 0.5_dp*(Dm_(i,j)+Dm_(i,j-1))  ! D cara
       aa(k) = -dtl*Gm/(dym*dyc)/rho(i,j)                             ! [Ec.47] a_k
       cc(k) = -dtl*Gp/(dyp*dyc)/rho(i,j)                             ! [Ec.47] c_k
       bb(k) = 1.0_dp - aa(k) - cc(k)                                 ! [Ec.47] b_k
       dd(k) = Ysp(i,j,s)                                             ! [Ec.46] f*
    end do
    aa(n) = -1.0_dp;  bb(n) = 1.0_dp;  cc(n) = 0.0_dp;  dd(n) = 0.0_dp  ! Neumann arriba
    call thomas(aa,bb,cc,dd,n,xx)
    do k = 1, n
       Ysp(i,NYS+k-1,s) = xx(k)
    end do
  end subroutine dif_y_impl_sp

  subroutine quimica()                                                ! [Ec.27-29]
    implicit none
    integer  :: i, j, n_sub
    real(dp) :: Yl(N_ESP), Tl, rl, Y_rad
    n_sub_max  = 0
    n_sub_tot  = 0
    n_celdas_q = 0
    !$omp parallel do default(shared) private(i,j,n_sub,Yl,Tl,rl,Y_rad) &
    !$omp   reduction(max:n_sub_max) reduction(+:n_sub_tot,n_celdas_q) &
    !$omp   schedule(dynamic,8)
    do j = NYS, NY
       do i = 1, NX
          Y_rad = Ysp(i,j,I_H) + Ysp(i,j,I_O) + Ysp(i,j,I_OH)         ! suma radicales
          if (Tg(i,j) < T_QUIM_MIN .and. Y_rad < Y_TRAZA) cycle       ! criterio salto (§18)
          Yl = Ysp(i,j,:)
          Tl = Tg(i,j)
          rl = rho(i,j)
          call paso_quimico(dt, rl, Yl, Tl, n_sub)                    ! [Ec.27-29]
          Ysp(i,j,:) = Yl
          Tg(i,j)    = Tl
          n_sub_max  = max(n_sub_max, n_sub)
          n_sub_tot  = n_sub_tot + n_sub
          n_celdas_q = n_celdas_q + 1
       end do
    end do
    !$omp end parallel do
  end subroutine quimica

  subroutine recuperar_T()                                            ! [Ec.7]
    implicit none
    integer  :: i, j
    real(dp) :: Tl
    !$omp parallel do default(shared) private(i,j,Tl)
    do j = NYS, NY
       do i = 1, NX
          Tl = Tg(i,j)
          call temperatura_de_entalpia(hmix(i,j), Ysp(i,j,:), Tl)     ! [Ec.7] h->T
          Tg(i,j) = Tl
       end do
    end do
    !$omp end parallel do
  end subroutine recuperar_T

  subroutine fuente_radiacion()                                       ! [Ec.19]
    implicit none
    integer  :: i, j
    real(dp) :: a_k, Xc, Xw, sumYW, tau, dyl, q, Tl
    qvol = 0.0_dp
    qrad = 0.0_dp
    if (.not. usar_radiacion) then
       call fuente_chispa()
       return
    end if
    do i = 1, NX
       q = 0.0_dp;  tau = 0.0_dp
       do j = NYS, NY-1
          Tl    = Tg(i,j)
          sumYW = sum(Ysp(i,j,:)/W_ESP)
          if (sumYW < 1.0e-12_dp) sumYW = 1.0_dp/W_ESP(I_O2)
          Xc = (Ysp(i,j,I_CO2)/W_ESP(I_CO2))/sumYW                    ! X_CO2
          Xw = (Ysp(i,j,I_H2O)/W_ESP(I_H2O))/sumYW                    ! X_H2O
          a_k = C_abs*(Xc + Xw)                                        ! kappa gris
          dyl = yg(i,j+1)-yg(i,j)                                      ! dy
          qvol(i,j) = qvol(i,j) - q_radiacion(Tl, Ysp(i,j,:), T_amb)  ! [Ec.19] qvol
          q   = q + a_k*SIGMA_SB*Tl**4*exp(-tau)*dyl                  ! flujo piso
          tau = tau + a_k*dyl                                          ! atenuacion
       end do
       qrad(i) = q
    end do
    call fuente_chispa()
  end subroutine fuente_radiacion

  subroutine fuente_chispa()                                          ! chispa de ignicion
    implicit none
    integer :: i, j
    if (tiempo < t_chispa_ini .or. tiempo > t_chispa_fin) return
    select case (modo_chispa)
    case (2)
       do i = 1, NX                                                   ! caja fija
          if (xg(i) < x_chispa_0 .or. xg(i) > x_chispa_1) cycle
          do j = NYS+1, NY
             if (yg(i,j) < y_chispa_0) cycle                          ! bajo piso
             if (yg(i,j) > y_chispa_1) exit                           ! sobre techo
             qvol(i,j) = qvol(i,j) + P_chispa                         ! fuente
          end do
       end do
    case default
       do i = 1, NX                                                   ! pegada a pared
          if (xg(i) < x_chispa_0 .or. xg(i) > x_chispa_1) cycle
          do j = NYS+1, NY
             if (yg(i,j) - esp(i) > dy_chispa) exit
             qvol(i,j) = qvol(i,j) + P_chispa
          end do
       end do
    end select
  end subroutine fuente_chispa

  subroutine balance_superficie()                                     ! [Ec.40-42]
    implicit none
    integer  :: i, k, s
    real(dp) :: Ta, Tb, Tm, Fa, Fb, Fm, dg, ds, h2, Yl(N_ESP)
    real(dp) :: rD, sumY, cpm, hw
    !$omp parallel do default(shared) &
    !$omp   private(i,k,s,Ta,Tb,Tm,Fa,Fb,Fm,dg,ds,h2,Yl,rD,sumY,cpm,hw)
    do i = 1, NX
       if (agotado(i)) then
          Tw(i) = Tg(i,NYS);  rdot(i) = 0.0_dp;  mdot(i) = 0.0_dp;  cycle  ! agotada
       end if
       dg = yg(i,NYS+1)-yg(i,NYS)                                     ! delta_g [Ec.41]
       ds = yg(i,NYS)-yg(i,NYS-1)                                     ! delta_s [Ec.41]
       h2 = hmix(i,NYS+1)                                             ! h celda 2 [Ec.41]
       Yl = Ysp(i,NYS,:)
       Ta = 280.0_dp;  Tb = 1400.0_dp                                 ! bracket
       Fa = f_sup(i, Ta, dg, ds, h2, Yl)                              ! [Ec.41] F(Ta)
       Fb = f_sup(i, Tb, dg, ds, h2, Yl)                              ! [Ec.41] F(Tb)
       if (Fa < 0.0_dp) then
          Tw(i) = Ta                                                  ! sin balance
       else
          do k = 1, 6
             if (Fb <= 0.0_dp) exit                                   ! bracket ok
             Tb = Tb + 400.0_dp                                       ! ensanchar
             Fb = f_sup(i, Tb, dg, ds, h2, Yl)
          end do
          do k = 1, 45                                                ! biseccion
             Tm = 0.5_dp*(Ta+Tb)
             Fm = f_sup(i, Tm, dg, ds, h2, Yl)
             if (Fa*Fm <= 0.0_dp) then
                Tb = Tm
             else
                Ta = Tm;  Fa = Fm
             end if
             if (Tb-Ta < 1.0e-4_dp) exit
          end do
          Tw(i) = 0.5_dp*(Ta+Tb)                                      ! Tw solucion
       end if
       if (es_grano(i)) then
          rdot(i) = A_pyr*exp(-Ea_pyr/(2.0_dp*RU_MOLAR*Tw(i)))        ! [Ec.40] r
       else
          rdot(i) = 0.0_dp                                            ! liner: r=0
       end if
       mdot(i)  = rho_s*rdot(i)                                       ! [Ec.40] m''
       qconv(i) = lam_(i,NYS)/cp_(i,NYS)*(h2 - hmix(i,NYS))/dg        ! [Ec.41] q_conv

       rD   = Dm_(i,NYS)/dg                                           ! [Ec.42] r_D
       sumY = 0.0_dp
       do s = 1, N_ESP
          if (s == I_C2H4) then
             Ysp(i,NYS,s) = (rD*Ysp(i,NYS+1,s) + mdot(i))/(mdot(i) + rD)  ! [Ec.42] con fuente
          else
             Ysp(i,NYS,s) =  rD*Ysp(i,NYS+1,s)            /(mdot(i) + rD)  ! [Ec.42] sin fuente
          end if
          if (Ysp(i,NYS,s) < 0.0_dp) Ysp(i,NYS,s) = 0.0_dp            ! clip Y
          sumY = sumY + Ysp(i,NYS,s)
       end do
       Ysp(i,NYS,:) = Ysp(i,NYS,:)/max(sumY,1.0e-12_dp)               ! normalizacion
       Tg(i,NYS) = Tw(i)
       call mezcla_cp_h(Ysp(i,NYS,:), Tw(i), cpm, hw)                 ! [Ec.6]
       hmix(i,NYS) = hw
    end do
    !$omp end parallel do
  end subroutine balance_superficie

  function f_sup(i, Twl, dg, ds, h2, Yl) result(F)                    ! [Ec.41] F(Tw)=0
    implicit none
    integer,  intent(in) :: i
    real(dp), intent(in) :: Twl, dg, ds, h2, Yl(N_ESP)
    real(dp) :: F, cpm, hw, r, qg, qr
    call mezcla_cp_h(Yl, Twl, cpm, hw)                                ! [Ec.6]
    qg = lam_(i,NYS)/cp_(i,NYS)*(h2 - hw)/dg                          ! [Ec.41] q_conv
    qr = 0.0_dp
    if (usar_radiacion) qr = emis_w*(qrad(i) - SIGMA_SB*Twl**4)       ! [Ec.41] q_rad
    if (tiempo <= t_ign_fin .and. xg(i) >= x_ign_0 .and. xg(i) <= x_ign_1) &
         qr = qr + q_ign                                              ! encendedor [Ec.41]
    if (es_grano(i)) then
       r = A_pyr*exp(-Ea_pyr/(2.0_dp*RU_MOLAR*Twl))                   ! [Ec.40] r
    else
       r = 0.0_dp
    end if
    F  = qg + qr - rho_s*r*dh_pyr - lam_s*(Twl - Ts(i,NYS-1))/ds      ! [Ec.41] balance
  end function f_sup

  subroutine solido_1d(dtl)                                           ! [Ec.43,44,46,47]
    implicit none
    real(dp), intent(in) :: dtl
    integer  :: i, j, n
    real(dp) :: aa(NYS), bb(NYS), cc(NYS), dd(NYS), xx(NYS)
    real(dp) :: dyp, dym, dyc, cpa, wj, conv, ysl(NYS)
    n = NYS
    !$omp parallel do default(shared) &
    !$omp   private(i,j,aa,bb,cc,dd,xx,dyp,dym,dyc,cpa,wj,conv,ysl)
    do i = 1, NX
       if (agotado(i)) cycle
       do j = 1, NYS
          ysl(j) = esp(i)*fsol(j)                                     ! y columna
       end do
       aa(1) = 0.0_dp;  bb(1) = 1.0_dp;  cc(1) = -1.0_dp;  dd(1) = 0.0_dp  ! fondo adiabatico
       do j = 2, n-1
          dyp = ysl(j+1)-ysl(j);  dym = ysl(j)-ysl(j-1)               ! pasos y
          dyc = 0.5_dp*(dyp+dym)
          cpa = cp_ap(Ts(i,j))                                        ! [Ec.44] cp ap
          aa(j) = -dtl*lam_s/(dym*dyc)/(rho_s*cpa)                    ! [Ec.47] a_k
          cc(j) = -dtl*lam_s/(dyp*dyc)/(rho_s*cpa)                    ! [Ec.47] c_k
          bb(j) = 1.0_dp - aa(j) - cc(j)                              ! [Ec.47] b_k
          wj   = -fsol(j)*rdot(i)                                     ! [Ec.43] w_j
          conv = wj*(Ts(i,j+1)-Ts(i,j-1))/(dyp+dym)                   ! [Ec.43] adv malla
          dd(j) = Ts(i,j) + dtl*conv                                  ! [Ec.43] fuente
       end do
       aa(n) = 0.0_dp;  bb(n) = 1.0_dp;  cc(n) = 0.0_dp;  dd(n) = Tw(i)  ! Dirichlet Tw
       call thomas(aa,bb,cc,dd,n,xx)                                  ! [Ec.46] resolver
       do j = 1, n
          Ts(i,j) = min(max(xx(j), 200.0_dp), 1500.0_dp)              ! clip T
       end do
    end do
    !$omp end parallel do
  end subroutine solido_1d

  pure function cp_ap(T) result(c)                                    ! [Ec.44]
    implicit none
    real(dp), intent(in) :: T
    real(dp) :: c
    c = cp_s
    if (abs(T - T_fus) < dT_fus) c = cp_s + L_fus/(2.0_dp*dT_fus)     ! [Ec.44] fusion
  end function cp_ap

  subroutine mover_superficie(dtl)                                    ! [Ec.45,37]
    implicit none
    real(dp), intent(in) :: dtl
    integer :: i
    do i = 1, NX
       if (agotado(i) .or. .not. es_grano(i)) cycle                   ! solo grano
       esp(i) = esp(i) - rdot(i)*dtl                                  ! [Ec.45] e^n+1
       if (esp(i) <= e_min) then
          esp(i) = e_min                                              ! queda liner
          agotado(i) = .true.
          rdot(i) = 0.0_dp;  mdot(i) = 0.0_dp
       end if
    end do
    call construir_malla()                                            ! [Ec.37] nueva y
  end subroutine mover_superficie

  subroutine gradientes_viscosos()                                    ! [Ec.33]
    implicit none
    integer :: i, j
    !$omp parallel do default(shared) private(i,j)
    do j = NYS, NY
       do i = 1, NX
          tau_a(i,j) = mu_(i,j)*d_dx_fis(u,i,j)                       ! [Ec.33] mu du/dx
          tau_b(i,j) = mu_(i,j)*d_dy_g (u,i,j)                        ! [Ec.33] mu du/dy
          tau_c(i,j) = mu_(i,j)*d_dx_fis(v,i,j)                       ! [Ec.33] mu dv/dx
          tau_d(i,j) = mu_(i,j)*d_dy_g (v,i,j)                        ! [Ec.33] mu dv/dy
          mu43(i,j)  = (4.0_dp/3.0_dp)*mu_(i,j)                       ! [Ec.33] 4/3 mu
       end do
    end do
    !$omp end parallel do
  end subroutine gradientes_viscosos

  subroutine filtro_sierra(fx, fy)                                    ! filtro par-impar
    implicit none
    real(dp), intent(inout) :: fx(NX,NY), fy(NX,NY)
    real(dp) :: ax(NX,NY), ay(NX,NY)
    integer  :: i, j
    if (EPS4 <= 0.0_dp) return
    ax = 0.0_dp;  ay = 0.0_dp
    !$omp parallel do default(shared) private(i,j)
    do j = NYS+1, NY-1
       do i = 3, NX-2
          ax(i,j) = fx(i+2,j) - 4.0_dp*fx(i+1,j) + 6.0_dp*fx(i,j) &
                  - 4.0_dp*fx(i-1,j) + fx(i-2,j)                      ! 4a dif x
          ay(i,j) = fy(i+2,j) - 4.0_dp*fy(i+1,j) + 6.0_dp*fy(i,j) &
                  - 4.0_dp*fy(i-1,j) + fy(i-2,j)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do default(shared) private(i,j)
    do j = NYS+2, NY-2
       do i = 2, NX-1
          ax(i,j) = ax(i,j) + fx(i,j+2) - 4.0_dp*fx(i,j+1) + 6.0_dp*fx(i,j) &
                  - 4.0_dp*fx(i,j-1) + fx(i,j-2)                      ! 4a dif y
          ay(i,j) = ay(i,j) + fy(i,j+2) - 4.0_dp*fy(i,j+1) + 6.0_dp*fy(i,j) &
                  - 4.0_dp*fy(i,j-1) + fy(i,j-2)
       end do
    end do
    !$omp end parallel do
    !$omp parallel do default(shared) private(i,j)
    do j = NYS+1, NY-1
       do i = 2, NX-1
          fx(i,j) = fx(i,j) - EPS4*ax(i,j)                            ! aplicar filtro
          fy(i,j) = fy(i,j) - EPS4*ay(i,j)
       end do
    end do
    !$omp end parallel do
  end subroutine filtro_sierra

  subroutine momento()                                                ! [Ec.32,33]
    implicit none
    integer  :: i, j
    real(dp) :: uu, vr, advu, advv, difu, difv, flot, v_pared
    real(dp), parameter :: c43 = 4.0_dp/3.0_dp, c23 = 2.0_dp/3.0_dp
    if (tensor_completo) call gradientes_viscosos()                   ! [Ec.33]
    if (MOMENTO_CONSISTENTE) then
       rus = rho*u;  rvs = rho*v                                      ! [Ec.32] rho^n+1 u^n
    else
       rus = ru_;  rvs = rv_
    end if
    !$omp parallel do default(shared) private(i,j,uu,vr,advu,advv,difu,difv,flot)
    do j = NYS+1, NY-1
       do i = 2, NX-1
          uu = u(i,j)
          vr = v(i,j) - u(i,j)*dYdx(i,j)                              ! [Ec.39] v_r
          advu = rho(i,j)*adv_up(u, i, j, uu, vr)                     ! [Ec.32] adv u
          advv = rho(i,j)*adv_up(v, i, j, uu, vr)                     ! [Ec.32] adv v
          if (tensor_completo) then
             difu = c43*difx(u, mu_, i, j) - c23*d_dx_fis(tau_d, i, j) &
                  + d_dy_g(tau_c, i, j)                               ! [Ec.33] dtau_xx + dtau_xy
             difv = difx(v, mu_, i, j) + d_dx_fis(tau_b, i, j) &
                  - c23*d_dy_g(tau_a, i, j)                           ! [Ec.33] dtau_yx + dtau_yy
          else
             difu = difx(u, mu_, i, j)                                ! reducido
             difv = difx(v, mu_, i, j)
          end if
          flot = (rho(i,j) - rho_in)*g_y                              ! [Ec.32] flotacion
          rus(i,j) = rus(i,j) + dt*(-advu + difu)                     ! [Ec.32] predictor u
          rvs(i,j) = rvs(i,j) + dt*(-advv + difv + flot)              ! [Ec.32] predictor v
       end do
    end do
    !$omp end parallel do
    !$omp parallel do default(shared) private(i,j,v_pared)
    do i = 2, NX-1
       do j = NYS, NY
          u(i,j) = rus(i,j)/rho(i,j)                                  ! [Ec.32] u = rho u/rho
          v(i,j) = rvs(i,j)/rho(i,j)
       end do
       u(i,NYS) = 0.0_dp                                              ! no slip
       v(i,NYS) = mdot(i)/max(rho(i,NYS),1.0e-6_dp)                   ! soplado [Ec.42]
       v_pared  = v(i,NYS)                                            ! copia local
       call dif_y_impl(i, u, mu_, 0.0_dp, .false., 0.0_dp, dt)        ! [Ec.46] dif y u
       if (tensor_completo) then
          call dif_y_impl(i, v, mu43, v_pared, .false., 0.0_dp, dt)   ! [Ec.46] dif y v (4/3)
       else
          call dif_y_impl(i, v, mu_,  v_pared, .false., 0.0_dp, dt)
       end if
       do j = NYS, NY
          rus(i,j) = rho(i,j)*u(i,j)
          rvs(i,j) = rho(i,j)*v(i,j)
       end do
    end do
    !$omp end parallel do
    call filtro_sierra(rus, rvs)                                      ! par-impar
    call contorno_velocidad(rus, rvs)
  end subroutine momento

  subroutine contorno_velocidad(mu_x, mv_y)                           ! contornos u,v
    implicit none
    real(dp), intent(inout) :: mu_x(NX,NY), mv_y(NX,NY)
    integer :: i, j
    do j = NYS, NY
       if (entrada_abierta(j)) then
          mu_x(1,j) = rho(1,j)*U_ent*cos_ent                          ! entrada u
          mv_y(1,j) = rho(1,j)*U_ent*sin_ent                          ! entrada v
       else
          mu_x(1,j) = 0.0_dp                                          ! pared ciega
          mv_y(1,j) = 0.0_dp
       end if
       mu_x(NX,j) = mu_x(NX-1,j)                                      ! salida d/dn=0
       mv_y(NX,j) = mv_y(NX-1,j)
       if (SIN_REFLUJO) mu_x(NX,j) = max(mu_x(NX,j), 0.0_dp)          ! anti-reflujo
    end do
    mu_x(1,NYS) = 0.0_dp;  mv_y(1,NYS) = 0.0_dp
    do i = 1, NX
       mu_x(i,NY)  = 0.0_dp                                           ! techo no slip
       mv_y(i,NY)  = 0.0_dp
       mu_x(i,NYS) = 0.0_dp                                           ! pared no slip
       mv_y(i,NYS) = mdot(i)                                          ! pared soplado
    end do
  end subroutine contorno_velocidad

  subroutine contorno_escalares()                                     ! contornos Y,h,T
    implicit none
    integer  :: i, j
    real(dp) :: cpm, hm
    do j = NYS, NY
       if (entrada_abierta(j)) then
          Ysp(1,j,:) = 0.0_dp;  Ysp(1,j,I_O2) = 1.0_dp                ! entrada O2
          Tg(1,j) = T_in
          call mezcla_cp_h(Ysp(1,j,:), T_in, cpm, hm)                 ! [Ec.6]
          hmix(1,j) = hm
       else
          Ysp(1,j,:) = Ysp(2,j,:)                                     ! pared: d/dn=0
          hmix(1,j)  = hmix(2,j)
          Tg(1,j)    = Tg(2,j)
       end if
       Ysp(NX,j,:) = Ysp(NX-1,j,:)                                    ! salida d/dn=0
       hmix(NX,j)  = hmix(NX-1,j)
       Tg(NX,j)    = Tg(NX-1,j)
    end do
    do i = 1, NX
       Ysp(i,NY,:) = Ysp(i,NY-1,:)                                    ! techo d/dn=0
       hmix(i,NY)  = hmix(i,NY-1)
       Tg(i,NY)    = Tg(i,NY-1)
    end do
  end subroutine contorno_escalares

  subroutine poisson_correccion()                                     ! [Ec.48,49,46]
    implicit none
    integer  :: i, j, k, n
    real(dp) :: rhs(NX,NY), res, resmax
    real(dp) :: aa(NYG), bb(NYG), cc(NYG), dd(NYG), xx(NYG)
    real(dp) :: dxp, dxm, dxc, dyp, dym, dyc, ae, aw, an, as_, ap
    real(dp) :: dpx, dpy, pmax, Mi, fac, cruz
    do j = NYS, NY
       do i = 1, NX
          div_(i,j) = d_dx_fis(rus,i,j) + d_dy_g(rvs,i,j)             ! [Ec.48] div rho u*
       end do
    end do
    do j = NYS, NY
       do i = 1, NX
          rhs(i,j) = (div_(i,j) + (rho(i,j)-rho_o(i,j))/dt)/dt        ! [Ec.48] RHS
       end do
    end do
    diag_cont_ant = maxval(abs(rhs(2:NX-1, NYS+1:NY-1)))*dt
    pmax = max(1.0_dp, maxval(abs(pr(:,NYS:NY))))                      ! escala residuo
    n = NYG
    do k = 1, it_poisson
       resmax = 0.0_dp
       do i = 2, NX-1                                                 ! barrido linea y
          dxp = xg(i+1)-xg(i);  dxm = xg(i)-xg(i-1);  dxc = 0.5_dp*(dxp+dxm)
          ae = 1.0_dp/(dxp*dxc);  aw = 1.0_dp/(dxm*dxc)               ! coef x [Ec.49]
          aa(1) = 0.0_dp;  bb(1) = 1.0_dp;  cc(1) = -1.0_dp;  dd(1) = 0.0_dp  ! Neumann pared
          do j = NYS+1, NY-1
             dyp = yg(i,j+1)-yg(i,j);  dym = yg(i,j)-yg(i,j-1)
             dyc = 0.5_dp*(dyp+dym)
             Mi  = dYdx(i,j)                                          ! [Ec.49] M
             fac = 1.0_dp + Mi*Mi                                     ! [Ec.49] 1+M^2
             an = fac/(dyp*dyc);  as_ = fac/(dym*dyc)                 ! [Ec.49] a_n, a_s
             ap = ae + aw + an + as_                                  ! [Ec.49] a_p
             cruz = -2.0_dp*Mi*(d_dy_g(pr,i+1,j) - d_dy_g(pr,i-1,j)) &
                    /(xg(i+1) - xg(i-1))                              ! [Ec.49] cruzado
             aa(j-NYS+1) = -as_                                       ! [Ec.47] a_k
             cc(j-NYS+1) = -an                                        ! [Ec.47] c_k
             bb(j-NYS+1) =  ap                                        ! [Ec.47] b_k
             dd(j-NYS+1) = -rhs(i,j) + ae*pr(i+1,j) + aw*pr(i-1,j) + cruz  ! [Ec.48] RHS
          end do
          aa(n) = -1.0_dp;  bb(n) = 1.0_dp;  cc(n) = 0.0_dp;  dd(n) = 0.0_dp  ! Neumann techo
          call thomas(aa,bb,cc,dd,n,xx)                               ! [Ec.46] resolver
          do j = NYS, NY
             res = abs(xx(j-NYS+1) - pr(i,j))
             if (res > resmax) resmax = res
             pr(i,j) = xx(j-NYS+1)
          end do
       end do
       do j = NYS, NY
          pr(1,j)  = pr(2,j)                                          ! entrada dp/dn=0
          pr(NX,j) = 0.0_dp                                           ! salida p fija
       end do
       if (resmax < tol_poisson*pmax) exit
    end do
    diag_k_poisson  = min(k, it_poisson)
    diag_poisson_ok = (k <= it_poisson)
    diag_res_p      = resmax/pmax
    do j = NYS+1, NY-1
       do i = 2, NX-1
          dpx = d_dx_fis(pr,i,j)
          dpy = d_dy_g(pr,i,j)
          ru_(i,j) = rus(i,j) - dt*dpx                                ! [Ec.48] correccion u
          rv_(i,j) = rvs(i,j) - dt*dpy                                ! [Ec.48] correccion v
       end do
    end do
    call contorno_velocidad(ru_, rv_)
    do j = NYS, NY
       do i = 1, NX
          u(i,j) = ru_(i,j)/max(rho(i,j),1.0e-6_dp)                   ! u = rho u/rho
          v(i,j) = rv_(i,j)/max(rho(i,j),1.0e-6_dp)
       end do
    end do
    do i = 1, NX
       u(i,NYS) = 0.0_dp;  v(i,NYS) = mdot(i)/max(rho(i,NYS),1.0e-6_dp)
       u(i,NY)  = 0.0_dp;  v(i,NY)  = 0.0_dp
    end do
    diag_cont_des = 0.0_dp;  diag_i_cont = 0;  diag_j_cont = 0
    do j = NYS+1, NY-1                                                ! diagnostico post
       do i = 2, NX-1
          res = abs(d_dx_fis(ru_,i,j) + d_dy_g(rv_,i,j) + (rho(i,j)-rho_o(i,j))/dt)
          if (res > diag_cont_des) then
             diag_cont_des = res;  diag_i_cont = i;  diag_j_cont = j
          end if
       end do
    end do
  end subroutine poisson_correccion

  subroutine diagnostico(itl, tgl)                                    ! log iteracion
    implicit none
    integer,  intent(in) :: itl
    real(dp), intent(in) :: tgl
    integer  :: i, j, im, jm, iconv
    real(dp) :: sp, umx
    umx = -1.0_dp;  im = 1;  jm = NYS
    do j = NYS, NY
       do i = 1, NX
          sp = sqrt(u(i,j)**2 + v(i,j)**2)                            ! |U|
          if (sp /= sp) sp = huge(1.0_dp)                             ! NaN al tope
          if (sp > umx) then
             umx = sp;  im = i;  jm = j
          end if
       end do
    end do
    diag_umax = umx;  diag_i_umax = im;  diag_j_umax = jm
    if (.not. diag_log_abierto) then
       open(newunit=un_diag, file='diag_divergencia.dat', status='replace')
       diag_log_abierto = .true.
       write(un_diag,'(A)') '# it  t_sol[s]  t_gas[s]  dt[s]  Umax[m/s]  i  j  '// &
            'T[K]  rho[kg/m3]  kP  conv  res_p  cont_antes  cont_despues  ic  jc  '// &
            'Tw[K]  mdot[kg/m2s]'
    end if
    if (umx > U_AVISO .and. diag_armado) then
       diag_armado   = .false.
       diag_n_avisos = diag_n_avisos + 1
       write(*,'(A,I3,A,I8,A,ES11.4,A,F9.1,A,I4,A,I4,A)') ' AVISO #', diag_n_avisos, &
            '  it=', itl, '  t_gas=', tgl, ' s  |U|max=', umx, ' m/s en i=', im, &
            ' j=', jm, '  -> volcado frame 9998'
       flush(6)
       call escribir_salidas(9998)                                    ! volcado aviso
    end if
    if (.not. diag_armado .and. umx < 0.5_dp*U_AVISO) diag_armado = .true.
    if (mod(itl, N_DIAG) == 0 .or. umx > U_AVISO) then
       iconv = 0;  if (diag_poisson_ok) iconv = 1
       write(un_diag,'(I9,3(1X,ES12.5),1X,ES12.5,2(1X,I4),2(1X,ES12.5),1X,I4,1X,I2,'// &
            '3(1X,ES11.4),2(1X,I4),2(1X,ES12.5))') &
            itl, tiempo, tgl, dt, umx, im, jm, Tg(im,jm), rho(im,jm), &
            diag_k_poisson, iconv, diag_res_p, diag_cont_ant, diag_cont_des, &
            diag_i_cont, diag_j_cont, Tw(im), mdot(im)
       flush(un_diag)
    end if
  end subroutine diagnostico

  subroutine verificar(itl)                                           ! chequeo NaN/absurdo
    implicit none
    integer, intent(in) :: itl
    integer :: i, j
    do j = NYS, NY
       do i = 1, NX
          if (u(i,j) /= u(i,j) .or. abs(u(i,j)) > 1.0e4_dp .or. &     ! NaN o |u|>1e4
              Tg(i,j) /= Tg(i,j) .or. Tg(i,j) > 6000.0_dp) then       ! NaN o T>6000
             write(*,'(A,I8,A,I4,A,I4,A,ES12.4,A,ES12.4)') &
               ' DIVERGE en it=',itl,'  i=',i,' j=',j,'  u=',u(i,j),' T=',Tg(i,j)
             write(*,'(A,ES12.4,A,I4,A,I4,A)') '   |U| maximo=', diag_umax, &
               '  en i=', diag_i_umax, ' j=', diag_j_umax, &
               '   (historial en diag_divergencia.dat)'
             if (diag_log_abierto) flush(un_diag)
             call escribir_salidas(9999)
             stop 2
          end if
       end do
    end do
  end subroutine verificar

  pure function zz(vv) result(r)                                      ! limpia subnormales
    implicit none
    real(dp), intent(in) :: vv
    real(dp) :: r
    r = vv
    if (abs(r) < 1.0e-30_dp) r = 0.0_dp
  end function zz

  subroutine escribir_malla_fija()                                    ! malla inicial
    implicit none
    integer :: i, j, un, n_sw, n_se, n_ne, n_nw
    open(newunit=un, file='malla_nodos.dat', status='replace')
    do j = 1, NY
       do i = 1, NX
          write(un,'(2(ES14.6,1X))') zz(xg(i)), zz(yg(i,j))
       end do
    end do
    close(un)
    open(newunit=un, file='malla_triangulos.dat', status='replace')
    do j = 1, NY-1
       do i = 1, NX-1
          n_sw = (j-1)*NX + i
          n_se = n_sw + 1
          n_ne = j*NX + i + 1
          n_nw = j*NX + i
          write(un,'(3(I8,1X))') n_sw-1, n_se-1, n_ne-1               ! triangulo 1
          write(un,'(3(I8,1X))') n_sw-1, n_ne-1, n_nw-1               ! triangulo 2
       end do
    end do
    close(un)
  end subroutine escribir_malla_fija

  subroutine escribir_salidas(frm)                                    ! volcado de campos
    implicit none
    integer, intent(in) :: frm
    integer :: i, j, un, s
    real(dp) :: sp
    character(len=40) :: nom
    write(nom,'(A,I4.4,A)') 'campo_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    do j = 1, NY
       do i = 1, NX
          if (j < NYS) then
             write(un,'(4(ES14.6,1X))') 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp  ! solido: 0
          else
             sp = sqrt(u(i,j)**2 + v(i,j)**2)
             write(un,'(4(ES14.6,1X))') zz(u(i,j)), zz(v(i,j)), zz(pr(i,j)), zz(sp)
          end if
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'temperatura_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    do j = 1, NY
       do i = 1, NX
          if (j < NYS) then
             write(un,'(3(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), zz(Ts(i,j))  ! T solido
          else
             write(un,'(3(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), zz(Tg(i,j))  ! T gas
          end if
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'especies_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    write(un,'(A)') '# x  y  C2H4  O2  CO2  H2O  CO  H2  H  O  OH'
    do j = 1, NY
       do i = 1, NX
          if (j < NYS) then
             write(un,'(11(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), &
                  (0.0_dp, s=1,N_ESP)                                 ! solido: 0
          else
             write(un,'(11(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), &
                  (zz(Ysp(i,j,s)), s=1,N_ESP)                         ! Y gas
          end if
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'obstaculo_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    do j = 1, NY
       do i = 1, NX
          if (j <= NYS .and. .not. agotado(i)) then
             write(un,'(2(ES14.6,1X),I3)') zz(xg(i)), zz(yg(i,j)), 1  ! obstaculo
          else
             write(un,'(2(ES14.6,1X),I3)') zz(xg(i)), zz(yg(i,j)), 0
          end if
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'malla_nodos_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    do j = 1, NY
       do i = 1, NX
          write(un,'(2(ES14.6,1X))') zz(xg(i)), zz(yg(i,j))           ! malla actual
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'densidad_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    do j = 1, NY
       do i = 1, NX
          if (j < NYS) then
             write(un,'(3(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), zz(rho_s)
          else
             write(un,'(3(ES14.6,1X))') zz(xg(i)), zz(yg(i,j)), zz(rho(i,j))  ! [Ec.9]
          end if
       end do
    end do
    close(un)
    write(nom,'(A,I4.4,A)') 'superficie_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    write(un,'(A)') '# x  espesor  Tw  rdot_mm_s  q_conv  q_rad'
    do i = 1, NX
       write(un,'(6(ES14.6,1X))') zz(xg(i)), zz(esp(i)), zz(Tw(i)), &
                                  zz(rdot(i)*1.0e3_dp), zz(qconv(i)), zz(qrad(i))
    end do
    close(un)
    call escribir_capa_limite(frm)
    write(nom,'(A,I4.4,A)') 'tiempo_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    write(un,'(ES16.8)') tiempo
    close(un)
  end subroutine escribir_salidas

  subroutine escribir_capa_limite(frm)                                ! delta99, y_llama, Re_x
    implicit none
    integer, intent(in) :: frm
    integer  :: i, j, un, j_fl, j_d
    real(dp) :: umax, d99, y_fl, ohmax, Re_x, xrel, frac
    real(dp) :: rho_e, u_e, mu_e, Tcol
    character(len=40) :: nom
    write(nom,'(A,I4.4,A)') 'capa_limite_', frm, '.dat'
    open(newunit=un, file=nom, status='replace')
    write(un,'(A)') '# x[m]  delta99[m]  y_llama[m]  frac  Tmax[K]  '// &
                    'u_e[m/s]  Re_x  OH_max'
    do i = 1, NX
       umax = 0.0_dp
       do j = NYS, NY
          if (u(i,j) > umax) umax = u(i,j)                            ! u_max columna
       end do
       d99 = yg(i,NY) - yg(i,NYS)                                     ! default puerto
       j_d = NY
       if (umax > 1.0e-6_dp) then
          do j = NYS, NY
             if (u(i,j) >= 0.99_dp*umax) then                         ! delta99
                d99 = yg(i,j) - yg(i,NYS)
                j_d = j
                exit
             end if
          end do
       end if
       if (d99 < 1.0e-9_dp) d99 = 1.0e-9_dp
       ohmax = 0.0_dp;  j_fl = NYS                                   ! y_llama: max OH
       do j = NYS, NY
          if (Ysp(i,j,I_OH) > ohmax) then
             ohmax = Ysp(i,j,I_OH);  j_fl = j
          end if
       end do
       y_fl = yg(i,j_fl) - yg(i,NYS)
       frac = y_fl / d99                                              ! ~0.33 Jens
       Tcol = 0.0_dp
       do j = NYS, NY
          if (Tg(i,j) > Tcol) Tcol = Tg(i,j)                          ! T_max columna
       end do
       rho_e = rho(i,j_d);  u_e = u(i,j_d);  mu_e = mu_(i,j_d)       ! Re local
       xrel  = max(xg(i) - x_ini_gr, 0.0_dp)
       Re_x  = 0.0_dp
       if (mu_e > 1.0e-12_dp) Re_x = rho_e*abs(u_e)*xrel/mu_e
       write(un,'(8(ES14.6,1X))') zz(xg(i)), zz(d99), zz(y_fl), zz(frac), &
                                  zz(Tcol), zz(u_e), zz(Re_x), zz(ohmax)
    end do
    close(un)
  end subroutine escribir_capa_limite

end module mod_solver