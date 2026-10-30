module mod_cinetica
  use mod_precision
  use mod_termo
  implicit none

  integer, parameter :: N_REAC = 15                                   ! 1 global + 14 elementales
  integer, parameter :: TIPO_ELEM    = 0                              ! elemental
  integer, parameter :: TIPO_TERCERO = 1                              ! tercer cuerpo M
  integer, parameter :: TIPO_LINDEM  = 2                              ! [Ec.24] Lindemann
  integer, parameter :: TIPO_GLOBAL  = 3                              ! global irreversible

  integer,  save :: tipo_reac(N_REAC)                                 ! tipo reaccion
  logical,  save :: es_reversible(N_REAC)                             ! con inversa
  real(dp), save :: nu_p (N_ESP, N_REAC)                              ! nu' reactivos [Ec.25]
  real(dp), save :: nu_pp(N_ESP, N_REAC)                              ! nu'' productos [Ec.25]
  real(dp), save :: nu_net(N_ESP, N_REAC)                             ! nu''-nu' [Ec.26]
  real(dp), save :: orden_f(N_ESP, N_REAC)                            ! orden directo [Ec.25]
  real(dp), save :: delta_n(N_REAC)                                   ! sum nu_net [Ec.23]
  real(dp), save :: A_f(N_REAC), b_f(N_REAC), Ea_f(N_REAC)            ! Arrhenius dir [Ec.20]
  real(dp), save :: A_0(N_REAC), b_0(N_REAC), Ea_0(N_REAC)            ! baja p [Ec.24]
  real(dp), save :: eff_M(N_ESP, N_REAC)                              ! eficiencias M [Ec.24,25]

  integer,  parameter :: MAX_SUBPASOS = 4000                         ! tope subpasos [Ec.28]
  real(dp), parameter :: DY_MAX_SUB   = 0.01_dp                       ! max|dY| subpaso [Ec.28]
  real(dp), parameter :: DT_MAX_SUB   = 25.0_dp                       ! max|dT| subpaso [K]
  real(dp), parameter :: FRAC_MIN_SUB = 1.0e-6_dp                     ! dt_sub minimo [Ec.28]
  real(dp), parameter :: T_QUIM_MIN   = 600.0_dp                      ! T piso (§18)
  real(dp), parameter :: Y_TRAZA      = 1.0e-14_dp                    ! traza (§18)
  real(dp), parameter :: Y_PISO       = 1.0e-20_dp                    ! piso Y [Ec.26]
  real(dp), parameter :: C_FRAC       = 1.0e-12_dp                    ! C piso [kmol/m3] [Ec.25]

  real(dp), save :: kf_tab(N_REAC, N_TABLA)                           ! kf tabulada [Ec.20]
  real(dp), save :: kb_tab(N_REAC, N_TABLA)                           ! kb tabulada [Ec.23]
  real(dp), save :: k0_tab(N_REAC, N_TABLA)                           ! k0 tabulada [Ec.24]
  real(dp), parameter :: DY_JAC = 1.0e-16_dp                          ! pert jacobiano [Ec.28]

contains

  subroutine init_cinetica()                                          ! carga R1-R15 (§16)
    implicit none
    integer :: j, k

    tipo_reac(:)     = TIPO_ELEM                                    ! default
    es_reversible(:) = .true.                                       ! default
    nu_p(:,:)        = 0.0_dp                                       ! default
    nu_pp(:,:)       = 0.0_dp
    orden_f(:,:)     = 0.0_dp
    eff_M(:,:)       = 1.0_dp
    A_0(:)  = 0.0_dp;  b_0(:) = 0.0_dp;  Ea_0(:) = 0.0_dp

    j = 1                                                            ! R1: C2H4+O2->2CO+2H2 (W&D)
    tipo_reac(j)     = TIPO_GLOBAL
    es_reversible(j) = .false.                                       ! global no se invierte
    nu_p (I_C2H4, j) = 1.0_dp
    nu_p (I_O2  , j) = 1.0_dp
    nu_pp(I_CO  , j) = 2.0_dp
    nu_pp(I_H2  , j) = 2.0_dp
    orden_f(I_C2H4,j)= 0.1_dp                                        ! orden frac (W&D) [Ec.25]
    orden_f(I_O2  ,j)= 1.65_dp                                       ! orden frac (W&D) [Ec.25]
    A_f(j)  = 1.1246826503806981e10_dp                               ! [Ec.20]
    b_f(j)  = 0.0_dp                                                 ! [Ec.20]
    Ea_f(j) = 1.2552e8_dp                                            ! [Ec.20]

    j = 2                                                            ! R2: 2O+M<=>O2+M
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_O ,j) = 2.0_dp
    nu_pp(I_O2,j) = 1.0_dp
    A_f(j) = 1.2e11_dp;  b_f(j) = -1.0_dp;  Ea_f(j) = 0.0_dp        ! [Ec.20]
    eff_M(I_CO ,j) = 1.75_dp                                         ! [Ec.25]
    eff_M(I_CO2,j) = 3.6_dp
    eff_M(I_H2 ,j) = 2.4_dp
    eff_M(I_H2O,j) = 15.4_dp

    j = 3                                                            ! R3: H+O+M<=>OH+M
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H ,j) = 1.0_dp
    nu_p (I_O ,j) = 1.0_dp
    nu_pp(I_OH,j) = 1.0_dp
    A_f(j) = 5.0e11_dp;  b_f(j) = -1.0_dp;  Ea_f(j) = 0.0_dp        ! [Ec.20]
    eff_M(I_CO ,j) = 1.5_dp                                          ! [Ec.25]
    eff_M(I_CO2,j) = 2.0_dp
    eff_M(I_H2 ,j) = 2.0_dp
    eff_M(I_H2O,j) = 6.0_dp

    j = 4                                                            ! R4: H2+O<=>H+OH
    nu_p (I_H2,j) = 1.0_dp
    nu_p (I_O ,j) = 1.0_dp
    nu_pp(I_H ,j) = 1.0_dp
    nu_pp(I_OH,j) = 1.0_dp
    A_f(j) = 3.87e1_dp;  b_f(j) = 2.7_dp;  Ea_f(j) = 2.619184e7_dp   ! [Ec.20]

    j = 5                                                            ! R5: CO+O(+M)<=>CO2(+M) Lindemann
    tipo_reac(j) = TIPO_LINDEM                                       ! [Ec.24]
    nu_p (I_CO ,j) = 1.0_dp
    nu_p (I_O  ,j) = 1.0_dp
    nu_pp(I_CO2,j) = 1.0_dp
    A_f(j) = 1.8e7_dp;   b_f(j) = 0.0_dp;  Ea_f(j) = 9.97884e6_dp    ! k_inf [Ec.20,24]
    A_0(j) = 6.02e8_dp;  b_0(j) = 0.0_dp;  Ea_0(j) = 1.2552e7_dp     ! k0   [Ec.24]
    eff_M(I_CO ,j) = 1.5_dp                                          ! [Ec.25]
    eff_M(I_CO2,j) = 3.5_dp
    eff_M(I_H2 ,j) = 2.0_dp
    eff_M(I_H2O,j) = 6.0_dp
    eff_M(I_O2 ,j) = 6.0_dp

    j = 6                                                            ! R6: CO+O2<=>CO2+O
    nu_p (I_CO ,j) = 1.0_dp
    nu_p (I_O2 ,j) = 1.0_dp
    nu_pp(I_CO2,j) = 1.0_dp
    nu_pp(I_O  ,j) = 1.0_dp
    A_f(j) = 2.5e9_dp;  b_f(j) = 0.0_dp;  Ea_f(j) = 1.999952e8_dp    ! [Ec.20]

    j = 7                                                            ! R7: H+O2<=>O+OH (ramificacion)
    nu_p (I_H ,j) = 1.0_dp
    nu_p (I_O2,j) = 1.0_dp
    nu_pp(I_O ,j) = 1.0_dp
    nu_pp(I_OH,j) = 1.0_dp
    A_f(j) = 2.65e13_dp;  b_f(j) = -0.6707_dp;  Ea_f(j) = 7.1299544e7_dp  ! [Ec.20]

    j = 8                                                            ! R8: 2H+M<=>H2+M
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H ,j) = 2.0_dp
    nu_pp(I_H2,j) = 1.0_dp
    A_f(j) = 1.0e12_dp;  b_f(j) = -1.0_dp;  Ea_f(j) = 0.0_dp         ! [Ec.20]
    eff_M(I_CO2,j) = 0.0_dp                                          ! tratado en R9-R11
    eff_M(I_H2 ,j) = 0.0_dp
    eff_M(I_H2O,j) = 0.0_dp

    j = 9                                                            ! R9: 2H+H2<=>2H2
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H ,j) = 2.0_dp
    nu_pp(I_H2,j) = 1.0_dp
    A_f(j) = 9.0e10_dp;  b_f(j) = -0.6_dp;  Ea_f(j) = 0.0_dp         ! [Ec.20]
    eff_M(:,   j) = 0.0_dp                                           ! solo H2
    eff_M(I_H2,j) = 1.0_dp

    j = 10                                                           ! R10: 2H+H2O<=>H2+H2O
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H ,j) = 2.0_dp
    nu_pp(I_H2,j) = 1.0_dp
    A_f(j) = 6.0e13_dp;  b_f(j) = -1.25_dp;  Ea_f(j) = 0.0_dp        ! [Ec.20]
    eff_M(:,    j) = 0.0_dp                                          ! solo H2O
    eff_M(I_H2O,j) = 1.0_dp

    j = 11                                                           ! R11: 2H+CO2<=>H2+CO2
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H ,j) = 2.0_dp
    nu_pp(I_H2,j) = 1.0_dp
    A_f(j) = 5.5e14_dp;  b_f(j) = -2.0_dp;  Ea_f(j) = 0.0_dp         ! [Ec.20]
    eff_M(:,    j) = 0.0_dp                                          ! solo CO2
    eff_M(I_CO2,j) = 1.0_dp

    j = 12                                                           ! R12: H+OH+M<=>H2O+M
    tipo_reac(j) = TIPO_TERCERO
    nu_p (I_H  ,j) = 1.0_dp
    nu_p (I_OH ,j) = 1.0_dp
    nu_pp(I_H2O,j) = 1.0_dp
    A_f(j) = 2.2e16_dp;  b_f(j) = -2.0_dp;  Ea_f(j) = 0.0_dp         ! [Ec.20]
    eff_M(I_H2 ,j) = 0.73_dp                                         ! [Ec.25]
    eff_M(I_H2O,j) = 3.65_dp

    j = 13                                                           ! R13: H2+OH<=>H+H2O
    nu_p (I_H2 ,j) = 1.0_dp
    nu_p (I_OH ,j) = 1.0_dp
    nu_pp(I_H  ,j) = 1.0_dp
    nu_pp(I_H2O,j) = 1.0_dp
    A_f(j) = 2.16e5_dp;  b_f(j) = 1.51_dp;  Ea_f(j) = 1.435112e7_dp   ! [Ec.20]

    j = 14                                                           ! R14: 2OH<=>H2O+O (Ea<0 GRI)
    nu_p (I_OH ,j) = 2.0_dp
    nu_pp(I_H2O,j) = 1.0_dp
    nu_pp(I_O  ,j) = 1.0_dp
    A_f(j) = 3.57e1_dp;  b_f(j) = 2.4_dp;  Ea_f(j) = -8.82824e6_dp    ! [Ec.20]

    j = 15                                                           ! R15: CO+OH<=>CO2+H
    nu_p (I_CO ,j) = 1.0_dp
    nu_p (I_OH ,j) = 1.0_dp
    nu_pp(I_CO2,j) = 1.0_dp
    nu_pp(I_H  ,j) = 1.0_dp
    A_f(j) = 4.76e4_dp;  b_f(j) = 1.228_dp;  Ea_f(j) = 2.9288e5_dp    ! [Ec.20]

    do j = 1, N_REAC                                                ! derivados
      do k = 1, N_ESP
        nu_net(k,j) = nu_pp(k,j) - nu_p(k,j)                        ! [Ec.26]
        if (tipo_reac(j) /= TIPO_GLOBAL) orden_f(k,j) = nu_p(k,j)   ! elemental: orden=nu' [Ec.25]
      end do
      delta_n(j) = sum(nu_net(:,j))                                 ! [Ec.23]
    end do
    print *, 'Mecanismo cargado: ', N_REAC, ' reacciones, ', N_ESP, ' especies'
  end subroutine init_cinetica

  subroutine pretabular_k()                                          ! [Ec.20,23,24]
    implicit none
    integer  :: it, j
    real(dp) :: Tl, cp_k(N_ESP), h_k(N_ESP), s_k(N_ESP)
    real(dp) :: kf(N_REAC), kb(N_REAC), expo, lnT, invT

    do it = 1, N_TABLA
      Tl = T_tab(it)
      call props_especies(Tl, cp_k, h_k, s_k)                       ! [Ec.5]
      call constantes_velocidad(Tl, h_k, s_k, kf, kb)               ! [Ec.20,23]
      kf_tab(:,it) = kf(:)
      kb_tab(:,it) = kb(:)
      lnT  = log(Tl)
      invT = 1.0_dp / Tl
      do j = 1, N_REAC
        if (tipo_reac(j) /= TIPO_LINDEM) then
          k0_tab(j,it) = 0.0_dp                                     ! sin k0
          cycle
        end if
        expo = b_0(j)*lnT - Ea_0(j)*invT/RU_KMOL                    ! [Ec.24] k0
        if (expo >  200.0_dp) expo =  200.0_dp
        if (expo < -200.0_dp) expo = -200.0_dp
        k0_tab(j,it) = A_0(j) * exp(expo)                           ! [Ec.24]
      end do
    end do
    print *, 'Constantes de velocidad pretabuladas: ', N_REAC, ' x ', N_TABLA
  end subroutine pretabular_k

  subroutine k_interpoladas(T, kf, kb, k0)                          ! [Ec.20,23,24]
    implicit none
    real(dp), intent(in)  :: T
    real(dp), intent(out) :: kf(N_REAC), kb(N_REAC), k0(N_REAC)
    integer  :: lo, j
    real(dp) :: frac

    call indice_T(T, lo, frac)                                      ! [Ec.4]
    do j = 1, N_REAC
      kf(j) = kf_tab(j,lo) + frac*(kf_tab(j,lo+1) - kf_tab(j,lo))   ! [Ec.5]
      kb(j) = kb_tab(j,lo) + frac*(kb_tab(j,lo+1) - kb_tab(j,lo))   ! [Ec.5]
      k0(j) = k0_tab(j,lo) + frac*(k0_tab(j,lo+1) - k0_tab(j,lo))   ! [Ec.5]
    end do
  end subroutine k_interpoladas

  subroutine constantes_velocidad(T, h_k, s_k, kf, kb)              ! [Ec.20,21,22,23]
    implicit none
    real(dp), intent(in)  :: T, h_k(N_ESP), s_k(N_ESP)
    real(dp), intent(out) :: kf(N_REAC), kb(N_REAC)
    integer  :: j, k
    real(dp) :: dS_Ru, dH_RuT, Kp, Kc, expo, lnT, invT

    lnT  = log(T)
    invT = 1.0_dp / T
    do j = 1, N_REAC
      expo = b_f(j)*lnT - Ea_f(j)*invT/RU_KMOL                     ! [Ec.20] exponent
      if (expo >  200.0_dp) expo =  200.0_dp                       ! anti-overflow
      if (expo < -200.0_dp) expo = -200.0_dp                       ! anti-underflow
      kf(j) = A_f(j) * exp(expo)                                    ! [Ec.20] kf

      if (.not. es_reversible(j)) then
        kb(j) = 0.0_dp                                              ! irreversible
        cycle
      end if

      dS_Ru  = 0.0_dp                                               ! [Ec.22]
      dH_RuT = 0.0_dp
      do k = 1, N_ESP
        if (nu_net(k,j) == 0.0_dp) cycle
        dS_Ru  = dS_Ru  + nu_net(k,j) * s_k(k) * W_ESP(k)           ! [Ec.22] dS°
        dH_RuT = dH_RuT + nu_net(k,j) * h_k(k) * W_ESP(k)           ! [Ec.22] dH°
      end do
      dS_Ru  = dS_Ru  / RU_KMOL
      dH_RuT = dH_RuT / (RU_KMOL * T)

      expo = dS_Ru - dH_RuT                                         ! [Ec.21] exponent
      if (expo >  200.0_dp) expo =  200.0_dp
      if (expo < -200.0_dp) expo = -200.0_dp
      Kp = exp(expo)                                                ! [Ec.21] Kp

      Kc = Kp * (P_REF / (RU_KMOL * T))**delta_n(j)                 ! [Ec.23] Kc
      if (Kc > 1.0e-30_dp) then
        kb(j) = kf(j) / Kc                                          ! [Ec.23] kb=kf/Kc
      else
        kb(j) = 0.0_dp                                              ! equilibrio desplazado
      end if
    end do
  end subroutine constantes_velocidad

  subroutine tasas_produccion(rho, Y, kf, kb, prod, dest)           ! [Ec.25,26]
    implicit none
    real(dp), intent(in)  :: rho, Y(N_ESP), kf(N_REAC), kb(N_REAC)
    real(dp), intent(out) :: prod(N_ESP), dest(N_ESP)

    real(dp) :: C(N_ESP), qf, qb, Cm
    integer  :: j, k

    do k = 1, N_ESP
      C(k) = rho * max(Y(k), 0.0_dp) / W_ESP(k)                     ! [Ec.25] Ck
    end do
    prod(:) = 0.0_dp                                                ! inicializa
    dest(:) = 0.0_dp

    do j = 1, N_REAC
      Cm = 1.0_dp                                                   ! sin tercer cuerpo
      if (tipo_reac(j) == TIPO_TERCERO) then
        Cm = 0.0_dp                                                 ! [Ec.25] M efectivo
        do k = 1, N_ESP
          Cm = Cm + eff_M(k,j) * C(k)                               ! [Ec.25] sum eff*C
        end do
      end if

      qf = kf(j)                                                    ! [Ec.25] qf=kf prod C^nu'
      do k = 1, N_ESP
        if (orden_f(k,j) <= 0.0_dp) cycle
        if (orden_f(k,j) == 1.0_dp) then
          qf = qf * C(k)                                            ! orden 1
        else if (orden_f(k,j) == 2.0_dp) then
          qf = qf * C(k) * C(k)                                     ! orden 2
        else
          if (C(k) <= C_FRAC) then                                  ! evitar (0)^0.1
            qf = 0.0_dp
            exit
          end if
          qf = qf * C(k)**orden_f(k,j)                              ! fraccionario
        end if
      end do

      qb = 0.0_dp
      if (es_reversible(j)) then
        qb = kb(j)                                                  ! [Ec.25] qb=kb prod C^nu''
        do k = 1, N_ESP
          if (nu_pp(k,j) <= 0.0_dp) cycle
          if (nu_pp(k,j) == 1.0_dp) then
            qb = qb * C(k)
          else if (nu_pp(k,j) == 2.0_dp) then
            qb = qb * C(k) * C(k)
          else
            if (C(k) <= C_FRAC) then
              qb = 0.0_dp
              exit
            end if
            qb = qb * C(k)**nu_pp(k,j)
          end if
        end do
      end if
      qf = qf * Cm                                                  ! [Ec.25] qf*M
      qb = qb * Cm                                                  ! [Ec.25] qb*M

      do k = 1, N_ESP                                               ! [Ec.26] reparto
        if (nu_pp(k,j) > 0.0_dp) then
          prod(k) = prod(k) + nu_pp(k,j) * qf * W_ESP(k)            ! [Ec.26] prod
          dest(k) = dest(k) + nu_pp(k,j) * qb * W_ESP(k)            ! [Ec.26] dest inversa
        end if
        if (nu_p(k,j) > 0.0_dp) then
          dest(k) = dest(k) + nu_p(k,j) * qf * W_ESP(k)             ! [Ec.26] dest directa
          prod(k) = prod(k) + nu_p(k,j) * qb * W_ESP(k)             ! [Ec.26] prod inversa
        end if
      end do
    end do
  end subroutine tasas_produccion

  subroutine k_lindemann(k0v, Cm, kf, kb)                           ! [Ec.24]
    implicit none
    real(dp), intent(in)    :: k0v(N_REAC), Cm(N_REAC)
    real(dp), intent(inout) :: kf(N_REAC), kb(N_REAC)
    integer  :: j
    real(dp) :: k0, Pr, factor

    do j = 1, N_REAC
      if (tipo_reac(j) /= TIPO_LINDEM) cycle                        ! solo Lindemann
      k0 = k0v(j)
      if (kf(j) > 1.0e-300_dp) then
        Pr     = k0 * Cm(j) / kf(j)                                 ! [Ec.24] Pr
        factor = Pr / (1.0_dp + Pr)                                 ! [Ec.24] factor
      else
        factor = 0.0_dp
      end if
      kf(j) = kf(j) * factor                                        ! [Ec.24] k_ef
      kb(j) = kb(j) * factor                                        ! [Ec.24] kb escalado
    end do
  end subroutine k_lindemann

  subroutine tercer_cuerpo(rho, Y, Cm)                              ! [Ec.25] M efectivo
    implicit none
    real(dp), intent(in)  :: rho, Y(N_ESP)
    real(dp), intent(out) :: Cm(N_REAC)
    integer  :: j, k
    real(dp) :: C(N_ESP)

    do k = 1, N_ESP
      C(k) = rho * max(Y(k), 0.0_dp) / W_ESP(k)                     ! [Ec.25] Ck
    end do
    Cm(:) = 1.0_dp
    do j = 1, N_REAC
      if (tipo_reac(j) /= TIPO_TERCERO .and. tipo_reac(j) /= TIPO_LINDEM) cycle
      Cm(j) = 0.0_dp
      do k = 1, N_ESP
        Cm(j) = Cm(j) + eff_M(k,j) * C(k)                           ! [Ec.25]
      end do
    end do
  end subroutine tercer_cuerpo

  subroutine terminos_qss(rho, Y, kf, kb, q, p)                     ! dY/dt = q - p Y
    implicit none
    real(dp), intent(in)  :: rho, Y(N_ESP), kf(N_REAC), kb(N_REAC)
    real(dp), intent(out) :: q(N_ESP), p(N_ESP)
    real(dp) :: prod(N_ESP), dest(N_ESP)
    integer  :: k

    call tasas_produccion(rho, Y, kf, kb, prod, dest)               ! [Ec.25,26]
    do k = 1, N_ESP
      q(k) = prod(k) / rho                                          ! q_k = P_k/rho
      p(k) = dest(k) / (rho * max(Y(k), Y_PISO))                    ! p_k = D_k/(rho Y)
    end do
  end subroutine terminos_qss

  pure function alfa_qss(pdt) result(a)                             ! funcion alfa Mott
    implicit none
    real(dp), intent(in) :: pdt
    real(dp) :: a, r, r2, r3

    if (pdt < 1.0e-12_dp) then;  a = 0.5_dp;  return;  end if       ! no rigido
    r  = 1.0_dp / pdt
    if (r > 1.0e8_dp) then;      a = 0.5_dp;  return;  end if
    r2 = r*r
    r3 = r2*r
    a  = (180.0_dp*r3 + 60.0_dp*r2 + 11.0_dp*r + 1.0_dp) / &
         (360.0_dp*r3 + 60.0_dp*r2 + 12.0_dp*r + 1.0_dp)
  end function alfa_qss

  subroutine paso_quimico(dt_total, rho, Y, T, n_sub)               ! [Ec.27,28,29]
    implicit none
    real(dp), intent(in)    :: dt_total, rho
    real(dp), intent(inout) :: Y(N_ESP), T
    integer,  intent(out)   :: n_sub

    real(dp) :: h_obj, t_hecho, dt_sub, dt_min
    real(dp) :: cp_k(N_ESP), h_k(N_ESP), s_k(N_ESP)
    real(dp) :: kf(N_REAC), kb(N_REAC), k0v(N_REAC), Cm(N_REAC)
    real(dp) :: prod(N_ESP), dest(N_ESP)
    real(dp) :: f0(N_ESP), fp(N_ESP), Jac(N_ESP,N_ESP)
    real(dp) :: Amat(N_ESP,N_ESP), rhs(N_ESP), dY(N_ESP), Yp(N_ESP)
    real(dp) :: dT_est, delta, suma, cp_mix, Y_neg
    integer  :: k, l, info
    logical  :: aceptado

    n_sub   = 0
    t_hecho = 0.0_dp
    dt_min  = dt_total * FRAC_MIN_SUB
    dt_sub  = dt_total
    h_obj = entalpia_mezcla(T, Y)                                   ! [Ec.6] h conservada

    do while (t_hecho < dt_total * (1.0_dp - 1.0e-12_dp))           ! bucle subpasos
      if (n_sub >= MAX_SUBPASOS) exit

      call props_especies(T, cp_k, h_k, s_k)                        ! [Ec.5]
      call k_interpoladas(T, kf, kb, k0v)                           ! [Ec.20,23,24]
      call tercer_cuerpo(rho, Y, Cm)                                ! [Ec.25]
      call k_lindemann(k0v, Cm, kf, kb)                             ! [Ec.24]

      cp_mix = 0.0_dp
      do k = 1, N_ESP
        cp_mix = cp_mix + Y(k) * cp_k(k)                            ! [Ec.6] cp_mix
      end do
      if (cp_mix < 1.0e-8_dp) cp_mix = 1000.0_dp

      call tasas_produccion(rho, Y, kf, kb, prod, dest)             ! [Ec.25,26]
      f0(:) = (prod(:) - dest(:)) / rho                             ! [Ec.27] f(Y)

      do l = 1, N_ESP                                               ! [Ec.28] jacobiano numerico
        delta = max(DY_JAC, 1.0e-6_dp * abs(Y(l)))
        Yp(:) = Y(:)
        Yp(l) = Yp(l) + delta
        call tasas_produccion(rho, Yp, kf, kb, prod, dest)
        fp(:) = (prod(:) - dest(:)) / rho
        do k = 1, N_ESP
          Jac(k,l) = (fp(k) - f0(k)) / delta                        ! [Ec.28] J
        end do
      end do

      aceptado = .false.                                            ! control de paso
      do while (.not. aceptado)
        do l = 1, N_ESP
          do k = 1, N_ESP
            Amat(k,l) = -Jac(k,l)                                   ! [Ec.28] -J
          end do
        end do
        do k = 1, N_ESP
          Amat(k,k) = Amat(k,k) + 1.0_dp / dt_sub                   ! [Ec.28] I/dt
          rhs(k)    = f0(k)                                         ! [Ec.28] f(Y^n)
        end do
        call resolver_lineal(Amat, rhs, dY, N_ESP, info)            ! [Ec.28] dY
        if (info /= 0) then
          dY(:) = f0(:) * dt_sub                                    ! fallback Euler
        end if

        dT_est = 0.0_dp
        Y_neg  = 0.0_dp
        do k = 1, N_ESP
          dT_est = dT_est - h_k(k) * dY(k)                          ! [Ec.6] dT a h=cte
          Y_neg  = min(Y_neg, Y(k) + dY(k))                         ! test positividad
        end do
        dT_est = dT_est / cp_mix

        if ((maxval(abs(dY)) > DY_MAX_SUB .or. abs(dT_est) > DT_MAX_SUB .or. Y_neg < -1.0e-12_dp) .and. dt_sub > dt_min) then
          dt_sub = dt_sub * 0.5_dp                                  ! reducir paso
        else
          aceptado = .true.                                         ! aceptar
        end if
      end do

      do k = 1, N_ESP
        Y(k) = Y(k) + dY(k)                                         ! aplicar dY
        if (Y(k) < 0.0_dp) Y(k) = 0.0_dp
        if (Y(k) > 1.0_dp) Y(k) = 1.0_dp
      end do

      call temperatura_de_entalpia(h_obj, Y, T)                     ! [Ec.7] h->T

      t_hecho = t_hecho + dt_sub
      n_sub   = n_sub + 1

      dt_sub = min(dt_sub * 1.6_dp, dt_total - t_hecho)             ! crecer paso
      if (dt_sub <= 0.0_dp) exit
    end do

    suma = sum(Y)                                                   ! normalizacion
    if (suma > 1.0e-12_dp) then
      Y(:) = Y(:) / suma
    else
      Y(:) = 0.0_dp;  Y(I_O2) = 1.0_dp                              ! celda vacia: O2
    end if
    call temperatura_de_entalpia(h_obj, Y, T)                       ! [Ec.7]
  end subroutine paso_quimico

  subroutine resolver_lineal(A, b, x, n, info)                      ! Gauss con pivoteo
    implicit none
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: A(n,n), b(n)
    real(dp), intent(out)   :: x(n)
    integer,  intent(out)   :: info
    integer  :: i, j, k, piv
    real(dp) :: maxv, factor, tmp

    info = 0
    do k = 1, n-1                                                   ! triangularizacion
      piv  = k
      maxv = abs(A(k,k))
      do i = k+1, n
        if (abs(A(i,k)) > maxv) then
          maxv = abs(A(i,k));  piv = i                              ! pivote parcial
        end if
      end do
      if (maxv < 1.0e-300_dp) then
        info = 1;  x(:) = 0.0_dp;  return                           ! singular
      end if
      if (piv /= k) then                                            ! intercambio
        do j = k, n
          tmp = A(k,j);  A(k,j) = A(piv,j);  A(piv,j) = tmp
        end do
        tmp = b(k);  b(k) = b(piv);  b(piv) = tmp
      end if
      do i = k+1, n
        factor = A(i,k) / A(k,k)
        if (factor == 0.0_dp) cycle
        do j = k, n
          A(i,j) = A(i,j) - factor * A(k,j)                         ! eliminacion
        end do
        b(i) = b(i) - factor * b(k)
      end do
    end do
    if (abs(A(n,n)) < 1.0e-300_dp) then
      info = 1;  x(:) = 0.0_dp;  return
    end if
    do i = n, 1, -1                                                 ! sustitucion atras
      tmp = b(i)
      do j = i+1, n
        tmp = tmp - A(i,j) * x(j)
      end do
      x(i) = tmp / A(i,i)
    end do
  end subroutine resolver_lineal

end module mod_cinetica