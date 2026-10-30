module mod_malla
  use mod_precision
  use mod_termo
  implicit none

  integer, parameter :: NX  = 241                                     ! nodos x
  integer, parameter :: NYS = 41                                      ! nodos solido; j=NYS: superficie
  integer, parameter :: NYG = 101                                     ! nodos gas (desde j=NYS)
  integer, parameter :: NY  = NYS+NYG-1                               ! filas totales
  integer, parameter :: NN  = NX*NY                                   ! nodos totales
  integer, parameter :: NTRI= (NX-1)*(NY-1)*2                         ! triangulos (dibujo)

  real(dp), parameter :: Lx  = 0.10_dp                                ! largo camara [m]
  real(dp), parameter :: e0  = 0.02_dp                                ! espesor inicial [m]

  real(dp), parameter :: L_precam = 0.025_dp                          ! precamara [m]
  real(dp), parameter :: L_poscam = 0.020_dp                          ! poscamara [m]
  real(dp), parameter :: x_ini_gr = L_precam                          ! inicio grano [m]
  real(dp), parameter :: x_fin_gr = Lx - L_poscam                     ! fin grano [m]
  real(dp), parameter :: L_grano  = x_fin_gr - x_ini_gr               ! largo grano [m]
  real(dp), parameter :: e_pared  = 2.0e-3_dp                         ! liner inerte [m]

  real(dp), parameter :: L_RAMPA_IN  = 2.0e-3_dp                      ! rampa entrada [m] [Ec.36]
  real(dp), parameter :: L_RAMPA_OUT = 10.0e-3_dp                     ! rampa salida [m]  [Ec.36]
  real(dp), parameter :: L_RAMPA = L_RAMPA_IN                         ! alias

  real(dp), parameter :: Hp0 = 20.0e-3_dp                             ! puerto inicial [m]
  real(dp), parameter :: Ly  = e0 + Hp0                               ! alto camara [m]
  real(dp), parameter :: Hp_max = 36.0e-3_dp                          ! corte por puerto [m]
  real(dp), parameter :: e_min  = e_pared                             ! agotamiento [m] [Ec.45]

  real(dp), parameter :: pot_x = 1.25_dp                              ! concentracion x (pre Ec.37)
  real(dp), parameter :: b_gas = 2.5_dp                               ! concentracion gas-pared (pre Ec.37)
  real(dp), parameter :: b_sol = 5.0_dp                               ! concentracion solido-pared (pre Ec.37)

  integer,  parameter :: modo_entrada = 3                             ! 1=mdot, 2=U, 3=G_ox
  real(dp), parameter :: mdot_ox   = 2.9e-3_dp                        ! mdot O2 total [kg/s]
  real(dp), parameter :: ancho_z   = 0.020_dp                         ! profundidad [m]
  real(dp), parameter :: U_in_fijo = 8.0_dp                           ! U entrada [m/s]
  real(dp), parameter :: G_ox_fijo = 14.0_dp                          ! G_ox impuesto [kg/m2s]

  logical,  parameter :: usar_ranura = .true.                         ! ranura baja
  real(dp), parameter :: y_ranura    = 0.5_dp*e0                      ! techo ranura [m]
  real(dp), parameter :: ang_entrada = 0.0_dp                         ! inclinacion [grados]
  real(dp), parameter :: PI_M = 3.14159265358979323846_dp             ! pi

  real(dp), parameter :: T_in  = 300.0_dp                             ! T entrada [K]
  real(dp), parameter :: T_amb = 300.0_dp                             ! T ambiente [K]
  real(dp), parameter :: g_y   = -9.81_dp                             ! gravedad [m/s2] [Ec.32]

  real(dp), parameter :: rho_s  = 960.0_dp                            ! densidad [kg/m3] [Ec.40]
  real(dp), parameter :: cp_s   = 1255.2_dp                           ! cp solido [J/kgK] [Ec.43]
  real(dp), parameter :: lam_s  = 0.1549_dp                           ! lambda solido [W/mK] [Ec.43]
  real(dp), parameter :: T_fus  = 403.0_dp                            ! fusion [K] [Ec.44]
  real(dp), parameter :: dT_fus = 15.0_dp                             ! semiancho fusion [K] [Ec.44]
  real(dp), parameter :: L_fus  = 2.10e5_dp                           ! latente [J/kg] [Ec.44]
  real(dp), parameter :: dh_pyr = 2.72e6_dp                           ! calor pirolisis [J/kg] [Ec.41]
  real(dp), parameter :: A_pyr  = 4.78e3_dp                           ! A pirolisis [m/s] [Ec.40]
  real(dp), parameter :: Ea_pyr = 2.5104e5_dp                         ! Ea [J/mol] [Ec.40]
  real(dp), parameter :: emis_w = 0.90_dp                             ! emisividad [Ec.41]

  integer,  parameter :: modo_chispa = 2                              ! 1=pegada, 2=caja
  real(dp), parameter :: t_chispa_ini = 0.02_dp                       ! [s]
  real(dp), parameter :: t_chispa_fin = 0.05_dp                       ! [s]
  real(dp), parameter :: x_chispa_0 = x_ini_gr - 0.002_dp             ! banda x [m]
  real(dp), parameter :: x_chispa_1 = x_ini_gr + 0.002_dp
  real(dp), parameter :: dy_chispa  = 0.003_dp                        ! alto sobre pared [m]
  real(dp), parameter :: y_chispa_0 = 0.004_dp                        ! piso caja [m]
  real(dp), parameter :: y_chispa_1 = 0.008_dp                        ! techo caja [m]
  real(dp), parameter :: P_chispa   = 2.0e9_dp                        ! potencia [W/m3]
  real(dp), parameter :: q_ign     = 1.5e6_dp                         ! flujo pared [W/m2] [Ec.41]
  real(dp), parameter :: t_ign_fin = 1.50_dp                          ! [s]
  real(dp), parameter :: x_ign_0 = x_ini_gr                           ! banda x [m]
  real(dp), parameter :: x_ign_1 = x_ini_gr + 0.023_dp

  logical,  parameter :: tensor_completo = .true.                     ! tau_ij completo [Ec.33]
  logical,  parameter :: usar_radiacion = .true.                      ! [Ec.19] activo
  real(dp), parameter :: C_abs = 0.60_dp                              ! absorcion [1/m] [Ec.19]
  real(dp), parameter :: RU_MOLAR = 8.31446261815324_dp               ! R [J/(mol K)] [Ec.40]

  real(dp), parameter :: CFL     = 0.35_dp                            ! Courant [Ec.50]
  real(dp), parameter :: FOU     = 0.20_dp                            ! Fourier [Ec.51]
  real(dp), parameter :: dt_max  = 2.0e-5_dp                          ! tope [s] [Ec.52]
  real(dp), parameter :: t_final = 120.0_dp                           ! t final [s]
  integer,  parameter :: iter_max = 200000000  !aqui los 2 min        ! iteraciones top
  integer,  parameter :: n_salida = 2000                              ! pasos entre volcados

  real(dp), parameter :: FAC_SOLIDO = 20.0_dp                         ! dt_sol = FAC*dt_gas
  real(dp), parameter :: t_acel_ini = 0.0_dp                          ! antes FAC=1

  integer,  parameter :: it_poisson  = 60                             ! barridos linea [Ec.48]
  real(dp), parameter :: tol_poisson = 1.0e-4_dp                      ! tol residuo [Ec.48]

  real(dp), save :: xg(NX), yg(NX,NY), dYdx(NX,NY)                    ! malla [Ec.37]; metrica [Ec.38]
  real(dp), save :: fsol(NYS), fgas(NYG)                              ! perfiles y [Ec.37]
  real(dp), save :: u(NX,NY), v(NX,NY), pr(NX,NY)                     ! vel; pr [Ec.48]
  real(dp), save :: ru_(NX,NY), rv_(NX,NY)                            ! rho u n+1 [Ec.32]
  real(dp), save :: rus(NX,NY), rvs(NX,NY)                            ! rho u* [Ec.32]
  real(dp), save :: Tg(NX,NY), hmix(NX,NY)                            ! T, h_mix [Ec.6]
  real(dp), save :: rho(NX,NY), rho_o(NX,NY)                          ! rho; previo [Ec.9]
  real(dp), save :: Ysp(NX,NY,N_ESP)                                  ! Y_k [Ec.34]
  real(dp), save :: Ynew(NX,NY,N_ESP), hnew(NX,NY)                    ! buffers
  real(dp), save :: mu_(NX,NY), lam_(NX,NY)                           ! mu, lam [Ec.13,14]
  real(dp), save :: cp_(NX,NY), Dm_(NX,NY)                            ! cp; rho D [Ec.6]
  real(dp), save :: qvol(NX,NY)                                       ! fuente rad [Ec.19]
  real(dp), save :: tau_a(NX,NY), tau_b(NX,NY)                        ! mu du/dx, dy [Ec.33]
  real(dp), save :: tau_c(NX,NY), tau_d(NX,NY)                        ! mu dv/dx, dy [Ec.33]
  real(dp), save :: mu43(NX,NY)                                       ! (4/3) mu [Ec.33]
  real(dp), save :: div_(NX,NY)                                       ! div rho u [Ec.48]

  real(dp), save :: Ts(NX,NYS)                                        ! T solido [Ec.43]
  real(dp), save :: esp(NX)                                           ! e(x,t) [Ec.36,45]
  real(dp), save :: Tw(NX)                                            ! T pared [Ec.41]
  real(dp), save :: rdot(NX)                                          ! r [Ec.40]
  real(dp), save :: mdot(NX)                                          ! m'' [Ec.40]
  real(dp), save :: qconv(NX)                                         ! q conv [Ec.41]
  real(dp), save :: qrad(NX)                                          ! q rad [Ec.41]
  logical,  save :: agotado(NX)                                       ! agotada [Ec.45]
  logical,  save :: es_grano(NX)                                      ! con combustible

  real(dp), save :: tiempo, dt                                        ! t, paso [s]
  real(dp), save :: U_ent, G_ox                                       ! entrada U, G_ox
  real(dp), save :: cos_ent = 1.0_dp, sin_ent = 0.0_dp                ! cos/sin ang_entrada
  real(dp), save :: rho_in                                            ! rho O2 [Ec.9]
  real(dp), save :: wk1(N_ESP,N_ESP), wk2(N_ESP,N_ESP)                ! Wilke [Ec.13]

  integer, save :: n_sub_max, n_sub_tot, n_celdas_q                   ! diag quimica [Ec.27]

contains

  pure function perfil_espesor(x) result(e)                          ! [Ec.36]
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: e, t, s1, s2
    t  = (x - x_ini_gr) / L_RAMPA_IN
    t  = min(max(t, 0.0_dp), 1.0_dp)
    s1 = t*t*(3.0_dp - 2.0_dp*t)                                      ! smoothstep [Ec.36]
    t  = (x_fin_gr - x) / L_RAMPA_OUT
    t  = min(max(t, 0.0_dp), 1.0_dp)
    s2 = t*t*(3.0_dp - 2.0_dp*t)                                      ! smoothstep [Ec.36]
    e = e_pared + (e0 - e_pared)*s1*s2                                ! [Ec.36]
  end function perfil_espesor

  subroutine init_wilke()                                              ! [Ec.13]
    implicit none
    integer :: i, j
    do j = 1, N_ESP
      do i = 1, N_ESP
        wk1(i,j) = (W_ESP(j)/W_ESP(i))**0.25_dp                       ! [Ec.13]
        wk2(i,j) = 1.0_dp / sqrt(8.0_dp*(1.0_dp + W_ESP(i)/W_ESP(j))) ! [Ec.13]
      end do
    end do
  end subroutine init_wilke

end module mod_malla