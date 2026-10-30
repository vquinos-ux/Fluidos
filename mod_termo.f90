module mod_termo
  use mod_precision
  implicit none

  integer, parameter :: N_ESP = 9                                     ! numero de especies
  integer, parameter :: I_C2H4 = 1                                    ! C2H4 (combustible)
  integer, parameter :: I_O2   = 2                                    ! O2 (oxidante)
  integer, parameter :: I_CO2  = 3                                    ! CO2 (producto)
  integer, parameter :: I_H2O  = 4                                    ! H2O (producto)
  integer, parameter :: I_CO   = 5                                    ! CO (intermedio)
  integer, parameter :: I_H2   = 6                                    ! H2 (intermedio)
  integer, parameter :: I_H    = 7                                    ! H (radical)
  integer, parameter :: I_O    = 8                                    ! O (radical)
  integer, parameter :: I_OH   = 9                                    ! OH (radical)

  character(len=4), parameter :: NOMBRE_ESP(N_ESP) = &               ! nombres
       [ 'C2H4', 'O2  ', 'CO2 ', 'H2O ', 'CO  ', 'H2  ', 'H   ', 'O   ', 'OH  ' ]

  real(dp), parameter :: W_ESP(N_ESP) = [ &                          ! W [kg/kmol]
       28.05316_dp, 31.99880_dp, 44.00950_dp, 18.01528_dp, 28.01010_dp, &
        2.01588_dp,  1.00794_dp, 15.99940_dp, 17.00734_dp ]

  real(dp), parameter :: RU_KMOL = 8314.46261815324_dp               ! R [J/(kmol K)]
  real(dp), parameter :: P_REF   = 1.0e5_dp                          ! p_ref [Pa] [Ec.23]

  integer,  parameter :: N_TABLA = 9601                              ! filas tabla [Ec.4]
  real(dp), parameter :: T_TAB_MIN = 200.0_dp                        ! T_min [K]
  real(dp), parameter :: T_TAB_MAX = 5000.0_dp                       ! T_max [K]
  real(dp), parameter :: T_TAB_D   = 0.5_dp                          ! dT [K] [Ec.4]

  real(dp), save :: T_tab(N_TABLA)                                   ! grilla T [K]
  real(dp), save :: cp_tab(N_TABLA, N_ESP)                           ! cp [J/kgK] [Ec.1]
  real(dp), save :: h_tab (N_TABLA, N_ESP)                           ! h abs [J/kg] [Ec.2]
  real(dp), save :: s_tab (N_TABLA, N_ESP)                           ! s° [J/kgK] [Ec.3]

  real(dp), parameter :: P_OPER = 101325.0_dp                        ! p operacion [Pa] [Ec.9]
  real(dp), parameter :: SIGMA_SB = 5.670374419e-8_dp                ! Stefan-Boltzmann [Ec.19]

  real(dp), parameter :: LJ_SIG(N_ESP) = [ &                         ! sigma LJ [A] [Ec.10]
       3.971_dp, 3.458_dp, 3.763_dp, 2.605_dp, 3.650_dp, &
       2.920_dp, 2.050_dp, 2.750_dp, 2.750_dp ]
  real(dp), parameter :: LJ_EPS(N_ESP) = [ &                         ! eps/k [K] [Ec.10]
       280.800_dp, 107.400_dp, 244.000_dp, 572.400_dp, 98.100_dp, &
        38.000_dp, 145.000_dp,  80.000_dp,  80.000_dp ]

  integer,  parameter :: N_PAR = N_ESP*(N_ESP-1)/2                   ! 36 pares binarios
  real(dp), save :: mu_tab (N_TABLA, N_ESP)                          ! mu [Pa s] [Ec.10]
  real(dp), save :: lam_tab(N_TABLA, N_ESP)                          ! lam [W/mK] [Ec.12]
  real(dp), save :: Dij_tab(N_TABLA, N_PAR)                          ! D bin [m2/s] [Ec.15]
  integer,  save :: idx_par(N_ESP, N_ESP)                            ! mapa (k,j)->par

  real(dp), save :: wilke_A(N_ESP,N_ESP)                             ! (Wj/Wk)^(1/4) [Ec.13]
  real(dp), save :: wilke_B(N_ESP,N_ESP)                             ! 1/sqrt(8(1+Wk/Wj)) [Ec.13]

contains

  subroutine cargar_tablas()                                          ! lee NASA [Ec.1-3]
    implicit none
    integer  :: k, i, unidad, io
    character(len=200) :: linea
    character(len=16)  :: archivo
    real(dp) :: d1, d2, d3

    do k = 1, N_ESP
      archivo = trim(adjustl(NOMBRE_ESP(k)))//'.txt'                 ! <ESPECIE>.txt
      open(newunit=unidad, file=archivo, status='old', action='read', iostat=io)
      if (io /= 0) then
        print *, 'ERROR: no se pudo abrir la tabla ', trim(archivo)
        print *, 'Genera las tablas primero con generar_tablas.f90'
        stop
      end if
      do                                                              ! saltar '#'
        read(unidad,'(A)', iostat=io) linea
        if (io /= 0) exit
        if (len_trim(linea) > 0) then
          if (linea(1:1) /= '#') exit
        end if
      end do
      backspace(unidad)                                               ! ultima leida era dato
      do i = 1, N_TABLA                                               ! leer filas
        read(unidad, *, iostat=io) T_tab(i), cp_tab(i,k), h_tab(i,k), s_tab(i,k), d1, d2, d3
        if (io /= 0) then
          print *, 'ERROR: tabla incompleta en ', trim(archivo), ' fila ', i
          stop
        end if
      end do
      close(unidad)
    end do
    print *, 'Tablas termodinamicas cargadas: ', N_ESP, ' especies x ', N_TABLA, ' puntos'
    call init_transporte()
  end subroutine cargar_tablas

  pure function omega_mu(T_red) result(om)                            ! [Ec.11]
    implicit none
    real(dp), intent(in) :: T_red
    real(dp) :: om, Tr
    Tr = max(T_red, 0.3_dp)
    om = 1.16145_dp*Tr**(-0.14874_dp) + 0.52487_dp*exp(-0.77320_dp*Tr) &  ! [Ec.11]
       + 2.16178_dp*exp(-2.43787_dp*Tr)
  end function omega_mu

  pure function omega_D(T_red) result(om)                             ! [Ec.17]
    implicit none
    real(dp), intent(in) :: T_red
    real(dp) :: om, Tr
    Tr = max(T_red, 0.3_dp)
    om = 1.06036_dp/Tr**0.15610_dp + 0.19300_dp/exp(0.47635_dp*Tr) &  ! [Ec.17]
       + 1.03587_dp/exp(1.52996_dp*Tr) + 1.76474_dp/exp(3.89411_dp*Tr)
  end function omega_D

  subroutine init_transporte()                                        ! tabula mu, lam, D
    implicit none
    integer  :: i, k, j, ipar
    real(dp) :: T, Tred, om, sig_kj, eps_kj, Mfac, Patm

    Patm = P_OPER / 101325.0_dp                                     ! p [atm] [Ec.15]

    do k = 1, N_ESP
      do i = 1, N_TABLA
        T    = T_tab(i)
        Tred = T / LJ_EPS(k)                                        ! T* [Ec.10]
        om   = omega_mu(Tred)
        mu_tab(i,k)  = 2.6693e-6_dp * sqrt(W_ESP(k)*T) / (LJ_SIG(k)**2 * om)  ! [Ec.10]
        lam_tab(i,k) = mu_tab(i,k) * (1.32_dp*cp_tab(i,k) + 0.45_dp*RU_KMOL/W_ESP(k))  ! [Ec.12]
      end do
    end do

    idx_par(:,:) = 0
    ipar = 0
    do k = 1, N_ESP-1
      do j = k+1, N_ESP
        ipar = ipar + 1
        idx_par(k,j) = ipar
        idx_par(j,k) = ipar
        sig_kj = 0.5_dp * (LJ_SIG(k) + LJ_SIG(j))                   ! [Ec.16]
        eps_kj = sqrt(LJ_EPS(k) * LJ_EPS(j))                        ! [Ec.16]
        Mfac   = 1.0_dp/W_ESP(k) + 1.0_dp/W_ESP(j)                  ! [Ec.15]
        do i = 1, N_TABLA
          T  = T_tab(i)
          om = omega_D(T/eps_kj)
          Dij_tab(i,ipar) = 1.0e-4_dp * 0.0018583_dp * sqrt(T**3 * Mfac) &  ! [Ec.15] (1e-4: cm2/s->m2/s)
                          / (Patm * sig_kj**2 * om)
        end do
      end do
    end do

    do k = 1, N_ESP
      do j = 1, N_ESP
        wilke_A(k,j) = (W_ESP(j)/W_ESP(k))**0.25_dp                 ! [Ec.13]
        wilke_B(k,j) = 1.0_dp / sqrt(8.0_dp*(1.0_dp + W_ESP(k)/W_ESP(j)))  ! [Ec.13]
      end do
    end do
    print *, 'Tablas de transporte generadas: mu, lambda y ', N_PAR, ' pares binarios'
  end subroutine init_transporte

  subroutine propiedades_mezcla(T, Y, mu_mix, lam_mix, cp_mix, D_k)  ! [Ec.13,14,18]
    implicit none
    real(dp), intent(in)  :: T, Y(N_ESP)
    real(dp), intent(out) :: mu_mix, lam_mix, cp_mix, D_k(N_ESP)

    real(dp) :: X(N_ESP), mu_k(N_ESP), lam_k(N_ESP), cp_k(N_ESP), Dij(N_PAR)
    real(dp) :: W_mix, frac, suma, den, sarit, sarmo, phi
    integer  :: lo, k, j, ipar

    call indice_T(T, lo, frac)                                      ! [Ec.4]
    do k = 1, N_ESP
      mu_k (k) = mu_tab (lo,k) + frac*(mu_tab (lo+1,k) - mu_tab (lo,k))  ! [Ec.5]
      lam_k(k) = lam_tab(lo,k) + frac*(lam_tab(lo+1,k) - lam_tab(lo,k))  ! [Ec.5]
      cp_k (k) = cp_tab (lo,k) + frac*(cp_tab (lo+1,k) - cp_tab (lo,k))  ! [Ec.5]
    end do
    do ipar = 1, N_PAR
      Dij(ipar) = Dij_tab(lo,ipar) + frac*(Dij_tab(lo+1,ipar) - Dij_tab(lo,ipar))  ! [Ec.5]
    end do

    W_mix  = peso_molecular_mezcla(Y)                               ! [Ec.8]
    cp_mix = 0.0_dp
    do k = 1, N_ESP
      X(k)   = max(Y(k), 0.0_dp) * W_mix / W_ESP(k)                 ! [Ec.8]
      cp_mix = cp_mix + Y(k)*cp_k(k)                                ! [Ec.6]
    end do
    if (cp_mix < 1.0e-8_dp) cp_mix = 1000.0_dp

    mu_mix = 0.0_dp                                                 ! [Ec.13]
    do k = 1, N_ESP
      if (X(k) < 1.0e-12_dp) cycle
      den = 0.0_dp
      do j = 1, N_ESP
        if (X(j) < 1.0e-12_dp) cycle
        phi = (1.0_dp + sqrt(mu_k(k)/mu_k(j))*wilke_A(k,j))**2 * wilke_B(k,j)  ! [Ec.13]
        den = den + X(j)*phi
      end do
      if (den > 1.0e-30_dp) mu_mix = mu_mix + X(k)*mu_k(k)/den      ! [Ec.13]
    end do
    if (mu_mix < 1.0e-12_dp) mu_mix = mu_k(I_O2)

    sarit = 0.0_dp
    sarmo = 0.0_dp
    do k = 1, N_ESP
      if (X(k) < 1.0e-12_dp) cycle
      sarit = sarit + X(k)*lam_k(k)                                 ! [Ec.14]
      sarmo = sarmo + X(k)/lam_k(k)                                 ! [Ec.14]
    end do
    if (sarmo > 1.0e-30_dp) then
      lam_mix = 0.5_dp*(sarit + 1.0_dp/sarmo)                       ! [Ec.14]
    else
      lam_mix = lam_k(I_O2)
    end if

    do k = 1, N_ESP                                                 ! [Ec.18]
      suma = 0.0_dp
      do j = 1, N_ESP
        if (j == k) cycle
        if (X(j) < 1.0e-12_dp) cycle
        suma = suma + X(j)/Dij(idx_par(k,j))                        ! [Ec.18]
      end do
      if (suma > 1.0e-30_dp) then
        D_k(k) = (1.0_dp - min(Y(k),0.999_dp)) / suma               ! [Ec.18]
      else
        if (k == I_O2) then
          D_k(k) = Dij(idx_par(I_O2,I_H2O))                         ! fallback
        else
          D_k(k) = Dij(idx_par(k,I_O2))                             ! fallback
        end if
      end if
    end do
  end subroutine propiedades_mezcla

  function q_radiacion(T, Y, T_inf) result(q)                        ! [Ec.19]
    implicit none
    real(dp), intent(in) :: T, Y(N_ESP), T_inf
    real(dp) :: q
    real(dp) :: th, ap_CO2, ap_H2O, ap_CO, kap, W_mix, Patm
    real(dp) :: X_CO2, X_H2O, X_CO

    if (T <= T_inf) then
      q = 0.0_dp
      return
    end if
    th = 1000.0_dp / T                                              ! theta = 1000/T

    ap_CO2 =  18.741_dp - 121.310_dp*th + 273.500_dp*th**2 &        ! RADCAL CO2
            - 194.050_dp*th**3 + 56.310_dp*th**4 - 5.8169_dp*th**5
    ap_H2O = -0.23093_dp - 1.12390_dp*th + 9.41530_dp*th**2 &       ! RADCAL H2O
            - 2.99880_dp*th**3 + 0.51382_dp*th**4 - 1.86840e-5_dp*th**5
    if (T < 750.0_dp) then
      ap_CO = 4.7869_dp - 0.06953_dp*T + 2.95775e-4_dp*T**2 &       ! RADCAL CO (bajo)
            - 4.25732e-7_dp*T**3 + 2.02894e-10_dp*T**4
    else
      ap_CO = 10.09_dp - 0.01183_dp*T + 4.7753e-6_dp*T**2 &         ! RADCAL CO (alto)
            - 5.87209e-10_dp*T**3 - 2.5334e-14_dp*T**4
    end if
    ap_CO2 = max(ap_CO2, 0.0_dp)                                    ! cortar negativos
    ap_H2O = max(ap_H2O, 0.0_dp)
    ap_CO  = max(ap_CO , 0.0_dp)

    W_mix = peso_molecular_mezcla(Y)                                ! [Ec.8]
    X_CO2 = max(Y(I_CO2),0.0_dp) * W_mix / W_ESP(I_CO2)             ! [Ec.8]
    X_H2O = max(Y(I_H2O),0.0_dp) * W_mix / W_ESP(I_H2O)             ! [Ec.8]
    X_CO  = max(Y(I_CO ),0.0_dp) * W_mix / W_ESP(I_CO)              ! [Ec.8]
    Patm = P_OPER / 101325.0_dp
    kap  = (X_CO2*ap_CO2 + X_H2O*ap_H2O + X_CO*ap_CO) * Patm        ! [Ec.19]: kappa_p
    q = 4.0_dp * SIGMA_SB * kap * (T**4 - T_inf**4)                 ! [Ec.19]
  end function q_radiacion

  subroutine indice_T(T, lo, frac)                                  ! [Ec.4]
    implicit none
    real(dp), intent(in)  :: T
    integer,  intent(out) :: lo
    real(dp), intent(out) :: frac
    real(dp) :: pos

    pos = (T - T_TAB_MIN) / T_TAB_D + 1.0_dp                       ! [Ec.4]: pi
    lo  = int(pos)                                                  ! [Ec.4]: l = [pi]
    if (lo < 1) then                                                ! saturacion abajo
      lo = 1;            frac = 0.0_dp
    else if (lo >= N_TABLA) then                                    ! saturacion arriba
      lo = N_TABLA - 1;  frac = 1.0_dp
    else
      frac = pos - real(lo, dp)                                     ! [Ec.4]: theta
    end if
  end subroutine indice_T

  subroutine props_especies(T, cp_k, h_k, s_k)                      ! [Ec.5] en cp,h,s
    implicit none
    real(dp), intent(in)  :: T
    real(dp), intent(out) :: cp_k(N_ESP), h_k(N_ESP), s_k(N_ESP)
    integer  :: lo, k
    real(dp) :: frac

    call indice_T(T, lo, frac)                                      ! [Ec.4]
    do k = 1, N_ESP
      cp_k(k) = cp_tab(lo,k) + frac * (cp_tab(lo+1,k) - cp_tab(lo,k))  ! [Ec.5]
      h_k (k) = h_tab (lo,k) + frac * (h_tab (lo+1,k) - h_tab (lo,k))  ! [Ec.5]
      s_k (k) = s_tab (lo,k) + frac * (s_tab (lo+1,k) - s_tab (lo,k))  ! [Ec.5]
    end do
  end subroutine props_especies

  function entalpia_mezcla(T, Y) result(h)                          ! [Ec.6]
    implicit none
    real(dp), intent(in) :: T, Y(N_ESP)
    real(dp) :: h
    integer  :: lo, k
    real(dp) :: frac, hk

    call indice_T(T, lo, frac)                                      ! [Ec.4]
    h = 0.0_dp
    do k = 1, N_ESP
      hk = h_tab(lo,k) + frac * (h_tab(lo+1,k) - h_tab(lo,k))       ! [Ec.5]
      h  = h + Y(k) * hk                                            ! [Ec.6]
    end do
  end function entalpia_mezcla

  function cp_mezcla(T, Y) result(cp)                               ! [Ec.6]
    implicit none
    real(dp), intent(in) :: T, Y(N_ESP)
    real(dp) :: cp
    integer  :: lo, k
    real(dp) :: frac, cpk

    call indice_T(T, lo, frac)                                      ! [Ec.4]
    cp = 0.0_dp
    do k = 1, N_ESP
      cpk = cp_tab(lo,k) + frac * (cp_tab(lo+1,k) - cp_tab(lo,k))   ! [Ec.5]
      cp  = cp + Y(k) * cpk                                         ! [Ec.6]
    end do
    if (cp < 1.0e-8_dp) cp = 1000.0_dp                              ! proteccion
  end function cp_mezcla

  subroutine temperatura_de_entalpia(h_obj, Y, T)                    ! [Ec.7] Newton
    implicit none
    real(dp), intent(in)    :: h_obj, Y(N_ESP)
    real(dp), intent(inout) :: T
    integer  :: it
    real(dp) :: h_act, cp_act, dT, h_lo, h_hi

    h_lo = entalpia_mezcla(T_TAB_MIN, Y)                            ! [Ec.6] saturacion
    h_hi = entalpia_mezcla(T_TAB_MAX, Y)
    if (h_obj <= h_lo) then;  T = T_TAB_MIN;  return;  end if
    if (h_obj >= h_hi) then;  T = T_TAB_MAX;  return;  end if
    if (T < T_TAB_MIN .or. T > T_TAB_MAX .or. T /= T) T = 1000.0_dp  ! semilla

    do it = 1, 30
      h_act  = entalpia_mezcla(T, Y)                                ! [Ec.6]
      cp_act = cp_mezcla(T, Y)                                      ! [Ec.6]
      dT     = (h_obj - h_act) / cp_act                             ! [Ec.7]
      if (dT >  400.0_dp) dT =  400.0_dp                            ! [Ec.7]
      if (dT < -400.0_dp) dT = -400.0_dp                            ! [Ec.7]
      T = T + dT
      if (T < T_TAB_MIN) T = T_TAB_MIN
      if (T > T_TAB_MAX) T = T_TAB_MAX
      if (abs(dT) < 1.0e-8_dp) exit
    end do
  end subroutine temperatura_de_entalpia

  function peso_molecular_mezcla(Y) result(W_mix)                   ! [Ec.8]
    implicit none
    real(dp), intent(in) :: Y(N_ESP)
    real(dp) :: W_mix, suma
    integer  :: k

    suma = 0.0_dp
    do k = 1, N_ESP
      suma = suma + Y(k) / W_ESP(k)                                 ! [Ec.8]
    end do
    if (suma < 1.0e-12_dp) suma = 1.0_dp / W_ESP(I_O2)              ! fallback
    W_mix = 1.0_dp / suma                                           ! [Ec.8]
  end function peso_molecular_mezcla

end module mod_termo