program verificar_entrada
  use mod_precision
  use mod_termo
  use mod_malla
  use mod_solver
  implicit none

  integer  :: i, j, n_abierto, n_ciego, j_tope                 ! contadores
  real(dp) :: t, h_abierta, A_abierta, mdot_tot, y_lo, y_hi    ! locales
  real(dp) :: x_impacto                                        ! impacto chorro-techo

  do j = 1, NYS                                               ! reconstruye fsol [Ec.37]
     t = real(j-1,dp)/real(NYS-1,dp)
     fsol(j) = 1.0_dp - (exp(b_sol*(1.0_dp-t))-1.0_dp)/(exp(b_sol)-1.0_dp)
  end do
  fsol(1) = 0.0_dp;  fsol(NYS) = 1.0_dp                       ! extremos exactos
  do j = 1, NYG                                               ! reconstruye fgas [Ec.37]
     t = real(j-1,dp)/real(NYG-1,dp)
     fgas(j) = (exp(b_gas*t)-1.0_dp)/(exp(b_gas)-1.0_dp)
  end do
  fgas(1) = 0.0_dp;  fgas(NYG) = 1.0_dp                       ! extremos exactos

  do i = 1, NX                                                ! x concentrado + e(x)
     t = real(i-1,dp)/real(NX-1,dp)
     xg(i)  = Lx*t**pot_x                                     ! mapeo pre-Ec.37
     esp(i) = perfil_espesor(xg(i))                           ! [Ec.36] e(x)
  end do

  call construir_malla()                                      ! [Ec.37,38]
  call actualizar_entrada()                                   ! U_ent, G_ox, cos/sin

  n_abierto = 0;  n_ciego = 0;  j_tope = NYS                  ! recorre i=1
  y_lo = yg(1,NYS);  y_hi = yg(1,NYS)                         ! extiende y
  do j = NYS, NY
     if (entrada_abierta(j)) then                             ! ranura: y<=y_ranura
        n_abierto = n_abierto + 1
        j_tope = j
        y_hi = yg(1,j)
     else
        n_ciego = n_ciego + 1                                 ! pared ciega
     end if
  end do

  h_abierta = y_hi - y_lo                                     ! alto de ranura [m]
  A_abierta = h_abierta*ancho_z                               ! area de paso [m2]
  mdot_tot  = rho_in*U_ent*cos_ent*A_abierta                  ! [Ec.9] rho*U*A

  write(*,'(A)') '====================================================='
  write(*,'(A)') ' VERIFICACION DE LA CONDICION DE ENTRADA'
  write(*,'(A)') '====================================================='
  write(*,'(A,L2)')      ' usar_ranura                  : ', usar_ranura
  write(*,'(A,F8.2,A)')  ' ang_entrada                  : ', ang_entrada, ' grados'
  write(*,'(A,I2)')      ' modo_entrada                 : ', modo_entrada
  write(*,'(A)') ' -----------------------------------------------------'
  write(*,'(A,F8.4)')    ' cos(ang)                     : ', cos_ent
  write(*,'(A,F8.4)')    ' sin(ang)                     : ', sin_ent
  write(*,'(A,F8.3,A)')  ' U_ent (modulo del chorro)    : ', U_ent, ' m/s'
  write(*,'(A,F8.3,A)')  '   componente axial    u      : ', U_ent*cos_ent, ' m/s'
  write(*,'(A,F8.3,A)')  '   componente vertical v      : ', U_ent*sin_ent, ' m/s'
  write(*,'(A,F8.3,A)')  ' rho_in                       : ', rho_in, ' kg/m3'
  write(*,'(A,F8.2,A)')  ' G_ox (reportado por el codigo): ', G_ox, ' kg/(m2 s)'
  write(*,'(A)') ' -----------------------------------------------------'
  write(*,'(A,I4,A,I4)') ' nodos ABIERTOS en i=1        : ', n_abierto, ' de ', NYG
  write(*,'(A,I4)')      ' nodos de PARED CIEGA en i=1  : ', n_ciego
  write(*,'(A,F8.3,A)')  ' piso de la entrada  y        : ', y_lo*1000.0_dp, ' mm'
  write(*,'(A,F8.3,A)')  ' techo de la entrada y        : ', y_hi*1000.0_dp, ' mm'
  write(*,'(A,F8.3,A)')  ' alto abierto al flujo        : ', h_abierta*1000.0_dp, ' mm'
  write(*,'(A,F10.6,A)') ' area de paso                 : ', A_abierta*1.0e6_dp, ' mm2'
  write(*,'(A,F8.3,A)')  ' caudal masico total          : ', mdot_tot*1000.0_dp, ' g/s'
  write(*,'(A)') ' -----------------------------------------------------'
  write(*,'(A,F8.3,A)')  ' alto del grano e0            : ', e0*1000.0_dp, ' mm'
  write(*,'(A,F8.3,A)')  ' techo de la camara Ly        : ', Ly*1000.0_dp, ' mm'
  if (sin_ent > 1.0e-6_dp) then                               ! chorro inclinado
     x_impacto = (Ly - y_hi)*cos_ent/sin_ent                  ! geometria: x toca techo
     write(*,'(A,F8.2,A)') ' el chorro toca el techo en x : ', x_impacto*1000.0_dp, ' mm'
     if (x_impacto < x_ini_gr) then                           ! antes del grano
        write(*,'(A)') '   -> ANTES del grano (x_ini_gr = 25 mm): el combustible'
        write(*,'(A)') '      recibe flujo ya rebotado, no el chorro directo.'
     else if (x_impacto < x_fin_gr) then                      ! sobre el grano
        write(*,'(A)') '   -> SOBRE el grano: el chorro impacta el techo en'
        write(*,'(A)') '      pleno tramo de combustible.'
     else                                                     ! poscamara
        write(*,'(A)') '   -> despues del grano, ya en la poscamara.'
     end if
  else
     write(*,'(A)') ' chorro horizontal: no toca el techo.'
  end if
  write(*,'(A,F8.4)')    ' Mach de entrada ~            : ', U_ent/430.0_dp  ! [Ec.30] Ma

  write(*,'(A)') ' -----------------------------------------------------'
  write(*,'(A,I2)')     ' modo_chispa                  : ', modo_chispa
  write(*,'(A,F8.2,A,F8.2,A)') ' banda x de la chispa         : ', &
       x_chispa_0*1000.0_dp, ' a ', x_chispa_1*1000.0_dp, ' mm'
  if (modo_chispa == 2) then                                  ! caja fija
     write(*,'(A,F8.2,A,F8.2,A)') ' banda y de la chispa (caja)  : ', &
          y_chispa_0*1000.0_dp, ' a ', y_chispa_1*1000.0_dp, ' mm'
  else                                                        ! pegada a pared
     write(*,'(A,F8.2,A)') ' alto sobre la superficie     : ', &
          dy_chispa*1000.0_dp, ' mm'
  end if

  n_abierto = 0;  y_lo = 1.0e9_dp;  y_hi = -1.0e9_dp         ! recorre caja chispa
  do i = 1, NX
     if (xg(i) < x_chispa_0 .or. xg(i) > x_chispa_1) cycle   ! fuera de banda x
     do j = NYS+1, NY
        if (modo_chispa == 2) then                            ! caja absoluta
           if (yg(i,j) < y_chispa_0) cycle                    ! bajo piso
           if (yg(i,j) > y_chispa_1) exit                     ! sobre techo
        else                                                  ! pegada
           if (yg(i,j) - esp(i) > dy_chispa) exit             ! alto agotado
        end if
        n_abierto = n_abierto + 1                             ! nodo de gas encendido
        y_lo = min(y_lo, yg(i,j));  y_hi = max(y_hi, yg(i,j))
     end do
  end do
  write(*,'(A,I6)')     ' nodos de GAS que enciende    : ', n_abierto
  if (n_abierto > 0) then
     write(*,'(A,F8.2,A,F8.2,A)') ' altura real que cubre        : ', &
          y_lo*1000.0_dp, ' a ', y_hi*1000.0_dp, ' mm'
     if (usar_ranura) then                                    ! chequeo coherencia
        if (y_lo > y_ranura) then                             ! chispa sobre ranura
           write(*,'(A)') '   -> AVISO: la chispa queda POR ENCIMA de la ranura.'
           write(*,'(A)') '      El chorro de O2 no la atraviesa.'
        else                                                  ! chispa en chorro
           write(*,'(A)') '   -> la chispa cae DENTRO del chorro de la ranura.'
        end if
     end if
  else                                                        ! fallo
     write(*,'(A)') '   -> ERROR: la chispa no toca ningun nodo de gas.'
     write(*,'(A)') '      Con esos limites no va a encender nunca.'
  end if
  write(*,'(A,F8.2,A,F8.2,A)') ' encendedor de superficie x   : ', &
       x_ign_0*1000.0_dp, ' a ', x_ign_1*1000.0_dp, ' mm'
  write(*,'(A)') '====================================================='
end program verificar_entrada