program combustion_pe
  use mod_precision
  use mod_termo
  use mod_cinetica
  use mod_malla
  use mod_solver
  implicit none

  integer  :: it, frame, i, n_gr                                  ! indices
  real(dp) :: dts, masa, masa0, Tmax, t_gas, Hp_medio             ! acumuladores

  call cargar_tablas()                                             ! tablas NASA [Ec.1-3]
  call init_cinetica()                                             ! mecanismo [Ec.20-24]
  call pretabular_k()                                              ! tabula kf,kb,k0 [Ec.20,23,24]
  call init_wilke()                                                ! constantes Wilke [Ec.13]
  call inicializar()                                               ! malla + estado [Ec.36-38]
  call escribir_malla_fija()                                       ! malla inicial

  frame  = 0                                                       ! contador volcados
  tiempo = 0.0_dp                                                  ! t solido [s]
  t_gas  = 0.0_dp                                                  ! t gas [s]
  masa0  = masa_combustible()                                      ! masa inicial [g]

  write(*,'(A,I0,A,I0,A,I0,A)') ' Malla: ', NN, ' nodos (', NX, ' x ', NY, ')'
  write(*,'(A,I0,A)')           ' Gas: ', NX*NYG, ' nodos ; solido: 1D por columna'
  write(*,'(A,ES10.3,A)')       ' Primer dy junto a la pared: ', &
       yg(1,NYS+1)-yg(1,NYS), ' m'
  write(*,'(A,F7.2,A,F8.2,A,F7.1,A)') ' Entrada: U=', U_ent, ' m/s   G_ox=', &
       G_ox, ' kg/(m2 s)   masa de combustible=', masa0, ' g'
  write(*,'(A,F6.2,A,F6.2,A,F7.4)') ' Puerto: Hp0=', Hp0*1000.0_dp, &
       ' mm -> Hp_max=', Hp_max*1000.0_dp, ' mm   Mach de entrada ~', U_ent/430.0_dp  ! [Ec.30] Ma
  write(*,'(A,F6.1,A,F6.1,A,F6.1,A)') ' Geometria: precamara ', L_precam*1000.0_dp, &
       ' mm | grano ', L_grano*1000.0_dp, ' mm | poscamara ', L_poscam*1000.0_dp, ' mm'
  write(*,'(A,F6.2,A,F6.2,A,I0,A,I0)') ' Liner inerte ', e_pared*1000.0_dp, &
       ' mm   rampas de ', L_RAMPA*1000.0_dp, ' mm   columnas con grano: ', &
       count(es_grano), ' de ', NX
  write(*,'(A,F7.2,A,F7.2)') ' Pendiente de la rampa: ', &
       1.5_dp*(e0-e_pared)/L_RAMPA, '   max|dY/dx| en la malla: ', maxval(abs(dYdx))  ! [Ec.38] M
  if (maxval(abs(dYdx)) > 5.0_dp) then                             ! aviso metrica
     write(*,'(A)') ' AVISO: la metrica dY/dx es grande. El termino u*dY/dx de la'
     write(*,'(A)') '        velocidad relativa a la malla va a recortar el paso de'
     write(*,'(A)') '        tiempo por Courant. Alargar L_RAMPA si se vuelve inviable.'
  end if
  if (U_ent/430.0_dp > 0.3_dp) then                                ! [Ec.30] Ma>0.3
     write(*,'(A)') ' AVISO: Mach > 0.3. La formulacion de bajo Mach deja de ser'
     write(*,'(A)') '        valida y el paso de tiempo se hara muy pequeno.'
     write(*,'(A)') '        Revisa Hp0, mdot_ox y G_ox: son incompatibles entre si.'
  end if
  write(*,'(A)') ' ---------------------------------------------------------------'
  flush(6)

  do it = 1, iter_max           !aqui podemos poner el limite                ! bucle principal

     call actualizar_entrada()                                     ! U_ent, G_ox segun modo
     call paso_tiempo()                                            ! [Ec.50-52]

     dts = dt                                                      ! dt solido
     if (tiempo >= t_acel_ini) dts = dt*FAC_SOLIDO                 ! desacople: FAC*dt_gas

     rho_o = rho                                                   ! guarda para [Ec.48]

     call fuente_radiacion()                                       ! [Ec.19]
     call transporte_escalares()                                   ! [Ec.34,35]
     call quimica()                                                ! [Ec.27-29]
     call recuperar_T()                                            ! [Ec.7]
     call contorno_escalares()                                     ! contornos Y,h,T
     call propiedades()                                            ! [Ec.9-14]

     call balance_superficie()                                     ! [Ec.40-42]
     call solido_1d(dts)                                           ! [Ec.43,44]
     call mover_superficie(dts)                                    ! [Ec.45]
     call momento()                                                ! [Ec.32,33]
     call poisson_correccion()                                     ! [Ec.48,49]

     tiempo = tiempo + dts                                         ! avanza t solido
     t_gas  = t_gas  + dt                                          ! avanza t gas

     call diagnostico(it, t_gas)                                   ! log iteracion
     if (mod(it,20) == 0) call verificar(it)                       ! chequeo NaN

     if (mod(it,500) == 0) then                                    ! latido progreso
        write(*,'(A,I8,A,ES10.3,A,ES10.3,A,ES9.2,A,F7.1,A,F7.1,A,I6,A,F8.1,A,I3,A)') &
             '   it=', it, '  t=', tiempo, ' s (gas ', t_gas, ')  dt=', dt, &
             '  Tmax=', maxval(Tg(:,NYS:NY)), ' K  Tw=', maxval(Tw), &
             ' K  celdas_q=', n_celdas_q, '  |U|max=', diag_umax, &
             '  kP=', diag_k_poisson, merge('  ', ' *', diag_poisson_ok)
        flush(6)
     end if

     if (mod(it,n_salida) == 0) then                               ! volcado periodico
        frame = frame + 1
        call escribir_salidas(frame)
        masa = masa_combustible()                                  ! masa actual [g]
        Tmax = maxval(Tg(:,NYS:NY))                                ! T max gas
        write(*,'(A,I8,A,F8.3,A,F7.1,A,F6.3,A,F6.3,A,F6.2,A,F5.1,A,F5.2,A)') &
             ' it=',it,' t=',tiempo,'s Tmax=',Tmax,'K Tw=', &
             maxval(Tw)/1000.0_dp,'kK r=',maxval(rdot)*1.0e3_dp, &
             'mm/s masa=',masa,'g (',100.0_dp*masa/masa0,'%) Gox=',G_ox, &
             ' kg/m2s'
        write(*,'(A,F7.4,A,F7.4,A,F5.1,A)') '          rho: min=', &
             minval(rho(:,NYS:NY)),'  max=',maxval(rho(:,NYS:NY)), &
             ' kg/m3  (relacion ',maxval(rho(:,NYS:NY))/ &
             minval(rho(:,NYS:NY)),'x)'
        write(*,'(A,3(A,ES10.3))') '          radicales:', &
             '  H=', maxval(Ysp(:,NYS:NY,I_H)), &
             '  O=', maxval(Ysp(:,NYS:NY,I_O)), &
             '  OH=',maxval(Ysp(:,NYS:NY,I_OH))
        write(*,'(A,I8,A,I5)') '          quimica: celdas=', n_celdas_q, &
             '  subpasos max=', n_sub_max
        flush(6)
     end if

     if (tiempo >= t_final) exit                                   ! t final

     Hp_medio = 0.0_dp;  n_gr = 0                                  ! puerto medio sobre grano
     do i = 1, NX
        if (.not. es_grano(i)) cycle
        Hp_medio = Hp_medio + (Ly - esp(i))                        ! Hp(i) = Ly - e
        n_gr     = n_gr + 1
     end do
     if (n_gr > 0) Hp_medio = Hp_medio/real(n_gr, dp)

     if (Hp_medio >= Hp_max) then                                  ! puerto agotado
        write(*,'(A,F7.3,A,F6.2,A)') ' Puerto sobre el grano abierto hasta '// &
             'Hp_max en t = ', tiempo, ' s (', Hp_max*1000.0_dp, ' mm)'
        exit
     end if

     if (all(agotado .or. .not. es_grano)) then                    ! combustible agotado
        write(*,'(A,F8.3,A)') ' Combustible agotado en t = ', tiempo, ' s'
        exit
     end if
  end do

  call escribir_salidas(frame+1)                                   ! volcado final
  write(*,'(A)') ' Simulacion terminada.'

end program combustion_pe