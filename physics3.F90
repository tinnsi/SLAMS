SUBROUTINE density(p)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION :: rho, porosity, tvol
   DOUBLE PRECISION :: orgC, calc, opal, clay, dens, TEPC
   DOUBLE PRECISION :: rho_solid
   INTEGER :: i, j 

   CALL fractal_dimension(p)
   CALL particle_volume(p, tvol)

! moles of solid stuff
   i = 1
   orgC = 0
   TEPC = 0
   DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   ENDDO
   calc = p%mineral(1) + p%mineral(2)
   opal = p%mineral(3)
   clay = p%mineral(4)
  
   rho_solid = (orgC*mw_orgC + &
               calc*mw_co + &
               opal*mw_si + &
               clay*mw_li + &
               TEPC*mw_TEP) / tvol

   IF (p%af .eq. 0 ) THEN
      porosity = 1.0 - p%Nn*(p%pr/p%r)**3
   ELSEIF (p%af .eq. 1 ) THEN
      porosity = 0.5
   ENDIF

   IF (porosity .lt. 0 .or. porosity .gt. 1) THEN 
   print*, 'porosity less than zero!', porosity, p%frac, p%s, p%Nn, &
           p%Nn*(p%pr/p%r)**3, TEPC, tvol
      porosity = 0.1
   ENDIF

   IF (p%Nn .eq. 1) THEN
      dens = rho_solid
      porosity = 0
   ELSEIF (p%Nn .gt. 1) THEN
      dens = rho_solid*(1.0-porosity) + rho_sw*porosity
   ENDIF
   
   p%por = porosity
   p%rho = dens

   IF (p%por .lt. 0 .or. p%por .gt. 1) THEN
      IF (p%Nn .gt. 1) THEN
         print*, 'r lt pr', p%Nn, p%r, tvol
         print*, rho_solid
         print*, p
      ENDIF
      p%r = p%pr
   ENDIF

   IF (p%r .eq. 0 .or. p%rho .le. 0) THEN
      print*, 'in density', p%id, p%n, p%Nn
      print*, p%rho, p%r, porosity
      print*, p%orgC(1,1), p%orgC(2,1), p%orgC(3,1)
      print*, vorgC, orgC, tvol
   ENDIF

   IF (dens .le. 0) THEN
      print*, 'rho - problem', porosity, dens, tvol
   ENDIF
 
   IF (p%Nn .gt. 1 .and. p%r .eq. p%pr) THEN
      print*, 'NOT GOOD', p%id, p%Nn
   ENDIF

   IF (TEP .lt. 0) print*, 'ERROR: density(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE density
!=====================================================================
SUBROUTINE fractal_dimension(p)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p

      p%frac = 2.0

END SUBROUTINE fractal_dimension
!=====================================================================
SUBROUTINE dissolve(p)
   USE the_info
! the aggregate is so tiny that it is considered dissolved

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION :: orgC, calc, arag, opal, clay, dens, TEPC
   INTEGER :: i, j

   CALL find_depth(p%z,j)

   IF (p%r .lt. 0.2 .and. p%b .eq. 0) THEN
      i = 1
      orgC = 0
      TEPC = 0
      DO WHILE (i .le. 10) 
         orgC = orgC + p%orgC(i,1) 
         TEPC = TEPC + p%TEP(i,1) 
         i = i + 1
      ENDDO
      calc = p%mineral(1) 
      arag = p%mineral(2)
      opal = p%mineral(3)
      clay = p%mineral(4)

      p%orgC = 0
      p%mineral = 0
      p%TEP = 0
  
      p%b = 1
      inventory(1) = inventory(1) - orgC*p%n - TEPC*p%n
      inventory(2) = inventory(2) - (arag+calc)*p%n
      inventory(3) = inventory(3) - opal*p%n
      inventory(4) = inventory(4) - clay*p%n

      lost(1) = lost(1) + orgC*p%n 
      lost(2) = lost(2) + (arag+calc)*p%n
      lost(3) = lost(3) + opal*p%n
      lost(4) = lost(4) + clay*p%n
      lost(5) = lost(5) + TEPC*p%n

      orgC_b(j) = orgC_b(j) + orgC*p%n + TEPC*p%n
      calcC_t(j) = calcC_t(j) + calc*p%n
      calcA_t(j) = calcA_t(j) + arag*p%n 
      opal_t(j) = opal_t(j) + opal*p%n 
   ENDIF

   i = 1
   DO WHILE (i .le. 10)
      IF (p%orgC(i,1) .lt. 1e-50) THEN
         p%orgC(i,1) = 0.0
      ENDIF
      i = i + 1
   ENDDO

   IF (TEPC .lt. 0) print*, 'ERROR: dissolve(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE dissolve
!====================================================================
SUBROUTINE disintegrate(p, harvest)
! the particle falls apart if it is not sticky enough to be 
! bound together 
! or, if there is too much shear acting on that particle
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest
   DOUBLE PRECISION :: new_Nn, old_Nn

   IF (p%s .lt. 0.02 .and. p%Nn .gt. 1) THEN
      CALL random_number(harvest)
      IF (harvest .gt. 0.99) THEN
         CALL break_p(p, harvest)
      ENDIF
   ENDIF
   IF (p%r .gt. 1e4) THEN
      CALL random_number(harvest)
      IF (timestep/86400.0*p%r/1.0e6 .gt. harvest) THEN
         CALL break_p(p,harvest)
      ENDIF
   ENDIF

END SUBROUTINE disintegrate
!====================================================================
SUBROUTINE particle_volume(p, tvol)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: tvol
   DOUBLE PRECISION :: ar, avol, orgC, TEPC, calc, opal, clay
   DOUBLE PRECISION :: vorgC, vTEPC, vcalc, vopal, vclay, pvol

! moles of solid stuff
   i = 1
   orgC = 0
   TEPC = 0
   DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   ENDDO
   calc = p%mineral(1) + p%mineral(2)
   opal = p%mineral(3)
   clay = p%mineral(4)
  
! change moles to volume
   vorgC = orgC * mw_orgC / rho_orgC
   vTEPC = TEPC * mw_TEP / rho_TEP
   vcalc = calc * mw_co / rho_co
   vopal = opal * mw_si / rho_si
   vclay = clay * mw_li / rho_li

   tvol = vorgC + vcalc + vopal + vclay + vTEPC

   IF (TEPC .lt. 0) print*, 'ERROR: partVol(TEP)', TEPC, p%TEP(1,1), testVar, p%id

END SUBROUTINE particle_volume
!====================================================================
SUBROUTINE particle_mass(p, mass)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: mass
   DOUBLE PRECISION :: orgC, TEPC, calc, opal, clay

! units of mass [g]
! moles of solid stuff
   i = 1
   orgC = 0
   TEPC = 0
   DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   ENDDO
   calc = p%mineral(1) + p%mineral(2)
   opal = p%mineral(3)
   clay = p%mineral(4)

   mass = orgC * mw_orgC + TEPC * mw_TEP + calc * mw_co + &
          opal * mw_si + clay * mw_li
   
   IF (TEPC .lt. 0) print*, 'ERROR: partMass(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE particle_mass
!====================================================================
SUBROUTINE stickiness(p, MLdepth, testVar)
! stickiness is a function of volume_TEP over vol_solid
   USE the_info
  
   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: MLdepth, testVar
   DOUBLE PRECISION :: tvol, TEPC, vTEPC, orgC, vorgC
   INTEGER :: i

   CALL particle_volume(p, tvol)

   orgC = 0
   TEPC = 0
   i = 1
   DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      i = i + 1
   ENDDO
   i = 1
   DO WHILE (i .le. 10) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   ENDDO
   vTEPC = TEPC * mw_TEP / rho_TEP
   vorgC = orgC * mw_orgC / rho_orgC

   p%s = vTEPC / tvol !* testVar

   IF (p%b .eq. 0) THEN
      IF (p%s .lt. 0 .or. TEPC .lt. 0) THEN
         print*, 'stickiness II', p
      ENDIF
   ENDIF

   IF (TEPC .lt. 0) print*, 'ERROR: stickiness(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE stickiness
!====================================================================
SUBROUTINE velocity(p, paraT, MLdepth)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(IN) :: paraT, MLdepth
   DOUBLE PRECISION :: d_rho, one, u, fu, du
   DOUBLE PRECISION :: new_u, d, nu, pm_sign
   INTEGER :: i

   CALL viscosity(p%z, nu, paraT, MLdepth)
   one = 1
   i = 1
   d_rho = ABS(p%rho - rho_sw)
   IF (p%rho .gt. rho_sw) THEN
      pm_sign = 1
   ELSEIF (p%rho .lt. rho_sw) THEN
      pm_sign = -1
   ENDIF
   d = 2*p%r/10000   !now d is in cm
   u = g*d_rho*d**2/(18*rho_sw*nu)
   DO WHILE (one .gt. 0.001 .and. i .lt. 10)
      i = i + 1
      fu = 24*nu*u/d + 6*u**1.5/(1+(d/nu)**0.5) + 0.4*u**2 &
      - 4*g*d_rho*d/(3*rho_sw)
      du = 24*nu/d + 9*u**0.5/(1+(d/nu)**0.5) + 0.8*u
      new_u = u - fu/du
      IF (new_u .eq. u) THEN 
         print*, 'problem in velocity', d, p%id
         print*, p
      ENDIF
      one = ABS(new_u - u)
      u = new_u
   END DO
   p%w = pm_sign * u ! cm/s

   IF (p%w .gt. 1e4) print*, 'velocity troubles', p%w, d_rho
   IF (p%w .lt. -1e4) print*, 'velocity oops', p%w, d_rho

END SUBROUTINE velocity
!====================================================================
SUBROUTINE check_depth(pa, pz, depth, i, j, z_max)

   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: pa
   DOUBLE PRECISION, INTENT(IN) :: depth
   INTEGER, INTENT(INOUT) :: pz, j, z_max
   INTEGER, INTENT(IN) :: i
   DOUBLE PRECISION :: d

   d = depth + dz/2.0
   IF (ABS(pa%z - d) .le. dz/2.0) THEN
      pz = i 
      z_max = j
      j = j + 1
   ENDIF
    
END SUBROUTINE check_depth
!=====================================================================
SUBROUTINE find_depth(z, i)
! z is depth in cm, i is the depth-box
   USE the_info

   DOUBLE PRECISION, INTENT(IN) :: z
   INTEGER, INTENT(INOUT) :: i
   INTEGER :: aha

   i = 1
   aha = 0
   DO WHILE (aha .eq. 0) 
      IF (z .ge. (i-1)*dz .and. z .lt. (i)*dz) THEN
         aha = 1
      ELSE 
         i = i + 1
      ENDIF
   ENDDO

END SUBROUTINE find_depth
!=====================================================================
SUBROUTINE find_size(r, i)
! r is depth in um, i is the radius box
   USE the_info

   INTEGER, INTENT(INOUT) :: i
   DOUBLE PRECISION, INTENT(IN) :: r
   INTEGER :: aha

   i = 1
   aha = 0
   DO WHILE (i .lt. n_size-1 .and. aha .eq. 0)
      IF (r .lt. radi(i)) THEN
         aha = 1
      ELSE
         i = i + 1
      ENDIF
   ENDDO

   IF (i .gt. 60) print*, 'find_size:', i, r

END SUBROUTINE find_size
!=====================================================================
SUBROUTINE find_velo(w, i)
! r is depth in um, i is the radius box
   USE the_info

   INTEGER, INTENT(INOUT) :: i 
   DOUBLE PRECISION, INTENT(IN) :: w
   INTEGER :: aha

   i = 1
   aha = 0
   DO WHILE (i .lt. n_size .and. aha .eq. 0)
      IF (w*864 .lt. radi(i)) THEN
         aha = 1
      ELSE
         i = i + 1
      ENDIF
   ENDDO

   IF (i .gt. 60) print*, 'find_velo:', i, w

END SUBROUTINE find_velo
!=====================================================================
SUBROUTINE find_dens(p, i)
! r is depth in um, i is the radius box
   USE the_info

   INTEGER, INTENT(INOUT) :: i 
   TYPE(agg_part), INTENT(INOUT) :: p  
   DOUBLE PRECISION :: density, orgC, TEPC, calc, opal, clay
   DOUBLE PRECISION :: rho_solid, tvol
   INTEGER :: aha

   i = 1
   aha = 0
   density = 0.8

   CALL particle_volume(p, tvol)
! moles of solid stuff
   i = 1
   orgC = 0
   TEPC = 0
   DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   ENDDO
   calc = p%mineral(1) + p%mineral(2)
   opal = p%mineral(3)
   clay = p%mineral(4)
  
   rho_solid = (orgC*mw_orgC + &
               calc*mw_co + &
               opal*mw_si + &
               clay*mw_li + &
               TEPC*mw_TEP) / tvol

   i = 1
   DO WHILE (i .le. n_size .and. aha .eq. 0)
      IF (rho_solid .lt. density) THEN
         aha = 1
      ELSE
         density = density + 0.04
         i = i + 1
      ENDIF
   ENDDO

   IF (i .gt. n_size) THEN
      print*, 'find_dens:', p%r, p%rho
      i = n_size
   ENDIF

   IF (TEPC .lt. 0) print*, 'ERROR: find_dens(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE find_dens
!=====================================================================
SUBROUTINE bottom(p, j)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   INTEGER, INTENT(IN) :: j
   DOUBLE PRECISION :: organic, TEPC
   INTEGER :: i, l, m

   IF (p%z .gt. seafloor) THEN
      p%b = 1
      i = 1
      organic = 0
      TEPC = 0
      DO WHILE (i .le. 10)
         organic = organic + p%orgC(i,1)*p%n
         TEPC = TEPC + p%TEP(i,1)*p%n
         seaBed(i) = seaBed(i) + p%orgC(i,1)*p%n + p%TEP(i,1)
         i = i + 1
      ENDDO
      sf(1) = sf(1) + organic 
      sf(2) = sf(2) + (p%mineral(1)+p%mineral(2))*p%n
      sf(3) = sf(3) + p%mineral(3)*p%n
      sf(4) = sf(4) + p%mineral(4)*p%n
      sf(5) = sf(5) + TEPC
      seaBed(11) = seaBed(11) + p%mineral(1)*p%n
      seaBed(12) = seaBed(12) + p%mineral(2)*p%n
      seaBed(13) = seaBed(13) + p%mineral(3)*p%n
      seaBed(14) = seaBed(14) + p%mineral(4)*p%n

      inorganic = p%mineral(1)*p%n + p%mineral(2)*p%n

      IF (p%af .gt. 0) THEN
         seaBed(15) = seaBed(15) + organic
         seaBed(16) = seaBed(16) + inorganic
         seaBed(17) = seaBed(17) + p%mineral(3)*p%n
         seaBed(18) = seaBed(18) + p%mineral(4)*p%n
      ENDIF
 
      sizeBed(j) = sizeBed(j) + organic

      CALL find_velo(p%w, l)
      veloBed(l) = veloBed(l) + organic

      m = 0
      CALL find_dens(p, m)
      densBed(m) = densBed(m) + organic
   ENDIF

   IF (TEPC .lt. 0) print*, 'ERROR: bottom(TEP)', TEP, p%TEP(1,1), testVar, p%id

END SUBROUTINE bottom
!=====================================================================
SUBROUTINE pool(p1, p2)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p1, p2
   INTEGER :: i

   i = 1
   DO WHILE (i .le. 10) 
      p1%orgC(i,1) = (p1%n*p1%orgC(i,1)+p2%n*p2%orgC(i,1))/(p1%n+p2%n)
      p1%TEP(i,1) = (p1%n*p1%TEP(i,1)+p2%n*p2%TEP(i,1))/(p1%n+p2%n)
      i = i + 1
   ENDDO

   i = 1
   DO WHILE (i .le. 4) 
      p1%mineral(i) = (p1%n*p1%mineral(i)+p2%n*p2%mineral(i))/(p1%n+p2%n)
      i = i + 1
   ENDDO

   p1%n = p1%n + p2%n
   p1%p = p1%p + p2%p
   p1%z = (p1%z*p1%n + p2%z*p2%n)/(p1%n+p2%n)

   p2%b = 1

END SUBROUTINE pool
!=====================================================================
SUBROUTINE radius(p)
   USE the_info
   ! this subroutine is called by stick() and break_p and respiration
   ! calculates new radius after coagulation or fragmentation

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION :: tvol, pvol, vol_fp

   CALL fractal_dimension(p)
   CALL particle_volume(p, tvol)

   pvol = tvol/p%Nn

   p%pr = ((3.0*pvol/(4.0*pi))**(1.0/3.0)*1e4)
   
   IF (p%af .eq. 0) THEN !if particle is an aggregate:
      p%r = p%pr*p%Nn**(1.0/p%frac)
   ELSEIF(p%af .eq. 1) THEN !particle is fecal pellet:
      vol_fp = tvol * 2.0
      p%r = ((3.0*vol_fp/(4.0*pi))**(1.0/3.0)*1e4)
   ELSE 
      print*, p%af, 'radius!!'
      p%r = p%pr*p%Nn**(1.0/p%frac)
   ENDIF

END SUBROUTINE radius
!=====================================================================
SUBROUTINE collide4(p1, p2, harvest, paraT, MLdepth, testVar)
   
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p1, p2
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, paraT, MLdepth, testVar
   DOUBLE PRECISION :: test1, test2, endofdt, semiTime
   DOUBLE PRECISION :: q, dt, p, muu, nu, epsi
   DOUBLE PRECISION :: beta, betaBr, betaSh, betaDS
   DOUBLE PRECISION :: r1, r2, w1, w2, prob, alpha
   DOUBLE PRECISION :: organic, inorganic, dwTotal, s, dw_orgC
   DOUBLE PRECISION :: r, mean, meanp
   INTEGER :: q_int

   r1 = p1%r/1e6  ! r1 in m
   r2 = p2%r/1e6  ! r2 in m
   w1 = p1%w/100  ! w1 in m/s
   w2 = p2%w/100  ! w1 in m/s
 
   p = MIN(r1,r2)/MAX(r1,r2)
   
   epsi = 1e-4 !turbulent energy dissipation rate [cm2/s3]
   CALL viscosity(p1%z, muu, paraT, MLdepth)  !muu dynamic viscosity
   nu = muu /rho_sw                  !nu kinematic viscisity [cm2/s]

test1 = p1%p + p2%p

   endofdt = timestep
   dt = timestep
   semiTime = 0
   test1 = p1%p + p2%p
   DO WHILE (semiTime .lt. endofdt .and. MAX(p1%Nn,p2%Nn) .lt. 1e9)
      alpha = MAX(p1%s,p2%s)
!Brownian motion:
      betaBr = 8*k*T/(6*mu)*(r1+r2)**2/(r1*r2)   ! in m^3/s
!Shear:
!      IF (p1%z .le. 10000) THEN
!         betaSh = 4/3.0*gamma*(r1+r2)**3              ! in m^3/s Rectil.
!      ELSEIF (p1%z .gt. 10000) THEN
!         betaSh = 4/3.0*(gamma/100)*(r1+r2)**3          ! in m^3/s
!      ENDIF
      betaSh = 9.8*(p**2/(1+2*p**2))*(epsi/nu)**0.5*(r1+r2)**3  ! in m^3/s C
!Differential settling
!      betaDS = pi*(r1+r2)**2*ABS(w2-w1)          ! rectilinear
      betaDS = 0.5*pi*MIN(r1,r2)**2*ABS(w2-w1)          ! in m^3/s

      beta = betaBr + betaSh + betaDS                   ! in m^3/s

      prob = alpha*beta*dt
      DO WHILE (prob .gt. 1.0) 
         dt = dt/10.0
         prob = alpha*beta*dt
      ENDDO

! find number of j-particles to collide with i-particles
      q = MAX(p1%n,p2%n)/MIN(p1%n,p2%n)
      mean = q * (1 - (1-prob)**MIN(p1%n,p2%n))
      meanp = mean - FLOOR(mean)
      CALL random_number(harvest)
      IF (meanp .gt. harvest) THEN
         q_int = CEILING(mean)
      ELSEIF(meanp .lt. harvest) THEN
         q_int = FLOOR(mean)
      ENDIF

      IF (q_int .lt. 0) THEN
         print*, 'ERROR: collide: q_int lt 0', q, mean, q_int, sigma
         q_int = 0
      ENDIF

      IF (q_int .gt. 0) THEN
         CALL stick2(p1,p2,q_int, MLdepth, testVar)
      ENDIF

      semiTime = semiTime + dt
   ENDDO

   test2 = p1%p + p2%p

   IF (ABS(test1-test2) .gt. 10) THEN
      print*, 'ERROR: collide4', test1, test2, p1%id, p2%id, p1%Nn, p2%Nn, p1%n, p2%n
   ENDIF

END SUBROUTINE collide4
!========================================================================
SUBROUTINE stick2(p1, p2, intRate, MLdepth, testVar)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p1, p2
   INTEGER, INTENT(INOUT) :: intRate
   DOUBLE PRECISION, INTENT(INOUT) :: MLdepth, testVar
   INTEGER :: i, whereAmI
   DOUBLE PRECISION :: test1, test2, x, bak_n1, bak_n2

   test1 = p1%orgC(1,1)*p1%n + p2%orgC(1,1)*p2%n
   x = 2.0  !this is to get out of problems with TYPE and MOD


   bak_n1 = p1%n
   bak_n2 = p2%n
   IF (p1%n .gt. p2%n .and. intRate*p2%n .lt. p1%n) THEN
      whereAmI = 1
      p2%Nn = p1%Nn*intRate + p2%Nn
      p2%af = 0
      i = 1
      DO WHILE (i .le. 10)
         p2%orgC(i,1) = p2%orgC(i,1) + intRate*p1%orgC(i,1)
         p2%TEP(i,1) = p2%TEP(i,1) + intRate*p1%TEP(i,1)
         i = i + 1
      ENDDO
      i = 1
      DO WHILE (i .le. 4)
          p2%mineral(i) = p2%mineral(i) + intRate*p1%mineral(i)
          i = i + 1
      ENDDO
      p1%n = p1%n - intRate*p2%n
      CALL radius(p2)
      CALL density(p2)
      CALL stickiness(p2, MLdepth, testVar)
   ELSEIF (p1%n .lt. p2%n .and. intRate*p1%n .lt. p2%n) THEN
      whereAmI = 2
      p1%Nn = p1%Nn + p2%Nn*intRate
      p1%af = 0
      i = 1
      DO WHILE (i .le. 10)
         p1%orgC(i,1) = p1%orgC(i,1) + intRate*p2%orgC(i,1)
         p1%TEP(i,1) = p1%TEP(i,1) + intRate*p2%TEP(i,1)
         i = i + 1
      ENDDO
      i = 1
      DO WHILE (i .le. 4)
         p1%mineral(i) = p1%mineral(i) + intRate*p2%mineral(i)
         i = i + 1
      ENDDO
      p2%n = p2%n - intRate*p1%n
      CALL radius(p1)
      CALL density(p1)
      CALL stickiness(p1, MLdepth, testVar)
   ELSEIF (p1%n .eq. p2%n) THEN               
      whereAmI = 3
      p1%Nn = p1%Nn + p2%Nn
      p2%Nn = p1%Nn
      p1%af = 0
      p2%af = 0
      i = 1
      DO WHILE (i .le. 10)
         p1%orgC(i,1) = p1%orgC(i,1) + p2%orgC(i,1) 
         p2%orgC(i,1) = p1%orgC(i,1)
         p1%TEP(i,1) = p1%TEP(i,1) + p2%TEP(i,1) 
         p2%TEP(i,1) = p1%TEP(i,1)
         i = i + 1
      ENDDO
      i = 1
      DO WHILE (i .le. 4)
         p1%mineral(i) = p1%mineral(i) + p2%mineral(i)  
         p2%mineral(i) = p1%mineral(i)
         i = i + 1
      ENDDO
      IF (MOD (p1%n,x) .eq. 0.0) THEN
         p1%n = p1%n/2.0
         p2%n = p1%n
      ELSE
         p1%n = (p1%n+1)/2.0
         p2%n = (p2%n-1)/2.0
      ENDIF
      CALL radius(p1)
      CALL radius(p2)
      CALL density(p1)
      CALL density(p2)
      CALL stickiness(p1, MLdepth, testVar)
      CALL stickiness(p2, MLdepth, testVar)
   ENDIF
   p1%p = p1%Nn * p1%n
   p2%p = p2%Nn * p2%n

   test2 = p1%orgC(1,1)*p1%n + p2%orgC(1,1)*p2%n

   IF (ABS(test1-test2) .gt. 1)  print*, 'stick2: failed conservation of mass', test1, test2, rate

   IF (p1%n .lt. 0) print*, 'p1 less than 0', p1%Nn, p2%Nn, WhereAmI, &
      bak_n1, bak_n2, p1%n, p2%n, intRate, p1%id, p2%id
   IF (p2%n .lt. 0) print*, 'p2 less than 0', p1%Nn, p2%Nn, WhereAmI, &
      bak_n1, bak_n2, p1%n, p2%n, intRate, p1%id, p2%id
   IF (p1%orgC(1,1) .lt. 0 .or. p2%orgC(1,1) .lt. 0) print*, 'orgC!!!', &
      WhereAmI, p1%orgC(1,1), p1%Nn, p1%n, p2%Nn, p2%n, intRate
   IF (p1%Nn .lt. 0) print*, 'N1 lt 0', WhereAmI, rate, bak_n1, bak_n2
   IF (p2%Nn .lt. 0) print*, 'N2 lt 0', WhereAmI, rate, bak_n1, bak_n2
!   IF (rate .gt. 1000) print*, WhereAmI, rate, p1%n, p2%n, p1%Nn, p2%Nn, 'o.o'

END SUBROUTINE stick2
!========================================================================
SUBROUTINE self_collide(p, harvest, MLdepth, testVar)
   
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, MLdepth, testVar
   DOUBLE PRECISION :: beta, betaBr, betaSh, betaDS
   DOUBLE PRECISION :: r1, w1,prob, alpha
   DOUBLE PRECISION :: test1, test2

   test1 = p%p

   r1 = p%r/1e6  ! r1 in m

   IF (p%n .gt. 2) THEN
      alpha = p%s
      betaBr = 8*k*T/(6*mu)*(2*r1)**2/(r1**2)   ! in m^3/s
      IF (p%z .le. 10000) THEN
         betaSh = 4/3.0*gamma*(r1+r2)**3              ! in m^3/s
      ELSEIF (p%z .gt. 10000) THEN
         betaSh = 4/3.0*(gamma/100)*(r1+r2)**3          ! in m^3/s
      ENDIF
      betaDS = 0.0
      beta = (betaBr + betaSh + betaDS)
!      prob = 1 - exp(-alpha*beta*sf1*timestep)
      prob = alpha*beta*p%n*timestep
      CALL random_number(harvest)
      
      IF (prob .gt. harvest) THEN
         CALL self_stick(p, MLdepth, testVar)
      ENDIF
   ENDIF

   IF (prob .gt. 1) THEN
      print*, 'self_collide', p%id, p%Nn, p%n
   ENDIF
 
   test2 = p%p
   IF (ABS(test1-test2) .gt. 1) THEN 
      print*, 'self', p%id, p%Nn
   ENDIF

END SUBROUTINE self_collide
!========================================================================
SUBROUTINE self_stick(p, MLdepth, testVar)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: MLdepth, testVar
   INTEGER :: i
   DOUBLE PRECISION :: test1, test2

   test1 = p%orgC(1,1)*p%n
   p%n = p%n / 2.0
   p%Nn = p%Nn * 2.0
   p%af = 0
   
   i = 1
   DO WHILE (i .le. 4)
      p%mineral(i) = p%mineral(i) * 2.0
      i = i + 1
   ENDDO

   i = 1 
   DO WHILE (i .le. 10) 
      p%orgC(i,1) = p%orgC(i,1) * 2.0
      p%TEP(i,1) = p%TEP(i,1) * 2.0
      i = i + 1
   ENDDO

   CALL radius(p)
   CALL density(p)
   CALL stickiness(p, MLdepth, testVar)

   test2 = p%orgC(1,1)*p%n

   IF (ABS(test1-test2) .gt. 1) THEN 
      print*, 'self_stick', p%Nn, p%n
   ENDIF
print*, 'SELF-STICK', p%id, p%TEP(1,1)

END SUBROUTINE self_stick
!========================================================================
SUBROUTINE break_p(p, harvest)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest
   DOUBLE PRECISION :: x, new, old, daughters, probability
   DOUBLE PRECISION :: Nn, new_Nn, new_n, extra, add, old_Nn, old_p, new_p
   DOUBLE PRECISION :: before, after, TEP_t
   INTEGER :: i

   daughters = 2
   CALL random_number(harvest)
   probability = 0.91*daughters**-1.56
   DO WHILE(probability .lt. harvest)
      daughters = daughters + 1
      probability = probability + 0.91*daughters**-1.56
   ENDDO
   TEP_t = 0
   i = 1
   DO WHILE (i .le. 10)
      TEP_t = TEP_t + p%TEP(i,1)*p%n
      i = i + 1
   ENDDO
   before = TEP_t
   IF (p%Nn .gt. daughters) THEN
      Nn = p%Nn
      old_Nn = p%Nn
      extra = MOD(Nn,daughters)
      Nn = p%Nn - extra

      new_Nn = Nn/daughters
      add = extra*p%n/new_Nn
      new_n = p%n*daughters + add
      new_p = new_Nn * new_n

      p%Nn = new_Nn
      p%n = new_n
      p%p = new_p
      p%af = 0
      
      i = 1
      DO WHILE (i .le. 10) 
         p%orgC(i,1) = p%orgC(i,1) * new_Nn/old_Nn
         p%TEP(i,1) = p%TEP(i,1) * new_Nn/old_Nn
         i = i + 1
      ENDDO
      i = 1
      DO WHILE (i .le. 4) 
         p%mineral(i) = p%mineral(i) * new_Nn/old_Nn
         i = i + 1
      ENDDO
      CALL radius(p)
   ENDIF 
   TEP_t = 0
   i = 1
   DO WHILE (i .le. 10)
      TEP_t = TEP_t + p%TEP(i,1)*p%n
      i = i + 1
   ENDDO
   after = TEP_t

   IF (ABS(before-after) .gt. before) print*, 'TEP - break', before, after

END SUBROUTINE break_p
!========================================================================
SUBROUTINE fecal_pellet(p, harvest)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest
   DOUBLE PRECISION :: test1, test2, oldN
   INTEGER :: i

   oldN = p%Nn
   test1 = p%orgC(1,1)*p%n
   CALL random_number(harvest)
   IF (p%af .eq. 1) THEN
      p%Nn = FLOOR(1000*harvest)
      IF (p%Nn .eq. 0) p%Nn = 100
   ENDIF

   p%n = p%p/p%Nn
   i = 1
   DO WHILE (i .le. 10) 
      p%orgC(i,1) = p%orgC(i,1)*p%Nn/oldN
      p%TEP(i,1) = p%TEP(i,1)*p%Nn/oldN
      i = i + 1
   ENDDO
   i = 1
   DO WHILE (i .le. 4) 
      p%mineral(i) = p%mineral(i)*p%Nn/oldN
      i = i + 1
   ENDDO
   test2 = p%orgC(1,1)*p%n
   CALL radius(p)
   CALL density(p)

   IF (ABS(test1-test2) .gt. 1) THEN
      print*, 'fecal pellet', test1, test2, '!'
      print*, p%Nn, p%n, oldN
   ENDIF
   IF (p%Nn .lt. 1 .or. p%n .lt. 0) THEN
      print*, 'fecal pellet' , p%Nn, p%n
   ENDIF

END SUBROUTINE fecal_pellet
!========================================================================
SUBROUTINE sink(p)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION :: z
   INTEGER :: i

   z = p%z + p%w*timestep

   i = 1
   DO WHILE (i .le. 9) 
      IF (p%z .lt. i*1e3 .and. z .gt. i*1e3) THEN 
         CALL sediment_trap(p,i)    ![i=1..10]
      ENDIF
      i = i + 1
   ENDDO
   DO WHILE (i .le. 20) 
      j = i - 9
      IF (p%z .lt. j*1e4 .and. z .gt. j*1e4) THEN 
         CALL sediment_trap(p,i)    ![i=11..20]
      ENDIF
      i = i + 1
   ENDDO
   j = 3
   DO WHILE (j .le. 8)
      IF (p%z .lt. j*5e4 .and. z .gt. j*5e4) THEN 
         CALL sediment_trap(p,i)     ![i=21..26]
      ENDIF
      i = i + 1
      j = j + 1
   ENDDO   

   p%z = p%z + p%w*timestep
   IF (p%z .lt. 0) p%z = 100

END SUBROUTINE sink
!========================================================================
SUBROUTINE respiration(p, bacZoo, lostOrg, paraT, MLdepth, testVar) 
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   INTEGER, INTENT(IN) :: bacZoo  
!  bacZoo = { 1(=bac), 2(=zoo) }
   DOUBLE PRECISION, INTENT(IN) :: paraT, MLdepth
   DOUBLE PRECISION, INTENT(INOUT) :: testVar
   DOUBLE PRECISION, INTENT(OUT) :: lostOrg
   DOUBLE PRECISION, DIMENSION(10) :: coef, coefTEP
   DOUBLE PRECISION :: temp, lostO, organic, temp_effect, TEPC, lostTEP
   INTEGER :: i, x

   CALL temperature(temp, p%z, paraT, MLdepth)
   CALL find_depth(p%z, x)
   
! total amount of organic carbon and TEP
   organic = 0
   TEPC = 0
   i = 1
   DO WHILE (i .le. 10)
      organic = organic + p%orgC(i,1)
      TEPC = TEPC + p%TEP(i,1)
      i = i + 1
   ENDDO
  
! bacZoo coefficients for bacterial (1) and zooplankton (2) respiration
! of "regular" organic carbon and TEP
   i = 1
   temp_effect = exp(DLOG(2.0)*(temp-30.0)/10.0)
   DO WHILE (i .le. 10)
      IF (bacZoo .eq. 1) THEN
         coef(i) = 0.1/p%orgC(i,2) * temp_effect !* testVar
         coefTEP(i) = 0.1/p%TEP(i,2) * temp_effect 
      ELSEIF (bacZoo .eq. 2) THEN
         IF (p%af .eq. 1) THEN      !Micro
            coef(i) = (1-(2.0**i/2.0**11))/timestep*0.9*temp_effect 
            coefTEP(i) = (1-(2.0**i/2.0**11))/timestep*0.9*temp_effect
         ENDIF
         IF (coef(i)*timestep .gt. 1) THEN
print*, 'SOS respiration', temp_effect, coef(i), i
            coef(i) = 1.0/timestep
         ENDIF
         IF (coefTEP(i)*timestep .gt. 1) THEN
print*, 'SOS - TEP', temp_effect, coef(i), i
            coefTEP(i) = 1.0/timestep
         ENDIF
      ENDIF
      i = i + 1
   ENDDO

   lostO = 0d00
   lostOrg = 0d00

IF (organic .gt. 0) THEN
   i = 1
   DO WHILE (i .le. 10) 
      lostO = p%orgC(i,1) * coef(i)*timestep
      p%orgC(i,1) = p%orgC(i,1) - lostO
      inventory(1) = inventory(1) - lostO*p%n
      lost(1) = lost(1) + lostO*p%n
      IF (bacZoo .eq. 1) THEN
         orgC_b(x) = orgC_b(x) + lostO*p%n
      ENDIF
      i = i + 1
      lostOrg = lostOrg + lostO
   ENDDO
ENDIF

IF (TEPC .gt. 0) THEN
   i = 1
   DO WHILE (i .le. 10)
      lostTEP = p%TEP(i,1) * coefTEP(i)*timestep
      p%TEP(i,1) = p%TEP(i,1) - lostTEP
      inventory(1) = inventory(1) - lostTEP*p%n
      lost(5) = lost(5) + lostTEP*p%n
      i = i + 1
      lostOrg = lostOrg + lostTEP
   ENDDO
ENDIF
   IF (lostOrg .gt. organic+TEPC) THEN
      print*, 'RESPIRATION - EMERGENCY', p%id
      print*, lostOrg, organic, coef
   ENDIF
   IF (p%TEP(1,1) .lt. 0) print*, 'respiration(TEP)', p%id

END SUBROUTINE respiration
!========================================================================
SUBROUTINE microzoop(p, zoop, harvest, paraT, paraD, MLdepth, testVar)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: zoop   !molC/m^3
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, paraT, paraD, MLdepth, testVar
   DOUBLE PRECISION :: zoo, dv, dv1, v, P_enc, orgCfrac, food, fact
   DOUBLE PRECISION :: P_break, apetite, microZ, depth, depth2
   INTEGER :: i, flag

   microZ = zoop
   depth = p%z/100.0
   flag = 0

   depth2 = depth/15.0
!   fact = 0.1*(depth2**0.5*exp(-depth2/80.0) - 3.9)
   fact = -4.5 
   P_enc = 10**fact * (-1/log(microZ)) * timestep !* testVar
   IF (-1/log(microZ) .gt. 1) THEN
      print*, 'ERROR: microzoop:', P_enc, microZ, depth
      P_enc = 0.9
   ENDIF

   CALL random_number(harvest)
   IF (P_enc .gt. harvest) THEN

      P_break = DATAN(p%r/1e4) !* 0.1 * testVar
      CALL random_number(harvest)
      harvest2 = harvest
      CALL random_number(harvest)
      CALL orgC_fraction(p, orgCfrac)
      apetite = orgCfrac 
      IF (P_break .gt. harvest .and. p%Nn .gt. 1) THEN
         CALL break_p(p, harvest)
      ELSEIF (apetite .gt. harvest2) THEN
         CALL ingestion(p, 1, paraT, paraD, harvest, MLdepth, testVar)
      ENDIF
   ENDIF
 
END SUBROUTINE microzoop
!========================================================================
SUBROUTINE ingestion(p, flagFP, paraT, paraD, harvest,MLdepth, testVar)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: paraT, paraD, MLdepth, testVar
   INTEGER, INTENT(IN) :: flagFP !FP: 1=micro 
   TYPE(agg_part)  :: p_old
   DOUBLE PRECISION :: organic, lostO, lostC, lostA, orgLost
   DOUBLE PRECISION :: calcite, aragonite, orgAssim, calcDiss, aragDiss
   DOUBLE PRECISION :: TEP, lostTEP!, S_c, S_a!, k_size
   INTEGER :: i, x

   p_old = p
   organic = 0
   i = 1
   TEP=0
   DO WHILE (i .le. 10)
      organic = organic + p%orgC(i,1)
      TEP = TEP + p%TEP(i,1)
      i = i + 1
   ENDDO
   IF (organic .lt. 0) print*, 'ERROR: ingestion:', organic, p%orgC(1,1)
   IF (TEP .lt. 0) print*, 'ERROR: ingestion(TEP)', TEP, p%TEP
   calcite = p%mineral(1)
   aragonite = p%mineral(2)
   
   CALL find_depth(p%z, x)

   calcDiss = 0.02*paraD
   aragDiss = calcDiss/5 !* paraD

! zooplankton organic carbon and TEP respiration
   IF (organic .gt. 0) THEN
      IF (flagFP .eq. 1) THEN     !micro
         p%af = 1
         CALL respiration(p, 2, lostO, paraT, MLdepth, testVar)
         orgC_G1(x) = orgC_G1(x) + lostO*p%n
      ENDIF 
 !zooplankton calcite dissolution
      IF (calcite .gt. 0) THEN
         lostC = lostO/(organic+TEP)*calcite*calcDiss
         IF (lostC .gt. calcite) THEN
            lostC = calcite
         ENDIF
         p%mineral(1) = p%mineral(1) - lostC
         lost(2) = lost(2) + lostC*p%n
         inventory(2) = inventory(2) - lostC*p%n
         calc_G1(x) = calc_G1(x) + lostC*p%n
      ENDIF
 !zooplankton aragonite dissolution
      IF (aragonite .gt. 0) THEN
         k_size = (p%mineral(2)/(p%mineral(2)+1e-10))
         S_a = 0.1
         lostA = lostO/(organic+TEP)*aragonite*aragDiss
!         lostA = lostO*aragDiss * (aragonite/(aragonite+1e-10))&
!                 *aragonite*timestep
!         lostA = k_size*k_caco3_a*S_a*p%mineral(1)*timestep
         IF (lostA .ge. aragonite) THEN
            lostA = aragonite
         ENDIF
         p%mineral(2) = p%mineral(2) - lostA
         lost(2) = lost(2) + lostA*p%n
         inventory(2) = inventory(2) - lostA*p%n
         calc_G1(x) = calc_G1(x) + lostA*p%n
      ENDIF
   ENDIF
CALL fecal_pellet(p, harvest)

i = 1
DO WHILE (i .lt. 10) 
   IF (p%orgC(i,1) .lt. 0) THEN
      p%orgC(i,1) = 0.0
   ENDIF
   IF (p%orgC(i,1) .lt. 1d-50 .and. p%orgC(i,1) .gt. 0.0) THEN
!      print*, 'tiny', p%orgC(i,1), i
      inventory(1) = inventory(1) - p%orgC(i,1)*p%n 
      lost(1) = lost(1) + p%orgC(i,1)*p%n 
      p%orgC(i,1) = 0.0
   ENDIF
   i = i + 1
ENDDO

END SUBROUTINE ingestion
!========================================================================
SUBROUTINE photolysis(p)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   
   IF (p%rho .lt. rho_sw) THEN
      p%r = 0.1
      CALL dissolve(p)
   ENDIF
   
END SUBROUTINE photolysis
!========================================================================
SUBROUTINE temperature(temp, depth, paraT, MLdepth)
   USE the_info

! depth is in cm, temp in C

   DOUBLE PRECISION, INTENT(OUT) :: temp
   DOUBLE PRECISION, INTENT(IN) :: depth, paraT, MLdepth
   DOUBLE PRECISION, DIMENSION(2,3) :: iTemp   
   DOUBLE PRECISION :: z, MLD
   INTEGER :: i

   MLD = MLdepth / 100.0  ! MLdepth is in [cm], MLD is in [m]

   iTemp(1,:) = [paraT,4,2]  
   iTemp(2,:) = [MLD,1000,4000]

   z = depth/100  !m
   i = 1
   DO WHILE ( z .ge. iTEMP(2,i) ) 
      i = i + 1
   ENDDO

   IF (depth .lt. MLdepth) THEN
     temp = paraT
   ELSE
      CALL interpolate(iTemp(1,i-1),iTemp(1,i),iTemp(2,i-1),iTemp(2,i),&
        temp,z)
   ENDIF

   IF (temp .gt. 35) THEN
      print*, 'temperature greater than 35'
      print*, iTEMP(1,1), iTemp(2,1), temp
   ENDIF

END SUBROUTINE temperature
!========================================================================
SUBROUTINE aging(p)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION :: older
   INTEGER :: i

   i = 9
   older = 0

   DO WHILE (i .ge. 1) 
      IF (i .eq. 1) THEN
         older = p%orgC(i,1) * timestep/(p%orgC(i,2))
      ELSE
         older = p%orgC(i,1) * timestep/(p%orgC(i,2)-p%orgC(i-1,2))
      ENDIF
      IF (older .gt. p%orgC(i,1)) THEN
         p%orgC(i+1,1) = p%orgC(i+1,1) + p%orgC(i,1)
         p%orgC(i,1) = 0.0
      ELSE
         p%orgC(i,1) = p%orgC(i,1) - older
         p%orgC(i+1,1) = p%orgC(i+1,1) + older
      ENDIF
      i = i - 1
   ENDDO
   
   i = 9
   older = 0
   DO WHILE (i .ge. 1) 
      IF (i .eq. 1) THEN
         older = p%TEP(i,1) * timestep/(p%TEP(i,2))
      ELSE
         older = p%TEP(i,1) * timestep/(p%TEP(i,2)-p%TEP(i-1,2))
      ENDIF
      IF (older .gt. p%TEP(i,1)) THEN
         p%TEP(i+1,1) = p%TEP(i+1,1) + p%TEP(i,1)
         p%TEP(i,1) = 0.0
      ELSE
         p%TEP(i,1) = p%TEP(i,1) - older
         p%TEP(i+1,1) = p%TEP(i+1,1) + older
      ENDIF
      i = i - 1
   ENDDO

END SUBROUTINE aging
!========================================================================
SUBROUTINE orgC_fraction(p, dwOrgCfrac)
! this subroutine calculates the organic carbon dry weight as a fraction
! of total weight

   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: dwOrgCfrac
   DOUBLE PRECISION :: organic, dw_mineral, dw_orgC, dw_total
   DOUBLE PRECISION ::  TEP, dw_TEP

   organic = 0
   TEP = 0
   dw_mineral = (p%mineral(1)+p%mineral(2))*mw_co + &
             p%mineral(3)*mw_si + p%mineral(4)*mw_li
   i = 1
   DO WHILE (i .le. 10) 
      organic = organic + p%orgC(i,1)
      TEP = TEP + p%TEP(i,1)
      i = i + 1
   ENDDO
   dw_orgC = organic*mw_orgC
   dw_TEP = TEP*mw_TEP

   dw_total = dw_orgC + dw_mineral + dw_TEP

   dwOrgCfrac = dw_orgC / dw_total

END SUBROUTINE orgC_fraction
!========================================================================
SUBROUTINE dissolution(p, paraC, paraT, MLdepth)

   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: paraC, paraT, MLdepth
   INTEGER :: x
   DOUBLE PRECISION :: R, Rd, Ta, temp, S_c, S_a
   DOUBLE PRECISION :: lostC, lostS, lostA, k_size

   lostA = 0
   lostC = 0
   lostS = 0
   S_c = 0
   S_a = 0

   IF (p%mineral(1)+p%mineral(2)+p%mineral(3) .gt. 0) THEN
      CALL temperature(temp, p%z, paraT, MLdepth)
      CALL deltaCO3(p%z, S_c, S_a, paraC)

      Ta = temp+273.15
      R =  1.32e16*exp(-11481/Ta)
      Rd = R/86400

      CALL find_depth(p%z, x)

! CaCO3 thermodynamics
 ! if amount of CaCO3 in agg less than in 1/2 coccolithoph
 ! then the whole thing dissolves, otherwise according to thermodyn.
      IF (S_c .gt. 0) THEN
         IF (p%mineral(1) .gt. 2.5e-12) THEN
            k_size = (p%mineral(1)/(p%mineral(1)+1e-10)) 
            lostC = k_caco3_c*S_c*p%mineral(1)*timestep
!            lostC = k_size*k_caco3_c*S_c*p%mineral(1)*timestep
!            IF (p%mineral(1) .lt. 2.5e-12 .or. p%mineral(1) .lt. lostC) THEN  
!               lostC = p%mineral(1)
!            ENDIF
            p%mineral(1) = p%mineral(1) - lostC
            inventory(2) = inventory(2) - lostC*p%n
            lost(2) = lost(2) + lostC*p%n
            calcC_t(x) = calcC_t(x) + lostC*p%n
         ENDIF
      ENDIF
! aragonite dissolution       
      IF (S_a .gt. 0) THEN
         IF (p%mineral(2) .gt. 2.5e-12) THEN
            k_size = (p%mineral(2)/(p%mineral(2)+1e-10)) 
            lostA = k_caco3_a*S_a*p%mineral(2)*timestep
!            lostA = k_size*k_caco3_a*S_a*p%mineral(2)*timestep
            IF (lostA .gt. p%mineral(2)) THEN
               print*, p%mineral(2), lostA, S_a
               print*, '-------mineral(2)-----------', p%z/100
            ENDIF
!            IF (p%mineral(2) .lt. 2.5e-12 .or. lostA .gt. p%mineral(2)) THEN  
!               lostA = p%mineral(2)
!            ENDIF
            p%mineral(2) = p%mineral(2) - lostA
            inventory(2) = inventory(2) - lostA*p%n
            lost(2) = lost(2) + lostA*p%n
            calcA_t(x) = calcA_t(x) + lostA*p%n
         ENDIF
      ENDIF

! opal dissolution
      IF (p%mineral(3) .gt. 0) THEN
         k_size = (p%mineral(3)/(p%mineral(3)+1e-10))   !DEFAULT
!          k_size = (p%mineral(3)/(p%mineral(3)+1e-9))  !ALTERNATE
         lostS = k_size*p%mineral(3)*Rd*timestep
         p%mineral(3) = p%mineral(3) - lostS
         inventory(3) = inventory(3) - lostS*p%n
         lost(3) = lost(3) + lostS*p%n
         opal_t(x) = opal_t(x) + lostS*p%n
      ENDIF


! volume of solid stuff
      i = 1
      orgC = 0d00
      DO WHILE (i .le. 10) 
         orgC = orgC + p%orgC(i,1) 
         i = i + 1
      ENDDO
      calc = p%mineral(1) + p%mineral(2)
      opal = p%mineral(3)
      clay = p%mineral(4)
  
! change moles to volume
      vorgC = orgC * mw_orgC / rho_orgC
      vcalc = calc * mw_co / rho_co
      vopal = opal * mw_si / rho_si
      vclay = clay * mw_li / rho_li

      tvol = vorgC + vcalc + vopal + vclay
   
      IF (lostS+lostC+lostA .gt. 0) THEN
         CALL radius(p)
      ENDIF
   ENDIF

IF (p%mineral(1) .lt. 0) print*, 'dissol', p%mineral(2), p%TEP(1,1)

END SUBROUTINE dissolution
!=======================================================================
SUBROUTINE deltaCO3(depth, S_c, S_a, paraC)
   USE the_info

   DOUBLE PRECISION, INTENT(IN) :: depth, paraC
   DOUBLE PRECISION, INTENT(INOUT) :: S_c, S_a
   DOUBLE PRECISION :: co3ion, co3ion1, co3ion2, co3ion3, co3iond
   DOUBLE PRECISION :: co3satC, co3satA, z

   z = depth/100  !m

   co3ion1 = paraC
   co3ion2 = paraC - 140
   co3ion3 = paraC - 140 !paraC/5 + 40
   co3iond = 1000
!   co3ion1 = 210 !PACIFIC (umol/kg)
!   co3ion2 = 65 
!   co3ion3 = 80
!   co3iond = 1000

!   co3ion1 = 250 !NATLANTIC(umol/kg)
!   co3ion2 = 150 
!   co3ion3 = 100
!   co3iond = 1000

!   co3ion1 = 230 !SOUTHERN(umol/kg)
!   co3ion2 = 100 
!   co3ion3 = 90
!   co3iond = 1000

   IF (z .lt. 2500) THEN
      co3satC = 0.01104*z + 41.8
      co3satA = 0.0168*z + 65
   ELSE
      co3satC = 0.01692*z + 27.1
      co3satA = 0.024*(z - 2500) + 107
   ENDIF

   IF (z .lt. co3iond) THEN
      slope = (z - co3iond)/((-co3iond)*2/pi)
      co3ion = co3ion1 + (co3ion2-co3ion1)*cos(slope)
   ELSE
      co3ion = co3ion2+(co3ion3-co3ion2)*tan(0.86*(z-co3iond)/(seafloor/100))
   ENDIF
   
   IF (co3ion .gt. co3satC) THEN
      S_c = 0
   ELSE
      S_c = (1 - co3ion/co3satC)**eta_c
   ENDIF

   IF (co3ion .gt. co3satA) THEN
      S_a = 0
   ELSE
      S_a = (1 - co3ion/co3satA)**eta_a
   ENDIF

END SUBROUTINE deltaCO3
!========================================================================
SUBROUTINE viscosity(depth, nu, paraT, MLdepth)
! dynamic viscosity in units g/cms
   USE the_info

   DOUBLE PRECISION, INTENT(IN) :: depth, paraT, MLdepth
   DOUBLE PRECISION, INTENT(OUT) :: nu
   DOUBLE PRECISION :: temp

   CALL temperature(temp, depth, paraT, MLdepth)

   nu = 5e-6 * exp(2250/(temp+273.15))
   IF (nu .lt. 1e-3 .or. nu .gt. 1) THEN
      print*, 'viscosity'
      print*, depth, temp, nu
   ENDIF

END SUBROUTINE
!========================================================================
SUBROUTINE mixedLayer(p, harvest, MLdepth)

   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, MLdepth
 
   CALL random_number(harvest)
 
   p%z = harvest * MLdepth
   
END SUBROUTINE
