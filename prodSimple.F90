SUBROUTINE production(pa, harvest, loop, time, prod, paraT)

! input is in terms of moles organic carbon.
! n is #/m3
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: pa
   INTEGER, INTENT(IN) :: loop, time
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, prod, paraT
   DOUBLE PRECISION :: a, b, c, d, x, y, z, w, f, dw, tw, org_ratio
   DOUBLE PRECISION :: labi, refr, calc, arag, opal, clay
   DOUBLE PRECISION :: vorgC, vcalc, vopal, vclay, tvol, age
   
   ! diatoms: 15mmol orgC, 5mmBSi
   ! coccoli: 7mmol orgC, 5mmCaCO3
   ! dust: 1mmol kaolinite
! this is 1000 mgC/m2d, then we multiply by x to get regional 
   pa%orgC = 0
   pa%mineral = 0
   IF (harvest .le. 0.04) THEN
      pa%orgC(1,1) = 7d-12
      pa%mineral(1) = 1.5d-12 
   ELSEIF (harvest .gt. 0.04 .and. harvest .le. 0.06) THEN
      pa%orgC(1,1) = 7d-12
      pa%mineral(2) = 1.5d-12 
   ELSEIF (harvest .gt. 0.06 .and. harvest .le. 0.3) THEN
      pa%orgC(1,1) = 15d-12 
      pa%mineral(3) = 5d-12 
   ELSEIF (harvest .gt. 0.3 .and. harvest .le. 0.99) THEN
     pa%orgC(1,1) = 1d-12
   ELSEIF (harvest .gt. 0.99 .and. harvest .le. 1) THEN
      pa%mineral(4) = 2.85d-13 
      pa%n = 3.39e5*timestep/loop  
   ENDIF
   IF (pa%orgC(1,1) .gt. 0) THEN 
      pa%n = prod/86400*timestep/loop/pa%orgC(1,1) 
   ENDIF
   pa%p = pa%n
   pa%TEP = 0d0
   pa%s = 0
   pa%Nn = 1
   pa%frac = 3d0
   CALL random_number(harvest)
   pa%z = 10000*exp(-5*harvest)
!   pa%z = 20000*exp(-4*harvest)
   pa%b = 0
   pa%af = 0
! age vector for orgC and TEP
   i = 1
   age = 1d00
   DO WHILE (i .le. 10)
      pa%orgC(i,2) = age*86400
      pa%TEP(i,2) = age*86400
      age = age * 2d00
      i = i + 1
   ENDDO

! volume of solid stuff
   i = 1
   orgC = 0
   DO WHILE (i .le. 10) 
      orgC = orgC + pa%orgC(i,1) 
      i = i + 1
   ENDDO
   calc = pa%mineral(1) + pa%mineral(2)
   opal = pa%mineral(3)
   clay = pa%mineral(4)
  
! change moles to volume
   vorgC = orgC * mw_orgC / rho_orgC
   vcalc = calc * mw_co / rho_co
   vopal = opal * mw_si / rho_si
   vclay = clay * mw_li / rho_li

   tvol = vorgC + vcalc + vopal + vclay

   pa%pr = (tvol*3/(4*pi))**(1/3.0) * 1e4
   pa%r = (tvol*3/(4*pi))**(1/3.0) * 1e4
   IF (pa%pr .gt. 1e4) THEN
      print*, 'production', tvol, pa%id
   ENDIF
   CALL fractal_dimension(pa)
   CALL density(pa)
   CALL velocity(pa, paraT)
   pa%Nn = 1

   npp(1) = npp(1) + pa%orgC(1,1)*pa%n
   npp(2) = npp(2) + pa%mineral(1)*pa%n+pa%mineral(2)*pa%n
   npp(3) = npp(3) + pa%mineral(3)*pa%n
   npp(4) = npp(4) + pa%mineral(4)*pa%n

   inventory(1) = inventory(1) + pa%orgC(1,1)*pa%n
   inventory(2) = inventory(2) + (pa%mineral(1)+pa%mineral(2))*pa%n
   inventory(3) = inventory(3) + pa%mineral(3)*pa%n
   inventory(4) = inventory(4) + pa%mineral(4)*pa%n

END SUBROUTINE production
!====================================================================
SUBROUTINE TEP_prod(p, organic, harvest, TEP_f, paraT)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, organic, TEP_f, paraT
   DOUBLE PRECISION :: tvol, age

   p%orgC = 0
   p%mineral = 0
   p%TEP = 0
   p%TEP(1,1) = 1d-12
   p%Nn = 1
   p%s = 1
   p%frac = 3d0
   p%n = TEP_f * organic/p%TEP(1,1)
   p%p = p%n
!   p%z = 20000*exp(-4*harvest)
   p%z = 10000*exp(-5*harvest)
   p%b = 0

! change moles to volume
   tvol = p%TEP(1,1) * mw_orgC / rho_orgC

   p%pr = (tvol*3/(4*pi))**(1/3.0) * 1e4
   p%r = (tvol*3/(4*pi))**(1/3.0) * 1e4
   i = 1
   age = 1d00
   DO WHILE (i .le. 10)
      p%orgC(i,2) = age*86400
      p%TEP(i,2) = age*86400
      age = age * 2d00
      i = i + 1
   ENDDO

   IF (p%pr .gt. 1e4) THEN
      print*, 'production', tvol, p%id
   ENDIF
   
   npp(5) = npp(5) + p%TEP(1,1)*p%n
   inventory(1) = inventory(1) + p%TEP(1,1)*p%n

   CALL fractal_dimension(p)
   CALL density(p)
   CALL velocity(p, paraT)

END SUBROUTINE TEP_prod
