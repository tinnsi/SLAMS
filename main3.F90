#define respi
#define dissol
#define microzoo
#define agin
#define sinking
#define dissolved
#define disintegr
#define housekeeping
#define selfCollision
#define photo
#define writing
#define array_size

MODULE the_info
   IMPLICIT NONE
   SAVE

   TYPE agg_part
      INTEGER :: id, b, af !identification number, bottom, agg or fecal pell
! b=0: in watercolumn, b=1: on bottom, b=2: has been reset
! af=0: aggregate
! af=1: don't do this. bad idea.
      DOUBLE PRECISION :: r, rho, por, s, z, w, n, p, Nn 
      DOUBLE PRECISION :: pr, frac
      DOUBLE PRECISION, DIMENSION(4) :: mineral
      DOUBLE PRECISION, DIMENSION(10,2) :: orgC, TEP
      !radius, density, porosity, stickiness, depth, settl. vel
      ![um,g/cm3,-,-,cm, cm/s, #,#_total,#/n
      ! n = number of aggregates this aggregate represents in m3
      ! p = number of primary particles in all the aggregages
      ! Nn = number of primary particles in each aggregate
      ! orgC : [molC/aggregate]
      ! minerals(calc, arag, opal, clay) 
   END TYPE

   DOUBLE PRECISION, DIMENSION(4) :: inventory = [0,0,0,0]
   DOUBLE PRECISION, DIMENSION(5) :: sf = [0,0,0,0,0] !seafloor
   DOUBLE PRECISION, DIMENSION(5) :: npp = [0,0,0,0,0] !net primary prod
   DOUBLE PRECISION, DIMENSION(5) :: lost = [0,0,0,0,0]   ! stuff lost in wc: org,calc,opal,clay,TEP
   DOUBLE PRECISION, DIMENSION(26,14) :: sedTrap
   DOUBLE PRECISION, DIMENSION(19) :: seaBed

! timestep should not be greater than 1 day (eg. disintegrate)
   INTEGER, PARAMETER :: timestep = 3600*8   !s
   INTEGER, PARAMETER :: endoftime = 1095*18
   INTEGER, PARAMETER :: seafloor = 400000    !cm
   INTEGER, PARAMETER :: dz = 1000         !cm
   INTEGER, PARAMETER :: n_db = 400  !number of depth bins in the watercolumn
   INTEGER, PARAMETER :: n_size = 60 !number of size bins for aggregates

! how much orgC and CaCO3 is cons by bact and zoop at each depth 
   DOUBLE PRECISION,DIMENSION(n_db) :: orgC_G1, orgC_b, opal_t
   DOUBLE PRECISION,DIMENSION(n_db) :: calcC_t, calcA_t, calc_G1
   DOUBLE PRECISION, DIMENSION(n_size) :: radi
   DOUBLE PRECISION, DIMENSION(n_size) :: sizeBed, veloBed, densBed
   DOUBLE PRECISION, DIMENSION(2,n_size) :: scaling
! physical constants
   DOUBLE PRECISION, PARAMETER :: a_number=200000
   DOUBLE PRECISION, PARAMETER :: k = 1.38d-23
   DOUBLE PRECISION, PARAMETER :: T = 285d00 !K
   DOUBLE PRECISION, PARAMETER :: mu = 8.9d-4, gamma = 0.01
   DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d00
   DOUBLE PRECISION, PARAMETER :: rho_sw = 1.02695d00
   DOUBLE PRECISION, PARAMETER :: rho_orgC = 1.06d00, mw_orgC = 28.8d00
   DOUBLE PRECISION, PARAMETER :: rho_co = 2.8275d00, mw_co = 100.1d00
   DOUBLE PRECISION, PARAMETER :: rho_si = 1.65d00, mw_si = 67.3d00
   DOUBLE PRECISION, PARAMETER :: rho_li = 2.65d00, mw_li = 60d00
   DOUBLE PRECISION, PARAMETER :: rho_TEP = 0.8, mw_TEP = 154.8d00 ! 154.8
   DOUBLE PRECISION, PARAMETER :: g = 981 !cm/s
   DOUBLE PRECISION, PARAMETER :: k_caco3_c = 5.0/86400
   DOUBLE PRECISION, PARAMETER :: k_caco3_a = 3.2/86400
   DOUBLE PRECISION, PARAMETER :: eta_c = 4.5
   DOUBLE PRECISION, PARAMETER :: eta_a = 4.2

END MODULE the_info
!====================================================================

PROGRAM main

   USE the_info
   IMPLICIT NONE

   102 FORMAT (E15.4E2, 1X, E15.4E2, 1X, I)
   103 FORMAT (E15.4E2, 1X, I5, 1X, E15.4E2)
   110 FORMAT (I4)
   111 FORMAT (E15.4E2)
   120 FORMAT (E15.6E2, 3X, F14.4, 1X, I)
   125 FORMAT (E15.6E2, 3X, F14.4, 1X, I, 1X, I)
   146 FORMAT (F, 1X, 15(E15.6E2, 1X))
   220 FORMAT (I5, 1X, 4E, 1X, I5)
   240 FORMAT (5E)
   250 FORMAT (I, 1X, 4E)
   
   INTEGER :: ka, time, i, j, l, agg_max, z_max, temp1, temp2, tmp1, tmp2
   INTEGER :: intT, intC, intD, intP, intS, intTEP, intTEP2, tsV
   INTEGER :: intFrac, yaff, resol
   INTEGER :: UNIT, free_count, loop
   INTEGER, DIMENSION(1) :: old, seed
   INTEGER, DIMENSION(a_number), TARGET :: agg_z
   INTEGER, DIMENSION(10), TARGET :: agg_pool
   INTEGER, DIMENSION(a_number) :: free_agg
   INTEGER, POINTER :: pz
   DOUBLE PRECISION, DIMENSION(1) :: harvest
   DOUBLE PRECISION :: depth, temp, dCO3c, dCO3a, max_r, organic, dummy
   DOUBLE PRECISION :: paraT, paraC, paraD, paraP, paraS
   DOUBLE PRECISION :: paraTEP
   DOUBLE PRECISION :: orgC1, orgC4, caco31, caco34, food
   DOUBLE PRECISION :: Fsize, resolution
   DOUBLE PRECISION :: before, after, TEP_frac, year, season
   DOUBLE PRECISION :: mass, day
   DOUBLE PRECISION, DIMENSION(n_db) :: orgCtotal, calctotal, rainratio
   DOUBLE PRECISION, DIMENSION(n_db) :: orgCMean, calcMean, rainMean
   DOUBLE PRECISION, DIMENSION(n_db) :: opalMean, clayMean
   DOUBLE PRECISION, DIMENSION(n_db) :: mass_sum
   DOUBLE PRECISION, DIMENSION(5) :: mass_i, mass_f
   DOUBLE PRECISION, DIMENSION(15,n_db) :: data_stuff
   DOUBLE PRECISION, DIMENSION(n_size,n_db) :: spec
   DOUBLE PRECISION, DIMENSION(n_size,n_db,2,10) :: data_mega
   DOUBLE PRECISION, DIMENSION(4), TARGET :: data_flux
   DOUBLE PRECISION, POINTER :: p_df(:) !pointer to data_flux
   DOUBLE PRECISION, DIMENSION(n_db) :: mic
   TYPE(agg_part), DIMENSION(a_number), TARGET :: agg
   TYPE(agg_part), POINTER :: pa, pa1, pa2
   CHARACTER (len=14) :: ctime, parP, parC, parT, parS
   CHARACTER (len=30) :: file_name, spec_name

   READ(5,*) intT,intC,intD,intP,intS,intTEP
   
   paraT = intT    !SST
   paraC = intC    !CO3=
   paraD = intD    !zoopl. dissol.
   paraP = intP    !primary prod.
   paraS = intS    !seasonality
   paraTEP = intTEP!fraction of pp that is TEP

#ifdef writing
   CALL open_files()
#endif

! random numbers
   ka = 1
   old = 1
   seed(1) = 1834
   CALL random_seed
   CALL random_seed(size=ka)
   CALL random_seed(GET=old(1:ka))
   CALL random_number(harvest(1))

! Initialize vectors and variables
   CALL initialize(time, day, year, agg_max, free_count, free_agg, tsV, mic, mass_i)

! Fraction of organic carbon that is TEP
   TEP_frac = paraTEP/100.0 !0.2

   DO WHILE (time .le. endoftime) 
      i = 1
      loop = FLOOR(20*1095/paraS) 
      organic = 0
      IF (time .eq. 1095*(year+1)) THEN
         year = year + 1
      ENDIF
      IF (time-year*1095 .lt. paraS) THEN
         season = paraP/12.0e3 * 1095/paraS   !prod per day [molC/d]
      ELSE
         season = 0
      ENDIF
! primary production 
      DO WHILE (i .le. loop .and. season .gt. 0)
         IF (free_count .gt. 1) THEN
            pa => agg(free_agg(free_count-1))
            pa%id = free_agg(free_count-1)
            free_count = free_count-1
         ELSEIF (free_count .eq. 1) THEN
            pa => agg(agg_max)
            pa%id = agg_max
            agg_max = agg_max + 1
         ENDIF
         IF (pa%id .gt. a_number) print*, 'END OF LIST', pa%id
         CALL random_number(harvest)
         CALL production(pa, harvest, loop, time, season, paraT)
                        
         mass_i(1) = mass_i(1) + pa%orgC(1,1)*pa%n
         mass_i(2) = mass_i(2) + (pa%mineral(1)+pa%mineral(2))*pa%n
         mass_i(3) = mass_i(3) + pa%mineral(3)*pa%n
         mass_i(4) = mass_i(4) + pa%mineral(4)*pa%n
         organic = organic + pa%orgC(1,1)*pa%n
         i = i + 1
      ENDDO 
! TEP production
      IF (season .gt. 0) THEN
         IF (free_count .gt. 1) THEN
            pa => agg(free_agg(free_count-1))
            pa%id = free_agg(free_count-1)
            free_count = free_count-1
         ELSEIF (free_count .eq. 1) THEN
            pa => agg(agg_max)
            pa%id = agg_max
            agg_max = agg_max + 1
         ENDIF
         CALL random_number(harvest)
         CALL TEP_prod(pa,organic,harvest, TEP_frac, paraT)
         mass_i(5) = mass_i(5) + pa%TEP(1,1)*pa%n
      ENDIF

! find agg from same depth and put them into an array
      depth = 0
      DO WHILE (depth .lt. seafloor) 
         z_max = 0
         agg_z = 0
         i = 1
         j = 1
         DO WHILE(i .lt. agg_max) 
            IF (agg(i)%b .eq. 0) THEN
               pz => agg_z(j)
               pa => agg(i)
               CALL check_depth(pa, pz, depth, i, j, z_max) 
            ENDIF
            i = i + 1
         ENDDO 
! make aggregates at similar depth aggregate
         i = 1
         j = 1
         DO WHILE (i .le. z_max) 
            DO WHILE (j .le. z_max) 
            CALL random_number(harvest)
            tmp1 = agg_z(FLOOR(harvest(1) * z_max))
            IF (FLOOR(harvest(1)*z_max) .eq. 0) tmp1 = agg_z(z_max)
            CALL random_number(harvest)
            tmp2 = agg_z(FLOOR(harvest(1) * z_max))
            IF (FLOOR(harvest(1)*z_max) .eq. 0) tmp2 = agg_z(z_max)
            pa1 => agg(tmp1)
            pa2 => agg(tmp2)
            IF (pa1%b .eq. 0 .and. pa2%b .eq. 0) THEN
               IF (tmp1 .eq. tmp2) THEN
#ifdef selfColl 
                  CALL self_collide(pa1, harvest)
#endif
               ELSEIF (tmp1 .ne. tmp2) THEN
                  CALL collide4(pa1, pa2, harvest, paraT)
               ENDIF
            ENDIF
            j = j + 1
         ENDDO 
         i = i + 1
         j = i + 1
      ENDDO 

! make zooplankton array using the agg_z() created above
         i = 1
         food = 0
         DO WHILE (i .le. z_max) 
            temp1 = agg_z(i)
            pa1 => agg(temp1)
            IF (pa1%af .eq. 0) THEN
               j = 1
               DO WHILE (j .le. 5) 
                  food = food + pa1%orgC(j,1)*pa1%n
                  j = j + 1
               ENDDO 
            ENDIF
            i = i + 1
         ENDDO 

         CALL find_depth(depth, i)

         food = food * 1/(dz/100.0)
         mic(i) = food                 !molC/m3

#ifdef housekeeping
         i = 1
         j = 1
         agg_pool = 0
         DO WHILE (i .le. z_max) 
            temp1 = agg_z(i)
            pa1 => agg(temp1)
            IF (pa1%Nn .eq. 1 .and. pa1%s .lt. 0.01 .and. &
                pa1%mineral(4) .lt. 2.8d-13) THEN
              agg_pool(j) = temp1
              j = j + 1
            ENDIF
            i = i + 1
         ENDDO 
         IF (j .gt. 3) THEN
            i = j - 2
            DO WHILE (i .ge. 1) 
               temp1 = agg_pool(j-2)
               temp2 = agg_pool(j-1)
               CALL pool(agg(temp1), agg(temp2))
               i = i - 2
               j = j - 2
            ENDDO 
         ENDIF
#endif
         depth = depth + dz
      ENDDO 
#ifdef writing
! write out the zooplankton array
         IF (MOD(time,1000) .eq. 0) THEN 
print*, 'inside', mic(4), size(mic)
            CALL write_zooplankton(time, mic)
print*, 'after write', size(mic)
         ENDIF
#endif
         IF (MOD(time,90) .eq. 0) THEN
            CALL write_npp(time, agg_max)
         ENDIF
         IF (MOD(time,1095) .eq. 0) THEN
print*, 'hello', time
            CALL write_sedTrap(time, day)
         ENDIF

! update density, velocity, ... for all particles
      i = 1 !aggegate
      j = 1 !depth
      l = 1 !sizebin
      DO WHILE (i .lt. agg_max) 
         IF (agg(i)%b .eq. 0) THEN
            pa => agg(i)
            CALL find_depth(pa%z, j)
            CALL find_size(pa%r, l, n_size)
            CALL stickiness(pa)
            CALL density(pa)
#ifdef microzoo
            IF (mic(j) .gt. 0) THEN
               CALL microzoop(pa,mic(j),harvest,paraT,paraD)
            ENDIF
#endif
#ifdef respi
            CALL respiration(pa, 1, dummy, paraT)
#endif 
#ifdef dissol
            CALL dissolution(pa,paraC,paraT)
#endif
#ifdef disintegr
            CALL disintegrate(pa, harvest)
#endif
#ifdef agin
            CALL aging(pa)
#endif
#ifdef sinking
            CALL fractal_dimension(pa)
            CALL radius(pa)
            CALL velocity(pa, paraT)
            CALL sink(pa, time)
            CALL bottom(pa,l)
#endif
#ifdef dissolved
            CALL dissolve(pa)
#endif
#ifdef photo
            IF (pa%z .le. 100 .and. pa%r .gt. 1e4) THEN
               CALL photolysis(pa)
            ENDIF
#endif
            CALL stickiness(pa)
         ENDIF
         IF (agg(i)%b .eq. 1) THEN
            free_agg(free_count) = i
            free_count = free_count + 1
            pa => agg(i)
            CALL reset(pa)
         ENDIF
         i = i + 1
      ENDDO 
! figure out how much stuff there is at each depth
      i = 1
      data_stuff = 0
      spec = 0
      DO WHILE (i .le. n_db) 
         data_stuff(1,i) = i*dz
         i = i + 1
      ENDDO 
      i = 1
      DO WHILE (i .lt. agg_max) 
         pa => agg(i)
         IF (pa%z .lt. seafloor) THEN
            CALL collect(pa, data_stuff, spec)
         ENDIF
         i = i + 1
      ENDDO  


! averaging the stuff in the watercolumn of the last tsV timesteps 
      IF(time .gt. endoftime - tsV) then
         CALL eotAverage(data_stuff, orgCMean, calcMean, opalMean, clayMean, rainMean)
      ENDIF

! calculate all the properties to make pictures
#ifdef writing
      IF (time .eq. endoftime ) THEN
         data_mega(:,:,:,:) = 0
         i = 1
         DO WHILE (i .lt. agg_max) 
            pa => agg(i)
            IF (pa%z .lt. seafloor) THEN 
               CALL collect2(pa, data_mega)
               IF (pa%b .eq. 1) THEN 
                  print*, 'oppsie!', pa%id 
               ENDIF
            ENDIF
            i = i + 1
         ENDDO  

         CALL write_propPic(data_mega)
         CALL write_seaFloor(seaBed, agg(1))
      ENDIF
#endif
      
! write the whole particle spectrum picture
#ifdef writing
      IF (MOD(time,1000) .eq. 0) THEN
print*, 'write_spec', spec(1,1)
         CALL write_specPic(spec, time)
      ENDIF
#endif

#ifdef writing
      i = 1
      DO WHILE (i .le. agg_max)
         IF (ABS(agg(i)%z - 100000) .le. dz/2.0) THEN
            pa => agg(i)
            CALL calculate_flux(pa, data_flux)
         ENDIF
         i = i + 1
      ENDDO
     
      p_df => data_flux
      CALL write_dataFlux(p_df, time)
#endif

#ifdef writing
      IF (time .gt. 50 .and. MOD(time,100) .eq. 0) THEN
         CALL write_sf(time)
      ENDIF
      CALL write_singlePart(agg(1), time, 48)
      CALL write_singlePart(agg(2), time, 49)
      CALL write_singlePart(agg(3), time, 50)
      CALL write_singlePart(agg(11), time, 51)
      
#endif
      time = time + 1

      IF (time .eq. endoftime - 365*3) THEN
         sizeBed(:) = 0
         veloBed(:) = 0
         densBed(:) = 0
         orgC_G1(:) = 0
         orgC_b(:) = 0
         calc_G1(:) = 0
         calcC_t(:) = 0
         calcA_t(:) = 0
         opal_t(:) = 0
         scaling(:,:) = 0
      ENDIF
   ENDDO
print*, 'hello - 1'
      
!write out the mean components of the component of the watercolumn.
   CALL write_eotA(orgCMean, calcMean, opalMean, clayMean, rainMean, tsV)
print*, 'after write eotA' 
   i = 1
   DO WHILE (i .lt. agg_max) 
      j = 1
      DO WHILE (j .le. 10) 
         mass_f(1) = mass_f(1) + agg(i)%orgC(j,1)*agg(i)%n
         mass_f(5) = mass_f(5) + agg(i)%TEP(j,1)*agg(i)%n
         j = j + 1
      ENDDO
      mass_f(2) = mass_f(2) + (agg(i)%mineral(1)+agg(i)%mineral(2))*agg(i)%n
      mass_f(3) = mass_f(3) + agg(i)%mineral(3)*agg(i)%n
      mass_f(4) = mass_f(4) + agg(i)%mineral(4)*agg(i)%n
      CALL find_size(agg(i)%r, l, n_size)
      scaling(1,l) = scaling(1,l) + 1
      scaling(2,l) = scaling(2,l) + agg(i)%n
      i = i + 1
   ENDDO

   print*, npp(1), sf(1), 'orgC prod - sf'
   print*, npp(2), sf(2), 'calc prod - sf'
   print*, npp(3), sf(3), 'opal prod - sf'
   print*, npp(4), sf(4), 'clay prod - sf'
   print*, npp(5), sf(5), 'TEP prod - sf'

   print*, mass_i(1), lost(1), mass_f(1) + lost(1) + sf(1) 
   print*, mass_i(2), lost(2), mass_f(2) + lost(2) + sf(2) 
   print*, mass_i(3), lost(3), mass_f(3) + lost(3) + sf(3) 
   print*, mass_i(4), lost(4), mass_f(4) + lost(4) + sf(4) 
   print*, mass_i(5), lost(5), mass_f(5) + lost(5) + sf(5) 
print*, 'before write end'
   CALL write_end()

END PROGRAM
