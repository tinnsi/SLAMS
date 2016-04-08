SUBROUTINE initialize(time, day, year, agg_max, free_count, free_agg, tsV, mic, mass_i)
   USE the_info

   INTEGER, INTENT(INOUT) :: time, agg_max, tsV, free_count
   INTEGER, DIMENSION(a_number), INTENT(INOUT) :: free_agg
   DOUBLE PRECISION, INTENT(INOUT) :: day, year
   DOUBLE PRECISION, DIMENSION(n_db), INTENT(INOUT) :: mic
   DOUBLE PRECISION, DIMENSION(5), INTENT(INOUT) :: mass_i
   INTEGER :: i
   DOUBLE PRECISION :: max_r

   i = 1
   max_r = 1d00
   DO i = 1, n_size, 1
      radi(i) = max_r
      max_r = max_r * 1.3d00
   END DO

   time = 1
   day = 365
   year = 0

   agg_max = 1 !the next available spot
   free_count = 1
   free_agg = 0
   tsV = 90

   sedTrap(:,:) = 0
   mic(:) = 0.0

   mass_i(:) = 0.0

END SUBROUTINE initialize
!=============================================================================
SUBROUTINE collect(pa, data_mega, spec)
   USE the_info

   TYPE(agg_part), INTENT(IN) :: pa
   DOUBLE PRECISION, DIMENSION(15,n_db), INTENT(INOUT) :: data_mega
   DOUBLE PRECISION, DIMENSION(n_size,n_db), INTENT(INOUT) :: spec
   DOUBLE PRECISION :: organic 
   INTEGER :: i, j, l, aha

   ! find the depth at which the aggregates are
   CALL find_depth(pa%z, i)

   data_mega(2,i) = data_mega(2,i) + pa%orgC(1,1)*pa%n
   data_mega(3,i) = data_mega(3,i) + pa%orgC(2,1)*pa%n
   data_mega(4,i) = data_mega(4,i) + pa%orgC(3,1)*pa%n
   data_mega(5,i) = data_mega(5,i) + pa%orgC(4,1)*pa%n
   data_mega(6,i) = data_mega(6,i) + pa%orgC(5,1)*pa%n
   data_mega(7,i) = data_mega(7,i) + pa%orgC(6,1)*pa%n
   data_mega(8,i) = data_mega(8,i) + pa%orgC(7,1)*pa%n
   data_mega(9,i) = data_mega(9,i) + pa%orgC(8,1)*pa%n
   data_mega(10,i) = data_mega(10,i) + pa%orgC(9,1)*pa%n
   data_mega(11,i) = data_mega(11,i) + pa%orgC(10,1)*pa%n
   data_mega(12,i) = data_mega(12,i) + pa%mineral(1)*pa%n
   data_mega(13,i) = data_mega(13,i) + pa%mineral(2)*pa%n
   data_mega(14,i) = data_mega(14,i) + pa%mineral(3)*pa%n
   data_mega(15,i) = data_mega(15,i) + pa%mineral(4)*pa%n

   ! find the sizebin of the aggregates
   aha = 0
   l = 1
   DO WHILE (l .lt. n_size-1 .and. aha .eq. 0)
      IF (pa%r .lt. radi(l)) THEN
         aha = 1
      ELSE
         l = l + 1
      ENDIF
   ENDDO

!   spec(l,i) = spec(l,i) + pa%n
   organic = pa%orgC(1,1) + pa%orgC(2,1) + pa%orgC(3,1) + &
             pa%orgC(4,1) + pa%orgC(5,1) + pa%orgC(6,1) + &
             pa%orgC(7,1) + pa%orgC(8,1) + pa%orgC(9,1) + &
             pa%orgC(10,1) 

   spec(l,i) = spec(l,i) + pa%n


   if (spec(l+1,i) .gt. 1e18) THEN
      print*, 'spec:', spec(l,i), i
   endif
   if (l .eq. n_size) THEN
      print*, 'l eq n_size', pa%r, pa%s
   endif

END SUBROUTINE collect
!=========================================================================
SUBROUTINE reset(pa)
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: pa

   pa%b = 2
   pa%af = 0
   pa%r = 0
   pa%rho = 0 
   pa%por = 0
   pa%s = 0
   pa%z = 0
   pa%w = 0
   pa%n = 0
   pa%p = 0
   pa%Nn = 0
   pa%pr = 0
   pa%frac = 0
   i = 1
   DO WHILE (i .le. 10)
      pa%orgC(i,1) = 0
      pa%orgC(i,2) = 0
      pa%TEP(i,1) = 0
      pa%TEP(i,2) = 0
      i = i + 1
   ENDDO
   i = 1
   DO WHILE (i .le. 4) 
      pa%mineral(i) = 0
      i = i + 1
   ENDDO

END SUBROUTINE reset
!==========================================================================
SUBROUTINE calculate_flux(pa, data_flux)
   USE the_info

   TYPE(agg_part), INTENT(IN) :: pa
   DOUBLE PRECISION, DIMENSION(4), INTENT(INOUT) :: data_flux
   DOUBLE PRECISION :: organic, inorganic
   INTEGER :: i 
   
   i = 1
   organic = 0
   inorganic = 0
   DO WHILE (i .le. 10) 
     organic = organic + pa%orgC(i,1)*pa%n
     i = i + 1
   ENDDO
   inorganic = (pa%mineral(1)+pa%mineral(2))*pa%n

   data_flux(1) = data_flux(1) + organic*pa%w/dz*timestep
   data_flux(2) = data_flux(2) + inorganic*pa%w/dz*timestep
   data_flux(3) = data_flux(3) + pa%mineral(3)*pa%w/dz*timestep*pa%n
   data_flux(4) = data_flux(4) + pa%mineral(4)*pa%w/dz*timestep*pa%n

END SUBROUTINE calculate_flux
!==========================================================================
SUBROUTINE collect2(pa, data_mega)
   USE the_info

   TYPE(agg_part), INTENT(IN) :: pa
   DOUBLE PRECISION, DIMENSION(n_size,n_db,2,10), INTENT(INOUT) :: data_mega
   DOUBLE PRECISION :: organic, TEP, calc, rainratio
   INTEGER :: max_i, max_l
   INTEGER :: i, j, l, aha, x, y

   aha = 0
   ! find the depth at which the aggregates are
   CALL find_depth(pa%z, i)

   ! find the sizebin of the aggregates
   aha = 0
   l = 1
   DO WHILE (l .lt. n_size-1 .and. aha .eq. 0)
      IF (pa%r .lt. radi(l)) THEN
         aha = 1
      ELSE
         l = l + 1
      ENDIF
   ENDDO

   m = 1
   organic = 0
   TEP = 0
   DO WHILE (m .le. 10) 
     organic = organic + pa%orgC(m,1)*pa%n
     TEP = TEP + pa%TEP(m,1)*pa%n
     m = m + 1
   ENDDO
   
   calc = pa%mineral(1)*pa%n + pa%mineral(2)*pa%n
   IF (calc .gt. 0) THEN
      rainratio = organic/calc
   ENDIF

   ! check if the indexes are within the array
#ifdef array_size
   max_i = UBOUND(data_mega,2)
   max_l = UBOUND(data_mega,1)
   CALL check_index(i, max_i, x)
   CALL check_index(l, max_l, y)
   IF (x .gt. 0 .or. y .gt. 0) THEN
      print*, 'OVERWRITE ARRAY in collect2'
      print*, i, max_i, l, max_l
      print*, pa%id, pa%z, pa%r
   ENDIF
#endif

   data_mega(l,i,1,1) = data_mega(l,i,1,1) + pa%rho*pa%n
   data_mega(l,i,2,1) = data_mega(l,i,2,1) + pa%n
   data_mega(l,i,1,2) = data_mega(l,i,1,2) + pa%w*pa%n
   data_mega(l,i,2,2) = data_mega(l,i,2,2) + pa%n
   data_mega(l,i,1,3) = data_mega(l,i,1,3) + pa%s*pa%n
   data_mega(l,i,2,3) = data_mega(l,i,2,3) + pa%n
   data_mega(l,i,1,4) = data_mega(l,i,1,4) + organic
   data_mega(l,i,1,5) = data_mega(l,i,1,5) + calc
   data_mega(l,i,1,6) = data_mega(l,i,1,6) + pa%mineral(3)*pa%n
   data_mega(l,i,1,7) = data_mega(l,i,1,7) + pa%mineral(4)*pa%n
   data_mega(l,i,1,8) = data_mega(l,i,1,8) + rainratio * pa%n
   data_mega(l,i,2,8) = data_mega(l,i,2,8) + pa%n
   data_mega(l,i,1,9) = data_mega(l,i,1,9) + pa%por*pa%n
   data_mega(l,i,2,9) = data_mega(l,i,2,9) + pa%n
   data_mega(l,i,1,10) = data_mega(l,i,1,10) + TEP

END SUBROUTINE collect2
!=========================================================================
SUBROUTINE sediment_trap(p, z)
! this subroutine collects the flux that passes a certain depth
   USE the_info

   TYPE(agg_part), INTENT(IN) :: p
   DOUBLE PRECISION :: TEP
   INTEGER, INTENT(IN) :: z
   INTEGER :: i

   i = 1
   TEP = 0
   DO WHILE (i .le. 10)
      sedTrap(z,i) = sedTrap(z,i) + p%orgC(i,1)*p%n 
      TEP = TEP + p%TEP(i,1)
      i = i + 1
   ENDDO

   sedTrap(z,11) = sedTrap(z,11) + p%mineral(1)*p%n + p%mineral(2)*p%n
   sedTrap(z,12) = sedTrap(z,12) + p%mineral(3)*p%n
   sedTrap(z,13) = sedTrap(z,13) + p%mineral(4)*p%n
   sedTrap(z,14) = sedTrap(z,14) + TEP*p%n
   
   IF (p%mineral(1) .lt. 0 .or. p%mineral(3) .lt. 0) THEN
      print*, 'sedTrap:', p%mineral(1), p%mineral(2)
   ENDIF
   IF (sedTrap(z,11) .lt. 0) print*, 'sedTrap problem', &
      p%Nn, p%n, p%id

END SUBROUTINE sediment_trap
!=========================================================================
SUBROUTINE open_files()
   USE the_info

   OPEN (UNIT=12, FILE="orgC", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=13, FILE="calc", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=14, FILE="opal", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=15, FILE="clay", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=16, FILE="rratio", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=18, FILE="seafloor", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=19, FILE="primaryProd", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=110, FILE="seafloor1yr", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=20, FILE="density.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=21, FILE="velocity.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=22, FILE="stickiness.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=23, FILE="organic.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=24, FILE="inorganic.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=25, FILE="opal.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=26, FILE="clay.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=27, FILE="rainratio.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=55, FILE="TEP.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=46, FILE="porosity.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=28, FILE="orgFloor", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=29, FILE="mineralFloor", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=31, FILE="ratioFloor", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=34, FILE="orgFlux1km", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=37, FILE="calcFlux1km", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=40, FILE="opalFlux1km", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=43, FILE="ratioFlux1km", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=45, FILE="fluxSeafloor", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=48, FILE="one.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=49, FILE="two.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=50, FILE="three.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=51, FILE="four.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=104, FILE="seaFloor.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=111, FILE="veloFloor.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=112, FILE="densFloor.out", STATUS='NEW', POSITION='APPEND')
   OPEN (UNIT=113, FILE="scaling.out", STATUS='NEW', POSITION='APPEND')

END SUBROUTINE open_files   
!=========================================================================
SUBROUTINE write_zooplankton(time, mic)
   USE the_info
   
   101 FORMAT (E, 1X, I5)

   INTEGER, INTENT(IN) :: time
   DOUBLE PRECISION, DIMENSION(n_db), INTENT(IN) :: mic
   INTEGER :: i, UNIT
   CHARACTER (len=30) :: file_name
   CHARACTER (len=14) :: ctime

print*, 'write zooplankton', file_name
! write out the zooplankton array
      WRITE(ctime,*) time
      ctime = adjustl(ctime)
      WRITE(file_name,"(a30)")(trim(ctime)//".zoo")
      OPEN(UNIT=56, FILE=file_name, STATUS='NEW', POSITION='APPEND')
      i = 1
      DO WHILE (i .le. n_db) 
         WRITE(56,101) mic(i), i
         i = i + 1
      ENDDO
      print*, ' zoop'

END SUBROUTINE write_zooplankton
!=========================================================================
SUBROUTINE write_npp(time, agg_max)
   USE the_info
   
   119 FORMAT (10F18.6)

   INTEGER, INTENT(IN) :: time, agg_max

   print*, '--------------', agg_max, '--------'
   print*, 'npp:', npp, 'sf:', sf, 'time:', time
   print*, '-----------------*-----------------'

   WRITE(110,119) sf(1), sf(2), sf(3), sf(4), sf(5), npp(1), &
                  npp(2), npp(3), npp(4), npp(5)
   sf(:) = 0
   npp(:) = 0

END SUBROUTINE write_npp
!=========================================================================
SUBROUTINE write_sedTrap(time, day)
   USE the_info
   
   200 FORMAT (13E, 1X, I5, 1X, I5, 1X, 3E)

   INTEGER, INTENT(IN) :: time
   DOUBLE PRECISION, INTENT(IN) :: day
   INTEGER :: i, j, UNIT
   CHARACTER (len=30) :: file_name
   CHARACTER (len=14) :: ctime

! parameters for file-names
            WRITE(ctime,*) time
            ctime = adjustl(ctime)

! open the files to write flux into 
            WRITE(file_name,"(a14)")("trap"//trim(ctime))
            OPEN(UNIT=53, FILE=file_name, STATUS='NEW', POSITION='APPEND')

! we average over a year, and return the flux in [molC/m2d]
            i = 1
            DO WHILE (i .le. 9)
               WRITE(53,200) sedTrap(i,1)/day, sedTrap(i,2)/day, sedTrap(i,3)/day,&
               sedTrap(i,4)/day, sedTrap(i,5)/day, sedTrap(i,6)/day, sedTrap(i,7)/day,&
               sedTrap(i,8)/day, sedTrap(i,9)/day, sedTrap(i,10)/day, sedTrap(i,11)/day,&
               sedTrap(i,12)/day, sedTrap(i,13)/day, time*timestep/86400, i*10, &
               (sedTrap(i,1)+sedTrap(i,2)+sedTrap(i,3)+sedTrap(i,4)+sedTrap(i,5)+&
               sedTrap(i,6)+sedTrap(i,7)+sedTrap(i,8)+sedTrap(i,9)+sedTrap(i,10))/day, &
               sedTrap(i,14)/day
               i = i + 1
            ENDDO
print*, 'sedTrap', i, time
            DO WHILE (i .le. 20)
               WRITE(53,200) sedTrap(i,1)/day, sedTrap(i,2)/day, sedTrap(i,3)/day,&
               sedTrap(i,4)/day, sedTrap(i,5)/day, sedTrap(i,6)/day, sedTrap(i,7)/day,&
               sedTrap(i,8)/day, sedTrap(i,9)/day, sedTrap(i,10)/day, sedTrap(i,11)/day,&
               sedTrap(i,12)/day, sedTrap(i,13)/day, time*timestep/86400, (i-9)*100, &
               (sedTrap(i,1)+sedTrap(i,2)+sedTrap(i,3)+sedTrap(i,4)+sedTrap(i,5)+&
               sedTrap(i,6)+sedTrap(i,7)+sedTrap(i,8)+sedTrap(i,9)+sedTrap(i,10))/day, &
               sedTrap(i,14)/day
               i = i + 1
            ENDDO
print*, 'sedTrap', i, time, j
            j = 15
            DO WHILE (i .le. 26)
               WRITE(53,200) sedTrap(i,1)/day, sedTrap(i,2)/day, sedTrap(i,3)/day,&
               sedTrap(i,4)/day, sedTrap(i,5)/day, sedTrap(i,6)/day, sedTrap(i,7)/day,&
               sedTrap(i,8)/day, sedTrap(i,9)/day, sedTrap(i,10)/day, sedTrap(i,11)/day,&
               sedTrap(i,12)/day, sedTrap(i,13)/day, time*timestep/86400, j*100, &
               (sedTrap(i,1)+sedTrap(i,2)+sedTrap(i,3)+sedTrap(i,4)+sedTrap(i,5)+&
               sedTrap(i,6)+sedTrap(i,7)+sedTrap(i,8)+sedTrap(i,9)+sedTrap(i,10))/day, &
               sedTrap(i,14)/day
               i = i + 1
               j = j + 5
            ENDDO
            sedTrap(:,:) = 0

END SUBROUTINE write_sedTrap
!=========================================================================
SUBROUTINE write_propPic(data_mega)
   USE the_info

   120 FORMAT (E15.6E2, 3X, F14.4, 1X, I)

   DOUBLE PRECISION, DIMENSION(n_size,n_db,2,10), INTENT(INOUT) :: data_mega
   INTEGER :: i, j, UNIT



         i = 1
         j = 1
         DO WHILE (j .le. n_db)
            DO WHILE (i .lt. n_size)
               IF (data_mega(i,j,1,1) .gt. 0) THEN
                  data_mega(i,j,1,1) = data_mega(i,j,1,1)/data_mega(i,j,2,1)
               ENDIF
               WRITE(20,120) data_mega(i,j,1,1), radi(i), j !density

               IF (data_mega(i,j,1,2) .gt. 0) THEN
                  data_mega(i,j,1,2) = data_mega(i,j,1,2)/data_mega(i,j,2,2)
               ENDIF
               IF (data_mega(i,j,2,2) .le. 0) THEN
                  data_mega(i,j,1,2) = 0.0
               ENDIF
               WRITE(21,120) data_mega(i,j,1,2)*864, radi(i), j !velocity

               IF (data_mega(i,j,1,3) .gt. 0) THEN
                  data_mega(i,j,1,3) = data_mega(i,j,1,3)/data_mega(i,j,2,3)
               ENDIF
               IF (data_mega(i,j,2,3) .le. 0) THEN
                  data_mega(i,j,1,3) = 0
               ENDIF
               WRITE(22,120) data_mega(i,j,1,3), radi(i), j   !stickiness

               data_mega(i,j,1,4) = data_mega(i,j,1,4)!/(radius(i+1)-radius(i))
               WRITE(23,120) data_mega(i,j,1,4), radi(i), j       ! organic

               data_mega(i,j,1,5) = data_mega(i,j,1,5)!/(radius(i+1)-radius(i))
               WRITE(24,120) data_mega(i,j,1,5), radi(i), j       ! inorganic

               data_mega(i,j,1,6) = data_mega(i,j,1,6)!/(radius(i+1)-radius(i))
               WRITE(25,120) data_mega(i,j,1,6), radi(i), j       ! opal

               data_mega(i,j,1,7) = data_mega(i,j,1,7)!/(radius(i+1)-radius(i))
               WRITE(26,120) data_mega(i,j,1,7), radi(i), j       ! clay

               IF (data_mega(i,j,1,8) .gt. 0) THEN   ! rain ratio
                  data_mega(i,j,1,8) = data_mega(i,j,1,8)/data_mega(i,j,2,8)
               ENDIF
               WRITE(27,120) data_mega(i,j,1,8), radi(i), j

               IF (data_mega(i,j,1,9) .gt. 0) THEN   ! porosity
                  data_mega(i,j,1,9) = data_mega(i,j,1,9)/data_mega(i,j,2,9)
               ENDIF
               WRITE(46,120) data_mega(i,j,1,9), radi(i), j

               data_mega(i,j,1,10) = data_mega(i,j,1,10)!/(radius(i+1)-radius(i))
               WRITE(55,120) data_mega(i,j,1,10), radi(i), j       ! TEP

               i = i + 1
            ENDDO
            j = j + 1
            i = 1
         ENDDO



END SUBROUTINE write_propPic
!=========================================================================
SUBROUTINE write_seaFloor(seaBed, oneAgg)
   USE the_info

   122 FORMAT (2(E))

   DOUBLE PRECISION, INTENT(IN) :: seaBed(:)
   TYPE(agg_part), INTENT(IN) :: oneAgg

         i = 1
         DO WHILE (i .le. 10)
            IF (i .eq. 1) THEN
               WRITE(28,122) seaBed(i), seaBed(i)/oneAgg%orgC(i,2)
            ELSEIF (i .gt. 1) THEN
               WRITE(28,122) seaBed(i), seaBed(i)/(oneAgg%orgC(i,2)-oneAgg%orgC(i-1,2))
            ENDIF
            i = i + 1
         ENDDO

         DO WHILE (i .le. 19)
            WRITE(29,122) seaBed(i)
            i = i + 1
         ENDDO

         WRITE(31,122) (seaBed(1)+seaBed(2)+seaBed(3)+seaBed(4)+seaBed(5)+&
                       seaBed(6)+seaBed(7)+seaBed(8)+seaBed(9)+seaBed(10))/&
                       (seaBed(11)+seaBed(12))

END SUBROUTINE write_seaFloor
!=========================================================================
SUBROUTINE write_specPic(spec, time)
   USE the_info

   102 FORMAT (E15.4E2, 1X, E15.4E2, 1X, I)

!   DOUBLE PRECISION, INTENT(IN) :: spec(:,:)
   DOUBLE PRECISION, DIMENSION(n_size,n_db), INTENT(IN) :: spec
   INTEGER, INTENT(IN) :: time
   INTEGER :: i, j, UNIT
   CHARACTER (len=30) :: file_name
   CHARACTER (len=14) :: ctime

   WRITE(ctime,*) time
   ctime = adjustl(ctime)
print*, ctime
   WRITE(file_name,"(a12)")("s"//trim(ctime)//".dat")
   OPEN(UNIT=17, FILE=file_name, STATUS='NEW', POSITION='APPEND')

   j = 1
   i = 1
   DO WHILE (j .le. n_db)
      DO WHILE (i .lt. n_size)
         IF (i .eq. 1) THEN
            WRITE(17,102) spec(1,j)/radi(1), radi(i), j*dz/100
         ELSEIF (i .gt. 1) THEN
            WRITE(17,102) spec(i,j)/(radi(i+1)-radi(i)), radi(i),&
                          j*dz/100
         ENDIF
         i = i + 1
      ENDDO
      j = j + 1
      i = 1
   ENDDO

END SUBROUTINE write_specPic
!=========================================================================
SUBROUTINE write_dataFlux(pflux, time)
   USE the_info

   121 FORMAT (5(E, 1X))


   DOUBLE PRECISION, INTENT(IN) :: pflux(4)
   INTEGER, INTENT(IN) :: time
   INTEGER :: UNIT

   WRITE(34,121) pflux, time*timestep/86400.0
   WRITE(37,121) pflux, time*timestep/86400.0
   WRITE(40,121) pflux, time*timestep/86400.0
   WRITE(43,121) pflux/pflux, time*timestep/86400.0

END SUBROUTINE write_dataFlux
!=========================================================================
SUBROUTINE write_singlePart(pa, time, u)
   USE the_info

   147 FORMAT (5(F,1X),13(E14.4E3, 1X))

   TYPE(agg_part), INTENT(IN) :: pa
   INTEGER, INTENT(IN) :: time, u 

   CALL find_depth(pa%z, i)

   WRITE(u,147) time*timestep/86400.0, pa%r, pa%rho, pa%w*864, &
                pa%s, pa%mineral(1)+pa%mineral(2), pa%mineral(3), &
                pa%mineral(4),pa%orgC(1,1)+pa%orgC(2,1)+pa%orgC(3,1), &
                pa%orgC(4,1)+pa%orgC(5,1)+pa%orgC(6,1), &
                pa%orgC(7,1)+pa%orgC(8,1)+pa%orgC(9,1)+pa%orgC(10,1), &
                pa%z, pa%Nn, pa%pr, pa%frac, pa%n, pa%por, &
                pa%TEP(1,1)+pa%TEP(2,1)+pa%TEP(3,1)+pa%TEP(4,1)+&
                pa%TEP(5,1)+pa%TEP(6,1)+pa%TEP(7,1)+pa%TEP(8,1)+&
                pa%TEP(9,1)+pa%TEP(10,1)

END SUBROUTINE write_singlePart
!=========================================================================
SUBROUTINE eotAverage(dataArray, orgCM, calcM, opalM, clayM, rainM)
   USE the_info

   147 FORMAT (5(F,1X),13(E14.4E3, 1X))

   DOUBLE PRECISION, DIMENSION(15,n_db), INTENT(INOUT) :: dataArray
   DOUBLE PRECISION, DIMENSION(n_db), INTENT(INOUT) :: orgCM, calcM, rainM
   DOUBLE PRECISION, DIMENSION(n_db), INTENT(INOUT):: opalM, clayM
   DOUBLE PRECISION, DIMENSION(n_db) :: orgCtotal, calctotal
   INTEGER :: i

   DO i = 1, n_db, 1
      orgCtotal(i) = dataArray(2,i)+dataArray(3,i)+dataArray(4,i)+&
                     dataArray(5,i)+dataArray(6,i)+dataArray(7,i)+&
                     dataArray(8,i)+dataArray(9,i)+dataArray(10,i)+&
                     dataArray(11,i)
      calctotal(i) = dataArray(12,i)+dataArray(13,i)
      orgCM(i) = orgCM(i) + orgCtotal(i)
      calcM(i) = calcM(i) + calctotal(i)
      opalM(i) = opalM(i) + dataArray(14,i)
      clayM(i) = clayM(i) + dataArray(15,i)
      rainM(i) = orgCM(i)/calcM(i) !that is averaged over to get orgCMean etc.`
   END DO

END SUBROUTINE eotAverage
!=========================================================================
SUBROUTINE write_eotA(orgCM, calcM, opalM, clayM, rainM, tsV)
   USE the_info

   102 FORMAT (E15.4E2, 1X, I)

   DOUBLE PRECISION, DIMENSION(n_db), INTENT(IN) :: orgCM, calcM, rainM
   DOUBLE PRECISION, DIMENSION(n_db), INTENT(IN):: opalM, clayM
   INTEGER, INTENT(IN) :: tsV
   INTEGER :: i

   DO i = 1, n_db, 1
      WRITE(12,102) orgCM(i)/tsV, i*dz
      WRITE(13,102) calcM(i)/tsV, i*dz
      WRITE(14,102) opalM(i)/tsV, i*dz
      WRITE(15,102) clayM(i)/tsV, i*dz
      WRITE(16,102) rainM(i), i*dz
   END DO

END SUBROUTINE write_eotA
!=========================================================================
SUBROUTINE write_sf(time)
   USE the_info

   145 FORMAT (E15.4E2, 1X, E15.4E2, 1X, E15.4E2, 1X, I4)

   INTEGER, INTENT(IN) :: time

   WRITE(45, 145) sf(1)/(time*timestep/86400), &
                  sf(2)/(time*timestep/86400), &
                  sf(1)/sf(2), time*timestep/86400


END SUBROUTINE write_sf
!=========================================================================
SUBROUTINE write_end()
   USE the_info

   104 FORMAT (E15.4E2, 1X, E15.4E2, 1X, E15.4E2)
   118 FORMAT (5(F12.6, 1X))
   245 FORMAT (6E)
   260 FORMAT (I5, 1X, 6E)

   INTEGER :: i, UNIT
   CHARACTER (len=30) :: file_name

print*, 'write_end'
   WRITE(18,118) sf(1), sf(2), sf(3), sf(4), sf(5)
   WRITE(19,118) npp(1), npp(2), npp(3), npp(4), npp(5)
   i = 1
   DO WHILE (i .le. n_size)
      IF (i .eq. 1) THEN
         WRITE(104,104) sizeBed(i), radi(i), sizeBed(i)/radi(i)
         WRITE(111,104) veloBed(i), radi(i), veloBed(i)/(radi(i)/864)
         WRITE(112,104) densBed(i), 0.8+0.03*(i-1), densBed(i)/0.03
         WRITE(113,245) scaling(1,i), scaling(2,i)/radi(i), radi(i)
      ELSE
         WRITE(104,104) sizeBed(i), radi(i), sizeBed(i)/(radi(i)-radi(i-1))
         WRITE(111,104) veloBed(i), radi(i), veloBed(i)/((radi(i)-radi(i-1))/864)
         WRITE(112,104) densBed(i), 0.8+0.03*(i-1), densBed(i)/0.03
         WRITE(113,245) scaling(1,i), scaling(2,i)/(radi(i)-radi(i-1)), radi(i)
      ENDIF
      i = i + 1
   ENDDO

print*, 'write_end', i
   WRITE(file_name,"(a30)")("consu")
   OPEN(UNIT=210, FILE=file_name, STATUS='NEW', POSITION='APPEND')
   i = 1
   DO WHILE (i .le. n_db)
      orgC_G1(i) = orgC_G1(i)/(dz/100)
      orgC_b(i)  = orgC_b(i)/(dz/100)
      calc_G1(i) = calc_G1(i)/(dz/100)
      calcC_t(i) = calcC_t(i)/(dz/100)
      calcA_t(i) = calcA_t(i)/(dz/100)
      opal_t(i)  = opal_t(i)/(dz/100)
      WRITE(210,260) i*dz/100, orgC_G1(i), orgC_b(i), &
      calc_G1(i), calcC_t(i), calcA_t(i), opal_t(i)
      i = i + 1
   ENDDO



END SUBROUTINE write_end
!=========================================================================
