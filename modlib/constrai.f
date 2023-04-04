      subroutine Constraints(xopt,nopt,ibound,bounds,inv)

*     CALLED BY

*     - the chosen objective function (e.g. LSQ)
*     - MAIN

*     SUBROUTINES NEEDED

*     None

*     PURPOSE

*     If inv equals 1 this routine transforms the free variables
*     to bounded ones according to ibound & bounds.  If inv equals -1
*     an inverse transformation is done. This is needed for calculation
*     of the initial values of the free variables - the user gives
*     initial values in units that have physical meaning i.e. initial
*     values for the bounded variables.


*     ARGUMENTS

      implicit none
      integer*4 nopt                              ! n of parameters and
      real*8    xopt(nopt)                        ! parameters to be optimized
      integer*4 ibound(nopt)                      ! see common_blocks.inc
      real*8    bounds(2,nopt)                    ! bounds for parameters
      integer*4 inv                               ! see PURPOSE

*     LOCAL VARIABLES

      integer*4 lohi                              ! ibound
      integer*4 i,j                               ! loop indices
      real*8    a,b,aplusb,bminusa                ! lower and upper bounds

      if( inv .eq. -1) then

          do i=1,nopt

              lohi = ibound(i)                    ! -1 for lower bounds and
                                                  ! +1 for upper bounds
              j = 1.5 + lohi*0.5                  ! j = 1 or 2 for bounds

              if(abs(lohi) .eq. 1) then

                 xopt(i) = dsqrt((-lohi)*(xopt(i) - bounds(j,i)))

              elseif(lohi .eq. 2) then

                 a       = bounds(1,i)
                 b       = bounds(2,i)
                 aplusb  = (a+b)/2.0d0
                 bminusa = (b-a)/2.0d0
                 a       = (xopt(i) - aplusb)/bminusa
                 if (a.gt.1.0d0) a = 1.0d0
                 if (a.lt.-1.0d0) a = -1.0d0
                 xopt(i) = dasin(a)

              elseif(lohi .eq. 0) then

                 continue

              else
                 write(*,*)
     &                'type of constraint must be -1,0,1 or 2!'
                 stop

              endif

          enddo

      elseif( inv .eq. 1) then

          do i=1,nopt

              lohi = ibound(i)                    ! -1 for lower bounds and
                                                  ! +1 for upper bounds
              j = 1.5 + lohi*0.5                  ! j = 1 or 2 for bounds

              if(abs(lohi) .eq. 1) then

                 xopt(i) = (-lohi)*( xopt(i) )**2 + bounds(j,i)

              elseif(lohi .eq. 2) then

                 a = bounds(1,i)
                 b = bounds(2,i)
                 a       = bounds(1,i)
                 b       = bounds(2,i)
                 aplusb  = (a+b)/2.0d0
                 bminusa = (b-a)/2.0d0
                 xopt(i) = aplusb + bminusa*dsin(xopt(i))

              elseif(lohi .eq. 0) then

                 continue

              else
                 write(*,*)
     &                'type of constraint must be -1,0,1 or 2!'
                 stop

              endif

          enddo

      else

          write(*,*) 'inv must be -1 or 1!'
          stop

      endif

      return

      end

