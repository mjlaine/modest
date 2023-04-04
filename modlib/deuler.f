
          subroutine deuler(nstates,t0aux,s,Fodeaux,tj,dstp,workeul)
          integer*4  nstates
          real*8     s(nstates),t0,tj,dstp
          real*8     dstpaux,taux,t0aux
          real*8     workeul(2*nstates)
          real*8     k1,k2,k3,k4,t1
          external   Fodeaux

            taux = t0aux

           do while(taux. lt. tj)
             call Fodeaux(nstates, taux, s, workeul)
             do i = 1,nstates
                s(i) = s(i) + dstp*workeul(i)
             enddo
             t0aux = taux
             taux  = taux + dstp
           enddo

           dstpaux = tj - t0aux
           call Fodeaux(nstates, tj, s, workeul)
           do i = 1,nstates
              s(i) = s(i) + dstpaux*workeul(i)
           enddo

c           do while(taux. lt. tj)
c             call Fodeaux(nstates, taux, s, workeul)
c             do i = 1,nstates
c                workeul(nstates+i) = s(i)
c                k1   = s(i) + dstp*workeul(i)
c                s(i) = s(i) + 0.5d0*k1
c             enddo
c             t1 = taux + 0.5d0*dstp
c             call Fodeaux(nstates,t1,s,workeul)
c             do i = 1,nstates
c                s(i) = workeul(nstates+i) + dstp*workeul(i)
c             enddo
c
c             t0aux = taux
c             taux  = taux + dstp
c            enddo
c
c             dstpaux = tj - t0aux
c             call Fodeaux(nstates, taux, s, workeul)
c             do i = 1,nstates
c                workeul(nstates+i) = s(i)
c                k1   = s(i) + dstpaux*workeul(i)
c                s(i) = s(i) + 0.5d0*k1
c             enddo
c             t1 = taux + 0.5d0*dstpaux
c             call Fodeaux(nstates,t1,s,workeul)
c             do i = 1,nstates
c                s(i) = workeul(nstates+i) + dstpaux*workeul(i)
c             enddo

          return
          end

