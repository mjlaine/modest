      subroutine dlsvrr(nra, nca, a, lda, ipath, info, tol, irank,
     &    s, u, ldu, v, ldv, wrk, work, e)

      implicit none
      integer nra, nca, lda, ipath, info, irank, ldu, ldv
      real*8  a(lda, nca), tol, s(lda), u(ldu, nra), v(ldv, nca),
     &        wrk(lda,nca), work(lda), e(nca)

      integer i, j, idim
      real*8  abstol, anorm, sumrow
      external dcopy, dsvdc

      do 10 i = 1, nca
          call dcopy(nra, a(1, i), 1, wrk(1, i), 1)
 10   continue

      call dsvdc(wrk, lda, nra, nca, s, e, u, ldu, v, ldv,
     &    work, ipath, info)

      idim = nra + 1
      if (nca.lt.idim) idim = nca
      irank = 0

      if (tol.lt.0) then
          anorm = 0
          do 20 j = 1, nra
              sumrow = 0
              do 30 i = 1, nca
                  sumrow = sumrow + dabs(a(j, i))
 30           continue
              if (anorm.lt.sumrow) anorm = sumrow
 20       continue
          abstol = dabs(tol)*anorm
      else
          abstol = tol
      end if

      do 40 i = 1, idim
          if (s(i).gt.abstol) irank = irank + 1
 40   continue

      return
      end

