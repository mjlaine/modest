      subroutine dscolr(ar,dar,ai,dai,
     &                  nrx, ncx, x, ldx, icomp, iordr, iret,
     &                  nkey, indkey, iperm, ngroup, ni)

*     ARGUMENTS
      integer*4 dar,dai
      real*8    ar(dar)
      integer*4 ai(dai)


      integer nrx, ncx, ldx, icomp, iordr, iret, nkey, ngroup
      double precision x(ldx, ncx)
      integer indkey(nkey), iperm(ncx), ni(ncx)

      integer i
*      double precision y(maxnopt + 1, maxnopt + 1)
      integer eqcol
      external dcopy

      include 'common3.inc'

      do 10040  i = 1, ncx
          iperm(i) = i
          ni(i) = 0
10040 continue

      call qcksrt(nrx, ncx, x, ldx, icomp, iordr,
     &            nkey, indkey, iperm)

      ngroup = 1
      ni(1) = 1
      do 10020  i = 1, ncx - 1
          if (eqcol(x(1, iperm(i)), x(1, iperm(i + 1)), nrx, icomp,
     &              nkey, indkey).eq.0)
     &        ngroup = ngroup + 1
          ni(ngroup) = ni(ngroup) + 1
10020 continue

      if (iret.eq.0) then
          do 10000  i = 1, ncx
*              call dcopy(nrx, x(1, iperm(i)), 1, y(1, i), 1)
              call dcopy(nrx,x(1,iperm(i)),1,
     &                   ar(pwrkdsco+(i-1)*(dnesopsi+1)),1)
10000     continue
          do 10010  i = 1, ncx
*              call dcopy(nrx, y(1, i), 1, x(1, i), 1)
              call dcopy(nrx,ar(pwrkdsco+(i-1)*(dnesopsi+1)),1,
     &                   x(1,i),1)
10010     continue
      end if

      return
      end


      integer function lecol(y1, y2, n, icomp, iordr, nkey, indkey)

      integer n, icomp, iordr, nkey
      double precision y1(n), y2(n)
      integer indkey(nkey)

      integer j
      double precision difference

      lecol = 1
      j = 1
10050 if (j.le.nkey) then
          if (icomp.eq.1) then
              difference = dabs(y1(indkey(j))) - dabs(y2(indkey(j)))
          else
              difference = y1(indkey(j)) - y2(indkey(j))
          end if
          if (iordr.eq.1) then
              if (difference.gt.0) then
                  return
              else if (difference.lt.0) then
                  lecol = 0
                  return
              end if
          else
              if (difference.lt.0) then
                  return
              else if (difference.gt.0) then
                  lecol = 0
                  return
              end if
          end if
          j = j + 1
          goto 10050
      end if

      return
      end


      integer function eqcol(y1, y2, n, icomp, nkey, indkey)

      integer n, icomp, nkey
      double precision y1(n), y2(n)
      integer indkey(nkey)

      integer j

      eqcol = 1
      j = 1
10050 if (j.le.nkey) then
          if (icomp.eq.1) then
              if (dabs(y1(indkey(j))).ne.dabs(y2(indkey(j)))) then
                  eqcol = 0
                  return
              end if
          else
              if (y1(indkey(j)).ne.y2(indkey(j))) then
                  eqcol = 0
                  return
              end if
          end if
          j = j + 1
          goto 10050
      end if

      return
      end


      subroutine qcksrt(nrx, ncx, x, ldx, icomp, iordr,
     &                  nkey, indkey, iperm)

      integer nrx, ncx, ldx, icomp, iordr, nkey
      double precision x(ldx, ncx)
      integer indkey(nkey), iperm(ncx)

      integer lecol
      integer ia

      parameter (m=7,nstack=50,fm=7875.,fa=211.,fc=1663.
     *    ,fmi=1.2698413e-4)
      dimension istack(nstack)
      jstack=0
      l=1
      ir=ncx
      fx=0.
10    if(ir-l.lt.m)then
        do 13 j=l+1,ir
          ia = iperm(j)
          do 11 i=j-1,1,-1
            if (lecol(x(1, iperm(i)), x(1, ia), nrx, icomp, iordr,
     &                nkey, indkey).eq.1) goto 12
            iperm(i+1)=iperm(i)
11        continue
          i=0
12        iperm(i+1)=ia
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        i=l
        j=ir
        fx=mod(fx*fa+fc,fm)
        iq=l+(ir-l+1)*(fx*fmi)
        ia=iperm(iq)
        iperm(iq)=iperm(l)
20      continue
21        if (j.gt.0) then
              if (lecol(x(1, iperm(j)), x(1, ia), nrx, icomp, iordr,
     &                  nkey, indkey).eq.0) then
                  j=j-1
                  goto 21
              end if
          end if
          if(j.le.i)then
            iperm(i)=ia
            goto 30
          end if
          iperm(i)=iperm(j)
          i=i+1
22        if (i.le.ncx) then
              if (lecol(x(1, ia), x(1, iperm(i)), nrx, icomp, iordr,
     &                  nkey, indkey).eq.0) then
                  i=i+1
                  goto 22
              end if
          end if
          if(j.le.i)then
            iperm(j)=ia
            i=j
            goto 30
          end if
          iperm(j)=iperm(i)
          j=j-1
        goto 20
30      jstack=jstack+2
        if(jstack.gt.nstack)pause 'nstack must be made larger.'
        if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
        else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
        end if
      end if
      goto 10
      end

