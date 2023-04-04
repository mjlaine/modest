c ===================================================================
c smbl.f:
c 23.10.94
c 24.10.94 - makeincfile-rutiini osin uusiksi
c 25.10.94 - tmp-tiedostot deletoidaan
c 25.10.94 - stringit pieniksi ...
c 08.03.95 - boundsien tarkistus ... (KESKENER[INEN) 
c 18.04.95 - Puolipiste modelvar-stringeihin ts. kopiointimerkki.
c 24.04.95 - Korjaus (lis{ys rajojen tarkistukseen bounds_check)
c  5.5.95  - korjattu cmpfiles-rutiini
c 21.6.95  - Korjattu tarkistusta rutiinissa uppd_lb
c            (nlpar -> mxnlpar,nstates0,mxnstat)
c 22.11.95 - bounds: sbounds comm. /HH
c 22.11.95 - discfile /HH
c ===================================================================
       subroutine smbl(task,model,modelvar,target,estvar,nsets,   
     *                 gpar,lpar,states0,s0file,lsim,lest,lexp,lsen,
     *                 lopt,bound,nbound,ngpar,nlpar,nstates0,nsim,nest,
     *                 nexp,nsen,nopt,nmlfile,error_number,mxnstat,
     *                 mxnlpar,mxchars,mxwords,mxmsent,mxnlaux,type, 
     *                 charcl,nglobal,nlocal,ncontrol,ninitval,laux,  
     *                 word,sentence,name,ncodir,codir,
     *                 evalfile,discfile,nfilepoint,ndiscpoint,
     *                 ngrid,igrid,codirvar,gridvar,maxnobs,
     *                 t0,tf,ntf,vardim,nx,nsbound,sbound,optcrit,
     *                 objfun,projectname,mustcmp,ftype,inttype,iperiod)

!   task          ! taski (exp,sim ...)          input
!   model         ! model (ode,...)              input
!   modelvar      ! input-str for modelvar-part  input
!   target        ! input-str for target-part    input
!   estvar        ! input-str for estvar-part    input
!   nsets         ! number of datasets           input
!   gpar          ! array of global vars         output
!   lpar          ! array of local vars          output
!   states0       ! array of initval vars        output
!   s0file        ! array of 0-1 vars            output
!   lsim          ! array of simul. vars         output
!   lest          ! array of est. vars           output
!   lexp          ! array of exp. vars           output
!   lsen          ! array of sens. vars          output
!   lopt          ! array of opt vars            output
!   bound         ! array of bounds              output
!   nbound        ! number of bounds             output
!   ngpar         ! number of global vars        output
!   nlpar         ! number of local vars         output
!   nstates0      ! number of initval. vars      output
!   nsim          ! number of sim. vars          output
!   nest          ! number of est. vars          output
!   nexp          ! number of exp. vars          output
!   nsen          ! number of sens. vars         output
!   nopt          ! number of opt. vars          output
!   nmlfile       ! name of nml-file             input
!   error_number  ! errornumber                  input/output
!   mxnstat       ! max of initval. vars         input
!   mxnlpar       ! max of local vars            input
!   mxchars       ! max length of strs           input
!   mxwords       ! max number of words of strs  input
!   mxmsent       ! max number of sentences      input
!   mxnlaux       ! max number of target vars    input
!   type          ! array of types of vars       work
!   charcl        ! array of types of chars      work
!   nglobal       ! ind. of glob. vars           work                         
!   nlocal        ! ind. of local vars           work
!   ncontrol      ! ind. of contr. vars          work
!   ninitval      ! ind. of initv. vars          work
!   laux          ! array of target-part vars    work
!   word          ! array of words of str        work
!   sentence      ! array of words of one sent   work
!   name          ! array of names of vars       work     
!   ncodir        ! number of codir-vars         output
!   codir         ! codir-array                  output
!   evalfile      ! filename with nfilepoints    input
!   discfile      ! filename with ndiscpoints    input
!   nfilepoint    ! number of filepoints         output
!   ndiscpoint    ! number of discpoints         output
!   ngrid         ! number of grid vars          output
!   codirvar      ! names of codir vars          input
!   gridvar       ! names of grid vars           input
!   maxnobs       ! max numb. of obs.            input
!   t0            ! array of t0-values           output
!   tf            ! array of tf-values           output
!   ntf           ! number of given tf-values    output
!   vardim        ! dimension of the variables   work
!   nx            ! number of control vars       output            
!   nsbound       ! number of sbounds            output      !!!! sbound
!   sbound        ! array of bounds              output      !!!! sbound
!   optcrit       ! criterion of the optim.      input       !!!! sbound
!   projectname   ! name of the project          input       
!   mustcmp       ! ind. of recompiling          output
!   ftype         ! the type of the fortran-fil. input
!   inttype       ! ind. of the integer-var.     output

       implicit none

       integer mxnstat,mxnlpar,mxchars,mxwords,mxmsent
c
c args:
       character*3 task,model
       character*(*) modelvar,target,estvar

       integer nsets
       double precision gpar(*),lpar(mxnlpar,*),
     *        states0(mxnstat,*)
       integer s0file(mxnstat,*),lsim(4,*),lest(4,*),
     *         lexp(4,*),lsen(4,*),lopt(4,*),nsim,nest,
     *         nexp,nsen,nopt,error_number,mxnlaux

       integer ncodir,nfilepoint,ndiscpoint,igrid(2,*)
       character*(*) codirvar,gridvar,evalfile,discfile
       integer codir(*),ngrid,maxnobs,ntf(*)
       double precision t0(*),tf(maxnobs,*)

       double precision bound(2,*)
       integer ngpar,nlpar,nstates0                        
       character*(*) nmlfile                                    
       integer charcl(0:mxchars+1)
       integer nglobal(mxmsent),nlocal(mxmsent),
     *         ncontrol(mxmsent),ninitval(mxmsent)
       integer laux(4,mxnlaux)
       character*(*) word(mxwords),sentence(mxwords)      
       character*(*) name(mxmsent)
       character*7 type(mxmsent)

       integer vardim(mxmsent,2)      
       integer nx                      

       integer nsbound
       real*8 sbound(2,*)
       character*3 optcrit,objfun

       integer nbound
       integer mustcmp
       character*3 ftype
       character*(*) projectname           
       integer inttype(*),iperiod

c
c locals:
       integer aascii(0:127)
       character*7 keyword
       character*8 part
       integer lenstr,nword,ithsentence,
     *         nsentence,nkey,nnumber,i,j,
     *         nmodelsentence,naux,
     *         lodevar/0/                                     !!!! 4.8.94 

       character*1 ch
       integer orignsentence !matrix-input
!25.10.94       character*3 loptcrit                                   !!!!!! 
       logical linttype

c ---------------------------------------------------------------
c
c initialize output-arguments:                                !!!! 16.9.94
         nbound       = 0
         ngpar        = 0
         nlpar        = 0
         nstates0     = 0
         nsim         = 0
         nest         = 0
         nexp         = 0
         nest         = 0
         nexp         = 0
         nsen         = 0
         nopt         = 0
         ncodir       = 0
         nfilepoint   = 0
         ngrid        = 0
         nsbound      = 0       
         error_number = 0
c
c 
!!!!
!!!! all (but not filenames) input-strings to lowercase: 
!!!! 25.10.94
       call case(task,len(task),1)       
       call case(model,len(model),1)
       call case(modelvar,len(modelvar),1)
       call case(target,len(target),1)
       call case(estvar,len(estvar),1)
       call case(evalfile,len(evalfile),1)
       call case(discfile,len(discfile),1)
       call case(codirvar,len(codirvar),1)
       call case(gridvar,len(gridvar),1)
       call case(optcrit,len(optcrit),1)
       call case(objfun,len(objfun),1)
       call case(projectname,len(projectname),1)
!!! 25.10.94      loptcrit = optcrit               !!!!!!
!!!       call case(loptcrit,3,1)          !!!!!! save optcrit with lower case.
c
c
       do i=1,mxmsent                   !!!! 4.8.94
         nglobal(i)  = 0                !!!! 4.8.94
         nlocal(i)   = 0                !!!! 4.8.94
         ncontrol(i) = 0                !!!! 4.8.94
         ninitval(i) = 0                !!!! 4.8.94
       enddo                            !!!! 4.8.94

       do i=1,mxmsent
         inttype(i) = 0
       enddo

c
c
         naux = 0

c define acceptable characters and write the
c aascii-array:
       call defchars(aascii)


************************modelvar-part:*******************************

       part = 'modelvar'
       lenstr  = len(modelvar)       
       ithsentence = 0
       nglobal(1)  = 0
       nlocal(1)   = 0
       ncontrol(1) = 0
       ninitval(1) = 0

c
c check the characters and write the 
c charcl-table: 
       call checkc(modelvar,lenstr,aascii,charcl,
     *      error_number)
       if(error_number.gt.0) then
         write(6,'(a)') ' ****************************************'
         write(6,'(a)') ' wrong character in the modelvar-string:'
         write(6,'(a)') modelvar(error_number:error_number)
         write(6,'(a)') ' ****************************************'
         stop
       endif

c
c low chars to upper:
       !25.10.94 call case(modelvar,lenstr,2)

c
c find the words from modelvar-string and write
c the word-array:
       call getwords(modelvar,lenstr,charcl,word,nword)


c loop://////////////////////////////////////////////////////////////
       dowhile(nword.gt.0)
           ithsentence = ithsentence + 1

           !find the next sentence:
           call next_sentence(word,nword,part,aascii,sentence,
     *          nsentence,keyword,nkey,error_number,linttype)

           if(error_number.gt.0) then
             write(6,'(a)') ' **************************************'
             write(6,'(a)') ' wrong syntax in the modelvar-string.'
             write(6,'(a)') ' **************************************'
             stop
           endif

           !***** matrix-input *************************************
           vardim(ithsentence,1) = 0
           vardim(ithsentence,2) = 0 
           orignsentence = nsentence
           if (sentence(6).eq.',') then
             call msentence(sentence,nsentence,ithsentence,mxmsent,
     *          vardim)        
             if(nkey.ne.11) then 
               write(6,'(a)') ' wrong syntax in the matrix-input'
               stop
             else
               nkey = nkey - 4  
             endif
           endif 
           !**** matrix-input****************************************

           !
           !update the nametable:
           call nametable(sentence,nsentence,ithsentence,
     *          keyword,nkey,name,type,nnumber,nglobal,nlocal,
     *          ncontrol,ninitval,error_number,linttype,
     *          inttype)

           if(error_number.gt.0) then
             write(6,'(a)') ' ****************************************'
             write(6,'(a)') ' wrong syntax in the modelvar-string.'
             write(6,'(a)') ' ****************************************'
             stop
           endif


           !update the gpar, lpar, states0 and s0file-array:       
           if(keyword.eq.'global') then
             call uppd_gpar(gpar,ngpar,sentence,nsentence,nkey,
     *        nnumber,error_number)
           elseif(keyword.eq.'local') then 
             call uppd_lpar(lpar,mxnlpar,nsets,s0file,sentence,
     *            nsentence,nkey,nlpar,keyword,error_number)
           elseif(keyword.eq.'initval') then
             call uppd_lpar(states0,mxnstat,nsets,s0file,
     *            sentence,nsentence,nkey,nstates0,keyword,
     *            error_number)
           endif

           if(error_number.gt.0) then
             write(6,'(a)') ' ****************************************'
             write(6,'(a)') ' wrong syntax in the modelvar-string.'
             write(6,'(a)') ' ****************************************'
             stop
           endif

           !write the t0- and tf-arrays:
           if(model.eq.'ode'.and.keyword.eq.'odevar') then
             call write_tf(sentence,nsentence,maxnobs,nsets,
     *            t0,tf,ntf)

             !*** matrix-input **** 
             do i = ithsentence,2,-1
               vardim(i,1) = vardim(i-1,1)
               vardim(i,2) = vardim(i-1,2)
             enddo          
             vardim(1,1) = 0
             vardim(1,2) = 0
             !*** matrix-input ****
             !*** odevar & inttype **** 27.10.94
             do i = ithsentence,2,-1
               inttype(i) = inttype(i-1)
             enddo
             inttype(1) = 0
             !*** odevar & inttype **** 27.10.94

             lodevar = 1    !!!! 4.8.94
           endif

           !check if we are in the end of the modelvar-part:
           nsentence = orignsentence        !matrix-input:
           nword = nword - nsentence
           if(nword.gt.0) then 
             do i=1,nword
               word(i) = word(nsentence + i)
             enddo
           endif
         enddo
c loop_end: //////////////////////////////////////////////////////////

       nmodelsentence = ithsentence             !!!! number of sentences
       nx = ncontrol(nmodelsentence)            !!!! 4.8.94
       if(model.eq.'ode'.and.lodevar.eq.0) then !!!! 4.8.94
         nx = nx + 1                            !!!! 4.8.94 
       endif                                    !!!! 4.8.94

c 
c write the dec and (possible) fortran-files:
       call makeincfile(nmlfile,nmodelsentence,
     * type,nglobal,nlocal,ncontrol,name,
     * mxmsent,vardim,                          !!! matrix-input: 
     * projectname,                             !!! 19-09
     * mustcmp,                                 !!! 19-09
     * model,                                   !!! 19-09
     * ftype,                                   !!! 19-09
     * inttype,
     * iperiod,
     * task)

       if(target.eq.' ') then                   !!!! 24.8.94
         !no target-string?                     !!!! 24.8.94
         return                                 !!!! 24.8.94
       endif                                    !!!! 24.8.94



************************target-part:***********************************

       part = 'target'
       lenstr = len(target)

c
c check the characters and write the 
c charcl-table: 
       call checkc(target,lenstr,aascii,charcl,
     *      error_number)

       if(error_number.gt.0) then
         write(6,'(a)') ' ****************************************'
         write(6,'(a)') ' wrong character in the target-string:'
         write(6,'(a)') target(error_number:error_number)
         write(6,'(a)') ' ****************************************'
         stop
       endif

c
c low chars to upper:
       !25.10.94 call case(target,lenstr,2)

c
c find the words from target-string and write
c the word-array:
       call getwords(target,lenstr,charcl,word,nword)

       if(error_number.gt.0) then
          write(6,'(a)') ' ****************************************'
          write(6,'(a)') ' wrong syntax in the target-string.'
          write(6,'(a)') ' ****************************************'
          stop
       endif

       ithsentence = 0

c loop://////////////////////////////////////////////////////////////
       dowhile(nword.gt.0)
           ithsentence = ithsentence + 1

           !find the next sentence:
           call next_sentence(word,nword,part,aascii,sentence,
     *          nsentence,keyword,nkey,error_number,linttype)

           if(error_number.gt.0) then
             write(6,'(a)') ' ****************************************'
             write(6,'(a)') ' wrong syntax in the target-string.'
             write(6,'(a)') ' ****************************************'
             stop
           endif

           !update the laux- and the bound-tables:
           call uppd_lb(part,sentence,nsentence,nmodelsentence,
     *              laux,naux,bound,name,type,nglobal,nlocal,
     *              ncontrol,ninitval,error_number,
********************** bounds_checking: *************************
     *              gpar,mxnlpar,lpar,mxnstat,states0,
     *              s0file)                              !Korjaus 24-04-95
           if(error_number.eq.9) stop
**********************                  *************************

           if(error_number.gt.0) then
             write(*,*) ' ***************************** '
             write(*,*) 'error_number:',error_number
             write(*,*) ' syntax error in target-string:'
             write(*,*) ' ***************************** '
             stop
           endif

 
           !check if we are in the end of the target-part:
           nword = nword - nsentence
           
           if(nword.gt.0) then
             do i = 1,nword
               word(i) = word(nsentence + i)
             enddo
           endif
       enddo
c loop_end: ///////////////////////////////////////////////////////

c task with uppercase:
       !25.10.94 call case(task,3,2)

c
c put the laux-table to the right table (lsim,lest,lexp, ...)
       if(task.eq.'sim') then
         do i=1,4
           do j=1,naux
             lsim(i,j) = laux(i,j)
           enddo
         enddo
         nsim = naux
       elseif(task.eq.'est') then
         do i=1,4
           do j=1,naux
             lest(i,j) = laux(i,j)
           enddo
         enddo
         nest = naux
       elseif(task.eq.'exp') then
         do i=1,4
           do j=1,naux
             lexp(i,j) = laux(i,j)
           enddo
         enddo
         nexp = naux
       elseif(task.eq.'sen') then
         do i=1,4
           do j=1,naux
             lsen(i,j) = laux(i,j)
           enddo
         enddo
         nsen = naux
       elseif(task.eq.'opt') then
         do i=1,4
           do j=1,naux
             lopt(i,j) = laux(i,j)
           enddo
         enddo
         nopt = naux
       endif          

       nbound = naux                             !!!!! 9.8.94
       if(task.eq.'sim'.or.task.eq.'sen') then   !!!! sbound
         do i=1,2
         do j=1,naux                             !!!! sbound
           sbound(i,j) = bound(i,j)              !!!! sbound
         enddo                                   !!!! sbound
         enddo
         nsbound = naux                          !!!! sbound
       endif                                     !!!! sbound
           

**********************codivar-part:*********************************

       if(codirvar.ne.' ') then      !!!! huom. kutsuvaan ohjelm.
                                     !!!! codirvar-stringi ja sille
                                     !!!! alustus = ' '.
         lenstr  = len(codirvar)       
         call checkc(codirvar,lenstr,aascii,charcl,
     *        error_number)

         if(error_number.gt.0) then
           write(6,'(a)') ' ****************************************'
           write(6,'(a)') ' wrong character in the codirvar-string:'
           write(6,'(a)') codirvar(error_number:error_number)
           write(6,'(a)') ' ****************************************'
           stop
         endif


         !low chars to upper chars:
         !25.10.94 call case(codirvar,lenstr,2)

         !find the words from string and write
         !the word-table:
         call getwords(codirvar,lenstr,charcl,word,nword)

c loop://////////////////////////////////////////////////////////////
         dowhile(nword.gt.0)
             ithsentence = ithsentence + 1

             !find the next sentence:
             call next_sentence(word,nword,part,aascii,sentence,
     *            nsentence,keyword,nkey,error_number,linttype)

             if(error_number.gt.0) then
               write(6,'(a)') ' ************************************'
               write(6,'(a)') ' wrong syntax in the codirvar-string.'
               write(6,'(a)') ' ************************************'
               stop
             endif

             call ccodir(sentence,nsentence,name,type,nmodelsentence,
     *            ncodir,codir,naux,laux,error_number,nglobal,nlocal,
     *            ncontrol,ninitval)


             if(error_number.gt.0) then
               write(6,'(a)') ' **************************************'
               write(6,'(a)') ' something wrong in the codirvar-array?'
               write(6,'(a)') ' **************************************'
               stop
             endif

             !check if we are in the end of the codirvar:
             nword = nword - nsentence
             if(nword.gt.0) then
               do i = 1,nword
                 word(i) = word(nsentence + i)
               enddo
             endif
         enddo
c loop_end: ///////////////////////////////////////////////////////
       endif


**********************gridvar-part:*********************************

       if(gridvar.ne.' ') then       !!!! huom. kutsuvassa ohjelmassa
                                     !!!! tulee olla gridvar-stringilla
                                     !!!! alustus = ' '.
         lenstr  = len(gridvar)       
         call checkc(gridvar,lenstr,aascii,charcl,
     *        error_number)

         if(error_number.gt.0) then
           write(6,'(a)') ' ***************************************'
           write(6,'(a)') ' wrong character in the gridvar-string:'
           write(6,'(a)') gridvar(error_number:error_number)
           write(6,'(a)') ' ***************************************'
           stop
         endif


         !low chars to upper chars:
         !25.10.94 call case(gridvar,lenstr,2)

         !find the words from string and write
         !the word-table:
         call getwords(gridvar,lenstr,charcl,word,nword)

c loop://////////////////////////////////////////////////////////////
         dowhile(nword.gt.0)
             ithsentence = ithsentence + 1

             !find the next sentence:
             call next_sentence(word,nword,part,aascii,sentence,
     *            nsentence,keyword,nkey,error_number,linttype)

             if(error_number.gt.0) then
               write(6,'(a)') '****************************************'
               write(6,'(a)') ' wrong syntax in the gridvar-string.'
               write(6,'(a)') '****************************************'
               stop
             endif

             call ccodir(sentence,nsentence,name,type,nmodelsentence,
     *            ngrid,igrid,naux,laux,error_number,nglobal,nlocal,
     *            ncontrol,ninitval)

             if(error_number.gt.0) then
               write(6,'(a)') ' **************************************'
               write(6,'(a)') ' something wrong in the gridvar-array?'
               write(6,'(a)') ' **************************************'
               stop
             endif

             !check if we are in the end of the gridvar:
             nword = nword - nsentence
             if(nword.gt.0) then
               do i = 1,nword
                 word(i) = word(nsentence + i)
               enddo
             endif
         enddo
c loop_end: ///////////////////////////////////////////////////////
       endif

       ngrid = ngrid / 2                    !!!! obs


*********************** nfilepoint:********************************
 
       if(evalfile.ne.'dummy.dat') then
         nfilepoint = 0
         open(3,file=evalfile,status='old')
111      continue
         read(3,'(a)',end=222) ch
         nfilepoint = nfilepoint + 1
         goto 111
222      continue
         close(3)
       endif

*********************** ndiscpoint:********************************
 
       if(discfile.ne.'dummy.dat') then
         ndiscpoint = 0
         open(3,file=discfile,status='old')
112      continue
         read(3,'(a)',end=223) ch
         ndiscpoint = ndiscpoint + 1
         goto 112
223      continue
         close(3)
       endif


*********************** estvar-part:********************************
 
       if(estvar.eq.' ') then
         !no estvar-string:
         !25.10.94 call case(task,3,1)
         return
       endif         

c       if((task.eq.'exp'.and.optcrit.eq.'glo').or.
c     *    (task.eq.'sen'.and.objfun.eq.'exp'.and.
c     *     optcrit.eq.'glo')) then
c          part = 'target'     !!! jotta rajat kirjautuisivat...
c       else
          part = 'estvar'
c       endif

       naux = 0
       lenstr  = len(estvar)       
c
c
c
c check the characters and write the 
c charcl-table: 
       call checkc(estvar,lenstr,aascii,charcl,
     *      error_number)

       if(error_number.gt.0) then
         write(6,'(a)') ' ***************************************'
         write(6,'(a)') ' wrong character in the estvar-string:'
         write(6,'(a)') estvar(error_number:error_number)
         write(6,'(a)') ' ***************************************'
         stop
       endif

c
c low chars to upper:
       !25.10.94 call case(estvar,lenstr,2)

c
c find the words from target-string and write
c the word-array:
       call getwords(estvar,lenstr,charcl,word,nword)

       if(error_number.gt.0) then
          write(6,'(a)') '****************************************'
          write(6,'(a)') ' wrong syntax in the estvar-string.'
          write(6,'(a)') '****************************************'
          stop
       endif

       ithsentence = 0

c loop://////////////////////////////////////////////////////////////
       dowhile(nword.gt.0)
           ithsentence = ithsentence + 1

           !find the next sentence:
           call next_sentence(word,nword,part,aascii,sentence,
     *          nsentence,keyword,nkey,error_number,linttype)

           if(error_number.gt.0) then
             write(6,'(a)') '****************************************'
             write(6,'(a)') ' wrong syntax in the estvar-string.'
             write(6,'(a)') '****************************************'
             stop
           endif

           !update the laux- and the bound-tables:
c           if(task.eq.'exp'.and.optcrit.eq.'glo'.or.
c     *     task.eq.'sen'.and.objfun.eq.'exp'.and.optcrit.eq.'glo') then  !!
c             call uppd_lb(part,sentence,nsentence,nmodelsentence,         !!
c     *              lest,nest,sbound,name,type,nglobal,nlocal,            !!
c     *              ncontrol,ninitval,error_number,
********************** bounds_checking: *************************
c     *              gpar,mxnlpar,lpar,mxnstat,states0,
c     *              s0file)                              !Korjaus 24-04-95
c           if(error_number.eq.9) stop
**********************                  *************************
c           else                                                           !!
             call uppd_lb(part,sentence,nsentence,nmodelsentence,         !!
     *              lest,nest,bound,name,type,nglobal,nlocal,             !!
     *              ncontrol,ninitval,error_number,
********************** bounds_checking: *************************
     *              gpar,mxnlpar,lpar,mxnstat,states0,
     *              s0file)                              !Korjaus 24-04-95
             if(error_number.eq.9) stop
c           endif                                                          !!
**********************                  *************************
           nsbound = nest

           if(error_number.gt.0) then
             write(*,*) ' ***************************** '
             write(*,*) ' syntax error in estvar-string:'
             write(*,*) ' ***************************** '
             stop
           endif

           !check if we are in the end of the estvar-part:
           nword = nword - nsentence
           if(nword.gt.0) then
             do i = 1,nword
               word(i) = word(nsentence + i)
             enddo
           endif
       enddo
c loop_end: ///////////////////////////////////////////////////////

c stop the estvar (and the whole subroutine)
       !25.10.94 call case(task,3,1)  !!! task-string palautetaan pienina kirj.
       return
       end


c================================================================
c   checkc:
c   checks the characters of the string and write the 
c   charcl-table.
c================================================================
       subroutine checkc(str,lenstr,aascii,charcl,iwrong)
       implicit none
c
c
c args:
c   str    - input-string                           - input
c   lenstr - length of the string                   - input
c   aascii - acceptable chars                       - input
c   charcl - class of the char(s)                   - output
c   iwrong - first non-accept. char (or 0)          - output
c            
c
c       values of the charcl:
c       charcl = -1            the character is non-accept.
c              =  0                             separating
c              =  1                             low char
c              =  2                             upper char
c              =  3                             number (or + - .)
c              =  9                             special:  ( ) ; :
c
c
c args:
       character*(*) str
       integer lenstr,aascii(0:127),charcl(0:lenstr+1),iwrong
c
c locals:
       integer i
c
c
       iwrong = 0
       charcl(0) = 0
       do 10 i = 1,lenstr
         charcl(i) = aascii(ichar(str(i:i)))
         if (charcl(i).lt.0) then
           iwrong = i
           return
         endif      
10     continue  
       charcl(lenstr + 1) = 0
       return
       end


c================================================================
c   defchars:
c   defines the acceptable chars and make the classification of the
c   chars: aascii
c modifioin tata nyt siten, etta pilkku on tastedes vain
c erikoismerkki, (joten sita ei voi kayttaa erottimena)
c================================================================
       subroutine defchars(aascii)
c lisasin alaviivan nimissas hyvaksyttavien
c merkkien joukkoon...
c args:
       integer aascii(0:127)
c
c locals:
       integer i
c
c non-accept. chars:
       do 10 i=1,127
         aascii(i) = -1
10     continue
c
c separating chars:    (blanko and comma)
       aascii(32) = 0
c
c
c low chars:
       do 20 i=97,122
         aascii(i) = 1
20     continue
c
c
c upper chars: 
       do 30 i=65,90
         aascii(i) = 2
30     continue
       aascii(ichar('_')) = 2
c
c
c numbers:
       aascii(ichar('.')) = 3
       aascii(ichar('+')) = 3
       aascii(ichar('-')) = 3
       do 40 i=48,57
         aascii(i) = 3
40     continue
c
c
c specials:
       aascii(ichar(',')) = 9    !matrix-input
       aascii(ichar(':')) = 9
       aascii(ichar(';')) = 9
       aascii(ichar('=')) = 9
       aascii(ichar('(')) = 9
       aascii(ichar(')')) = 9
c
       return
       end

c================================================================
c   getwords:
c   separates the words from the string.
c================================================================
       subroutine getwords(str,lenstr,charcl,word,nword)
c
c
c args:
c   str    - input-string                           - input
c   lenstr - stringlength                           - input
c   charcl - class of the characters                - input
c   word  - the word table                         - output
c   nword - the size of the word-table            - input/output
c            
c
c args:
       integer lenstr,nword,charcl(0:lenstr+1)
       character*(*) str,word(*)
c
c locals:
       integer inew,iold,ibegin,iend,i
c
       nword = 0
       inew = charcl(0)
       do 10 i = 1,lenstr+1
         iold = inew
         inew = charcl(i)
         if(iold.ge.1.and.iold.le.3) iold = 1
         if(inew.ge.1.and.inew.le.3) inew = 1
         if    (iold.eq.0.and.inew.eq.1) then
           nword = nword + 1
           ibegin = i
         elseif(iold.eq.0.and.inew.eq.9) then
           nword = nword + 1
           ibegin = i
           iend   = i
           word(nword) = str(ibegin:iend)
         elseif(iold.eq.1.and.inew.eq.0) then
           iend = i - 1
           word(nword) = str(ibegin:iend)
         elseif(iold.eq.1.and.inew.eq.9) then
           iend = i - 1
           word(nword) = str(ibegin:iend)
           nword = nword + 1
           ibegin = i
           iend   = i
           word(nword) = str(ibegin:iend)
         elseif(iold.eq.9.and.inew.eq.1) then
           nword = nword + 1
           ibegin = i
         elseif(iold.eq.9.and.inew.eq.9) then
           nword = nword + 1
           ibegin = i
           iend   = i
           word(nword) = str(ibegin:iend)
         endif
10       continue
         return
         end

c ========================================================================
c nametable:
c updates the name,type and the nglobal, nlocal-...tables
c ========================================================================
c
c
c !!! odevar-muuttujatyypin kasittely
c !!! 15-06-94
c !!!
        subroutine nametable(sentence,nsentence,ithsentence,
     *             keyword,nkey,name,type,nnumber,nglobal,nlocal,
     *             ncontrol,ninitval,ierror,
     *             linttype,inttype)
c
c args:
       character*(*) sentence(*)                           !input
       integer nsentence,ithsentence                       !input
       character*7 keyword                                 !input
       integer nkey                                        !input
       character*(*) name(*)                               !output
       character*7 type(*)                                 !output
       integer nnumber                                     !output
       integer nglobal(*),nlocal(*),ncontrol(*),ninitval(*),!input/output
     *         ierror                                      !output
       logical linttype                                    !input
       integer inttype(*)                                  !output
c
c locals:
       integer ilow,istep,iupp,i
c
c
       name(ithsentence) = sentence(1)
      
       type(ithsentence) = keyword
c  
c
       if(linttype) then
         inttype(ithsentence) = 1
       else
         inttype(ithsentence) = 0
       endif
c
c
       istep = 1
       if(sentence(2).ne.'(') then
c        single:
         ilow  = 1
         iupp  = 1              
       else
c        vector:
c        check the syntax of vector-sentence: 
c        ' name ( ilow : iupp ) keyword ... '
         if(sentence(4).ne.':'.and.sentence(6).ne.':') then
           ierror = 1   
           return
         else
           call readint(sentence(3),ilow)
           if(ilow.ne.1) then 
             ierror = 2 
             return
           endif
           call readint(sentence(5),iupp)
         endif
       endif
c
c number of numbers of the ith sentence: nnumber
       nnumber = (iupp - ilow) / istep + 1
c       
       if(ithsentence.eq.1) then
         nglobal(1)  = 0
         nlocal(1)   = 0
         ncontrol(1) = 0
         ninitval(1) = 0
       else
         nglobal(ithsentence)  = nglobal(ithsentence - 1)
         nlocal(ithsentence)   = nlocal(ithsentence - 1)
         ncontrol(ithsentence) = ncontrol(ithsentence - 1)
         ninitval(ithsentence) = ninitval(ithsentence - 1)
       endif
c
       if(keyword.eq.'global') then
         nglobal(ithsentence) = nglobal(ithsentence) + nnumber
       elseif(keyword.eq.'local') then
         nlocal(ithsentence)  = nlocal(ithsentence) + nnumber
       elseif(keyword.eq.'control') then
         ncontrol(ithsentence) = ncontrol(ithsentence) + nnumber
       elseif(keyword.eq.'initval') then
         ninitval(ithsentence) = ninitval(ithsentence) + nnumber
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c 15-06-94
       elseif(keyword.eq.'odevar') then
         do i=ithsentence,2,-1
           nglobal(i)  = nglobal(i-1)
           nlocal(i)   = nlocal(i-1)
           ncontrol(i) = ncontrol(i-1) + 1
           ninitval(i) = ninitval(i-1)
           type(i)     = type(i-1)
           name(i)     = name(i-1)
         enddo         
         nglobal(1)  = 0
         nlocal(1)   = 0
         ncontrol(1) = 1   
         ninitval(1) = 0
         type(1)     = 'control'
         name(1)     = sentence(1)
       endif
c
       return
c666    continue
c       ierror = 3
c       return
       end

c================================================================
c next_sentence:
c finds the (next) sentence from the word-arrray: 
c
c================================================================
       subroutine next_sentence(word,nword,part,aascii,sentence,
     *            nsentence,keyword,nkey,ierror,linttype)
c
c
c  args:
c    word      -    the table of the (remaining) 
c                   words                         - input
c    nword     -    number of words               - input
c    part      -    the part of the process: 
c                   (modelvar/target/estvar)      - input
c    aascii    -    character-class               - input
c    sentence  -    the table of the words        - input/output
c    nsentence -    number of the words of the
c                   sentence                      - input/output
c    keyword   -    keyword of the sentence       - output
c    nkey      -    index of the keyword          - output
c    ierror    -    error_number                  - input/output     
c    linttype  -    ind. of integer-variabletype  - output

c args:
       character*(*) word(*)
       integer nword
       character*8 part
       character*(*) sentence(*)
       character*7 keyword
       integer nkey,ierror,nsentence
       integer aascii(0:127)
       logical linttype
c
c locals:
       integer i
c max of the length of the words is 80.
       character*80 str
       character*1 first
c
       do i=1,nword
         sentence(i) = ' '
       enddo
c
c
c find the keyword:
       if(part.eq.'target'.or.part.eq.'estvar') then
         keyword = part
c...     nkey = 0
         nkey = 1
         i = nkey
         goto 20
       endif
       i = 0
10     continue
       i = i + 1
       if(i.gt.nword) then
c        no keywords!!!!!!
         ierror = 1
         return
       else
         str = word(i)
       endif
c
       if(str.ne.'global'.and.str.ne.'local'.and.
     *    str.ne.'control'.and.str.ne.'initval'.and.
     *    str.ne.'odevar'.and.str.ne.'globali'.and.
     *    str.ne.'locali') then 
         goto 10
       else
         nkey    = i
         keyword = word(nkey)
       endif
c
c
       if(keyword.eq.'globali') then
         keyword = 'global'
         linttype = .true.
       elseif(keyword.eq.'locali') then
         keyword = 'local'
         linttype = .true.
       else
         linttype = .false.
       endif
c
c find the next variable-name:
       i = nkey
20     continue
       i = i + 1
       if(i.gt.nword) then
         nsentence = nword
       else
         str = word(i)
         first = str(1:1)
!25.10.94         if(aascii(ichar(first)).ne.2.or.str.eq.'file') then
         if(aascii(ichar(first)).ne.1.or.str.eq.'file') then 
           goto 20
         else
           nsentence = i - 1
         endif
       endif
c
       do 30 i = 1,nsentence
         sentence(i) = word(i)
30     continue
       return
       end                   

c================================================================
c   uppd_gpar:
c   updates the gpar-table:
c
c  Modified 18-04-95: The colon_character:;
c
c================================================================
       subroutine uppd_gpar(gpar,ngpar,sentence,nsentence,nkey,
     *            nnumber,ierror)
c
c  args:
c      gpar      -     gpar-table                     - input/output
c      ngpar     -     size of table                  - input/output
c      sentence  -     sentence (table of words)      - input
c      nsentence -     size of sentence               - input
c      nkey      -     the index of the keyword       - input
c      nnumber   -     the number of adding numbers   - input
c      ierror    -     error-number                   - output
c
c args:
       double precision gpar(*)
       integer ngpar,nkey,nnumber,ierror
       character*(*) sentence(*)
c
c locals:
       integer nadd,i,nsentence
c
c
********************** The Colon-character: ***********************
       if(sentence(nsentence).eq.';') then
         nadd = nsentence - nkey - 1
         if(nadd.gt.nnumber.or.nadd.le.0) then
c          wrong number of numbers in the sentence ...
           ierror = 1
           return
         endif
         do i = 1,nadd
           call readdbl(sentence(nkey+i),gpar(ngpar+i))
         enddo
         do i = nadd+1,nnumber
           call readdbl(sentence(nkey+nadd),gpar(ngpar+i))
         enddo  
         ngpar = ngpar + nnumber
         return
       endif
*********************;;;;;;;;;;;;;;;;;;;;;;;;;************************
c
c
       nadd = nsentence - nkey
       if(nadd.eq.1) then 
         do 10 i = 1,nnumber 
           call readdbl(sentence(nkey+1),gpar(ngpar+i))
10       continue
       elseif(nadd.eq.nnumber) then
         do 20 i = 1,nnumber
           call readdbl(sentence(nkey+i),gpar(ngpar+i))
20       continue
       else
c        wrong number of adding numbers in the sentence...
         ierror = 1
         return
       endif
c
       ngpar = ngpar + nnumber
       return
c
c666    continue
c       ierror = 2
c       return
       end

c================================================================
c uppd_lb:
c updates the laux- and the bound-tables.
c 
c version 25-04-94
c version 24-04-95 bounds_checking
c================================================================
       subroutine uppd_lb(part,sentence,nsentence,nmodelsentence,
     *            laux,naux,bound,name,type,nglobal,nlocal,
     *            ncontrol,ninitval,error_number,

********************** bounds_checking: *************************
     *            gpar,mxnlpar,lpar,mxnstat,states0,
     *            s0file                           !Korjaus 24-04-95
     *            )
**********************                  *************************

        

c
c
c  args:
       character*7 part                                     !input
       character*(*) sentence(*)                            !input
       integer nsentence,                                   !input
     *         nmodelsentence,                              !input
     *         laux(4,*),                                   !output
     *         naux                                         !input/output
       double precision bound(2,*)                          !output
       character*(*) name(*)                                !input
       character*7 type(*)                                  !input
       integer nglobal(*),nlocal(*),ncontrol(*),ninitval(*),!input
     *         error_number                                 !output

*************************** bounds_checking: **********************
       integer mxnlpar,mxnstat
       double precision gpar(*),lpar(mxnlpar,*),
     *        states0(mxnstat,*)
       integer s0file(mxnstat,*)                  !Korjaus 24-04-95
***************************                  **********************

c
c  locals:
c the max of the word length = 80
       character*80 variable
       character*7  keyword
       integer ilow,istep,iupp,nfirst,i,
     *         n0,jglob,jcont,jloc,j,
c25-04-94 ======================================================
     *         nsets,    !!! the number of number series
     *         n1        !!! the number of numbers for
                         !!! one dataset   = 3 for local,initval
                         !!! and one obs   = 4 for control-vars.
c25-04-94 =======================================================
c

************************ bounds_checking: ***********************
       logical check_bounds/.true./
       double precision value
       integer ij
*****************************************************************

       if(sentence(2).ne.'(') then
c single:
         variable = sentence(1)
         ilow = 1
         istep= 1
         iupp = 1
         nfirst = 2       !!! the first number
       elseif(sentence(6).eq.')') then
         variable = sentence(1)
         call readint(sentence(3),ilow)
         istep = 1
         call readint(sentence(5),iupp)
         nfirst = 7
       elseif(sentence(8).eq.')') then
         variable = sentence(1)
         call readint(sentence(3),ilow)
         call readint(sentence(5),istep)
         call readint(sentence(7),istep)
         nfirst = 9
       elseif(sentence(4).eq.')') then
         variable = sentence(1)
         call readint(sentence(3),ilow)
         istep = 1
         iupp = ilow
         nfirst = 5
       endif       

**************************** bounds_checking: *********************
       ij = ilow
       if(nfirst.eq.2) ij = 0
****************************                  *********************


c                
c get the name from the nametable:
       i = 0
10     continue
       i = i + 1
       if(i.gt.nmodelsentence) then
         error_number = 1
c        the variable-name not found!!!!
         return
       endif         
       if(variable.ne.name(i)) then
         goto 10
       endif
c
c the type of the variable:
       keyword = type(i)
c             
c25-04-94 ===========================================================
c check the syntax of the variablename:
c if(keyword.ne.'global') then
c it shoud be 'name' or 'name(i)' so
c ilow = iupp and istep = 1.                   
       if(keyword.ne.'global'.and.iupp.ne.ilow) then
         error_number = 7   
         return
       endif
c
c25-04-94 ============================================================ 
       if (i.gt.1) then
       if(keyword.eq.'global') then
         n0 = nglobal(i-1)
       elseif(keyword.eq.'local') then
         n0 = nlocal(i-1)
       elseif(keyword.eq.'control') then
         n0 = ncontrol(i-1)
       elseif(keyword.eq.'initval') then
         n0 = ninitval(i-1)
       endif           
       else
          n0=0
       endif
       j = i
c
c update the laux-table:    
       ilow = ilow + n0
       iupp = iupp + n0
c
       jglob = 0
       jcont = 0
       jloc  = 0
c
c25-04-94 ======================================================================
c n1    = the number of numbers for one vars 
c nsets = the number of number-series 
       if(keyword.ne.'global') then
         if(keyword.eq.'control') then
           n1 = 4
         else
           n1 = 3
         endif
         nsets = (nsentence - nfirst + 1) / n1
c check that we total number of numbers is n1*nsets:
         i = nsentence - nfirst + 1
         if(i.ne.nsets*n1) then
           error_number = 8    
           return
         endif  
         iupp = ilow + nsets - 1
       endif
c25-04-94 ===============================================================================
c
       do 20 i=ilow,iupp,istep
         naux = naux + 1
         if(keyword.eq.'global') then
           laux(1,naux) = 1
         elseif(keyword.eq.'local') then
           laux(1,naux) = 2
         elseif(keyword.eq.'control') then
           laux(1,naux) = 3
         elseif(keyword.eq.'initval') then
           laux(1,naux) = 4
         endif
c25-04-94 ==============================================================================
         if(keyword.eq.'global') then
           laux(2,naux) = i
         else
           laux(2,naux) = ilow
         endif
c25-04-94 ==============================================================================
c
         if(keyword.eq.'global'.or.keyword.eq.'local'.or.
     *      keyword.eq.'initval') then
           laux(3,naux) = 0
         elseif(keyword.eq.'control') then
           call readint(sentence(nfirst+jcont),laux(3,naux))
         endif
         if(keyword.eq.'global') then
           laux(4,naux) = 0
         elseif(keyword.eq.'local') then
           call readint(sentence(nfirst+jloc),laux(4,naux))
         elseif(keyword.eq.'control') then
           call readint(sentence(nfirst+1+jcont),laux(4,naux))
         elseif(keyword.eq.'initval') then
           if(sentence(nfirst + jloc).eq.'file') then
             laux(4,naux) = 0
           else
             call readint(sentence(nfirst+jloc),laux(4,naux))
           endif
         endif
c
c bounds:
         if(part.eq.'estvar') return        !!!! LIS[TTY 5.6.95
         if(part.ne.'estvar') then
           if(keyword.eq.'global') then
             call readdbl(sentence(nfirst + jglob),bound(1,naux))
             call readdbl(sentence(nfirst + 1 + jglob),bound(2,naux))
           elseif(keyword.eq.'local'.or.keyword.eq.'initval') then
             call readdbl(sentence(nfirst + 1 + jloc), bound(1,naux))
             call readdbl(sentence(nfirst + 2 + jloc),bound(2,naux))
           elseif(keyword.eq.'control') then
             call readdbl(sentence(nfirst + 2 + jcont),bound(1,naux))
             call readdbl(sentence(nfirst + 3 + jcont),bound(2,naux))
           endif
         endif

********************** bounds_checking: ******************************
         if    (keyword.eq.'global') then
           value = gpar(laux(2,naux))       
         elseif(keyword.eq.'local') then
           value = lpar(laux(2,naux),laux(4,naux))
*KORJAUS 24-04 *************************************************        
***      elseif(keyword.eq.'in1[Aitval'.and.laux(4,naux).ne.0) then
         elseif(keyword.eq.'initval'.and.
     *          s0file(laux(2,naux),laux(4,naux)).ne.1) then
*KORJAUS 24-04 *************************************************
           value = states0(laux(2,naux),laux(4,naux)) 
         else
           check_bounds = .false.           
         endif
         if(check_bounds) then
           if(bound(1,naux).le.value.and.bound(2,naux).ge.value) then
             continue     !!! OK
           else
             !!! ERROR
             error_number = 9
             write(*,*) '*************************************'

             write(*,*) 'Error:'
             write(*,*) 'The initvalue for the variable: '
             if(ij.eq.0) then
               write(*,*) '   ',variable
             else
               len = index(variable,' ')
             write(*,'(a,a,a,i3,a)') '  ', variable(1:len),'(',ij,')'
             endif
             write(*,*) 'is not in the interval (',bound(1,naux),
     1                  bound(2,naux),')'

             write(*,*) '*************************************'
           endif   
         endif  
************************** end of bounds_checking *************************

         jglob = jglob + 2
         jloc  = jloc  + 3
         jcont = jcont + 4
20     continue
c
c
c ===============================================================
       return
       end
c
c

c================================================================
c uppd_lpar
c updates the lpar- and/or s0file-arrays.
c 
c================================================================
       subroutine uppd_lpar(lpar,mxnlpar,nsets,s0file,sentence,
     *            nsentence,nkey,nlpar,keyword,
     *            ierror)
c
c  args:
c 25-04-94 =======================================================
c      lpar      -     lpar-table (or states0-table)  - input/output
c      mxnlpar   -     size of table                  - input/output
c      nsets     -     the number of datasets         - input
c      s0file    -     'file'-matrix                  - output
c      sentence  -     sentence (table of words)      - input
c      nsentence -     size of sentence               - input
c      nkey      -     the index of the keyword       - input
c      nlpar     -     the number of lpar-vars        - input/output
c      keyword   -     the type of the var            - input
c      ierror    -     error-number                   - output
c
c args:
       integer nsets,mxnlpar
       double precision lpar(mxnlpar,*)
       integer s0file(mxnlpar,*)
       integer nkey,ierror,nlpar
       character*(*) sentence(*)
       character*7 keyword
c
c locals:
       integer nsentence,iinit,nlpar0,hm,ilow,iupp,
     *         idata,iword,ii,jj
       double precision xdouble
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      max of the word-length is 80.
       character*80 word
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c how many number-series (=hm):
        if(sentence(2).ne.'(') then
          hm = 1
        else
c         the syntax must be 'name ( ilow : iupp ) nkey'
          if(sentence(2).ne.'('.or.sentence(4).ne.':'.or.
     *       sentence(6).ne.')') then
c            wrong syntax in vector-name
             ierror = 1
             return
          else
            call readint(sentence(3),ilow)
            if(ilow.ne.1) then
c             wrong syntax in vector-name. lower bound must be 1.
              ierror = 2               
              return
            endif
            call readint(sentence(5),iupp)
          endif
          hm = iupp
        endif
c
c so there are hm vector-elements.
c
c       fix the nlpar-value (= nlpar0)
        nlpar0 = nlpar
        iword  = nkey                   !index of the word of the sentence.
        idata  = 0                      !index of the dataset
        nlpar  = nlpar + 1              !is this right???
10      continue
        iword = iword + 1
        if(iword.gt.nsentence) goto 666
        word = sentence(iword)
        if(word.eq.';') then
c  we must set the current value(s) for all the
c  remaining datasets.
          do jj = idata+1,nsets
            lpar(nlpar,jj) = xdouble
            if(keyword.eq.'initval') s0file(nlpar,jj) = iinit
          enddo
          idata = nsets
        elseif(word.eq.'file') then
          if(keyword.ne.'initval') then
c do'nt use files in local-sentence       
            ierror = 3
            return
          else
            idata = idata + 1
            if(idata.gt.nsets) then
              idata = idata - nsets
              nlpar = nlpar + 1
            endif
            xdouble = 0.
            iinit   = 1
            lpar(nlpar,idata) = xdouble
            s0file(nlpar,idata) = iinit
          endif
        else
          call readdbl(word,xdouble)
          idata = idata + 1
          if(idata.gt.nsets) then
            idata = idata - nsets
            nlpar = nlpar + 1
          endif
          lpar(nlpar,idata) = xdouble
          if(keyword.eq.'initval') then
            iinit = 0
            s0file(nlpar,idata) = iinit
          endif
        endif
        goto 10
c
c
666     continue
c
c now idata shoud be 0 or nsets.
        if(idata.ne.0.and.idata.ne.nsets) then
          ierror = 4
          return
        endif
c
c set the last lpar (s0file)-row for all the 
c remaining vector-elements:
        do ii=nlpar+1,nlpar0 + hm
          do jj=1,nsets
            lpar(ii,jj) = lpar(nlpar,jj)
            if(keyword.eq.'initval') 
     *        s0file(ii,jj) = s0file(nlpar,jj)   
          enddo
        enddo
        nlpar = nlpar0 + hm   !!! update the nlpar-value
c
c
        return
        end


c================================================================
c case:
c lower to upper (icase = 2)
c upper to lower (icase = 1)
c   :
c================================================================
        subroutine case(str,lenstr,into)
c
c  args: str      !input-string                -input/output
c        lenstr   !the length of the string    -input
c        into     !to what                     -input
c
c args:
       character*(*) str
       integer lenstr
       integer into
c
c locals:
       integer iadd,ilow,iupp,iasc,i  
c
c
        if(into.eq.2) then
c pienet isoiksi:
          iadd = -32
          ilow   = 97
          iupp   = 122
        elseif(into.eq.1) then
c isot pieniksi:
          iadd = 32
          ilow   = 65
          iupp   = 90
        else
          return
        endif
c
c
        do 10 i=1,lenstr
          iasc = ichar(str(i:i))
          if(iasc.ge.ilow.and.iasc.le.iupp) then
            iasc = iasc + iadd
            str(i:i) = char(iasc)
          endif
10      continue
        return
        end

c================================================================
c readint:
c reads one integer from string with format.
c muutin sen verran, etta luvussa (= stringissa saa)
c olla piste lopussa.
c 26-04
c================================================================
       subroutine readint(str,inumber)
c
c arsg:
       character*(*) str
       integer inumber
c
c locals:
       integer i,lenstr,ierror
       character*10 iformat
c!!!!!!!!!!!!!!!!!!!!!!!!!!
       iformat = '(i20)'
c!!!!!!!!!!!!!!!!!!!!!!!!!!
       lenstr = len(str)
       i = lenstr + 1
10     continue
       i = i - 1
       if(i.le.0) then
         ierror = 1
         return
       elseif(str(i:i).eq.' ') then
         goto 10
       endif
       lenstr = i
c
       if(lenstr.gt.20) then
         ierror = 2
       else
         if(str(lenstr:lenstr).eq.'.') lenstr = lenstr-1
         read(str(1:lenstr),iformat,err=666) inumber           
       endif
       return
c
666    continue
       ierror = 3
       return
       end

c================================================================
c readdbl:
c reads one double from string with format.
c now it can read input-numbers without period.
!!! modifioitu 4.10.94
!!! 5d12 -> 5.d12
!!!
c================================================================
       subroutine readdbl(str,xdouble)
c
c arsg:
       character*(*) str
       character*10 dformat
       double precision xdouble
c
c locals:
       integer i,lenstr,ierror,ie,ii
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       dformat = '(d20.16)'
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       lenstr = len(str)
       i = lenstr + 1
10     continue
       i = i - 1
       if(i.le.0) then
         ierror = 1
         return
       elseif(str(i:i).eq.' ') then
         goto 10
       endif
       lenstr = i
c
       if(lenstr.gt.20) then
         ierror = 2
       else
c        check if there is a period ('.')
         if(index(str(1:lenstr),'.').le.0) then
           if(index(str(1:lenstr),'e').gt.0) then
c put . before e:      
             ie = index(str(1:lenstr),'e')
             do ii = lenstr+1,ie+1,-1
               str(ii:ii) = str(ii-1:ii-1)
             enddo
             str(ie:ie) = '.'
             lenstr = lenstr + 1
           elseif(index(str(1:lenstr),'d').gt.0) then   !!! 4.10
c put . before d:                                       !!! 4.10
             ie = index(str(1:lenstr),'d')              !!! 4.10
             do ii = lenstr+1,ie+1,-1                   !!! 4.10
               str(ii:ii) = str(ii-1:ii-1)              !!! 4.10
             enddo                                      !!! 4.10
             str(ie:ie) = '.'                           !!! 4.10
             lenstr = lenstr + 1                        !!! 4.10
           else                                      
c put . to last char.:
             lenstr = lenstr+1
             str(lenstr:lenstr) = '.'
           endif
         endif
         read(str(1:lenstr),dformat,err=666) xdouble           
       endif
       return
c
666    continue
       ierror = 3
       return
       end

c ======================================
c noblancs:
c
c removes all blancs from string
c ======================================
c
       subroutine noblancs(str,lenstr)
c
c
       character*(*) str
       integer lenstr,i,j
c
c
       lenstr = len(str)
       j = 0
       i = 0
10     continue
       i = i + 1       
       if(i.gt.lenstr) goto 666
       if(str(i:i).eq.' ') goto 10
       j = j + 1
       str(j:j) = str(i:i)
       goto 10
666    continue
       lenstr = j
       return
       end

c ====================================================
c ccodir:
c writes the codir-arrays from the codirvar-string.
c ====================================================
        subroutine ccodir(sentence,nsentence,name,type,
     *  nmodelsentence,ncodir,codir,naux,laux,ierror,
     *  nglobal,nlocal,ncontrol,ninitval)
c
c
c args:
        character*(*) sentence(*)
        integer nsentence
        character*(*) name(*)
        character*7 type(*)
        integer nmodelsentence
        integer ncodir
        integer codir(*),ierror,nnumber
        integer naux
        integer laux(4,naux),nglobal(*),nlocal(*),
     *          ncontrol(*),ninitval(*)        
c
c
c locals:
        character*20 varname
        character*7 typevar
        integer i,typenum,typeind,
     *  obsnum,dataset !!! 3.8.94   
c
c
        typeind = 0
        typenum = 0
        obsnum  = 0
        dataset = 0
c
c
        varname = sentence(1)
!!!        if(nsentence.eq.1) then
!!!          nnumber = 1
!!!        elseif(sentence(2).eq.'(') then
!!!          call readint(sentence(3),nnumber)
!!!! 22.8.94          nnumber = nnumber - 1
!!!        endif
        nnumber = 1                            !!! 9.9.94
        if(sentence(2).eq.'(' .and. nsentence.gt.1) then !!! 12.9.94
          call readint(sentence(3),nnumber)    !!! 9.9.94
        endif                                  !!! 9.9.94
c
c
c muuttujan tyyppi:
        i = 1
        dowhile(varname.ne.name(i).and.
     *          i.le.nmodelsentence)
          i = i + 1
        enddo
        if(i.gt.nmodelsentence) then
c vaaran niminen muuttuja...
          ierror = 1
          return
        endif       
c tyyppi:
!!!! 22.8.94:
!!!! aikaisemmin typenum = nglobal(i) (nlocal(i) yms.
!!!! nyt typenum = nglobal(i-1)
        typenum=0
        typevar = type(i)
        if(typevar.eq.'global') then
          typeind = 1
          if (i.gt.1) typenum = nglobal(i-1)
        elseif(typevar.eq.'local') then
          typeind = 2
          if (i.gt.1) typenum = nlocal(i-1)
        elseif(typevar.eq.'control') then
          typeind = 3
          if (i.gt.1) typenum = ncontrol(i-1)
        elseif(typevar.eq.'initval') then
          typeind = 4
          if (i.gt.1) typenum = ninitval(i-1)
        endif

        typenum = typenum + nnumber          !!!!! 22.8.94

!!! 3.8.94        typenum = typenum + nnumber - 1
!!!        if(typevar.eq.'global') then      !!! 3.8.94
!!!          typenum = typenum + nnumber - 1 !!! 3.8.94
!!!        endif                             !!! 3.8.94
c
c
!!! 3.8.94 *****************************************
        if(typevar.eq.'local'.or.typevar.eq.'initval') then
          call readint(sentence(2),dataset)
        elseif(typevar.eq.'control') then
          call readint(sentence(2),obsnum)
          call readint(sentence(3),dataset)
        endif    
!!!        ******************************************
c
c
        i = 1
        dowhile((laux(1,i).ne.typeind.or.laux(2,i).ne.     !!! 3.8.94
     *           typenum.or.laux(3,i).ne.obsnum.or.        !!! 3.8.94
     *           laux(4,i).ne.dataset).and.i.le.naux)      !!! 3.8.94
          i = i + 1            
        enddo
c
c
        if(i.gt.naux) then
c         jokin on pielessa mita enimmin...
          ierror = 2
          return
        endif        
        ncodir = ncodir + 1
        codir(ncodir) = i
        return
        end

c ===============================================================
c  makeincfile:
c  
c  luo <projectname>.inc-tiedoston (.tmp-tiedoston), vertaa sita
c  vanhaan <>.inc-tiedostoon ja jos eroavat, niin tasta
c  .tmp-tiedostosta tulee uusi <>.inc-tiedosto. edelleen 
c  talloin rutiini kirjoittaa uudet <projectname>*.f-tiedostot.
c ===============================================================

       subroutine makeincfile(
     * nmlfile,                 ! name of the nml-file (input)
     * nmodelsentence,          ! number of sent. of modelvar-str (input)
     * type,                    ! array of types of vars (input)
     * nglobal,                 ! array of ind. of glob. vars (input)
     * nlocal,                  ! array of ind. of loc. vars (input)
     * ncontrol,                ! array of ind. of contr. vars (input)
     * name,                    ! array of names of vars (input)
     * mxmsent,                 ! max numb of modelsentence (input)
     * vardim,                  ! array of number of matrix-dims (input)
     * projectname,             ! name of the project (input)     !!! 19-09  
     * mustcmp,                 ! (= 0, not-changed, 1 = changed, 2 = new )(output) !!! 19-09  
     * model,                   ! model-string ('ode' or 'alg')   !!! 19-09
     * ftype,                   ! type of the fortran-files       !!! 19-09
     * inttype,                 ! ind. of integer-value-type
     * iperiod,                 ! hh
     * task)                    ! hh
c
c

       implicit none

c subroutine writes the inc-file.
c obs: to chanel 2...
c
c
c args:
       character*(*) nmlfile
       integer nmodelsentence
       character*7 type(*)
       integer nglobal(*),nlocal(*),ncontrol(*)
       character*(*) name(*)
       integer mxmsent,vardim(mxmsent,2)  !matrix-input
       character*7 projectname
       integer mustcmp
       character*3 model
       character*3 ftype
       integer inttype(*)
       character*3 task
c
c
c locals:
       character*80 varname           !!!!! obs.
       character*6 iobs     !',iobs'
       character*80 apu,empty/'  '/
       integer lenapu
       character*5 tblname
       integer iperiod,ij,lenvar,kkk,ilow,iupp,ilen,
     *         nn,mxlenvar
       character*17 typestr
       integer typelen
c
c
       if(projectname.eq.'dummy.t') then
         write(*,*) '===================================='
         write(*,*) 'error: you must give the projectname'
         write(*,*) 'in your nml-file!'
         write(*,*) '===================================='
         stop
c         return
       endif
!!! 19-09      iperiod = index(nmlfile,'.')
!!! 19-09      if(iperiod.eq.0) then
!!! 19-09        iperiod = len(nmlfile)
!!! 19-09         dowhile(nmlfile(iperiod:iperiod).eq.' ')
!!! 19-09           iperiod = iperiod - 1
!!! 19-09         enddo
!!! 19-09       else
!!! 19-09         iperiod = iperiod - 1
!!! 19-09       endif
               iperiod = index(projectname,'.')              !!! 19-09
               if(iperiod.eq.0) then                         !!! 19-09
               iperiod = len(projectname)                    !!! 19-09
                 dowhile(projectname(iperiod:iperiod).eq.' ')!!! 19-09
                   iperiod = iperiod - 1                     !!! 19-09
                 enddo                                       !!! 19-09
               else                                          !!! 19-09
                 iperiod = iperiod - 1                       !!! 19-09
               endif                                         !!! 19-09
c
       open(2,file=projectname(1:iperiod) // '.tmp',status='unknown')    !!! 19-09
       write(2,'(a)') '       integer*4 ii,jj'

******** max of the length of variablesnames, mxlenvar: ********************
       mxlenvar = 0
       do ij=1,nmodelsentence
         varname = name(ij)
         lenvar  = len(varname)
         kkk     = lenvar
         dowhile(varname(kkk:kkk).eq.' ')
           kkk = kkk - 1
         enddo  
         lenvar = kkk
         mxlenvar = max(mxlenvar,lenvar)
       enddo

 
******** write the variable declarations: *******************************
       do ij=1,nmodelsentence
         if(type(ij).ne.'initval') then
           varname = name(ij)
           lenvar  = len(varname)        !!! = 80
           !25.10.94 call case(varname,lenvar,1)   !!! var.name with low.case
           kkk = lenvar
           dowhile(varname(kkk:kkk).eq.' ')
             kkk = kkk - 1
           enddo
           lenvar = kkk
           if(ij.eq.1) then 
             ilow = 0
           else
             if(type(ij).eq.'global') then
               ilow = nglobal(ij-1)
             elseif(type(ij).eq.'local') then
               ilow = nlocal(ij-1)
             elseif(type(ij).eq.'control') then
               ilow = ncontrol(ij-1)
             endif
           endif
c
           apu = ' '
           iobs      = ' ,iobs'
           if(type(ij).eq.'global') then
               iupp    = nglobal(ij)
               tblname = 'gpar'
               ilen    = 1
           elseif(type(ij).eq.'local') then
               iupp = nlocal(ij)
               tblname = 'lpar'
               ilen    = 1
           elseif(type(ij).eq.'control') then
               iupp = ncontrol(ij)
               tblname = 'xdata'
               ilen    = 6
           endif
           nn = iupp - ilow
           if(inttype(ij).eq.1) then
             typestr = '       integer*4 '
             typelen = 17
           else
             typestr = '       real*8 '
             typelen = 14
           endif

           if(nn.le.1) then
              write(2,'(a,a)')        typestr(1:typelen),
     *              varname(1:lenvar)
!matrix-input:           else
           elseif(vardim(ij,1).eq.0) then                  !matrix-input 
               write(apu,'(a,i5,a)') '(',iupp-ilow,')'
               call noblancs(apu,lenapu)
               write(2,'(a,a,a)') typestr(1:typelen),
     *               varname(1:lenvar),apu(1:lenapu)

c ******************matrix-input:*********************
           else
               write(apu,'(a,i5,a,i5,a)') '(',
     *               vardim(ij,1),',',vardim(ij,2),')'
               call noblancs(apu,lenapu)
               write(2,'(a,a,a)') typestr(1:typelen),
     *               varname(1:lenvar),apu(1:lenapu)
c ******************matrix-input********************** 

           endif
         endif
       enddo
c 
c
       write(2,'(a)') 'c '
       write(2,'(a)') 'c '

************* write the assignment statesments *************************
       do ij=1,nmodelsentence
         if(type(ij).ne.'initval') then
           varname = name(ij)
           lenvar  = len(varname)        !!! = 80
           !25.10.94 call case(varname,lenvar,1)   !!! var.name with low.case
           kkk = lenvar
           dowhile(varname(kkk:kkk).eq.' ')
             kkk = kkk - 1
           enddo
           lenvar = kkk
           if(ij.eq.1) then 
             ilow = 0
           else
             if(type(ij).eq.'global') then
               ilow = nglobal(ij-1)
             elseif(type(ij).eq.'local') then
               ilow = nlocal(ij-1)
             elseif(type(ij).eq.'control') then
               ilow = ncontrol(ij-1)
             endif
           endif
c
           iobs      = ' ,iobs'
           if(type(ij).eq.'global') then
               iupp    = nglobal(ij)
               tblname = 'gpar'
               ilen    = 1
           elseif(type(ij).eq.'local') then
               iupp = nlocal(ij)
               tblname = 'lpar'
               ilen    = 1
           elseif(type(ij).eq.'control') then
               iupp = ncontrol(ij)
               tblname = 'xdata'
               ilen    = 6
           endif
           nn = iupp - ilow
           if(nn.le.1) then
             write(apu,'(a,a,i5,a,a)') tblname,
     *             '(',iupp,iobs(1:ilen),')'
             call noblancs(apu,lenapu)
             if(lenvar.lt.mxlenvar) then
               write(2,'(a,a,a,a,a)') '       ',varname(1:lenvar),
     *           empty(1:mxlenvar-lenvar),' = ',apu(1:lenapu)
             else
               write(2,'(a,a,a,a)')   '       ',varname(1:lenvar),
     *                                    ' = ',apu(1:lenapu)
             endif
!matrix-input else
           elseif(vardim(ij,1).eq.0) then !matrix-input
             write(apu,'(i5)') nn
             call noblancs(apu,lenapu)
             write(2,'(a,a)') '       do ii=1,',
     *             apu(1:lenapu)
             if(ilow.gt.0) then
               write(apu,'(a,a,i5,a,a)') tblname,'(ii + ',ilow,
     *               iobs(1:ilen),')'  !!!
             else
               write(apu,'(a,a,a,a)') tblname,'(ii',iobs(1:ilen),')'  !!!!
             endif
             call noblancs(apu,lenapu)
             write(2,'(a,a,a,a)') '         ',varname(1:lenvar),  !!!!
     *             '(ii) = ',apu(1:lenapu)
             write(2,'(a)') '       enddo'
           else
c ********************matrix-input*********
             write(apu,'(i5)') vardim(ij,1)
             call noblancs(apu,lenapu)
             write(2,'(a,a)') '       do ii=1,',
     *             apu(1:lenapu)
             write(apu,'(i5)') vardim(ij,2)
             call noblancs(apu,lenapu)
             write(2,'(a,a)') '       do jj=1,',
     *             apu(1:lenapu)
             if(ilow.gt.0) then
               write(apu,'(a,a,i5,a,i5,a,a)') tblname,'((ii-1)*',
     *               vardim(ij,2),'+jj+',ilow,iobs(1:ilen),')'  !!!
             else
               write(apu,'(a,a,i5,a,a,a)') tblname,'((ii-1)*',
     *                vardim(ij,2),'+jj',iobs(1:ilen),')'       !!!
             endif
             call noblancs(apu,lenapu)
             write(2,'(a,a,a,a)') '         ',varname(1:lenvar),  !!!!
     *             '(ii,jj) = ',apu(1:lenapu)
             write(2,'(a)') '       enddo'
             write(2,'(a)') '       enddo'
c ******************matrix-input*************
           endif
         endif
       enddo

***************** close the tmp-file *********************************
       close(2)
c
c 
********************makedecfile:**************************************
      call makedecfile(projectname,iperiod,mustcmp)
c
c
********************makefortranfiles:*********************************
      call makefortranfiles(model,projectname,iperiod,ftype,task)
c
c
       return
       end


c =====================================================
c makedecfile:
c
c compares the (new) tmp- and the (old) dec-file 
c and 
c - if the dec-file not exist then
c   the tmp-file is the new dec-file
c   mustcmp = 2 
c 
c - the dec-file <> tmp-file then
c   the tmp-file is the new dec-file
c   mustcmp = 1
c
c - the dec-file = tmp-file then
c   mustcmp = 0    
c =====================================================
       subroutine makedecfile(projectname,iperiod,mustcmp)
c
c
c args:
       character*(*) projectname
       integer iperiod,mustcmp
c
c
c locals:
       character*80 decfile,tmpfile
       logical filexist,ldiffer,deletetmp/.true./
c
c
       decfile = projectname(1:iperiod) // '.inc'
       tmpfile = projectname(1:iperiod) // '.tmp'
       inquire(file=decfile,exist=filexist)
       if(.not.filexist) then
         call copyfile(tmpfile,decfile,deletetmp)
         mustcmp = 2
       else
         call cmpfiles(tmpfile,decfile,ldiffer)
         if(ldiffer) then
           call copyfile(tmpfile,decfile,deletetmp)
           mustcmp = 1
         else
           mustcmp = 0
         endif
       endif                  
c
c
       return
       end

c =================================================================
c makefortranfiles:
c writes the fortranfiles <projectname>*.<ftype>:
c =================================================================
      subroutine makefortranfiles(model,projectname,iperiod,ftype,task)
c
c
c  args:
       character*(*) model,projectname,ftype,task     !inputs
       integer iperiod                                !input
c
c
c  locals:
       logical filexist
       character*80 filename
c
c

******************* write <project>m.for: *********************************
       filename = projectname(1:iperiod) // 'm.' // ftype 
       inquire(file=filename,exist=filexist)
       if(.not.filexist) then
         if(model(1:3).eq.'ode') then
           open(2,file=filename,
     *          status='unknown')
           write(2,'(a)') '      subroutine fode(ns,t,s,ds,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           nsaux,nstatea,'
           write(2,'(a)') '     &           states0,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns,nsaux,nstatea   ',
     *        ' !n of state variables'
           write(2,'(a)') '       real*8    t           !time'
           write(2,'(a)') 
     *          '       real*8    s(nstatea) !state variables'
           write(2,'(a)') '       real*8    ds(ns)      !derivatives'
           write(2,'(a)') '       integer*4 nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    states0(nstatea)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') 'c     local variables ',
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '      include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') 'c     user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '      return'
           write(2,'(a)') '      end'
           close(2,status='keep') 
         elseif (model(1:3).eq.'alg') then
           open(2,file=filename,
     *         status='unknown')
           write(2,'(a)') '      subroutine falg(ns,s,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns       ',
     *           '  !n of state variables'
           write(2,'(a)') 
     *           '       real*8    s(ns)       !state variables'
           write(2,'(a)') ' '     
           write(2,'(a)') '       integer*4 nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') '* local variables ', 
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '       include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') '* user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '       return'
           write(2,'(a)') '       end'
           close(2,status='keep')
         elseif (model(1:3).eq.'imp') then
           open(2,file=filename,
     *         status='unknown')
           write(2,'(a)') '      subroutine Fimp(ns,s,f,iflg,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           nsaux,nstatea,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns,nsaux,nstatea ',
     *           '  !n of state variables'
           write(2,'(a)')
     *           '       real*8    s(nstatea)  !state variables'
           write(2,'(a)')
     *           '       real*8    f(ns)       !state variables'
           write(2,'(a)') ' '     
           write(2,'(a)') '       integer*4 iflg,nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') '* local variables ', 
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '       include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') '* user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '       return'
           write(2,'(a)') '       end'
           close(2,status='keep')
         endif
       endif

******************* write <project>o.for: *********************************
       filename = projectname(1:iperiod) // 'o.' // ftype
       inquire(file=filename,exist=filexist)
       if(.not.filexist) then
         if(model.ne.'alg') then
           open(2,file=filename,
     *          status='unknown')
           write(2,'(a)') 
     *           '      subroutine observations(s,ns,yest,ny,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           nsaux,nstatea,'
           write(2,'(a)') '     &           states0,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns,nsaux,nstatea   ',
     *       ' !n of state variables'
           write(2,'(a)') 
     *           '       real*8    s(nstatea) !state variables'
           write(2,'(a)') 
     *           '       integer*4 ny          !n of obs. vars'
           write(2,'(a,a)') '       real*8    yest(ny)   ',
     *        ' !observed variables'
           write(2,'(a)') ' '
           write(2,'(a)') '       integer*4 nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    states0(nstatea)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') '* local variables ', 
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '       include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') 'c      user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '       return'
           write(2,'(a)') '       end'
           close(2,status='keep') 
         else
           open(2,file=filename,
     *          status='unknown')
           write(2,'(a)') 
     *           '      subroutine observations(s,ns,yest,ny,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           nsaux,nstatea,'
           write(2,'(a)') '     &           states0,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns,nsaux,nstatea  ',
     *        '  !n of state variables'
           write(2,'(a)') 
     *           '       real*8    s(nstatea) !state variables'
           write(2,'(a)') '       integer*4 ny'
           write(2,'(a,a)') '       real*8    yest(ny)  ',
     *        '  !observed variables'
           write(2,'(a)') ' '
           write(2,'(a)') '       integer*4 nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    states0(nstatea)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') '* local variables ', 
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '       include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') '* user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '       return'
           write(2,'(a)') '       end'
           close(2,status='keep') 
         endif   
       endif

******************* write <project>i.for: *********************************
       filename = projectname(1:iperiod) // 'i.' // ftype
       inquire(file=filename,exist=filexist)
       if(.not.filexist) then
         if(model.eq.'ode') then
           open(2,file=filename,
     *          status='unknown')
           write(2,'(a)') '      subroutine inits0(ns,t,s,'
           write(2,'(a)') '     &           xdata,nx,nobs,'
           write(2,'(a)') '     &           nsaux,nstatea,'
           write(2,'(a)') '     &           states0,'
           write(2,'(a)') '     &           gpar,ngpar,'
           write(2,'(a)') '     &           lpar,nlpar,'
           write(2,'(a)') '     &           iobs,iset)'
           write(2,'(a)') ' '
           write(2,'(a)') '      implicit none'
           write(2,'(a)') ' '
           write(2,'(a)') 'c     arguments'
           write(2,'(a)') ' '
           write(2,'(a,a)') '       integer*4 ns,nsaux,nstatea   ',
     *           '!n of state variables'
           write(2,'(a)') '       real*8    t           !time'
           write(2,'(a)') 
     *           '       real*8    s(nstatea) !state variables'
           write(2,'(a)') ' '
           write(2,'(a)') '       integer*4 nx,nobs,ngpar,nlpar'
           write(2,'(a)') '       real*8    xdata(nx,nobs)'
           write(2,'(a)') '       real*8    states0(ns)'
           write(2,'(a)') '       real*8    gpar(ngpar)'
           write(2,'(a)') '       real*8    lpar(nlpar)'
           write(2,'(a)') '       integer*4 iobs,iset'
           write(2,'(a)') ' '
           write(2,'(a,a)') '* local variables ', 
     *       '(all user defined variables must be declared here!)'
           write(2,'(a)') ' '
           write(2,'(a,a,a)') '       include ''',
     *           projectname(1:iperiod),'.inc'''
           write(2,'(a)') ' '
           write(2,'(a)') '* user code:'
           write(2,'(a)') ' '
           write(2,'(a)') '       return'
           write(2,'(a)') '       end'
           close(2,status='keep') 
         endif
       endif
c
******************* write <project>u.for: *********************************
       filename = projectname(1:iperiod) // 'u.' // ftype
       inquire(file=filename,exist=filexist)
       if(.not.filexist) then
        if(task(1:3).eq.'opt') then
        open(2,file=filename,
     *  status='unknown')
        write(2,'(a)') '      real*8 function optfun(s,ns,yest,ny,'
        write(2,'(a)') '     &                 xdata,nx,nobs,'
        write(2,'(a)') '     &                 nsaux,nstatea,'
        write(2,'(a)') '     &                 ydata,weight,mostydat,'
        write(2,'(a)') '     &                 states0,'
        write(2,'(a)') '     &                 gpar,ngpar,'
        write(2,'(a)') '     &                 lpar,nlpar,'
        write(2,'(a)') '     &                 iset)'
        write(2,'(a)') ' '
        write(2,'(a)') '      implicit none'
        write(2,'(a)') ' '
        write(2,'(a)') 'c     arguments'
        write(2,'(a)') ' '
        write(2,'(a)') '      integer*4 ns,ny,nx,nobs,nsaux,nstatea'
        write(2,'(a)') '      integer*4 mostydat,ngpar,nlpar,iset'
        write(2,'(a)') '      integer*4 ny(nobs)'
        write(2,'(a)') '      real*8    s(nstatea,nobs)        '
        write(2,'(a)') '      real*8    yest(mostydat,nobs)   '
        write(2,'(a)') '      real*8    xdata(nx,nobs)        '
        write(2,'(a)') '      real*8    ydata(mostydat,nobs)  '
        write(2,'(a)') '      real*8    weight(mostydat,nobs) '
        write(2,'(a)') '      real*8    states0(nstatea)      '
        write(2,'(a)') '      real*8    gpar(ngpar)           '
        write(2,'(a)') '      real*8    lpar(nlpar)           '

        write(2,'(a)') ' '
        write(2,'(a,a)') '* local variables ',
     *    '(all user defined variables must be declared here!)'
        write(2,'(a)') ' '
c        write(2,'(a,a,a)') '       include ''',
c     *        projectname(1:iperiod),'.inc'''
        write(2,'(a)') ' '
        write(2,'(a)') '* user code:'
        write(2,'(a)') ' '
        write(2,'(a)') '      return'
        write(2,'(a)') '      end'
        close(2,status='keep')
        endif
       endif

c
       return
       end

c =====================================================
c write_tf:
c writes tf- and t0-arrays.
c =====================================================
       subroutine write_tf(sentence,nsentence,maxnobs,
     *            nsets,t0,tf,ntf)
c
c
c args:
        character*(*) sentence(*)
        integer nsentence,maxnobs,nsets
        double precision t0(*),
     *  tf(maxnobs,nsets)
        integer ntf(*)
c
c
c locals:
        integer k,   ! index of dataset 
     *  nword,       ! index of word in sentence
     *  nnumber,     ! number of numbers in sentence
     *  ii,j         
c
c
        do k=1,nsets
          ntf(k) = 0
        enddo

        k     = 0
        nword = 0
10      continue        
        nnumber = 0
        k = k + 1
        if(k.eq.1) nword = 2
20      continue
        nnumber = nnumber + 1
        nword = nword + 1
        if(nword.gt.nsentence) goto 666
        if(sentence(nword).eq.';') then
          ntf(k) = nnumber - 2
          goto 10
        elseif(nnumber.eq.1) then
          call readdbl(sentence(nword),t0(k))
        elseif(nnumber.gt.1) then
          call readdbl(sentence(nword),tf(nnumber-1,k))
        endif
        goto 20
666     continue
        if(k.gt.1) then                  !!!!!! 16-06-94 
          do j = k,nsets                 !!!!!! obs obs obs
            t0(j)  = t0(k-1) 
            ntf(j) = ntf(k-1)
            do ii = 1,ntf(j)
              tf(ii,j) = tf(ii,k-1)
            enddo
          enddo
        endif                            !!!!!! 16-06-94
c
c
        return
        end             
c 
c 
c ======================================================
c msentence:
c konvertoi matriisityyppisen lauseen. 
c ======================================================
       subroutine msentence(sentence,     !input/output
     *                      nsentence,    !input/output
     *                      ithsentence,  !input
     *                      mxmsent,      !input
     *                      vardim)       !output
c
c
c args:
       character*(*) sentence(*)
       integer nsentence,ithsentence,mxmsent,
     *         vardim(mxmsent,2)
c
c
c locals:
       integer ii,jj,kk,i
c
c
       call readint(sentence(5),ii)
       call readint(sentence(9),jj)
       kk = ii*jj              
       write(sentence(5),'(i5)') kk
       do i=6,nsentence-4
         sentence(i) = sentence(i+4)
       enddo
       nsentence = nsentence - 4
c
c
       vardim(ithsentence,1) = ii
       vardim(ithsentence,2) = jj
       return
       end

c ============================================================
c cmpfiles:
c
c compares the files: file1 and file2. the logical variable
c ldiffer = .true. iff the file1 != file2. 
c obs: the size of the records is 80!!!!!
c chanels are 7, 8
c 
c Korjattu 5.5.95 kts. 5.5.95 AJN
c ============================================================
       subroutine cmpfiles(file1,file2,ldiffer)
c
c
c

       implicit none

c args:
       character*80 file1,file2
       logical ldiffer
c
c
c locals:
       character*80 rec1,rec2           !!!! obs
       integer i,lines1,lines2
       integer iunit1/7/,             !!!! obs 
     *         iunit2/8/
       character*1 ch
c
c
c number of lines (lines1, lines2):
       open(iunit1,file=file1,status='old')
       i = 0
10     continue
       read(iunit1,'(a)',end=20) ch
       i = i + 1
       goto 10          
20     continue       
       lines1 = i

       open(iunit2,file=file2,status='old')
       i = 0
30     continue
       read(iunit2,'(a)',end=40) ch
       i = i + 1
       goto 30
40     continue
       lines2 = i
c
c
       if(lines1.ne.lines2) then
         ldiffer = .true.
       else
         rewind(iunit1)
         rewind(iunit2)
         rec1 = ' '
         rec2 = ' '
         i    = 0
         do while ((rec1.eq.rec2).and.i.lt.lines1)
           read(iunit1,'(a)') rec1
           read(iunit2,'(a)') rec2
           i = i + 1
         enddo
*** old: if(i.lt.lines1) then  ********************************
         if(rec1.ne.rec2) then               !!! 5.5.95 AJN             
           ldiffer = .true.
         else
           ldiffer = .false.
         endif
       endif
c
c
       close(iunit1)
       close(iunit2)
       return
       end    
                  

c ==============================================================
c copyfile:
c
c copies the records from oldfile to newfile.
c obs: the size of the records is 80.!!!!!!!
c      units are 7,8
c ===============================================================
c
       subroutine copyfile(oldfile,newfile,deleteold)
c args:
       character*80 oldfile,newfile
       logical deleteold
c
c
c locals:
       character*80 rec               !!!! obs
       integer iold/7/,             !!!! obs
     *         inew/8/       
c
c
       open(iold,file=oldfile,status='old')
       open(inew,file=newfile,status='unknown')
c
c
10     continue
       read(iold,'(a)',end=20) rec
       write(inew,'(a)')       rec
       goto 10
20     continue
       if(deleteold) then
         close(iold,status='delete')
       else
         close(iold,status='keep')
       endif
       close(inew,status='keep')
c
c
       return
       end      
             










