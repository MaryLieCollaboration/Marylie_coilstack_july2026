************************************************************************
*  header                      AFRONT                                  *
*   Front end of MARYLIE                                               *
************************************************************************
*
************************************************************************
*      MARYLIE and all its subroutines are copyrighted (1987) by       *
*                        Alex J. Dragt.                                *
*      All rights to MARYLIE and its subroutines are reserved.         *
************************************************************************
*
*  Version Date: 5/25/98
*
***********************************************************************
*  All of MARYLIE 3.0 is believed to conform to Fortran 77 standards.
*  In addition, all of MARYLIE 3.0 is self contained save for
*
*     1) the time functions called in the MARYLIE
*        subroutine cputim
*
*        and
*
*     2) the 'exit' subroutine called (and only called) in the
*        MARYLIE subroutine myexit.
*
*  The two subroutines cputim and myexit and the subroutine mytime
* (which is the only one that calls cputim) are kept separately
*  under the heading "XTRA".
*********************************************************************
c
*********************************************************************
c
*  Main program for MARYLIE 3.0
*  Written originally by Rob Ryne ca 1984 and revised by
*  Petra Schuett November 1987
*
      program mary30
c
      include 'impli.inc'
      include 'files.inc'
c
      dimension p(24)
c
c initiate starting time
      p(1)=1.d0
      p(2)=12.d0
      p(3)=0.d0
      call mytime(p)
c
c write out copyright message at terminal
      write(jof,90)
   90 format(/,1h ,'***MARYLIE 3.0***',
     & /,1h ,'Prerelease Development Version 5/25/98'
     & /,1h ,'Copyright 1987 Alex J. Dragt',
     & /,1h ,'All rights reserved',/)
c
c initialize commons (other than those in block data's)
      call initia
c
c read master input file
      call dumpin
c
c debug output
c      call dump(jodf)
c      call dump(jof)
c
c initialize lie algebraic things
      call tables
      call tables5
      call binom5
      call init
c
c write out status at terminal
      write(jof,100)
  100 format(1h ,'Data input complete; going into #labor.')
c
c begin actual calculations
      call tran
c
c the end of the program is reached in MYEXIT
c
      end
c
***********************************************************************
      block data misc
c-----------------------------------------------------------------------
c initialize miscellaneous common variables
c Written by Petra Schuett, November 1987 based on earlier work of
c Rob Ryne and Liam Healy
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------------------
      include 'bmline.inc'
      include 'coment.inc'
      include 'core.inc'
      include 'files.inc'
      include 'mldex.inc'
      include 'items.inc'
      include 'lims.inc'
      include 'maxcat.inc'
      include 'infin.inc'
      include 'zeroes.inc'
c--------------------
      data latt/mlabor*'        '/
      data num/mlabor*1/
c--------------------
      data mline/maxcmt*' '/
c--------------------
      data na/0/,nb/0/,noble/0/,np/0/,mppu/0/
c--------------------
      data ordcat/4/,topcat/monoms/
      data ordlib/4/,toplib/monoms/
c--------------------
      data bottom/0,1, 7,28, 84,210,462, 924,1716,3003,5005, 8008,12376/
      data top   /0,6,27,83,209,461,923,1715,3002,5004,8007,12375,18563/
c--------------------
      data lf,icf,jfcf,mpi,mpo,jif,jof,jodf,ibrief,iquiet/
     #     11, 13,  14, 15, 16,  5,  6,  12,     2,     0/
c--------------------
      data inuse/maxlum*0/
c--------------------
      data lmade/itmno*0/
c--------------------
      data xinf,     yinf,    tinf,    ginf/
     #     1000.,    1000.,   1000.,   1000./
c--------------------
      data fzer,     detz/
     #     0.d0,     0.d0/
c
c  There is an additional data statement in proc.
c  It has approximately the following form.
c  Its exact contents may be found by looking in subroutine ?? in proc.
c
c      data key/'sm','bm','sf','bf','m(','s(','f(','u(','dz','tx','ty',
c     *'ts','cx','cy','qx','qy','hh','vv','tt','hv','ht','vt','ax','bx',
c     *'gx','ay','by','gy','at','bt','gt','ex','ey','et','wx','wy','wt'/
c      data kask/2*4,4*3,3*2,28*1/
c
      end
c
************************************************************************
c
      block data names
c-----------------------------------------------------------------------
c define type codes allowed in MARYLIE
c Written by Petra Schuett in November 1987 based on earlier
c work of Rob Ryne
c-----------------------------------------------------------------------
      include 'impli.inc'
c---------
c commons
c---------
      include 'codes.inc'
c-----------------------------------------------------------------------
c components of master input file
      data ling/'#comment',
     #          '#beam   ',
     #          '#menu   ',
     #          '#lines  ',
     #          '#lumps  ',
     #          '#loops  ',
     #          '#labor  ',
     #          '#xdata  '/
c
c menu entries
c
c 1: simple elements
c
      data (ltc(1,j),j=1,30)/
     #  'drft    ','nbnd    ','pbnd    ','gbnd    ','prot    ',
     #  'gbdy    ','frng    ','cfbd    ','quad    ','sext    ',
     #  'octm    ','octe    ','srfc    ','arot    ','twsm    ',
     #  'thlm    ','cplm    ','cfqd    ','dism    ','sol     ',
     #  'mark    ','jmap    ','dp      ','recm    ','spce    ',
     #  'cfrn    ','coil    ','intg    ','rmap    ','arc     '/
      data (nrp(1,j),j=1,30)/
     #    1,6,4,6,2,
     #    4,5,6,4,2,
     #    2,2,6,1,4,
     #    6,5,4,5,6,
     #    0,0,0,6,1,
     #    4,11,6,6,4/
c
c 2: user-supplied elements
c
      data (ltc(2,j),j=1,20)/
     #  'usr1    ','usr2    ','usr3    ','usr4    ','usr5    ',
     #  'usr6    ','usr7    ','usr8    ','usr9    ','usr10   ',
     #  'usr11   ','usr12   ','usr13   ','usr14   ','usr15   ',
     #  'usr16   ','usr17   ','usr18   ','usr19   ','usr20   '/
      data (nrp(2,j),j=1,20)/20*6/
c
c 3: parameter sets
c
      data (ltc(3,j),j=1,9)/
     #  'ps1     ','ps2     ','ps3     ','ps4     ','ps5     ',
     #  'ps6     ','ps7     ','ps8     ','ps9     '/
      data (nrp(3,i),i=1,9)/9*6/
c
c 4: random elements
c    Note: The j indices for these entries must be aligned with
c          the j indices for their "simple element" counterparts.
c          This is done by using "dummy" entries.
c
      data (ltc(4,j),j=1,24)/
     #  'rdrft   ','rnbnd   ','rpbnd   ','rgbnd   ','rprot   ',
     #  'rgbdy   ','rfrng   ','rcfbd   ','rquad   ','rsext   ',
     #  'roctm   ','rocte   ','rsrfc   ','rarot   ','rtwsm   ',
     #  'rthlm   ','rcplm   ','rcfqd   ','rdism   ','rsol    ',
     #  'dummark ','dumjmap ','dumdp   ','rrecm   '/
      data (nrp(4,i),i=1,24)/24*2/
c
c 5: random user-supplied elements
c
      data (ltc(5,j),j=1,9)/
     #  'rusr1   ','rusr2   ','rusr3   ','rusr4   ','rusr5   ',
     #  'rusr6   ','rusr7   ','rusr8   ','rusr9   '/
      data (nrp(5,i),i=1,9)/9*2/
c
c 6: random parameter sets
c
      data (ltc(6,j),j=1,9)/
     #  'rps1    ','rps2    ','rps3    ','rps4    ','rps5    ',
     #  'rps6    ','rps7    ','rps8    ','rps9    '/
      data (nrp(6,i),i=1,9)/9*2/
c
c 7: simple commands
c
      data (ltc(7,j),j=1,40)/
     #  'rt      ','sqr     ','symp    ','tmi     ','tmo     ',
     #  'pmif    ','circ    ','stm     ','gtm     ','end     ',
     #  'ptm     ','iden    ','whst    ','inv     ','tran    ',
     #  'revf    ','rev     ','mask    ','num     ','rapt    ',
     #  'eapt    ','of      ','cf      ','wnd     ','dwnd    ',
     #  'ftm     ','wps     ','time    ','cdf     ','bell    ',
     #  'wmrt    ','wcl     ','paws    ','inf     ','dims    ',
     #  'zer     ','wuca    ','tpol    ','dpol    ','cbm     '/
      data (nrp(7,i),i=1,40)/
     #   6,0,2,4,1,
     #   3,6,1,2,0,
     #   5,0,2,0,0,
     #   1,0,4,4,4,
     #   3,6,6,5,6,
     #   4,2,3,1,0,
     #   2,3,0,5,6,
     #   3,4,6,5,3/
c
c 8: advanced commands
c
      data (ltc(8,j),j=1,39)/
     #  'cod     ','amap    ','dia     ','dnor    ','exp     ',
     #  'pdnf    ','psnf    ','radm    ','rasm    ','sia     ',
     #  'snor    ','tadm    ','tasm    ','tbas    ','gbuf    ',
     #  'trsa    ','trda    ','smul    ','padd    ','pmul    ',
     #  'pb      ','pold    ','pval    ','fasm    ','fadm    ',
     #  'sq      ','wsq     ','ctr     ','asni    ','pnlp    ',
     #  'csym    ','psp     ','mn      ','bgen    ','tic     ',
     #  'ppa     ','moma    ','geom    ','fwa     '/
      data (nrp(8,i),i=1,39)/
     #   6,6,5,5,3,
     #   5,5,5,5,5,
     #   5,4,6,1,2,
     #   6,6,6,3,3,
     #   3,5,3,6,6,
     #   4,6,4,6,4,
     #   1,3,2,6,6,
     #   6,6,6,6/
c
c 9: procedures and fitting and optimization
c
      data (ltc(9,j),j=1,37)/
     #  'bip     ','bop     ','tip     ','top     ',
     #  'aim     ','vary    ','fit     ','opt     ',
     #  'con1    ','con2    ','con3    ','con4    ','con5    ',
     #  'mrt0    ',
     #  'mrt1    ','mrt2    ','mrt3    ','mrt4    ','mrt5    ',
     #  'fps     ',
     #  'cps1    ','cps2    ','cps3    ','cps4    ','cps5    ',
     #  'cps6    ','cps7    ','cps8    ','cps9    ',
     #  'dapt    ','grad    ','rset    ','flag    ','scan    ',
     #  'mss     ','gendip  ','taylor  '/
      data (nrp(9,i),i=1,37)/
     #   1,1,1,1,
     #   6,6,6,6,
     #   6,6,6,6,6,
     #   1,
     #   6,6,6,6,6,
     #   1,
     #   6,6,6,6,6,
     #   6,6,6,6,
     #   4,6,5,6,6,
     #   6,9,0/
c
      end
c
c**********************************************************************
c
      subroutine buffin(th,tmh,thsave,tmhsav,lumpno)
c-----------------------------------------------------------------------
c     stores a map and polynomial coeffs into a buffer
c  Written by J. Howard, Fall 1986
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'stack.inc'
c
      dimension th(monoms),tmh(6,6)
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c
      do 100 i = 1,6
      do 100 j = 1,6
      tmhsav(i,j,lumpno) = tmh(i,j)
100   continue
c
      do 200 i = 1,monoms
      thsave(i,lumpno) = th(i)
200   continue
      return
      end
************************************************************************
      subroutine bufout(thsave,tmhsav,th,tmh,lumpno)
c-----------------------------------------------------------------------
c     reads a map and polynomial coefficients out of a buffer
c  Written by J. Howard, Fall 1986
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'stack.inc'
c
      dimension th(monoms),tmh(6,6)
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c
      do 100 i = 1,6
      do 100 j = 1,6
      tmh(i,j) = tmhsav(i,j,lumpno)
100   continue
c
      do 200 i = 1,monoms
      th(i) = thsave(i,lumpno)
200   continue
c
      return
      end
************************************************************************
      subroutine cnumb(string,num,lnum)
c-----------------------------------------------------------------------
c  This routine finds out, whether string codes an integer number.
c  if so, it converts it to num
c
c  Input: string character*10 input string
c  Output:num    integer      correspondent number
c         lnum   logical      =.true. if string is a number
c
c  Author: Petra Schuett
c          October 19, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
c
      integer num
      character string*10
      logical lnum
c
      character*10 digits
      save digits
      data digits /'0123456789'/
c-----------------
c The first char may be minus or digit
c      write(jodf,*)'cnumb: string= ',string
      if(index(digits,string(1:1)).eq.0 .and. string(1:1).ne.'-') then
        lnum=.false.
c        write(jodf,*)'first char is no digit'
        return
      endif
c All other characters must be digits...
      do 1 k=2,10
      if(index(digits,string(k:k)).eq.0) then
c ...or trailing blanks
        if(string(k:10).ne.' ') then
c        write(jodf,*)'string(',k,':10) is not blank'
         lnum=.false.
         return
        else
         goto 11
        endif
      endif
   1  continue
c string is a number
  11  read(string,*,err=999) num
      lnum=.true.
      return
c-----------------
c error exit
 999  write(jof ,99) string
      write(jodf,99) string
 99   format(' ---> error in cnumb: string ',a10,' could not be'
     &     , ' converted to number')
      call myexit
      end
c
************************************************************************
c
      subroutine comwtm(j,jrep)
c-----------------------------------------------------------------------
c  routine to combine the jth lump with the total map n times
c  Written by J. Howard, Fall 1986
c  Modified by Alex Dragt, 15 June 1988, to save storage
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'map.inc'
      include 'core.inc'
c----------
c start
c----------
c make n positive in case lump has a negative repetition number
      n=iabs(jrep)
c
      do 100 i=1,n
      call concat(th,tmh,thl(1,j),tmhl(1,1,j),th,tmh)
  100 continue
      return
      end
c
*******************************************************************
c
      subroutine cqlate(icfile,norder,ntimes,nwrite,isend)
c
c-----------------------------------------------------------------------
c  circulate ntimes times; print every nwrite
c
c  input : icfile (integer) = filenumber for reading input rays
c          norder (integer) = NOT USED
c          ntimes (integer) = number of turns through actual line
c          nwrite (integer) = modulo for writing coordinates to jfcf
c          isend  (integer) = 1 output only to jof
c                           = 2 output only to jodf
c                           = 3 both
c  in common /files/: jfcf  < 0 full precision coordinates
c                           > 0 standard format coordinates
c
c  Written by Rob Ryne ca 1984
c  slightly changed by Petra Schuett (labels, if-then-else ...)
c                      October 30,1987
c------------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'items.inc'
      include 'elmnts.inc'
      include 'rays.inc'
      include 'talk.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'map.inc'
      include 'deriv.inc'
      include 'core.inc'
      include 'nturn.inc'
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
      character*8 string(5),str
      logical     ljof,ljodf
c-----------------------------------------------------------------------
c  start routine
c-----------------------------------------------------------------------
c  see if a loop exists
      if(nloop.le.0) then
      write(jof ,*) ' error from cqlate: no loop has been specified'
      write(jodf,*) ' error from cqlate: no loop has been specified'
      return
      endif
c
c  read initial conditions, if requested:
      if(icfile.gt.0) call raysin(icfile)
c--------------------
c  write items in loop, if requested:
      ljof  = isend.eq.1 .or. isend.eq.3
      ljodf = isend.eq.2 .or. isend.eq.3
      if (ljof .or. ljodf) then
c  write loop name
      if(ljof ) write(jof ,510) ilbl(nloop)
      if(ljodf) write(jodf,510) ilbl(nloop)
  510 format(1h ,'circulating through ',a8,' :')
c  write loop contents
      do 1 jtot=0,joy,5
       kmax=min(5,joy-jtot)
       do 2 k=1,kmax
        jk1=jtot+k
c element
        if(mim(jk1).lt.0) then
          string(k)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(k)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(k)=ilbl(inuse(mim(jk1)))
        endif
   2  continue
      if(ljof)  write(jof ,511)(string(k),k=1,kmax)
      if(ljodf) write(jodf,511)(string(k),k=1,kmax)
  511 format(' ',5(1x,a8))
   1  continue
      endif
c--------------------
c  circulation procedure
c
c  set up control indices
      ibrief=0
      kwrite=nwrite
      if(kwrite.eq.0)kwrite=1
c--------------------
c  circulation through loop
      do 1000 nturn=1,ntimes
c----------
c  handle each loop-element
      do 100 kk=1,joy
c  check to see if all particles lost
c  Go at least one turn...
      if(nturn.gt.1) then
      if(nlost.ge.nrays) then
       write(jof,*) 'all particles lost'
       return
      endif
      endif
c
c check to see if item is a lump, user routine, or element
      ip = mim(kk)
      if(ip.gt.5000)then
c procedure for a user routine
        nip=ip-5000
        if(nt2(nip).eq.1) call user1(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.2) call user2(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.3) call user3(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.4) call user4(pmenu(1+mpp(nip)))
        if(nt2(nip).eq.5) call user5(pmenu(1+mpp(nip)))
c
      else if(ip.lt.0) then
c procedure for turtling thru elements (which may be in lines):
        call lmnt(nt1(-ip),nt2(-ip),pmenu(1+mpp(-ip)),1)
c
      else
c procedure for a lump
       do 110 nn=1,nrays
c check to see if particle has already been lost
       if (istat(nn).eq.0) then
         do 112 l=1,6
  112    zi(l)=zblock(nn,l)
c call symplectic tracker
         jwarn=0
         call evalsr(tmhl(1,1,ip),zi,zf,dfl(1,1,ip),
     #               rdfl(1,1,ip),rrjacl(1,1,1,ip))
c check to see if particle was 'lost' by the symplectic ray tracer
         if (jwarn.ne.0) then
           istat(nn)=nturn
           nlost=nlost + 1
           ihist(nlost,1)=nturn
           ihist(nlost,2)=nn
         endif
c copy results of ray trace into storage array
         do 114 l=1,6
  114    zblock(nn,l)=zf(l)
       endif
  110  continue
      endif
  100 continue
c
c end of one turn
c----------
c check to see if results should be written out
      if(mod(nturn,kwrite).eq.0) then
      do 200 m1=1,nrays
      if(istat(m1).eq.0) then
c  procedure for writing out results of ray trace
c
c  procedure for standard format
       if(jfcf.gt.0) then
        write(jfcf,520)(zblock(m1,m2),m2=1,6)
  520   format(6(1x,1pe12.5))
c
c  procedure for full precision
       else if(jfcf.lt.0) then
        do 525 m2=1,6
        write(-jfcf,*) zblock(m1,m2)
  525   continue
       endif
c
      endif
  200 continue
      endif
c
 1000 continue
c end of loop thru turns
c-----------------------------------------------------------------------
      return
      end
c
***********************************************************************
c
      subroutine cread(kbeg,msegm,line,string,lfound)
c-----------------------------------------------------------------------
c  This routine searches line(kbeg:80) for the next string; strings are
c  delimited by ' ','*' or ','.
c
c  Input: line   character*80 input line
c         kbeg   integer      line is only searched behind kbeg
c                             if a string is found, kbeg is set to
c                             possible start of next string
c         msegm  integer      segment number. for msegm=1 (comment
c                             section), the length of a string is not
c                             checked.
c  output:string character*10 found string
c         lfound logical      =.true. if a string was found
c
c  Written by Rob Ryne ca 1984
c  Rewritten by Petra Schuett
c          October 19, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
c
      integer kbeg
      character line*80, string*10
      logical lfound
c
c look for first character of string
c
      do 1 k=kbeg,80
       if(     (line(k:k).ne.' ')
     &    .and.(line(k:k).ne.'*')
     &    .and.(line(k:k).ne.',')) then
c found:
         lfound=.true.
c string has max 10 characters
         lmax  = min(k+10,80)
c look for delimiter
         do 2 l=k,lmax
          if(    (line(l:l).eq.' ')
     &       .or.(line(l:l).eq.'*')
     &       .or.(line(l:l).eq.',')) then
c if found, string is known
            string=line(k:l-1)
            kbeg  =l+1
c done.
            return
          endif
  2      continue
c no delimiter behind string found...
         string=line(k:lmax)
c ...end of line
         if(lmax.eq.80) then
           kbeg=80
c ...or string too long (ignored for comment section)
         else if (msegm.ne.1) then
           write(jof,99) kbeg,line
  99       format(' ---> warning from cread:'/
     &            '      the following line has a very long string at ',
     &            'position ',i2,' :'/'      ',a80)
           kbeg=lmax+1
         endif
c ...anyway, its done.
         return
       endif
  1   continue
c No string was found after all.
      lfound=.false.
      return
      end
c
c************************************************
      subroutine dumpin
c-----------------------------------------------------------------------
c  This routine organizes the data input from file lf, the master input
c  file.
c  This file is divided into "components" beginning with a code "#..."
c  The available codes are given in common/sharp/.
c  In dumpin, they are numbered as they occur in that common. The
c  component (segment) currently being read has number "msegm".
c  The entries of the component "#menu" will be called "elements".
c  The term "item" is used for entries of the components "#lines",
c  "#lumps" and "#loops". Entries in "#labor" will be called "tasks".
c
c  Output is transferred via several commons. They are explained below.
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett October 21, 1987
c  Internal text data segment #xdata added by Tom Mottershead, Feb 2012
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
c
      include 'bmline.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'coment.inc'
      include 'mldex.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'cxdata.inc'
c-----------------------------------------------------------------------
c local variables:
      character*8 strarr(40),string
      character*80 line
      integer narr(40)
      logical leof,lcont
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
      leof = .false.
      rewind lf
      msegm = -1
      write(jof,13)
      write(jodf,13)
  13  format('Marylie 2014 version with #xdata')
c--------------------
c read first line of master input file (should set msegm)
c
  10  call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
c--------------------
c new component (segment) begins, branch to appropriate part of code
c
   1  goto(100,200,300,400,500,600,700,800),msegm
c
c error exit:
      write(jof ,99)
      write(jodf,99)
  99  format(1x,'problems at 1st goto of routine dumpin')
      call myexit
c--------------------
c  #comment
c
  100 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 1) goto 1
c
      np=np+1
      if(np.gt.maxcmt) then
        write(jof ,199) maxcmt
        write(jodf,199) maxcmt
 199    format(1x,'error in dumpin:',
     # ' too many comment lines (>= maxcmt = ',i6,') in #comment')
        call myexit
      endif
      mline(np)=line
      goto 100
c--------------------
c  #beam
c
  200 read(lf,*,err=290,end=1000)brho,gamm1,achg,sl
c  computation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      ts=sl/c
      goto 10
c
c error exit:
  290 write(jof ,299)
      write(jodf,299)
  299 format(1h ,'data input error detected by dumpin near #beam')
      call myexit
c--------------------
c  #menu
c
  300 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 3) goto 1
c
      na=na+1
      if(na.gt.mnumax) then
        write(jof ,399) mnumax
        write(jodf,399) mnumax
  399   format(1h ,'error in dumpin:'  ,
     & ' too many items (>= mnumax = ',i6,') elements in #menu')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,396) strarr(1)
        write(jodf,396) strarr(1)
  396   format(' error detected by dumpin in #menu: name ',
     & a8,' is doubly defined'/
     &         ' the second definition will be ignored')
        na = na-1
        goto 300
      endif
c check if element/command name is present ( F. Neri 4/14/89 ):
      if ( itot .lt. 2 ) then
        write(jof ,1396) strarr(1)
        write(jodf ,1396) strarr(1)
 1396   format(' error in #menu: ',a8,' has no type code!')
        call myexit
      endif
      if ( itot .gt. 2 ) then
        write(jof ,1397) strarr(1)
        write(jodf ,1397) strarr(1)
 1397   format(' error in #menu: ',a8,' has more than one type code!')
        call myexit
      endif
c new item in menu:
      lmnlbl(na)=strarr(1)
c string is name of element/command type, look up element indices
      string=strarr(2)
      do 325 m = 1,9
      do 325 n = 1,40
      if(string.eq.ltc(m,n)) then
        nt1(na) = m
        nt2(na) = n
c read parameters. Number of parameters is given in nrp
        imax=nrp(m,n)
        if(imax.eq.0)goto 300
        mpp(na) = mppu
        read(lf,*,err=390)(pmenu(i+mpp(na)),i=1,imax)
        mppu = mppu + imax
        goto 300
c normal end of a menu element/command
      endif
  325 continue
c
c error: unknown element/command name
      write(jof ,398)(strarr(j),j=1,2)
      write(jodf,398)(strarr(j),j=1,2)
  398 format(1h ,'dumpin error in ',a8,': type code ',a8,' not found.'/
     #       1h ,'this item will be ignored')
      na=na-1
      goto 300
c
c error in parameter input
  390 write(jof ,397)lmnlbl(na)
      write(jodf,397)lmnlbl(na)
  397 format(1h ,'dumpin data input error at item ',a8)
      call myexit
c--------------------
c  #lines,#lumps,#loops
c
  400 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 4) goto 1
      goto 410
  500 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 5) goto 1
      goto 410
  600 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 6) goto 1
c
  410 nb=nb+1
      if(nb.gt.itmno)then
        write(jof ,499) itmno
        write(jodf,499) itmno
  499   format(1h ,'error in dumpin:' ,
     & ' too many lines, lumps, and loops (sum >= itmno = ',i6, ')')
        call myexit
      endif
c check for doubly defined names
      call lookup(strarr(1),itype,indx)
      if(itype.ne.5) then
        write(jof ,497) ling(msegm),strarr(1)
        write(jodf,497) ling(msegm),strarr(1)
  497   format(1x,'dumpin error in ',
     &  a8,': name ',a8,' is doubly defined'/
     &         ' the second definition will be ignored')
        na = na-1
        if(msegm .eq.4) then
          goto 400
        else if(msegm.eq.5) then
          goto 500
        else if(msegm.eq.6) then
          goto 600
        endif
      endif
c new item
      ilbl(nb)=strarr(1)
c read components of item
      imin=0
c--
c repeat...
  420 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if(imin+itot.gt.itmln)then
        write(jof ,498) itmln,ilbl(nb)
        write(jodf,498) itmln,ilbl(nb)
  498   format(1h ,'error in dumpin:',
     & ' too many entries (> itmln = ',i6,') in ',a8)
        call myexit
      endif
c store names and rep rates of components
      do 425 i=1,itot
        icon(imin+i,nb)=strarr(i)
        irep(imin+i,nb)=narr(i)
  425 continue
      imin=imin+itot
      if(lcont) goto 420
c ... until no more continuation lines
c--
c now set length and type of element
      ilen(nb)=imin
      ityp(nb)=msegm-2
c go back to appropriate component (segment)
      if(msegm .eq.4) then
        goto 400
      else if(msegm.eq.5) then
        goto 500
      else if(msegm.eq.6) then
        goto 600
      endif
c--------------------
c  #labor
c
  700 call rearec(line,leof,msegm,strarr,narr,itot,lcont)
      if (leof) goto 1000
      if (msegm .ne. 7) goto 1
c
      noble=noble+1
      if(noble.gt.mlabor) then
        write(jof ,799) mlabor
        write(jodf,799) mlabor
  799   format(1h ,'error in dumpin:',
     & ' too many entries (>= mlabor = ',i6,') in #labor array')
        call myexit
      endif
c new task
      latt(noble)=strarr(1)
      num(noble)=narr(1)
      goto 700
ctm2014-----------------------
c  #data  internal text data segment
  800 continue
      jrec = 0
      nxd = 0
      jj = 0
 810  read(lf,117,end=1000) line
 117  format(a)
c
c     first check for new Marylie segment
c
      nseg = 0
      do 820 ii = 1,8
      nseg = index(line,ling(ii))
      if(nseg.ne.0) then
         msegm = ii
         go to 890
      endif
  820  continue
      jj = jj + 1
ctm   write(6,819) jj, line
 819  format('line',i3,a)
      kk = index(line,'ixd>:')
      if(kk.gt.0) then
         read(line(kk+6:),*) kurid
ctm      write(6,*) kk,'=kk',kurid,'=kurid'
         go to 810
      endif
      jrec = jrec + 1
      if(jrec.gt.maxrec) go to 890
      if(nxd.lt.jrec) nxd = jrec
      mltext(jrec) = line
      lunex(jrec) = kurid
      go to 810
      return
 890  continue
      if(nxd.le.0) go to 1
      write(jof,*) nxd,' lines read into mltext common.'
      write(jodf,*) nxd,' lines read into mltext common.'
      do 895 j = 1,nxd
         write(jof,897) j,lunex(j),mltext(j)
         write(jodf,897) j,lunex(j),mltext(j)
 897     format(i3,': id=',i3,2x,a)
 895  continue
      go to 1
c--------------------
c normal return at end of file
 1000 continue
      return
      end
c
************************************************************************
      subroutine dump(iu)
c-----------------------------------------------------------------------
c  this subroutine dumps the status of the common-blocks
c  as they are filled by dumpin
c  input: iu   output-file
c
c  written by Petra Schuett
c             October 26, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'bmline.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'coment.inc'
      include 'mldex.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c start routine
c-----------------------------------------------------------------------
c comments
      if(np.ne.0) then
        write(iu,500) ling(1)
        write(iu,*) (mline(i),i=1,np)
      endif
c--------------------
c  beam
      write(iu,500) ling(2)
      write(iu,*) brho
      write(iu,*) gamm1
      write(iu,*) achg
      write(iu,*) sl
c--------------------
c  menu
      write(iu,500) ling(3)
      do 10 k=1,na
       write(iu,520) lmnlbl(k),ltc(nt1(k),nt2(k))
       imax=nrp(nt1(k),nt2(k))
       if(imax.ne.0) then
        write(iu,522)(pmenu(i+mpp(k)),i=1,imax)
       endif
   10 continue
c--------------------
c  lines,lumps,loops
      if(nb.ne.0) then
       do 40 ii=2,4
       write(iu,500) ling(ii+2)
        do 40 k=1,nb
        if(ityp(k).eq.ii) then
         write(iu,530) ilbl(k)
         write(iu,532)(irep(l,k),icon(l,k),l=1,ilen(k))
        endif
   40  continue
      endif
c--------------------
c labor
      if(noble.ne.0) then
      write(iu,500) ling(7)
      do 100 j=1,noble
       write(iu,540) num(j),latt(j)
  100 continue
      endif
      return
c-----------------------------------------------------------------------
c format
c-----------------------------------------------------------------------
  500 format(1h ,a8)
  510 format(1h ,a80)
  520 format(1h ,1x,a8,1x,a8)
  522 format((1h ,3(1x,1pg22.15)))
  530 format(1h ,1x,a8)
  532 format((1h ,1x,5(i5,'*',a8),1x,:'&'))
  540 format(1h ,1x,i4,'*',a8)
      end
************************************************************************
      subroutine initia
c-----------------------------------------------------------------------
c  This routine initializes some constants in common blocks, which are
c  not initialized in block data
c
c  Petra Schuett  10/30/87
c  Alex Dragt 6/20/88
c-----------------------------------------------------------------------
c
      include 'impli.inc'
      include 'pie.inc'
      include 'parm.inc'
      include 'param.inc'
      include 'frnt.inc'
c
c set up constants
c
c      pi    = 3.14159265358979323846264d0
      pi = 4.d0*atan(1.d0)
      pi180 = pi/180.d0
      twopi = 2.d0*pi
      c     = 2.99792458d+08
c
c set up default values for cfbd fringe field parameters
c
      cfbgap=0.d0
      cfblk1=.5
      cfbtk1=.5
c
      return
      end
c
***********************************************************************
      subroutine initst(string,nrept)
c-----------------------------------------------------------------------
c initialize stack, putting element 'string' with rep. factor nrept
c on top
c
c input: string character*8 initial stack element (loop,line or lump!)
c        nrept  integer     rep. factor
c
c  Petra Schuett, October 30,1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'items.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      character*8 string
c-------
c start
c-------
      np = 1
      lstac(1) = string
      loop (1) = nrept
      call lookup(string,ntype,ith)
      if(ntype.eq.1 .or. ntype.eq.5) then
       write(jof,510) string
  510  format(' >>>>warning from initst: stack element',a8,' is an'
     &     ,  ' element or an unused label')
      else
       nslot(1) = newsl(1,ith)
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
c-----------------------------------------------------------------------
      return
      end
c
************************************************************************
c
      subroutine lmnt(nt1,nt2,prms,ntrk)
c-----------------------------------------------------------------------
c this routine actually switches to the routines that handle single
c elements, commands, etc.
c
c input: nt1  = group type code of element
c        nt2  = index of element in its group
c        prms = array of params for this element
c        ntrk = 0 if accumulated transfer map is to be constructed
c             = 1 if tracking
c
c  Written by Rob Ryne ca 1984
c  modified by J. Howard 7-87
c              P. Schuett 11-87
c              A. Dragt 6/15/88
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------------------
      include 'map.inc'
      include 'parm.inc'
      include 'parset.inc'
      include 'rays.inc'
      include 'files.inc'
      include 'codes.inc'
      include 'pie.inc'
      include 'infin.inc'
      include 'zeroes.inc'
      include 'frnt.inc'
c--------------------------------
c parameter types and dimensions
c--------------------------------
      integer nt1,nt2,ntrk
      dimension prms(*), reftmp(6)
c-----------------
c local variables
c-----------------
c maps
      double precision h(monoms),ht1(monoms)
      double precision mh(6,6),mht1(6,6)
c
c leading .. and trailing fringe field: lfrn (tfrn) .ne.0, if it is
c to be taken into account
      integer lfrn,tfrn
c Input args nt1,nt2,prms must not be changed. work with kt1,kt2,p.
      dimension p(24)
      character*3 kynd
      equivalence (p1,p(1)),(p2,p(2)),(p3,p(3)),(p4,p(4)),(p5,p(5)),
     #            (p6,p(6))
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c
c Preliminary procedure for random elements or commands
c
      if(nt1.eq.4 .or. nt1.eq.5 .or. nt1.eq.6) then
        nfile = nint(prms(1))
        iecho = nint(prms(2))
        kt1 = nt1 - 3
        kt2 = nt2
        call randin(nfile,kt1,kt2,p)
c
c Echo back if requested
      if( iecho.eq.1 .or. iecho.eq.3) then
        write(jof, 7005) nfile,kt1,kt2
        write(jof ,*)(p(i),i=1,nrp(kt1,kt2))
      endif
      if( iecho.eq.2 .or. iecho.eq.3) then
        write(jodf,7005) nfile,kt1,kt2
        write(jodf,*)(p(i),i=1,nrp(kt1,kt2))
      endif
 7005   format(' random lmnt check:nfile,kt1,kt2=',3i5/
     &         ' parameters found:')
c
c
      else
c
c Preliminary procedure for all other type codes
c
        kt1 = nt1
        kt2 = nt2
      maxp = nrp(nt1,nt2)
      if(maxp.gt.24) maxp = 24
ctm12 write(6,*) 'lmnt: maxp=', maxp
        do 5 i=1,maxp
          p(i)=prms(i)
  5     continue
      endif
c
c Select appropriate action depending on value of kt1,kt2
c
      go to (11,12,13,14,15,16,17,18,19), kt1
c
c     1: simple elements *************************************
c
11    continue
c          'drft    ','nbnd    ','pbnd    ','gbnd    ','prot    ',
      go to(101,       102,       103,       104,       105,
c          'gbdy    ','frng    ','cfbd    ','quad    ','sext    ',
     &      106,       107,       108,       109,       110,
c          'octm    ','octe    ','srfc    ','arot    ','twsm    ',
     &      111,       112,       113,       114,       115,
c          'thlm    ','cplm    ','cfqd    ','dism    ','sol     ',
     &      116,       117,       118,       119,       120,
c          'mark    ','jmap    ','dp      ','recm    ','spce    ',
     &      121,       122,       123,       124,       125,
c          'cfrn    ','coil    ','intg    ','rmap    ','arc     '/
     &      126,       127,       128,       129,       130),kt2
c
c 'drft    ': drift
c
101   call drift(p1,h,mh)
      goto 2000
c
c 'nbnd    ': normal entry bend
c the map for this element is of the form
c gfrngg*nbend*gfrngg with the leading and trailing fringe field
c maps optional
c
102   continue
c     angdeg=p1
      gap=p2
      xk1=p3
      rho=brho/p4
      lfrn=nint(p5)
      tfrn=nint(p6)
c compute nbend
      call nbend(rho,p1,h,mh)
c put on leading fringe field
      if(lfrn.ne.0) then
        call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
        call concat(ht1,mht1,h,mh,h,mh)
      endif
c put on trailing fringe field
      if(tfrn.ne.0) then
        call gfrngg(0.d0,rho,2,ht1,mht1,gap,xk1)
        call concat(h,mh,ht1,mht1,h,mh)
      endif
      goto 2000
c
c  'pbnd    ':
c  parallel faced bend, including leading and trailing
c  pole face rotations and fringe fields.
c  Symmetric bends only are permitted under this option.
c  For parallel-faced magnet with asymmetric entry and exit,
c  use the general bending magnet.
c  The map for the parallel faced bend is of the form
c  prot*gfrng*pbend*gfrng*prot
c
103   continue
      psideg=p1/2.d0
      gap=p2
      xk1=p3
      rho=brho/p4
c compute pbend
      call pbend(rho,p1,h,mh)
c put on the leading fringe field
      call gfrngg(psideg,rho,1,ht1,mht1,gap,xk1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
      call prot(psideg,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on trailing fringe field
      call gfrngg(psideg,rho,2,ht1,mht1,gap,xk1)
      call concat(h,mh,ht1,mht1,h,mh)
c put om trailing prot
      call prot(psideg,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      goto 2000
c
c 'gbnd    ': general bending magnet
c
c  The map for the general bend is of the form
c  prot*gfrng*gbend*gfrng*prot
c
104   continue
      gap=p4
      xk1=p5
      rho=brho/p6
c compute gbend
      call gbend(rho,p1,p2,p3,h,mh)
c put on the leading fringe field
      call gfrngg(p2,rho,1,ht1,mht1,gap,xk1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on leading prot
      call prot(p2,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on trailing fringe field
      call gfrngg(p3,rho,2,ht1,mht1,gap,xk1)
      call concat(h,mh,ht1,mht1,h,mh)
c put on trailing prot
      call prot(p3,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      goto 2000
c
c 'prot    ': rotation of reference plane
c
105   continue
      kind=nint(p2)
      call prot(p1,kind,h,mh)
      goto 2000
c
c 'gbdy    ': body of a general bending magnet
c
106   continue
      rho=brho/p4
      call gbend(rho,p1,p2,p3,h,mh)
      goto 2000
c
c 'frng    ': hard edge dipole fringe fields
c
107   continue
      rho=brho/p4
      iedge=nint(p5)
      gap = p2
      xk1 = p3
      call gfrngg(p1,rho,iedge,h,mh,gap,xk1)
      goto 2000
c
c 'cfbd    ': combined function bend (normal entry and exit)
c
c The map for the combined function bend is of the form
c [(arot*frquad*arotinv)frquad*gfrngg]*
c cfbend*
c [gfrngg*frquad*(arot*frquad*arotinv)]
c with the leading and trailing fringe fields optional.
c The gap size cfbgap and normalized leading and trailing
c field integrals cfblk1 and cfbtk1 for the dipole
c are taken from the table in block common frnt.
c These entries are initialized to 0, .5, .5, respectively
c at the beginning of a Marylie run, and can be changed using
c the type code cfrn.
c The factors of the form (arot*frquad*arotinv) give skew
c quad fringe fields.
c
108   continue
      rho=brho/p2
      lfrn=nint(p3)
      tfrn=nint(p4)
      ijopt=nint(p5)
      iopt=mod(ijopt,10)
c      write(6,*) 'iopt=',iopt
      if ((iopt.lt.1) .or. (iopt.gt.3)) then
      write(jof,*) 'WARNING: parameter iopt outside allowed range in',
     # ' element with type code cfbd'
      endif
      ipset=nint(p6)
c compute cfbend
      call cfbend(p,h,mh)
c put on leading dipole and quad fringe fields
      if(lfrn.ne.0) then
c compute and put on leading dipole fringe field
      gap=cfbgap
      xk1=cfblk1
      call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
c      call nfrng(rho,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c compute and put on leading quad fringe fields
      if((ipset.gt.0) .and. (ipset.le.maxpst)) then
      if(iopt.eq.1) then
      bqd=pst(1,ipset)
      aqd=pst(2,ipset)
      endif
      if(iopt.eq.2) then
      bqd=pst(1,ipset)
      aqd=pst(2,ipset)
      endif
      if(iopt.eq.3) then
      bqd=brho*pst(1,ipset)
      aqd=brho*pst(2,ipset)
      endif
c compute and put on normal quad fringe field
      call frquad(bqd,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c compute and put on skew quad fringe field
      angr=-pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      call frquad(aqd,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      endif
      endif
c put on trailing dipole and quad fringe fields
      if(tfrn.ne.0) then
c compute and put on trailing dipole fringe field
      gap=cfbgap
      xk1=cfbtk1
      call gfrngg(0.d0,rho,1,ht1,mht1,gap,xk1)
c      call nfrng(rho,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c compute and put on trailing quad fringe fields
      if((ipset.gt.0) .and. (ipset.le.maxpst)) then
      if(iopt.eq.1) then
      bqd=pst(1,ipset)
      aqd=pst(2,ipset)
      endif
      if(iopt.eq.2) then
      bqd=pst(1,ipset)
      aqd=pst(2,ipset)
      endif
      if(iopt.eq.3) then
      bqd=brho*pst(1,ipset)
      aqd=brho*pst(2,ipset)
      endif
c compute and put on normal quad fringe field
      call frquad(bqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c compute and put on skew quad fringe field
      angr=pi/4.d0
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      call frquad(aqd,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      angr=-angr
      call arot(angr,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      endif
      endif
      goto 2000
c
c 'quad    ': quadrupole
c
c the map for this element is of the form
c frquad*quad*frquad
c with the leading and trailing fringe field maps optional
c
109   continue
      lfrn=nint(p3)
      tfrn=nint(p4)
c compute the map for the quad
      if(p2.lt.0.d0) call dquad(p1,-p2,h,mh)
      if(p2.eq.0.d0) call drift(p1,h,mh)
      if(p2.gt.0.d0) call fquad(p1,p2,h,mh)
c put on leading fringe field
      if(lfrn.ne.0) then
      call frquad(p2,-1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
      endif
c put on trailing fringe field
      if(tfrn.ne.0) then
      call frquad(p2,1,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
      endif
      goto 2000
c
c 'sext    ': sextupole
c
110   call sext(p1,p2,h,mh)
      goto 2000
c
c 'octm    ': mag. octupole
c
111   call octm(p1,p2,h,mh)
      goto 2000
c
c 'octe    ': elec. octupole
c
112   call octe(p1,p2,h,mh)
      goto 2000
c
c 'srfc    ': short rf cavity
c
113   omega=twopi*p2
      call srfc(p,omega,h,mh)
      goto 2000
c
c 'arot    ': axial rotation
c
114   p1rad=pi180*p1
      call arot(p1rad,h,mh)
      goto 2000
c
c 'twsm    ': linear matrix via twiss parameters
c
115   iplane=nint(p1)
      call twsm(iplane,p2,p3,p4,h,mh)
      goto 2000
c
c 'thlm    ': thin lens low order multipole
c
116   call thnl(p1,p2,p3,p4,p5,p6,h,mh)
      goto 2000
c
c 'cplm    ': "compressed" low order multipole
c
117   call cplm (p,h,mh)
      goto 2000
c
c 'cfqd    ': combined function quadrupole
c
118   call cfdrvr(p,h,mh)
      goto 2000
c
c dispersion matrix (dism)
c
119   call dism(p,h,mh)
      go to 2000
c
c 'sol     ': solenoid
c
120   call gensol(p,h,mh)
      go to 2000
c
c 'mark    ': marker
c
121   continue
      return
c
c 'jmap    ': j mapping
c
122   call jmap(h,mh)
      go to 2000
c
c 'dp      ': data point
c
123   continue
      return
c
c  REC multiplet (recm)
c
124   continue
      call gnrec3(p,h,mh)
      go to 2000
c
c 'spce    ': space
c
125   continue
      return
c
c 'cfrn    ': change or write out fringe field parameters
c for combined function dipole
c
126   continue
      mode=nint(p1)
      if (mode .eq. 0) then
c change fringe field parameters
      cfbgap=p2
      cfblk1=p3
      cfbtk1=p4
      return
      endif
c write out values of fringe field parameters
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of fringe field parameters gap, ennfi, exnfi are:'
      write(jof,*) cfbgap,cfblk1,cfbtk1
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of fringe field parameters gap, ennfi, exnfi are:'
      write(jodf,*) cfbgap,cfblk1,cfbtk1
      endif
      return
c
c     'coil'
c
127   continue
      call coil(prms)
      return
c
c     'intg'
c
128   continue
      call integ(p,h,mh)
      go to 2000
c
c     'rmap'
c
129   continue
      call rmap(p,h,mh)
      go to 2000
c
c     'arc'
c
130   continue
      return
c
c     2: user-supplied elements ************************************
c
c           'usr1    ','usr2    ','usr3    ','usr4    ','usr5    ',
12    go to (201,       202,       203,       204,       205,
c           'usr6    ','usr7    ','usr8    ','usr9    ','usr10   ',
     &       206,       207,       208,       209,       210,
c           'usr11   ','usr12   ','usr13   ','usr14   ','usr15   ',
     &       211,       212,       213,       214,       215,
c           'usr16   ','usr17   ','usr18   ','usr19   ','usr20   ',
     &       216,       217,       218,       219,       220),kt2
c
201   call user1(p)
      return
202   call user2(p)
      return
203   call user3(p)
      return
204   call user4(p)
      return
205   call user5(p)
      return
206   call user6(p,th,tmh)
      return
207   call user7(p,th,tmh)
      return
208   call user8(p,th,tmh)
      return
209   call user9(p,th,tmh)
      return
210   call user10(p,th,tmh)
      return
211   call user11(p,th,tmh)
      return
212   call user12(p,th,tmh)
      return
213   call user13(p,th,tmh)
      return
214   call user14(p,th,tmh)
      return
215   call user15(p,th,tmh)
      return
216   call user16(p,th,tmh)
      return
217   call user17(p,th,tmh)
      return
218   call user18(p,th,tmh)
      return
219   call user19(p,th,tmh)
      return
220   call user20(p,th,tmh)
      return
c
c     3: parameter sets **********************************************
c
c           'ps1     ','ps2     ','ps3     ','ps4     ','ps5     ',
13    go to (301,       302,       303,       304,       305,
c           'ps6     ','ps7     ','ps8     ','ps9     '/
     #       306,       307,       308,       309),kt2
c
301   call pset(p,1)
      return
302   call pset(p,2)
      return
303   call pset(p,3)
      return
304   call pset(p,4)
      return
305   call pset(p,5)
      return
306   call pset(p,6)
      return
307   call pset(p,7)
      return
308   call pset(p,8)
      return
309   call pset(p,9)
      return
c
c     4-6: random elements
c
14    continue
15    continue
16    continue
      write(jof ,9014)
      write(jodf,9014)
 9014 format(' error in lmnt: random element reached')
      call myexit
c
c     7: simple commands *************************************
c
c           'rt      ','sqr     ','symp    ','tmi     ','tmo     ',
17    go to (701,       702,       703,       704,       705,
c           'pmif    ','circ    ','stm     ','gtm     ','end     ',
     &       706,       707,       708,       709,       710,
c           'ptm     ','iden    ','whst    ','inv     ','tran    ',
     &       711,       712,       713,       714,       715,
c           'revf    ','rev     ','mask    ','num     ','rapt    ',
     &       716,       717,       718,       719,       720,
c           'eapt    ','of      ','cf      ','wnd     ','dwnd    ',
     &       721,       722,       723,       724,       725,
c           'ftm     ','wps     ','time    ','cdf     ','bell    ',
     &       726,       727,       728,       729,       730,
c           'wmrt    ','wcl     ','paws    ','inf     ','dims    ',
     &       731,       732,       733,       734,       735,
c           'zer     ','wuca    ','tpol    ','dpol    ','cbm     '/
     &       736,       737,       738,       739,       740),kt2
c
c 'rt      ': ray trace
c
701   icfile=nint(p1)
      nfcfle=nint(p2)
      norder=nint(p3)
      ntrace=nint(p4)
      nwrite=nint(p5)
      iotemp=jof
      jfctmp=jfcf
      jfcf=nfcfle
      if(p6.lt.0.)jof=jodf
      ibrief=iabs(nint(p6))
      call trace(icfile,norder,ntrace,nwrite,th,tmh)
      jof=iotemp
      jfcf=jfctmp
      return
c
c     square the existing map:
c
702   call concat(th,tmh,th,tmh,th,tmh)
      write(jof,5702)
 5702 format(1h ,'existing map concatenated with itself')
      if(ntrk.eq.1)write(jof,9702)
 9702 format(1h ,'warning: map squared in turtle mode')
      return
c
c     symplectify matrix in transfer map
c
703   continue
      iopt=nint(p1)
      kind=nint(p2)
      call sympl(iopt,kind,th,tmh)
      return
c
c     input transfer map from an external file:
c
704   continue
      iopt=nint(p1)
      ifile=nint(p2)
      nopt=nint(p3)
      nskp=nint(p4)
c    rewind only option
      if (nopt.eq.1 .and. nskp.eq.-1) then
       rewind ifile
       write(jof,5704) ifile
 5704  format(1x,'file unit ',i3,' rewound')
       return
      endif
c
c    other options
c
      mpitmp=mpi
      mpi=ifile
      call mapin(nopt,nskp,h,mh)
      mpi=mpitmp
c  option when tracking (procedure at end)
      if(ntrk.eq.1) goto 2000
c  options when not tracking
      if(iopt.eq.1) call concat(th,tmh,h,mh,th,tmh)
      if(iopt.eq.2) call mapmap(h,mh,th,tmh)
      return
c
c     output transfer map to an external file (tmo):
c
705   ifile=nint(p1)
      mpotmp=mpo
      mpo=ifile
      call mapout(0,th,tmh)
      mpo=mpotmp
      return
c
c     print contents of file master input file:
c
706   itype=nint(p1)
      ifile=nint(p2)
      isend=nint(p3)
      jtmp=jodf
      jodf=ifile
      if(isend.eq.1.or.isend.eq.3)call pmif(jof,itype)
      if(isend.eq.2.or.isend.eq.3)call pmif(jodf,itype)
      jodf=jtmp
      return
c
c     Note: the program should never get here.
c     (cqlate is called directly from tran)
c
707   write(jof ,9707)
      write(jodf,9707)
 9707 format(' error: reached element "circ" in routine lmnt')
      call myexit
c
c     store the existing transfer map
c
708   continue
      kynd='stm'
      nmap=nint(p1)
      if ((nmap.gt.5).or.(nmap.lt.1)) then
        write(jof,9708) nmap
 9708   format(1x,'nmap=',i3,1x,'trouble with stm:nmap < 1 or > 5')
        call myexit
      else
        call strget(kynd,nmap,th,tmh)
        return
      endif
c
c     get transfer map from storage
c
709   continue
      kynd='gtm'
      iopt=nint(p1)
      nmap=nint(p2)
c
      if ((nmap.gt.5).or.(nmap.lt.1)) then
        write(jof,9709) nmap
 9709   format(1x,'nmap=',i3,1x,'trouble with gtm:nmap < 1 or > 5')
        call myexit
        return
      endif
c
      call strget(kynd,nmap,h,mh)
c option when tracking(procedure at end)
      if(ntrk .eq. 1) goto 2000
c options when not tracking
      if(iopt.eq.1) call concat(th,tmh,h,mh,th,tmh)
      if(iopt.eq.2) call mapmap(h,mh,th,tmh)
      return
c
c     end of job:
c
710   write(jof,5710)
 5710 format(/1x,'end of MARYLIE run')
      call myexit
      return
c
c     print transfer map:
c
711   continue
      n1=nint(p1)
      n2=nint(p2)
      n3=nint(p3)
      n4=nint(p4)
      n5=nint(p5)
      if(n5.eq.1)call pcmap(n1,n2,n3,n4,th,tmh)
      if(n5.eq.2)call psrmap(n1,n2,th,tmh)
      if(n5.eq.3)call pdrmap(n1,n2,th,tmh)
      return
c
c     identity mapping:
c
712   call ident(th,tmh)
      return
c
c     write history of beam loss
c
713   call whst(p)
      goto 2000
c
c     inverse:
c
714   call inv(th,tmh)
      return
c
c     transpose:
c
715   call mtran(tmh)
      return
c
c     reverse factorization:
c
716   iord=nint(p1)
      call revf(iord,th,tmh)
      return
c
c     Dragt's reversal
c
717   call rev(th,tmh)
      return
c
c     mask off selected portions of transfer map:
c
718   call mask(p,th,tmh)
      return
c
c     number lines in a file
c
719   call num(p)
      return
c
c     aperture particle distribution
c
720   call rapt(p)
      return
c
721   continue
      mode=nint(p1)
      if (mode .eq. 1) call eapt(p)
      return
c
c     open files
c
722   call of(p)
      return
c
c     close files
c
723   call cf(p)
      return
c
c     window particle distribution
c
724   call wnd(p)
      return
725   call dwnd(p)
      return
c
c     filter transfer map
c
726   call ftm(p,th,tmh)
      return
c
c     write parameter set
c
727   ipset = nint(p(1))
      isend = nint(p(2))
      call wps(ipset,isend)
      return
c
c     write time
c
728   call mytime(p)
      return
c
c     change output drop file
c
729   jodf = nint(p(1))
      return
c
c     ring bell
c
730   continue
      call bell
      return
c
c     write value of merit function
c
731   ifn = nint(p(1))
      isend = nint(p(2))
      call wmrt(ifn,isend)
      return
c
c     write contents of loop
c
732   call wcl(p)
      return
c
c     pause (paws)
c
733   continue
      write(jof,*) ' press return to continue'
      read(5,7330) iwxyz
7330  format(a1)
      return
c
c     change or write out infinities (inf)
c
734   continue
      mode=nint(p1)
      if (mode .eq. 0) then
c change infinities
      xinf=p2
      yinf=p3
      tinf=p4
      ginf=p5
      return
      endif
c write out values of infinities
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of infinities xinf, yinf, tinf, ginf are:',
     & xinf,yinf,tinf,ginf
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of infinities xinf, yinf, tinf, ginf are:',
     & xinf,yinf,tinf,ginf
      endif
      return
c
c get dimensions (dims)
c
735   continue
      return
c
c     change or write out values of zeroes (zer)
c
736   continue
      mode=nint(p1)
      if (mode .eq. 0) then
c change values of zeroes
      fzer=p2
      detz=p3
      return
      endif
c write out values of zeroes
      isend=mode
      if( isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*)
     & 'values of zeroes fzero, detzero are:',
     & fzer,detz
      endif
      if( isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*)
     & 'values of zeroes fzero, detzero are:',
     & fzer,detz
      endif
      return
c
ctm Feb 2012      dump ucalc buffer
c
737   continue
      call wuca(p)
      return
c
c     twiss polynomial (tpol)
c
738   continue
      call tpol(p,th,tmh)
      return
c
c     dispersion polynomial (dpol)
c
739   continue
      call dpol(p,th,tmh)
      return
c
c     change or write out beam parameters (cbm)
c
740   continue
      job=nint(p1)
      if (job .eq. -1) then
      prms(1)=0.d0
      prms(2)=brho
      prms(3)=gamm1
      endif
      if (job .eq. 0) then
      brho=p2
      gamm1=p3
c  recomputation of relativistic beta and gamma factors:
      gamma=gamm1+1.d0
      stuff2=gamm1*(gamma+1.d0)
      stuff1=sqrt(stuff2)
      beta=stuff1/gamma
      endif
      if ((job .eq. 1) .or. (job .eq. 3)) then
      write(jof,*) ' beam parameters are: ', brho,gamm1
      endif
      if ((job .eq. 2) .or. (job .eq. 3)) then
      write(jodf,*) ' beam parameters are: ', brho,gamm1
      endif
      return
c
c
c     8: advanced commands ********************************
c
c           'cod     ','amap    ','dia     ','dnor    ','exp     ',
18    go to (801,       802,       803,       804,       805,
c           'pdnf    ','psnf    ','radm    ','rasm    ','sia     ',
     &       806,       807,       808,       809,       810,
c           'snor    ','tadm    ','tasm    ','tbas    ','gbuf    ',
     &       811,       812,       813,       814,       815,
c           'trsa    ','trda    ','smul    ','padd    ','pmul    ',
     &       816,       817,       818,       819,       820,
c           'pb      ','pold    ','pval    ','fasm    ','fadm    ',
     &       821,       822,       823,       824,       825,
c           'sq      ','wsq     ','ctr     ','asni    ','pnlp    ',
     &       826,       827,       828,       829,       830,
c           'csym    ','psp    ','mn       ','bgen    ','tic     ',
     &       831,       832,      833,        834,       835,
c           'ppa     ','moma   ','geom     ','fwa     '/
     &       836,       837,      838,        839),kt2
c
c     off-momentum closed orbit analysis
c
801   call cod(p,th,tmh)
      return
c
c     apply map to a function or moments
c
802   call amap(p,th,tmh)
      return
c
c     dynamic invariant analysis
c
803   call dia(p,th,tmh)
      return
c
c     dynamic normal form analysis
c
804   call dnor(p,th,tmh)
      return
c
c     compute exponential
c
805   call cex(p,th,tmh)
      return
c
c     compute power of dynamic normal form
c
806   call pdnf(p,th,tmh)
      return
c
c     compute power of static normal form
c
807   call psnf(p,th,tmh)
      return
c
c     resonance analyze dynamic map
c
808   call radm(p,th,tmh)
      return
c
c     resonance analyze static map
c
809   call rasm(p,th,tmh)
      return
c
c     static invariant ayalysis
c
810   call sia(p,th,tmh)
      return
c
c     static normal form analysis
c
811   call snor(p,th,tmh)
      return
c
c     twiss analyze dynamic map
c
812   call tadm(p,th,tmh)
      return
c
c     twiss analyze static map
c
813   call tasm(p,th,tmh)
      return
c
c     translate basis
c
814   call tbas(p,th,tmh)
      return
c
c     get buffer contents
c
815   continue
c
c test control parameters
c
      nmap=nint(p2)
      if (nmap.gt.5 .or. nmap.lt.1) then
        write(jof,9815) nmap
 9815   format(1x,'nmap=',i3,1x,'trouble with gbuf:nmap < 1 or > 5')
        call myexit
      endif
c option when tracking (procedure at end)
      if(ntrk .eq. 1) then
        p1temp=p(1)
        p(1)=2
        call gbuf(p,h,mh)
        p(1)=p1temp
        goto 2000
      endif
c options when not tracking
      call gbuf(p,th,tmh)
        return
c
c     transport static (script) A
c
816   call trsa(p,th,tmh)
      return
c
c    transport dynamic script A
c
817   call trda(p,th,tmh)
      return
c
c     multiply polynomial by a scalar
c
818   call smul(p,th,tmh)
      return
c
c     add two polynomials
c
819   call padd(p,th,tmh)
      return
c
c     multiply two polynomials
c
820   call pmul(p,th,tmh)
      return
c
c     Poisson bracket two polynomials
c
821   call pbpol(p,th,tmh)
      return
c
c     polar decompose matrix portion of transfer map
c
822   call pold(p,th,tmh)
      return
c
c     evaluate a polynomial
c
823   call pval(p,th,tmh)
      return
c
c     fourier analyze static map
c
824   call fasm(p,th,tmh)
      return
c
c     fourier analyze dynamic map
c
825   call fadm(p,th,tmh)
      return
c
c     select quantities
c
826   call sq(p)
      return
c
c     write selected quantities
c
827   call wsq(p)
      return
c
c     change tune ranges
c
828   continue
      call subctr(p)
      return
c
c     apply script N inverse
c
829   continue
      call asni(p)
      return
c
c     compute power of nonlinear part
c
830   continue
      call pnlp(p,th,tmh)
      return
c
c     check for symplecticity
c
831   continue
      isend=nint(p1)
      call csym(isend,tmh,ans)
      return
c
c     (psp) compute scalar product of two polynomials
c
832   call psp(p,th,tmh)
      return
c
c     (mn) compute matrix norm
c
833   call submn(p,th,tmh)
      return
c
c     (bgen) generate a beam
c
834   call bgen(p)
      return
c
c     (tic) translate (move) initial conditions
c
835   call tic(p)
      return
c
c     (ppa) principal planes analysis
c
836   call ppa(p,th,tmh)
      return
c
c     (moma) moment and map analysis
c
837   call moma(p)
      return
c
c     (geom) compute geometry of a loop
c
838   call geom(p)
      return
c
c     (fwa) copy file to working array
c
839   call fwa(p)
      return
c
c     9: procedures and fitting and optimization *************************
c
c           'bip     ','bop     ','tip     ','top     ',
19    go to (901,       902,       903,       904,
c           'aim     ','vary    ','fit     ','opt     ',
     &       905,       906,       907,       908,
c           'con1    ','con2    ','con3    ','con4    ','con5    ',
     &       909,       910,       911,       912,       913,
c           'mrt0    ',
     &       914,
c           'mrt1    ','mrt2    ','mrt3    ','mrt4    ','mrt5    ',
     &       915,       916,       917,       918,       919,
c           'fps     ',
     &       920,
c           'cps1    ','cps2    ','cps3    ','cps4    ','cps5    ',
     &       921,       922,       923,       924,       925,
c           'cps6    ','cps7    ','cps8    ','cps9    ',
     &       926,       927,       928,       929,
c           'dapt    ','grad    ','rset    ','flag    ','scan    ',
     &       930,       931,       932,       933,       934,
c           'mss     ','gendip  ','taylor  '/
     &       935,       936,       937),kt2
c     begin procedures
c
901   call bip(p)
      return
902   call bop(p)
      return
c
c     end procedures
c
903   call subtip(p)
      return
904   call subtop(p)
      return
c
c     specify aims
c
905   call aim(p)
      return
c
c     specify quantities to be varied
c
906   call vary(p)
      return
c
c     fit to achieve aims
c
907   call fit(p)
      return
c
c     optimize
c
908   call opt(p)
      return
c
c     constraints
c
909   call con1(p)
      return
910   call con2(p)
      return
911   call con3(p)
      return
912   call con4(p)
      return
913   call con5(p)
      return
c
c     merit functions
c
c     least squares merit function
c
914   call mrt0
      return
c
c     user supplied merit functions
c
915   call mrt1(p)
      return
916   call mrt2(p)
      return
917   call mrt3(p)
      return
918   call mrt4(p)
      return
919   call mrt5(p)
      return
c
c     free parameter sets
c
920   ipset=nint(p1)
      call fps(ipset)
      return
c
c     capture parameter sets
c
921   ipset=1
      call cps(prms,p,ipset)
      return
922   ipset=2
      call cps(prms,p,ipset)
      return
923   ipset=3
      call cps(prms,p,ipset)
      return
924   ipset=4
      call cps(prms,p,ipset)
      return
925   ipset=5
      call cps(prms,p,ipset)
      return
926   ipset=6
      call cps(prms,p,ipset)
      return
927   ipset=7
      call cps(prms,p,ipset)
      return
928   ipset=8
      call cps(prms,p,ipset)
      return
929   ipset=9
      call cps(prms,p,ipset)
      return
c
c     compute dynamic aperture (dapt)
c
930   continue
      call dapt(p)
      return
c
c     gradient (grad)
c
931   continue
      call grad(p)
      return
c
c     rset
c
932   continue
      call rset(p)
      return
c
c     flag
c
933   continue
      call flag(p)
      return
c
c     scan
c
934   continue
      call scan(p)
      return
c
c     mss
c
935   continue
      call mss(p)
      return
c
c spare1 is now gendip (ctm Feb 2012)
c
936   continue
      call gendip(p,h,mh,reftmp)
c put on leading prot
      call prot(p2,1,ht1,mht1)
      call concat(ht1,mht1,h,mh,h,mh)
c put on trailing prot
      p3neg=p3
      call prot(p3neg,2,ht1,mht1)
      call concat(h,mh,ht1,mht1,h,mh)
c
cryne hardwired includedrifts:
      includedrifts=1
      if(includedrifts.eq.1)then
c put on leading negative drift
        z1neg=-p(5)*p(6)/abs(cos(p(2)*pi180))
cryne   if(idproc.eq.0.and.iprintmsgs.eq.1)write(6,*)'z1neg=',z1neg
        call drift(z1neg,ht1,mht1)   !drift3 in ML/I
        call concat(ht1,mht1,h,mh,h,mh)
c put on trailing negative drift
        z2neg=-p(5)*p(6)/abs(cos(p(3)*pi180))
cryne   if(idproc.eq.0.and.iprintmsgs.eq.1)write(6,*)'z2neg=',z2neg
        call drift(z2neg,ht1,mht1)   !drift3 in ML/I
        call concat(h,mh,ht1,mht1,h,mh)
      endif
      goto 2000
c
c     spare2 is now taylor. ctm feb2012
c
937   continue
      call tugen(th,tmh)
      return
c
c  ......... concatenate or track before exiting ........
c
 2000 if(ntrk.eq.0)then
        call concat(th,tmh,h,mh,th,tmh)
        return
      else
        call trace(0,5,1,0,h,mh)
        return
      endif
      end
c
************************************************************************
c
      subroutine lookup(string,itype,index)
c-----------------------------------------------------------------------
c  this subroutine determines whether the input string 'string'
c  is an element,line,lump,loop, or unused label.
c
c  input:  string  character*8   item name
c  output: itype   integer       =1, if string is an element
c                                =2, .... line
c                                =3, .... lump
c                                =4, .... loop
c                                =5, .... unused label
c          index   integer       index of 'string' in its array
c
c  Written by Rob Ryne ca 1984
c  rewritten by Petra Schuett
c             October 30, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'elmnts.inc'
      include 'items.inc'
      include 'mldex.inc'
c-----------------------------------------------------------------------
c parameter types
c-----------------------------------------------------------------------
      character*8 string
      integer itype,index
c-----------------------------------------------------------------------
c start routine
c-----------------------------------------------------------------------
c element?
      do 10 n=1,na
      if(string.eq.lmnlbl(n)) then
       itype = 1
       index = n
       return
      endif
   10 continue
c item?
      do 20 n=1,nb
      if(string.eq.ilbl(n)) then
       itype = ityp(n)
       index = n
       return
      endif
   20 continue
c not found:
      itype=5
      return
      end
c
***********************************************************************
c
      subroutine low(line)
c  Converts all uppercase characters to lowercase.
c  Written by Liam Healy, Feb. 28, 1985.
c
      character*80 line
c
      character*26 lower,upper
      character*1 blank
      save lower,upper,blank
      data lower/'abcdefghijklmnopqrstuvwxyz'/
      data upper/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      data blank/' '/
c
      lnbc=1
      do 100 i=1,80
        loc=index(upper,line(i:i))
        if(loc.gt.0) then
          line(i:i)=lower(loc:loc)
        endif
        if (line(i:i).ne.blank) lnbc=i
  100 continue
      return
      end
c
***********************************************************************
      subroutine lumpit(mth)
c-----------------------------------------------------------------------
c  translate map of mth lump to special format and save it in core
c
c  input mth (integer) : index of lump in /items/
c
c  Written by Rob Ryne ca 1984
c  changed by Petra Schuett
c             October 30,1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'items.inc'
      include 'map.inc'
      include 'deriv.inc'
      include 'core.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
c inew keeps track of the lump to be deleted next, if core is full
      integer inew
      save inew
      data inew /0/
c-----------------------------------------------------------------------
c  start routine
c-----------------------------------------------------------------------
c find vacant space in core
      icore = 0
      do 1 i=maxlum,1,-1
       if(inuse(i).eq.0) icore=i
   1  continue
c if core is full, destroy oldest lump
      if(icore.eq.0) then
       inew=inew+1
       if(inew.eq.maxlum)inew=1
       icore=inew
       lmade(inuse(icore))=0
       write(jof,510) ilbl(inuse(icore))
 510   format(1h ,'lump ',a8,' deleted;')
      endif
c--------------------
c calculate rjac,df,..
      call canx(tmh,th)
c ..rrjac and rdf
      call rearr
c--------------------
c now store all the info
      do 10 n1=1,6
      do 10 n2=1,83
      dfl(n1,n2,icore)=df(n1,n2)
   10 continue
      do 20 n1=1,3
      do 20 n2=1,84
   20 rdfl(n1,n2,icore)=rdf(n1,n2)
      do 30 n1=1,3
      do 30 n2=1,3
      do 30 n3=1,28
   30 rrjacl(n1,n2,n3,icore)=rrjac(n1,n2,n3)
      do 40 n1=1,6
      do 40 n2=1,6
   40 tmhl(n1,n2,icore)=tmh(n1,n2)
      do 50 n1=1,monoms
   50 thl(n1,icore)=th(n1)
c--------------------
c set pointers
      inuse(icore) = mth
      lmade(mth)   = icore
c--------------------
      write(jof,520) ilbl(mth),icore
      write(jodf,520)  ilbl(mth),icore
  520 format(1h ,'lump ',a8,' constructed and stored.','(',i2,')')
c--------------------
      return
      end
c
************************************************************************
c
      function newsl(mp,kth)
c-----------------------------------------------------------------------
c newsl is the first icon to be treated. Depending on the sign of the
c repetition factor loop(mp), it is either the first or the last one.
c
c Petra Schuett, November 6,1987
c-----------------------------------------------------------------------
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'items.inc'
c-------
c start
c-------
      if(loop(mp).ge.0) then
        newsl = 1
      else
        newsl = ilen(kth)
      endif
      return
      end
************************************************************************
      function nsign(number)
c Petra Schuett, November 6,1987
c
      if (number .ge. 0) then
        nsign = 1
      else
        nsign = -1
      endif
      return
      end
************************************************************************
      subroutine pop(lempty)
c-----------------------------------------------------------------------
c  pop stack
c
c  output: lempty = .true. if stack is empty
c
c  Petra Schuett  October 30,1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'items.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      logical lempty
c-------
c start
c-------
      np = np - 1
      if(np .le. 0) then
       lempty = .true.
      else
       lempty = .false.
       call lookup(lstac(np),ntype,ith)
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
      return
      end
************************************************************************
      subroutine push
c-----------------------------------------------------------------------
c  push stack: add actual icon on top of it
c
c  Petra Schuett  October 30,1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'items.inc'
      include 'files.inc'
c-------
c start
c-------
      np = np + 1
      if(np .gt. mstack) then
       write(jof ,510)
       write(jodf,510)
 510   format(' error in push: stack overflow')
       call myexit
      else
       lstac(np)   = icon(nslot(np-1),ith)
       loop(np)    = irep(nslot(np-1),ith) * nsign(loop(np-1))
       nslot(np)   = newsl(np,jth)
       nslot(np-1) = nslot(np-1) + nsign(loop(np-1))
       ith         = jth
       ntype       = mtype
       call lookup(icon(nslot(np),ith),mtype,jth)
      endif
      return
      end
************************************************************************
      subroutine readin(line,leof)
c-----------------------------------------------------------------------
c  Reads a line from file lf, puts it in the character variable 'line'.
c  Designed so that filters may put on, such as the routine to convert
c  all uppercase characters to lower case.
c  Written by Liam Healy, May 1, 1985.
c
c  Changed by Petra Schuett, October 19, 1987 :
c    use logical: leof=.true. , if end of file is encountered
c-----------------------------------------------------------------------
      include 'files.inc'
      character*80 line
      logical leof
c----Routine----
      read(lf,800,end=100) line
  800 format(a)
c      call low(line)
      return
  100 leof=.true.
      return
      end
c
************************************************************************
      subroutine rearec(line,leof,msegm,strarr,narr,itot,lcont)
c-----------------------------------------------------------------------
c  This routine reads a new line and interprets it partly:
c   - if a line starts with '!' (comment line), it is skipped and
c                  another line is read
c   - if a #... code is read, it sets msegm, indicating the start of
c                  a new input component (segment).
c   - lines in the comment component (segment) are not interpreted.
c   - all other lines are cut into stings and numbers assuming a form:
c        n1*str1 n2*str2 ... n3*str3 n4*str4,n5*str5  (etc)
c        the numbers n1,n2,... are optional
c        the strings and numbers are separated by blanks and/or commas
c   - a line ending with & is continued on the next line (lcont=true)
c   - anything after a '!' is ignored
c
c  Input: msegm      integer      input component (segment) currently
c                                 being read
c
c  Output:line       character*80 line read
c         leof       logical      =.true. if end of file is encountered
c         msegm      integer      future component (segment) to be read
c         strarr(40) character*8  array of strings str1,str2,...
c         narr(40)   integer      array of numbers n1,n2,...
c         itot       integer      number of strings found
c         lcont      logical      =.true. if line to be continued
c
c  Author: Petra Schuett
c          October 19, 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
c
      include 'files.inc'
      include 'codes.inc'
c
      integer msegm,narr(40),itot
      character*8  strarr(40)
      character*80 line
      logical leof,lcont
c
      character*10 string
      logical lnum,lfound
c-----------------
c init
      lcont=.false.
c read new line
  1   call readin(line,leof)
c end of file
      if (leof) return
c convert to lower case, except for comment component (segment)
      if (msegm.ne.1) call low(line)
c find first string
      kbeg=1
      call cread(kbeg,msegm,line,string,lfound)
c empty line is ignored, except in comment component (segment)
      if(.not.lfound) then
        if(msegm.eq.1) then
         return
        else
         goto 1
        endif
      endif
c new component (segment) starts...
      if(string(1:1).eq.'#')then
c in case previous component (segment) was comment, conv to lower case
        call low(line)
c ...which one?
        do 2 k=1,8
          if(string(1:8).eq.ling(k)) then
            msegm = k
c            write(jodf,*) msegm
            return
          endif
  2     continue
c ... no match:
c default is #beam after #comment
        if(msegm.eq.1) then
          msegm = 2
          return
        else
c but in all other cases, this should not happen!
          write(jof,99) string
  99      format(' ---> warning from rearec:'/
     &           '      user name ',a10,' begins with #')
        endif
      endif
c if #comment line is read, no interpretation
      if(msegm.eq.1) return
c comment line
      if (string(1:1).eq.'!') goto 1
c.........................................................
c now interpret line
c first init itot
      itot = 0
      do 10 i=1,40
c end of line
      if((.not.lfound).or.(string(1:1).eq.'!')) return
c line to be continued
      if(string(1:1).eq.'&') then
        lcont = .true.
        return
      endif
c so we found another string
      itot=itot+1
c is it a number?
      call cnumb(string,num,lnum)
c      write(jodf,*)string,'=',num,'lnum=',lnum
      if((.not.lnum).and.(string(1:1).ne.'-')) then
c.. this must be "string"
        narr(i)=1
        strarr(i)=string(1:8)
      else if(.not.lnum) then
c.. it is "-string"
        narr(i)=-1
        strarr(i)=string(2:9)
      else
c.. it is the number of "n*string"
        narr(i)=num
        call cread(kbeg,msegm,line,string,lfound)
c        write(jodf,*)kbeg,string,lfound
        if (.not.lfound) then
c.. error
          write(jof ,98) line
          write(jodf,98) line
  98      format(' ---> warning from rearec:',/,
     &           '      the following line contains a number which is'
     &          ,' not followed by a string:',/,'      ',a80/)
          call myexit
        endif
c.. normal way
        strarr(i)=string(1:8)
      endif
c find next string and start over again
      call cread(kbeg,msegm,line,string,lfound)
c      write(jodf,*)kbeg,string,lfound
  10  continue
c
c more than 40 strings, which are separated by single characters,
c cannot occur in a line of 80 characters.
      end
c
c***********************************************************************
c
      subroutine tran
c-----------------------------------------------------------------------
c organizes translation of input into work to be done
c
c Written by Rob Ryne ca 1984
c adapted to new version 9 Nov 87 by
c Petra Schuett
c Modified 31 Aug 88 by Alex Dragt
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c-----------------------------------------------------------------------
c common blocks
c-----------------------------------------------------------------------
      include 'map.inc'
      include 'bmline.inc'
      include 'mldex.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'core.inc'
      include 'loop.inc'
      include 'files.inc'
      include 'labpnt.inc'
c-----------------------------------------------------------------------
c local variables
c-----------------------------------------------------------------------
      dimension rh(monoms),rmh(6,6)
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c initialize
      call ident(th,tmh)
      jicnt=0
      jocnt=0
c main loop through latt (labor)
ctm 3/01  EXPLICIT DO LOOP REMOVED for other compilers
ctm 3/01       do 5 lp=1,noble
      lp = 1
  10  continue
        call lookup(latt(lp),ntype,ith)
        if (ntype.eq.1) then
c  command 'circ' ...
          if(nt1(ith).eq.7 .and. nt2(ith).eq.7 ) then
            icfile=nint(pmenu(1+mpp(ith)))
            nfcfle=nint(pmenu(2+mpp(ith)))
            norder=nint(pmenu(3+mpp(ith)))
            ntimes=nint(pmenu(4+mpp(ith)))
            nwrite=nint(pmenu(5+mpp(ith)))
            isend =nint(pmenu(6+mpp(ith)))
            jfctmp=jfcf
            jfcf  =nfcfle
            call cqlate(icfile,norder,ntimes,nwrite,isend)
            jfcj=jfctmp
          else
c  ... or other element/command
            call trlmnt(ith,num(lp))
          endif
        else if (ntype.eq.2) then
c  line
          call trobj(latt(lp),num(lp),0)
        else if (ntype.eq.3) then
c  lump
c
c  ignore or destroy lump with zero repetition number
        if( num(lp).eq.0) then
c  if lump is unmade, ignore it
        if( lmade(ith).eq.0 ) then
            write(jof, 510) ilbl(ith)
            write(jodf,510) ilbl(ith)
 510        format(1x,'unmade lump ',a8,
     #      ' in #labor with rep no. = 0 ignored')
        else
c  if lump is made, destroy it
          do 30 k = 1,maxlum
           if(inuse(k).eq.ith) then
            inuse(k) = 0
            write(jof, 520) ilbl(ith)
            write(jodf,520) ilbl(ith)
 520        format(1x,'lump ',a8,
     #      ' in #labor with rep n. = 0 destroyed')
           endif
 30       continue
          lmade(ith) = 0
        endif
        endif
c
c   otherwise
        if( num(lp).ne.0) then
          call trobj(latt(lp),num(lp),0)
        endif
c
        else if (ntype.eq.4) then
c  loop
          if( num(lp).ne.0) nloop=ith
          call trloop(ilbl(ith),num(lp))
c  store unmade lumps and replace lump-number by number in core
          do 40 i=1,joy
            if(mim(i).ge.0.and.mim(i).le.5000) then
              jth=mim(i)
              if(lmade(jth).eq.0) then
                call mapmap(th,tmh,rh,rmh)
                call ident(th,tmh)
                call trobj(ilbl(jth),1,0)
c               call lumpit(jth)
                call mapmap(rh,rmh,th,tmh)
              endif
              mim(i)=lmade(jth)
            endif
  40      continue
        else
c unused label
          write(jof,610) latt(lp)
  610     format(1h ,'warning from tran: ',a8,' not found.')
        endif
   5  continue
ctm 3/01  remove explicit do loop
      lp = lp + 1
      if(lp.le.noble) go to 10
      return
      end
c
c***********************************************************************
c
      subroutine trlmnt(ith,irep)
c-----------------------------------------------------------------------
c handle single item (element/command) in the menu
c
c input: ith  index of the item in menu
c        irep repetition factor
c
c Petra Schuett, Nov.9,1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c
      include 'elmnts.inc'
c------ Put ith in inmenu common variable ------------
      inmenu = ith
c----------
      do 10 i = 1,iabs(irep)
        call lmnt(nt1(ith),nt2(ith),pmenu(1+mpp(ith)),0)
 10   continue
      return
      end
c
c***********************************************************************
c
      subroutine trloop(string,nrept)
c-----------------------------------------------------------------------
c     interprets loops
c     Written by Rob Ryne ca 1984
c     adapted to use of strings and slightly simplified in logic
c     by Petra Schuett  11-10-87
c     Fixed by Filippo Neri and Alex Dragt 8-29-88
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'files.inc'
      include 'loop.inc'
c---------------
c parameter type
c---------------
      character*8 string
c-----------------
c local variables
c-----------------
c lempty = .true. , when stack is empty
      logical lempty
      save lempty
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c???????????????????????????????????????????????????????????????????????
c       printout counters: prints out first nm icons
c
c      nm = 100
c      nc = 0
c      write(jof, 700)
c      write(jodf,700)
c 700  format(/' entering trloop:'/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c procedure when nrept=0
c
      if(nrept.eq.0) then
          write(jodf,508) string
          write(jof, 508) string
 508      format(1x,'loop ',a8,' with rep no. = 0 has been ignored')
      return
      endif
c
c procedure when nrept .ne. 0
c
      lnrept=nsign(nrept)
      joy=1
c  initialize the stacks
      call initst(string,lnrept)
c
c-----------------------------------------------------------------------
c here we start with a stack element, either a new one or one that
c just has been popped
c
 1000 continue
c
      if(ntype .eq. 1) then
c menu entry should never occur here ...
        write(jodf,910) lstac(np)
        write(jof ,910) lstac(np)
  910   format(1x,'error in trloop: menu entry ',a8,
     &            ' found at start of routine.')
        call myexit
      else if (ntype.eq.3) then
c ... neither should a lump
        write(jodf,920) lstac(np)
        write(jof ,920) lstac(np)
  920   format(1x,'error in trloop: lump ',a8,
     &            ' found at start of routine.')
        call myexit
      else if (ntype.eq.5) then
c unknown label is ignored
        write(jodf,930) lstac(np)
        write(jof, 930) lstac(np)
  930   format(1x,'warning from trloop: object ',a8,' not found.')
c
c main part is loops or lines (should be the only part used)
c
      else if (ntype.eq.2 .or. ntype.eq.4) then
c
c here we begin a computational loop, handling simple
c entries in lines or loops
c
 2000   continue
c
        if((ntype.eq.2 .or. ntype.eq.4) .and.
     &     (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a line or loop
          loop(np) = loop(np) - nsign(loop(np))
          if(loop(np).ne.0) then
            nslot(np) = newsl(np,ith)
            call lookup(icon(nslot(np),ith),mtype,jth)
            goto 1000
          endif
        else if (mtype.eq.2 .and. irep(nslot(np),ith).eq.0) then
c if  a line and rep no. = 0 ignore line
          write(jodf,509) icon(nslot(np),ith)
          write(jof, 509) icon(nslot(np),ith)
 509      format(1x,'line ',a8,' with rep no. = 0 has been ignored')
          nslot(np) = nslot(np) + nsign(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.4 .and. irep(nslot(np),ith).eq.0) then
c if a loop and rep no. = 0 ignore loop
          write(jodf,510) icon(nslot(np),ith)
          write(jof, 510) icon(nslot(np),ith)
 510      format(1x,'loop ',a8,' with rep no. = 0 has been ignored')
          nslot(np) = nslot(np) + nsign(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.5) then
c ignore unknown label
          write(jodf,930) icon(nslot(np),ith)
          write(jof, 930) icon(nslot(np),ith)
          nslot(np) = nslot(np) + nsign(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.1 .or. mtype.eq.3) then
c line or loop points to a menu entry or lump
          do 10 i=1,iabs(irep(nslot(np),ith))
            if(mtype.eq.1 .and. nt1(jth).eq.2) then
c user supplied element
              if(nt2(jth).gt.5) then
                write(jof ,940) icon(nslot(np),ith)
                write(jodf,940) icon(nslot(np),ith)
 940            format(/' warning from trloop: ',a8,' is a usern',
     &                  ' with n>5')
              endif
              mim(joy) = 5000+jth
            elseif (mtype.eq.1) then
c other menu entry
              mim(joy)=-jth
            elseif (mtype.eq.3) then
c lump
              mim(joy)=jth
            endif
            joy=joy+1
  10      continue
          if(joy.gt.joymax)then
            write(jodf,950) joymax
            write(jof ,950) joymax
 950        format(1x,'error: array length >= joymax (',
     #      i6,') in trloop')
            call myexit
          endif
c - next item:
          nslot(np) = nslot(np) + nsign(loop(np))
          if (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith)) goto 2000
          call lookup(icon(nslot(np),ith),mtype,jth)
          goto 2000
        else if (mtype.eq.2 .or. mtype.eq.4) then
c line or loop points to line or loop
          call push
          goto 1000
        endif
c
c end of interpretation of a line/loop
c
      endif
c end of main way through stack-element
c-----------------------------------------------------------------------
c     pop stack; see what's there. this is the normal way to return
c
      call pop(lempty)
      if (lempty) then
        joy = joy-1
        return
      endif
      goto 1000
c
      end
c
c***********************************************************************
c
      subroutine trobj(string,nrept,ntrk)
c-----------------------------------------------------------------------
c     interprets lines and lumps
c  Written by Rob Ryne ca 1984
c        Modified to store unmade lumps (tro5)
c
c        New features:
c          1. Stores and retrieves unmade lumps in lines
c          2. Stores and retrieves nested lumps
c          3. Ignores lumps with nrept = 0
c          4. Ignores lines with nrept = 0
c
c        Jim Howard   CDG   7-7-87
c
c    adapted to use of strings and slightly simplified in logic
c        Petra Schuett  11-9-87
c
c    modified by Alex Dragt 31 August 1988
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
c--------
c commons
c--------
      include 'stack.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'core.inc'
      include 'map.inc'
      include 'files.inc'
c---------------
c parameter type
c---------------
      character*8 string
c-----------------
c local variables
c-----------------
      dimension thsave(monoms,mstack),tmhsav(6,6,mstack)
c making(np) = 1  while this lump is under construction
      dimension making(mstack)
c lempty = .true. , when stack is empty
      logical lempty
      save thsave,tmhsav,making,lempty
c-----------------------------------------------------------------------
c start
c-----------------------------------------------------------------------
c  initialize arrays
c
      do 5 k = 1,mstack
        making(k) = 0
  5   continue
      call initst(string,nrept)
c???????????????????????????????????????????????????????????????????????
c       printout counters: prints out first nm calls to lmnt
c
c      nm = 100
c      nc = 0
c      write(jof, 700)
c      write(jodf,700)
c 700  format(/' entering trobj:'/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c ignore incoming line with zero repetition number
c
      if(ntype.eq.2.and.nrept.eq.0) then
        write(jof , 510) ilbl(ith)
        write(jodf, 510) ilbl(ith)
 510    format(1x,'line ',a8,' with rep no. = 0 has been ignored')
        return
      endif
c
c-----------------------------------------------------------------------
c here we start with a stack element, either a new one or one that
c just has been popped
c
 1000 continue
c
      if(ntype .eq. 1) then
c menu element should never occur here ...
        write(jodf,910) lstac(np)
        write(jof ,910) lstac(np)
  910   format(1x,'error in trobj: menu element ',a8,
     # ' found at beginning of routine.')
        call myexit
      else if (ntype.eq.4) then
c ... neither should a loop
        write(jodf,920) lstac(np)
        write(jof ,920) lstac(np)
  920   format(1x,'error in trobj: loop ',a8,
     # ' found at beginning of routine.')
        call myexit
      else if (ntype.eq.5) then
c unknown label is ignored
        write(jodf,930) lstac(np)
        write(jof, 930) lstac(np)
  930   format(1x,'warning from trobj: object ',a8,' not found.')
c
c main part is lumps or lines (should be the only part used)
c
      else if (ntype.eq.3 .and. loop(np).eq.0) then
c ignore lump with rep factor 0
          write(jof, 520) ilbl(ith)
          write(jodf,520) ilbl(ith)
 520      format(1x,'lump ',a8,' with rep no. = 0 has been ignored')
      else if (ntype.eq.3 .and. lmade(ith).ne.0) then
c a lump which is already in the core
c???????????????????????????????????????????????????????????????????????
c        write(jof, 710) ith
c        write(jodf,710) ith
c 710    format(/' picking up old lump, no.',i4)
c        if(nc.lt.nm) write(jodf,715) jth,lmade(jth),loop(np),np
c        if(nc.lt.nm) write(jof, 715) jth,lmade(jth),loop(np),np
c 715    format(/5x,'jth =',i4,'    lmade(jth) =',i4,
c     &          5x,'lumprep =',i3,'  np =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - combine previously made lump with the total map
        call comwtm(lmade(ith),loop(np))
c
      else if (ntype.eq.2 .or. ntype.eq.3) then
c all other lines and lumps
c
c first, lumps need special handling, if they are new
        if (ntype.eq.3 .and. making(np).ne.1) then
          making(np) = 1
c???????????????????????????????????????????????????????????????????????
c          write(jodf,720) ith,np,loop(np)
c          write(jof, 720) ith,np,loop(np)
c 720      format(/' starting lump no. ',i3,
c     1            '  np =',i3,'  loop(np) =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - save current map and initialize a new lump
          call buffin(th,tmh,thsave,tmhsav,np)
          call ident(th,tmh)
        endif
c
c here we begin a computational loop, handling simple elements
c in lines or lumps
c
 2000   continue
c
        if(ntype.eq.2 .and.
     &     (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a line
c????????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,730)
c          if(nc.lt.nm) write(jof, 730)
c  730     format(/' *********************** end of line')
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          loop(np) = loop(np) - nsign(loop(np))
          if(loop(np).ne.0) then
            nslot(np) = newsl(np,ith)
            call lookup(icon(nslot(np),ith),mtype,jth)
            goto 1000
          endif
        else if (ntype.eq.3 .and.
     &           (nslot(np).eq.0 .or. nslot(np).gt.ilen(ith))) then
c end of a lump
c????????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,740)
c          if(nc.lt.nm) write(jof, 740)
c  740     format(/' *********************** end of lump')
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c - store lump
          call lumpit(ith)
c - recall running total map
          call bufout(thsave,tmhsav,th,tmh,np)
c - combine new lump with total map:
c???????????????????????????????????????????????????????????????????????
c          write(jof, 750) loop(np)
c          write(jodf,750) loop(np)
c 750      format(/'  combine new lump with total map: lumprep =',i3/)
c          write(jof, 755) ith,lmade(ith)
c          write(jodf,755) ith,lmade(ith)
c 755      format(/' ith =',i3,'   lmade =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call comwtm(lmade(ith),loop(np))
          making(np) = 0
        else if (mtype.eq.4) then
c lump or line points to a loop; error
          write(jof, 940) lstac(np),icon(nslot(np),ith)
          write(jodf,940) lstac(np),icon(nslot(np),ith)
 940      format(1x,'error in trobj: ',a8,
     &    ' contains a loop (',a8,').')
          call myexit
        else if (mtype.eq.2 .and. irep(nslot(np),ith).eq.0) then
c if rep no. = 0 ignore line
          write(jof ,510) icon(nslot(np),ith)
          write(jodf,510) icon(nslot(np),ith)
          nslot(np) = nslot(np) + nsign(loop(np))
          call lookup(icon(nslot(np),ith),mtype,jth)
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,760) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 760) ith,mtype,jth
c 760      format(/' 0*line: ith, mtype,   jth:'/5x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.5) then
c ignore unknown label
          write(jodf,930) icon(nslot(np),ith)
          write(jof, 930) icon(nslot(np),ith)
          nslot(np) = nslot(np) + nsign(loop(np))
          call lookup(icon(nslot(np),ith),mtype,jth)
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,762) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 762) ith,mtype,jth
c 762      format(/' unkn. : ith, mtype,   jth:'/5x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.1) then
c line or lump points to an element
c???????????????????????????????????????????????????????????????????????
c          nc = nc + 1
c          if(nc.lt.nm) write(jodf,770) np
c          if(nc.lt.nm) write(jof, 770) np
c 770      format(/' element in line or lump: np =',i3/)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          call trlmnt(jth,irep(nslot(np),ith))
c - next element:
          nslot(np) = nslot(np) + nsign(loop(np))
          call lookup(icon(nslot(np),ith),mtype,jth)
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jodf,764) ith,mtype,jth
c          if(nc.lt.nm) write(jof, 764) ith,mtype,jth
c 764      format(/' next icon: ith, mtype,   jth:'/7x,3i7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 2000
        else if (mtype.eq.2 .or. mtype.eq.3) then
c line or lump points to line or lump
          call push
c???????????????????????????????????????????????????????????????????????
c          if(nc.lt.nm) write(jof, 780) np
c          if(nc.lt.nm) write(jodf,780) np
c780      format(/' *********************** pushing stack: new np =',i3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          goto 1000
        endif
c
c end of interpretation of a line/lump
c
      endif
c end of main way through stack-element
c-----------------------------------------------------------------------
c     pop stack; see what's there. this is the normal way to return
c
      call pop(lempty)
      if (lempty) return
c
c???????????????????????????????????????????????????????????????????????
c      write(jof, 790) np
c      write(jodf,790) np
c790  format(/' *************************** popping stack: new np =',i3)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
      goto 1000
      end
c
c end of file
