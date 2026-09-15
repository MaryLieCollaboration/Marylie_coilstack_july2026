************************************************************************
* header:                 INPUT AND OUTPUT                             *
* Input and output for maps, matrices, files, arrays, and parameter    *
* sets                                                                 *
************************************************************************
      subroutine bell
c This subroutine rings the bell at the terminal
c Written by Alex Dragt, 13 July 1988
c
      include 'files.inc'
c
      character*1 ding
      ding=char(7)
c
      write(jof,*) ding
c
      return
      end
c
*********************************************************************
c
      subroutine cf(p)
c  subroutine to close files
c  Written by Alex Dragt, Spring 1987
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),ip(6)
      include 'files.inc'
c
c set up control indices
      do 10 j=1,6
      ip(j)=nint(p(j))
   10 continue
c
c close indicated files
      do 20 j=1,6
      n=ip(j)
      if( n.gt.0) close(unit=n, err=30)
      go to 20
   30 write(jof,*) 'error in closing file unit ',n
   20 continue
c
      return
      end
c
******************************************************************
c
      subroutine flag(p)
c display and set flags in block commons status and files
c
      include 'status.inc'
      include 'files.inc'
c
c calling arrays
      dimension p(6)
c
c set up control parameters
      job=nint(p(1))
c
c change flags if job=0
      if (job .eq. 0) then
      continue
      endif
c
c write flags if job=1, 2, or 3
      if (job .eq. 1 .or. job .eq. 3) then
      write (jof,10) imbad,iquiet
   10 format (/,1h ,'imbad =',i3,3x,'iquiet =',i3)
      endif
      if (job .eq. 2 .or. job .eq. 3) then
      write (jodf,10) imbad,iquiet
      endif
c
      return
      end
c
************************************************************************
c
      subroutine fwa(pp)
c
c   Prototype routine to write results from a file into the array ucalc.
c   At the moment results are put into ucalc.  Eventually they should
c   be put in the working array wa(*).
c   Written by A. Dragt 8/28/92.
c
      include 'impli.inc'
      include 'param.inc'
c
c common blocks
c
      include 'elmnts.inc'
      include 'items.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
      include 'usrdat.inc'
c
      dimension pp(6)
c
c local variables
c
      dimension temp(6)
c
c  set up control indices
c
      ifile=nint(pp(1))
      jcol=nint(pp(2))
      item=nint(pp(3))
 
c
c  start routine
c
      read (ifile,*) temp(1),temp(2),temp(3),temp(4),temp(5),temp(6)
      ucalc(item) = temp(jcol)
      do 10 i=1,6
      wa(i)=temp(i)
  10  continue
c
      return
      end
c
************************************************************************
c
      subroutine mapin(nopt,nskp,h,mh)
c read nonzero matrix elements and monomials from file unit mpi
      include 'impli.inc'
      include 'param.inc'
      double precision h,mh
      include 'files.inc'
      dimension mh(6,6),h(monoms)
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
c
      if(nopt.eq.0)goto 5
      rewind mpi
      write(jof,210)mpi
    5 continue
c     initialize arrays:
      do 10 j=1,monoms
   10 h(j)=0.d0
      do 20 k=1,6
      do 20 l=1,6
   20 mh(k,l)=0.
c
c  skip nskp maps:
      ns=nskp
   55 if(ns.eq.0)goto 100
   60 read(mpi,*)i,j,temp
      if(i.eq.6.and.j.eq.6)goto 80
      goto 60
   80 read(mpi,*)k,temp
      if(k.eq.monoms)goto 90
      goto 80
   90 ns=ns-1
      goto 55
c
c  now read in the map:
  100 continue
  160 read(mpi,*)i,j,temp
      mh(i,j)=temp
      if(i.eq.6.and.j.eq.6)goto 180
      goto 160
  180 read(mpi,*)k,temp
      h(k)=temp
      if(k.eq.monoms)goto 200
      goto 180
  200 continue
      write(6,205) mpi,nskp
  205 format(1x,'map read in from file ',i3,'; ',i3,
     #' record(s) skipped')
  210 format(1x,'file unit ',i2,' rewound')
      return
      end
c
***********************************************************************
c
      subroutine mapout(nopt,h,mh)
c  output present matrix and polynomials (nonzero values)
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1984
c   modified Oct 89 by Tom Mottershead to work without requiring
c   a prior map file.
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
c
c calling arrays
c
      double precision h,mh
      dimension mh(6,6),h(monoms)
c
c check to see if file assignment makes sense
      if(mpo .le. 0) return
c
      nout=mpo
ctm:5/2013 NOTE change: gfortran thinks reading past EOF is a fatal error.
c So skip error assumptions and just write a new file, the normal case
c  the parameter nopt is not currently used
c
c  go to end of file. Actually, skip it.  CTM 5/2013
c
c  15 read(nout,*,end=17,err=19)dummy
c     goto 15
c  17 write(jof,18) nout
c  18 format(' adding current map to file on unit',i3)
c     go to 25
c
c      read error means map file does not exist, so start a new one
c
  19  write(jof,21) nout
  21  format(' starting new map file on unit',i3)
c
c  write out map
c
  25  continue
      do 30 i=1,6
      do 30 j=1,6
   30 if(mh(i,j).ne.0.)write(nout,*)i,j,mh(i,j)
      if(mh(6,6).eq.0.)write(nout,*)6,6,mh(6,6)
      do 40 k=1,monoms
   40 if(h(k).ne.0.)write(nout,*)k,h(k)
      if(h(monoms).eq.0.)write(nout,*)monoms,h(monoms)
c      write(jof,100) nout
c  100 format(1x,'map written on file ',i2)
      return
      end
c
***********************************************************************
c
      subroutine mapsnd(iopt,nmap,ta,tm,ha,hm)
c this is a subroutine for sending a map to some buffer or file
c  Written by Alex Dragt, Spring 1987
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'buffer.inc'
c
      dimension ta(monoms),ha(monoms)
      dimension tm(6,6),hm(6,6)
      character*3 kynd
c
c procedure when iopt = 0
      if (iopt.eq.0) then
      if (nmap.ge.1 .and. nmap.le.5) then
      kynd='stm'
      call strget(kynd,nmap,ta,tm)
      endif
c
      if (nmap.eq.0) call mapmap(ta,tm,ha,hm)
c
      if (nmap.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmap.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmap.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmap.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmap.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
      endif
c
c procedure when iopt .gt. 0
      if (iopt.gt.0) then
      mpot=mpo
      mpo=iopt
      call mapout(0,ta,tm)
      mpo=mpot
      endif
c
      return
      end
c
**************************************************************
c
      subroutine num(p)
c This is a subroutine for numbering lines in a file
c Written by Alex Dragt, Spring 1987
      include 'impli.inc'
      dimension p(6)
      include 'rays.inc'
c
c set up control indices
      iopt=nint(p(1))
      nfile=nint(p(2))
      ifirst=nint(p(3))
      istep=nint(p(4))
c
c procedure for writing out a file with numbered lines
      line=ifirst
      do 10 i=1,nrays
      write(nfile,100) line,
     # zblock(i,1), zblock(i,2), zblock(i,3),
     # zblock(i,4), zblock(i,5), zblock(i,6)
  100 format(1x,i5,6(1x,1pe11.4))
      line=line+istep
   10 continue
c
      return
      end
c
***********************************************************************
c
      subroutine of(p)
c  subroutine to open files
c  Written by Alex Dragt, Spring 1987
      include 'impli.inc'
      character*6 unit
      dimension p(6),ip(6)
      dimension unit(50)
      include 'files.inc'
c
c set up file names
c
      data (unit(i),i=1,50)/
     #'unit01','unit02','unit03','unit04','unit05',
     #'unit06','unit07','unit08','unit09','unit10',
     #'unit11','unit12','unit13','unit14','unit15',
     #'unit16','unit17','unit18','unit19','unit20',
     #'unit21','unit22','unit23','unit24','unit25',
     #'unit26','unit27','unit28','unit29','unit30',
     #'unit31','unit32','unit33','unit34','unit35',
     #'unit36','unit37','unit38','unit39','unit40',
     #'unit41','unit42','unit43','unit44','unit45',
     #'unit46','unit47','unit48','unit49','unit50'/
c
c set up control indices
      do 10 j=1,6
      ip(j)=nint(p(j))
   10 continue
c
c open indicated files
      do 20 j=1,6
      n=ip(j)
      if( n.gt.0 .and. n.le.50 .and. n.ne.lf .and. n.ne.jof
     # .and. n.ne.jodf) then
      open(unit=n, file=unit(n), status='unknown', err=30)
      endif
      go to 20
   30 write(jof,*) 'error in opening file unit ',n
   20 continue
c
      return
      end
c
***********************************************************************
c
      subroutine pcmap(n1,n2,n3,n4,fa,fm)
c  routine to print m,f3,f4 and t,u.
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
c
      include 'impli.inc'
      include 'param.inc'
      integer colme(6,0:6)
      include 'expon.inc'
      include 'pbkh.inc'
      include 'files.inc'
      include 'zeroes.inc'
c
      dimension fa(monoms),fm(6,6)
      dimension t(monoms),u(monoms),u2(monoms)
c
c  test for matrix write
      if(n1.eq.0) goto 20
c
c  procedure for writing out matrix
c  write matrix at terminal
      if(n1.eq.1.or.n1.eq.3)then
        write(jof,13)
   13   format(/1h ,'matrix for map is :'/)
        write(jof,15)((fm(k,i),i=1,6),k=1,6)
   15   format(6(1x,1pe12.5))
      endif
c  write matrix on file 12
      if(n1.eq.2.or.n1.eq.3)then
        write(jodf,13)
        write(jodf,15)((fm(k,i),i=1,6),k=1,6)
c        write(88,115)((fm(k,i),i=1,6),k=1,6)
 115  format(2x,2(1pg23.15))
      endif
c
c  test for polynomial write
   20 continue
      if(n2.eq.0)goto 30
c
c  procedure for writing out polynomial
c  write polynomial at terminal
      if(n2.eq.1.or.n2.eq.3)then
        write(jof,22)
   22   format(/1h ,'nonzero elements in generating polynomial are :'/)
        do 25 i=1,monoms
c  test to see if fa(i) is large enough to write
        if( abs(fa(i)) .le. fzer) goto 25
        write(jof,27)i,(expon(j,i),j=1,6),fa(i)
   27   format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',1pg21.14)
   25   continue
      endif
c  write polynomial on file 12
      if(n2.eq.2.or.n2.eq.3)then
        write(jodf,22)
        do 26 i=1,monoms
c  test to see if fa(i) is large enough to write
        if( abs(fa(i)) .le. fzer ) goto 26
        write(jodf,27)i,(expon(j,i),j=1,6),fa(i)
   26   continue
      endif
c
c  prepare for higher order matrix write if required
   30 continue
      if(n3.gt.0.or.n4.gt.0) call brkts(fa)
c
c  test for t matrix write
      if(n3.eq.0) goto 40
c
c  procedure for writing t matrix
c  write out heading
      if(n3.eq.1.or.n3.eq.3)write(jof,32)
      if(n3.eq.2.or.n3.eq.3)write(jodf,32)
   32   format(/1h ,'nonzero elements in second order matrix are :'/)
c  write out contents
        do 35 i=1,6
        call xform(pbh(1,i),2,fm,i-1,t)
        do 36 n=7,27
c  test to see if t(n) is large enough to write
        if( abs(t(n)) .le. fzer ) goto 36
        if(n3.eq.1.or.n3.eq.3)
     #  write(jof,38) i,n,i,(expon(j,n),j=1,6),t(n)
        if(n3.eq.2.or.n3.eq.3)
     #  write(jodf,38) i,n,i,(expon(j,n),j=1,6),t(n)
   38 format(2x,'t',i1,'(',i3,')','=t',i1,'( ',3(2i1,1x),')=',1pg21.14)
c        write(88,138) i,n,i,(expon(j,n),j=1,6),t(n)
c  138  format(2i5,7i3,1pg25.15)
   36   continue
   35   continue
c
 
c  test for u matrix write
   40 continue
      if(n4.eq.0) goto 50
c
c  procedure for writing u matrix
c  write out heading
      if(n4.eq.1.or.n4.eq.3)write(jof,42)
      if(n4.eq.2.or.n4.eq.3)write(jodf,42)
   42 format(/1h ,'nonzero elements in third order matrix are :'/)
c  write out contents
         do  44 i=1,6
         call xform(pbh(1,i),3,fm,i-1,u)
         call xform(pbh(1,i+6),3,fm,1,u2)
         do  45 n=28,83
         u(n)=u(n)+u2(n)/2.d0
c  test to see if u(n) is large enough to write
        if( abs(u(n)) .le. fzer ) goto 45
         if(n4.eq.1.or.n4.eq.3)
     #   write(jof,46) i,n,i,(expon(j,n),j=1,6),u(n)
         if(n4.eq.2.or.n4.eq.3)
     #   write(jodf,46) i,n,i,(expon(j,n),j=1,6),u(n)
   46 format(2x,'u',i1,'(',i3,')','=u',i1,'( ',3(2i1,1x),')=',1pg21.14)
c        write(88,138) i,n,i,(expon(j,n),j=1,6),u(n)
   45    continue
   44    continue
c
c  procedure if all nj's are zero or are faulty
   50 continue
c
      return
      end
c
***********************************************************************
c
      subroutine pdrmap(n1,n2,fa,fm)
c  routine to print mij in the cartesian basis
c  and the f's in the dynamic resonance basis.
c    F. Neri  6/3/1986
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'drl.inc'
      include 'zeroes.inc'
c
      dimension fa(monoms),fm(6,6)
c
c  beginning of routine
c
c  writing out matrix
      if(n1.eq.0)goto 20
      if(n1.eq.1.or.n1.eq.3)write(jof,33)
      if(n1.eq.2.or.n1.eq.3)write(jodf,33)
   33 format(/1h ,'matrix for map is :'/)
      if(n1.eq.1.or.n1.eq.3)write(jof,35)((fm(k,i),i=1,6),k=1,6)
      if(n1.eq.2.or.n1.eq.3)write(jodf,35)((fm(k,i),i=1,6),k=1,6)
c   35 format(6(1x,e12.5))
   35 format(6(1x,1pe12.5))
c
c  writing out f's
   20 if(n2.eq.0)goto 30
      if(n2.eq.1.or.n2.eq.3)write(jof,55)
      if(n2.eq.2.or.n2.eq.3)write(jodf,55)
   55 format(/1h ,'nonzero elements in generating polynomial in',/,
     #1x,'the dynamic resonance basis are :'/)
      do 22 i=1,monoms
c  test to see if fa(i) is large enough to write
        if( abs(fa(i)) .le. fzer ) goto 22
      if(n2.eq.1.or.n2.eq.3)write(jof,60)i,dln(i),fa(i)
      if(n2.eq.2.or.n2.eq.3)write(jodf,60)i,dln(i),fa(i)
   60 format(1x,'f(',i3,')=f( ',a7,' )=',1pg21.14)
   22 continue
c
   30 continue
      return
      end
c
*****************************************************************
c
      subroutine pmif(iu,itype)
c  subroutine to write out the master input file
c  written by Rob Ryne ca 1984
c
      include 'impli.inc'
      include 'param.inc'
      include 'bmline.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'coment.inc'
      include 'codes.inc'
      include 'mldex.inc'
      include 'files.inc'
      include 'parm.inc'
      include 'cxdata.inc'
c
      character*79 line
c
       save
c------------------
c start routine
c------------------
      goto(4,90)itype
c  echo exact contents of file lf:
      rewind lf
    1 read(lf,1002,end=3)line
 1002 format(a)
      write(iu,1003)line
 1003 format(1x,a)
      goto 1
    3 continue
      write(iu,*)
      return
c
c  comments
    4 if(np.eq.0)goto 5
      write(iu,530)ling(1)
      write(iu,6)(mline(i),i=1,np)
c mline is character*80, but is written out with an a79 format
c so that the last character (which should be a blank) does not
c cause carriage returns (and hence blank lines) on laser printers
    6 format(a79)
c    6 format(a80)
c  beam
    5 write(iu,530)ling(2)
      write(iu,*)brho
      write(iu,*)gamm1
      write(iu,*)achg
      write(iu,*)sl
c  menu
      write(iu,530)ling(3)
      do 10 k=1,na
         write(iu,600)lmnlbl(k),ltc(nt1(k),nt2(k))
         imax=nrp(nt1(k),nt2(k))
         if(imax.eq.0)goto 10
         write(iu,603)(pmenu(i+mpp(k)),i=1,imax)
   10 continue
c
ctm2014  new internal text files
c
      if(nxd.eq.0) go to 28
      write(iu,530)ling(8)
      do 20 jj = 1,nxd
         lux = lunex(jj)
         write(6,13) jj, lux, mltext(jj)
  13     format('jj,lux=',2i4,a)
         if(lux.ne.kurid) then
            kurid = lux
            write(iu,17) kurid
  17   format(6x,'ixd>:',i4)
         endif
         write(iu,*) mltext(jj)
  20  continue
c
c  lines,lumps,loops
c
  28  if(nb.eq.0)goto 90
      do 40 ii=2,4
      write(iu,530)ling(ii+2)
      do 30 k=1,nb
      if(ityp(k).ne.ii)goto 30
      write(iu,790)ilbl(k)
      write(iu,791)(irep(l,k),icon(l,k),l=1,ilen(k))
   30 continue
   40 continue
c  labor
   90 if(noble.eq.0)goto 999
      write(iu,530)ling(7)
      do 100 j=1,noble
  100 write(iu,800)num(j),latt(j)
c  529 format(1h ,'#comments')
  530 format(1h ,a8)
  600 format(1h ,1x,a8,1x,a8)
c  603 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  603 format((1h ,3(1x,1pg22.15)))
  790 format(1h ,1x,a8)
  791 format((1h ,1x,5(i5,'*',a8),1x,:'&'))
  800 format(1h ,1x,i4,'*',a8)
  999 continue
c      write(iu,*)
      return
      end
c
***********************************************************************
c
      subroutine pset(p,k)
c  subroutine to read in parameter set values
c  the integer k labels the parameter set
c  Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
      dimension p(6)
c
c test on k
      if (k.gt.maxpst .or. k.le.0) then
      write(6,*) 'improper attempt to store parameters in',
     # ' a parameter set with k=',k
      call myexit
      endif
c
c store parameter values
      do 10 i=1,6
   10 pst(i,k)=p(i)
c
      return
      end
c
***********************************************************************
c
      subroutine psrmap(n1,n2,fa,fm)
c  routine to print mij in the cartesian basis
c  and the f's in the static resonance basis.
c  Written by Alex Dragt, Fall 1986
c
      include 'impli.inc'
      include 'param.inc'
      include 'srl.inc'
      include 'files.inc'
      include 'zeroes.inc'
c
      dimension fa(monoms),fm(6,6)
c
c  beginning of routine
c
c  writing out matrix
      if(n1.eq.0)goto 20
      if(n1.eq.1.or.n1.eq.3)write(jof,33)
      if(n1.eq.2.or.n1.eq.3)write(jodf,33)
   33 format(/1h ,'matrix for map is :'/)
      if(n1.eq.1.or.n1.eq.3)write(jof,35)((fm(k,i),i=1,6),k=1,6)
      if(n1.eq.2.or.n1.eq.3)write(jodf,35)((fm(k,i),i=1,6),k=1,6)
c   35 format(6(1x,e12.5))
   35 format(6(1x,1pe12.5))
c
c writing out f's
   20 if(n2.eq.0)goto 30
      if(n2.eq.1.or.n2.eq.3)write(jof,55)
      if(n2.eq.2.or.n2.eq.3)write(jodf,55)
   55 format(/1h ,'nonzero elements in generating polynomial in',/,
     #1x,'the static resonance basis are :'/)
      do 22 i=1,monoms
c  test to see if fa(i) is large enough to write
      if( abs(fa(i)) .le. fzer ) goto 22
      if(n2.eq.1.or.n2.eq.3)write(jof,60)i,sln(i),fa(i)
      if(n2.eq.2.or.n2.eq.3)write(jodf,60)i,sln(i),fa(i)
   60 format(1x,'f(',i3,')=f( ',a6,' )=',1pg21.14)
   22 continue
c
   30 continue
      return
      end
c
******************************************************************
c
      subroutine randin(nfile,kt1,kt2,p)
c  This is a subroutine for reading in and checking the parameters
c  for random elements.
c  Written by Alex Dragt, ca 1985, and modified by F. Neri
c  and Alex Dragt on 29 January 1988
      include 'impli.inc'
      include 'param.inc'
      double precision p(6)
      include 'files.inc'
      include 'codes.inc'
      include 'parset.inc'
c
      if(nfile.gt.0) then
c  Read in parameters from file:
   2    read(nfile,*,end=4,err=6)(p(i),i=1,nrp(kt1,kt2))
        goto 8
   4    rewind nfile
        goto 2
   6    write(6,7) nfile,kt1,kt2
   7    format(1x,'nfile=',i4,' kt1,kt2=',2i4,' read error in randin')
        call myexit
   8    continue
      if(kt1.eq.3) return
      endif
      if(nfile.lt.0) then
c  Read in parameters from parameter sets:
        ipset = -nfile
        if(ipset.gt.maxpst) then
          write(jof,*) ' num of pset too large in rnd lmnt.'
          call myexit
        endif
        do 7008 i=1,6
          p(i) = pst(i,ipset)
 7008   continue
      endif
      if(nfile.eq.0) then
        write(jof,*) ' zero isource number in rnd lmnt.'
        call myexit
      endif
c
c  Check on parameter values corresponding to control indices:
c
c         'rdrft   ','rnbnd   ','rpbnd   ','rgbnd   ','rprot   ',
      goto(10,        20,        30,        40,        50,
c         'rgbdy   ','rfrng   ','rcfbd   ','rquad   ','rsext   ',
     &     60,        70,        80,        90,        100,
c         'roctm   ','rocte   ','rsrfc   ','rarot   ','rtwsm   ',
     &     110,       120,       130,       140,       150,
c         'rthlm   ','rcplm   ','rcfqd   ','riftm   ','rsol    ',
     &     160,       170,       180,       190,       200,
c         'dummark ','dumjmap ','dumdp   ','rrect   '/
     &     210,       220,       230,       240),kt2
c
c complain if get here
      write(jof,*) 'kt2=',kt2,
     & ' error in randin, kt2 out of range'
      call myexit
      return
c
c  drift:
  10  continue
      return
c
c  normal entry bend:
  20  call rcheck(0,1,p(3),'lfrn','nbnd',nfile)
      call rcheck(0,1,p(4),'tfrn','nbnd',nfile)
      return
c
c  parallel faced bend:
  30  continue
      return
c
c  general bending magnet:
  40  call rcheck(0,1,p(5),'lfrn','gbnd',nfile)
      call rcheck(0,1,p(6),'tfrn','gbnd',nfile)
      return
c
c  leading or trailing rotation for a parallel faced bend:
  50  call rcheck(0,1,p(2),'kind','prot',nfile)
      return
c
c  body of a general bending magnet:
  60  continue
      return
c
c  hard edge fringe fields of a normal entry bend:
  70  call rcheck(0,1,p(2),'iedg','frng',nfile)
      return
c
c  combined function bend:
  80  continue
      return
c
c  quadrupole:
  90  call rcheck(0,1,p(3),'lfrn','quad',nfile)
      call rcheck(0,1,p(4),'tfrn','quad',nfile)
      return
c
c  sextupole:
 100  continue
      return
c
c  mag. octupole:
 110  continue
      return
c
c  elec. octupole:
 120  continue
      return
c
c  short rf cavity
 130  continue
      return
c
c  axial rotation:
 140  continue
      return
c
c  linear matrix via twiss parameters:
 150  call rcheck(1,3,p(1),'ipla','twsm',nfile)
      return
c
c  thin lens low order multipole:
 160  continue
      return
c
c  "compressed" low order multipole:
 170  continue
      return
c
c  combined function quadrupole
 180  continue
      return
c
c  initial final transfer matrix
 190  continue
      return
c
c  solenoid
 200  continue
      return
c
c  dummy marker
 210  continue
      write(jof,*) ' error in randin, dummy marker type code reached'
      call myexit
      return
c
c  dummy jmap
 220  continue
      write(jof,*) ' error in randin, dummy jmap type code reached'
      call myexit
      return
c
c  dummy data point
 230  continue
      write(jof,*)
     & ' error in randin, dummy data point type code reached'
      call myexit
      return
c
c  REC quadrupole triplet
 240  continue
      return
c
      end
c
***************************************************************
c
      subroutine raysin(icfile)
c  routine to read in initial conditions of rays to be traced
c  and to initialize the arrays istat and ihist
c  Written by Robert Ryne ca 1984, and modified by Alex Dragt
      include 'impli.inc'
      include 'rays.inc'
      include 'files.inc'
      include 'parset.inc'
c
c procedure for reading a single set of initial conditions from
c a parameter set
c
      if(icfile.lt.0) then
      ipset=-icfile
      if(ipset.gt.maxpst) then
      write(jof,*) 'parameter icfile out of range'
      call myexit
      return
      endif
      nrays=1
      do 10 j=1,6
   10 zblock(1,j)=pst(j,ipset)
      goto 130
      endif
c
c procedure for reading initial conditions from a file
c
      nrays=0
      rewind icfile
      do 100 i=1,maxray
      read(icfile,*,end=110,err=120)(zblock(i,j),j=1,6)
  100 nrays=i
  110 continue
      if (nrays.eq.0) goto 120
      goto 125
  120 write(jof,121)nrays
  121 format(1x,'trouble in routine raysin following record #',i5)
      call myexit
      return
  125 continue
      write(jof,117)nrays,icfile
  117 format(1x,i5,' ray(s) read in from file ',i3)
c
  130 continue
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
      return
      end
c
***********************************************************************
c
      subroutine rcheck(min,max,arg,ivar,itype,nfile)
c  This is a subroutine for checking that control
c  parameters for random elements lie within the allowed range.
c  Written by Alex Dragt ca 1985
      double precision arg
      character*4 ivar,itype
      iarg=nint(arg)
      if (iarg.lt.min.or.iarg.gt.max) goto 10
      return
   10 write (6,100) ivar,iarg,itype,nfile
  100 format(1x,a4,'=',i4,1x,'trouble with random element',1x,a4,1x,
     #'read from file',1x,i4)
      call myexit
      return
      end
c
c**********************************************************************
c
      subroutine subctr(p)
c  subroutine to change tune range
c  Written by Alex Dragt, 18 August 1988
c
      include 'impli.inc'
c
      dimension p(6)
c
      write(6,*) 'in subroutine subctr'
c
      return
      end
c
************************************************************************
c
      subroutine wcl(pp)
c  subroutine to write out contents of a loop
c  Written by Alex Dragt, 23 August 1988
c  Based on the subroutines cqlate and pmif
c
      include 'impli.inc'
      include 'param.inc'
c
c common blocks
c
      include 'elmnts.inc'
      include 'items.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
c
      dimension pp(6)
c
c local variables
c
      character*8 string(5),str
      logical     ljof,ljodf
c
c  set up control indices
c
      iopt=nint(pp(1))
      ifile=nint(pp(2))
      isend=nint(pp(3))
c
c  start routine
c
c  see if a loop exists
      if(nloop.le.0) then
      write(jof ,*) ' error from wcl: no loop has been specified'
      write(jodf,*) ' error from wcl: no loop has been specified'
      return
      endif
c
c  procedure when iopt=1 (write only names of loop contents)
      if (iopt.eq.1) then
      ljof  = isend.eq.1 .or. isend.eq.3
      ljodf = isend.eq.2 .or. isend.eq.3
      if (ljof .or. ljodf) then
c  write loop name
      if(ljof ) write(jof ,510) ilbl(nloop)
      if(ljodf) write(ifile,510) ilbl(nloop)
  510 format(/,1h ,'contents of loop ',a8,' :')
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
      if(ljodf) write(ifile,511)(string(k),k=1,kmax)
  511 format(' ',5(1x,a8))
   1  continue
      endif
      endif
c
c  procedure when iopt=2 (write out names of loop contents with & signs)
c  code to be installed
      if(iopt .eq. 2) write(6,*) 'option 2 not yet installed for wcl'
c
c  procedure when iopt=3 (write out names of loop contents with % and & signs)
c  code to be installed
      if(iopt .eq. 3) write(6,*) 'option 3 not yet installed for wcl'
c
c  procedure when iopt=4 or iopt=5 (write out full loop contents)
      if (iopt.eq.4 .or. iopt.eq.5) then
      ljof  = isend.eq.1 .or. isend.eq.3
      ljodf = isend.eq.2 .or. isend.eq.3
      if (ljof .or. ljodf) then
c
      do 10 i=1,2
      iskip=0
      if(i.eq.1 .and. ljof) then
      iu=jof
      iskip=1
      endif
      if(i.eq.2 .and. ljodf) then
      iu=ifile
      iskip=1
      endif
      if (iskip.eq.0) go to 10
c  comments
      write(iu,530) ling(1)
c  write loop name
      write(iu,512) ilbl(nloop)
  512 format(1h ,' contents of loop ',a8)
c  beam
      write(iu,530) ling(2)
      write(iu,*) brho
      write(iu,*) gamm1
      write(iu,*) achg
      write(iu,*) sl
c  biglist heading
      write(iu,127)
  127 format(1h ,'#biglist')
   10 continue
c
c  contents of biglist
c  write loop contents
      do 137 jk1=1,joy
c element
        if(mim(jk1).lt.0) then
          string(1)=lmnlbl(-mim(jk1))
c user supplied element
        else if(mim(jk1).gt.5000) then
          string(1)=lmnlbl(mim(jk1)-5000)
c lump
        else
          string(1)=ilbl(inuse(mim(jk1)))
        endif
      call lookup(string(1),itype,item)
c      write(6,513) string(1)
c  513 format(1x,a8)
c      write(6,*) 'itype and item are ',itype, item
c procedure for a menu item
      if(itype.eq.1) then
      k=item
c case where iopt=4
      if (iopt.eq.4) then
      if(ljof) write(jof,600)lmnlbl(k),ltc(nt1(k),nt2(k))
      if(ljodf) write(ifile,600)lmnlbl(k),ltc(nt1(k),nt2(k))
  600 format(1h ,1x,a8,1x,a8)
         imax=nrp(nt1(k),nt2(k))
         if(imax.eq.0)goto 137
      if(ljof) write(jof,603)(pmenu(i+mpp(k)),i=1,imax)
      if(ljodf) write(ifile,603)(pmenu(i+mpp(k)),i=1,imax)
c  603 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  603 format((1h ,3(1x,1pg22.15)))
      endif
c case where iopt=5
      if (iopt.eq.5) then
      imax=nrp(nt1(k),nt2(k))
      if(ljof)
     # write(jof,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
      if(ljodf)
     # write(ifile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
  605 format(1h ,1x,a8,1x,a8,1x,i5,1x,i5,1x,i5)
         if(imax.eq.0)goto 137
      if(ljof) write(jof,607)(pmenu(i+mpp(k)),i=1,imax)
      if(ljodf) write(ifile,607)(pmenu(i+mpp(k)),i=1,imax)
c  607 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  607 format((1h ,3(1x,1pg22.15)))
      endif
      endif
c procedure for a lump
      if(itype.eq.3) then
c case where iopt=4
      if (iopt .eq. 4) then
      if(ljof) write(jof,514) string(1)
      if(ljodf) write(ifile,514) string(1)
  514 format(1x,1x,a8,1x,'lump')
      endif
c case where iopt=5
      if (iopt .eq. 5) then
      if(ljof) write(jof,515) string(1),mim(jk1)
      if(ljodf) write(ifile,515) string(1),mim(jk1)
  515 format(1x,1x,a8,1x,'lump',9x,'0',4x,'-1',2x,i4)
      endif
      endif
  137 continue
      endif
      endif
c
  530 format(1h ,a8)
c
      return
      end
c
************************************************************************
c
      subroutine whst(p)
c  subroutine to write out history of beam loss
c  Written by Alex Dragt, Fall 1986
      include 'impli.inc'
      dimension p(6)
      include 'rays.inc'
c  begin routine
      ifile=nint(p(1))
      job=nint(p(2))
c  determine what job is to be done
      if (job.eq.2) goto 200
c  procedure for writing out istat
c
      do 5 k=1,nrays
      write (ifile,50) k, istat(k)
   50 format (1h ,i10,i10)
    5 continue
      if(nrays.gt.0) write(6,*) 'istat written on file ',ifile
      return
c
c  procedure for writing out ihist
c
  200 continue
      do 10 k=1,nlost
      write (ifile,100) k, ihist(k,1), ihist(k,2)
  100 format (1h ,i10,i10,i10)
   10 continue
      write(6,*) 'nlost =',nlost
      if(nlost.eq.0) write(6,*) 'ihist not written out'
      if(nlost.gt.0) write(6,*) 'ihist written on file ',ifile
      return
      end
c
***********************************************************************
c
      subroutine wmrt(ifn,isend)
c subroutine for writing out value of a merit function
c Written by Alex Dragt, 25 July 1988
c
      include 'impli.inc'
      include 'files.inc'
      include 'parset.inc'
      include 'merit.inc'
c
c     write out value of merit function ifn
c
c     if(ifn.lt.0 .or. ifn.gt.5) return
      if(isend.eq.1 .or. isend.eq.3) then
      write(jof,*) 'value of merit function',ifn,' is',val(ifn)
      endif
      if(isend.eq.2 .or. isend.eq.3) then
      write(jodf,*) 'value of merit function',ifn,' is',val(ifn)
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine wps(ipset,isend)
c subroutine for writing out values in a parameter set
c Written by Alex Dragt, 30 January 1988
c
      include 'impli.inc'
      include 'files.inc'
      include 'parset.inc'
c
c Check to see that ipset is within range
      if ((ipset.lt.1) .or. (ipset.gt.maxpst)) then
        write (jof,*) 'WARNING: ipset out of range in command',
     #  ' with typecode wps'
        return
      endif
c
c Write out values of parameters
c
      if ((isend.eq.1) .or. (isend.eq.3)) then
      write (jof,*) 'values of parameters in the parameter set',ipset
      write (jof,100)( pst(j,ipset), j=1,6)
      endif
      if ((isend.eq.2) .or. (isend.eq.3)) then
      write (jodf,*) 'values of parameters in the parameter set',ipset
      write (jodf,100)( pst(j,ipset), j=1,6)
      endif
  100 format((1h ,3(1x,1pg22.15)))
c
      return
      end
c
c end of file
c
*******************************************
c
      subroutine wuca(p)
c subroutine for writing out contents of the array UCALC.
c Written by Alex Dragt, 4 March 1999.
c
      include 'impli.inc'
      include 'usrdat.inc'
      include 'mldex.inc'
      include 'files.inc'
c
      dimension p(6)
c
c set up control indices
c
      kmin=nint(p(1))
      kmax=nint(p(2))
      isend=nint(p(3))
      ifile=nint(p(4))
c
c           ucalc dump
c
c write headings
c
      if((isend .eq. 1) .or. (isend .eq. 3)) then
      write(jof,*) ' k and ucalc(k) for nonzero values of array'
      endif
c
      if((isend .eq. 2) .or. (isend .eq. 3)) then
      write(jodf,*) ' k and ucalc(k) for nonzero values of array'
      endif
c
c write out nonzero values
c
      do 100  k=kmin, kmax
c
      if((ifile .gt. 0) .and. (ucalc(k) .ne. 0.d0)) then
      write(ifile,*) k, ucalc(k)
      endif
c
      if((isend .eq. 1) .or. (isend .eq. 3)) then
      if (ucalc(k) .ne. 0.d0) then
      write(jof,*) k, ucalc(k)
      endif
      endif
c
      if((isend .eq. 2) .or. (isend .eq. 3)) then
      if (ucalc(k) .ne. 0.d0) then
      write(jodf,*) k, ucalc(k)
      endif
      endif
c
 100  continue
c
      return
      end
