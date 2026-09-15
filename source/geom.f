      subroutine geom(pp)
c
c   Routine to compute geometry of a loop.
c   Written by A. Dragt 8/27/92.
c   Modified 5/27/98 AJD.
c   Based on the subroutines cqlate and pmif
c    Outputs modified by CTM  8 June 98
c-------------------------------------------
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
      include 'pie.inc'
      include 'parset.inc'
c
      dimension pp(6)
c
c local variables
c
      character*8 string(5),str
      dimension ex(3), ey(3), ez(3)
      dimension exl(3), eyl(3), ezl(3)
      dimension exlt(3), eylt(3), ezlt(3)
      logical ldb
      ldb = .false.
c
c  set up control indices
c
      iopt=nint(pp(1))
      si=pp(2)
      ti=pp(3)
      ipset1=nint(pp(4))
      ipset2=nint(pp(5))
      isend=nint(pp(6))
c
      jflag=1
c
c get contents of psets
c
      if (ipset1.gt.0 .and. ipset1.le.maxpst) then
      xi=pst(1,ipset1)
      yi=pst(2,ipset1)
      zi=pst(3,ipset1)
      phid=pst(4,ipset1)
      thetad=pst(5,ipset1)
      psid=pst(6,ipset1)
      phir=phid*pi180
      thetar=thetad*pi180
      psir=psid*pi180
      endif
      if (ipset2.gt.0 .and. ipset2.le.maxpst) then
         ifile=nint(pst(1,ipset2))
         jfile=nint(pst(2,ipset2))
         kfile=nint(pst(3,ipset2))
         lfile=nint(pst(4,ipset2))
         if(lfile.gt.0) ldb = .true.
         mfile=nint(pst(5,ipset2))
         ndpts=nint(pst(6,ipset2))
      endif
c
c  set up variables and constants
c
      ss=si
      tt=ti
      x=xi
      y=yi
      z=zi
      do i=1,3
      ex(i)=0.d0
      ey(i)=0.d0
      ez(i)=0.d0
      end do
      ex(1)=1.d0
      ey(2)=1.d0
      ez(3)=1.d0
      if(ldb) write(lfile,*) ' thetar =', thetar
      call erotv(phir,thetar,psir,ex,exl)
      call erotv(phir,thetar,psir,ey,eyl)
      call erotv(phir,thetar,psir,ez,ezl)
      if(ldb) write(lfile,*) 'ezl =',ezl(1),ezl(2),ezl(3)
      vr=1.d0/(beta*c)
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
c scan the loop
c
ctm   write(6,*) '## geom big loop:',joy,'<=joy'
      do 137 jk1=1,joy
c record element category
        icat=0
        wide = 0.0
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
c      if(ldb) write(lfile,513) string(1)
c  513 format(1x,a8)
c      if(ldb) write(lfile,*) 'itype and item are ',itype, item
c procedure for a menu item
      if(itype.eq.1) then
      k=item
      imax=nrp(nt1(k),nt2(k))
c see if item is a simple command
      if (nt1(k) .eq. 7) then
c sapt
      if (nt2(k) .eq. 35) then
      xmin=pmenu(1+mpp(k))
      xmax=pmenu(2+mpp(k))
      ymin=pmenu(3+mpp(k))
      ymax=pmenu(4+mpp(k))
      taumin=pmenu(5+mpp(k))
      taumax=pmenu(6+mpp(k))
      endif
      endif
c see if item is a simple element
ctm   write(6,*) ' Simple element?',k,'=k',nt1(k),nt2(k),'=nt1,nt2'
      if (nt1(k) .eq. 1) then
c drift
      if (nt2(k) .eq. 1) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c spce
      if (nt2(k) .eq. 25) then
      icat=1
      aleng = pmenu(1+mpp(k))
      endif
c quad
      if (nt2(k) .eq. 9) then
      icat=1
      wide = 1.
      aleng = pmenu(1+mpp(k))
      endif
c cfqd
      if (nt2(k) .eq. 18) then
      icat=1
       wide = 1.
      aleng = pmenu(1+mpp(k))
      endif
c sext
      if (nt2(k) .eq. 10) then
      icat=1
        wide = 2.
      aleng = pmenu(1+mpp(k))
      endif
c octm
      if (nt2(k) .eq. 11) then
      icat=1
       wide = 3.
      aleng = pmenu(1+mpp(k))
      endif
c recm
      if (nt2(k) .eq. 24) then
      icat=1
         wide = 4.
      aleng = pmenu(2+mpp(k)) - pmenu(1+mpp(k))
      endif
c sol
      if (nt2(k) .eq. 20) then
      icat=1
       wide = -1.
      aleng = pmenu(2+mpp(k)) - pmenu(1+mpp(k))
      endif
c nbend
      if (nt2(k) .eq. 2) then
      icat=2
        wide = 5.
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c pbend
      if (nt2(k) .eq. 3) then
      icat=2
            wide = 5.
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c gbend
      if (nt2(k) .eq. 4) then
      icat=2
         wide = 5.
      angd=pmenu(1+mpp(k))
      b=pmenu(6+mpp(k))
      endif
c gbdy
      if (nt2(k) .eq. 6) then
      icat=2
         wide = 5.
      angd=pmenu(1+mpp(k))
      b=pmenu(4+mpp(k))
      endif
c cfbd
      if (nt2(k) .eq. 8) then
      icat=2
          wide = 5.
      angd=pmenu(1+mpp(k))
      b=pmenu(2+mpp(k))
      endif
c arot
      if (nt2(k) .eq. 14) then
      icat=3
          wide = 5.
      angd=pmenu(1+mpp(k))
      endif
c prot
      if (nt2(k) .eq. 5) then
c      angd=pmenu(1+mpp(k))
c      kind=nint(pmenu(2+mpp(k)))
      endif
c
c carry out computations
c
c compute and write out data points
      if(iopt .eq. 1) then
c procedure for a straight element
      if(icat .eq. 1) then
      if(ldb) write(lfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,
     * nt1(k),nt2(k)
  605 format(1h ,1x,a8,1x,a8,1x,i5,1x,i5,1x,i5)
         if(imax.eq.0)goto 137
      if(ldb) write(lfile,607)(pmenu(i+mpp(k)),i=1,imax)
c  607 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
  607 format((1h ,3(1x,1pg22.15)))
      if(ldb) write(lfile,*) ' length = ', aleng
      if(ldb) write(lfile,*) ' distance = ', ss
      xx = ss
      yy = xx
      if(ldb) write(lfile,611) ss, aleng, wide, xx, yy, lmnlbl(k)
  611 format(5f12.4,3x,a)
      ndptst=ndpts
      if(ndpts .le. 1) ndptst=2
      daleng=aleng/float(ndptst)
      k1 = nt1(k)
      k2 = nt2(k)
      do i=0,ndptst
      sst = ss + daleng*float(i)
      ttt = tt + daleng*float(i)*vr
      xt = x + daleng*float(i)*ezl(1)
      yt = y + daleng*float(i)*ezl(2)
      zt = z + daleng*float(i)*ezl(3)
      if(ldb) write(lfile,*) ' when i = ',i
      if(ldb) write(lfile,*) ' values of sst and ttt '
      if(ldb) write(lfile,*) sst,ttt
      if(ldb) write(lfile,*) ' values of xt,yt,zt '
      if(ldb) write(lfile,*) xt,yt,zt
      write(35,235) lmnlbl(k),ltc(k1,k2),jk1,icat,k1,k2,xt,yt,zt,wide
  235 format(2x,a,1x,a,i4,3i3,3f13.5,f6.1)
      if(ldb) write(lfile,*) ' vectors exl,eyl,ezl '
      if(ldb) write(lfile,*) exl(1), exl(2), exl(3)
      if(ldb) write(lfile,*) eyl(1), eyl(2), eyl(3)
      if(ldb) write(lfile,*) ezl(1), ezl(2), ezl(3)
      end do
      endif
c procedure for a bending (dipole) element
      if(icat .eq. 2) then
      if(ldb) write(lfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,
     * nt1(k),nt2(k)
         if(imax.eq.0)goto 137
      if(ldb) write(lfile,607)(pmenu(i+mpp(k)),i=1,imax)
      angr=angd*pi180
      rho=dabs(brho/b)
ctm   patch for negative bends
      aleng=dabs(rho*angr)
      if(ldb) write(lfile,*) ' length = ', aleng
      if(ldb) write(lfile,*) ' distance = ', ss
      xx = ss
      yy = xx
      if(ldb) write(lfile,611) ss, aleng, wide, xx, yy, lmnlbl(k)
      ndptst=ndpts
      if(ndpts .le. 1) ndptst=2
      daleng=aleng/float(ndptst)
      dangr=angr/float(ndptst)
      k1 = nt1(k)
      k2 = nt2(k)
      do i=0,ndptst
      sst = ss + daleng*float(i)
      ttt = tt + daleng*float(i)*vr
      angrt=dangr*float(i)
      sfact=rho*dsin(angrt)
      cfact=dabs(rho*(1.d0 - dcos(angrt)))
c Note: signs have been adjusted to take into account that
c the rotation axis is - eyl.
      xt = x -cfact*exl(1) + sfact*ezl(1)
      yt = y -cfact*exl(2) + sfact*ezl(2)
      zt = z -cfact*exl(3) + sfact*ezl(3)
      angrt=-angrt
      call rotv(eyl,angrt,exl,exlt)
      call rotv(eyl,angrt,ezl,ezlt)
      if(ldb) write(lfile,*) ' when i = ',i
      if(ldb) write(lfile,*) ' values of sst and ttt '
      if(ldb) write(lfile,*) sst,ttt
      if(ldb) write(lfile,*) ' values of xt,yt,zt '
      if(ldb) write(lfile,*) xt,yt,zt
      write(35,235) lmnlbl(k),ltc(k1,k2),jk1,icat,k1,k2,xt,yt,zt,wide
      if(ldb) write(lfile,*) ' vectors exl,eyl,ezl '
      if(ldb) write(lfile,*) exlt(1), exlt(2), exlt(3)
      if(ldb) write(lfile,*) eyl(1), eyl(2), eyl(3)
      if(ldb) write(lfile,*) ezlt(1), ezlt(2), ezlt(3)
      end do
      endif
c
c write out results
ctm   write(6,*) 'geom output:',jflag,ifile,kfile,'=jflag,ifile,kfile.'
      if (jflag .eq. 1) then
      write(ifile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,
     * nt1(k),nt2(k)
c  605 format(1h ,1x,a8,1x,a8,1x,i5,1x,i5,1x,i5)
         if(imax.eq.0)goto 136
      write(ifile,607)(pmenu(i+mpp(k)),i=1,imax)
c  607 format((1h ,3(1x,d22.15)))
c Output using Mottershead's favorite pg format
c  607 format((1h ,3(1x,1pg22.15)))
      if(ldb) write(lfile,*) ' length = ', aleng
      if(ldb) write(lfile,*) ' distance = ', ss
      xx = ss
      yy = xx
c      write(kfile,611) ss, aleng, wide, xx, yy, lmnlbl(k)
      k1 = nt1(k)
      k2 = nt2(k)
      write(kfile,621) jk1,k,k1,k2,ss,aleng,wide,lmnlbl(k),ltc(k1,k2)
 621  format(4i4,2f13.6,f6.1,3x,a,3x,a)
      endif
c
      endif
c
c update geometry
c
c procedure for a straight element
      if(icat .eq. 1) then
      ss = ss + aleng
      tt = tt + aleng*vr
      x = x + aleng*ezl(1)
      y = y + aleng*ezl(2)
      z = z + aleng*ezl(3)
      endif
c procedure for a bending (dipole) element
      if(icat .eq. 2) then
      angr=angd*pi180
      rho=brho/b
ctm          aleng=rho*angr
ctm   patch for negative bends
      aleng=dabs(rho*angr)
      ss = ss + aleng
      tt = tt + aleng*vr
      sfact=rho*dsin(angr)
      cfact=dabs(rho*(1.d0 - dcos(angr)))
c Note: signs have been adjusted to take into account that
c the rotation axis is - eyl.
      x = x -cfact*exl(1) + sfact*ezl(1)
      y = y -cfact*exl(2) + sfact*ezl(2)
      z = z -cfact*exl(3) + sfact*ezl(3)
      angr=-angr
      call rotv(eyl,angr,exl,exl)
      call rotv(eyl,angr,ezl,ezl)
      endif
c procedure for an arot
      if(icat .eq. 3) then
      angr=angd*pi180
c Note: sign of angr has been adjusted to ...
      angr=-angr
      call rotv(ezl,angr,exl,exl)
      call rotv(ezl,angr,eyl,eyl)
      endif
c
 136  continue
c
c look for a data point
      if (nt2(k) .eq. 23) then
      write(ifile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,nt1(k),nt2(k)
      if(ldb) write(lfile,605)lmnlbl(k),ltc(nt1(k),nt2(k)),imax,
     * nt1(k),nt2(k)
      if(ldb) write(lfile,*) ' distance, time = ', ss,tt
      if(ldb) write(lfile,*) ' xyz =', x,y,z
      if(ldb) write(lfile,*) ' exl = ', exl(1),exl(2),exl(3)
      if(ldb) write(lfile,*) ' eyl = ', eyl(1),eyl(2),eyl(3)
      if(ldb) write(lfile,*) ' ezl = ', ezl(1),ezl(2),ezl(3)
      if(ldb) write(lfile,*) ' xmin,xmax = ',xmin,xmax
      write (jfile,520) ss,x,y,z,tt,0.
  520 format(6(1x,1pe12.5))
      endif
      endif
      endif
c procedure for a lump
      if(itype.eq.3) then
      write(ifile,515) string(1),mim(jk1)
  515 format(1x,1x,a8,1x,'lump',9x,'0',4x,'-1',2x,i4)
      endif
  137 continue
c
      return
      end
c
************************************************************************
c
