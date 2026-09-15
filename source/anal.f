***********************************************************************
* header              ANALYSIS (advanced commands)                    *
*  Routines for advanced commands and advanced analysis               *
***********************************************************************
c
      subroutine amap(p,fa,fm)
c subroutine for applying a map to a function or a set of moments
c Written by Alex Dragt, Fall 1986
c 10/29/91 job 3 and 4 code by Filippo Neri, put in by Johannes van Zeijts
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6)
      dimension fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms5),gm(6,6)
      dimension ha(monoms5),hm(6,6)
      dimension t1a(monoms5)
      dimension t2a(monoms5)
      character*3  kynd
c
c set up control indices
      job=nint(p(1))
      isend=nint(p(2))
      ifile=nint(p(3))
      nopt=nint(p(4))
      nskip=nint(p(5))
      nmpo=nint(p(6))
c
c procedure for reading in function or moments when job= 0, 1, or 2:
      if (job.eq.0 .or. job.eq.1 .or. job.eq.2) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      call ident(ga,gm)
      kynd='gtm'
      call strget(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
c
c procedure for reading in moments when job = 3 or 4:
      if (job.eq.3 .or. job.eq.4) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      kynd='gtm'
	call strget5(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin5(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
      continue
c
c procedure for letting map act on a function:
      if(job.eq.0) then
c let map characterized by fa,fm act on ga
c the result is the array ha
c this amounts to computing ha = (Dtranspose)*ga
c where D = D(fa,fm)
      call ident(ha,hm)
      call fxform(fa,fm,ga,ha)
      endif
c
c procedure for letting map act on moments:
c letting the map act on moments amounts to computing ha = D*ga
c
c procedure for job = 1 or 2:
      if (job.eq.1 .or. job.eq.2) then
c
c clear arrays
      call ident(ha,hm)
      do 10 i=1,monoms
   10 t1a(i)=0.d0
c
c when job = 2, compute all transformed moments through 4'th moments
      imax=209
c when job = 1, compute transformed moments only through 2'nd moments
      if (job.eq.1) imax=27
c
c perform calculation
      do 20 i=7,imax
      t1a(i)=1.d0
      call fxform(fa,fm,t1a,t2a)
      t1a(i)=0.d0
      do 30 j=1,209
   30 ha(i)=ha(i)+t2a(j)*ga(j)
   20 continue
c put result in buffer 1
      call mapmap(ha,hm,buf1a,buf1m)
      endif
c
c procedure for job = 3 or 4:
      if (job.eq.3 .or. job.eq.4) then
c when job = 4, compute all transformed moments through 4'th moments
      imax=209
c when job = 3, compute transformed moments only through 2'nd moments
      if (job.eq.3) imax=27
c clear arrays
      call ident(ha,hm)
      do 15 i=1,monoms5
   15 t1a(i)=0.d0
c perform calculation
      do 25 i=7,imax
      t1a(i)=1.d0
      call fxfrm5(fa,fm,t1a,t2a)
      t1a(i)=0.d0
      do 35 j=1,923
   35 ha(i)=ha(i)+t2a(j)*ga(j)
c      if(ha(i).ne.0) write(5,*) i,ha(i)
   25 continue
      continue
c put result in buffer 1, but to 5'th order
      call mapmap5(ha,hm,buf1a,buf1m)
      endif
c
c write out results if desired
c code needs to be modified to write out transformed function or moments
      if (job.eq.0) then
      if (isend.eq.1.or.isend.eq.3) then
      write(6,*) 'ready to write out transformed function'
c
c      write (jof,*) 'xe2=',xemit2
c      write (jof,*) 'ye2=',yemit2
c      write (jof,*) 'te2=',temit2
c      write (jof,*) 'xee2=',wex
c      write (jof,*) 'yee2=',wey
c      write (jof,*) 'tee2=',wet
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write(6,*) 'ready to write out transformed function'
c      write (jodf,*) 'xe2=',xemit2
c      write (jodf,*) 'ye2=',yemit2
c      write (jodf,*) 'te2=',temit2
c      write (jodf,*) 'xee2=',wex
c      write (jodf,*) 'yee2=',wey
c      write (jodf,*) 'tee2=',wet
      endif
      endif
c
      if (job.ne.0) then
      if (isend.eq.1.or.isend.eq.3) then
      write(6,*) 'ready to write out transformed moments'
c
c      write (jof,*) 'xe2=',xemit2
c      write (jof,*) 'ye2=',yemit2
c      write (jof,*) 'te2=',temit2
c      write (jof,*) 'xee2=',wex
c      write (jof,*) 'yee2=',wey
c      write (jof,*) 'tee2=',wet
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write(6,*) 'ready to write out transformed moments'
c      write (jodf,*) 'xe2=',xemit2
c      write (jodf,*) 'ye2=',yemit2
c      write (jodf,*) 'te2=',temit2
c      write (jodf,*) 'xee2=',wex
c      write (jodf,*) 'yee2=',wey
c      write (jodf,*) 'tee2=',wet
      endif
      endif
c
c write out result ha if nmpo > 0
      mpot=mpo
      mpo=nmpo
      if (nmpo.gt.0) call mapout(0,ha,hm)
      mpo=mpot
c
      return
      end
c
***********************************************************************
c
      subroutine asni(p)
c  This subroutine applies powers of script N inverse to phase space data.
c  It does this using analytic formulas.
c
      include 'impli.inc'
      include 'param.inc'
      include 'rays.inc'
c
c  calling array
      dimension p(6)
c
c  local array
      dimension fa(monoms), fr(monoms)
      dimension fm(6,6)
      character*3  kynd
c
c  set up control parameters
      iopt=nint(p(1))
      nmap=nint(p(2))
      nfcf=nint(p(3))
      istart=nint(p(4))
      igroup=nint(p(5))
      nwrite=nint(p(6))
c
c  begin calculation
c
c  get script N from storage
      kynd='gtm'
      call strget(kynd,nmap,fa,fm)
c
c  procedure for a static map
c
      if( iopt.eq.1) then
c  compute linear phase advances and linear time of flight
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
      wt=fm(5,6)
c
c  transform nonlinear part of map to the static resonance basis
      call ctosr(fa,fr)
c
c  begin outer loop over the sets of particles
      nset=nint(float(nrays)/float(igroup))
      do 10 i=1,nset
c  begin inner loop over the particles within a set
      do 20 j=1,igroup
c
c  get phase-space coordinates
      iray=(i-1)*igroup+j
      do 30 k=1,6
   30 zi(k)=zblock(iray,k)
c
c  compute emittances
      ex2=zi(1)**2 + zi(2)**2
      ey2=zi(3)**2 + zi(4)**2
      pt=zi(6)
c
c  compute phase advances and time of flight terms
c  incorporate - sign needed for inverse in the definition of an
      an=-float(istart+(i-1)*nwrite)
c  compute x and y phase advances
      phix = wx - 2.*pt*fr(28) - 2.*pt*pt*fr(84)
     & - 4.*ex2*fr(87) - 2.*ey2*fr(89)
      phiy = wy - 2.*pt*fr(29) - 2.*pt*pt*fr(85)
     & - 4.*ey2*fr(88) - 2.*ex2*fr(89)
c  compute time-like drift terms
      drt = pt*wt - ex2*fr(28) - ey2*fr(29)
     & - 2.*pt*ex2*fr(84) -2.*pt*ey2*fr(85)
     & - 3.*pt*pt*fr(30) - 4.*pt*pt*pt*fr(86)
c
c  set up matrix quantities
      cx=cos(an*phix)
      sx=sin(an*phix)
      cy=cos(an*phiy)
      sy=sin(an*phiy)
      tof=an*drt
c
c  apply matrix to transverse coordinates
      zf(1)= cx*zi(1)+sx*zi(2)
      zf(2)=-sx*zi(1)+cx*zi(2)
      zf(3)= cy*zi(3)+sy*zi(4)
      zf(4)=-sy*zi(3)+cy*zi(4)
c  transform time deviation and energy deviation
      zf(5)= zi(5)+tof
      zf(6)= zi(6)
c
c  write out results
      write (nfcf,100) zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
  100 format(6(1x,1pe12.5))
c
   20 continue
   10 continue
c
      endif
c
c  procedure for a dynamic map
c
      if( iopt.eq.2) then
c  compute linear phase advances
      cwx=fm(1,1)
      swx=fm(1,2)
      wx=atan2(swx,cwx)
      cwy=fm(3,3)
      swy=fm(3,4)
      wy=atan2(swy,cwy)
      cwt=fm(5,5)
      swt=fm(5,6)
      wt=atan2(swt,cwt)
c
c  transform nonlinear part of map to the dynamic resonance basis
      call ctodr(fa,fr)
c
c  begin outer loop over the sets of particles
      nset=nint(float(nrays)/float(igroup))
      do 40 i=1,nset
c  begin inner loop over the particles within a set
      do 50 j=1,igroup
c
c  get phase-space coordinates
      iray=(i-1)*igroup+j
      do 60 k=1,6
   60 zi(k)=zblock(iray,k)
c
c  compute emittances
      ex2=zi(1)**2 + zi(2)**2
      ey2=zi(3)**2 + zi(4)**2
      et2=zi(5)**2 + zi(6)**2
c
c  compute phase advances
c  incorporate - sign needed for inverse in the definition of an
      an=-float(istart+(i-1)*nwrite)
c  this part of code not yet complete
      phix=wx
      phiy=wy
      phit=wt
c
c  set up matrix quantities
      cx=cos(an*phix)
      sx=sin(an*phix)
      cy=cos(an*phiy)
      sy=sin(an*phiy)
      ct=cos(an*phit)
      st=sin(an*phit)
c
c  apply matrix to coordinates
      zf(1)= cx*zi(1)+sx*zi(2)
      zf(2)=-sx*zi(1)+cx*zi(2)
      zf(3)= cy*zi(3)+sy*zi(4)
      zf(4)=-sy*zi(3)+cy*zi(4)
      zf(5)= ct*zi(5)+st*zi(6)
      zf(6)=-st*zi(5)+ct*zi(6)
c
c  write out results
      write (nfcf,100) zf(1),zf(2),zf(3),zf(4),zf(5),zf(6)
c
   50 continue
   40 continue
c
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine betmap(ana,anm,ba,bm)
c This is a subroutine for finding the betatron portion of a map.
c The map ana, anm is assumed to be a map about the fixed point.
c Written by Alex Dragt, 4 August 1986
      include 'impli.inc'
      include 'param.inc'
c
c Calling arrays
      dimension ana(monoms),anm(6,6)
      dimension ba(monoms),bm(6,6)
c
c Local arrays
      dimension tempa(monoms),tempm(6,6)
c
c Extraction of betatron term b.
c First pass: terms linear in pt.
      call clear(tempa,tempm)
      call matmat(anm,tempm)
      tempa(33)=ana(33)
      tempa(38)=ana(38)
      tempa(42)=ana(42)
      tempa(45)=ana(45)
      tempa(53)=ana(53)
      tempa(57)=ana(57)
      tempa(60)=ana(60)
      tempa(67)=ana(67)
      tempa(70)=ana(70)
      tempa(76)=ana(76)
      call mapmap(tempa,tempm,ba,bm)
c Second pass: terms quadratic in pt.
      call inv(ba,bm)
      call concat(ba,bm,ana,anm,ba,bm)
      tempa(104)=ba(104)
      tempa(119)=ba(119)
      tempa(129)=ba(129)
      tempa(135)=ba(135)
      tempa(154)=ba(154)
      tempa(164)=ba(164)
      tempa(170)=ba(170)
      tempa(184)=ba(184)
      tempa(190)=ba(190)
      tempa(200)=ba(200)
      call mapmap(tempa,tempm,ba,bm)
      return
      end
c
***********************************************************************
c
      subroutine cex(p,ga,gm)
c This routine computes exp(:f).
c The array f=fa in extalk is destroyed in the process.
c Written by Alex Dragt, Spring 1987
 
      include 'impli.inc'
      include 'param.inc'
c
      include 'extalk.inc'
      include 'hmflag.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),ga(monoms),gm(6,6)
c
c Local arrays
      dimension em(6,6),fm(6,6)
      dimension y(224)
      character*3  kynd
c
c set up power and control indices
      power=p(1)
      nmapf=nint(p(2))
      nmapg=nint(p(3))
c
c get the array fa
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
c
c perform calculation
c
c set up exponent
      call csmul(power,fa,fa)
c
c compute a scaling factor to bring exponent within range of a
c taylor expansion or GENMAP
      call matify(em,fa)
      call mnorm(em,res)
      kmax=1
      scale=.5d0
   10 continue
      test=res*scale
      if (test.lt..1d0) goto 20
      kmax=kmax+1
      scale=scale/2.d0
      go to 10
   20 continue
c
c select procedure
      itest=1
      do 30 i=28,monoms
      if (fa(i).ne.0.d0) itest=2
      if (itest.ne.1) go to 40
   30 continue
   40 continue
      if (itest.eq.1) go to 50
      if (itest.eq.2) go to 80
c
c procedure using taylor series
   50 continue
      write (6,*) 'exp(:f:) computed using taylor series'
      call clear(fa,fm)
c rescale em
      call smmult(scale,em,em)
c compute taylor series result fm=exp(scale*em)
      call exptay(em,fm)
c raise the result to the 2**kmax (= 1/scale) power
      do 60 i=1,kmax
      call mmult(fm,fm,fm)
   60 continue
      goto 200
c
c procedure using genmap
   80 continue
      write(6,*) 'exp(:f:) computed using GENMAP'
c rescale fa
      call csmul(scale,fa,fa)
c setup and initialize for GENMAP routines
      iflag=1
      t=0.d0
      ns=50.d0
      h=.02d0
      ne=224
      do 90 i=1,ne
   90 y(i)=0.d0
      do 100 i=1,6
      j=7*i
  100 y(j)=1.d0
c call GENMAP routines
      call adam11(h,ns,'start',t,y)
      call putmap(y,fa,fm)
c
c raise the result to the 2**kmax (= 1/scale) power
      do 110 i=1,kmax
      call concat(fa,fm,fa,fm,fa,fm)
  110 continue
      go to 200
c
c decide where to put results
c
  200 continue
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='stm'
      call strget(kynd,nmapg,fa,fm)
      endif
c
      if (nmapg.eq.0) call mapmap(fa,fm,ga,gm)
c
      if (nmapg.eq.-1) call mapmap(fa,fm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(fa,fm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(fa,fm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(fa,fm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(fa,fm,buf5a,buf5m)
c
      return
      end
c
*******************************************************************************
c
      subroutine chrexp(iopt,delta,ta,tm,am1,am2,am3)
c
c This subroutine computes the chromatic expansion of the map ta,tm
c for the case in which the f3 and f4 parts of ta contain only terms
c linear and quadratic in pt, respectively.
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
c Calling arrays
      dimension ta(monoms),tm(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
c Local arrays
      dimension t1a(monoms),t1m(6,6),t2m(6,6)
c--------
c procedure when IOPT = 1 (delta=pt).
      if (iopt.eq.1) then
c Calculation of matrix associated with pt*f2 terms in ta.
c Set up f2a in t1a.
      call clear(t1a,am1)
      t1a(7)=ta(33)
      t1a(8)=ta(38)
      t1a(9)=ta(42)
      t1a(10)=ta(45)
      t1a(13)=ta(53)
      t1a(14)=ta(57)
      t1a(15)=ta(60)
      t1a(18)=ta(67)
      t1a(19)=ta(70)
      t1a(22)=ta(76)
c Compute matrix t1m corresponding to :f2a:.
      call matify(t1m,t1a)
c      write(6,*) 'result from chrexp'
c      call pcmap(1,0,0,0,t1a,t1m)
c Compute am1=t1m*tm
      call mmult(t1m,tm,am1)
c Calculation of matrix associated with (pt**2)*f2 terms in ta.
c Set up f2a in t1a.
      call clear(t1a,am2)
      t1a(7)=ta(104)
      t1a(8)=ta(119)
      t1a(9)=ta(129)
      t1a(10)=ta(135)
      t1a(13)=ta(154)
      t1a(14)=ta(164)
      t1a(15)=ta(170)
      t1a(18)=ta(184)
      t1a(19)=ta(190)
      t1a(22)=ta(200)
c Compute matrix t2m corresponding to :f2a:.
      call matify(t2m,t1a)
c      write(6,*) 'result from chrexp'
c      call pcmap(1,0,0,0,t1a,t2m)
c Compute am3=t2m+t1m*t1m/2.
      call mmult(t1m,t1m,am3)
      call smmult(.5d0,am3,am3)
      call madd(t2m,am3,am3)
c Compute am2=(t2m+t1m*t1m/2.)*tm
      call mmult(am3,tm,am2)
c Compute am3=tm+delta*am1+delta2*am2 for specific value of delta.
      delta2=delta*delta
      call matmat(tm,am3)
      call smmult(delta,am1,t1m)
      call madd(am3,t1m,am3)
      call smmult(delta2,am2,t1m)
      call madd(am3,t1m,am3)
      endif
c
c procedure when IOPT = 2 (delta=dp/p0).
      if (iopt.eq.2) then
      continue
      endif
c
      return
      end
c
*****************************************************************
c
      subroutine cod(p,th,tmh)
c This is a subroutine for computing closed orbit data.
c Written by Alex Dragt, 6 November 1985
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'parm.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6)
      dimension th(monoms),tmh(6,6)
c
c Local arrays
      dimension ana(monoms),anm(6,6)
      dimension ta(monoms),tm(6,6)
      dimension ba(monoms),bm(6,6)
      dimension ca(monoms),cm(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
c Set up control indices:
      iopt=nint(p(1))
      delta=p(2)
      idata=nint(p(3))
      ipmaps=nint(p(4))
      isend=nint(p(5))
      iwmaps=nint(p(6))
c
c Write headings
      if (isend.eq.1 .or. isend.eq.3) write(jof,90)
      if (isend.eq.2 .or. isend.eq.3) write(jodf,90)
   90 format(/,1x,'closed orbit analysis for static map')
c
c Computation of fixed point and map about it.
      call fxpt(th,tmh,ana,anm,ta,tm)
c
c Procedure for output of closed orbit location.
      if(idata.eq.1 .or. idata.eq.3) then
c Procedure when IOPT = 1:
      if (iopt.eq.1) then
c Compute location of closed orbit for given delta value
      delta2=delta*delta
      delta3=delta*delta2
      xc=delta*tm(1,6)-delta2*ta(63)-delta3*ta(174)
      pxc=delta*tm(2,6)+delta2*ta(48)+delta3*ta(139)
      yc=delta*tm(3,6)-delta2*ta(79)-delta3*ta(204)
      pyc=delta*tm(4,6)+delta2*ta(73)+delta3*ta(194)
c Put results in common/fitdat/ array
      dz(1)=tm(1,6)
      dz(2)=tm(2,6)
      dz(3)=tm(3,6)
      dz(4)=tm(4,6)
c Write out results
      do 5 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 5
      write(ifile,100)
  100 format(/,1x,'closed orbit data for delta defined in terms of',
     # 1x,'P sub tau:')
      write(ifile,120)
  120 format(1x,'location of closed orbit (x,px,y,py)')
      write(ifile,130)
  130 format(/,1x,'terms linear in delta')
      write(ifile,140) tm(1,6),tm(2,6),tm(3,6),tm(4,6)
  140 format(1x,4(d15.8,2x))
      write(ifile,150)
  150 format(/,1x,'terms quadratic in delta')
      write(ifile,140) -ta(63),ta(48),-ta(79),ta(73)
      write(ifile,160)
  160 format(/,1x,'terms cubic in delta')
      write(ifile,140) -ta(174),ta(139),-ta(204),ta(194)
      write(ifile,110) delta
  110 format(/,1x,'location of closed orbit when delta = ',d15.8)
      write(ifile,140) xc,pxc,yc,pyc
    5 continue
      endif
c Procedure when IOPT = 2
      if (iopt.eq.2) then
c Compute location of closed orbit for given delta value
c
c the code below needs to be modified
      delta2=delta*delta
      delta3=delta*delta2
      xc=delta*tm(1,6)-delta2*ta(63)-delta3*ta(174)
      pxc=delta*tm(2,6)+delta2*ta(48)+delta3*ta(139)
      yc=delta*tm(3,6)-delta2*ta(79)-delta3*ta(204)
      pyc=delta*tm(4,6)+delta2*ta(73)+delta3*ta(194)
c Put results in common/fitdat/ array
      dz(1)=tm(1,6)
      dz(2)=tm(2,6)
      dz(3)=tm(3,6)
      dz(4)=tm(4,6)
c end of code to be modified
c
c Write out results
      do 7 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 7
      write(ifile,102)
  102 format(/,1x,'closed orbit data for delta defined in terms of',
     # 1x,'momentum deviation:')
c
c the code below needs to be modified
      write(ifile,120)
      write(ifile,130)
      write(ifile,140) tm(1,6),tm(2,6),tm(3,6),tm(4,6)
      write(ifile,150)
      write(ifile,140) -ta(63),ta(48),-ta(79),ta(73)
      write(ifile,160)
      write(ifile,140) -ta(174),ta(139),-ta(204),ta(194)
      write(ifile,110) delta
      write(ifile,140) xc,pxc,yc,pyc
c end of code to be modified
c
    7 continue
      write(6,*) 'IOPT = 2 case not yet installed completely'
      endif
      endif
c
c Factorization of map into betatron and remaining terms.
c Computation of betatron term script B.
      call betmap(ana,anm,ba,bm)
c Computation of remaining nonlinear correction map script C.
c Use buf1a,buf1m as a temporary array
      call mapmap(ba,bm,buf1a,buf1m)
      call inv(buf1a,buf1m)
      call concat(buf1a,buf1m,ana,anm,ca,cm)
c Computation of B for specific value of delta.
      call chrexp(iopt,delta,ba,bm,am1,am2,am3)
c
c Procedure for output of twiss matrix and corrections.
      if(idata.eq.2 .or. idata.eq.3) then
c Procedure when IOPT =1
      if (iopt.eq.1) then
c Print out matrices bm, am1, and am2.
      do 9 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 9
      write(ifile,198)
  198 format(//,1x,'twiss matrix for delta defined in terms of',
     #1x,'P sub tau:')
      write(ifile,200)
  200 format(/,1x,'on energy twiss matrix')
      call pcmap(i,0,0,0,ba,bm)
      write(ifile,300)
  300 format(//,1x,'matrix for delta correction')
      call pcmap(i,0,0,0,ba,am1)
      write(ifile,400)
  400 format(//,1x,'matrix for delta**2 correction')
      call pcmap(i,0,0,0,ba,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
  402 format(//,1x,'value of twiss matrix when delta= ',d15.8)
      call pcmap(i,0,0,0,ba,am3)
    9 continue
      endif
c Procedure when IOPT = 2
      if (iopt.eq.2) then
c Print out matrices bm, am1, and am2.
      do 11 i=1,2
      if (i.eq.1) then
      iflag=1
      ifile=jof
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 11
      write(ifile,199)
  199 format(//,1x,'twiss matrix for delta defined in terms of',
     #1x,'momentum deviation:')
      write(ifile,202)
  202 format(/,1x,'on momentum twiss matrix')
      call pcmap(i,0,0,0,ba,bm)
      write(ifile,300)
      call pcmap(i,0,0,0,ba,am1)
      write(ifile,400)
      call pcmap(i,0,0,0,ba,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
      call pcmap(i,0,0,0,ba,am3)
   11 continue
      endif
      endif
c
c Procedure for output of the maps BC, B, and C.
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      do 13 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 13
      write(ifile,600)
  600 format(//,1x,'total transfer map about the closed orbit')
      call pcmap(i,i,0,0,ana,anm)
      write(ifile,700)
  700 format(//,1x,'betatron factor of transfer map')
      call pcmap(i,i,0,0,ba,bm)
      write(ifile,800)
  800 format(//,1x,'nonlinear factor of transfer map')
      call pcmap(i,i,0,0,ca,cm)
   13 continue
      endif
c
c Procedure for output of script T.
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      do 15 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 15
      write(ifile,900)
  900 format(//,1x,'transfer map script T to the closed orbit')
      call pcmap(i,i,0,0,ta,tm)
   15 continue
      endif
c
c Put maps in buffers
      call mapmap(ana,anm,buf1a,buf1m)
      call mapmap(ba,bm,buf2a,buf2m)
      call mapmap(ca,cm,buf3a,buf3m)
      call clear(buf4a,buf4m)
      call mapmap(buf4a,am3,buf4a,buf4m)
      call mapmap(ta,tm,buf5a,buf5m)
c
c Procedure for writing of maps.
      if(iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ana,anm)
      call mapout(0,ba,bm)
      call mapout(0,ca,cm)
      call mapout(0,buf4a,buf4m)
      call mapout(0,ta,tm)
      mpo=mpot
      endif
      return
      end
c
***********************************************************************
c
      subroutine csym(isend,fm,ans)
c
c This subroutine checks the symplectic condition for the matrix fm
c Written by Alex Dragt, 4 October 1989
c
      include 'impli.inc'
      include 'files.inc'
c
c Calling arrays
      dimension fm(6,6)
c
c Temporary arrays
      dimension tm(6,6)
c
c-----Procedure
c
      call matmat(fm,tm)
      call minv(tm)
      call mmult(tm,fm,tm)
      do 10 i=1,6
   10 tm(i,i)=tm(i,i)-1.d0
      call mnorm(tm,ans)
c
c Write out results if desired
c
      if (isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' symplectic violation = ',ans
      endif
      if (isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' symplectic violation = ',ans
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine dia(p,fa,fm)
c this is a routine for computing invariants in the dynamic case
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'rays.inc'
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6)
c
c set up control indices
      iopt=nint(p(1))
      ipinv=nint(p(2))
      ipmaps=nint(p(3))
      isend=nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'dynamic invariant analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'dynamic invariant analysis'
      endif
c
c begin calculation
c
c find the transforming (conjugating) map script A
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm)
c remove offensive terms from f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put script A in buffer 1 and script N in buffer 2
      call mapmap(ta,tm,buf1a,buf1m)
      call mapmap(ga,gm,buf2a,buf2m)
c
c procedure for computing invariant
c clear buffers 3, 4, and 5
      call ident(buf3a,buf3m)
      call ident(buf4a,buf4m)
      call ident(buf5a,buf5m)
c invert script A
      call inv(ta,tm)
c
c computation of regular invariants
      if (iopt .eq. 0) then
c computation of x invariant
      call clear(t1a,t1m)
      t1a(7)=1.d0
      t1a(13)=1.d0
c put invariant in buffer 3
      call fxform(ta,tm,t1a,buf3a)
c computation of y invariant
      call clear(t1a,t1m)
      t1a(18)=1.d0
      t1a(22)=1.d0
c put invariant in buffer 4
      call fxform(ta,tm,t1a,buf4a)
c computation of t invariant
      call clear(t1a,t1m)
      t1a(25)=1.d0
      t1a(27)=1.d0
c put invariant in buffer 5
      call fxform(ta,tm,t1a,buf5a)
      endif
c
c computation of mixed invariant
c NOTE: THE CASE IOPT < 0 STILL HAS TO BE PROGRAMMED
      if(iopt .gt. 0) then
      read(iopt,*) anx,any,ant
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      call clear(t1a,t1m)
      t1a(7)=anx
      t1a(13)=anx
      t1a(18)=any
      t1a(22)=any
      t1a(25)=ant
      t1a(27)=ant
c put invariant in buffer 3
      call fxform(ta,tm,t1a,buf3a)
      endif
c
c procedure for putting out invariants
      if (ipinv .eq. 1) then
c putting out regular invariants
      if (iopt .eq. 0) then
      do 10 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 10
      write(ifile,*)
      write(ifile,*) 'x invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
      write(ifile,*)
      write(ifile,*) 'y invariant polynomial'
      call pcmap(0,j,0,0,buf4a,buf4m)
      write(ifile,*)
      write(ifile,*) 't invariant polynomial'
      call pcmap(0,j,0,0,buf5a,buf5m)
   10 continue
      endif
c putting out mixed invariant
      if (iopt .ne. 0) then
      do 15 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 15
      write(ifile,*)
      write(ifile,*) 'mixed invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
   15 continue
      endif
      endif
c
c procedure for printing out maps
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form map script N'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   20 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      if (iopt .eq. 0) then
      call mapout(0,buf3a,buf3m)
      call mapout(0,buf4a,buf4m)
      call mapout(0,buf5a,buf5m)
      endif
      if (iopt .ne. 0) then
      call mapout(0,buf3a,buf3m)
      endif
      mpo=mpot
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine dnor(p,fa,fm)
c this is a subroutine for normal form analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6)
c
c set up control indices
      keep=  nint(p(1))
      idata= nint(p(2))
      ipmaps=nint(p(3))
      isend= nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'dynamic normal form analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'dynamic normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
ctm04 call dpur2(fa,fm,ga,gm,ta,tm,t1m)
      call dpur2(fa,fm,ga,gm,ta,tm)
c remove offensive terms from f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(ga,gm,buf2a,buf2m)
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=gm(1,1)
      swx=gm(1,2)
      wx=atan2(swx,cwx)
      cwy=gm(3,3)
      swy=gm(3,4)
      wy=atan2(swy,cwy)
      cwt=gm(5,5)
      swt=gm(5,6)
      wt=atan2(swt,cwt)
c set up normal form for exponent
      do 10 i=1,27
   10 ga(i)=0.
      ga(7)=-wx/2.d0
      ga(13)=-wx/2.d0
      ga(18)=-wy/2.d0
      ga(22)=-wy/2.d0
      ga(25)=-wt/2.d0
      ga(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or. idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,ga,g1a)
c store results in buffer 3
      call mapmap(g1a,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,ga,gm)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
*******************************************************************
c
      subroutine fadm(p,fa,fm)
c  this subroutine fourier analyzes a dynamic transfer map
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'fadm not yet available'
      return
      end
c
*******************************************************************
c
      subroutine fasm(p,fa,fm)
c  this subroutine fourier analyzes a static transfer map
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'fasm not yet available'
      return
      end
c
*******************************************************************
c
      subroutine gbuf(p,fa,fm)
c  this subroutine gets a map from an auxiliary buffer and
c  either concatenates it with the map in the main buffer,
c  or uses it to replace the map in the main buffer.
c  Written by Alex Dragt, Spring 1987
c  Modified by Alex Dragt, 17 June 1988
c  Modified by Alex Dragt, 13 October 1988
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c  Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c set up and test control indices
      iopt=nint(p(1))
      i=nint(p(2))
      if (i.lt.1 .or. i.gt.5) then
      write(6,*) 'trouble with index nmap in gbuf'
      return
      endif
c
c if iopt=1, concatenate the existing map with the map in bufi
      if(iopt.eq.1) then
      if(i.eq.1) call concat(fa,fm,buf1a,buf1m,fa,fm)
      if(i.eq.2) call concat(fa,fm,buf2a,buf2m,fa,fm)
      if(i.eq.3) call concat(fa,fm,buf3a,buf3m,fa,fm)
      if(i.eq.4) call concat(fa,fm,buf4a,buf4m,fa,fm)
      if(i.eq.5) call concat(fa,fm,buf5a,buf5m,fa,fm)
      return
      endif
c
c if iopt=2, replace the existing map with the map in bufi
      if(iopt.eq.2) then
      if(i.eq.1) call mapmap(buf1a,buf1m,fa,fm)
      if(i.eq.2) call mapmap(buf2a,buf2m,fa,fm)
      if(i.eq.3) call mapmap(buf3a,buf3m,fa,fm)
      if(i.eq.4) call mapmap(buf4a,buf4m,fa,fm)
      if(i.eq.5) call mapmap(buf5a,buf5m,fa,fm)
      return
      endif
c
      write(6,*) 'trouble with index iopt in gbuf'
      return
      end
c
c*****************************************************************************
c
c      subroutine geom(pp)
c
c
************************************************************************
c
      subroutine hmltn1(h)
c This subroutine is used to specify a constant hamiltonian.
c Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'param.inc'
      include 'extalk.inc'
      include 'hmflag.inc'
      dimension h(monoms)
c
c  begin computation
      iflag=0
      do 10 i=1,monoms
   10 h(i)=-fa(i)
c
      return
      end
c
********************************************************************
c
      subroutine moma(p)
c this subroutine is for moment and map analysis
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6)
      dimension fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms5),gm(6,6)
      dimension ha(monoms5),hm(6,6)
      dimension t1a(monoms5)
      dimension t2a(monoms5)
      character*3  kynd
c
c set up control indices
      job=nint(p(1))
      isend=nint(p(2))
      ifile=nint(p(3))
      nopt=nint(p(4))
      nskip=nint(p(5))
      nmpo=nint(p(6))
c
c procedure for reading in function or moments when job = 1,2, or 3:
      if (job.eq.1 .or. job.eq.2 .or. job.eq.3) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      call ident(ga,gm)
      kynd='gtm'
      call strget(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
c
c procedure for reading in moments when job = 4,5, or 6:
      if (job.eq.4 .or. job.eq.5 .or. job.eq.6) then
c test for file read or internal map fetch
      if(ifile.lt.0) then
      nmap=-ifile
      kynd='gtm'
	call strget5(kynd,nmap,ga,gm)
      else
      mpit=mpi
      mpi=ifile
      call mapin5(nopt,nskip,ga,gm)
      mpi=mpit
      endif
      endif
      continue
c
c compute eigen emittances
      if (job.eq.1 .or. job.eq.4) call eigemt(2,ga)
      if (job.eq.2 .or. job.eq.5) call eigemt(4,ga)
      if (job.eq.3 .or. job.eq.6) call eigemt(6,ga)
c compute mean square eigen-emittances and put results in fitbuf
      wex=buf1a(7)**2
      wey=buf1a(18)**2
      wet=buf1a(25)**2
c
c write out results if desired
c code needs to be modified to write out eigen emittances
      if (isend.eq.1.or.isend.eq.3) then
      write (jof,*) 'xee2=',wex
      write (jof,*) 'yee2=',wey
      write (jof,*) 'tee2=',wet
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write (jodf,*) 'xee2=',wex
      write (jodf,*) 'yee2=',wey
      write (jodf,*) 'tee2=',wet
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine padd(p,ha,hm)
c this subroutine adds two polynomials
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      character*3  kynd
c
c set up control indices
      nmapf=nint(p(1))
      nmapg=nint(p(2))
      nmaph=nint(p(3))
c
c  clear arrays and get maps
      call ident(ta,tm)
      if (nmapf.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
      if (nmapg.eq.0) call mapmap(ha,hm,ga,gm)
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='gtm'
      call strget(kynd,nmapg,ga,gm)
      endif
c
c perform calculation
      call cpadd(fa,ga,ta)
c
c decide where to put results
c
      if (nmaph.ge.1 .and. nmaph.le.5) then
      kynd='stm'
      call strget(kynd,nmaph,ta,tm)
      endif
c
      if (nmaph.eq.0) call mapmap(ta,tm,ha,hm)
c
      if (nmaph.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmaph.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmaph.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmaph.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmaph.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
********************************************************************
c
      subroutine pbpol(p,fa,fm)
c this subroutine poisson brackets two polynomials
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write (6,*) 'in subroutine pbpol'
c
      return
      end
c
**********************************************************************
c
      subroutine pdnf(p,ha,hm)
c this subroutine computes powers of a dynamic normal form
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ta(monoms)
      dimension fm(6,6),tm(6,6)
      character*3  kynd
c
c set up control indices
      jinopt=nint(p(1))
      if (jinopt.eq.1) pow=p(2)
      if (jinopt.eq.2) npowf=nint(p(2))
      nmapi=nint(p(3))
      joutop=nint(p(4))
      nmapo=nint(p(5))
c
c clear array and get map
      call ident(ta,tm)
      if (nmapi.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapi.ge.1 .and. nmapi.le.5) then
      kynd='gtm'
      call strget(kynd,nmapi,fa,fm)
      endif
c
c perform calculation
c
c procedure when jinopt=1
      if (jinopt.eq.1) then
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
c
c procedure when jinopt=2
      if (jinopt.eq.2) then
c
c procedure when npowf .lt. 0
      ipset=-npowf
      if (ipset.gt.0 .and. ipset.le.maxpst) then
      do 10 k=1,6
      pow=pst(k,ipset)
      if (pow.ne.0.d0) then
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
   10 continue
      endif
c procedure when npowf .gt.0
      if (npowf.gt.0) then
   20 continue
      read(npowf,*,end=30) pow
      call cpdnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      goto 20
   30 continue
      endif
c
      endif
c
      return
      end
c
********************************************************************
c
      subroutine pmul(p,ha,hm)
c this subroutine multiplies two polynomials
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ga(monoms),ta(monoms)
      dimension fm(6,6),gm(6,6),tm(6,6)
      character*3  kynd
c
c set up control indices
      nmapf=nint(p(1))
      nmapg=nint(p(2))
      nmaph=nint(p(3))
c
c clear array and get maps
      call ident(ta,tm)
      if (nmapf.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
      if (nmapg.eq.0) call mapmap(ha,hm,ga,gm)
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='gtm'
      call strget(kynd,nmapg,ga,gm)
      endif
c
c perform calculation
      call cpmul(fa,ga,ta)
c
c decide where to put results
c
      if (nmaph.ge.1 .and. nmaph.le.5) then
      kynd='stm'
      call strget(kynd,nmaph,ta,tm)
      endif
c
      if (nmaph.eq.0) call mapmap(ta,tm,ha,hm)
c
      if (nmaph.eq.-1) call mapmap(ta,tm,buf1a,buf1m)
      if (nmaph.eq.-2) call mapmap(ta,tm,buf2a,buf2m)
      if (nmaph.eq.-3) call mapmap(ta,tm,buf3a,buf3m)
      if (nmaph.eq.-4) call mapmap(ta,tm,buf4a,buf4m)
      if (nmaph.eq.-5) call mapmap(ta,tm,buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine pnlp(p,ga,gm)
c this subroutine raises the nonlinear part of a map to a power
c Written by Alex Dragt, Fall 1988
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),ga(monoms),gm(6,6)
c
c Local arrays
      dimension fa(monoms)
      dimension fm(6,6)
      character*3  kynd
c
c set up scalar and control indices
      iopt=nint(p(1))
      power=p(2)
      nmapf=nint(p(3))
      nmapg=nint(p(4))
c
c get map
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
c
c perform calculation
      if(iopt.eq.0) then
      call mclear(fm)
      do 10 i=1,6
  10  fm(i,i)=1.d0
      endif
      call csmul(power,fa,fa)
c
c decide where to put results
c
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='stm'
      call strget(kynd,nmapg,fa,fm)
      endif
c
      if (nmapg.eq.0) call mapmap(fa,fm,ga,gm)
c
      if (nmapg.eq.-1) call mapmap(fa,fm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(fa,fm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(fa,fm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(fa,fm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(fa,fm,buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine pold(p,fa,fm)
c this is a subroutine for polar decomposition
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ta(monoms)
      dimension tm(6,6),rm(6,6),pdsm(6,6),revec(6,6)
      dimension reval(6)
c
c set up control indices
      write (6,*) ' in pold'
      mapin=nint(p(1))
      isend=nint(p(2))
      idata=nint(p(3))
      ipmaps=nint(p(4))
      iwmaps=nint(p(5))
c
c write heading(s)
c
c begin calculation
      call polr(fa,fm,rm,pdsm,reval,revec)
c put out data
      call ident(ta,tm)
c
c put out maps
c
c put maps in buffers
      call ident(ta,tm)
      call mapmap(ta,rm,buf1a,buf1m)
      call mapmap(ta,pdsm,buf2a,buf2m)
      call mapmap(fa,tm,buf3a,buf3m)
      do 10 i=1,6
      tm(i,i)=reval(i)
   10 continue
      call mapmap(ta,tm,buf4a,buf4m)
      call mapmap(ta,revec,buf5a,buf5m)
c
c write out maps
c
      return
      end
c
**********************************************************************
c
      subroutine ppa(p,fa,fm)
c
c     This routine computes focal lengths and principal planes for the
c     current map.
c         C. T. Mottershead  LANL AT-3 / 2 Oct 89
c---------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'fitdat.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6)
c
      isend = int(p(1))
      if(iquiet.eq.1) isend = 0
      eps = 1.e-9
c
c       x-plane
c
      finv = fm(2,1)
      if(abs(finv).lt.eps) finv = eps
      f = -1.0/finv
      z2 = f*(1.0 - fm(1,1))
      z1 = f*(1.0 - fm(2,2))
      fx = f
      xb = z1
      xa = z2
      xu = f - z1
      xd = f - z2
c
c       y-plane
c
      finv = fm(4,3)
      if(abs(finv).lt.eps) finv = eps
      f = -1.0/finv
      z2 = f*(1.0 - fm(3,3))
      z1 = f*(1.0 - fm(4,4))
      fy = f
      yb = z1
      ya = z2
      yu = f - z1
      yd = f - z2
c
c     print the matrix and focal lengths
c
      if(isend.lt.2) go to 400
      lun = jodf
 300  continue
      write(lun,17) fx,fy
  17  format(' * PPA  Focal lengths : fx =',1pg15.7,'  fy =',1pg15.7)
      write(lun,33)
  33  format(' * Principal Planes (before and after):')
      write(lun,37) xb,xa,yb,ya
  37  format(1x,' xb=',1pg15.7,' xa=',1pg15.7,' yb=',1pg15.7,' ya=',
     * 1pg15.7)
      write(lun,27) xu,yu
  27  format(' * Focal points (upstream): xu =',1pg15.7,
     * '  yu =',1pg15.7)
      write(lun,29) xd,yd
  29  format('              (downstream): xd =',1pg15.7,
     * '  yd =',1pg15.7)
 400  if((isend.eq.1).or.(isend.eq.3)) then
         lun = jof
         isend = 0
         go to 300
      endif
c
      return
      end
c
********************************************************************
c
      subroutine psnf(p,ha,hm)
c this subroutine computes powers of a static normal form
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ta(monoms)
      dimension fm(6,6),tm(6,6)
      character*3  kynd
c
c set up control indices
      jinopt=nint(p(1))
      if (jinopt.eq.1) pow=p(2)
      if (jinopt.eq.2) npowf=nint(p(2))
      nmapi=nint(p(3))
      joutop=nint(p(4))
      nmapo=nint(p(5))
c
c get map and clear arrays
      if (nmapi.eq.0) call mapmap(ha,hm,fa,fm)
      if (nmapi.ge.1 .and. nmapi.le.5) then
      kynd='gtm'
      call strget(kynd,nmapi,fa,fm)
      endif
      call ident(ta,tm)
c
c perform calculation
c
c procedure when jinopt=1
      if (jinopt.eq.1) then
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
c
c procedure when jinopt=2
      if (jinopt.eq.2) then
c
c procedure when npowf .lt. 0
      ipset=-npowf
      if (ipset.gt.0 .and. ipset.le.maxpst) then
      do 10 k=1,6
      pow = pst(k,ipset)
      if (pow.ne.0.d0) then
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      endif
   10 continue
      endif
c procedure when npowf .gt.0
      if (npowf.gt.0) then
   20 continue
      read(npowf,*,end=30) pow
      call cpsnf(pow,fa,fm,ta,tm)
      call mapsnd(joutop,nmapo,ta,tm,ha,hm)
      goto 20
   30 continue
      endif
c
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine psp(p,ha,hm)
c this subroutine computes the scalar product of two polynomials
c Written by Alex Dragt, 10/23/89
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension fa(monoms),ga(monoms)
      dimension tm(6,6)
      character*3  kynd
c
c set up control indices
      job=nint(p(1))
      nmapf=nint(p(2))
      nmapg=nint(p(3))
      isend=nint(p(4))
c
c  clear arrays and get maps
      if (nmapf.eq.0) call mapmap(ha,hm,fa,tm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,tm)
      endif
      if (nmapg.eq.0) call mapmap(ha,hm,ga,tm)
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='gtm'
      call strget(kynd,nmapg,ga,tm)
      endif
c
c perform calculation
      call cpsp(fa,ga,ans1,ans2,ans3,ans4)
c
c decide where to send and put results
c
      write(6,*) ' ans1=',ans1,' ans2=',ans2
      write(6,*) ' ans3=',ans3,' ans4=',ans4
c
      return
      end
c
********************************************************************
c
      subroutine pval(p,ga,gm)
c this is a routine for evaluating a polynomial
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'rays.inc'
      include 'files.inc'
      include 'parset.inc'
c
c Calling arrays
      dimension p(6),ga(monoms),gm(6,6)
c
c Local arrays
      dimension fa(monoms),fm(6,6)
      dimension zit(6)
      character*3  kynd
c
c set up control indices
      mapin=nint(p(1))
      idata=nint(p(2))
      iwnum=nint(p(3))
c
c get polynomial
c
      if (mapin.eq.0) call mapmap(ga,gm,fa,fm)
      if (mapin.ge.1 .and. mapin.le.5) then
      kynd='gtm'
      call strget(kynd,mapin,fa,fm)
      endif
c
c compute value(s) of polynomial and write them out
c
      if (iwnum.gt.0) then
      write(jof,*)
      write(jof,*) 'value(s) of polymomial written on file ',iwnum
c
c procedure when idata > 0
      if (idata.gt.0) then
      ipset=idata
c get phase space data from the parameter set ipset
      if(ipset.lt.1 .or. ipset.gt.maxpst) then
       do 50 i=1,6
  50   zit(i) = 0.0d0
      else
       do 60 i=1,6
  60   zit(i) = pst(i,ipset)
      endif
c
c carry out computation
      call evalf(zit,fa,val2,val3,val4)
      k=1
      write(iwnum,100) k,val2,val3,val4,0.,0.
  100 format(1x,i12,5(1x,1pe12.5))
      endif
c
c procedure when idata = 0
      if (idata.eq.0) then
      do 70 k=1,nrays
c check on status of k'th ray; skip this ray if its status is 'lost'
      if (istat(k).ne.0) goto 70
      do 80 j=1,6
   80 zit(j)=zblock(k,j)
      call evalf(zit,fa,val2,val3,val4)
      write(iwnum,100) k,val2,val3,val4,0.,0.
   70 continue
      endif
c
      endif
c
      return
      end
c
**********************************************************************
c
      subroutine radm(p,fa,fm)
c this is a subroutine for resonance analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 17 June 1988
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ta(monoms),t1a(monoms)
      dimension tm(6,6),t1m(6,6)
      dimension look(3),pmask(6)
c
c set up control indices
      iopt=nint(p(1))
      i2=nint(p(2))
      i3=nint(p(3))
      i4=nint(p(4))
      iwmaps=nint(p(5))
c
c compute isend
      do 10 j=1,3
      look(j)=0
      if (i2.eq.j .or. i3.eq.j .or. i4.eq.j) look(j)=1
   10 continue
      isend=0
      if (look(1).eq.1) isend=1
      if (look(2).eq.1) isend=2
      if (look(1).eq.1 .and. look(2).eq.1) isend=3
      if (look(3).eq.1) isend=3
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*)
      write(jof,*) 'resonance analysis of dynamic map'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'resonance analysis of dynamic map'
      endif
c
c beginning of calculation
c
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,buf2a,buf2m,ta,tm)
c
c procedure for removing third order terms
      if(iopt.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'third order terms removed'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'third order terms removed'
      endif
      call dpur3(buf2a,buf2m,buf3a,buf3m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,ta,tm)
c put map in proper place
      call mapmap(buf3a,buf3m,buf2a,buf2m)
      endif
c
c resonance decompose purified map:
      call matmat(buf2m,buf3m)
      call ctodr(buf2a,buf3a)
c
c procedure for writing resonance driving terms at terminal (file jof)
      if (isend.eq.1 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.1 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.1 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.1 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jof,*)
      write(jof,*) 'requested resonance driving terms written as a map'
      call pdrmap(0,1,buf4a,buf4m)
      endif
c
c procedure for writing resonance driving terms on external file
c (file jodf)
      if (isend.eq.2 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.2 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.2 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.2 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jodf,*)
      write(jodf,*) 'requested resonance driving terms written as a map'
      call pdrmap(0,2,buf4a,buf4m)
      endif
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ta,tm)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c put the transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c buffers 2 and 3 already contain the purified map in the cartesian
c and dynamic resonance bases, respectively
c clear the remaining buffers
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine rasm(p,fa,fm)
c this is a subroutine for resonance analysis of static maps
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 17 June 1988
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ta(monoms),t1a(monoms)
      dimension tm(6,6),t1m(6,6)
      dimension look(3),pmask(6)
 
c
c set up control indices
      iopt=nint(p(1))
      i2=nint(p(2))
      i3=nint(p(3))
      i4=nint(p(4))
      iwmaps=nint(p(5))
c
c compute isend
      do 10 j=1,3
      look(j)=0
      if (i2.eq.j .or. i3.eq.j .or. i4.eq.j) look(j)=1
   10 continue
      isend=0
      if (look(1).eq.1) isend=1
      if (look(2).eq.1) isend=2
      if (look(1).eq.1 .and. look(2).eq.1) isend=3
      if (look(3).eq.1) isend=3
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*)
      write(jof,*) 'resonance analysis of static map'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'resonance analysis of static map'
      endif
c
c beginning of calculation
c
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,buf2a,buf2m,ta,tm,t1m)
c
c procedure for removing third order terms
      if(iopt.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'third order terms removed'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'third order terms removed'
      endif
c remove offensive chromatic terms from f3 part of map:
      call scpur3(buf2a,buf2m,buf3a,buf3m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(buf3a,buf3m,buf2a,buf2m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,ta,tm)
      endif
c
c resonance decompose purified map:
      call matmat(buf2m,buf3m)
      call ctosr(buf2a,buf3a)
c
c procedure for writing resonance driving terms at terminal (file jof)
      if (isend.eq.1 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.1 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.1 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.1 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jof,*)
      write(jof,*) 'requested resonance driving terms written as a map'
      call psrmap(0,1,buf4a,buf4m)
      endif
c
c procedure for writing resonance driving terms on external file (file jodf)
      if (isend.eq.2 .or. isend.eq.3) then
      call mapmap(buf3a,buf3m,buf4a,buf4m)
c compute masking parameters pmask(j)
      pmask(1)=1.
      pmask(2)=0.
      pmask(3)=0.
      pmask(4)=0.
      if (i2.eq.2 .or. i2.eq.3) pmask(2)=1.
      if (i3.eq.2 .or. i3.eq.3) pmask(3)=1.
      if (i4.eq.2 .or. i4.eq.3) pmask(4)=1.
c mask of unwanted portions of buf4a
      call mask(pmask,buf4a,buf4m)
c display result in sr basis
      write(jodf,*)
      write(jodf,*) 'requested resonance driving terms written as a map'
      call psrmap(0,2,buf4a,buf4m)
      endif
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,ta,tm)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c put the transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c buffers 2 and 3 already contain the purified map in the cartesian
c and static resonance bases, respectively
c clear the remaining buffers
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine sia(p,fa,fm)
c this is a routine for computing invariants in the static case
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 17 June 1988
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'rays.inc'
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6)
c
c set up control indices
      iopt=nint(p(1))
      ipinv=nint(p(2))
      ipmaps=nint(p(3))
      isend=nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'static invariant analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'static invariant analysis'
      endif
c
c begin calculation
c
c find the transforming (conjugating) map script A
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,ga,gm,ta,tm,t1m)
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put script A in buffer 1 and script N in buffer 2
      call mapmap(ta,tm,buf1a,buf1m)
      call mapmap(g1a,g1m,buf2a,buf2m)
c
c procedure for computing invariants
c clear buffers 3, 4, and 5
      call ident(buf3a,buf3m)
      call ident(buf4a,buf4m)
      call ident(buf5a,buf5m)
c
c invert script A
      call inv(ta,tm)
c
c computation of regular invariants
      if(iopt.eq.0) then
c computation of x invariant
      call clear(t1a,t1m)
      t1a(7)=1.d0
      t1a(13)=1.d0
c put invariant in buffer 3
      call fxform(ta,tm,t1a,buf3a)
c computation of y invariant
      call clear(t1a,t1m)
      t1a(18)=1.d0
      t1a(22)=1.d0
c put invariant in buffer 4
      call fxform(ta,tm,t1a,buf4a)
      endif
c
c computation of mixed invariant
c NOTE: THE CASE IOPT < 0 STILL HAS TO BE PROGRAMMED
      if(iopt .gt. 0) then
      read(iopt,*) anx,any
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jof,*) 'parameters read in from file',iopt
      endif
      call clear(t1a,t1m)
      t1a(7)=anx
      t1a(13)=anx
      t1a(18)=any
      t1a(22)=any
c put invariant in buffer 3
      call fxform(ta,tm,t1a,buf3a)
      endif
c
c procedure for putting out invariants
      if (ipinv .eq. 1) then
c putting out regular invariants
      if (iopt.eq.0) then
      do 10 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 10
      write(ifile,*)
      write(ifile,*) 'x invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
      write(ifile,*)
      write(ifile,*) 'y invariant polynomial'
      call pcmap(0,j,0,0,buf4a,buf4m)
   10 continue
      endif
c putting out mixed invariant
      if (iopt .ne. 0) then
      do 15 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 15
      write(ifile,*)
      write(ifile,*) 'mixed invariant polynomial'
      call pcmap(0,j,0,0,buf3a,buf3m)
   15 continue
      endif
      endif
c
c procedure for printing out maps
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form map script N'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   20 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      if (iopt .eq. 0) then
      call mapout(0,buf3a,buf3m)
      call mapout(0,buf4a,buf4m)
      endif
      if (iopt .ne. 0) then
      call mapout(0,buf3a,buf3m)
      endif
      mpo=mpot
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine smul(p,ga,gm)
c this subroutine multiplies a polynomial by a scalar
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 17 June 1988
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
      include 'parset.inc'
c
c Calling arrays
      dimension p(6),ga(monoms),gm(6,6)
c
c Local arrays
      dimension fa(monoms)
      dimension fm(6,6)
      character*3  kynd
c
c set up scalar and control indices
      iopt=nint(p(1))
      scalar=p(2)
      nmapf=nint(p(3))
      nmapg=nint(p(4))
      ipset=nint(p(5))
      islot=nint(p(6))
      if (ipset .gt. 0) scalar=pst(islot,ipset)
c
c get map
      if (nmapf.eq.0) call mapmap(ga,gm,fa,fm)
      if (nmapf.ge.1 .and. nmapf.le.5) then
      kynd='gtm'
      call strget(kynd,nmapf,fa,fm)
      endif
c
c perform calculation
      if(iopt.eq.0) then
      call mclear(fm)
      do 10 i=1,6
  10  fm(i,i)=1.d0
      endif
      call csmul(scalar,fa,fa)
c
c decide where to put results
c
      if (nmapg.ge.1 .and. nmapg.le.5) then
      kynd='stm'
      call strget(kynd,nmapg,fa,fm)
      endif
c
      if (nmapg.eq.0) call mapmap(fa,fm,ga,gm)
c
      if (nmapg.eq.-1) call mapmap(fa,fm,buf1a,buf1m)
      if (nmapg.eq.-2) call mapmap(fa,fm,buf2a,buf2m)
      if (nmapg.eq.-3) call mapmap(fa,fm,buf3a,buf3m)
      if (nmapg.eq.-4) call mapmap(fa,fm,buf4a,buf4m)
      if (nmapg.eq.-5) call mapmap(fa,fm,buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine snor(p,fa,fm)
c this is a subroutine for normal form analysis of static maps
c Written by Alex Dragt, Spring 1987
c Modified 20 June 1988
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'files.inc'
      include 'buffer.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6)
      dimension tm(6,6),t1m(6,6)
c
c set up control indices
      keep=  nint(p(1))
      idata= nint(p(2))
      ipmaps=nint(p(3))
      isend= nint(p(4))
      iwmaps=nint(p(5))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write (jof,*)
      write (jof,*) 'static normal form analysis'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write (jodf,*)
      write (jodf,*) 'static normal form analysis'
      endif
c
c begin calculation
c
c remove offensive terms from matrix part of map:
      call spur2(fa,fm,ga,gm,ta,tm,t1m)
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c remove offensive terms from f4 part of map:
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map:
      call concat(t1a,t1m,ta,tm,ta,tm)
c put transforming map in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c put transformed map in buffer 2
      call mapmap(g1a,g1m,buf2a,buf2m)
c
c procedure for computing normal form exponent and pseudo hamiltonian
      call ident(buf3a,buf3m)
      if (idata.eq.1 .or. idata.eq.2 .or. idata.eq.3) then
c compute phase advances
      cwx=g1m(1,1)
      swx=g1m(1,2)
      wx=atan2(swx,cwx)
      cwy=g1m(3,3)
      swy=g1m(3,4)
      wy=atan2(swy,cwy)
c compute momentum compaction
      wt=g1m(5,6)
c set up normal form for exponent
      do 10 i=1,27
   10 g1a(i)=0.
      g1a(7)=-wx/2.d0
      g1a(13)=-wx/2.d0
      g1a(18)=-wy/2.d0
      g1a(22)=-wy/2.d0
      g1a(27)=-wt/2.d0
      endif
c transform exponent to get pseudo hamiltonian
      if (idata.eq.2 .or.idata.eq.3) then
      call inv(ta,tm)
      call fxform(ta,tm,g1a,ga)
c store results in buffer 3
      call mapmap(ga,buf3m,buf3a,buf3m)
      endif
c
c procedure for putting out data
      do 20 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 20
      if (idata.eq.1 .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'exponent for normal form'
      call pcmap(0,j,0,0,g1a,g1m)
      endif
      if (idata.eq.2. .or. idata.eq.3) then
      write(ifile,*)
      write(ifile,*) 'pseudo hamiltonian'
      call pcmap(0,j,0,0,buf3a,buf3m)
      endif
   20 continue
c
c procedure for printing out maps
      do 30 j=1,2
      ifile=0
      if (j.eq.1) then
      if (isend.eq.1 .or. isend.eq.3) ifile=jof
      endif
      if (j.eq.2) then
      if (isend.eq.2 .or. isend.eq.3) ifile=jodf
      endif
      if (ifile.eq.0) goto 30
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normalizing map script A'
      call pcmap(j,j,0,0,buf1a,buf1m)
      endif
      if (ipmaps.eq.2. .or. ipmaps.eq.3) then
      write(ifile,*)
      write(ifile,*) 'normal form script N for transfer map'
      call pcmap(j,j,0,0,buf2a,buf2m)
      endif
   30 continue
c
c procedure for writing out maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      mpo=mpot
      endif
c
c put maps in buffers
c buffers 1 and 2 already contain the transforming map script A
c and the normal form map script N, respectively.
c buffer 3 contains the map which has for its matrix the identity
c matrix  and for its array the pseudo hamiltonian.
c clear buffers 4 and 5
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine submn(p,ha,hm)
c This subroutine computes the norm of a matrix
c Written by Alex Dragt, 10/20/90
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
c
c Calling arrays
      dimension p(6),ha(monoms),hm(6,6)
c
c Local arrays
      dimension tm(6,6)
c
c set up control indices
      iopt=nint(p(1))
      isend=nint(p(2))
c
c copy matrix
      call matmat(hm,tm)
c
c perform calculation
c
c modify matrix if required
      if (iopt .eq. 1) then
      do 10 i=1,6
 10   tm(i,i) = tm(i,i) -1.d0
      endif
c
c compute norm
      call mnorm(tm,ans)
c
c decide where to send and put results
c
      write(6,*) ' matrix norm = ',ans
c
      return
      end
c
***********************************************************************
c
      subroutine tadm(p,fa,fm)
c this is a subroutine for twiss analysis of dynamic maps
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 20 June 1988
c this program will eventually have to be rewritten to improve the
c output data and its format
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms)
      dimension gm(6,6),g1m(6,6),tm(6,6),t1m(6,6)
c
c set up control indices
      idata=nint(p(1))
      ipmaps=nint(p(2))
      isend=nint(p(3))
      iwmaps=nint(p(4))
c
c write headings
      if (isend.eq.1 .or. isend.eq.3) then
      write(jof,*)
      write(jof,*) 'twiss analysis of dynamic map'
      endif
      if (isend.eq.2 .or. isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'twiss analysis of dynamic map'
      endif
c
c beginning of calculation
c remove offensive terms from matrix part of map:
      call dpur2(fa,fm,ga,gm,ta,tm)
c store purifying map script A2 in buffer 1
      call mapmap(ta,tm,buf1a,buf1m)
c
c compute tunes:
      cwx=gm(1,1)
      swx=gm(1,2)
      cwy=gm(3,3)
      swy=gm(3,4)
      cwt=gm(5,5)
      swt=gm(5,6)
      pi=4.*atan(1.d0)
      wx=atan2(swx,cwx)
      if (wx.lt.0.) wx=wx+2.*pi
      wy=atan2(swy,cwy)
      if (wy.lt.0.) wy=wy+2.*pi
      wt=atan2(swt,cwt)
c     if (wt.lt.0.) wt=wt+2.*pi
      tnx=wx/(2.*pi)
      tny=wy/(2.*pi)
      tnt=wt/(2.*pi)
c
c put results in commom/fitdat/ array
      tux=tnx
      tuy=tny
      tus=tnt
c
c remove f3 part of map:
      call dpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,ta,tm)
c resonance decompose purified map:
      call ctodr(g1a,ga)
c
c compute dependence of tune on betatron amplitude:
      hhi=ga(84)
      vvi=ga(85)
      tti=ga(86)
      hvi=ga(87)
      hti=ga(88)
      vti=ga(89)
      hhn=-2.d0*hhi/pi
      vvn=-2.d0*vvi/pi
      ttn=-2.d0*tti/pi
      hvn=-hvi/pi
      htn=-hti/pi
      vtn=-vti/pi
c
c put results in commom/fitdat/ array
      hh=hhn
      vv=vvn
      tt=ttn
      hv=hvn
      ht=htn
      vt=vtn
c
c remove f4 part of map
      call dpur4(g1a,g1m,ga,gm,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,ta,tm,ta,tm)
c
c compute envelopes and twiss parameters
c use the map in buffer 1
      call mapmap(buf1a,buf1m,t1a,t1m)
c compute envelopes
      exh2=t1m(1,1)**2+t1m(1,2)**2
      exv2=t1m(1,3)**2+t1m(1,4)**2
      ext2=t1m(1,5)**2+t1m(1,6)**2
      epxh2=t1m(2,1)**2+t1m(2,2)**2
      epxv2=t1m(2,3)**2+t1m(2,4)**2
      epxt2=t1m(2,5)**2+t1m(2,6)**2
      eyh2=t1m(3,1)**2+t1m(3,2)**2
      eyv2=t1m(3,3)**2+t1m(3,4)**2
      eyt2=t1m(3,5)**2+t1m(3,6)**2
      epyh2=t1m(4,1)**2+t1m(4,2)**2
      epyv2=t1m(4,3)**2+t1m(4,4)**2
      epyt2=t1m(4,5)**2+t1m(4,6)**2
      eth2=t1m(5,1)**2+t1m(5,2)**2
      etv2=t1m(5,3)**2+t1m(5,4)**2
      ett2=t1m(5,5)**2+t1m(5,6)**2
      epth2=t1m(6,1)**2+t1m(6,2)**2
      eptv2=t1m(6,3)**2+t1m(6,4)**2
      eptt2=t1m(6,5)**2+t1m(6,6)**2
      exh=sqrt(exh2)
      exv=sqrt(exv2)
      ext=sqrt(ext2)
      epxh=sqrt(epxh2)
      epxv=sqrt(epxv2)
      epxt=sqrt(epxt2)
      eyh=sqrt(eyh2)
      eyv=sqrt(eyv2)
      eyt=sqrt(eyt2)
      epyh=sqrt(epyh2)
      epyv=sqrt(epyv2)
      epyt=sqrt(epyt2)
      eth=sqrt(eth2)
      etv=sqrt(etv2)
      ett=sqrt(ett2)
      epth=sqrt(epth2)
      eptv=sqrt(eptv2)
      eptt=sqrt(eptt2)
c compute twiss parameters
      call inv(t1a,t1m)
c compute invariants using buffers 2 thru 5
c computation of x invariant
      call clear(buf2a,buf2m)
      buf2a(7)=1.d0
      buf2a(13)=1.d0
      call fxform(t1a,t1m,buf2a,buf3a)
c computation of y invariant
      call clear(buf2a,buf2m)
      buf2a(18)=1.d0
      buf2a(22)=1.d0
      call fxform(t1a,t1m,buf2a,buf4a)
c computation of t invariant
      call clear(buf2a,buf2m)
      buf2a(25)=1.d0
      buf2a(27)=1.d0
      call fxform(t1a,t1m,buf2a,buf5a)
c get twiss parameters from the invariants
c terms for horizontal (x) plane
c 'diagonal' terms
      ax=buf3a(8)/2.d0
      bx=buf3a(13)
      gx=buf3a(7)
c all terms: print all terms later as a twiss invariant
c terms for vertical (y) plane
c 'diagonal' terms
      ay=buf4a(19)/2.d0
      by=buf4a(22)
      gy=buf4a(18)
c all terms: print all terms later as a twiss invariant
c terms for temporal (t) plane
c 'diagonal' terms
      at=buf5a(26)/2.d0
      bt=buf5a(27)
      gt=buf5a(25)
c all terms: print all terms later as a twiss invariant
c
c procedure for writing out tunes and anharmonicities
      if (idata.eq.1 .or. idata.eq.12 .or.
     # idata.eq.13 .or. idata.eq.123) then
      do 10 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 10
      if (iflag.eq.0) goto 10
c write out tunes:
      write (ifile,*)
      write (ifile,*) 'horizontal tune =',tnx
      write (ifile,*) 'vertical tune =',tny
      write (ifile,*) 'temporal tune =',tnt
c write out normalized anharmonicities
      write(ifile,*)
      write(ifile,*) 'normalized anharmonicities'
      write(ifile,*) ' hhn=',hhn
      write(ifile,*) ' vvn=',vvn
      write(ifile,*) ' ttn=',ttn
      write(ifile,*) ' hvn=',hvn
      write(ifile,*) ' htn=',htn
      write(ifile,*) ' vtn=',vtn
   10 continue
      endif
c
c procedure for printing out twiss parameters and envelopes
      if (idata.eq.2 .or. idata.eq.12 .or.
     # idata.eq.23 .or. idata.eq.123) then
      do 20 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 20
      if (iflag.eq.0) goto 20
      write(ifile,*)
      write(ifile,*) 'horizontal twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) ax,bx,gx
      write(ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,buf3a,buf3m)
      write(ifile,*)
      write(ifile,*) 'vertical twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) ay,by,gy
      write(ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,buf4a,buf4m)
      write(ifile,*)
      write(ifile,*) 'temporal twiss parameters'
      write(ifile,*) 'diagonal terms (alpha,beta,gamma)'
      write(ifile,*) at,bt,gt
      write(ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,buf5a,buf5m)
      write(ifile,*)
      write(ifile,*) 'horizontal envelopes (exh,exv,ext;epxh,epxv,epxt)'
      write(ifile,*) exh,exv,ext
      write(ifile,*) epxh,epxv,epxt
      write(ifile,*)
      write(ifile,*) 'vertical envelopes (eyh,eyv,eyt;epyh,epyv,epyt)'
      write(ifile,*) eyh,eyv,eyt
      write(ifile,*) epyh,epyv,epyt
      write(ifile,*)
      write(ifile,*) 'temporal envelopes (eth,etv,ett;epth,eptv,eptt)'
      write(ifile,*) eth,etv,ett
      write(ifile,*) epth,eptv,eptt
   20 continue
      endif
c
c procedure for printing out eigenvectors
      if (idata.eq.3 .or. idata.eq.13 .or.
     # idata.eq.23 .or. idata.eq.123) then
      do 30 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 30
      if (iflag.eq.0) goto 30
c write out the matrix buf1m
      write(ifile,*)
      write(ifile,*) 'matrix of eigenvectors'
      call pcmap (i,0,0,0,buf1a,buf1m)
   30 continue
      endif
c
c procedure for printing of maps
c
c put out script A
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
      do 40 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 40
      if (iflag.eq.0) goto 40
      write(ifile,*)
      write(ifile,*) 'transforming map script A'
      call pcmap(i,i,0,0,ta,tm)
   40 continue
      endif
c
c put out script N
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
      do 50 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (isend.eq.0) goto 50
      if (iflag.eq.0) goto 50
      write(ifile,*)
      write(ifile,*) 'normal form map script N'
      call pcmap(i,i,0,0,ga,gm)
   50 continue
      endif
c
c procedure for writing of maps
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,ta,tm)
      call mapout(0,ga,gm)
      mpo=mpot
      endif
c
c put maps in buffers
c buffer 1 already contains script A2
      call mapmap(ta,tm,buf2a,buf2m)
      call mapmap(ga,gm,buf3a,buf3m)
      call clear(buf4a,buf4m)
      call clear(buf5a,buf5m)
c
      return
      end
c
***********************************************************************
c
      subroutine tasm(p,fa,fm)
c this is a subroutine for twiss analysis of static maps
c Written by Alex Dragt, Spring 1987
c Modified by Alex Dragt, 20 June 1988
c Again modified by Alex Dragt, 19 August 1988
c
      include 'impli.inc'
      include 'param.inc'
c
      include 'parm.inc'
      include 'files.inc'
      include 'buffer.inc'
      include 'fitdat.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),g1a(monoms)
      dimension ta(monoms),t1a(monoms),t2a(monoms)
      dimension gm(6,6),g1m(6,6)
      dimension tm(6,6),t1m(6,6),t2m(6,6)
      dimension am1(6,6),am2(6,6),am3(6,6)
c
c temporary local arrays
      dimension temp1a(monoms),temp2a(monoms),temp3a(monoms)
      dimension temp1m(6,6),temp2m(6,6),temp3m(6,6)
      dimension bm1(6,6),bm2(6,6),bm3(6,6)
      dimension rm1(6,6),rm2(6,6),rm3(6,6)
      dimension em3(6,6)
      dimension pm1(6,6),pm2(6,6),pm3(6,6)
      dimension um1(6,6)
      dimension dm1(6,6)
c
c set up control indices
      iopt=nint(p(1))
      delta=p(2)
      idata=nint(p(3))
      ipmaps=nint(p(4))
      isend=nint(p(5))
      iwmaps=nint(p(6))
c
c write headings
      if (isend.eq.1.or.isend.eq.3) then
      write(jof,*)
      write(jof,*) 'twiss analysis of static map'
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write(jodf,*)
      write(jodf,*) 'twiss analysis of static map'
      endif
c
c first compute closed orbit to get dispersion functions
c this routine does not put out dispersions (subroutine cod does)
c but does put them in common /fitdat/ array for possible fitting
c or plotting
      call fxpt(fa,fm,temp1a,temp1m,temp3a,temp3m)
      dz(1)=temp3m(1,6)
      dz(2)=temp3m(2,6)
      dz(3)=temp3m(3,6)
      dz(4)=temp3m(4,6)
c
c preparatory steps for starting main calculation
c one objective of this calculation is to find script Ac,
c the transforming map with respect to the closed orbit
c remove offensive terms from matrix part of map:
      call clear(buf5a,buf5m)
      call spur2(fa,fm,ga,gm,t1a,t1m,buf5m)
c temporarily save the transforming map associated with sa2 in
c buffer 5 for later use
c
c compute tunes:
      cwx=gm(1,1)
      swx=gm(1,2)
      cwy=gm(3,3)
      swy=gm(3,4)
      pi=4.*atan(1.d0)
      wx=atan2(swx,cwx)
      if (wx.lt.0.) wx=wx+2.*pi
      wy=atan2(swy,cwy)
      if (wy.lt.0.) wy=wy+2.*pi
      tx=wx/(2.*pi)
      ty=wy/(2.*pi)
c put results in commom/fitdat/ array
      tux=tx
      tuy=ty
c
c preparatory steps for continuing calculation
c remove offensive chromatic terms from f3 part of map:
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map
      call concat(t1a,t1m,buf5a,buf5m,ta,tm)
c resonance decompose purified map:
      call ctosr(g1a,t2a)
c
c compute chromaticities
c
c procedure when IOPT = 1
      if (iopt.eq.1) then
c compute first order chromaticities:
      chrox1=-(1.d0/pi)*t2a(28)
      chroy1=-(1.d0/pi)*t2a(29)
c compute second order chromaticities:
      chrox2=-(2.d0/pi)*t2a(84)
      chroy2=-(2.d0/pi)*t2a(85)
c put results in commom/fitdat/ array
      cx=chrox1
      cy=chroy1
      qx=chrox2
      qy=chroy2
c
c compute tunes about closed orbit
      delta2=delta*delta
      txc=tx+delta*chrox1+delta2*chrox2
      tyc=ty+delta*chroy1+delta2*chroy2
      tsc=txc-tyc
c
      endif
c
c procedure when IOPT = 2
      if (iopt.eq.2) then
c compute first order chromaticities:
      chrox1=(beta/pi)*t2a(28)
      chroy1=(beta/pi)*t2a(29)
c compute second order chromaticities:
      beta2=beta*beta
      beta3=beta*beta2
      chrox2=(t2a(28)*(beta-beta3)-2.d0*t2a(84)*beta2)/pi
      chroy2=(t2a(29)*(beta-beta3)-2.d0*t2a(85)*beta2)/pi
c put results in commom/fitdat/ array
      cx=chrox1
      cy=chroy1
      qx=chrox2
      qy=chroy2
c
c compute tune about closed orbit
      delta2=delta*delta
      txc=tx+delta*chrox1+delta2*chrox2
      tyc=ty+delta*chroy1+delta2*chroy2
      tsc=txc-tyc
c
      endif
c
c write out tunes and chromaticities
      if (idata.eq.1 .or. idata.eq.12 .or.
     # idata.eq.13 .or. idata.eq.123) then
      if(isend.eq.0) goto 11
      do 10 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 10
      write (ifile,*)
      if(iopt.eq.1) then
      write (ifile,*) 'tunes and chromaticities for delta defined in',
     #' terms of P sub tau:'
      endif
      if(iopt.eq.2) then
      write (ifile,*) 'tunes and chromaticities for delta defined in',
     #' terms of momentum deviation:'
      endif
      write (ifile,*)
      write (ifile,*) 'horizontal tune =',tx
      write (ifile,*) 'first order horizontal chromaticity =',chrox1
      write (ifile,*) 'second order horizontal chromaticity =',chrox2
      write (ifile,*) 'horizontal tune when delta =',delta
      write (ifile,*) txc
      write (ifile,*)
      write (ifile,*) 'vertical tune =',ty
      write (ifile,*) 'first order vertical chromaticity =',chroy1
      write (ifile,*) 'second order vertical chromaticity =',chroy2
      write (ifile,*) 'vertical tune when delta =',delta
      write (ifile,*) tyc
      write (ifile,*)
      write (ifile,*) 'tune separation when delta=',delta
      write (ifile,*) tsc
   10 continue
   11 continue
      endif
c
c preparatory steps for continuing calculation
c remove offensive geometric terms from f3 part of map:
      call sgpur3(g1a,g1m,ga,gm,t2a,t2m)
c accumulate transforming map
      call concat(t2a,t2m,ta,tm,t1a,t1m)
c resonance decompose purified map:
      call ctosr(ga,t2a)
c
c compute dependence of tune on betatron amplitude:
      hhi=t2a(87)
      vvi=t2a(88)
      hvi=t2a(89)
      hhn=-2.d0*hhi/pi
      vvn=-2.d0*vvi/pi
      hvn=-hvi/pi
c put results in commom/fitdat/ array
      hh=hhn
      vv=vvn
      hv=hvn
c
c write out normalized anharmonicities
      if (idata.eq.1 .or. idata.eq.12 .or.
     # idata.eq.13 .or. idata.eq.123) then
      if(isend.eq.0) goto 21
      do 20 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 20
      write(ifile,*)
      write(ifile,*) 'normalized anharmonicities'
      write(ifile,*) ' hhn=',hhn
      write(ifile,*) ' vvn=',vvn
      write(ifile,*) ' hvn=',hvn
   20 continue
   21 continue
      endif
c
c complete computation of script Ac and script N
c store script Ac in buffer 1 and script N in buffer 2
      call spur4(ga,gm,buf2a,buf2m,t2a,t2m)
c accumulate transforming map to get script Ac
      call concat(t2a,t2m,t1a,t1m,buf1a,buf1m)
c
c compute twiss parameter expansions (invariants) and envelopes
c
c preparatory steps for continuing calculation
c compute fixed point and map around it, and put the transforming
c map to the closed orbit in buffer 3
      call fxpt(fa,fm,g1a,g1m,buf3a,buf3m)
c extract the betatron factor of the map and store it in buffer 4
      call betmap(g1a,g1m,buf4a,buf4m)
c find the transforming (conjugating) map script Ab for the betatron factor
c use the map in buffer 5 to purify the f2 part of the betatron factor
      call sndwch(buf5a,buf5m,buf4a,buf4m,g1a,g1m)
c remove offensive chromatic terms from f3 part of betatron factor
      call scpur3(g1a,g1m,ga,gm,ta,tm)
c accumulate transforming map:
      call concat(ta,tm,buf5a,buf5m,t2a,t2m)
c remove offensive terms from f4 part of betatron factor
      call spur4(ga,gm,g1a,g1m,t1a,t1m)
c accumulate transforming map to get script Ab:
      call concat(t1a,t1m,t2a,t2m,ta,tm)
c store script Ab in buffer 5
      call mapmap(ta,tm,buf5a,buf5m)
c invert script Ab
      call inv(ta,tm)
c
c compute invariants
c computation of x invariant
      call clear(ga,gm)
      ga(7)=1.d0
      ga(13)=1.d0
      call fxform(ta,tm,ga,t1a)
c computation of y invariant
      call clear(ga,gm)
      ga(18)=1.d0
      ga(22)=1.d0
      call fxform(ta,tm,ga,t2a)
c
c preliminary calculations required for envelopes and eigenvalues
      if (idata.eq.2 .or. idata.eq.3 .or.
     # idata.eq.12 .or. idata.eq.13 .or.
     # idata.eq.23 .or. idata.eq.123) then
c put script Ab in ta,tm
      call mapmap(buf5a,buf5m,ta,tm)
c make chromatic expansion of ta,tm
c the result will be used to compute both envelopes and eigenvectors
      call chrexp(iopt,delta,ta,tm,am1,am2,am3)
      endif
c
c see if output of twiss parameters and envelopes is desired
      if (idata.eq.2 .or. idata.eq.12 .or.
     # idata.eq.23 .or. idata.eq.123) then
c
c continue with calculation
c
c terms for horizontal (x) plane
c 'diagonal terms'
      ax0=t1a(8)/2.d0
      bx0=t1a(13)
      gx0=t1a(7)
c all terms: later print t1a as a map
c terms for vertical (y) plane
c 'diagonal'terms
      ay0=t2a(19)/2.d0
      by0=t2a(22)
      gy0=t2a(18)
c all terms: later print t2a as a map
c
c put horizontal and vertical results in commom/fitdat/ array
      ax=ax0
      bx=bx0
      gx=gx0
      ay=ay0
      by=by0
      gy=gy0
c
c compute envelopes
      exhc2=am3(1,1)**2+am3(1,2)**2
      exvc2=am3(1,3)**2+am3(1,4)**2
      epxhc2=am3(2,1)**2+am3(2,2)**2
      epxvc2=am3(2,3)**2+am3(2,4)**2
      eyhc2=am3(3,1)**2+am3(3,2)**2
      eyvc2=am3(3,3)**2+am3(3,4)**2
      epyhc2=am3(4,1)**2+am3(4,2)**2
      epyvc2=am3(4,3)**2+am3(4,4)**2
      exhc=sqrt(exhc2)
      exvc=sqrt(exvc2)
      epxhc=sqrt(epxhc2)
      epxvc=sqrt(epxvc2)
      eyhc=sqrt(eyhc2)
      eyvc=sqrt(eyvc2)
      epyhc=sqrt(epyhc2)
      epyvc=sqrt(epyvc2)
c
c write out twiss functions invariants, and envelopes
      if(isend.eq.0) goto 31
      do 30 i=1,2
      if (i.eq.1) then
      ifile=jof
      iflag=1
      if (isend.eq.2) iflag=0
      endif
      if (i.eq.2) then
      ifile=jodf
      iflag=1
      if (isend.eq.1) iflag=0
      endif
      if (iflag.eq.0) goto 30
      write (ifile,*)
      write (ifile,*) 'twiss parameters, invariants, and envelopes'
c
c write twiss functions and invariants
      write (ifile,*)
      write (ifile,*) 'horizontal parameters'
      write (ifile,*) 'on energy diagonal terms (alpha,beta,gamma)'
      write (ifile,*) ax0,bx0,gx0
      write (ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,t1a,t1m)
      write (ifile,*)
      write (ifile,*) 'vertical parameters'
      write (ifile,*) 'on energy diagonal terms (alpha,beta,gamma)'
      write (ifile,*) ay0,by0,gy0
      write (ifile,*) 'full twiss invariant written as a map'
      call pcmap(0,i,0,0,t2a,t2m)
c
c write out envelopes
      if(iopt.eq.1) then
      write (ifile,*)
      write (ifile,*) 'envelopes for delta defined in terms of',
     #' P sub tau with delta =',delta
      endif
      if(iopt.eq.2) then
      write (ifile,*) 'envelopes for delta defined in terms of',
     #' momentum deviation with delta =',delta
      endif
      write (ifile,*)
      write (ifile,*) 'normalized horizontal envelope coefficients'
     #,' (exhc,exvc;epxhc,epxvc)'
      write (ifile,*) exhc,exvc
      write (ifile,*) epxhc,epxvc
      write (ifile,*)
      write (ifile,*) 'normalized vertical envelope coefficients'
     #,' (eyhc,eyvc;epyhc,epyvc)'
      write (ifile,*) eyhc,eyvc
      write (ifile,*) epyhc,epyvc
c
   30 continue
   31 continue
      endif
c
c Procedure for output of eigenvectors.
      if (idata.eq.3 .or. idata.eq.13 .or.
     # idata.eq.23 .or. idata.eq.123) then
c Print out matrices tm, am1, and am2.
      if(isend.eq.0) goto 41
      do 40 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 40
c
c procedure when IOPT = 1
      if (iopt.eq.1) then
      write(ifile,198)
  198 format(//,1x,'eigenvector expansion for delta defined',
     #1x,'in terms of P sub tau:')
      write(ifile,200)
  200 format(/,1x,'on energy matrix of eigenvectors')
      endif
c
c procedure when IOPT = 2
      if (iopt.eq.2) then
      write(ifile,199)
  199 format(//,1x,'eigenvector expansion for delta defined',
     #1x,'in terms of momentum deviation:')
      write(ifile,201)
  201 format(/,1x,'on momentum matrix of eigenvectors')
      endif
c
      call pcmap(i,0,0,0,ta,tm)
      write(ifile,300)
  300 format(//,1x,'delta correction')
      call pcmap(i,0,0,0,ta,am1)
      write(ifile,400)
  400 format(//,1x,'delta**2 correction')
      call pcmap(i,0,0,0,ta,am2)
c Print out value of twiss matrix
      write(ifile,402) delta
  402 format(//,1x,'matrix of eigenvectors when delta= ',d15.8)
      call pcmap(i,0,0,0,ta,am3)
c
c test results
c
      write(ifile,*) 'test results'
**************************************************************************
***************************************************************************
c compute matrix for betatron portion of map
      call chrexp(iopt,delta,buf4a,buf4m,bm1,bm2,bm3)
c compute tune matrix
      call spur2(fa,fm,ga,gm,t1a,t1m,t2m)
      call scpur3(ga,gm,g1a,g1m,t1a,t1m)
      call clear(temp1a,temp1m)
      call clear(temp3a,temp3m)
      call matmat(gm,temp1m)
      call ctosr(g1a,temp2a)
      temp3a(28)=temp2a(28)
      temp3a(29)=temp2a(29)
      temp3a(84)=temp2a(84)
      temp3a(85)=temp2a(85)
      call srtoc(temp3a,temp1a)
      call chrexp(iopt,delta,temp1a,temp1m,rm1,rm2,rm3)
c set up matrix of eigenvectors
      call matmat(am3,em3)
c form the product pm1=em3*rm3
      call mmult(em3,rm3,pm1)
c invert pm1
      call inv(t1a,pm1)
c form the product pm2=bm3*em3
      call mmult(bm3,em3,pm2)
c form the product pm3=(pm1 inverse)*pm2
      call mmult(pm1,pm2,pm3)
c form the negative identity matrix
      call ident(temp1a,um1)
      call smmult(-1.d0,um1,um1)
c form the difference dm1=pm3-um1
      call madd(pm3,um1,dm1)
c print the results bm3 and dm1
      call pcmap(i,0,0,0,ta,bm3)
      call pcmap(i,0,0,0,ta,dm1)
**************************************************************
**************************************************************
c
   40 continue
   41 continue
      endif
c
c Procedure for printing of maps.
c
      if (ipmaps.eq.1 .or. ipmaps.eq.3) then
c print out script Ac and script N
      if(isend.eq.0) goto 51
      do 50 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 50
      write(ifile,600)
  600 format(//,1x,'transforming map with respect to
     # the closed orbit')
      call pcmap(i,i,0,0,buf1a,buf1m)
      write(ifile,700)
  700 format(//,1x,'normal form for transfer map')
      call pcmap(i,i,0,0,buf2a,buf2m)
   50 continue
   51 continue
      endif
c
      if (ipmaps.eq.2 .or. ipmaps.eq.3) then
c print out betatron portion of map and script Ab
      if(isend.eq.0) goto 61
      do 60 i=1,2
      if (i.eq.1) then
      iflag=1
      if (isend.eq.2) iflag=0
      ifile=jof
      endif
      if (i.eq.2) then
      iflag=1
      if (isend.eq.1) iflag=0
      ifile=jodf
      endif
      if (iflag.eq.0) goto 60
      write(ifile,610)
  610 format(//,1x,'betatron factor of transfer map')
      call pcmap(i,i,0,0,buf4a,buf4m)
      write(ifile,710)
  710 format(//,1x,'transforming map for betatron factor')
      call pcmap(i,i,0,0,buf5a,buf5m)
   60 continue
   61 continue
      endif
c
c Procedure for writing of maps.
      if (iwmaps.gt.0) then
      mpot=mpo
      mpo=iwmaps
      call mapout(0,buf1a,buf1m)
      call mapout(0,buf2a,buf2m)
      call mapout(0,buf3a,buf3m)
      call mapout(0,buf4a,buf4m)
      call mapout(0,buf5a,buf5m)
      mpo=mpot
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine tbas(p,fa,fm)
c this routine translates bases
c Written by Alex Dragt, Spring 1987
c
      include 'impli.inc'
      include 'param.inc'
c
c Calling arrays
      dimension p(6),fa(monoms),fm(6,6)
c
c Local arrays
      dimension ga(monoms),gm(6,6)
c
      call mapmap(fa,fm,ga,gm)
      iopt=nint(p(1))
      if (iopt.eq.1) call ctosr(ga,fa)
      if (iopt.eq.2) call ctodr(ga,fa)
      if (iopt.eq.3) call srtoc(ga,fa)
      if (iopt.eq.4) call drtoc(ga,fa)
      return
      end
c
*******************************************************************
c
      subroutine trda(p,fa,fm)
c  this subroutine transports a dynamic script A
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'trda not yet available'
      return
      end
 
c
c
*******************************************************************
c
      subroutine trsa(p,fa,fm)
c  this subroutine transports a static script A
      include 'impli.inc'
      include 'param.inc'
      dimension p(6),fa(monoms),fm(6,6)
c
      write(6,*) 'trsa not yet available'
      return
      end
