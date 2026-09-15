c=====================================================
      subroutine user1(p)
c   Modified to plot beam envelopes  18 Dec 90
c   Overhauled some more to take various input beams
c          CTM   Oct 91
c
c          Ray tracing added  Mar 96  FRN (&CTM)
c          FOV calculation added  Dec 97 (CTM)
c          both moved to other USER routines May 01
c
c  This routine applies the current map in (fa,fm) to the initial sigma
c  matrix given in the parameter set p = (s11,s12,s22,s33,s34,s44), to
c  produce the output sigma matrix of the same form, which is saved in
c
c    control parameters
c    p(1) = job = -N, 0, +N :    -N (N<0) means silent running in
c                  a fit or scan loop, with no output files.
c           job = 0 (default) is to fill in /linbuf/
c                 N = 1 means also slice the elements and fill in the
c                       sincos buffer for future use in FOV and RAY.
c                 N = 2 means to compute the envelope from the initial
c                       moments read in by USER8. The resulting envelope
c                       excursions are stored in UCALC(33 - 36).
c                 N = 3 means to also generate the table of envelope
c                       excursions in each element. ISEND controls where
c                       (and if) it is written.
c    p(2) = menv = 0 for no output files. (used in fit loops, etc.)
c                = 1 to write the sincos buffer to the .MSC file (lun=34
c                = 2 to write the beam envelope to the .ENV file (lun=24
c                = 3 to write both.
c    p(3) = mincut = minimum no. of slices/element
c    p(4) = step: nominal thickness of one slice. nslice=L/step
c    p(5) = isend as usual: controls jof and jodf output.
c    p(6) = jtran = 0 for no translation files.
c                 = 1 to write ADLIB dump on unit 22
c                 = 2 to write TRACE3D format to .TRC file on unit 36
c                 = 3 to write old TRANSPORT format to .TRN on unit 38
c                 = 4 to write MAD format to .MAD file on unit 40
c
c    beam parameters: xx(nn) = xmax = sqrt(betax*xemit)
c                ax(nn) = alphax
c                px(nn) = thetax = sqrt(gammax*xemit)
c                yy(nn) = ymax = sqrt(betay*yemit)
c                ay(nn) = alphay
c                py(nn) = thetay = sqrt(gammay*yemit)
c     Tom Mottershead  LANL  31-Mar-89
c             revised Dec 97
c---------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'linbuf.inc'
      include 'sigbuf.inc'
      include 'usrdat.inc'
      include 'parset.inc'
      include 'sincos.inc'
      include 'taylor.inc'
      save
      dimension xc(2),pxc(2),xs(2),pxs(2),yc(2),pyc(2),ys(2),pys(2)
      dimension p(6), fm(6,6)
      logical  lterm, lfile, lenv, jenv, msc, nosig
c
c       decode input parameters
c
      job  = nint(p(1))
      menv = nint(p(2))
      mincut = nint(p(3))
      step = p(4)
      isend = nint(p(5))
      jtran = nint(p(6))
c
c       backward compatibility resets for no output if job < 0
c
      msc = .false.
      lenv = .false.
      jenv = .false.
      nosig = .true.
      if(job.lt.0) then
         iquiet=1
         job = -job
      endif
      if(job.gt.1) nosig = .false.
      if(job.eq.3) jenv = .true.
      if((menv.eq.1).or.(menv.eq.3)) msc = .true.
      if(menv.eq.2) lenv = .true.
      lterm = .false.
      lfile = .false.
      if((isend.eq.1).or.(isend.eq.3)) lterm = .true.
      if(isend.gt.1) lfile = .true.
c
c       temp setup of dump files
c
      ltrace = 1
      ltran = 0
      if(jtran.eq.1) ltran = 22
      if(jtran.eq.2) ltran = 36
      zero = 0.0d0
      one  = 1.0d0
      two  = 2.0d0
c
c      fillin linbuf from the loop arrays
c
      npsu = 9
      call filin
c
c     count the beamline
c
      ndrft = 0
      nquad = 0
      do 10 nn = nbgn,nend
      if(lintyp(nn).eq.0) then
         ndrft = ndrft+1
      else
         nquad = nquad+1
      endif
  10  continue
      if(lterm) write(jof,13) ndrft,nquad
      if(lfile) write(jodf,13) ndrft,nquad
  13  format('  Line has',i4,' drifts and',i4,' quads')
      if(job.eq.0) return
      nep = 1 + ndrft +3*nquad
c        type *, ' ACTPAR on FILIN exit:'
c        type *, aapar,bbpar,ccpar,ddpar,eepar,ffpar
      if(lterm) write(jof,*) ' Beam Line Summary'
      if(lfile) write(jodf,*) ' Beam Line Summary'
      total = zero
      gltot = zero
      zbeg = zero
c
c        write beamline list and element blocks to .ENV file
c
      ap = 0.0
      nn = 0
      if(lterm) write(jof,16)
      if(lfile) write(jodf,16)
  16  format(3x,'n',2x,'name',4x,' codes',2x,'aper(cm) length(cm)',4x,
     * 'Gradient',5x,'PoleTip(g)',3x,'path(cm)')
      do 20 nn = nbgn, nend
      ltyp = lintyp(nn)
      zbeg = total
      total = total + elong(nn)
      cmtot = 100.0*total
      cml = 100.0*elong(nn)
      glp = elong(nn)*strong(nn)
      gltot = gltot + abs(glp)
      rcm = 100.0*radius(nn)
      qupole = 1.e4*strong(nn)*radius(nn)
      sxpole = 1.e4*pssext(nn)*radius(nn)**2
      ocpole = 1.e4*psoctu(nn)*radius(nn)**3
      if(lterm) then
         write(jof,17) nn, cname(nn), mltyp1(nn), mltyp2(nn),
     *   lintyp(nn), rcm, cml, strong(nn), qupole, cmtot
         if(pssext(nn).ne.0.0) write(jof,18) 'sext',pssext(nn),sxpole
         if(psoctu(nn).ne.0.0) write(jof,18) 'oct ',psoctu(nn),ocpole
      endif
      if(lfile) then
         write(jodf,17) nn, cname(nn), mltyp1(nn), mltyp2(nn),
     *   lintyp(nn), rcm, cml, strong(nn), qupole, cmtot
         if(pssext(nn).ne.0.0) write(jodf,18) 'sext',pssext(nn),sxpole
         if(psoctu(nn).ne.0.0) write(jodf,18) 'oct ',psoctu(nn),ocpole
      endif
  17  format(i4,2x,a,3i2,f8.2,f10.3,f15.7,f15.5,f10.3)
  18  format(16x,a4,18x,f15.7,f15.5)
  20  continue
      ucalc(39) = total
      ucalc(40) = gltot
      nask = total/step
      if(lterm) write(jof,23) total,nask,step,gltot
      if(lfile) write(jodf,23) total,nask,step,gltot
  23  format(' total length = ',f15.5,' in',i5,' steps of',f9.4,
     * f12.6,'= Sum|GL|')
c
c       set chromatic correlation coefficients
c
      mu = 3
      call setcor(mu,uwx,uwy)
c
c       write blocks on top of .ENV file
c
      lun = 24
      if(lenv) call blocks(lun,nep,uwx,uwy)
c
c       initialize sin/cos variables
c
      jp = 1
      ni = 1
      zu = zero
      xc(ni)  = one
      pxc(ni) = zero
      xs(ni)  = zero
      pxs(ni) = one
      yc(ni)  = one
      pyc(ni) = zero
      ys(ni)  = zero
      pys(ni) = one
c
c     initialize sin,cosine file
c
      lun = 34
      if(msc) call blocks(lun,nep,uwx,uwy)
      rcm = 100.0*radius(1)
      if(msc) write(lun,147) zero,rcm,one,zero,one,zero
c
c     initialize element loop
c
      ni = 0
      nf = 0
      kka = 0
      kkb = 1
      xmax = 0.0
      ymax = 0.0
      xmin = 1.0e38
      ymin = 1.0e38
      if(jenv) then
         if(lterm) write(jof,88)
         if(lfile) write(jodf,88)
      endif
  88  format('  element  cuts',4x,'zxmax',6x,'xmax',6x,'xmin',5x,
     *'zymax',6x,'ymax',6x,'ymin')
c
c     initialize envelope points
c
      xx(1) = xxin
      ax(1) = axin
      px(1) = pxin
      yy(1) = yyin
      ay(1) = ayin
      py(1) = pyin
      if(lenv) write(24,27) zu,xx(1),px(1),yy(1),py(1),ax(1),ay(1)
  27  format(f10.4,4(1pe11.4),2(1pe12.4))
c
c     big loop over all elements
c
      do 200 nn = nbgn, nend
c
c       setup transfer matrix for one slice
c
      nslice = elong(nn)/step
      if(nslice.lt.mincut) nslice = mincut
      delta = elong(nn)/float(nslice)
      if(lterm) write(jof,193) nn, elong(nn), nslice, delta
      if(lfile) write(jodf,193) nn, elong(nn), nslice, delta
 193  format('Element',i4,' Length=',f10.5,i4,' slices of dL=',f10.4)
      grad = strong(nn)
      rcm = 100.0*radius(nn)
      call rmquad(delta,grad,fm)
c-------------------   x plane  ------------
      cfx = fm(1,1)
      sfx = fm(1,2)
      cpx = fm(2,1)
      spx = fm(2,2)
      txx = cfx**2
      tax = 2.0*cfx*sfx
      tpx = sfx**2
      txa = cfx*cpx
      taa = cfx*spx+sfx*cpx
      tpa = sfx*spx
      txp = cpx**2
      tap = 2.0*cpx*spx
      tpp = spx**2
c-------------------   y plane  ------------
      cfy = fm(3,3)
      sfy = fm(3,4)
      cpy = fm(4,3)
      spy = fm(4,4)
      ryy = cfy**2
      ray = 2.0*cfy*sfy
      rpy = sfy**2
      rya = cfy*cpy
      raa = cfy*spy+sfy*cpy
      rpa = sfy*spy
      ryp = cpy**2
      rap = 2.0*cpy*spy
      rpp = spy**2
c
c         propagate through all the slices
c
      exmax = 0.0
      eymax = 0.0
      exmin = 1.0e30
      eymin = 1.0e30
      zbgn = zu
      do 100 kk = 1, nslice
        ni = kka + 1
        nf = kkb + 1
        kka = 1- kka
        kkb = 1- kkb
        zu = zu + delta
c
c    collect Sinelike and Cosinelike rays for FOV calculation
c
      xc(nf)  = cfx*xc(ni)+sfx*pxc(ni)
      pxc(nf) = cpx*xc(ni)+spx*pxc(ni)
      xs(nf)  = cfx*xs(ni)+sfx*pxs(ni)
      pxs(nf) = cpx*xs(ni)+spx*pxs(ni)
      yc(nf)  = cfy*yc(ni)+sfy*pyc(ni)
      pyc(nf) = cpy*yc(ni)+spy*pyc(ni)
      ys(nf)  = cfy*ys(ni)+sfy*pys(ni)
      pys(nf) = cpy*ys(ni)+spy*pys(ni)
      cx(jp) = xc(nf)
      sx(jp) = xs(nf)
      cy(jp) = yc(nf)
      sy(jp) = ys(nf)
      za(jp) = zu
      if(msc) write(lun,147) za(jp),rcm,cx(jp),sx(jp),cy(jp),sy(jp)
  147 format(6f13.7)
      jp = jp+1
      if(nosig) go to 100
c
c      compute x-plane final beam ellipse
c
      xss = xx(ni)**2
      xa = -xx(ni)*ax(ni)*px(ni)/sqrt(1.0d0 + ax(ni)**2)
      xps = px(ni)**2
c      sigf1 = (cfx**2)*xss + 2.0*cfx*sfx*xa + (sfx**2)*xps
c      sigf2 = cfx*cpx*xss + (cfx*spx+sfx*cpx)*xa + sfx*spx*xps
c      sigf3 = (cpx**2)*xss + 2.0*cpx*spx*xa + (spx**2)*xps
      sigf1 = txx*xss + tax*xa + tpx*xps
      sigf2 = txa*xss + taa*xa + tpa*xps
      sigf3 = txp*xss + tap*xa + tpp*xps
      exsq = sigf1*sigf3 - sigf2**2
      if(exsq.le.0.0) exsq = 1.e-30
      ex = sqrt(exsq)
      xx(nf) = sqrt(sigf1)
      if(xx(nf).gt.xmax) then
         xmax = xx(nf)
         zxmax = zu
      endif
      if(xx(nf).gt.exmax) then
         exmax = xx(nf)
         zxm = zu
      endif
      if(xx(nf).lt.exmin) exmin = xx(nf)
      px(nf) = sqrt(sigf3)
      ax(nf) = -sigf2/ex
c
c     compute y-plane beam ellipse
c
      yss  = yy(ni)**2
      ya  = -yy(ni)*ay(ni)*py(ni)/sqrt(1.0d0 + ay(ni)**2)
      yps = py(ni)**2
c      sigf4 = (cfy**2)*yss + 2.0*cfy*sfy*ya + (sfy**2)*yps
c      sigf5 = cfy*cpy*yss + (cfy*spy+sfy*cpy)*ya + sfy*spy*yps
c      sigf6 = (cpy**2)*yss + 2.0*cpy*spy*ya + (spy**2)*yps
      sigf4 = ryy*yss + ray*ya + rpy*yps
      sigf5 = rya*yss + raa*ya + rpa*yps
      sigf6 = ryp*yss + rap*ya + rpp*yps
      eysq = sigf4*sigf6 - sigf5**2
      if(eysq.le.0.0) eysq = 1.e-30
      ey = sqrt(eysq)
      yy(nf) = sqrt(sigf4)
      if(yy(nf).gt.ymax) then
         ymax = yy(nf)
         zymax = zu
      endif
      if(yy(nf).gt.eymax) then
         eymax = yy(nf)
         zym = zu
      endif
      if(yy(nf).lt.eymin) eymin = yy(nf)
      py(nf) = sqrt(sigf6)
      ay(nf) = -sigf5/ey
      if(lenv) write(24,27) zu,xx(nf),px(nf),yy(nf),py(nf),ax(nf)
     * ,ay(nf)
 100  continue
      if(nosig) go to 200
      exmax = 100.0*exmax
      eymax = 100.0*eymax
      exmin = 100.0*exmin
      eymin = 100.0*eymin
      dzz = 100.0*delta
      zxm = 100.0*(zxm-zbgn)
      zym = 100.0*(zym-zbgn)
      if(jenv) then
      if(lterm) write(jof,117) cname(nn), nslice, zxm, exmax, exmin,
     * zym, eymax, eymin
      if(lfile) write(jodf,117) cname(nn), nslice, zxm, exmax, exmin,
     * zym, eymax, eymin
      endif
  117 format(2x,a,i4,6f10.3)
 200  continue
      if(nosig) go to 400
      if(jenv) then
         if(lterm) write(jof,203), nf, xmax, zxmax, ymax, zymax
         if(lfile) write(jodf,203), nf, xmax, zxmax, ymax, zymax
      endif
 203  format(i6,' points generated.  xmax =',f10.6,' at z =',f10.6,/
     * 26x,'ymax =',f10.6,' at z =',f10.6)
c
c      save envelope excursions and elliptic parameters in ucalc
c
      big = xmax
      if(ymax.gt.big) big = ymax
      ucalc(33) = xmax
      ucalc(34) = ymax
      ucalc(35) = xmax - ymax
      ucalc(36) = big
      ucalc(46) = sqrt(ymax/xmax)
      ucalc(47) = sqrt(xmax*ymax)
      ucalc(48) = ucalc(40)*(ucalc(47)**3)
  400 continue
c
c              optional beamline translation files
c
 600  if(jtran.eq.1) call adump(ltran)
      if(jtran.eq.2) call trcdmp(ltran,ltrace)
      if(jtran.eq.3) call trnspt(ltran)
      if(jtran.eq.4) call madump(ltran)
 777  return
      end
c=====================================================
      subroutine filin
c
c       filin's job is to fill in /linbuf/ from the contents of
c       a loop. Derived from wcl, written by Alex Dragt, 23 August 1988
c       Based on the SRs cqlate and pmif. WCL modified by F. Ner
c       Dec 90 to fill LINBUF. That function of wcl save is this routine
c       April 91, by T. Mottershead.
c
      include 'impli.inc'
      include 'param.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
c
c common blocks
c
      include 'elmnts.inc'
      include 'items.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'parset.inc'
      include 'files.inc'
      include 'loop.inc'
      include 'core.inc'
c
      dimension ltype(50)
      save
      data (ltype(j),j=1,25)/
     #   0,2,2,2,2,
     #   2,3,2,1,0,
     #   0,0,4,5,7,
     #  -1,7,1,7,6,
     #   7,7,7,7,7/
c
c    set parameter set to use (old=9, now pass in from caller CTM 5/00)
c
c     npsu = 9
c
c  start routine by checking to see if a loop exists
c
      if(nloop.le.0) then
         write(jodf,*)  '  error in filin: no loop in labor'
         write(jof,*)  '  error in filin: no loop in labor'
         return
      endif
c
c  initialize linbuf pointers:
c
      nbgn = 1
      nend = 0
      npst = 0
c
c  contents of biglist
c
c      write(jodf,*) joy,' = joy in filin'
      do 100 jj=1,joy
      mmv = mim(jj)
c       element
        if(mmv.lt.0) then
          ktyp = 1
          mnu = -mmv
c       user supplied element
        else if(mmv.gt.5000) then
          ktyp = 1
          mnu = mmv-5000
c       lump
        else
          ktyp = 3
        endif
c       write(jodf,*) '  filin: ',jj,'=jj',mmv,'=mim(jj)',mnu,'=mnu',
c     * ktyp,'=ktyp'
c procedure for a menu item
      if(ktyp.ne.1) go to 100
         imax=nrp(nt1(mnu),nt2(mnu))
         if(imax.eq.0)goto 100
c
c   reset default radius from pset(npsu) (T. Mottershead 23 Jan 91)
c   keep all values from the same pset for other purposes, e.g. FOV calc
c   (CTM 12/97)
c   CTM 11/99: save all instances of npsu in linbuf for future use
c
      if((nt1(mnu).eq.3).and.(nt2(mnu).eq.npsu)) then
         aapar = pmenu(1+mpp(mnu))
         bbpar = pmenu(2+mpp(mnu))
         ccpar = pmenu(3+mpp(mnu))
         ddpar = pmenu(4+mpp(mnu))
         eepar = pmenu(5+mpp(mnu))
         ffpar = pmenu(6+mpp(mnu))
         npst = npst + 1
         aap(npst) = aapar
         bbp(npst) = bbpar
         ccp(npst) = ccpar
         ddp(npst) = ddpar
         eep(npst) = eepar
         ffp(npst) = ffpar
      endif
c
c  Put simple element in linbuf (F. Neri 12/18/1990)
c
      if (nt1(mnu).eq.1) then
        lt = ltype(nt2(mnu))
        if(lt.eq.1.or.lt.eq.0) then
          nend = nend+1
          cname(nend) = lmnlbl(mnu)
          lintyp(nend) = lt
          mltyp1(nend) = nt1(mnu)
          mltyp2(nend) = nt2(mnu)
          elong(nend) = pmenu(1+mpp(mnu))
          radius(nend) = aapar
          strong(nend) = 0.d0
c
c        normal quad
c
          if(nt2(mnu).eq.9) strong(nend) = pmenu(2+mpp(mnu))
c
c   quick & dirty fix for new cfqd (think again about when psets are exe
c
          if(nt2(mnu).eq.18) then
               ids = nint(pmenu(2+mpp(mnu)))
               strong(nend) = pst(1,ids)
               pssext(nend) = pst(3,ids)
               psoctu(nend) = pst(5,ids)
          endif
        endif
      endif
 100  continue
      return
      end
c=====================================================
      subroutine setcor(kor,wx,wy)
      include 'impli.inc'
      include 'param.inc'
      include 'taylor.inc'
      include 'map.inc'
      include 'parset.inc'
      include 'usrdat.inc'
      save
      zero = 0.0d0
      two = 2.0d0
c
c     No corr: Straight in
c
      if(kor.eq.0) then
         wx = zero
         wy = zero
      endif
c
c     Fourier corr
c
      if(kor.eq.1) then
         wx = -tmh(1,1)/tmh(1,2)
         wy = -tmh(3,3)/tmh(3,4)
      endif
c
c     Chromatic corr for identity lens
c
      if(kor.eq.2) then
         wx = -th(38)/(two*th(53))
         wy = -th(70)/(two*th(76))
      endif
c
c     Chromatic corr in general from Taylor map
c
      if(kor.eq.3) then
         call tugen(th,tmh)
         wx = -tumat(12,1)/tumat(17,1)
         wy = -tumat(21,3)/tumat(24,3)
      endif
c
c     general correlation from ucalc(N) and N+1
c
      if(kor.lt.0) then
        nc = -kor
        if(nc.gt.250) return
        wx = ucalc(nc)
        wy = ucalc(nc+1)
      endif
      return
      end
c=====================================================
      subroutine blocks(lun,nep,wx,wy)
      include 'impli.inc'
      include 'linbuf.inc'
      include 'parm.inc'
      save
      zero = 0.0d0
      ap = zero
      zbeg = zero
      nflg = 0
      nn = 0
      write(lun,12) nep,nflg,wx,wy,brho,beta,gamm1
      write(lun,13) zbeg, ap, strong(1), nn,nn, cname(1)
  12  format(2i4,5f14.7)
  13  format(3f16.7,2i5,2x,a8)
c
c         write envelope file top section
c
      do 20 nn = nbgn, nend
        ltyp = lintyp(nn)
        zbeg = total
        total = total + elong(nn)
        cmtot = 100.0*total
        cml = 100.0*elong(nn)
        glp = elong(nn)*strong(nn)
        gltot = gltot + abs(glp)
        rcm = 100.0*radius(nn)
        kk = nn
        if(ltyp.ne.0) then
            ap = radius(nn)
            write(lun,13)  zbeg, ap, strong(nn), nn, ltyp, cname(nn)
            write(lun,13) total, ap, strong(nn), nn, ltyp, cname(nn)
            kk = nn + 1
        endif
        ap = 0.0
        write(lun,13) total, ap, zero, nn, ltyp, cname(kk)
  20  continue
      return
      end
c--------------------------------------
 
c=====================================================
      subroutine rmquad(zlen,grad,rm)
c
c     computes matrix rm for a quadrupole of length zlen meters
c     and with field gradient of grad tesla/meter
c     if grad>0, the quad is horizontally focusing, vertically defocusin
c     if grad<0, the quad is vertically focusing, horizontally defocusin
c     If grad=0, it is simply a drift.
c       T. Mottershead LANL AT-3, 10/89, copied from MARYLIE routines
c       FQUAD and DQUAD by D. Douglas, 1982.
c------------------------------------------------------------------
      implicit double precision (a-h,o-z)
c
c   beam and kinematic constants:
      include 'parm.inc'
      save
      double precision rm(6,6)
c
c     initialize rm to the identitiy matrix
c
      do 40 i=1,6
      do 20 j=1,6
      rm(i,j)=0.0
  20  continue
      rm(i,i)=+1.0d0
  40  continue
c
c     add drift terms to rm
c
      rm(1,2)= zlen
      rm(3,4)= zlen
      rm(5,6)= zlen/(gamma*beta)**2
c
c      drift return
c
      if(grad.eq.0.0) return
c
c    positive quad (horizontally focusing)
c
      square = grad/brho
      if(square.gt.0.0) then
         ifp = 1
         jfp = 2
         idp = 3
         jdp = 4
      endif
c
c    negative quad (horizontally defocusing)
c
      if(square.lt.0.0) then
         idp = 1
         jdp = 2
         ifp = 3
         jfp = 4
         square = -square
      endif
c
c     coefficients for transverse rm
c
      wavek=dsqrt(square)
      zk=zlen*wavek
      coshkz=(dexp(zk)+dexp(-zk))/(2.0d0)
      sinhkz=(dexp(zk)-dexp(-zk))/(2.0d0)
      coskz=dcos(zk)
      sinkz=dsin(zk)
c
c      matrix in the focusing plane
c
      rm(ifp,ifp) = coskz
      rm(ifp,jfp) = sinkz/wavek
      rm(jfp,jfp) = coskz
      rm(jfp,ifp) = -wavek*sinkz
c
c      matrix in the defocusing plane
c
      rm(idp,idp) = coshkz
      rm(idp,jdp) = sinhkz/wavek
      rm(jdp,jdp) = coshkz
      rm(jdp,idp) = wavek*sinhkz
      return
      end
c=====================================================
      subroutine adump(ldump)
      include 'impli.inc'
      include 'parm.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
      include 'sigbuf.inc'
      save
c      write(ldump,*) '  Beam:', beta, gamm1, brho
c      write(ldump,*) '  Xsig:', xxin, axin, pxin
c      write(ldump,*) '  Ysig:', yyin, ayin, pyin
c      write(ldump,*) '  Zsig:', zzin, azin, pzin
      write(ldump,*) beta, gamm1, brho
      write(ldump,*) xxin, axin, pxin
      write(ldump,*) yyin, ayin, pyin
      write(ldump,*) zzin, azin, pzin
      do 700 nn = nbgn, nend
      write(ldump,617) cname(nn),nn,lintyp(nn),mltyp1(nn),mltyp2(nn),
     * radius(nn), elong(nn), strong(nn)
 617  format(2x,a8,i4,3i3,2f12.7,1pg24.16)
 700  continue
      return
      end
c--------------------------
c=====================================================
      subroutine trcdmp(ltran,ltrace)
      include 'impli.inc'
      include 'parm.inc'
      include 'linbuf.inc'
      include 'actpar.inc'
      include 'sigbuf.inc'
      include 'usrdat.inc'
      save
      ntrc = ltrace-1
      write(ltran,*)  '$data'
      do 500 nn = nbgn, nend
      nnu = nn + ntrc
      tlong = 1000.0*elong(nn)
c
c     drift
c
      if(lintyp(nn).eq.0)  then
         ityp = 1
         write(ltran,488) nnu,cname(nn),nnu,ityp,nnu,tlong
      endif
c
c     quad
c
      if(lintyp(nn).eq.1)  then
         ityp = 3
         write(ltran,488) nnu,cname(nn),nnu,ityp,nnu,strong(nn),tlong
      endif
 488  format(1x,'cmt(',i3,')=''',a,''' nt(',i3,')=',i3,', a(1,',i3
     * ,')=',2(1pg17.9,','))
 500  continue
      write(ltran,*)  '$end'
      return
      end
c=====================================================
      subroutine trnspt(ltran)
      return
      end
c=====================================================
      subroutine madump(ltran)
      return
      end
c=====================================================
      subroutine user2(p)
c
c CTM 16 Apr 2001:
c New USER 2 routine generates a Taylor map from the current map in
c common/map/th(monoms),tmh(6,6) in map.inc.  The Taylor map generation
c is adapted from sunroutine pcmap that writes a taylor map in
c the cartesian basis.
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'map.inc'
      include 'expon.inc'
      include 'taylor.inc'
      include 'parset.inc'
      dimension p(6), iuse(6), zz(6), zvec(monoms)
      character*1 tu, tag(2)
      logical  lterm, lfile
c
c     LUT for rearranging Giorgelli indexing by aberration class
c
      dimension lut(209),kblok(35)
      data lut /1,3,7,9,18,28,30,39,64,84,86,95,120,
     * 175,2,4,8,10,14,19,29,31,35,40,54,65,85,87,91,
     * 96,110,121,155,176,13,15,22,34,36,43,50,55,68,90,92,
     * 99,106,111,124,145,156,179,49,51,58,74,105,107,114,130,141,
     * 146,159,185,140,142,149,165,195,6,12,21,33,42,67,89,98,
     * 123,178,17,24,38,45,57,70,94,101,113,126,158,181,53,60,
     * 76,109,116,132,148,161,187,144,151,167,197,27,48,73,104,129,
     * 184,63,79,119,135,164,190,154,170,200,83,139,194,174,204,209,
     * 5,11,20,32,41,66,88,97,122,177,16,23,37,44,56,69,
     * 93,100,112,125,157,180,52,59,75,108,115,131,147,160,186,143,
     * 150,166,196,26,47,72,103,128,183,62,78,118,134,163,189,153,
     * 169,199,82,138,193,173,203,208,25,46,71,102,127,182,61,77,
     * 117,133,162,188,152,168,198,81,137,192,172,202,207,80,136,191,
     * 171,201,206,205/
c
       data kblok /14,20,18,12,5,10,12,9,4,6,6,3,3,
     *  2,1,10,12,9,4,6,6,3,3,2,1,6,6,3,3,2,1,3,2,1,1/
c
c    User routine for computine and printing Taylor map coefficients
c
c     P(1) = iopt = 0 for compute only
c                 = 1 to also print
c     P(2) = ipq  = 0 for no T,U components
c                 = 1 to print Q components (x,y) only
c                 = 2 to print P components (Px,Py) only
c                 = 3 to print both Q and P components
c                 = 4 to print all components, including TOF
c     P(3) = lun  = file to write machine readback format (0 if none)
c     P(4) = isend for formatted writes (0=none, 1=jof, 2=jodf, 3=both)
c     P(5) = iref = pset number containing reference point to weight by
c                  (0=none)
c     P(6) = print threshhold for weighted coefficient (at ref pt)
c
      iopt = nint(p(1))
      ipq  = nint(p(2))
      lun  = nint(p(3))
      isend= nint(p(4))
      iref = nint(p(5))
      fmin = p(6)
      zero = 0.0d0
      lterm = .false.
      lfile = .false.
      if((isend.eq.1).or.(isend.eq.3)) lterm = .true.
      if((isend.eq.2).or.(isend.eq.3)) lfile = .true.
      zero = 0.0d0
c
c     set components list and reference point
c
      imax = 6
      do 10 j = 1,6
        iuse(j) = j
        zz(j) = pst(j,iref)
  10  continue
      if(ipq.eq.3) imax = 4
      if(ipq.eq.1) then
         imax = 2
         iuse(2) = 3
      endif
      if(ipq.eq.2) then
         imax = 2
         iuse(1) = 2
         iuse(2) = 4
      endif
c
c      fill z monomial vector
c
      call zmono(iopt,zz,zvec)
c
c     reordered dump of Lie polynomials
c
      if(lterm) write(jof,19) zz
      if(lfile) write(jodf,19) zz
  19  format(/,' Ref Z:',6(1pe11.2))
      if(lterm) write(jof,22)
      if(lfile) write(jodf,22)
   22 format(/1h ,'nonzero elements in generating polynomial are :'/)
      kmax = 35
      jmin = 1
      jmax = 0
      do 60 k = 1, kmax
         jmax = jmax + kblok(k)
         nfuse = 0
         do 50 j = jmin, jmax
            n = lut(j)
c
c    test to see if th(n) is large enough to write
c
         if( abs(th(n)) .le. zero) go to 50
         frp = th(n)*zvec(n)
         if(lterm) write(jof,27)n,(expon(m,n),m=1,6),th(n),th(n),frp
         if(lfile) write(jodf,27)n,(expon(m,n),m=1,6),th(n),th(n),frp
  27     format(2x,'f(',i3,')=f( ',3(2i1,1x),')=',f21.6,2(1pe12.2))
         nfuse = nfuse + 1
  50     continue
         if(nfuse.gt.0) then
             if(lterm) write(jof,*) ' '
             if(lfile) write(jodf,*) ' '
         endif
         jmin = jmax + 1
  60  continue
c
c     generate Taylor map
c
      call tugen(th,tmh)
c
c      write formatted dump of taylor map
c
      tag(1) = 't'
      tag(2) = 'u'
c
c          T-matrix elements
c
      if(lterm) write(jof,61) fmin
      if(lfile) write(jodf,61) fmin
  61  format(/,'  T-Matrix elements >',1pe10.3,/)
c
c     optional machine readback file of taylor map
c
      if(lun.gt.0) then
         ii = 0
         jj = 0
         write(lun,68) ii,jj,zero
  68     format(2i4,1pg23.14)
         do 75 ii = 1, 6
            do 70 jj = 1, 6
            if(tmh(ii,jj).ne.zero) write(lun,68) ii,jj,tmh(ii,jj)
  70        continue
  75     continue
         do 85 jj = 1, 6
           do 80 ii = 1, 83
           if(tumat(ii,jj).ne.zero) then
             write(lun,78) jj,ii,tumat(ii,jj),(expon(m,ii),m=1,6)
  78         format(2i4,1pg23.14,3x,6i2)
           endif
  80       continue
  85     continue
      endif
c
c     begin formatted writes
c
      do 90 ii = 1,imax
      i = iuse(ii)
      kmax = 35
      jmin = 1
      jmax = 0
      tu = tag(1)
      do 120 k = 1, kmax
         jmax = jmax + kblok(k)
         nfuse = 0
         do 100 j = jmin, jmax
            n = lut(j)
            if(n.gt.27) go to 100
            tval = tumat(n,i)
            if(dabs(tval).le.fmin) go to 100
            tref = tval*zvec(n)
      if(lterm) then
         write(jof,38) tu,i,n,tu,i,(expon(m,n),m=1,6),tval,tval,tref
      endif
      if(lfile) then
         write(jodf,38) tu,i,n,tu,i,(expon(m,n),m=1,6),tval,tval,tref
      endif
   38 format(2x,a1,i1,'(',i2,')',' = ',a1,i1,'( ',2i1,2(i2,i1),' ) =',
     * f21.6,' =',1pe10.2,'  ref:',1pe10.2)
c      if(lun.gt.0) write(lun,138) i,n,(expon(m,n),m=1,6),tumat(n,i)
  138 format(2i3,6i2,1pg21.13)
         nfuse = nfuse + 1
 100     continue
         jmin = jmax + 1
 120  continue
         if(lterm) write(jof,*) ' '
         if(lfile) write(jodf,*) ' '
  90  continue
c
c      print sorted U-matrix
c
      tu = tag(2)
      if(lterm) write(jof,91) fmin
      if(lfile) write(jodf,91) fmin
  91  format(/,'  U-Matrix elements >',1pe10.3,/)
c      if(lun.gt.0) write(lun,91) fmin
      kmax = 35
      jmin = 1
      jmax = 0
      do 180 k = 1, kmax
         jmax = jmax + kblok(k)
      do 190 ii = 1,imax
         i = iuse(ii)
         nfuse = 0
         do 170 j = jmin, jmax
            n = lut(j)
            if(n.lt.28) go to 170
            if(n.gt.83) go to 170
            tval = tumat(n,i)
            if(dabs(tval).le.fmin) go to 170
            tref = tval*zvec(n)
      if(lterm) then
       write(jof,38) tu,i,n,tu,i,(expon(m,n),m=1,6),tval,tval,tref
      endif
      if(lfile) then
       write(jodf,38) tu,i,n,tu,i,(expon(m,n),m=1,6),tval,tval,tref
      endif
c      if(lun.gt.0) write(lun,138) i,n,(expon(m,n),m=1,6),tumat(n,i)
         nfuse = nfuse + 1
 170     continue
         if(nfuse.gt.0) then
            if(lterm) write(jof,*) ' '
            if(lfile) write(jodf,*) ' '
         endif
 190  continue
         jmin = jmax + 1
 180  continue
      return
      end
c
c=====================================================
      subroutine zmono(nopt,zz,zvec)
      include 'impli.inc'
      include 'param.inc'
      include 'expon.inc'
      include 'vblist.inc'
      dimension zz(6), zvec(monoms)
c      do 20 n = 1, monoms
c      j = 1
c      kode = expon(j,n)
c  18  continue
c      j = j + 1
c      kode = 10*kode + expon(j,n)
c      if(j.lt.6) go to 18
c      write(89,19) n,kode,(expon(j,n),j=1,6),(vblist(k,n),k=1,4)
c  19  format(2i8,2x,6i2,'  vb:',4i2)
c  20  continue
      one = 1.0d0
      nord = 2
      do 100 n = 1, monoms
      if(n.eq.28) nord = 3
      if(n.eq.84) nord = 4
      prod = one
        do 50 j = 1,nord
           prod = prod*zz(vblist(j,n))
  50    continue
      zvec(n) = prod
 100  continue
      return
      end
c=====================================================
      subroutine tugen(fa,fm)
c  pcmap routine to print m,f3,f4 and t,u.
c Written by D. Douglas ca 1982 and modified by Rob Ryne
c and Alex Dragt ca 1986
ctm modified to just generate the taylor map in 'taylor.inc'
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'expon.inc'
      include 'pbkh.inc'
      include 'taylor.inc'
c
      dimension fa(monoms),fm(6,6),fsav(monoms)
      dimension t(monoms),u(monoms),u2(monoms)
c
c  prepare for higher order matrix generation
c
ctm   write(6,*) ' in TUGEN'
      do 10 jj = 1, monoms
      fsav(jj) = fa(jj)
  10  continue
      call brkts(fa)
c
c  procedure for generating t-matrix
c
        do 35 i=1,6
        call xform(pbh(1,i),2,fm,i-1,t)
        do 36 n=7,27
        tumat(n,i) = t(n)
   36   continue
   35   continue
c
c  procedure for generating U-matrix
c
      do  44 i=1,6
         call xform(pbh(1,i),3,fm,i-1,u)
         call xform(pbh(1,i+6),3,fm,1,u2)
         do  45 n=28,83
            u(n)=u(n)+u2(n)/2.d0
            tumat(n,i) = u(n)
   45    continue
   44 continue
      return
      end
c=====================================================
      subroutine user3(p)
c
c          user3 dumps the specified range of ucalc in full
c          precision on unit nfile
 
      include 'impli.inc'
      include 'usrdat.inc'
      include 'coment.inc'
      include 'mldex.inc'
      include 'param.inc'
      include 'map.inc'
      dimension p(6)
c
c     parameters
c
c     p(1) = nfile = unit to write on
c     p(2) = kmin = array min to write
c     p(3) = kmax = array max to write
c     p(4) = job: job = 0 to dump range of ucalc to unit nfile
c                 job = 1 to dump current matrix to unit nfile
c                 job = 2 to list type codes in group order (call listyp)
c                 job = 3 lists all type codes in alphabetical order(allist)
c                 job = 4 to edit comment block
c                 job = 5 to call Marylie dump of all commons (echo input format)
c                 job = 6 to write a LaTeX file of all type codes (texlis)
c     p(5) = kform = format flag for dumps = 0 for full * precision
c     p(6) = nextra = spare flag
c
      nfile = nint(p(1))
      if(nfile.le.0) return
      kmin = nint(p(2))
      kmax = nint(p(3))
      job = nint(p(4))
      kform = nint(p(5))
      nextra = nint(p(6))
c
c           ucalc dump is job 0
c
      if(job.eq.0) then
      write(nfile,*) '  UCALC DUMP'
         do 100  k=kmin, kmax
           if(ucalc(k).ne.0.0) write(nfile,*) k, ucalc(k)
 100     continue
      endif
c
c           matrix print is job 1
c
      if(job.eq.1) then
         write(nfile,*) ' --------   Current Matrix  ----------'
         write(nfile,130) (tmh(1,j),j=1,2),(tmh(3,j),j=3,4)
         write(nfile,130) (tmh(2,j),j=1,2),(tmh(4,j),j=3,4)
  130    format('  Rx:',2f15.9,'    Ry:',2f15.9)
      endif
c
c        list typecodes in this version of marylie
c
      if(job.eq.2) call listyp(nfile)
      if(job.eq.3) call allist(nfile)
c
c      comment editing
c
      if(job.eq.4) then
      write(jof,*)'  Edit comments (e.g. to show purpose of this run)'
          call caedit(maxcmt, np, mline)
      endif
c
c           marylie input dump
c
      if(job.eq.5) call dump(nfile)
c
c       TeX listing of type codes
c
      if(job.eq.6) call texlis(nfile)
      return
      end
c=====================================================
      subroutine listyp(lun)
c
c       listyp  lists the type codes in the current vesion of marylie
c       T. Mottershead   LANL  AT-3  17 Jan 91
c------------------------------------------------
      include 'codes.inc'
      character*16 tcode(9), word
      dimension numcmd(9), nord(9), jlo(2),jhi(2)
      data nord /1,4,7,8,9,2,3,5,6/
      jlo(1)=1
      jhi(1)=5
      jlo(2)=6
      jhi(2)=9
      ib = ichar(' ')
c
c     scan for null string to indicate number of commands
c
      ntot = 0
      do 20 j = 1, 9
         max = 0
         do 10 n=1,40
         word = ltc(j,n)
         iw = ichar(word(1:1))
         if(iw.le.ib) go to 15
         max = max + 1
  10     continue
  15     numcmd(j) = max
         ntot = ntot + max
  20  continue
c
c        report total
c
      write(lun,*) ' Marylie Type Codes '
      write(lun,*)
      write(lun,*) ntot,' type codes in this version'
c
c        writeout
c
      do 90  kk = 1,2
      jmin = jlo(kk)
      jmax = jhi(kk)
      nmax = 0
      do 30 j=jmin, jmax
      mm = numcmd(nord(j))
      if(mm.gt.nmax) nmax = mm
  30  continue
      write(lun,*)
      write(lun,33) (nord(k),k=jmin,jmax)
  33  format(' index ',2x,5('nrp (G',i1,')',4x))
      do 50 n = 1 ,nmax
         do 40 j=jmin, jmax
         tcode(j) = '        '
         word = ' '
         jju = nord(j)
         if(n.le.numcmd(jju)) word = ltc(jju,n)
         iw = ichar(word(1:1))
         if(iw.gt.ib) tcode(j) = word
  40  continue
      write(lun,43) n,(nrp(nord(j),n), tcode(j), j=jmin,jmax)
  43  format(i4,5x,5(i2,2x,a8))
  50  continue
      write(lun,31) (numcmd(nord(j)), j=jmin, jmax)
  31  format(' total:',i8,4i12)
  90  continue
      return
      end
c
c=====================================================
      subroutine allist(lun)
c
c       allist makes an alphabetical table of all
c       type codes in the current vesion of marylie
c       T. Mottershead   Hillsboro, OR The Ides of March 2011
c------------------------------------------------
      include 'codes.inc'
      character*16  tcode(360), word
      dimension nseq(360), npars(360), kind(360)
c
c     fill unsorted tcode array, stopping on null string for each group
c
      max = 0
      do 20 j = 1, 9
         do 10 n=1,40
         word = ' '
         word = ltc(j,n)
         icw = ichar(word(1:1))
c         write(66,17) j,n,nrp(j,n),ltc(j,n),word,icw
c  17     format(' j,n,nrp=',3i4,' ltc,word=: ',2a,' icw=',i3)
c         if(word.eq.' ') go to 10
         if(icw.lt.30) go to 10
         max = max + 1
         tcode(max) = word
         npars(max) = nrp(j,n)
         kind(max) = j
  10     continue
  20  continue
c
c      alphabetize and write Latex table
c
      write(6,*) max,'= max = number of type codes'
      call lexort(max,tcode,nseq)
      do 50  k = 1,max
      jj = nseq(k)
      write(lun,37) tcode(jj),npars(jj),kind(jj),k
  37  format(1x,a,10x,3i5)
  50  continue
      return
      end
c=====================================================
      subroutine texlis(lun)
c
c       texlis makes an alphabetical LaTeX table of all
c       type codes in the current vesion of marylie
c       T. Mottershead   LANL  AT-3  16 Jan 92
c------------------------------------------------
      include 'codes.inc'
      character*16 tcode(360), word
      dimension nseq(360), npars(360), kind(360)
c
c     fill unsorted tcode array, stopping on null string for each group
c
      long = 40
      max = 0
      do 20 j = 1, 9
         do 10 n=1,40
         word = ltc(j,n)
         if(ichar(word(1:1)).eq.0) go to 20
         max = max + 1
         tcode(max) = word
         npars(max) = nrp(j,n)
         kind(max) = j
  10     continue
  20  continue
c
c      alphabetize and write Latex table
c
      call lexort(max,tcode,nseq)
      call headoc(lun,'MaryLie Type Codes')
      call tabdef(lun,5)
      write(lun,27)
  27  format(' No. & Code  & Kind & NRP  & \\')
      write(lun,*) '  \hline'
      do 50  k = 1,max
      jj = nseq(k)
      write(lun,37) k, tcode(jj), kind(jj), npars(jj)
  37  format(i6,' & ',a,' &',i3,' &',i3,' &  \\')
      if(mod(k,long).eq.0) then
        call tabend(lun)
        write(lun,*) ' \newpage'
        call tabdef(lun,5)
        write(lun,27)
        write(lun,*) '  \hline'
      endif
  50  continue
      call tabend(lun)
      call endoc(lun)
      return
      end
c=====================================================
      subroutine headoc(lun,title)
c
c       writes LaTeX document header to LUN from fortran programs.
c
c      T. Mottershead  AT-3   16 Jan 91
c-------------------------------------------------
      character title*(*)
c
c           Write Document opening
c
      write(lun,*) ' \documentstyle[12pt]{article} '
      write(lun,*) ' \setlength{\textwidth}{6.5in}'
      write(lun,*) ' \setlength{\oddsidemargin}{0in} '
      write(lun,*) ' \setlength{\topmargin}{-0.5in}'
      write(lun,*) ' \setlength{\textheight}{9in} \pagestyle{empty}'
      write(lun,*) ' \begin{document}'
      write(lun,*) ' '
      write(lun,*) ' {\large '
      write(lun,14) title
  14  format(' \begin{center} {\bf ',a,' }  \\
     * \today    \end{center} }')
      return
      end
c=====================================================
      subroutine tabdef(lun,ncol)
c
c       writes LaTeX table headers to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c-------------------------------------------------
      write(lun,*), ' \begin{center}'
      write(lun,*), ' \begin{tabular}{',('|c',j=1,ncol),'|}'
      write(lun,*), ' \hline'
      return
      end
c=====================================================
      subroutine tabend(lun)
c
c       writes LaTeX table end to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c-------------------------------------------------
      write(lun,*), ' \hline \end{tabular} \end{center}'
      return
      end
c=====================================================
      subroutine endoc(lun)
c
c       writes LaTeX document end to LUN from fortran programs.
c      T. Mottershead  AT-3   16 Jan 91
c-------------------------------------------------
      write(lun,*), ' \end{document}'
      return
      end
c=====================================================
      subroutine caedit(max, nl, text)
c
c     caedit is a simple character array editing routine.
c     text(i), i=1,nl is the array of character strings to be edited.
c          The length of the strings is determined internally.
c          All strings are the same length.
c     nl = actual number of character strings, which may be increased
c          or decreased by this routine.
c     max = maximum number of character strings allowed by the calling
c           program.
c
c           T. Mottershead/ LANL  AT-3 /  24 Jan 91
c=========================================================
      character*(*) text(max)
      character*74  card
c
c         print array
c
  10  do 20 j=1, nl
      card = text(j)
      write(6,13) j,card
  13  format(i3,a)
  20  continue
c
c         command input
c
  30  write(6,32)
  32  format(' enter command: +N = insert new line N, 0 = quit,',
     * ' -N = delete line N')
      read(5,*) kl
      if(kl.eq.0) return
c
c     delete line N
c
      if(kl.lt.0) then
         nd = -kl
         if(nd.gt.nl) go to 30
         last = nd
         card = ' '
         write(6,37) last
  37     format(' enter last line to delete <',i3,'>:')
         read(5,64) card
         if(card.ne.' ') then
            read(card,*) last
         endif
  40     nl = nl-1
         do 50 n = nd,nl
         text(n) = text(n+1)
  50     continue
         if(last.le.nd) go to 10
         last = last-1
         go to 40
      endif
c
c       add new line N
c
      if(kl.gt.0) then
         if(kl.gt.nl) kl = nl+1
         write(6,*) ' Enter new line(s) <CR>=print'
  60     card = ' '
         read(5,64) card
  64     format(a)
         if(card.eq.' ') go to 10
         nl = nl+1
         do 90 n = nl-1, kl, -1
         text(n+1) = text(n)
  90     continue
         text(kl)=card
         kl = kl+1
         go to 60
      endif
      go to 10
      end
c-----------------------------------------------------
      subroutine user4(p)
cccc      call u4test
      return
      end
c=====================================================
      subroutine user5(pp)
      include 'impli.inc'
      include 'files.inc'
      include 'linbuf.inc'
      include 'usrdat.inc'
      parameter (maxa=200, maxpt=2000)
      dimension ktyp(maxa), path(maxa), bend(maxa), pp(*)
      dimension xx(maxpt),yy(maxpt),key(maxpt)
      logical ltty, ldsk
      ltty = .true.
      ldsk = .true.
c      write(jof,*) ' USER5 Parameter set processor called for arcs.'
      npsu = 7
      call filin
c      type *, npsu,'=npsu',npst,' of these are in loop:'
      if(npst.le.0) return
      zero = 0.0d0
      one = 1.0d0
      two = 2.0d0
      pi = 3.14159265358979d0
      radeg = pi/180.0d0
      rdef = 192.0d0/pi
      rru = rdef
c================================================
      iop = nint(pp(1))
      xbeg = pp(2)
      ybeg = pp(3)
      if(iop.eq.1) then
         xbeg = pp(2)*dcos(radeg*pp(3))
         ybeg = pp(2)*dsin(radeg*pp(3))
      endif
      angle = pp(4)
      iud = nint(pp(5))
      if(iud.lt.0) then
         iud = -iud
         ltty = .false.
         ldsk = .false.
      endif
      lun = nint(pp(6))
c
c        echo and process contents of loop
c
      arclen = zero
      do 20 j = 1, npst
        job = nint(aap(j))
        path(j) = bbp(j)
        bend(j) = ccp(j)
        theta = radeg*bend(j)
        if(job.eq.1) rru = dabs(path(j)/theta)
        if(job.eq.2) path(j) = dabs(rru*radeg*bend(j))
        if(job.eq.3) bend(j) = path(j)/(rru*radeg)
        ktyp(j) = nint(ddp(j))
        arclen = arclen + path(j)
        if(ltty) then
         write(jof,17) j,job,ktyp(j),path(j),bend(j),rru,bbp(j),ccp(j)
        endif
        if(ldsk) then
         write(jodf,17) j,job,ktyp(j),path(j),bend(j),rru,bbp(j),ccp(j)
        endif
  17    format(i6,2i3,5f13.6)
  20  continue
c
c      data in, draw arcs
c
      xx(1) = xbeg
      yy(1) = ybeg
      delta = two
      call arcs(npst,ktyp,path,bend,angle,delta,maxpt,xx,yy,key,nn)
      if(ltty) write(jof,32) arclen,nn
      if(ldsk) write(jodf,32) arclen,nn
  32  format(' Total path length =',f12.6,' in',i4,' points')
      if(ltty) write(jof,34) xx(nn),yy(nn),angle
      if(ldsk) write(jodf,34) xx(nn),yy(nn),angle
  34  format(' Endpoint:',f13.6,' = x',f13.6,' = y  at',f10.4,' deg.')
c
c       optional save to ucalc
c
      if((iud.gt.0).and.(iud.le.250)) then
        ucalc(iud) = arclen
        ucalc(iud+1) = angle
        ucalc(iud+2) = xx(nn)
        ucalc(iud+3) = yy(nn)
      endif
c
c       optional coordinate dump
c
      if(lun.le.0) return
      do 50 j = 1,nn
        write(lun,47) key(j),xx(j),yy(j)
  47    format(i5,2f12.6)
  50  continue
      return
      end
c=====================================================
      subroutine arcs(nseg,ktyp,path,bend,angle,delta,max,xx,yy,key,nn)
c
c        Generates (xx(i),yy(i),i=1,nn) for drawing a sequence of
c        connected arcs. The curve starts from x(1),y(1) at input angle,
c        and ends at x(nn),y(nn) headed in direction of final angle.
c        nseg = number of path segments k=1,nseg
c        path(k) = arclength of kth segment
c        bend(k) = bend angle in degrees of kth segment
c        ktyp(k) = ID or type code for kth segment
c        delta = degrees of bend per plot point (max allowed)
c        max = maximum number of output points allowed
c
c        C. T. Mottershead  24 Aug 99
c=====================================================
      implicit double precision (a-h,o-z)
      dimension path(*),bend(*),ktyp(*),xx(max),yy(max),key(max)
      pi = 3.14159265358979d0
      radeg = pi/180.0d0
      bmin = 0.001
      one = 1.0d0
      xu = xx(1)
      yu = yy(1)
      key(1) = ktyp(1)
      lout = 44
      nn = 0
c
c       loop over all path segments
c
      do 200 k = 1,nseg
        theta = radeg*angle
        angle = angle + bend(k)
        xs = path(k)*dcos(theta)
        ys = path(k)*dsin(theta)
        write(lout,207) k,ktyp(k),path(k),bend(k),angle,xs,ys
  207   format(' Arc:',2i4,5f12.4)
c
c     straight segment
c
        if(dabs(bend(k)).lt.bmin) then
          nn = nn+1
          if(nn.gt.max) return
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          nn = nn+1
          if(nn.gt.max) return
          xu = xu + xs
          yu = yu + ys
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          np = 1
c          write(lout,97) nn,key(nn),xx(nn),yy(nn),xs,ys
          go to 190
        endif
c
c      curved segment (starts with duplicate of previous point)
c
        nn = nn+1
        if(nn.gt.max) return
        xx(nn) = xu
        yy(nn) = yu
        key(nn) = ktyp(k)
        phi = radeg*bend(k)
        radius = path(k)/phi
        np = nint(bend(k)/delta)
        if(np.lt.0) np = -np
        if(np.eq.0) np = 1
        fnp = float(np)
        alpha = phi/fnp
        xs = xs/fnp
        ys = ys/fnp
        cosa = dcos(alpha)
        sina = dsin(alpha)
        sfac = sina/alpha
        cfac = (one-cosa)/alpha
        x0 = xs
        y0 = ys
c
c     update plot points in arc
c
        do 100 j = 1,np
          dx = sfac*xs - cfac*ys
          dy = sfac*ys + cfac*xs
          xu = xu + dx
          yu = yu + dy
          nn = nn+1
          if(nn.gt.max) return
          xx(nn) = xu
          yy(nn) = yu
          key(nn) = ktyp(k)
          xs = cosa*x0 - sina*y0
          ys = cosa*y0 + sina*x0
          x0 = xs
          y0 = ys
  100   continue
  190 continue
  200 continue
      return
      end
c
ccccccccccccccccc     Blocks     cccccccccccccccccccccccccccccc
c
c=====================================================
      subroutine user6(p,f,fm)
c
c     alternate sigma matrix input
c
c              p(1) = alphax
c              p(2) = betax
c              p(3) = epsx
c              p(4) = alphay
c              p(5) = betay
c              p(6) = epsy
      include 'impli.inc'
      include 'param.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'sigbuf.inc'
      dimension f(monoms)
      dimension fm(6,6)
      dimension p(6), psig(6)
      alphax = p(1)
      betax = p(2)
      epsx = p(3)
      if(epsx.le.0.0) epsx = 1.0
      alphay = p(4)
      betay = p(5)
      epsy = p(6)
      if(epsy.le.0.0) epsy = 1.0
      psig(1) = sqrt(epsx*betax)
      psig(2) = alphax
      gamx = (1.0 + alphax**2)/betax
      psig(3) = sqrt(epsx*gamx)
      psig(4) = sqrt(epsy*betay)
      psig(5) = alphay
      gamy = (1.0 + alphay**2)/betay
      psig(6) = sqrt(epsy*gamy)
      call user8(psig, f, fm)
      return
      end
c=====================================================
      subroutine user7(p,g,mh)
      include 'impli.inc'
      dimension p(6)
ccc      call u7test(p)
      return
      end
c=====================================================
      subroutine user8(p,f,rm)
c  This routine applies the current map in (f,rm) to the initial sigma
c  matrix given in the parameter set p = (s11,s12,s22,s33,s34,s44), to
c  produce the output sigma matrix of the same form, which is saved in
c  ucalc(1:6) in the same order.
c
c    parameters: p(1) = xmax = sqrt(betax*xemit)
c                p(2) = alphax
c                p(3) = thetax = sqrt(gammax*xemit)
c                p(4) = ymax = sqrt(betay*yemit)
c                p(5) = alphay
c                p(6) = thetay = sqrt(gammay*yemit)
c     Tom Mottershead  LANL  31-Mar-89
c     f Neri           LANL   5-Dec-90  Full 4D calculation.
c     Elliptical beam kicks.  8-Apr-92
c====================================================================
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'parm.inc'
      include 'sigbuf.inc'
      include 'parset.inc'
      dimension f(monoms)
      dimension rm(6,6)
      dimension p(6)
      logical  lterm, lfile
      isend = 3
      level = 1
      lout = 77
      lterm = .true.
      lfile = .true.
c
c     save input parameters to sigbuf for user1
c
c      load input beam from #menu
c
      if(p(1).gt.0.0) then
         xxin = p(1)
         axin = p(2)
         pxin = p(3)
         yyin = p(4)
         ayin = p(5)
         pyin = p(6)
         go to 100
      endif
c
c       optional input specifications
c
      job = nint(-p(1))
      if(job.gt.0) then
         inform = nint(p(2))
         lunin  = nint(p(3))
         lout  = nint(p(4))
c
c       output switches
c
         level  = nint(p(5))
         isend  = nint(p(6))
         lterm = .false.
         lfile = .false.
         if((isend.eq.1).or.(isend.eq.3)) lterm = .true.
         if(isend.ge.2) lfile = .true.
c
c      aberrations only case
c
         if(job.eq.5) go to 300
         if(inform.ne.0) go to 100
         if(lunin.gt.0) then
            rewind(lunin)
            read(lunin,*) xxin,axin,pxin,yyin,ayin,pyin
         endif
         if(lunin.lt.0) then
            kpu = -lunin
            if((kpu.lt.1).or.(kpu.gt.maxpst)) then
              write(jof,*)  ' USR8 Error:',kpu,' is not a valid pset'
              call myexit
            endif
            xxin = pst(1,kpu)
            axin = pst(2,kpu)
            pxin = pst(3,kpu)
            yyin = pst(4,kpu)
            ayin = pst(5,kpu)
            pyin = pst(6,kpu)
         endif
      endif
c
c         end of input options, map the beam
c
 100  continue
      if(iquiet.eq.1) then
         lterm = .false.
         lfile = .false.
      endif
c
c        clear sigma matrices
c
      do 150 i=1,4
      do 150 j=1,4
        sig0(i,j) = 0.d0
        sigf(i,j) = 0.d0
  150 continue
c
c        Initialize Sigma Matrix from input beam parameter block
c
      sig0(1,1) = xxin**2
      sig0(1,2) = -xxin*axin*pxin/sqrt(1.0d0 + axin**2)
      sig0(2,2) = pxin**2
      sig0(2,1) = sig0(1,2)
      sig0(3,3) = yyin**2
      sig0(3,4) = -yyin*ayin*pyin/sqrt(1.0d0 + ayin**2)
      sig0(4,4) = pyin**2
      sig0(4,3) = sig0(3,4)
c
c             Transport Ellipse ( 4D )
c
      do 250 i=1,4
      do 250 j=1,4
         do 200 k=1,4
         do 200 n=1,4
            sigf(i,j) = sigf(i,j) + rm(i,n)*sig0(n,k)*rm(j,k)
 200     continue
 250  continue
c      type *, ' sigf:', sigf
c
c      normal moments
c
      xxf = sigf(1,1)
      xpf = sigf(1,2)
      ppf = sigf(2,2)
      yyf = sigf(3,3)
      yqf = sigf(3,4)
      qqf = sigf(4,4)
c
c       cross moments
c
      xyf = sigf(1,3)
      ypf = sigf(2,3)
      xqf = sigf(1,4)
      detc = xxf*yyf-xyf**2
c      type *,' detc =',detc
      xfocus = (yyf*xpf-xyf*ypf)/detc
      yfocus = (xxf*yqf-xyf*xqf)/detc
      xcross = (xxf*ypf-xyf*xpf)/detc
      ycross = (yyf*xqf-xyf*yqf)/detc
c
c       normal emittances
c
      exsq = xxf*ppf - xpf**2
      if(exsq.le.0.0) exsq = 1.e-30
      ex = sqrt(exsq)
      betax = xxf/ex
      gammax = ppf/ex
c
      eysq = yyf*qqf - yqf**2
      if(eysq.le.0.0) eysq = 1.e-30
      ey = sqrt(eysq)
      betay = yyf/ey
c
      ucalc(1) = sqrt(xxf)
      ucalc(2) = -xpf/ex
      ucalc(3) = sqrt(ppf)
      ucalc(4) = sqrt(yyf)
      ucalc(5) = -yqf/ey
      ucalc(6) = sqrt(qqf)
      ucalc(7) = betax
      ucalc(8) = ex
      ucalc(9) = xfocus
      ucalc(10) = xcross
      ucalc(11) = betay
      ucalc(12) = ey
      ucalc(13) = yfocus
      ucalc(14) = ycross
      ucalc(50) = xxf-yyf
c
c        save output beam
c
c
      write(6,*) 'U8OUT: lout=',lout
      if (lout.gt.0) write(lout,*) (ucalc(j), j= 1,6)
      if (lout.lt.0) then
         kpu = -lout
         if((kpu.lt.1).or.(kpu.gt.maxpst)) then
            write(jof,*)  ' USR8 Error:',kpu,' is not a valid pset'
            call myexit
         endif
         do 400 j=1,6
         pst(j,kpu) = ucalc(j)
 400     continue
      endif
c
c     copy output back over input if job=3
c
      if(job.eq.3) then
         xxin = ucalc(1)
         axin = ucalc(2)
         pxin = ucalc(3)
         yyin = ucalc(4)
         ayin = ucalc(5)
         pyin = ucalc(6)
      endif
c
c  Change of Sigma with energy: first order ( f. Neri 12-5-1990 ).
c
      dsig(1,1) = -2*f(38)*sigf(1,1) - 4*f(53)*sigf(1,2) -
     -  2*f(57)*sigf(1,3) - 2*f(60)*sigf(1,4)
      dsig(1,2) = 2*f(33)*sigf(1,1) + f(42)*sigf(1,3) +
     -  f(45)*sigf(1,4) - 2*f(53)*sigf(2,2) -
     -  f(57)*sigf(2,3) - f(60)*sigf(2,4)
      dsig(2,2) = 4*f(33)*sigf(1,2) + 2*f(38)*sigf(2,2) +
     -  2*f(42)*sigf(2,3) + 2*f(45)*sigf(2,4)
      dsig(1,3) = -(f(45)*sigf(1,1)) - f(60)*sigf(1,2) -
     -  f(38)*sigf(1,3) - f(70)*sigf(1,3) -
     -  2*f(76)*sigf(1,4) - 2*f(53)*sigf(2,3) -
     -  f(57)*sigf(3,3) - f(60)*sigf(3,4)
      dsig(2,3) = -(f(45)*sigf(1,2)) + 2*f(33)*sigf(1,3) -
     -  f(60)*sigf(2,2) + f(38)*sigf(2,3) -
     -  f(70)*sigf(2,3) - 2*f(76)*sigf(2,4) +
     -  f(42)*sigf(3,3) + f(45)*sigf(3,4)
      dsig(3,3) = -2*f(45)*sigf(1,3) - 2*f(60)*sigf(2,3) -
     -  2*f(70)*sigf(3,3) - 4*f(76)*sigf(3,4)
      dsig(1,4) = f(42)*sigf(1,1) + f(57)*sigf(1,2) +
     -  2*f(67)*sigf(1,3) - f(38)*sigf(1,4) +
     -  f(70)*sigf(1,4) - 2*f(53)*sigf(2,4) -
     -  f(57)*sigf(3,4) - f(60)*sigf(4,4)
      dsig(2,4) = f(42)*sigf(1,2) + 2*f(33)*sigf(1,4) +
     -  f(57)*sigf(2,2) + 2*f(67)*sigf(2,3) +
     -  f(38)*sigf(2,4) + f(70)*sigf(2,4) +
     -  f(42)*sigf(3,4) + f(45)*sigf(4,4)
      dsig(3,4) = f(42)*sigf(1,3) - f(45)*sigf(1,4) +
     -  f(57)*sigf(2,3) - f(60)*sigf(2,4) +
     -  2*f(67)*sigf(3,3) - 2*f(76)*sigf(4,4)
      dsig(4,4) = 2*f(42)*sigf(1,4) + 2*f(57)*sigf(2,4) +
     -  4*f(67)*sigf(3,4) + 2*f(70)*sigf(4,4)
c
      ds11 = dsig(1,1)
      ds12 = dsig(1,2)
      ds22 = dsig(2,2)
      ds33 = dsig(3,3)
      ds34 = dsig(3,4)
      ds44 = dsig(4,4)
c
      ucalc(21) = ds11
      ucalc(22) = ds12
      ucalc(23) = ds22
      ucalc(24) = ds33
      ucalc(25) = ds34
      ucalc(26) = ds44
c
      dfxdp = -beta*(xxf*ds12 - xpf*ds11)/(xxf**2)
      dfydp = -beta*(yyf*ds34 - yqf*ds33)/(yyf**2)
      ucalc(15) = dfxdp
      ucalc(16) = dfydp
      ucalc(37) = dfxdp - dfydp
      ucalc(38) = dfxdp**2 + dfydp**2
c
c        rms cubic aberration amplitude
c
c      aa = f(84)
c      bb = f(95)
c      cc = f(175)
c      qq0 = 5.0*(aa**2 + cc**2) + 0.5*bb**2 + bb*(aa+cc)
c      abcube = sqrt(qq0)
c    Get Ellipticity from user1:
       rx = ucalc(46)
       if (rx .eq. 0.0 ) then
         write(6,*) ' user8: Ellipticity not set in user1, so take r=1'
         rx = 1.0
       endif
       ry = 1/rx
c
c    Mathematica formula for quadratic and cubic kicks
c
 300  continue
c      abquad=3.375*(f(28)**2 + f(64)**2) + 0.875*(f(30)**2 + f(39)**2)
c     * + 0.75*(f(30)*f(64) + f(28)*f(39))
c      abcube=5.*(f(84)**2 + f(175)**2) + f(95)*(f(175)+f(84)) +
c     * 0.5*f(95)**2 + 0.75*f(86)*f(120) + 0.875*(f(86)**2+f(120)**2)
c
      abquad=3.375*rx**4*f(28)**2 + 0.375*rx**4*f(30)**2 +
     -  0.5*rx**2*ry**2*f(30)**2 +
     -  0.75*rx**2*ry**2*f(28)*f(39) +
     -  0.5*rx**2*ry**2*f(39)**2 + 0.375*ry**4*f(39)**2 +
     -  0.75*rx**2*ry**2*f(30)*f(64) + 3.375*ry**4*f(64)**2
      abcube=5.*rx**6*f(84)**2 + 0.3125*rx**6*f(86)**2 +
     -  0.5625*rx**4*ry**2*f(86)**2 +
     -  rx**4*ry**2*f(84)*f(95) +
     -  0.25*rx**4*ry**2*f(95)**2 +
     -  0.25*rx**2*ry**4*f(95)**2 +
     -  0.375*rx**4*ry**2*f(86)*f(120) +
     -  0.375*rx**2*ry**4*f(86)*f(120) +
     -  0.5625*rx**2*ry**4*f(120)**2 +
     -  0.3125*ry**6*f(120)**2 + rx**2*ry**4*f(95)*f(175) +
     -  5.*ry**6*f(175)**2
c
c
      ptquad=3.375*(f(49)**2 + f(74)**2) + 0.875*(f(51)**2 + f(58)**2)
     * + 0.75*(f(51)*f(74) + f(49)*f(58))
c
      ptcube=5.*(f(140)**2 + f(195)**2) + f(149)*(f(195)+f(140)) +
     * 0.5*f(149)**2 + 0.75*f(142)*f(165) + 0.875*(f(142)**2+f(165)**2)
c
      abfour = 0.0
      abfive = 0.0
      ucalc(17) = sqrt(abquad)
      ucalc(18) = sqrt(abcube)
      ucalc(19) = sqrt(abfour)
      ucalc(20) = sqrt(abfive)
      ucalc(60) = sqrt(ptquad)
      ucalc(61) = sqrt(ptcube)
c
      cx = -2.0*beta*f(33)
      cy = -2.0*beta*f(67)
      ucalc(41) = cx
      ucalc(42) = cy
      ucalc(43) = cx**2 + cy**2
      ucalc(44) = cx - cy
      ucalc(45) = ucalc(1)*ucalc(4)
      cxim = -2.0*beta*f(53)
      ucalc(62) = cxim
      if(f(38).ne.0.0) xvs = -2.0*f(53)/f(38)
      ucalc(63) = xvs
      cyim = -2.0*beta*f(76)
      ucalc(64) = cyim
      if(f(70).ne.0.0) yvs = -2.0*f(76)/f(70)
      ucalc(65) = yvs
      if(lterm) call u8out(jof,level)
      if(lfile) call u8out(jodf,level)
      return
      end
c=====================================================
      subroutine u8out(lun,level)
      include 'impli.inc'
      include 'param.inc'
      include 'usrdat.inc'
      include 'sigbuf.inc'
      save
      character*3 tag
      write(lun,21)
  21  format(/' Output beam parameters:',/,'  x=u(1)    ',
     *'   alfx=u(2)    px=u(3)      y=u(4)       alfy=u(5) ',
     *'   py=u(6)   ')
      write(lun,25) (ucalc(j), j=1,6)
  25  format(6(1pe13.5))
      write(lun,22)
  22  format(6x,'beta',15x,'emittance',10x,'focus',14x,'cross-focus')
      tag = ' X:'
      write(lun,24), tag, (j,ucalc(j),j=7,idiv)
  24  format(a3,i3,'>',1pe14.7,3(i4,'>',1pe14.7))
      tag = ' Y:'
      write(lun,24), tag, (j,ucalc(j),j=11,14)
      if(level.lt.1) return
c
c      also write aberration coefficients for level 1
c
      write(lun,27)
  27  format(/,' Chromatic and RMS Geometric Aberration Coefficients:',
     */,'     dFx/dp=u(15)      dFy/dp=u(16)      K2(R/m^2)=u(17)',
     *'   K3(R/m^3)=u(18)')
      write(lun,29) (ucalc(j),j=15,18)
  29  format(4(1pg18.7))
      if(level.lt.2) return
c
c     also write dsigma/de for level 2
c
      write(lun,23)
  23  format(/' Energy dependence of sigma matrix:',/,'  ds11=u(21)',
     *'   ds12=u(22)   ds22=u(23)   ds33=u(24)   ds34=u(25)',
     *'   ds44=u(26)')
      write(lun,25) (ucalc(j),j=21,26)
c      write(lun,*) ' Sigma matrix:'
c      do 201 i=1,4
c      write(lun,26) (sigf(i,j),j=1,4)
c  201 continue
c 26  format(3x,4(1pe15.7))
c         write(lun,*) ' Energy derivative of Sigma matrix:'
c      do 202 i=1,4
c         write(lun,26) (dsig(i,j),j=1,4)
c  202 continue
      return
      end
c=====================================================
      subroutine user9(p,fa,fm)
c
c     This routine does arithmetic operations on the current map.
c         C. T. Mottershead  LANL AOT-1 / 4 April 96
ctm       Feb 2012  completed arithmetic (add,div,mult) for R,f,T,u in
ctm       most backward compatible way possible.
c--------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'taylor.inc'
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6), xv(30)
c
      job = int(p(1))
      jndx = int(p(2))
      kndx = int(p(3))
      nusr = int(p(4))
      aa = p(5)
      bb = p(6)
      low = 16
      idiv = 100
c      type *,' USR9:',job,jndx,kndx,nusr,aa,bb
      if(job.lt.0) then
         job = -job
         idxa = nint(p(5))
         if((idxa.gt.0).and.(idxa.le.250)) aa = ucalc(idxa)
         idxb = nint(p(6))
         if((idxb.gt.0).and.(idxb.le.250)) bb = ucalc(idxb)
      endif
      zero = 0.0d0
      value = zero
c
c       job = 0 is to set constants into ucalc
c
      if(job.eq.0) then
         if((jndx.gt.0).and.(jndx.le.250)) ucalc(jndx) = aa
         if((kndx.gt.0).and.(kndx.le.250)) ucalc(kndx) = bb
c         write(6,17) jndx,aa,kndx,bb
c  17   format(' in USR9, set u(',i3,')=',f12.6,'  u(',i3,')=',f12.6)
         return
      endif
c
c         a valid user location is required for job.ne.0
c
      if((nusr.lt.1).or.(nusr.gt.250)) then
         write(jof,117) nusr
 117     format(' USR9 ERROR:',i8,' is not a valid UCALC location')
         return
      endif
c
c        job = 1 is to store linear combinations of R-matrix elements
c
      if(job.eq.1) then
         idiv = 10
         if((jndx.gt.10).and.(jndx.lt.67)) then
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*fm(ii,jj)
         endif
         if((kndx.gt.10).and.(kndx.lt.67)) then
             ii = kndx/idiv
             jj = kndx - idiv*ii
             value = value + bb*fm(ii,jj)
         endif
      endif
c
c        job = 2 is to store linear combinations of Lie coefficients
c
      if(job.eq.2) then
         if((jndx.gt.0).and.(jndx.le.monoms)) value=value+aa*fa(jndx)
         if((kndx.gt.0).and.(kndx.le.monoms)) value=value+bb*fa(kndx)
      endif
c
c        job=3 is to store linear combinations of UCALC entries
c
      if(job.eq.3) then
         if((jndx.gt.0).and.(jndx.le.250)) value=value+aa*ucalc(jndx)
         if((kndx.gt.0).and.(kndx.le.250)) value=value+bb*ucalc(kndx)
      endif
c
c        job=4 is to store linear combinations of Taylor coefficients
c
      if(job.eq.4) then
         if((jndx.gt.low).and.(jndx.lt.683)) then
             idiv = 100
             if(jndx.lt.idiv) idiv = 10
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*tumat(jj,ii)
ctm   write(6,*) 'jndx: tumat(',jj,ii,')=',tumat(jj,ii)
         endif
         if((kndx.gt.low).and.(kndx.lt.683)) then
             idiv = 100
             if(kndx.lt.idiv) idiv = 10
             ii = kndx/idiv
             jj = kndx - idiv*ii
             value = value + bb*tumat(jj,ii)
ctm   write(6,*) 'kndx: tumat(',jj,ii,')=',tumat(jj,ii)
         endif
      endif
c
c        job=5 is to compute ratios of Lie coefficients
c
      if(job.eq.5) then
         if((jndx.gt.0).and.(jndx.le.monoms)) value=value+aa*fa(jndx)
         if((kndx.gt.0).and.(kndx.le.monoms)) then
              denom = fa(kndx) + bb
              if(denom.ne.zero) value = value/denom
         endif
      endif
c
c        job=6 is to compute ratios of Taylor coefficients
c
      if(job.eq.6) then
         if((jndx.gt.low).and.(jndx.lt.683)) then
             idiv = 100
             if(jndx.lt.idiv) idiv = 10
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*tumat(jj,ii)
         endif
         if((kndx.gt.low).and.(kndx.lt.683)) then
             idiv = 100
             if(kndx.lt.idiv) idiv = 10
             ii = kndx/idiv
             jj = kndx - idiv*ii
             denom = tumat(jj,ii) + bb
             if(denom.ne.zero) value = value/denom
         endif
      endif
c
c        job=7 is to compute ratios of ucalc entries
c
      if(job.eq.7) then
         if((jndx.gt.0).and.(jndx.le.250)) value=value+aa*ucalc(jndx)
         if((kndx.gt.0).and.(kndx.le.250)) then
             denom = ucalc(kndx) + bb
             if(denom.ne.zero) value = value/denom
         endif
      endif
c
c        job = 8 is to store ratios of R-matrix elements
c
      if(job.eq.8) then
         idiv = 10
         if((jndx.gt.10).and.(jndx.lt.67)) then
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*fm(ii,jj)
         endif
         if((kndx.gt.10).and.(kndx.lt.67)) then
             ii = kndx/idiv
             jj = kndx - idiv*ii
             denom = fm(ii,jj) + bb
             if(denom.ne.zero) value = value/denom
         endif
      endif
c
c        job=9 is to store products of Lie coefficients
c
      if(job.eq.9) then
         if((jndx.gt.0).and.(jndx.le.monoms)) value=value+aa*fa(jndx)
         if((kndx.gt.0).and.(kndx.le.monoms)) then
              value = value*(fa(kndx) + bb)
         endif
      endif
c
c        job=10 is to store products of Taylor coefficients
c
      if(job.eq.10) then
         if((jndx.gt.low).and.(jndx.lt.683)) then
             idiv = 100
             if(jndx.lt.idiv) idiv = 10
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*tumat(jj,ii)
         endif
         if((kndx.gt.low).and.(kndx.lt.683)) then
             idiv = 100
             if(kndx.lt.idiv) idiv = 10
             ii = kndx/idiv
             jj = kndx - idiv*ii
             value = value*(tumat(jj,ii) + bb)
         endif
      endif
c
c        job=11 is to store products of ucalc entries
c
      if(job.eq.11) then
         if((jndx.gt.0).and.(jndx.le.250)) value=value+aa*ucalc(jndx)
         if((kndx.gt.0).and.(kndx.le.250)) then
             value = value*(ucalc(kndx) + bb)
         endif
      endif
c
c        job = 12 is to store products of R-matrix elements
c
      if(job.eq.12) then
         idiv = 10
         if((jndx.gt.10).and.(jndx.lt.67)) then
             ii = jndx/idiv
             jj = jndx - idiv*ii
             value = value + aa*fm(ii,jj)
         endif
         if((kndx.gt.10).and.(kndx.lt.67)) then
             ii = kndx/idiv
             jj = kndx - idiv*ii
             value = value*(fm(ii,jj) + bb)
         endif
      endif
c
c       job=13 is to compute the weighted RMS of the two ucalc entries:
c
      if(job.eq.13) then
         if((jndx.gt.0).and.(jndx.le.250)) then
             value=value+aa*ucalc(jndx)**2
         endif
         if((kndx.gt.0).and.(kndx.le.250)) then
             value=value+bb*ucalc(kndx)**2
         endif
         value = dsqrt(value)
      endif
c
c        job=14 is to compute reciprocals of selected variables (loon = kndx)
c
      if(job.eq.14) then
         if((kndx.gt.0).and.(kndx.lt.4)) then
             call getvar(kndx,nv,xv)
             if((jndx.le.nv).and.(jndx.gt.0)) then
                xpar =  xv(jndx)
                if(xpar.ne.0.0d0) value = aa/xpar
             endif
         endif
      endif
      ucalc(nusr) = value
      return
      end
