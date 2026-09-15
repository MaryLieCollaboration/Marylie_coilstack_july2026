c      New version of proc using amdii in fit  D-Day June 6 1994
c-------------------------------------------
      subroutine aim(pp)
c
c      aim allows user selection of map elements to fit
c     C. T. Mottershead /AT-6  & L. B. Schweitzer /UCB   Dec 86
c      New character controled version  July 87 & July 88
c
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'labpnt.inc'
      include 'files.inc'
      include 'aimdef.inc'
      include 'keyset.inc'
      include 'fitbuf.inc'
      include 'map.inc'
c
      dimension pp(6)
      character*1 yn
      character*8 word
      character*128 card
      dimension bufr(4)
      save
c-----!----------------------------------------------------------------!
      if(jicnt.ne.0) return
      keyrms = 11
      ier = 0
      iaim = nint(pp(1))
      infile=nint(pp(2))
      if(infile.eq.0) infile = jif
      logf = nint(pp(3))
      iquiet = nint(pp(4))
      loon = nint(pp(5))
      if((loon.lt.1).or.(loon.gt.3)) loon = 1
      ksq(loon) = iaim
      do 3 j=1,3
      lsq(j,loon) = 0
  3   continue
      isend = nint(pp(6))
      jeof = 0
  5   nf = 0
      do 10 j=1,maxf
      wts(j,loon) = 1.0
  10  continue
      go to 45
c
c       help section
c
  20  if(infile.eq.jif) then
      write(jof,*) ' '
      write(jof,*) ' Symbols are defined for the following categories:'
c      write(jof,*) '     (Enter the catagory number to get help)'
      write(jof,*) '  1. The current map, and various stored maps and wo
     $rk buffers in map format.'
      write(jof,*) '  2. Standard quantities derived from the current ma
     $p, including tunes,'
      write(jof,*) '     chromaticities, anharmonicities, dispersions, a
     $nd twiss parameters.'
      write(jof,*) '  3. Beam parameters.'
      write(jof,*) '  4. User calculated quantities.'
      write(jof,*) '  -----Select help catagory-----'
      ihelp=0
      card = ' '
ctm14 read(infile,17,end=45) card
      call nexlin(infile,jeof,card)
      if(jeof.lt.0) go to 45
      if(card.eq.' ') go to 45
      call txtnum(card,20,1,ng,bufr)
      if(ng.gt.0) ihelp = bufr(1)
      if((ihelp.gt.0).and.(ihelp.le.4)) go to (25,30,35,40), ihelp
      go to 45
  25  continue
      write(jof,*) '                 ***  Quantities in Map Format  ***'
      write(jof,*) ' '
      write(jof,*) '  The current transfer map:'
      write(jof,*) '     r(i,j) = first order (linear) matrix elements i
     $,j=1,6.'
      write(jof,*) '     f(i) = Lie polynomial coefficients, i=1,209.'
      write(jof,*) ' '
      write(jof,*) '  Stored maps in locations 1 through 5:'
      write(jof,*) '     sm1(i,j) ... sm5(i,j) = first order (linear) ma
     $trix elements i,j=1,6.'
      write(jof,*) '     sf1(i) ... sf5(i) = Lie polynomial coefficients
     $, i=1,209.'
      write(jof,*) ' '
      write(jof,*) '  Result arrays in buffers 1 through 5:'
      write(jof,*) '     bm1(i,j) ... bm5(i,j) = first order (linear) ma
     $trix elements i,j=1,6.'
      write(jof,*) '     bf1(i) ... bf5(i) = Lie polynomial coefficients
     $, i=1,209.'
      write(jof,*) '   Taylor Map: tm(i,j), i=1,6 is component;  j=1,83
     * is Giorgelli index'
      go to 45
  30  continue
      write(jof,*) '      ***  Standard quantities computed by COD, TASM
     $, and TADM ***'
      write(jof,*) ' '
      write(jof,*) ' Tunes:'
      write(jof,*) '    tx, ty = horizontal and vertical tunes.'
      write(jof,*) '        ts = synchrotron (temporal) tune for a dynam
     $ic map,'
      write(jof,*) '             or temporal dispersion (eta) for a stat
     $ic map.'
      write(jof,*) ' '
      write(jof,*) ' Chromaticities:'
      write(jof,*) '    cx, cy = horizontal and vertical 1st order chrom
     $aticities.'
      write(jof,*) '    qx, qy = horizontal and vertical 2nd order (quad
     $ratic) chromaticities.'
      write(jof,*) ' '
      write(jof,*) ' Anharmonicities:'
      write(jof,*) '    hh, vv, tt, hv, ht, vt = second order anharmonic
     $ities.'
      write(jof,*) ' '
      write(jof,*) ' Dispersions:'
      write(jof,*) '    dz1, dz2, dz3, dz4 = first order dispersions for
     $ the'
      write(jof,*) '                         transverse phase space coor
     $dinates.'
      write(jof,*) ' '
      write(jof,*) ' Twiss parameters (zeroth order):'
      write(jof,*) '    ax,bx,gx = horizontal alpha, beta, gamma.'
      write(jof,*) '    ay,by,gy = vertical alpha, beta, gamma.'
      write(jof,*) '    at,bt,gt = temporal alpha, beta, gamma.'
      go to 45
  35  continue
      write(jof,*) '      ***  Beam Parameters  ***'
      write(jof,*) ' '
      write(jof,*) ' Particle coordinates from tracking:'
c-----!----------------------------------------------------------------!
      write(jof,*) '    z(i,j) = jth component of ith particle in',
     * ' tracking buffer.'
      write(jof,*) ' '
      write(jof,*) ' Moments (from AMAP):'
      write(jof,*) '    s(i,j) = 2nd moments of the beam <zi,zj>.'
c-----!----------------------------------------------------------------!
      write(jof,*) '    bf1(i) = moments of the beam in standard Lie mon
     *omial sequence.'
      write(jof,*) ' '
      write(jof,*) ' Emittances (from AMAP):'
      write(jof,*) '    ex,ey,et = rms emittances for the x, y, and t pl
     $anes.'
      write(jof,*) '    wx,wy,wt = eigen emittances for coupled systems.
     $'
      go to 45
  40  continue
c-----!----------------------------------------------------------------!
      write(jof,*) '   ***  User defined quantities computed by USER',
     $ ' and MERIT routines  ***'
      write(jof,*) ' '
      write(jof,*) '    u(i) = ucalc(i) for i=1 to 250, stored in common
     $/usrdat/'
      write(jof,*) ' '
      write(jof,*) '    umi = val(i) for i=1 to 5, computed by user ',
     * 'merit functions','          MRT1 through MRT5 and stored in',
     * ' common/merit/'
      endif
c
c           prompt for and read input line
c
  45  if(infile.eq.jif) then
      write(jof,*) ' '
      if(iaim.eq.1) write(jof,*) ' Select quantities by entering their s
     $ymbols'
      if(iaim.eq.2) write(jof,*) ' Enter aim(s) in the form:  symbol = t
     $arget value'
      if(iaim.eq.3) write(jof,*) ' Enter aim(s) in the form:  symbol = t
     $arget value, weight'
      write(jof,*) ' (Enter a # sign when finished, or type help to revi
     $ew the defined symbols)'
      endif
  50  card = ' '
ctm14 15  read(infile,17,end=200) card(2: )
      call nexlin(infile,jeof,card)
      if(jeof.lt.0) go to 200
  17  format(a)
      if(card.eq.' ') go to 50
      call low(card)
      ifin = index(card,'#')
c
c         scan input line for keywords and numeric values
c
      do 150 ku=1,nkeys
      nask = kask(ku)
      nget = 0
      ntgt = nask
      if(iaim.eq.1) nask = nask-1
      if(iaim.eq.3) nask = nask+1
      ju = ku
      if(ju.gt.13) ju = 13
      jmax = maxj(ju)
  60  call keynum(card,key(ku),nask,nget,nloc,bufr)
      if(nloc.eq.0) go to 150
c
c     found key(ku) ; check and interpret input
c
      if(nget.ne.nask) go to 60
c
c      scalar variables
c
      idu = 0
      jdu = 0
      nsu = 0
      if(ku.gt.13) go to 72
c
c      stored maps
c
      nj = 1
      if(ku.le.4) then
         nsu = bufr(1)
         nj = 2
      endif
      idu = bufr(nj)
      imax = jmax
      if (ku.eq.7) imax = 999
      if((idu.lt.1).or.(idu.gt.imax)) go to 60
      if(ku.lt.9) then
         jdu = bufr(nj+1)
         if((jdu.lt.1).or.(jdu.gt.jmax)) go to 60
      endif
c
c     accept and report the entry
c
  72  if(nf.ge.maxa) then
          write(jof,*) maxa,' = maximum aims allowed'
          write(jodf,*) maxa,' = maximum aims allowed'
          go to 200
      endif
c
c       trap and flag rmserr requests
c
      if(ku.eq.keyrms) then
          lsq(idu,loon) = 1
          write(jof,74) idu
  74      format(' RMS',i1,' enabled')
          go to 60
      endif
      nf=nf+1
      if(iaim.gt.1) target(nf,loon) = bufr(ntgt)
      if(iaim.eq.3) wts(nf,loon) = bufr(nget)
      kyf(nf,loon) = ku
      nsf(nf,loon) = nsu
      idf(nf,loon) = idu
      jdf(nf,loon) = jdu
ctm   write(6,*) ' AIM got:',ku,nsu,'=ku,nsu',idu,jdu,'=idu,jdu'
c
c     echo accepted aim and see if there are any more
c
      if(ku.gt.12) then
         write(word,42) key(ku)
  42     format(6x,a2)
         go to 100
      endif
      go to (75,75,80,80,85,85,88,78,90,90,95,95,95), ku
  75  write(word,76) key(ku), nsu, idu, jdu
  76  format(a2,i1,'(',i1,',',i1,')')
      go to 100
  78  write(word,79) key(ku), idu, jdu
  79  format(a2,'(',i1,',',i2,')')
      go to 100
  80  write(word,81) key(ku), nsu, idu
  81  format(a2,i1,'(',i3,')')
      go to 100
  85  write(word,86) key(ku), idu, jdu
  86  format(2x,a2,i1,',',i1,')')
      go to 100
  88  write(word,89) key(ku), idu, jdu
  89  format(a2,i3,',',i1,')')
      go to 100
  90  write(word,91) key(ku), idu
  91  format(2x,a2,i3,')')
      go to 100
  95  write(word,96) key(ku), idu
  96  format(5x,a2,i1)
 100  write(jof,101) nf, word
cccccccccccc                  target(nf,loon)
 101  format('AIM:  accept',i3,':    ',a)
ccccccccccccc                     ,' = ',1pg22.14)
      qname(nf,loon) = word
      go to 60
 150  continue
c
c       check for pleas for help
c
      if(index(card,'help').ne.0) go to 20
      if(index(card,'HELP').ne.0) go to 20
      if(ifin.eq.0) go to 50
c
c      End of input, echo variables and target values
c
 200  continue
      nsq(loon) = nf
      if((isend.eq.1).or.(isend.eq.3)) call wrtaim(jof,iaim,loon)
      if((isend.eq.2).or.(isend.eq.3)) call wrtaim(jodf,iaim,loon)
c
      if(infile.eq.jif) then
         write(jof,*) ' '
         write(jof,*) ' Do you wish to start over? (y/n) <n>:'
         yn='n'
         read(jif,177) yn
 177     format(a)
         if ((yn .eq.'y').or.(yn .eq.'Y')) go to 5
      endif
c
c        write log file if requested
c
      if(logf.gt.0) then
      call wrtaim(logf,4,loon)
      write(logf,*) ' #'
      endif
      return
      end
c
ccccccccccccccccccccccccccccc   AMDII  ccccccccccccccccccccccccc
c
      subroutine amdii(iter,nv,fv,xv,ef,em,step,info)
c
c     Adaptive MultiDimensional Inverse Interpolation algorithm for
c     solving simultaneous non-linear equations.
c
c     C. T. Mottershead  LANL AT-3   Jan 93
c     email: motters@atdiv.lanl.gov, phone (505) 667-9730
c
c    Parameters:
c     iter = external iteration counter
c     for iter = 0, initialize routine by passing in:
c        fv(i) = target(i), i=1,nv
c        em = errtol = quitting tolerance for Max(fv(i)-target(i)) on
c             normal running calls with iter>nv
c        step = delta = feeler step size for first nv iterations.
c        info = maxcut = maximum allowed cuts in the reach from initial
c               point to final target. Reach is defined as the fraction
c               of the distance between the best point found so far and
c               the ultimate target. Full reach = 1.0, in which case we
c               try to step xv(i) to fit fv(i) = target(i). If this fail
c               and maxcut>0, we cut the reach in half. Thus maxcut=N
c               means attempts as short as 2**-N of the original distanc
c               are allowed. If the iteration succeeds twice with a reac
c               less than 1.0, the reach is doubled (but never exceeds 1
c       iter > 0: normal running. ic = internal iteration counter. Norma
c               self incrementing on each call, but may be set back to 0
c               cause new orthogonal feeler steps if the matrix becomes
c               degenerate.
c   On call:
c       nv = number of variables (= dimension of fit, maximum of 40)
c       fv(i), i=1,nv = computed or measured function values
c              at the current point xv(i), i=1,nv.
c   On return:
c       xv(i), i=1,nv = recommended new value for the independent
c            variables. The function values fv(i), i=1,nv at this new
c            point should be acquired and passed in on the next iteratio
c       em = Max(fv(i)-target(i)) = error measure of new incoming point.
c            Used internally for branching and quitting decisions.
c       step = length of step = cartesian distance between incoming old
c              xv(i) and returned new xv(i), i=1,nv
c       info = control flag. Quit if info > 0.
c            = -N for normal return with reach cut by 2**N. Keep going.
c            =  0 for normal return at full reach toward final target.
c                 Keep going.
c            =  1 Converged: Quit because new point is below error
c                 tolerance: em<errtol
c            =  2 Quit because error failed to decrease at full reach wh
c                 below ANTS tolerance threshold. In this case the xv(i)
c                 returned will be the best found so far.
c            =  3 Quit because step < slo at full reach.
c                  (But slo=0.0 in this version, so this should never ha
c            =  4 Quit because the new point is the worst yet, and no mo
c                 reach cuts are allowed, i.e. we hit MAXCUT without win
c            =  5 Abort because of complete degeneracy of points, making
c                 internal matrix = 0. This version should not allow thi
c                 to happen.
c---------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
      include 'amdiip.inc'
      dimension xv(*),fv(*)
ctm   save stpmin, iwin, ymax, ymin
      save
c-----------------------------------------------------------
c
c         dimension check
c
      if(nv.gt.maxv) then
        write(jof,*), ' AMDII error: ',nv,'=nv > maxv =',maxv
        info = 6
        return
      endif
c
c       initialization:  store targets and flags for iter = 0
c
      if(iter.le.0) then
         do 5 i = 1,nv
            fultgt(i) = fv(i)
            partgt(i) = fv(i)
   5     continue
         ncut = 0
         ic = 0
         level = 0
         reach = 1.0
         detol = 1.0e-9
         errtol = ef
         auxtol = em
         antmin = 1.0e-6
c        antol = 1.e4*errtol
         antol = auxtol
         if(antol.le.0.0) antol = antmin
         delta = step
         stpmin = 0.01*delta
         slo = 0.0
         maxcut = info
c   undocumented function limits test
         ymax = 1.0e+38
         ymin = -ymax
         if(info.lt.0) then
            maxcut = -info
            ymax = xv(1)
            ymin = xv(2)
         endif
         return
      endif
c
c        running for iter>0
c
      info = -ncut
      last = nv+1
      ic = ic + 1
c
c       compute error measure of new point and test for convergence
c
      em = difmax(nv,fv,partgt)
      ef = difmax(nv,fv,fultgt)
      if(ef.le.errtol) then
         info = 1
         step = 0.0
         return
      endif
c
c     save new xv,fv and em and increment xv on first  nv  calls:
c
  10  continue
      if(ic.le.nv) then
         do 20 i = 1,nv
            ff(i,ic) = fv(i)
            xx(i,ic) = xv(i)
  20     continue
         err(ic) = em
         do 30 i = 1,nv
           xv(i) = xx(i,1)
  30     continue
         step = delta*xv(ic)
         if(abs(step).lt.stpmin) step=stpmin
         xv(ic) = xv(ic) + step
         return
      endif
c
c      at ic = nv+1, we have just enough to take a MDII step, so save
c      the incoming point and do it:
c
      if(ic.eq.last) then
         isav = last
         go to 100
      endif
c
c      converged enough to extend reach for next MDII step
c
      if((ncut.gt.0).and.(em.lt.antol)) then
         iwin = 1-iwin
         if(iwin.eq.0) ncut = ncut - 1
         call mdicut(nv)
         level = -1
         call mdipik(nv,ibest,iworst,best,worst)
         antol = best*antmin
         if(antol.lt.antmin) antol = antmin
         em = difmax(nv,fv,partgt)
c        write(jof,*), ' New em =',em,'    antol=',antol
         info = -ncut
         go to 60
      endif
c----------------------------------------------------------------
c     When ic > nv+1, search for best and worst previous points:
c
      call mdipik(nv,ibest,iworst,best,worst)
      antol = best*antmin
      if(antol.lt.antmin) antol = antmin
c
c     return to best case and quit if error fails to decrease at full
c     reach while below antol
c
      if((ncut.eq.0).and.(em.ge.best).and.(em.lt.antol)) then
        info = 2
        do 40 i = 1,nv
          xv(i) = xx(i,ibest)
          fv(i) = ff(i,ibest)
  40    continue
        em = difmax(nv,fv,partgt)
        ef = difmax(nv,fv,fultgt)
        return
      endif
c
c   New point is the worst yet. Cut reach if allowed, otherwise restore
c   best available point and quit:
c
      if(em.ge.worst) then
         if(ncut.lt.maxcut) then
            ncut = ncut + 1
            call mdicut(nv)
            level = 1
            call mdipik(nv,ibest,iworst,best,worst)
            antol = best*antmin
            if(antol.lt.antmin) antol = antmin
            em = difmax(nv,fv,partgt)
c        write(jof,*), ' New em =',em,'    antol=',antol
            info = -ncut
            go to 60
         else
            info = 4
            do 50 i = 1,nv
            xv(i) = xx(i,ibest)
  50        continue
            return
         endif
      endif
c--------------------------------------------------------------------
c
c     arrange data pool for MDII step
c
  60  continue
      isav = 0
      if(em.lt.best) isav = last
      if((em.ge.best).and.(em.lt.worst)) isav = iworst
c
c     Overwrite the previous worst case with the last one, in preparatio
c     saving the new one in last place.
c
      if(isav.eq.last) then
         do 70 i = 1,nv
         ff(i,iworst) = ff(i,last)
         xx(i,iworst) = xx(i,last)
  70     continue
         err(iworst) = err(last)
      endif
c
c    if best is not last, swap it there
c
      if((isav.ne.last).and.(ibest.ne.last)) then
      write(jof,*), ' Swapping',ibest,'=ibest  for last'
         do 90 i=1,nv
         temp = ff(i,last)
         ff(i,last) = ff(i,ibest)
         ff(i,ibest) = temp
         temp = xx(i,last)
         xx(i,last) = xx(i,ibest)
         xx(i,ibest) = temp
  90     continue
         temp = err(last)
         err(last) = err(ibest)
         err(ibest) = temp
      endif
c
c     save new xv,fv and em in column isav of xx,ff, and err arrays
c
 100  continue
      if(isav.gt.0) then
         do 120 i = 1,nv
         ff(i,isav) = fv(i)
         xx(i,isav) = xv(i)
 120     continue
         err(isav) = em
      endif
c-------------------------------------------------------------------
c     compute xv using the amdii algorithm for ic > nv:
c
 200  continue
c
c     compute ga matrix for linear solver, with the best rhs vector
c     as the nv+1st column.
c
      gav = 0.0
      do 230 i = 1,nv
        kolum = nv*(i-1)
        do 210 j = 1,nv
          ga(j+kolum) = ff(j,i)-ff(j,last)
          gav = gav + abs(ga(j+kolum))
 210    continue
 230  continue
c
c      check for null matrix (serious error)
c
      if(gav.eq.0.0) then
        info = 5
        return
      endif
c
c     ok, ga matrix not 0, so normalize to make mean element = 1
c     (to prevent underflow or overflow)
c
      nvsq = nv**2
      gav = gav/float(nvsq)
      gfac = 1.0/gav
      do 260 i = 1,nv
      kolum = nv*(i-1)
      do 250 j = 1,nv
      ga(j+kolum) = gfac*ga(j+kolum)
 250  continue
      ga(i+nvsq) = gfac*(ff(i,last) - partgt(i))
 260  continue
c
c     Solve for correction vector yy:
c
 265  call solveh(ga,nv,1,det)
c
c     redo feeler steps if determinant is bad
c
c      if(abs(det).lt.detol) then
c         write(jof,*), ' Restart from best point because determinant = ',det
c         do 270 i = 1,nv
c            xv(i) = xx(i,ibest)
c            fv(i) = ff(i,ibest)
c 270     continue
c         em = difmax(nv,fv,partgt)
c         ef = difmax(nv,fv,fultgt)
cc         delta = 2.0*delta
c         ic = 1
c         go to 10
c      endif
c
c       extract yy from solveh return
c
      ytot = 0.0
      do 280 j = 1,nv
      yy(j) = ga(j+nvsq)
      ytot = ytot + yy(j)
 280  continue
c
c     Set new parameter vector for next round:
c
      do 340 i = 1,nv
        xbar(i) = 0.0
        do 320 j = 1,nv
          xbar(i) = xbar(i) + yy(j)*xx(i,j)
 320    continue
 340  continue
c
c       add the step to xn+1 to get xv
c
      step = 0.0
 400  do 450 i = 1,nv
      yy(i) = ytot*xx(i,last) - xbar(i)
      if(yy(i).gt.ymax) yy(i) = ymax
      if(yy(i).lt.ymin) yy(i) = ymin
      step = step + yy(i)**2
      xv(i) = xx(i,last) + yy(i)
 450  continue
      step = sqrt(step)
      if((step.lt.slo).and.(ncut.eq.0)) info = 3
      return
      end
c**************************************************************
      subroutine bip(p)
c  subroutine for beginning an inner procedure
c  Written by Alex Dragt on 8-8-88.
      include 'impli.inc'
      dimension p(6)
      include 'labpnt.inc'
      nitms=nint(p(1))
      jbipnt=lp
c      write(6,*) 'lp=',lp
      jicnt = 0
      return
      end
c
************************************************************************
c
      subroutine bop(p)
c  subroutine for beginning an outer procedure
c  Written by Alex Dragt on 8-8-88.
      include 'impli.inc'
      dimension p(6)
      include 'labpnt.inc'
      notms=nint(p(1))
      jbopnt=lp
c      write(6,*) 'lp=',lp
      jocnt = 0
      return
      end
c
************************************************************************
c
      subroutine chekit(loon,jiter,maxit)
c
c       returns current loop counter (jiter), and
c       maximum number of iterations (maxit) for specified loop.
c       jiter = current loop counter
c       maxit = maximum number of iterations
c       loon = loop number = 1 for inner loop
c                          = 2 for outer loop
c
c      T. Mottershead LANL AT-3  7 Oct 91
c----------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
      include 'labpnt.inc'
      if(loon.eq.2) then
         jiter = jocnt
         maxit = notms
      else
         jiter = jicnt
         maxit = nitms
      endif
      iquiet = 1
      return
      end
c
c********************************************************************
c
      subroutine cps(prms,p,ipset)
c  subroutine for capturing a parameter set
c  Written by Alex Dragt, 28 July 1988
c
      include 'impli.inc'
      include 'parset.inc'
      include 'psflag.inc'
c
      dimension prms(6),p(6)
c
c  test to see that ipset is in range
c
      if(ipset.le.0 .or. ipset.gt.maxpst) then
c print some message
      return
      endif
c
c  capture a parameter set and flag it
c
      if(icapt(ipset).eq.0) then
      do 10 i=1,6
   10 prms(i)=pst(i,ipset)
      icapt(ipset)=1
      return
      endif
c
c  transfer p to a captured parameter set
c
      do 20 i=1,6
   20 pst(i,ipset)=p(i)
      return
      end
c
***********************************************************************
c
      subroutine dapt(p)
c  This subroutine finds the dynamic aperture in one dimension.
c
      include 'impli.inc'
      include 'param.inc'
      include 'rays.inc'
      include 'parset.inc'
      include 'usrdat.inc'
c
c  calling array
      dimension p(6)
c
c  saved quantities
      save
ctm   save ient, gut, bad
c
c  set up control parameters
      ipset = nint(p(2))
      prec=p(3)
      gutic=0.
      badic=1.e5
c
c  see if this is initial entry into this subroutine
      if(ient .eq. 0) then
      gut=gutic
      bad=badic
      ient=1
      endif
c
c  examine result of tracking, and proceed accordingly
c
c  procedure if particle not lost
      if(nlost .eq. 0) then
      gut=pst(1,ipset)
c
      if(bad .eq. badic) then
      pst(1,ipset)=2.*gut
      endif
c
      if(bad .ne. badic) then
      pst(1,ipset)=gut+(bad-gut)/2.
      endif
c
      acc=(bad-gut)/bad
      write(6,*) gut, pst(1,ipset), bad, acc
      return
      endif
c
c  procedure if particle lost
      if(nlost .ne. 0) then
      bad=pst(1,ipset)
c
      if(gut .eq. gutic) then
      pst(1,ipset)=bad/2.
      endif
c
      if(gut .ne. gutic) then
      pst(1,ipset)=gut+(bad-gut)/2.
      endif
c
      acc=(bad-gut)/bad
      write(6,*) gut, pst(1,ipset), bad, acc
      endif
c
      return
      end
c
**********************************************************************
c
      double precision function difmax(n,aa,bb)
      implicit double precision (a-h,o-z)
      dimension aa(*), bb(*)
c
c         evaluate error of incoming functions
c
      difmax = 0.0
      do 18  i=1,n
        del = abs(aa(i)-bb(i))
        difmax = dmax1(del,difmax)
  18  continue
      return
      end
c****************************************************************
      subroutine fit(pp)
c
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'status.inc'
      include 'mldex.inc'
c     include 'fitbuf.inc'
      include 'aimdef.inc'
c
      dimension pp(6)
c
      integer iter, jiter
      character*40 quit(6)
      parameter (maxf=100)
      dimension  xv(maxf), aims(maxf), lout(2)
      save
ctm   save kut, info, em, ef, step, iter, jiter, kfit, nrpt
ctm   write(6,*) 'FIT entry params:',(pp(ip), ip=1,6)
c
c         reason for quiting messsages
c
      quit(1) = 'Converged: error < tolerance'
      quit(2) = 'Best effort: Hit bottom at full reach'
      quit(3) = 'step size below limit'
      quit(4) = '**FAILURE: error increase at full cut'
      quit(5) = '**NULL MATRIX: Nothing is Happening!!! **'
      quit(6) = '**AMDII OVERFLOW: Too many variables **'
c
c     check and reset error flag from eig4 or eig6 to ok value
      if (imbad .ne. 0) imbad=0
c
c     check iteration counter for specified loop
c
      loon =1
      nfit  = nint(pp(1))
      if(nfit.lt.0) then
          nfit = -nfit
          loon = 2
      endif
      call chekit(loon,jiter,maxit)
ctm   write(6,117),loon,jiter,maxit,info,kfit
c117  format('Iteration check: loon,jiter,maxit,info,kfit=',5i5)
c
c       pick up parameters only on first pass (jiter=0)
c
      if(jiter.eq.0) then
        write(6,*) 'Init FIT: set',kfit,'=kfit when',jiter,'=jiter'
        kfit = 0 ! one more time flag
        auxtol = pp(2)
        errtol = pp(3)
        delta = pp(4)
        mode = nint(pp(5))
        nrpt = 0
        if(mode.lt.0) nrpt = -mode
c        write(jof,*), mode,'=MODE',nrpt,'=NRPT'
        isend = nint(pp(6))
c
c       check for x damping
c
        relax = 0.5
        nrep = 0
          if(nfit.lt.0) then
             nrep = -nfit
             if(auxtol.ne.0.0) relax = auxtol
          endif
        ibsave = ibrief
        ibrief = -7
ctm   write(6,*) '>> set kfit = 0  on first pass'
c
c           load fixed data into amdii
c
        nf = nsq(loon)
        do 20 j=1,nf
        aims(j) = target(j,loon)
        write(jof,*), ' target',j,' =',aims(j)
  20    continue
        call amdii(jiter,nf,aims,xv,errtol,auxtol,delta,nfit)
      endif
c
c     set up output levels and routing
c
      lout(1) = jof
      lout(2) = jodf
      ja = 1
      if(iabs(isend).eq.2) ja = 2
      jb = 1
      if(iabs(isend).gt.1) jb = 2
      if(isend.eq.0) ja=3
c
c      check the last go around flag to restore final values everywhere
c
ctm   write(6,*) kfit,'=kfit with',info,'=info (stop if>0) at go around'
      if(kfit.eq.1) go to 700
c--------------------------------------------------------------
c     Have another go at it: Collect the current x-parameter vector,
c     and computed function values:
c
      call getaim(loon,nf,aims)
      call getvar(loon,nx,xv)
      if(nx.ne.nf) then
         write(jof,31) nx,nf,loon
         write(jof,31) nx,nf,loon
  31     format(' **FATAL FIT ERROR:',i3,' variables with',i3,
     *' aims. Abort loop',i3)
         go to 1000
      endif
c
c     compute the next estimate of the x-parameters, and put them
c     back in the parameter buffer:
c
      iter=jiter+1
ctm   write(6,*) iter,'=iter before amdii',kfit,'=kfit',info,'=info'
      call amdii(iter,nf,aims,xv,ef,em,step,info)
ctm   write(6,*) iter,'=iter  AFTER amdii',kfit,'=kfit',info,'=info'
      call putvar(loon,xv)
c
c      check the MDII error flag. If quit, do the reset pass.
c
ctm   write(6,*) info,'=info',kfit,'=kfit. If info>0, set kfit=1'
      if(info.gt.0) then
         kfit = 1
         return
      endif
c
c          iteration reports
c
      if(ja.gt.jb) go to 100
      do 60 j = ja,jb
      lun = lout(j)
c
c          optional reports every nrpt iterations
c
      if((nrpt.gt.0).and.(mod(iter,nrpt).eq.0)) then
       write(6,*) '>> Report for',nrpt,'=nrpt when',iter,'=iter'
        write(lun,53) (aims(k), k=1,nf)
  53    format(' AIMS: ',4(1pg18.9))
        write(lun,54) (xv(k), k=1,nx)
  54    format(' XNEW: ',4(1pg18.9))
      endif
      kut = 2**(-info)
      write(lun,57) iter,ef,step,em,kut
  57  format(' Iter',i4,' Error=',1pe11.4,',  Step=',1pe11.4,
     * ',  SubErr=',1pe11.4,' @cut=',i5)
      write(lun,58)
  58  format(1x,55('-'))
  60  continue
 100  continue
      return
c
c     Converged at full reach, or failed. Display final results.
c
 700  continue
       write(6,*) 'Final fit report for iteration',iter
      if(isend.eq.0) go to 1000
      ja = 1
      jb = 2
      do 900 j=ja,jb
      lun = lout(j)
      if((info.lt.1).or.(info.gt.6)) then
         info = 6
         write(6,*) 'info out of bounds, set to 6 for quit message'
      endif
      write(6,*) lun,' = lun for final fit report'
      write(lun,*) iter,'= iter final value.'
      write(lun,710) iter,info,quit(info),kut
  710 format(' Quit on iteration',i6,' for reason',i2,': ',a,/
     *' Final values with reach =',i5,' are:')
      call wrtaim(lun,2,loon)
c
c     Display new values for parameters.
c
      call wrtvar(loon,lun)
c
c     Display maximum percent error.
c
      write(lun,*)' '
      write(lun,840) em
 840   format(2x,'Maximum error is    ',1pg13.6)
      write(lun,841) errtol
 841   format(2x,'Maximum allowed was ',1pg13.6)
 900  continue
 1000 ibrief = ibsave
      call stopit(loon)
      return
      end
c
**********************************************************************
c
      subroutine fps(ipset)
c
c  this subroutine frees a parameter set so that it can be recaptured
c  if desired
c  Written by Alex Dragt, 28 July 1988
c
      include 'impli.inc'
      include 'parset.inc'
      include 'psflag.inc'
c
c  test to see that ipset is within range
c
      if(ipset.le.0 .or. ipset.gt.maxpst) then
c write error message
      return
      endif
c
c  free the indicated parameter set
c
      icapt(ipset)=0
c
      return
      end
c
*************************  GETAIM   ********************************
c
      subroutine getaim(loon,naim,aims)
c
c      getaim returns the number of aims selected (nf,loon), and
c      their current values (aims(j), j=1,nf) from aim block loon
c
c      C. T. Mottershead  LANL   AT-3    16 Aug 90
c           modified  29 May 91 to double aim blocks, selected by loon.
c--------------------------------------------------
      include 'impli.inc'
      include 'aimdef.inc'
      dimension aims(*)
c
      naim = nsq(loon)
      do 20 nu = 1, naim
      ku = kyf(nu,loon)
      nsu = nsf(nu,loon)
      idu = idf(nu,loon)
      jdu = jdf(nu,loon)
      aims(nu) = valuat(ku,nsu,idu,jdu)
 20   continue
      return
      end
c
c*************************   GETVAR    ***************************
c
      subroutine getvar(loon,nx,xv)
c
c     collect the current x-parameter vector
c     mm(i,kvs) is the i-th marylie menu item selected in vary
c     idx(i,kvs) is the index (1 to 6) of the i-th variable parameter
c       T. Mottershead   LANL  AT-3   23 July 90
c-------------------------------------------------------
      include 'impli.inc'
      include 'elmnts.inc'
      include 'xvary.inc'
      dimension xv(*)
      kvs = loon
      if((kvs.lt.1).or.(kvs.gt.3)) kvs = 1
      nx = nva(kvs)
      do 210 i = 1,nx
      mnu = mm(i,kvs)
      xv(i) = pmenu(idx(i,kvs)+mpp(mnu))
 210  continue
      return
      end
c
c***********************************  GRAD  ***************************
c
      subroutine grad(pp)
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'aimdef.inc'
      include 'xvary.inc'
c
      parameter (maxg=50)
      common/gradat/nx,ny,dxx,dyy,xcen(maxg),ycen(maxg)
      dimension pp(6)
      dimension yval(maxg), xval(maxg)
      logical ltty, lfile
      character*12 vname
      character*8 yname(maxg)
c
c       check loop and option flag
c
      loon = 1
      job = nint(pp(1))
      if(job.lt.0) then
         job = -job
         loon = 2
      endif
      call chekit(loon,iter,maxit)
c
c     initialization for iter = 0 (ic=1)
c
      if(iter.eq.0) then
         isides = nint(pp(2))
         lun = nint(pp(3))
         delta = pp(4)
         scale = pp(5)
         isend = nint(pp(6))
         ltty = .false.
         lfile = .false.
         if((isend.eq.1).or.(isend.eq.3)) ltty = .true.
         if((isend.eq.2).or.(isend.eq.3)) lfile = .true.
c
c     collect the central x-parameter vector
c
         call getvar(loon,nx,xcen)
c
c      collect central y-values
c
         ny = 0
         if((job.eq.1).or.(job.eq.3)) then
            call getaim(loon,ny,ycen)
            do 10 k=1,ny
            yname(k) = qname(k,loon)
  10        continue
         endif
c
c     compute any rms errors selected in current loon
c
         if(job.ge.2) then
           do 15 ku = 1,3
           if(lsq(ku,loon).ne.1) go to 15
           ny = ny + 1
           ycen(ny) = rmserr(ku)
           write(yname(ny),14) ku
  14       format(2x,'rms',i1,2x)
  15       continue
         endif
c
c      increment the first x-variable on the first call
c
         do 20 i=1,nx
         xval(i) = xcen(i)
  20     continue
         xval(1) = xcen(1) + delta
         call putvar(loon,xval)
         return
      endif
c
c       all done, jump out of labor loop
c
      if(iter.gt.nx) then
          write(jif,*) ' gradient loop is finished'
          call stopit(loon)
          return
      endif
c
c     normal running:  compute and save numerical gradients for iter=1,n
c
      call getvar(loon,nx,xval)
      ny = 0
      if((job.eq.1).or.(job.eq.3)) call getaim(loon,ny,yval)
c
c     compute any rms errors selected in current loon
c
      if(job.ge.2) then
        do 515 ku = 1,3
        if(lsq(ku,loon).ne.1) go to 515
        ny = ny + 1
        yval(ny) = rmserr(ku)
 515    continue
      endif
      delx = xval(iter) - xcen(iter)
      do 600 n = 1,ny
      dydx = (yval(n) - ycen(n))/delx
      if(dydx.eq.0.0) go to 600
      qsq = dydx/delx
      qs = 0.0
      if(qsq.gt.0.0) qs = scale*sqrt(qsq)
      vname = varnam(iter,loon)
      if(ltty) write(jof,577) n,iter,dydx,qs,yname(n),vname
      if(lfile) write(jodf,577) n,iter,dydx,qs,yname(n),vname
      if(lun.gt.0) write(lun,577) n,iter,dydx,qs,yname(n),vname
  577 format(2i5,1pg22.14,1pg16.8,2x,a,' / ',a)
  600 continue
c
c    increment parameter vector on first nx-1 calls, and restore
c    central values on last call:
c
      if(iter.lt.nx) then
         do 520 i=1,nx
         xval(i) = xcen(i)
 520     continue
         xval(iter+1) = xcen(iter+1) + delta
         call putvar(loon,xval)
      else
         call putvar(loon,xcen)
      endif
      return
      end
c
**********************************************************************
c
      subroutine keynum(text,word,nask,nget,nloc,bufr)
c
c      keynum scans the character string 'text' for one
c      keyword, followed by a list of associated numbers.
c
c     parameters:
c          text  - the input character string to be searched
c          word  - the single keyword
c          nask  - maximum length of number list
c          nget  - actual number of numbers found
c          nloc  - location in text string of keyword found
c          bufr  - output array of real numbers
c
c     T. Mottershead  Los Alamos  Aug 1985
c----------------------------------------------------------
      implicit double precision (a-h,o-z)
      character text *(*), word*(*)
      dimension bufr(1)
      max=len(text)
      nch=len(word)
      nloc=index(text,word)
      if(nloc.eq.0) go to 40
      jj=nloc+nch-1
      text(nloc:jj)=' '
      nrem=max-nloc+1
      if(nask.gt.0) call txtnum(text(nloc: ),nrem,nask,nget,bufr)
c      write(6,17) nrem,nask,nget,(bufr(k), k=1,nask)
c  17  format(' in keynum:',3i3,3(1pg14.6))
  40  return
      end
c
************************************************************************
c
      subroutine lexort(num,namlst,indx)
c
c      General character array sorting subroutine.
c        parameters -
c        num    - dimension of the arrays
c        namlst - input array of character variables to be
c                 indexed in alphabetical order. This input
c                 array is not altered.
c        indx   - output sorted index array, generated so that
c                 namlst(indx(k)), k=1,num is in alphabetical order.
c      C. T. Mottershead/ AT-3  2 Feb 88
c----------------------------------------
      dimension indx(*)
      character*(*) namlst(*)
      character*8 first
c
c     initialize the index array
c
      do 5 j = 1, num
      indx(j) = j
   5  continue
c
c        find smallest remaining value and corresponding index
c
      j=1
  10  jin = indx(j)
      first=namlst(jin)
      mv=j
c
c     find first remaining name, and location of its pointer in the
c     index array (i.e namlst(indx(mv)) is the lowest remaining name)
c
      do 30 k=j,num
      kin = indx(k)
      if(lge(namlst(kin),first)) go to 30
      first=namlst(kin)
      mv=k
   30 continue
c
c      swap the index of the lowest name (in location mv) with the
c      first index (in location j) in the remaining list
c
      isav = indx(j)
      indx(j) = indx(mv)
      indx(mv) = isav
c
c        move starting location down one and repeat if necessary
c
      j=j+1
      if(j.lt.num) go to 10
   77 return
      end
c
************************************************************************
c
      subroutine mdicut(nv)
c
c       reset partial target (PARTGT) and rejudge the pool (XX,FF) relat
c       to the new target
c
c         C. T. Mottershead  AT-3  6 July 92
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'amdiip.inc'
      include 'files.inc'
c-----------------------------------------------------------
      last = nv+1
      npower = 2**ncut
      write(jof,*), npower,'<<----Reach cut factor'
      reach = 1.0/float(npower)
      do 20 k = 1, nv
      partgt(k) = ff(k,last) + reach*(fultgt(k) - ff(k,last))
  20  continue
c
c         reevaluate errors of data pool relative to new target
c
      do 60 n=1,last
      em = 0.0
         do 40 i=1,nv
         df = abs(ff(i,n)-partgt(i))
         em = dmax1(df,em)
  40     continue
      err(n) = em
  60  continue
      return
      end
c**************************************************************
      subroutine mdipik(nv,ibest,iworst,best,worst)
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'amdiip.inc'
c-----------------------------------------------------------
      last = nv + 1
      best=1.e30
      worst=-best
      ibest=1
      iworst=1
      do 60 i = 1,last
      if(err(i).lt.best) then
         best = err(i)
         ibest=i
      endif
      if(err(i).gt.worst) then
         worst = err(i)
         iworst = i
      endif
  60  continue
c     write(6,67) ibest,best,iworst,worst
c  67  format(' MDIPIK:',i3,' is best=',1pg15.7,i8,' is worst=',1pg15.7)
      return
      end
c*************************************************************
c      subroutine mrt0
c
c   This is the least squares merit function defined by aims for loop 1
c
c     T. Mottershead, LANL  Feb 89
c-----------------------------------------------------
c      include 'impli.inc'
c      include 'merit.inc'
c
c        gather the computed function values
c
c      loon = 1
c      write(6,*) 'in mrt0'
c      fval = rmserr(loon)
c      val(0) = fval
c      return
c      end
c
*****************************************************************
c
      subroutine mss(pp)
c
c     optimization routine to least squares fit the selected aims to the
c     specified targets.
c     T. Mottershead  LANL  AT-3  24 Nov 92
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'status.inc'
      include 'mldex.inc'
      include 'aimdef.inc'
c
      parameter (maxf = 30)
      dimension pp(6), tol(4)
      dimension  aims(maxf), xv(maxf), fv(maxf), lout(2)
      save fval, nx
c
c     check and reset error flag from eig4 or eig6 to ok value
      if (imbad .ne. 0) imbad=0
c
c     check iteration counter for specified loop
c
      nrpt = 1
      loon =1
      nopt  = nint(pp(1))
      if(nopt.lt.0) then
          nopt = -nopt
          loon = 2
      endif
      call chekit(loon,jiter,maxit)
      iter = jiter + 1
      if(jiter.eq.0) then
c
c       pick up parameters only on first pass (jicnt=0)
c
        ftol = pp(2)
        gtol = pp(3)
        delta = pp(4)
        xtol = pp(5)
        tol(1) = gtol
        tol(2) = xtol
        tol(3) = ftol
        tol(4) = delta
        isend = nint(pp(6))
        info = 0
        ibsave = ibrief
        ibrief = -7
        iscale = 0
        mode = 0
c
c     set up output levels and routing
c
        lout(1) = jof
        lout(2) = jodf
        italk = 0
        if(isend.lt.0) italk = -isend
        if(isend.eq.4) italk = 2
        if(isend.eq.-4) italk = 1
        ja = 1
        if(iabs(isend).eq.2) ja = 2
        jb = 1
        if(iabs(isend).gt.1) jb = 2
      write(jof,*) ' MSS is initialized'
      endif
c
c     progress report on previous iteration;
c     compute fv and chi-squared
c
      call getaim(loon,nf,aims)
      do 40 jj = 1,nf
      fv(jj) = aims(jj) - target(jj,loon)
40    continue
      fold = fval
      fval = rmserr(loon)
      change = fval - fold
c      if(abs(change).lt.ftol) go to 700
      if(ja.gt.jb) go to 100
      do 60 j = ja,jb
      lun = lout(j)
      write(lun,*)' '
      write(lun,51) iter, fval, change
  51  format(2x,'MSS iteration',i5,' F=',1pg18.9,'  DF=',1pg14.6)
      if((nrpt.gt.0).and.(mod(iter,nrpt).eq.0)) then
        write(lun,54) (xv(k), k=1,nx)
  54    format(' XNEW: ',4(1pg18.9))
      endif
      write(lun,58)
  58  format(1x,55('-'))
  60  continue
 100  continue
c
c     Collect the current x-parameter vector.
c
      call getvar(loon,nx,xv)
c
c     compute the next estimate of the minimum,
c
      if(nopt.eq.1) then
         call nls(mode,iter,nx,xv,nf,fv,iscale,wts(1,loon),tol,info)
      else
c
c       use QSO on rms error merit function
c
        call qso(mode,iter,nx,xv,fval,iscale,tol,info)
      endif
      write(jof,*)'  optimizer return:',info,'=info'
      write(jof,*)'  XV=',(xv(ij),ij=1,nx)
c
c     scatter the new x-parameter estimate, and let the labor continue.
c
      if(info.eq.0) then
         call putvar(loon,xv)
         return
      endif
c
c        Termination. Display final results.
c
 700  continue
      if((info.lt.4).or.(info.gt.6)) call putvar(loon,xv)
      if(isend.eq.0) go to 1000
      do 900 j=ja,jb
      lun = lout(j)
      call message(info,lun)
      write(lun,710) iter
  710 format(' After ',i4,' total iterations, final values are:')
      call wrtaim(lun,2,loon)
      call wrtvar(loon,lun)
      write(lun,*)' '
      write(lun,840) change
 840  format(2x,'Final change was    ',1pg13.6)
c      write(lun,841) dftol
c 841  format(2x,'Tolerance was       ',1pg13.6)
 900  continue
 1000 call stopit(loon)
      ibrief = ibsave
      return
      end
c
c******************************************************************
c
      subroutine nexlin(lun,jeof,card)
c
c    Get next line of text from file or internal text block.
c    lun = logical unit number of file or block.
c          if lun > 0 read from disk file
c          if lun < 0 copy from internal data block
c    jeof = i/o flag returned for this read:
c         = +N ---> nexlin is returning line N of file or block.
c         =  0 ---> ERROR: no data returned. Do not use card.
c         = -N ---> Normal EOF. File had only N records, do not use card.
c     card = The character*80 line of text from the file/block.
c
c    C.T.Mottershead Jan 2014
c----------------------------------------------------------
      character*(*) card
      include 'cxdata.inc'
      data kurec /0/, kurid /0/
c
c      the normal external file read for lun > 0
c
      if(lun.gt.0) then
         read(lun,111,end=190) card
 111     format(a)
         kurec = kurec + 1
         jeof = kurec
         return
 190     jeof = -kurec
         return
      endif
c
c   error return for lun = 0
c
      if(lun.eq.0) then
         write(6,*) 'NEXLIN error:',lun,'=lun'
         jeof = 0
         return
      endif
c
c      internal block copy for lun < 0
c
      ilun = -lun
      if(ilun.eq.kurid) go to 50
c
c    otherwise hunt for first instance of ilun
c
      mm = 1
  20  continue
      if(lunex(mm).eq.ilun) go to 30
      mm = mm + 1
      if(mm.le.maxrec) go to 20
c
c         error return: lun not found
c
      write(jof,28) lun
  28  format('NEXLIN error: lun',i3,' not found')
      jeof = lun
      return
c
c      found lun, reset kurid = lun
c
  30  continue
      kurid = ilun
      kzero = mm - 1
      kurec = 0
c
c       normal internal read of next record in same lun
c
  50  continue
      kurec = kurec + 1
      jrec = kzero + kurec
      if(lunex(jrec).ne.ilun) go to 70
      card = mltext(jrec)
      return
c
c     EOF on internal block
c
  70  continue
      jeof = -kurec
      kurec = 0
      write(6,73) ilun, jeof
  73  format('NEXLIN EOF: unit',i3,' has only',i3,' records')
      return
      end
c
**********************************************************************
c
      subroutine opt(pp)
c
c     optimization routine
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'status.inc'
      include 'mldex.inc'
c                   include 'labpnt.inc'
      include 'merit.inc'
      parameter (maxf = 30)
      dimension  aims(maxf), xv(maxf), xold(maxf), lout(2)
c
      dimension pp(6), tol(4)
      save fval
c
c
c     check and reset error flag from eig4 or eig6 to ok value
      if (imbad .ne. 0) imbad=0
c
c     check iteration counter for specified loop
c
      nrpt = 1
      loon =1
      job = nint(pp(1))
      if(job.lt.0) then
          job = -job
          loon = 2
      endif
      call chekit(loon,jiter,maxit)
      iter = jiter + 1
      if(jiter.eq.0) then
c
c       pick up parameters only on first pass (jicnt=0)
c
        ftol = pp(2)
        gtol = pp(3)
        delta = pp(4)
        xtol = pp(5)
        tol(1) = gtol
        tol(2) = xtol
        tol(3) = ftol
        tol(4) = delta
        isend = nint(pp(6))
        info = 0
        ibsave = ibrief
        ibrief = -7
        iscale = 0
        mode = 0
c
c     set up output levels and routing
c
        lout(1) = jof
        lout(2) = jodf
        ja = 1
        if(iabs(isend).eq.2) ja = 2
        jb = 1
        if(iabs(isend).gt.1) jb = 2
      write(jof,*) ' OPT is initialized'
      endif
c
c     Collect the current x-parameter vector.
c
      call getvar(loon,nx,xold)
      do 20 j=1,nx
      xv(j) = xold(j)
  20  continue
c
c          collect selected merit function
c
      fold = fval
      if(job.eq.1) then
         call getaim(loon,nf,aims)
         fval = aims(1)
      else
         fval = rmserr(loon)
      endif
      change = fval - fold
c
c       use QSO on selected merit function
c
      call qso(mode,iter,nx,xv,fval,iscale,tol,info)
      dxsq = 0.0
      do 40 j=1,nx
      dxsq = (xv(j) - xold(j))**2
 40   continue
      step = sqrt(dxsq)
c
c     progress report
c
      if(ja.gt.jb) go to 600
      do 60 j = ja,jb
      lun = lout(j)
      write(lun,*)' '
      write(lun,51) iter, fval, change, step
  51  format(i5,'=iter   F=',1pg20.12,'  DF=',1pg14.6,
     * '  DX=',1pg14.6)
      if((nrpt.gt.0).and.(mod(iter,nrpt).eq.0)) then
        write(lun,54) (xv(k), k=1,nx)
  54    format(' XNEW: ',4(1pg18.10))
c        write(lun,57) (dx(k), k=1,nx)
c  57    format(' DX: ',6(1pg12.4))
      endif
      write(lun,58)
  58  format(1x,55('-'))
  60  continue
 600  continue
c
c     scatter the new x-parameter estimate, and let the labor continue.
c
      if(info.eq.0) then
         call putvar(loon,xv)
         return
      endif
c
c        Termination. Display final results.
c
 700  continue
      if((info.lt.4).or.(info.gt.6)) call putvar(loon,xv)
      if(isend.eq.0) go to 1000
      do 900 j=ja,jb
      lun = lout(j)
      call message(info,lun)
      write(lun,710) iter
  710 format(' After ',i4,' total iterations, final values are:')
      if(job.eq.1) then
         call wrtaim(lun,1,loon)
      else
         final = rmserr(loon)
         write(lun,721) loon, final
  721 format(' Final RMS error for loop',i2,' is',1pg18.10)
      endif
      call wrtvar(loon,lun)
      write(lun,*)' '
      write(lun,840) change, step
 840  format(2x,'Final change was ',1pg13.6,2x,'final step was ',
     * 1pg13.6)
      write(lun,841) ftol, xtol
 841  format(2x,'Tolerances: ftol=',1pg13.6,12x,'xtol=',1pg13.6)
 900  continue
 1000 call stopit(loon)
      ibrief = ibsave
      return
      end
c
**********************************************************************
c
      subroutine putvar(loon,xnew)
c
c     scatter the new x-parameter vector back to the parameter blocks
c     including reseting the dependent variables.
c     mm(i,kvs) is the i-th marylie menu item selected in vary
c     idx(i,kvs) is the index (1 to 6) of the i-th variable parameter
c        The first nv of these specify the independent variables being
c        adjusted by the fit or optimization routines.
c        The next nd entries specify parameters linearly dependent on
c        one of the independent variables in the form
c            base_value + slope*variable
c        idv(j,kvs) specifies which independent variable (1 to nv)
c        the j-th dependent parameter is tied to.
c        xbase(j) and xslope(j) are the relevant base value and
c        slope specified in vary.
c
c           T. Mottershead  LANL AT-3     23 July 90
c-------------------------------------------------------------
      include 'impli.inc'
      include 'elmnts.inc'
      include 'xvary.inc'
      dimension xnew(*)
c
c       simple copy back for the independent varables
c
      kvs = loon
      if(kvs.ne.2) kvs = 1
      nv = nva(kvs)
      nd = nda(kvs)
      do 20 i=1,nv
      mnu = mm(i,kvs)
      pmenu(idx(i,kvs)+mpp(mnu)) = xnew(i)
 20   continue
      if (nd.le.0) return
c
c       compute new values for linearly dependendent variables
c
      do 50 j=1,nd
      idu = idx(j+nv,kvs)
      mnu = mm(j+nv,kvs)
      depvar = xbase(j,kvs) + xslope(j,kvs)*xnew(idv(j,kvs))
      pmenu(idu+mpp(mnu)) = depvar
  50  continue
      return
      end
c
ccccccccccccccccccccccc    RMSERR    cccccccccccccccccccccccc
c
      double precision function rmserr(loon)
c
c   This function calculates the RMS error of the specified loop
c
c     T. Mottershead, LANL AT-3   8 Dec 92
c-----------------------------------------------------
      include 'impli.inc'
      include 'aimdef.inc'
c
      rmserr = 0.0
      nf = nsq(loon)
      if(nf.le.0) return
      fnum = float(nf)
      chisq = 0.0
      do 20 nu = 1, nf
      ku = kyf(nu,loon)
      nsu = nsf(nu,loon)
      idu = idf(nu,loon)
      jdu = jdf(nu,loon)
      aimiss = valuat(ku,nsu,idu,jdu) - target(nu,loon)
      chisq = chisq + wts(nu,loon)*aimiss**2
 20   continue
      rmserr = sqrt(chisq/fnum)
      return
      end
c
***************************************************************
c
      subroutine rset(pp)
c
c  this subroutine resets items in the #menu component interactively
c  made from parts of 'vary', 19 Aug 91, T. Mottershead, LANL
c
      include 'impli.inc'
      include 'files.inc'
      include 'codes.inc'
      include 'elmnts.inc'
      character*80 card
      dimension pp(6), menitm(20), indvar(20)
      logical ltty, lfile
c
c      write(6,*) ' in subroutine rset'
c-----------------------------------------------------
      job = nint(pp(1))
      infile=nint(pp(2))
      if(infile.eq.0) infile = jif
      logf = nint(pp(3))
      lun = nint(pp(4))
      isend = nint(pp(5))
      lfile = .false.
      ltty = .false.
      if((isend.eq.1).or.(isend.eq.3)) ltty = .true.
      if((isend.eq.2).or.(isend.eq.3)) lfile = .true.
      ktot = 0
      maxrs = 20
      jeof = 0
  10  if(ltty) call varlis(jof)
c
c     select variable parameters
c
  20  if(ltty) write(jof,32)
  32  format(/' To RESET #menu items, enter name and (optional)',
     * ' parameter index.',/' Type * to relist, or a # sign',
     * ' when finished.')
c
c       read a new input line
c
ctm14  40  read(infile,41,end=200) card
 40   continue
      call nexlin(infile,jeof,card)
      if(jeof.lt.0) go to 200
  41  format(a)
      if(card.eq.' ') go to 40
      call low(card)
      ifin = index(card,'#')
      list = index(card,'*')
      call vreset(infile,card,maxrs,ktot,menitm,indvar)
      if(list.ne.0) go to 10
      if(ifin.eq.0) go to 20
 200  continue
c
c     Display reset selections
c
      if(ktot.eq.0) return
      if(ltty) write(jof,203) ktot
      if(lfile) write(jodf,203) ktot
      if(lun.gt.0) write(lun,203) ktot
 203  format(/' The following',i3,' #menu item(s) have been reset:')
      if(ltty) write(jof,175)
      if(lfile) write(jodf,175)
      if(lun.gt.0) write(lun,175)
 175  format(' No.   Item   ',4x,'Type',5x,'Parameter',3x,'New value'
     * /60('-'))
      do 170 i = 1,ktot
      idu = indvar(i)
      mnu = menitm(i)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      if(ltty) write(jof, 176) i, lmnlbl(mnu), ltc(mt1,mt2), idu,
     *  pmenu(idu+mpp(mnu))
      if(lfile) write(jodf, 176) i, lmnlbl(mnu), ltc(mt1,mt2),
     * idu, pmenu(idu+mpp(mnu))
      if(lun.gt.0) write(lun, 176) i, lmnlbl(mnu), ltc(mt1,mt2),
     * idu, pmenu(idu+mpp(mnu))
 176  format(i3,4x,a8,3x,a8,i6,5x,1pg21.14)
 170  continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine scan(pp)
c
c       routine for scanning parameters
c       T. Mottershead   LANL AT-3     29-MAR-91
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
c
      parameter (maxg=50)
      common/gradat/nx,nz,dxx,dyy,xcen(maxg),ycen(maxg)
      dimension pp(6)
      dimension xval(maxg)
      logical ltty, lfile
      save
ctm   save ltty, lfile, ncyc, deltax
      loon =1
      nopt  = nint(pp(1))
      if(nopt.lt.0) then
          nopt = -nopt
          loon = 2
      endif
      call chekit(loon,iter,maxit)
c
c     initialization for iter = 0 (ic=1)
c
      if(iter.eq.0) then
        ndx = nint(pp(2))
        ndy = nint(pp(3))
        ncyc = ndx + 1
        if(ndx.lt.0) ncyc = 1 - 2*ndx
        if(ncyc.ge.maxit) ncyc = maxit - 1
        deltax = pp(4)
        deltay = pp(5)
        isend = nint(pp(6))
        ltty = .false.
        lfile = .false.
        if((isend.eq.1).or.(isend.eq.3)) ltty = .true.
        if((isend.eq.2).or.(isend.eq.3)) lfile = .true.
c
c     save the original x-parameter vector and set up loop
c
        call getvar(loon,nx,xcen)
        xx = xcen(1)
        if(ndx.lt.0) xx = xx + float(ndx)*deltax
c
c      increment the first x-variable
c
        do 20 i=1,nx
        xval(i) = xcen(i)
  20    continue
        xval(1) = xx
        call putvar(loon,xval)
        return
      endif
c
c     normal running:  print selected quantities for iter=1,nx:
c
      if(iter.gt.ncyc) go to 900
      call getvar(loon,nx,xval)
      xx = xval(1)
ctm   yy = xval(2)
      yy = 0.0d0
      if(ltty) write(jof,433) iter, xx, yy
      if(lfile) write(jodf,433) iter, xx, yy
 433  format(' ** Scan Step',i4,' xx, yy =',2(1pg15.7))
c
c    increment parameter vector on first ncyc-1  calls:
c
      if(iter.lt.ncyc) then
         do 520 i=1,nx
         xval(i) = xcen(i)
 520     continue
         xx = xx + deltax
         xval(1) = xx
         call putvar(loon,xval)
      else
c        restore original values on last pass
         if(nopt.eq.1) call putvar(loon,xcen)
      endif
      return
c
c       all done, put it back the way it was
c       and jump out of labor loop
c
 900  continue
      if(ltty) write(jof,*) '  Scan loop is finished'
      if(lfile) write(jodf,*) '  Scan loop is finished'
      call stopit(loon)
      return
      end
c
*****************************************************************
c
      subroutine solveh(augmat,nrow,nca,det)
c  Solves the linear equations
c       m*a = b
c  where m is  nrow by nrow
c        a is  nrow by nca
c        b is  nrow by nca
c  On entry :   augmat = m augmented by b
c               nca = number of columns in a or b
c  On return:   augmat = identity augmented by a
c               det = determinant of m
      double precision augmat(nrow,nrow+nca), det
c
c  Coded from a routine originally written for the HP41C calculator
c  (Dearing, p.46).  Written by Liam Healy, February 14, 1985.
c  routine and some variables renamed by Tom Mottershead, 15-Nov-88.
c----Variables----
c  nrow = dimension of matrix = total number of rows.
c  nca  = number of columns augmented
c  ncol = total number of columns = nrow + nca.
      integer nrow, nca, ncol
c
c  irow,jcol,ir,jc,roff,mr=row and column numbers, row and column indice
c  row offset, row number for finding maxc.
      integer irow,jcol,ir,jc,roff,mr
c
c  elem, mrow, hold = matrix element, its row number, and held value.
c  const = constant used in multiplication
      double precision elem,hold,const
      integer mrow
c
c----Routine----
      det=1.
      ncol=nrow+nca
      jcol=0
      do 100 irow=1,nrow
 300    jcol=jcol+1
        if (jcol.le.ncol) then
c               find max of abs of mat elts in this col, and its row num
          elem=0.
          mrow=0
          do 120 mr=irow,nrow
            if (abs(augmat(mr,jcol)).ge.abs(elem)) then
              elem=augmat(mr,jcol)
              mrow=mr
            endif
  120     continue
          det=det*elem
          if (elem.eq.0.) goto 300
          do 140 jc=1,ncol
            augmat(mrow,jc)=augmat(mrow,jc)/elem
  140     continue
          if (mrow.ne.irow) then
c               swap the rows
            do 160 jc=1,ncol
              hold=augmat(mrow,jc)
              augmat(mrow,jc)=augmat(irow,jc)
              augmat(irow,jc)=hold
  160       continue
            det=-det
          endif
          if (irow.lt.nrow) then
            do 190 ir=1,nrow
              if (ir.ne.irow) then
c                 multiply row row by const & subtract from row ir
                const=augmat(ir,jcol)
                do 180 jc=1,ncol
                  augmat(ir,jc)=augmat(ir,jc)-augmat(irow,jc)*const
  180           continue
              endif
  190       continue
          endif
        endif
  100 continue
c
c  Matrix is now in upper triangular form.
c  To solve equation, must get it to the identity.
      do 200 roff=nrow,1,-1
        do 200 irow=roff-1,1,-1
          const=augmat(irow,roff)
          do 200 jc=irow,ncol
            augmat(irow,jc)=augmat(irow,jc)-const*augmat(roff,jc)
  200 continue
      return
      end
c
*********************************************************************
c
      subroutine sq(p)
c
c  This subroutine selects quantities to be written out by a command
c  with type code wsq (write selected quantities).
c  It it essentially a special way of calling the subroutine aim.
c  Written by Alex Dragt 6/15/89
c  Revised (interim) by Tom Mottershead, 30 May 91
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
      dimension p(6), pp(6)
c
c set up control parameters for pass thru to aim
c
      pp(1) = 1.
      pp(2) = p(1)
      pp(3) = p(2)
      pp(4) = 0.
      pp(5) = p(3)
      pp(6) = p(4)
c
      write(jof,*) ' '
      write(jof,*) ' In subroutine sq'
      call aim(pp)
c
      return
      end
c
c********************************************************************
c
      subroutine stopit(loon)
c
c       stops loop number loon = 1 for inner loop
c                              = 2 for outer loop
c
c      T. Mottershead LANL AT-3  7 Oct 91
c----------------------------------------------------
      include 'impli.inc'
      include 'files.inc'
      include 'labpnt.inc'
       write(6,*) ' STOP loon =',loon
      if(loon.eq.2) then
         jocnt = notms
      else
         jicnt = nitms
      endif
      iquiet = 0
      return
      end
c
************************************************************************
c
      subroutine subtip(p)
c  subroutine for terminating an inner procedure
c  Written by Alex Dragt on 8-8-88. Modified 9-3-90 AJD
c
      include 'impli.inc'
      include 'labpnt.inc'
      include 'files.inc'
c
      dimension p(6)
      character*1 yn
c
      iopt=nint(p(1))
      if (iopt.eq.1) then
c         jif = 5
         write(jof,*) ' '
         write(jof,*)
     #   ' Do you wish to jump out of inner procedure? (y/n) <n>:'
         yn='n'
         read(jif,177) yn
 177     format(a)
         if ((yn .eq.'y').or.(yn .eq.'Y')) then
          jicnt=0
          return
         endif
      endif
c
      jicnt=jicnt+1
      if (jicnt.lt.nitms) then
c      write(6,*) 'lp in epro is',lp
       lp=jbipnt
       return
      else
       jicnt=0
       iquiet=0
       return
      endif
c
      end
c
*********************************************************************
c
      subroutine subtop(p)
c  subroutine for terminating an outer procedure
c  Written by Alex Dragt on 8-8-88. Modified 9-3-90 AJD
c
      include 'impli.inc'
      include 'labpnt.inc'
      include 'files.inc'
c
      dimension p(6)
      character*1 yn
c
      iopt=nint(p(1))
      if (iopt.eq.1) then
c         jif = 5
         write(jof,*) ' '
         write(jof,*)
     #   ' Do you wish to jump out of outer procedure? (y/n) <n>:'
         yn='n'
         read(jif,177) yn
 177     format(a)
         if ((yn .eq.'y').or.(yn .eq.'Y')) then
          jocnt=0
          return
         endif
      endif
c
      jocnt=jocnt+1
      if (jocnt.lt.notms) then
c      write(6,*) 'lp in epro is',lp
       lp=jbopnt
       return
      else
       jocnt=0
       iquiet = 0
       return
      endif
c
      end
c
***************************************************************
c
      subroutine txtnum(str,nch,nask,nget,bufr)
c
c      txtnum scans the first nch characters of the string 'str' for
c      real numbers in any format.  all non-numeric characters
c      serve as delimiters.
c
c     parameters:
c          str   - the input character string to be searched
c          nch   - length of character string to be searched
c          nask  - maximum length of number list
c          nget  - actual number of numbers found
c          bufr  - output array of real numbers
c
c     T. Mottershead  Los Alamos  June 1985
c----------------------------------------------------------
      implicit double precision(a-h,o-z)
c-----------------------------------------------------------------------
      dimension bufr(*)
      dimension numc(4) ,ntyp(17)
      character*1  ic
      character*17 numric
      character numdat*30, card*30,  str*(*)
      data ntyp /10*1,2*2,3,4*4/
      numric='0123456789+-.EeDd'
      maxch=len(str)
      if(maxch.gt.nch) maxch=nch
      indx=0
      maxc=30
      nget=0
c     write(6,17)
c 17  format(' err nget -----numc-----',16x,'extracted text',5x
c    $,'value')
c
c     scan for numeric characters
c
   90 do 95 j=1,4
  95  numc(j)=0
      idc=0
  100 continue
      indx=indx+1
      if(indx.le.maxch) go to 110
      if(idc.gt.0) go to 200
      go to 599
  110 ic=str(indx:indx)
      nn=index(numric,ic)
      if(nn.eq.0) go to 200
c
c     valid numeric character found, save in decode buffer:
c
      kk=ntyp(nn)
      numc(kk)=numc(kk)+1
      idc=idc+1
      numdat(idc:idc)=ic
      go to 100
c
c     delimiter character, check for valid string
c
  200 continue
      if(idc.eq.0) go to 100
      if((numc(1).eq.0).or.(numc(2).gt.2)) go to 90
      if((numc(3).gt.1).or.(numc(4).gt.1)) go to 90
c
c     possible string, right justify and decode
c
      if(idc.gt.maxc) idc=maxc
      m=maxc-idc+1
      card=' '
      value=-7777777.
      card(m:maxc)=numdat(1:idc)
      read(card,301,iostat=numerr) value
  301 format(f30.0)
c      write(6,317) numerr,nget+1,numc,card,value
c 317  format(6i4,a30,1pg16.7)
c
c      save the good ones in the output buffer
c
      if(numerr.eq.0) then
         nget=nget+1
         bufr(nget)=value
         if(nget.eq.nask) go to 599
      endif
      go to 90
  599 return
      end
c
ccccccccccccccccccccccc    valuat   ccccccccccccccccccccccccc
c
      function valuat(ku,nsu,idu,jdu)
      include 'impli.inc'
      include 'param.inc'
      include 'fitdat.inc'
      include 'buffer.inc'
      include 'stmap.inc'
      include 'usrdat.inc'
      include 'map.inc'
      include 'rays.inc'
      include 'merit.inc'
      include 'taylor.inc'
      dimension matsig(6,6)
      data matsig /7,8,9,10,11,12,8,13,14,15,16,17,9,14,18,19,20,21,
     * 10,15,19,22,23,24,11,16,20,23,25,26,12,17,21,24,26,27/
c      if(ku.ge.12) then
      if(ku.gt.13) then
         valuat = fitval(ku - 9)
         return
      endif
      go to (10,20,30,40,50,60,70,80,90,100,110,120,130), ku
c
c     stored maps
c
  10  continue
      valuat = sma(idu,jdu,nsu)
      return
c
c       map buffers
c
  20  continue
      valuat = bma(idu,jdu,nsu)
      return
c
c     stored polynomials
c
  30  continue
      valuat = sfa(idu,nsu)
      return
c
c       polynomial buffers
c
  40  continue
      valuat = bfa(idu,nsu)
      return
c
c      map matrix
c
  50  continue
      valuat = tmh(idu,jdu)
      return
c
c       beam sigma matrix
c
  60  continue
      isig = matsig(idu,jdu)
      valuat = bfa(isig,1)
      return
c
c        ray coordinate
c
  70  continue
      valuat = zblock(idu,jdu)
      return
c
c      Taylor map components
c
  80  continue
      valuat = tumat(jdu,idu)
ctm   write(6,*) 'Taylor Map: TM(',idu,jdu,') =',valuat
      return
c
c      current polynomial
c
  90  continue
      valuat = th(idu)
      return
c
c      user data
c
 100  continue
      valuat = ucalc(idu)
      return
c
c      merit function
c
 110  continue
      valuat = val(idu)
      return
c
c      least squares is just flagged
c
 120  continue
      valuat = -777.0
      return
c
c       dispersions
c
 130  continue
      valuat = dz(idu)
      return
      end
c
cccccccccccccccccccccccc   VARLIS   cccccccccccccccccccccccccc
c
      subroutine varlis(jof)
c
c     display list of menu item labels
c
      include 'impli.inc'
      include 'mldex.inc'
      include 'elmnts.inc'
c
      dimension kseq(8), lseq(mnumax)
      character*8 nwrd(8)
c
      call lexort(na,lmnlbl,lseq)
      nnw = (na+7)/8
      write(jof,11)
  11  format(' MARYLIE #menu entries available to be varied:',
     */72('-'))
      do 20 j = 1,nnw
         do 15 nn=1,8
            kseq(nn) = lseq(j+(nn-1)*nnw)
  15     continue
         do 18 i=1,8
            if (kseq(i).eq.0) nwrd(i) = '        '
            if (kseq(i).ne.0) nwrd(i) = lmnlbl(kseq(i))
  18     continue
      write(jof,19) (nwrd(i),i=1,8)
  19  format(8(1x,a8))
  20  continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vary(pp)
c
c     VARY selects the #menu parameters to vary in fitting procedures
c
c      C. T. Mottershead /LANL  & L. B. Schweitzer/UCB
c-----------------------------------------------------
c
      include 'impli.inc'
      include 'param.inc'
      include 'aimdef.inc'
      include 'xvary.inc'
      include 'mldex.inc'
      include 'labpnt.inc'
      include 'files.inc'
      include 'elmnts.inc'
      include 'codes.inc'
c
c lmnlbl(i) = name of ith element (numbered as they occur in #menu)
c pmenu(k+mpp(i))    = kth parameter of ith element
c nt1(i)    = group index of the elements type (see /monics/)
c nt2(i)    = type index  .....
c
      dimension pp(6), menitm(20), indvar(20)
c
c ltc(i,k) = name of element k in group i
c nrp(i,k) = number of parameters
c both are constants, set in block data names
c
      character*128 card, retext
      character*12  vword
      character*1 yn
      logical nodep
      save
c
      loon = nint(pp(5))
      if((loon.lt.1).or.(loon.gt.3)) loon = 1
      call chekit(loon,iter,maxit)
      if(iter.ne.0) return
      ivar=nint(pp(1))
      nodep = .true.
      if(ivar.lt.0) then
         nodep = .false.
         ivar = -ivar
      endif
c      jif = 5
      infile=nint(pp(2))
      if(infile.eq.0) infile = jif
      logf = nint(pp(3))
      iscale = nint(pp(4))
      if(iscale.ne.0) write(jof,*)' VARY:scale and bounds not ready yet'
      isend = nint(pp(6))
      maxrs = 20
      ktot = 0
      jeof = 0
      nmax = nsq(loon)
      if(ivar.gt.2) nmax = maxv
   5  nrem = nmax
      nv = 0
      ibgn = 1
      do 8 j=1,maxv
      idx(j,loon) = 0
      mm(j,loon) = 0
   8  continue
      if(infile.ne.jif) go to 40
  10  call varlis(jof)
c
c     select variable parameters
c
  30  if(infile.ne.jif) go to 40
      write(jof,32)
  32  format(/' Enter name and (optional) parameter index for',
     *' #menu element(s) to be varied.',/' Elements named following a
     * $ may be reset only.  Type * to relist')
      if(ivar.eq.1) then
         write(jof,36) nrem
  36     format(' Selection of',i4,' more elements is required')
      else
         write(jof,34) nrem
  34  format(' Select up to',i4,' element(s): (Enter a # sign when finis
     $hed)')
      endif
c
c       read a new input line
c
  40  continue
      call nexlin(infile,jeof,card)
      if(jeof.lt.0) go to 200
ctm14  40  read(infile,41,end=200) card
  41  format(a)
c     write(jof,41) card
      if(card.eq.' ') go to 40
      call low(card)
      ifin = index(card,'#')
      list = index(card,'*')
      idol = index(card,'$')
      retext = ' '
c
c       check for reset commands
c
      if(idol.gt.0) then
         retext = card(idol+1:80)
         card(idol+1:80) = ' '
         call vreset(infile,retext,maxrs,ktot,menitm,indvar)
      endif
c
c     look for normal variable selections
c
      call vscan(card,kount,menitm,indvar)
      write(jof,53) card,kount
  53  format(/,'VARY card:',a,/,i3,' items found:')
c
c      finish the incomplete specifications
c
      do  100  kr = 1, kount
      mnu  = menitm(kr)
      idu  = indvar(kr)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      npar=nrp(mt1,mt2)
      if(npar.eq.1) idu = 1
      nv = nv + 1
c      write(jof,*) ' card',nv,'=nv',kr,'=kr',mnu,'=mnu',npar,'=npar'
      mm(nv,loon) = mnu
c
c      ask for new selections for out of bounds parameters
c
  80  if((idu.gt.0).and.(idu.le.npar)) go to 85
      write(jof,81) nv, lmnlbl(mnu), ltc(mt1,mt2)
  81  format(' No.',i3,' is ',2a8,'.  Select parameter to vary:')
      write(jof,83) (pmenu(k+mpp(mnu)),k = 1,npar)
  83  format(' p=',6(1pg12.4))
      read(jif,*) idu
      go to 80
c
c       save valid parameter index and show what is picked so far
c
  85  idx(nv,loon) = idu
  90  continue
c 90  write(jof,91) nv, lmnlbl(mnu), ltc(mt1,mt2), idu, npar
c 91  format(' No.',i3,' is ',a8,1x,a8,'.  Parameter',i2,' out of',i2,
c    *' selected.')
      nch = index(lmnlbl(mnu),' ')
c      write(jof,92) nch, mnu, lmnlbl(mnu)
c  92  format('  nch,mnu=',2i6,'  name = ',a)
      if(nch.eq.0) nch=9
      write(vword,93) idu
  93  format(9x,'(',i1,')')
      vword(11-nch:9) = lmnlbl(mnu)(1:nch-1)
      varnam(nv,loon) = vword
      write(jof,95) varnam(nv,loon),pmenu(idu+mpp(mnu))
  95  format(1x,a,' = ',1pg21.14)
 100  continue
c
c     See if we are finished. If not, go back to read another card
c
      nrem = nmax - nv
ctm   write(jof,*) ' VARY end check:',nrem,ifin,list
      if(nrem.le.0) go to 200
      if(ifin.ne.0) go to 200
      if(list.ne.0) go to 10
      go to 30
c
c     Display final selections
c
 200  continue
      kprt = 0
      if((isend.eq.1).or.(isend.eq.3)) lun = jof
      if(isend.eq.2) lun = jodf
 201  write(lun,203)
 203  format(/'   Variable #menu elements selected:')
      write(lun,175)
 175  format(' No.  Element',4x,'Type',5x,'Parameter',3x,'Present value'
     * /60('-'))
      do 170 i = 1,nv
      idu = idx(i,loon)
      mnu = mm(i,loon)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      write(lun,176) i,lmnlbl(mnu),ltc(mt1,mt2),idu,pmenu(idu+mpp(mnu))
 176  format(i3,4x,a8,3x,a8,i6,5x,1pg21.14)
 170  continue
      if((isend.eq.3).and.(kprt.eq.0)) then
          kprt = 1
          lun = jodf
          go to 201
      endif
c
      if(infile.eq.jif) then
         write(jof,*) ' '
         write(jof,*) ' Do you wish to start over? (y/n) <n>:'
         yn='n'
         read(jif,41) yn
         if ((yn .eq.'y').or.(yn .eq.'Y')) go to 5
      endif
c
c     select dependent parameters
c
      nd = 0
      if(nodep) go to 400
 300  continue
      nd = 0
      nmax = maxv - nv
      nrem = nmax
 430  if(infile.ne.jif) go to 340
      write(jof,332) nrem
 332  format(/' Define up to',i4,' dependent variables V by entering the
     $ir names',/,' and (optional) parameter indices. (Enter a # sign wh
     $en finished)')
c
c       read a new input line
c
ctm14  340  read(infile,41,end=600) card
 340  call nexlin(infile,jeof,card)
      if(jeof.lt.0) go to 600
c     write(jof,41) card
      if(card.eq.' ') go to 340
      call low(card)
      ifin = index(card,'#')
      list = index(card,'*')
c
c     look for normal variable selections
c
         call vscan(card,kount,menitm,indvar)
c
c      finish the incomplete specifications
c
      if(kount.eq.0) go to 530
      do  500  kr = 1, kount
      mnu  = menitm(kr)
      idu  = indvar(kr)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      npar=nrp(mt1,mt2)
      if(npar.eq.1) idu = 1
      nd = nd + 1
      mm(nd + nv,loon) = mnu
c
c      ask for new selections for out of bounds parameters
c
 380  if((idu.gt.0).and.(idu.le.npar)) go to 355
      write(jof,81) nd, lmnlbl(mnu), ltc(mt1,mt2)
      write(jof,83) (pmenu(k+mpp(mnu)),k = 1,npar)
      read(jif,*) idu
      go to 380
c
c           read new parameter value
c
 355  idx(nd + nv,loon) = idu
c     write(jof,91) nd, lmnlbl(mnu), ltc(mt1,mt2), idu, npar
      if(infile.eq.jif) write(jof,357) nv
 357  format(' Enter ID (1 to',i4,') of independent variable x , and the
     * derivative dV/dx:')
ctm14    read(infile,*) ivn, deriv
      call nexlin(infile,jeof,card)
      read(card,*) ivn, deriv
      idv(nd,loon) = ivn
      xslope(nd,loon) = deriv
      pval = pmenu(idx(ivn,loon)+mpp(mm(ivn,loon)))
      xbase(nd,loon) = pmenu(idu+mpp(mnu)) - xslope(nd,loon)*pval
 500  continue
c
c     See if we are finished. If not, go back to read another card
c
 530  nrem = nmax - nd
      if(nrem.le.0) go to 600
      if(ifin.ne.0) go to 600
      if(list.ne.0) go to 200
      go to 430
c
c     Display final selections
c
 600  continue
      if(nd.le.0) go to 650
      kprt = 0
      if((isend.eq.1).or.(isend.eq.3)) lun = jof
      if(isend.eq.2) lun = jodf
 301  write(lun,303)
 303  format(/'   Dependent #menu elements selected:')
      write(lun,375)
 375  format(' No.  Element',4x,'Type',5x,'Parameter',3x,'Present value'
     *,6x,' IDV  Slope', /72('-'))
      do 570 i = 1,nd
      idu = idx(i+nv,loon)
      mnu = mm(i+nv,loon)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      write(lun,376) i,lmnlbl(mnu),ltc(mt1,mt2),idu,pmenu(idu+mpp(mnu)),
     * idv(i,loon), xslope(i,loon)
 376  format(i3,3x,a8,3x,a8,i4,5x,1pg21.14,i5,0pf10.5)
 570  continue
      if((isend.eq.3).and.(kprt.eq.0)) then
          kprt = 1
          lun = jodf
          go to 301
      endif
c
 650  if(infile.eq.jif) then
         write(jof,*) ' '
         write(jof,*) ' Do you wish to start over? (y/n) <n>:'
         yn='n'
         read(jif,41) yn
         if ((yn .eq.'y').or.(yn .eq.'Y')) go to 300
      endif
c
c       write final selections to log file
c
 400  if(logf.gt.0) then
         do 410 i = 1,nv
         write(logf,407)  lmnlbl(mm(i,loon)), idx(i,loon)
 407     format(2x,a8,2x,i2)
 410     continue
         if(nd.gt.0) then
            if(ivar.gt.1) write(logf,*) '  #'
            do 450 j = 1,nd
            i = j + nv
            write(logf,407)  lmnlbl(mm(i,loon)), idx(i,loon)
            write(logf,*) idv(j,loon), xslope(j,loon)
 450        continue
         endif
         write(logf,*) '  #'
      endif
      nva(loon) = nv
      nda(loon) = nd
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vreset(infile,retext,maxrs,ktot,mrset,irset)
      include 'impli.inc'
      include 'files.inc'
      include 'elmnts.inc'
      include 'codes.inc'
      dimension menitm(20), indvar(20), mrset(*), irset(*)
      character*(*) retext
      jout = jof
      if(infile.ne.jif) jout = jodf
c
c lmnlbl(i) = name of ith element (numbered as they occur in #menu)
c pmenu(k+mpp(i))    = kth parameter of ith element
c nt1(i)    = group index of the elements type (see /monics/)
c nt2(i)    = type index  .....
c ltc(i,k) = name of element k in group i
c nrp(i,k) = number of parameters
c both are constants, set in block data names
c
         call vscan(retext,kount,menitm,indvar)
       write(jof,*) '  after vscan:',kount,'=kount'
       write(jof,*) ' menitm:',(menitm(j),j=1,kount)
       write(jof,*) ' indvar:',(indvar(j),j=1,kount)
c
c     complete the resets
c
         do  60  kr = 1, kount
         mnu  = menitm(kr)
         idu  = indvar(kr)
         mt1 = nt1(mnu)
         mt2 = nt2(mnu)
         npar=nrp(mt1,mt2)
         if(npar.eq.0) go to 60
       write(jof,*) ' vreset',kr,'=kr',mnu,'=mnu',idu,'=idu',mt1,'=mt1'
c     *,mt2,'=mt2', npar,'=npar'
         if(npar.eq.1) idu = 1
c
c      ask for new selections for out of bounds parameters
c
  50   if((idu.gt.0).and.(idu.le.npar)) go to 55
       write(jout,51) lmnlbl(mnu), ltc(mt1,mt2)
  51   format(' Select parameter to reset for ',2a8)
       write(jout,53) (pmenu(k+mpp(mnu)),k = 1,npar)
  53   format(1x,6(1pe13.5))
       read(infile,*) idu
       go to 50
c
c           valid index, read new parameter value
c
  55     pval = pmenu(idu+mpp(mnu))
         write(jout,57) idu, lmnlbl(mnu), pval
  57     format(' Enter new value for parameter',i2,' of ',a8,
     *   ' <',1pg16.8,'>:')
         read(infile,*) pval
         pmenu(idu+mpp(mnu)) = pval
         ktot = ktot + 1
         if(ktot.le.maxrs) then
            mrset(ktot) = mnu
            irset(ktot) = idu
         endif
  60     continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine vscan(card,kount,menitm,indvar)
c
c          VSCAN scans the character array card for valid #menu
c          items (menitm) and parameter selections (indvar)
c          C. T. Mottershead /LANL AT-3  18 July 90
c-----------------------------------------------------
      include 'impli.inc'
      include 'elmnts.inc'
      include 'codes.inc'
c
c nt1(i)    = group index of the elements type (see /monics/)
c nt2(i)    = type index  .....
c nrp(i,k) = number of parameters
c both are constants, set in block data names
c
      dimension menitm(*), indvar(*)
      character*(*) card
      character*10 string
      logical lfound, lnum, leftel
      save
c
c       parse the input card
c
      kbeg = 1
      leftel = .false.
      kount = 0
      npar = 0
  50  continue
      msegm=2
      call cread(kbeg,msegm,card,string,lfound)
ctm   write(6,*) 'VSCAN cread:',kount,'=kount',kbeg,string
      if(.not.lfound) return
c
c      string found, see if it is an element
c
      call lookup(string,itype,mnu)
      if(itype.ne.1) go to 60
c
c        it is an element
c
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      npar=nrp(mt1,mt2)
ctm   write(6,43) mnu,mt1,mt2,npar
c 43  format('VSCAN element: mnu=',4i4,'=npar')
c
c       reject if no parameters
c
      if(npar.eq.0) then
          leftel = .false.
          go to 50
      endif
c
c      accept it because it has parameters
c
      kount = kount + 1
      menitm(kount) = mnu
      indvar(kount) = 0
      leftel = .true.
      go to 50
c
c      it isn't an element, see if it is a number
c
  60  call cnumb(string,num,lnum)
c
c     if it is a number, and a preceeding element is defined, and it is
c     in bounds, interpret it as the selected parameter for that element
c
      if(lnum.and.leftel) then
         if((num.ge.1).and.(num.le.npar)) indvar(kount) = num
      endif
      leftel = .false.
      go to 50
      end
c
cccccccccccccccccccccc   WRTAIM  ccccccccccccccccccccccccc
c
      subroutine wrtaim(lun,kaim,loon)
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      include 'aimdef.inc'
      include 'fitbuf.inc'
      include 'files.inc'
      dimension aimbuf(3)
      if(kaim.le.3) write(lun,*)' '
      if(kaim.le.3) write(lun,*)' Aims selected : '
      if(kaim.eq.1) write(lun,103)
 103  format(' No.     item',6x,'present value',/35('-'))
      if(kaim.eq.2) write(lun,105)
 105  format(' No.     item',8x,'present value',8x,'target value',
     * /53('-'))
      if(kaim.eq.3) write(lun,107)
 107  format(' No.     item',8x,'present value',8x,'target value',
     * 8x,'weights', /70('-'))
      jmin = 1
      jmax = kaim
      if(kaim.gt.3) then
         jmin = 2
         jmax = 3
      endif
      nf = nsq(loon)
ctm   write(6,*) ' WRTAIM:',nf,'=nf',kaim,'=kaim',lun,'=lun'
      do 100 nu = 1, nf
      ku = kyf(nu,loon)
      nsu = nsf(nu,loon)
      idu = idf(nu,loon)
      jdu = jdf(nu,loon)
ctm   write(6,97) loon,nu,ku,nsu,idu,jdu
c 97  format('in WRTAIM:loon,nu=',2i3,'  ku,nsu=',2i3,'  idu,jdu=',2i3)
      aimbuf(1) = valuat(ku,nsu,idu,jdu)
      aimbuf(2) = target(nu,loon)
      aimbuf(3) = wts(nu,loon)
      write(lun,96) nu, qname(nu,loon), (aimbuf(j),j=jmin,jmax)
  96  format(i3,'  ',a,' = ',2(1pg20.9),0pf14.4)
 100  continue
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wrtsq(lun,kform,nf,names,values)
c
c  Output selected quantities array to unit lun in format kform.
c
c      kform = format selection flag
c            = 0 to write names only
c            = 1 to write a line of values (6 digit accuracy)
c            = 2 to write 3 names and values per line (8 digit accuracy)
c            = 3 to write a single value per line, full precision
c      Tom Mottershead, 31 May 91
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      dimension values(*)
      character*(*) names(*)
c       write(lun,*) ' *WRTSQ:',nf,' = number items'
c
c      format 0: labels only
c
      if(kform.eq.0) then
        write(lun,63) (names(i), i=1,nf)
  63    format(10(1x,a12))
      endif
c
c      format 1: 6 digit columns
c
      if(kform.eq.1) then
ctm   write(6,*) 'WRTSQ >>> ',kform,'=kform',lun,'=lun  <<<'
        write(lun,67) (values(i), i=1,nf)
  67    format(10(1pe13.5))
      endif
c
c      format 2: 3 on a line
c
      if(kform.eq.2) then
c         write(lun,43) ((names(i),values(i)),i=1,nf)
         write(lun,43) (names(i),values(i),i=1,nf)
  43     format(3(a12,'=',1pg13.6))
      endif
c
c      format 3: full precision
c
      if(kform.eq.3) then
         do 25 nu = 1, nf
         write(lun,27) nu, names(nu), values(nu)
  27     format(' SQ',i3,':    ',a,' = ', 1pg23.15)
  25     continue
      endif
      return
      end
c
c*******************************************************************
c
      subroutine wrtvar(kpik,lun)
c
c     Display new values for parameters.
c
      include 'impli.inc'
      include 'elmnts.inc'
      include 'codes.inc'
      include 'xvary.inc'
      loon = kpik
      if((loon.lt.1).or.(loon.gt.3)) loon = 1
      nv = nva(loon)
      nd = nda(loon)
      write(lun,*)' '
      write(lun,*)' New values for parameters:'
      write(lun,65)
  65  format(' No.  Element',4x,'Type',3x,'Parameter',3x,'Present value'
     *,8x,' IDV  Slope', /72('-'))
      do 50 i = 1,nv
      idu = idx(i,loon)
      mnu = mm(i,loon)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      write(lun,48) i,lmnlbl(mnu),ltc(mt1,mt2),idu,pmenu(idu+mpp(mnu))
  48  format(i3,3x,a8,3x,a8,i4,5x,1pg21.14)
  50  continue
      do 70 i = 1,nd
      j = i + nv
      idu = idx(j,loon)
      mnu = mm(j,loon)
      mt1 = nt1(mnu)
      mt2 = nt2(mnu)
      write(lun,68) j,lmnlbl(mnu),ltc(mt1,mt2),idu,pmenu(idu+mpp(mnu)),
     * idv(i,loon), xslope(i,loon)
  68  format(i3,3x,a8,3x,a8,i4,5x,1pg21.14,i5,0pf10.4)
  70  continue
      return
      end
c
ccccccccccccccccccc    wsq   cccccccccccccccccccccccccc
c
      subroutine wsq(pp)
c subroutine for writing out selected quantities
c Written by Alex Dragt, 21 August 1988
c Filled in by Tom Mottershead, 30 May 91
c
c-----!----------------------------------------------------------------!
      include 'impli.inc'
      include 'files.inc'
      include 'aimdef.inc'
      include 'elmnts.inc'
      include 'xvary.inc'
      include 'usrdat.inc'
      dimension pp(6), qbuf(30)
      character*12 names(30), ctemp
      ctemp = ' '
      job = nint(pp(1))
      loon = nint(pp(2))
      if((loon.lt.1).or.(loon.gt.3)) loon = 1
      lun  = nint(pp(3))
      kform = nint(pp(4))
      jform = nint(pp(5))
      isend = nint(pp(6))
      nv = 0
      nf = 0
      numq = 0
ctm   write(jof,*)job,loon,lun,kform,jform,isend,' = WSQ (PP)'
      if(job.eq.4) go to 100
c
c       copy variables to print buffers
c
      if(job.eq.1) go to 50
      nv = nva(loon)
      numq = nv
      if(numq.gt.30) numq = 30
      do 40 i = 1,numq
      mnu = mm(i,loon)
      qbuf(i) = pmenu(idx(i,loon)+mpp(mnu))
      names(i) = varnam(i,loon)
  40  continue
c
c        copy selected quantities (aims) to print buffers
c
  50  if(job.eq.2) go to 100
      nf = nsq(loon)
      numq = nf + nv
      if(numq.gt.30) then
         numq = 30
         nf = 30 - nv
      endif
      do 60 nu = 1, nf
      ku = kyf(nu,loon)
      nsu = nsf(nu,loon)
      idu = idf(nu,loon)
      jdu = jdf(nu,loon)
      qbuf(nu+nv) = valuat(ku,nsu,idu,jdu)
      ctemp(3:) = qname(nu,loon)
      names(nu+nv) = ctemp
  60  continue
c
c     print any rms errors selected in current loon
c
 100  continue
      do 150 ku = 1,3
      lu = lsq(ku,loon)
      if(lu.ne.1) go to 150
      numq = numq+1
      fval = rmserr(ku)
      qbuf(numq) = fval
      write(names(numq),141) lu
 141  format(3x,'rms',i1,5x)
 150  continue
c
c     print buffers
c
c      write(jof,*)' WSQ:  nva =', nva
c      write(jof,*)' WSQ:  nsq =', nsq
c      write(jof,*)' * WSQ:',nv,'=nv',nf,'=nf',numq,'=numq total'
      if(lun.gt.0) call wrtsq(lun,kform,numq,names,qbuf)
      if(isend.ge.2) call wrtsq(jodf,jform,numq,names,qbuf)
      if((isend.eq.1).or.(isend.eq.3))
     *     call wrtsq(jof,jform,numq,names,qbuf)
c
c       copy to ucalc if lun<0  22 Aug 00  CTM from AJD's version
c
      if(lun.lt.0) then
         ibgn = -lun
         do 200 j = 1, numq
           ucalc(ibgn+j-1) = qbuf(j)
  200    continue
      endif
      return
      end
