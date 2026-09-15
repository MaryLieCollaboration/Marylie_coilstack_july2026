************************************************************************
* header              GENREC (GENMAP for a pattern of REC quads)       *
*  All routines needed for this special GENMAP                         *
************************************************************************
c
      subroutine chkdet(y)
      include 'impli.inc'
ctm2014   changes for GFORTRAN. It seems to dislike y(7) arguement.
ctm2014 This was subr errchk -> renamed chkdet
      dimension y(*)
      mm = 6
      s1=1.d0 -(y(1+mm)*y(8+mm)-y(2+mm)*y(7+mm))
      s2=1.d0 -(y(15+mm)*y(22+mm)-y(16+mm)*y(21+mm))
      write(6,100)s1,s2
  100 format(1x,'1. - 2 x 2 determinants = ',e14.7,1x,e14.7)
      return
      end
c
************************************************************************
c
      subroutine g012(z,g,gz,gzz)
c  This routine computes g, dg/dz, and d/dz(dg/dz) for
c  a REC Quadrupole Multiplet.
c  Written by Alex Dragt, Fall 1986, and based on work of
c  Rob Ryne and F. Neri
c  This version by T. Mottershead allows up to 9(?) different
c  REC Quadrupoles.
      include 'impli.inc'
      include 'quadpn.inc'
c----------------------------------------
      v(z,r)=1./sqrt(1.+(z/r)**2)
      f(z,a,b)=0.5-
     #.0625*z*(1./a+1./b)*v(z,a)**2*v(z,b)**2/(v(z,a)+v(z,b))*
     #(v(z,a)**2+v(z,a)*v(z,b)+v(z,b)**2+4.+8./(v(z,a)*v(z,b)))
      fz(z,a,b)=-.1875*(1./a+1./b)*v(z,a)**2*v(z,b)**2*
     #(v(z,a)**3+v(z,b)**3 + v(z,a)**2*v(z,b)**2/(v(z,a)+v(z,b)))
      fzz(z,a,b)= .1875*(1./a+1./b)*z*v(z,a)**2*v(z,b)**2 *(
     #5.*(v(z,a)**5/a**2+v(z,b)**5/b**2)
     #+v(z,a)**2*v(z,b)**2 *(
     # 2.*(v(z,b)/a**2+v(z,a)/b**2)
     #+4.*(v(z,a)**2/a**2+v(z,b)**2/b**2)/(v(z,a)+v(z,b))
     #-(v(z,a)**3/a**2+v(z,b)**3/b**2)/(v(z,a)+v(z,b))**2 ))
c----------------------------------------
c
      g=0
      gz=0
      gzz=0
      kuad = 0
      za = di
      do 60 i=1,ncyc
      do 50 k=1,nqt
      gk = ga(k)
      wk = wd(k)
      dk = dr(k)
      aa = ra(k)
      bb = rb(k)
      do 40 j=1,nr(k)
      zb = za + wk
      g  = g  +gk*(  f(z-zb,aa,bb)-f(z-za,aa,bb))
      gz = gz +gk*( fz(z-zb,aa,bb)-fz(z-za,aa,bb))
      gzz= gzz+gk*(fzz(z-zb,aa,bb)-fzz(z-za,aa,bb))
      kuad = kuad + 1
      if(kuad.ge.maxq) return
      za = zb + dk
  40  continue
  50  continue
  60  continue
      return
      end
c
************************************************************************
c
      subroutine gnrec3(p,fa,fm)
c
c This is a subroutine for computing the map for patterns of REC quads.
c It is based on Rob Ryne's GENMAP as modified by Alex Dragt 12/18/86.
c
c Further modified 4 June 87 by Tom Mottershead to allow any number
c of cycles through a pattern of up to 5 REC quads, with arbitrary
c numbers of repeats.  Actually probably more general than this. AJD
c
c Modified by Filippo Neri Jan. 16 1989 to include multipoles.
c Further modified by Alex Dragt on Bastille Day 1989 to slightly
c change input and output format.
c
c The input parameters come from recn, and the
c auxilary parameter NPQUAD, MULTIPOLES, plus NQT pset defining the quads.
c
c  The parameters of recm are:
c  1.    zi  Initial integration point.
c  2.    zf  Final integration point.
c  3.    NS  Number of integration steps.
c  4.    IFILE (profile file number).
c  5.    MULTIPOLES index of pset containin multipole information
c                   for this integration region. If MULTIPOLES is 0, then
c                   the value of all multipoles is set to zero.
c  6.    NPQUAD     index of pset containing quad pattern information.
c---------------------------------------------------------------
c   The pset MULTIPOLES has the same format as used in the cfq element:
c  1.    Normal sextupole ( Tesla/meter^2 ).
c  2.    Skew   sextupole ( Tesla/meter^2 ).
c  3.    Normal octupole  ( Tesla/meter^3 ).
c  4.    Skew   octupole  ( Tesla/Meter^3 ).
c  5.    UNUSED ( must be there!)
c  6.    UNUSED ( nust be there!)
c---------------------------------------------------------------
c  The parameter set NPQUAD defines the set of Halbach quad units as follows:
c  1.    Di   Initial drift to first quad ( starting from Z = 0.)
c  2.    NQT  Number of psets used for quads.
c  3.    IPS  Index of inital pset used for quads. The following
c             psets in order define the successive quads.
c  4.    MAXQ Maximum number of quads actually used, including multiple
c             cycles thru pattern (but not repeats).
c  5.    NCYC Number of cycles thru pattern, except that the sequence
c             stops when MAXQ is reached.
c  6.    ISEND     output control:  0 = quiet running
c                                   1 = print on terminal (jof)
c                                   2 = print on std. output file (jodf)
c                                   3 = print on both
c---------------------------------------------------------------
c     The parameter set type codes are used to define the Halbach quad unit
c     cells in a manner parallel to the normal quad type codes:
c      ps(1) = ql, the quad length in meters.
c      ps(2) = fg, the field gradient in Tesla/meter (if the quad were
c                  infinitely long).
c      ps(3) = ra, the inner radius.   (These radii control the shape of the
c      ps(4) = rb, the outer radius.       fringe field)
c      ps(5) = td, the trailing drift to the front face of the next quad in
c                  the pattern (meters). Note: the trailing drift assigned to
c                  the last quad in the whole pattern is ignored; the
c                  integration proceeds to the final point zf anyway.
c      ps(6) = nr, the number of consecutive repetitions, or multiplicity,
c                  of this quad unit cell in the basic pattern.
c-------------------------------------------------------------------
c
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'quadpn.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'recmul.inc'
c
c  calling arrays
      dimension p(6)
      dimension fa(monoms), fm(6,6)
c
c  local arrays
      dimension pb(6)
      dimension y(monoms+15)
c
c timing variables
c      real ttaa, ttbb
c
c use equivalence statement to make the various parameter sets pstj
c available as if they were in a two dimensional array
c
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c
c  get interval and number of steps from GENREC parameters
c
      zi = p(1)
      zf = p(2)
      ns = nint(p(3))
      ifile = nint(p(4))
      mpole = nint(p(5))
      npquad = nint(p(6))
c
c  get multipole values from the parameter set mpole
c  Note that if mpole is zero, then the multipoles are set to zero.
c
      if (mpole.lt.1 .or. mpole.gt.maxpst) then
        do 500 i=1,6
  500   pb(i) = 0.0d0
      else
        do 600 i=1,6
  600   pb(i) = pst(i,mpole)
      endif
c
c compute multipole coefficients
c
      bsex=pb(1)
      asex=pb(2)
      boct=pb(3)
      aoct=pb(4)
      fsxnr=bsex/(3.d0*brho)
      fsxsk=asex/(3.d0*brho)
      focnr=boct/(4.d0*brho)
      focsk=aoct/(4.d0*brho)
c
c  get REC pattern information from pset npquad
c
      di = pst(1,npquad)
      nqt = nint(pst(2,npquad))
      if (nqt .gt. 6) then
      write(jof,*) ' input error: nqt=',nqt
      write(jof,*) ' this value is > 6 and therefore too large'
      call myexit
      endif
      ips = nint(pst(3,npquad))
      maxq = nint(pst(4,npquad))
      ncyc = nint(pst(5,npquad))
      isend = nint(pst(6,npquad))
      jtty = 0
      jdsk = 0
      if((isend.eq.1).or.(isend.eq.3)) jtty = 1
      if((isend.eq.2).or.(isend.eq.3)) jdsk = 1
c
      h=(zf-zi)/float(ns)
c
c     echo input parameters
c
      if(jtty.eq.1) write( jof,13) ns,zi,zf
      if(jdsk.eq.1) write(jodf,13) ns,zi,zf
  13  format(/' Integrating in ',i6,' steps from',f10.5,'=zi to',
     *f10.5,'=zf')
      if (mpole.lt.1 .or. mpole.gt.maxpst) then
      if(jtty.eq.1) write( jof,*)
     #' mpole=0, all multipoles are zero'
      if(jdsk.eq.1) write( jodf,*)
     #' mpole=0, all multipoles are zero'
      endif
      if (mpole.ge.1 .and. mpole.le.maxpst) then
      if(jtty.eq.1) write( jof,14) mpole,bsex,asex,boct,aoct
      if(jdsk.eq.1) write(jodf,14) mpole,bsex,asex,boct,aoct
  14  format(1x,'multipole strengths from pset',i2,':',/,
     *1x,' Sextupole:',2(1pg12.4),/,
     *1x,' Octupole: ',2(1pg12.4))
      endif
      if(jtty.eq.1) write( jof,15) npquad,ncyc,nqt,maxq,di
      if(jdsk.eq.1) write(jodf,15) npquad,ncyc,nqt,maxq,di
  15  format(' Pattern from pset',i2,':',/,
     *i4,' cycle(s) of',i3,' type(s) of section(s) with a maximum of'
     *,i3,' section(s)',/,'  di=',f10.5)
      if(jtty.eq.1) write( jof,*) ' types of section(s) are:'
      if(jdsk.eq.1) write(jodf,*) ' types of section(s) are:'
      if(jtty.eq.1) write( jof,16)
      if(jdsk.eq.1) write(jodf,16)
  16  format(1x,'pset',2x,'length',4x,'strength',3x,
     *'radii: inner',6x,'outer',6x,'tdrift',4x,'number')
c
c      inititialize the nqt quad types from the first nqt parameter sets
c
      do 20 n = ips, nqt+ips-1
      kn = n-ips+1
      wd(kn) = pst(1,n)
      ga(kn) = pst(2,n)
      ra(kn) = pst(3,n)
      rb(kn) = pst(4,n)
      dr(kn) = pst(5,n)
      nr(kn) = nint(pst(6,n))
      if(jtty.eq.1) write(jof,17)
     # n,wd(kn),ga(kn),ra(kn),rb(kn),dr(kn),nr(kn)
      if(jdsk.eq.1) write(jodf,17)
     # n,wd(kn),ga(kn),ra(kn),rb(kn),dr(kn),nr(kn)
  17  format(i4,5f12.6,i6)
  20  continue
c
c     call VAX system routine for timing report
c
c      ttaa = secnds(0.0)
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=0.d0
      y(2)=0.d0
      y(3)=0.d0
      y(4)=0.d0
      y(5)=0.d0
      y(6)=0.d0
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c      cmp2=(1.d0/(beta*gamma))**2
c
c  initialize map to the identity map:
      ne=224
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c  Set up multipoles:
c
c     compute useful numbers
c
      sl2=sl*sl
      sl3=sl*sl2
      bet2=beta*beta
      bet3=beta*bet2
      gam2=gamma*gamma
c
c  write out gradient and derivatives on file ifile:
      ipflag=0
      if (ifile .ne. 0) then
      if (ifile .lt.0 ) then
      ipflag=1
      ifile=-ifile
      endif
      do 100 i=0,ns
      s=zi + float(i)*h
      call g012(s,g,gz,gzz)
  100 write(ifile,101)s,g,gz,gzz,0.,0.
  101 format(6(1x,1pg12.5))
      write(jof,*) ' '
      write(jof,*) ' profile written on file ',ifile
      endif
c
c  return identity map if ifile was < 0
      if (ipflag .eq. 1) then
      call ident(fa,fm)
      return
      endif
c
c  do the computation:
      t=zi
      iflag = 2
      call adam11(h,ns,'start',t,y)
      call chkdet(y)
      call putmap(y,fa,fm)
      call csym(1,fm,ans)
c
c     call VAX system routine for timing report
c
c      ttbb = secnds(ttaa)
c      if(jtty.eq.1) write( jof,567) ttbb
c      if(jdsk.eq.1) write(jodf,567) ttbb
c 567  format(' GENREC integration time = ',f12.2,' sec.')
c
      return
      end
c
**********************************************************************
c
      subroutine hmltn2(t,y,h)
c
c  this routine is used to specify h(z) for a REC quad multiplet
c  Written by Alex Dragt, Fall 1986, and based on work of
c  Rob Ryne and F. Neri
c  Modified by F. Neri, Jan 11 1989 to include
c  sextupole and octupole fields added to the
c  REC quadrupole.
c
      include 'impli.inc'
      include 'param.inc'
      include 'parm.inc'
      include 'combs.inc'
      include 'recmul.inc'
      dimension h(monoms),y(*)
c
c  begin calculation
c
c  compute gradients
      call g012(t,g,gz,gzz)
c
c  initialization
      do 10 i=7,monoms
   10 h(i)=0.d0
c
c 2nd order
      h(7)=0.5d0*g*sl*qbyp
      h(13)=0.5d0/sl
      h(18)=-h(7)
      h(22)=h(13)
      h(27)=0.5d0/sl*(ptg**2-1.d0)
c 3rd order
c with sextupoles
      h(28)=fsxnr*sl2
      h(30)=-3.d0*fsxsk*sl2
      h(39)=-3.d0*fsxnr*sl2
      h(53)=-0.5d0*ptg/sl
      h(64)=fsxsk*sl2
      h(76)=h(53)
      h(83)=0.5d0/sl*(ptg-ptg**3)
c 4th order
      h(84)=-gzz*sl**3*qbyp/12.d0
      h(175)=-h(84)
      h(85)=-0.25d0*gz*sl**2*qbyp
      h(110)=-h(85)
      h(96)=h(85)
      h(176)=-h(85)
      h(140)=0.125d0/sl
      h(195)=h(140)
      h(149)=2.d0*h(140)
      h(154)=0.25d0*(3.d0*ptg**2-1.d0)/sl
      h(200)=h(154)
      h(209)=0.125d0*(ptg**2-1.d0)*(5.d0*ptg**2-1.d0)/sl
c
c  add octupoles
c
      h(84)=h(84)+focnr*sl3
      h(86)=-4.0d0*focsk*sl3
      h(95)=-6.d0*focnr*sl3
      h(120)=4.d0*focsk*sl3
      h(175)=h(175)+focnr*sl3
c
      return
      end
