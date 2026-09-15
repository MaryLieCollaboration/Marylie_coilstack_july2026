************************************************************************
* header:                 COMMANDS (simple)
*  These modules deal with simple commands that affect maps and phase
*  space data.  Ray trace related commands are found in TRAC, and input
*  output related commands are found in INPU.
************************************************************************
c
      subroutine cmom(amom)
c this subroutine computes the moments of the particle distribution
c stored in zblock.  The moments are put in the array amom.
c Written by Alex Dragt, 1 July 1991.
c
      include 'impli.inc'
      include 'param.inc'
      include 'rays.inc'
c
c calling arrays
      dimension amom(monoms)
c
c working arrays
      dimension vmon(monoms)
      dimension z(6)
c
c clear amom
      do 5 i=1,monoms
    5 amom(i)=0.d0
c
c computational loop
      do 10 k=1,nrays
c get a ray
      do 20 i=1,6
   20 z(i)=zblock(k,i)
c compute moments for the ray
      call evalm(z,vmon)
c add results to moment sum
      do 30 i=1,monoms
   30 amom(i)=amom(i)+vmon(i)
   10 continue
c
c normalize by the number of particles
      fact=1.d0/float(nrays)
      do 40 i=1,monoms
   40 amom(i)=fact*amom(i)
c
      return
      end
c
      subroutine cmom5(amom)
c this subroutine computes the moments of the particle distribution
c stored in zblock.  The moments are put in the array amom.
c Written by Alex Dragt, 1 July 1991.
c Modified to generate <5> and <6>, by Johannes van Zeijts, 11 November 1991.
c
      include 'impli.inc'
      include 'param.inc'
      include 'rays.inc'
c
c calling arrays
      dimension amom(monoms5)
c
c working arrays
      dimension vmon(monoms5)
      dimension z(6)
c
c clear amom
      do 5 i=1,monoms5
    5 amom(i)=0.d0
c
c computational loop
      do 10 k=1,nrays
c get a ray
      do 20 i=1,6
   20 z(i)=zblock(k,i)
c compute moments for the ray
      call evalm5(z,vmon)
c add results to moment sum
      do 30 i=1,monoms5
   30 amom(i)=amom(i)+vmon(i)
   10 continue
c
c normalize by the number of particles
      fact=1.d0/float(nrays)
      do 40 i=1,monoms5
   40 amom(i)=fact*amom(i)
c
      return
      end
c
******************************************************************
c
      subroutine dpol(p,fa,fm)
c This is a subroutine for setting up a quadratic polynomial
c described in terms of dispersion parameters.
c Written by Alex Dragt, 29 March 1991
c
      include 'impli.inc'
      include 'param.inc'
c
      dimension fa(monoms)
      dimension fm(6,6)
      dimension p(6)
c
c set map to the identity
      call ident(fa,fm)
c
c set up polynomial
      fa(12)=p(2)
      fa(17)=-p(1)
      fa(21)=p(4)
      fa(24)=-p(3)
      fa(27)=-p(5)/2.d0
c
      return
      end
c
***********************************************************************
c
      subroutine dwnd(p)
c  this is a subroutine for windowing tracking results from a dynamic map
c  Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'rays.inc'
      dimension p(6)
c
      iturn=iturn+1
      do 100 k=1,nrays
      if (istat(k).ne.0) goto 100
c  examine particle location
      ax=abs(zblock(k,1))
      apx=abs(zblock(k,2))
      ay=abs(zblock(k,3))
      apy=abs(zblock(k,4))
      at=abs(zblock(k,5))
      apt=abs(zblock(k,6))
      if (ax.gt.p(1) .or. apx.gt.p(2)
     # .or. ay.gt.p(3) .or. apy.gt.p(4)
     # .or. at.gt.p(5) .or. apt.gt.p(6))
     # then
c  particle outside window
      istat(k)=iturn
      nlost=nlost+1
      ihist(nlost,1)=iturn
      ihist(nlost,2)=k
      endif
 100  continue
      return
      end
c
*********************************************************************
c
      subroutine eapt(p)
c  this is a subroutine for an elliptic aperture
c  Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'rays.inc'
      dimension p(6)
c
      smas=p(2)
      raxs=p(3)
      iturn=iturn+1
      do 100 k=1,nrays
      if (istat(k).ne.0) goto 100
      x=zblock(k,1)
      y=zblock(k,3)
      x2=x*x
      y2=y*y
      test=x2+raxs*y2
c  examine particle location
      if (test.ge.smas) goto 10
c  particle within aperture
      goto 100
  10  continue
c  particle outside aperture
      istat(k)=iturn
      nlost=nlost+1
      ihist(nlost,1)=iturn
      ihist(nlost,2)=k
 100  continue
      return
      end
c
*********************************************************************
c
      subroutine ftm(p,fa,fm)
c  this subroutine filters transfer maps
c  Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      dimension p(6),fa(monoms),fm(6,6)
      dimension ta(monoms),tm(6,6)
c
      ifile=nint(p(1))
      nopt=nint(p(2))
      nskp=nint(p(3))
      kind=nint(p(4))
      mpitmp=mpi
      mpi=ifile
      call mapin(nopt,nskp,ta,tm)
      mpi=mpitmp
c  determine procedure
      if(kind.eq.0) goto 5
      if(kind.eq.1) goto 15
c  procedure for normal filter
   5  continue
      do 10 i=1,monoms
      if(ta(i).eq.0.) fa(i)=0.
  10  continue
      return
c  procedure for a reversed filter
  15  continue
      do 20 i=1,monoms
      if(ta(i).gt.0.) fa(i)=0.
  20  continue
      return
      end
 
c
***********************************************************************
c
      subroutine ident(h,xmh)
c Written by D. Douglas, ca 1982
      include 'impli.inc'
      include 'param.inc'
      dimension h(monoms),xmh(6,6)
      do 10 i=1,6
      do 10 j=1,6
   10 xmh(i,j)=0.
      do 15 i=1,6
   15 xmh(i,i)=1.
      do 20 i=1,monoms
   20 h(i)=0.
      return
      end
c
***********************************************************************
c
      subroutine inv(h,mh)
c  Returns the inverse of the map represented by h, mh.
c  Written by Liam Healy, April 30, 1985.
c     implicit none
      include 'param.inc'
      include 'symp.inc'
      double precision h(*),mh(6,6)
      double precision temp(6,6)
      integer iord,ind,j,k,l
c----Routine----
      call mclear(temp)
c  Reverse the factorization:
      iord=0
      call revf(iord,h,mh)
c  Invert matrix:
      do 120 j=1,6
      do 120 k=1,6
      do 120 l=1,6
  120  temp(j,k)=temp(j,k)+mh(l,j)*jm(l,k)
c  Clear out mh to fill it with its inverse:
      call mclear(mh)
c  Calculate inverse of matrix mh:
      do 140 j=1,6
      do 140 l=1,6
      do 140 k=1,6
  140  mh(j,k)=mh(j,k)+jm(l,j)*temp(l,k)
c  Invert polynomials:
      do 100 ind=1,monoms
  100   h(ind)=-h(ind)
      return
      end
c
**********************************************************************
c
      subroutine mask(wipe,h,mh)
c
c  Sets higher order monomial coefficients to zero as specified by wipe
c  Written by Liam Healy, Spring 1985
c
c     implicit none
c  Variables
c  wipe = array of flags: if nth element is less than 0.5, nth order is
c         set to the identity map
      double precision wipe(*)
c  h, mh = polynomial and matrix (input and output)
      double precision h(*),mh(6,6)
c  ord = order of polynomial
      integer ord
c
      include 'maxcat.inc'
      include 'lims.inc'
c
c----Routine----
      if (wipe(1).le.0.5) then
        do 100 i=1,6
          do 100 j=1,6
  100       mh(i,j)=0.
        do 120 i=1,6
  120     mh(i,i)=1.
      endif
      do 200 ord=2,ordcat
        if (wipe(ord).le.0.5) then
          do 140 i=bottom(ord),top(ord)
  140       h(i)=0.
        endif
  200 continue
      return
      end
c
***********************************************************************
c
      subroutine mtran(mh)
c  Takes the transpose of the matrix mh.
c  Written by Liam Healy, April 16, 1985.
c     implicit none
      double precision mh(6,6),hold
      integer i,j
c
c----Routine----
      do 100 i=1,6
        do 100 j=1,i-1
          hold=mh(j,i)
          mh(j,i)=mh(i,j)
          mh(i,j)=hold
  100 continue
      return
      end
c
***********************************************************************
c
      subroutine rapt(p)
c  this is a subroutine for a rectangular aperture
c  Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'rays.inc'
      dimension p(6)
c
      xmin=p(1)
      xmax=p(2)
      ymin=p(3)
      ymax=p(4)
      iturn=iturn+1
      do 100 k=1,nrays
      if (istat(k).ne.0) goto 100
      x=zblock(k,1)
      y=zblock(k,3)
c  examine particle location
      if (x.le.xmin) goto 10
      if (x.ge.xmax) goto 10
      if (y.le.ymin) goto 10
      if (y.ge.ymax) goto 10
c  particle within aperture
      goto 100
  10  continue
c  particle outside aperture
      istat(k)=iturn
      nlost=nlost+1
      ihist(nlost,1)=iturn
      ihist(nlost,2)=k
 100  continue
      return
      end
c
***********************************************************************
c
      subroutine rev(h,mh)
c This is a subroutine for reversing a map.
c Written by Alex Dragt on Friday, 13 Sept 1985.
c     implicit none
      include 'param.inc'
      include 'expon.inc'
c
      double precision h(*),mh(6,6)
      double precision temp(6,6)
      double precision r(6,6)
      save r
c Define the reversing matrix r by a set of data statements:
      data (r(1,j),j=1,6)/ 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/
      data (r(2,j),j=1,6)/ 0.d0,-1.d0, 0.d0, 0.d0, 0.d0, 0.d0/
      data (r(3,j),j=1,6)/ 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0/
      data (r(4,j),j=1,6)/ 0.d0, 0.d0, 0.d0,-1.d0, 0.d0, 0.d0/
      data (r(5,j),j=1,6)/ 0.d0, 0.d0, 0.d0, 0.d0,-1.d0, 0.d0/
      data (r(6,j),j=1,6)/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/
c----Routine----
c Compute inverse map:
      call inv(h,mh)
c Reverse matrix portion of map:
      call mmult (r,mh,temp)
      call mmult (temp,r,mh)
c Reverse polynomial portion of map:
c For now, work only with n in the interval [28,209].
c Change range of n later when f1's are implemented.
      do 10 n=28,209
        n2 = expon(2,n)
        n4 = expon(4,n)
        n5 = expon(5,n)
        m = n2+n4+n5
        mm2 = mod(m,2)
        if (mm2.eq.0) h(n)=-h(n)
   10 continue
      return
      end
c
***********************************************************************
c
      subroutine revf(iord,h,mh)
c  Changes the order of factorization.
c  If iord=0, it is assumed that the map specified by h,mh is in the
c  standard MARYLIE order exp(:f2:)exp(:f3:)exp(:f4). This
c  routine will then return g3,g4 corresponding to the factorization
c  exp(:g4:)exp(:g3:)exp(:f2:).
c  If iord=1, it is assumed that the map specified by h,mh is in
c  reversed order.  This routine then returns the map in the standard
c  MARYLIE order.
c  Written by Liam Healy, April 16, 1985.
c
c----Variables----
c     implicit none
      include 'param.inc'
      double precision h(*),mh(6,6)
      integer iord
c
      double precision wipe(6)
      double precision ftmp1(monoms),fmtmp1(6,6)
      double precision ftmp2(monoms),fmtmp2(6,6)
      integer colme(6,0:6)
c
c----Routine----
      if(iord.gt.0) go to 200
c----Procedure for going from MARYLIE order to reversed order:
      do 100 i=3,4
        call xform(h,i,mh,i-3,ftmp1)
  100 continue
      do 120 ind=1,monoms
  120   h(ind)=ftmp1(ind)
      return
c----Procedure for going from reversed order to MARYLIE order:
  200 continue
      call mapmap(h,mh,ftmp1,fmtmp1)
      call mapmap(h,mh,ftmp2,fmtmp2)
      wipe(1)=0.
      wipe(2)=1.
      wipe(3)=1.
      wipe(4)=1.
      call mask(wipe,ftmp1,fmtmp1)
      wipe(1)=1.
      wipe(2)=0.
      wipe(3)=0.
      wipe(4)=0.
      call mask(wipe,ftmp2,fmtmp2)
      call concat(ftmp1,fmtmp1,ftmp2,fmtmp2,h,mh)
      return
      end
c
***********************************************************************
c
      subroutine strget(kynd,nmap,fa,fm)
c This is a subroutine for storing and getting maps.
c A total of 5 maps can be stored and retrieved.
c The incoming and outgoing maps are represented by fa,fm.
c In the store mode, the maps are stored in sf1,sm1 to sf5,sm5.
c In the retrieve mode, the map is gotten from sf1,sm1 to sf5,sm5.
c The maps sf1,sm1 to sf5,sm5 are stored in the block common stmap.
c Written by Alex Dragt, Spring 1987. Modified 10/13/88 AJD.
c
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'stmap.inc'
c
c Calling arrays
      dimension fa(monoms),fm(6,6)
      character*3 kynd
c
      if (kynd.eq.'gtm') goto 100
c Procedure for storing maps:
      write(jof,400) nmap
  400 format(1x,'map stored in location',2x,i2)
      goto(10,20,30,40,50),nmap
   10 call mapmap(fa,fm,sf1,sm1)
      return
   20 call mapmap(fa,fm,sf2,sm2)
      return
   30 call mapmap(fa,fm,sf3,sm3)
      return
   40 call mapmap(fa,fm,sf4,sm4)
      return
   50 call mapmap(fa,fm,sf5,sm5)
      return
c
c Procedure for getting maps:
  100 continue
      write(jof,450) nmap
  450 format(1x,'map gotten from location',2x,i2)
      goto(110,120,130,140,150),nmap
  110 call mapmap(sf1,sm1,fa,fm)
      return
  120 call mapmap(sf2,sm2,fa,fm)
      return
  130 call mapmap(sf3,sm3,fa,fm)
      return
  140 call mapmap(sf4,sm4,fa,fm)
      return
  150 call mapmap(sf5,sm5,fa,fm)
      return
c
      end
c
*********************************************************************
c
      subroutine sympl(iopt,kind,fa,fm)
c  This is a symplectification subroutine
c Written by Alex Dragt, Spring 1987.
      include 'impli.inc'
      include 'param.inc'
      dimension fa(monoms),fm(6,6)
c
      if (kind.eq.1) call sympl1(iopt,fm)
      if (kind.eq.2) call sympl2(iopt,fm)
      if (kind.eq.3) call sympl3(iopt,fm)
c
      return
      end
c
***********************************************************************
c
      subroutine tpol(p,fa,fm)
c This is a subroutine for computing a quadratic polynomial
c described in terms of twiss parameters.
c Written by Alex Dragt, 21 December 1990
      include 'impli.inc'
      include 'param.inc'
      dimension p(6)
      dimension fa(monoms),fm(6,6)
c-----
c set map to the identity
      call ident(fa,fm)
c set up twiss parameters
      ax=p(1)
      bx=p(2)
      ay=p(3)
      by=p(4)
      at=p(5)
      bt=p(6)
      gx=(1.+ax*ax)/bx
      gy=(1.+ay*ay)/by
      gt=(1.+at*at)/bt
      if (bx .lt. 0.) then
      ax=0.
      bx=0.
      gx=0.
      endif
      if (by .lt. 0.) then
      ay=0.
      by=0.
      gy=0.
      endif
      if (bt .lt. 0.) then
      at=0.
      bt=0.
      gt=0.
      endif
c set up polynomial
      fa(7)=gx
      fa(8)=2.d0*ax
      fa(13)=bx
      fa(18)=gy
      fa(19)=2.d0*ay
      fa(22)=by
      fa(25)=gt
      fa(26)=2.d0*at
      fa(27)=bt
      return
      end
c
c**********************************************************************
c
      subroutine wnd(pp)
c  this is a subroutine for windowing tracking results in
c  all six variables.
c  Written by Alex Dragt, 12 January 1991.
      include 'impli.inc'
      include 'rays.inc'
      dimension pp(6)
c
c set up control parameters
      iplane=nint(pp(1))
      if ((iplane .lt. 1) .or. (iplane .gt. 3)) then
      write(6,*) 'iplane out of range in wnd'
      return
      endif
      icon=2*(iplane - 1)
      qmin=pp(2)
      qmax=pp(3)
      pmin=pp(4)
      pmax=pp(5)
c
      iturn=iturn+1
      do 100 k=1,nrays
      if (istat(k).ne.0) goto 100
c  examine particle
      q=zblock(k,1+icon)
      p=zblock(k,2+icon)
      if   (q.lt.qmin .or. q.gt.qmax
     # .or. p.lt.pmin .or. p.gt.pmax)
     # then
c  particle outside window
      istat(k)=iturn
      nlost=nlost+1
      ihist(nlost,1)=iturn
      ihist(nlost,2)=k
      endif
 100  continue
      return
      end
c
c end of file
