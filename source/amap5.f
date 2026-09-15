      subroutine amap5(p,fa,fm)
c subroutine for applying a map to a function or a set of moments
c THIS ROUTINE IS MERGED INTO AMAP (anal.f)
c
      implicit double precision (a-h,o-z)
      dimension p(6)
      dimension fa(923),fm(6,6)
ctm2014 dimension changed per gfortran error message.
      dimension ga(923),gm(6,6)
      dimension ha(923),hm(6,6)
      dimension t1a(923)
      dimension t2a(923)
      include 'files.inc'
      character*3 kynd
c
c set up control indices
      mode=nint(p(1))
      isend=nint(p(2))
      ifile=nint(p(3))
      nopt=nint(p(4))
      nskip=nint(p(5))
      nmpo=nint(p(6))
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
c test for mode
c compute ha = D*ga
c clear arrays
      do 10 i=1,923
      t1a(i)=0.d0
   10 ha(i)=0.d0
c perform calculation
      do 20 i=7,27
      t1a(i)=1.d0
      call fxfrm5(fa,fm,t1a,t2a)
      t1a(i)=0.d0
      do 30 j=1,923
   30 ha(i)=ha(i)+t2a(j)*ga(j)
      if(ha(i).ne.0) write(5,*) i,ha(i)
   20 continue
c compute squares of rms emittances
      xemit2=ha(7)*ha(13)-ha(8)*ha(8)
      yemit2=ha(18)*ha(22)-ha(19)*ha(19)
      if (isend.eq.1.or.isend.eq.3) then
      write (6,*) 'xemit2=',xemit2
      write (6,*) 'yemit2=',yemit2
      endif
      if (isend.eq.2.or.isend.eq.3) then
      write (jodf,*) 'xemit2=',xemit2
      write (jodf,*) 'yemit2=',yemit2
      endif
c write out result ha if nmpo > 0
      mpot=mpo
      mpo=nmpo
      if (nmpo.gt.0)  write (6,*) 'mapout5 not implemented'
c      call mapout5(0,ha,hm)
      mpo=mpot
c write results into the array ucalc
c first clear the array
c      do 40 i=1,250
c   40 ucalc(i)=0.d0
c complete the task
c      nuvar=211
c      do 50 i=1,209
c   50 ucalc(i)=ha(i)
c      if (mode.eq.2.or.mode.eq.3) then
c      ucalc(210)=xemit2
c      ucalc(211)=yemit2
c      endif
      return
      end
c
      subroutine fxfrm5(ga,gm,fa,ha)
c this subroutine transforms a function f.
c exp(:ga2:)exp(:ga3:)...exp(:ga6:) is applied to f.
      implicit double precision (a-h,o-z)
      dimension ga(923),fa(923),ha(923),t1a(923),t2a(923)
      dimension gm(6,6)
      do 100 i=1,923
        ha(i) = 0.0d0
        t1a(i) = fa(i)
        t2a(i) = 0.0d0
  100 continue
      do 200 ideg = 4,3,-1
        call exphf5(ga,ideg,t1a,6,t2a)
        do 222 i = 1,923
          t1a(i) = t2a(i)
  222   continue
  200 continue
      call xform5(t1a,6,gm,0,ha)
      do 400 ideg = 5,2,-1
        call xform5(t1a,ideg,gm,1,ha)
  400 continue
      return
      end
c
      subroutine exphf5(h,ideg,f,maxf,trf)
c  Applies Exp(:h:) on polynomial f.
c  h is a polynomial of degree ideg.
c  f has terms from 1 thru maxf.
c  The result is trf, which has terms 1-maxf.
c
c   Written By F. Neri 9/26/86.
c
      double precision f(923),h(923),trf(923)
      double precision tmpf1(923),tmpf2(923)
c  maxf has to be .le. 6.
      integer maxf
      integer maxpow, ifact, maxord
      include 'lims.inc'
      do 1 i=1,top(maxf)
        tmpf1(i) = f(i)
   1  continue
      do 3 i = 1,top(maxf)
        trf(i) = f(i)
    3 continue
      maxpow = int((maxf-1)/(ideg-2))
      ifact = 1
      do 10 n=1,maxpow
        ifact =  ifact * (-n)
        maxord = int(maxf - (ideg-2) )
        do 11 i=1,top(maxf)
          tmpf2(i) = 0.d0
   11   continue
        do 20 iord=maxord,1,-1
          call pbkt5(tmpf1,iord,h,ideg,tmpf2)
          call pmadd5(tmpf2,iord+ideg-2,(1.d0/ifact),trf)
  20    continue
        do 30 i=1,top(maxf)
          tmpf1(i) = tmpf2(i)
   30   continue
  10  continue
      return
      end
c
      subroutine mapin5(nopt,nskp,h,mh)
c read nonzero matrix elements and monomials from file unit mpi
      implicit double precision (a-h,o-z)
      double precision h,mh
      include 'files.inc'
      dimension mh(6,6),h(923)
c
      if(nopt.eq.0)goto 5
      rewind mpi
      write(jof,210)mpi
    5 continue
c     initialize arrays:
      do 10 j=1,923
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
      if(k.eq.923)goto 90
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
      if(k.eq.923)goto 200
      goto 180
  200 continue
      write(6,205) mpi,nskp
  205 format(1x,'map read in from file ',i3,'; ',i3,
     #' record(s) skipped')
  210 format(1x,'file unit ',i2,' rewound')
      return
      end
c
      subroutine pmadd5(f,n,coeff,h)
      implicit double precision (a-h,o-z)
      dimension f(923),h(923)
	include 'lims.inc'
      if(coeff.eq.1.d0) goto 20
      do 10 i=bottom(n),top(n)
        h(i) = h(i) + f(i)*coeff
 10   continue
      return
 20   continue
      do 30 i =bottom(n),top(n)
        h(i) = h(i) + f(i)
 30   continue
      return
      end
c
c
      subroutine pbkt5(f,ordf,g,ordg,pb)
c  Calculates the Poisson Bracket of the order ordf part of f with
c  the order ordg part of g, leaving the result in pb.
c
      include 'param.inc'
c----Variables----
c  f,g,pb = two arrays and their Poisson bracket
      double precision f(*),g(*),pb(*)
c  ordf,ordg = order desired from f and g
      integer ordf,ordg
c  indf,indg,indpb = index in f,g and pb
      integer indf,indg,indpb
c  pr = array of exponents for the product of f(indf) and g(indg)
c  prd =copy of pr
c  prx, prp = exponents of particular x, p variables, from pr
      integer pr(6),prd(6),prx,prp
c  xvb, pvb, vbl = x variable, p variable (pvb=xvb+1), variable number
      integer xvb,pvb,vbl
c  expon = table of exponents
      include 'expon5.inc'
c  bottom, top = lowest and highest monomial number for each order
      include 'lims.inc'
      include 'prodex5.inc'
      integer vbf,vbg
      integer conj(6),sigj(6)
      data conj /2,1,4,3,6,5/
      data sigj /1,-1,1,-1,1,-1/
	external iprod
	integer iprod
c
c----Routine----
c  initialize array pb
      do 80 indpb=bottom(ordf+ordg-2),top(ordf+ordg-2)
   80   pb(indpb)=0.
c  pick individual indf and indg, find what element of pb it affects,
c  and calculate the new value.  Loop for all indeces in the specified orders.
      do 160 vbf = 1,6
        vbg = conj(vbf)
        sign = sigj(vbf)
        do 100 indf=bottom(ordf-1),top(ordf-1)
          if(f(prodex5(vbf,indf)).eq.0.0d0 ) goto 100
          do 110 indg=bottom(ordg-1),top(ordg-1)
            indpb = iprod(indg,indf)
            pb(indpb) = pb(indpb) +
     #      sign * f(prodex5(vbf,indf)) * g(prodex5(vbg,indg))
     #      * (expon5(vbf,indf)+1) * (expon5(vbg,indg)+1)
  110     continue
  100   continue
  160 continue
      return
      end
c
      subroutine xform5(f,no,m,nu,l)
ctm04   nu is a dummy argument for consistancy with the calling sequence
ctm04   I don't know what it is supposed to do.
c
c     transforms arguements of polynomial of degree n
c     by the lie transformation whose
c     matrix representation is m.  the coefficients
c     of the resultant polynomial are stored in the
c     the array l thus contains the coefficients of f(m*z).
c
      implicit double precision (a-h,o-z)
      double precision l,m,f
      dimension f(923),l(923)
      double precision temp(923)
      dimension m(6,6)
	  include 'expon5.inc'
	  include 'prodex5.inc'
	  include 'vblist5.inc'
	  include 'lims.inc'
c
c     initialise arrays
c
      do 10 kp=1,top(no)
      l(kp)=0.0d0
      temp(kp) = 0.0d0
   10 continue
      do 100 n= bottom(no), top(no)
        if(f(n).eq.0.0d0) goto 100
        do 101 kp=1, top(no-1)
          temp(kp) = 0.0d0
  101   continue
        do 102 k=1,6
          if(m(vblist5(1,n),k).eq.0.0d0) goto 102
          temp(k) = f(n)*m(vblist5(1,n),k)
  102    continue
        do 110 ior=1,no-1
          k1 = vblist5(ior+1,n)
          n1 = bottom(ior)
          do 112 k=1,6
            if(m(k1,k).eq.0.0d0) goto 112
            xm = m(k1,k)
c
            do 111 nn=n1, top(ior)
            temp(prodex5(k,nn)) = temp(prodex5(k,nn)) + temp(nn)*xm
  111       continue
  112     continue
  110   continue
  100 continue
      do 200 nn= bottom(no), top(no)
        l(nn) = l(nn)+temp(nn)
  200 continue
      return
      end
c
      subroutine ident5(h,xmh)
      include 'impli.inc'
      include 'param.inc'
      dimension h(monoms5),xmh(6,6)
      do 10 i=1,6
      do 10 j=1,6
   10 xmh(i,j)=0.
      do 15 i=1,6
   15 xmh(i,i)=1.
      do 20 i=1,monoms5
   20 h(i)=0.
      return
      end
c
      subroutine mapmap5(rh,rmh,th,tmh)
c  Written by Rob Ryne, ca 1982
      include 'impli.inc'
      include 'param.inc'
      dimension rh(monoms5),th(monoms5),rmh(6,6),tmh(6,6)
      do 10 i=1,6
      do 10 j=1,6
   10 tmh(i,j)=rmh(i,j)
      do 20 i=1,monoms5
   20 th(i)=rh(i)
      return
      end
c
      subroutine mapout5(nopt,h,mh)
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
      dimension mh(6,6),h(monoms5)
c
c check to see if file assignment makes sense
      if(mpo .le. 0) return
c
      nout=mpo
c
c  the parameter nopt is not currently used
c
c  go to end of file
c
   15 read(nout,*,end=17,err=19)dummy
      goto 15
   17 write(jof,18) nout
   18 format(' adding current map to file on unit',i3)
      go to 25
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
      do 40 k=1,monoms5
   40 if(h(k).ne.0.)write(nout,*)k,h(k)
      if(h(monoms5).eq.0.)write(nout,*)monoms5,h(monoms5)
c      write(jof,100) nout
c  100 format(1x,'map written on file ',i2)
      return
      end
c
      subroutine strget5(kynd,nmap,fa,fm)
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
      dimension fa(monoms5),fm(6,6)
      character*3 kynd
c
      if (kynd.eq.'gtm') goto 100
c Procedure for storing maps:
      write(jof,400) nmap
  400 format(1x,'map stored in location',2x,i2)
      goto(10,20,30,40,50),nmap
   10 call mapmap5(fa,fm,sf1,sm1)
      return
   20 call mapmap5(fa,fm,sf2,sm2)
      return
   30 call mapmap5(fa,fm,sf3,sm3)
      return
   40 call mapmap5(fa,fm,sf4,sm4)
      return
   50 call mapmap5(fa,fm,sf5,sm5)
      return
c
c Procedure for getting maps:
  100 continue
      write(jof,450) nmap
  450 format(1x,'map gotten from location',2x,i2)
      goto(110,120,130,140,150),nmap
  110 call mapmap5(sf1,sm1,fa,fm)
      return
  120 call mapmap5(sf2,sm2,fa,fm)
      return
  130 call mapmap5(sf3,sm3,fa,fm)
      return
  140 call mapmap5(sf4,sm4,fa,fm)
      return
  150 call mapmap5(sf5,sm5,fa,fm)
      return
c
      end
