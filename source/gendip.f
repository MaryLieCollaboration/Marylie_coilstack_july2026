************************************************************************
* header              GENDIP (GENMAP for a "parallel" face dipole      *
*                     magnet with soft fringe fields                   *
*  All routines needed for this special GENMAP                         *
************************************************************************
c
      subroutine bderivs(z,y,b)
c  This routine computes b(z) = Int(-Inf to z) By(z') dz
c  and the second and fourth derivative b2(z) and b4(z) for
c  (a dipole magnet.)
      include 'impli.inc'
      include 'parm.inc'
      include 'param.inc'
      include 'dip.inc'
      include 'vecpot.inc'   !cryne 18 March 2011
c----------------------------------------
      double precision z, y(*), b(0:6,0:6)
c----------------------------------------
      external amyf,f1,f2,f4,f6
c----------------------------------------
cryne 9 April 2011 included factor of pi
      idproc = 0
      pi=2.d0*asin(1.d0)
      epsz2=gap/(sl*pi)
ccc   epsz2 = gap/sl
      bb = (By*sl)
c     write(6,*)'gap,epsz2,bb=',gap,epsz2,bb
      do 10 i=0,6
        do 10 j=0,6
          b(i,j) = 0.0d0
 10   continue
c
c      b = (amyf(z-za,epsz2) - amyf(z-zb,epsz2)) * bb
      b(0,0)= (f1(z-za,epsz2) -f1(z-zb,epsz2)) * bb
      b(0,1)= (f2(z-za,epsz2) -f2(z-zb,epsz2)) * bb
      b(0,3)= (f4(z-za,epsz2) -f4(z-zb,epsz2)) * bb
c      b6= (f6(z-za,epsz2)-f6(z-zb,epsz2)) * bb
      Ax(0)=0.d0  !hardwired for the case y=0.d0
      Ay(0)=0.d0  !hardwired for the case y=0.d0
      Ax(1)=0.d0  !hardwired for the case y=0.d0
      Ay(1)=0.d0  !hardwired for the case y=0.d0
      Ax(3)=0.d0  !hardwired for the case y=0.d0
      Ay(3)=0.d0  !hardwired for the case y=0.d0
      Az(1)=-b(0,0)/brho !hardwired...
      Az(3)=0.d0   !hardwired...
c     write(38,*),z,b(0,0)
c     write(6,*) z, b(0,0),' = z,b'
      return
      end
c
************************************************************************
      function amyf(y,eps)
      double precision amyf,y,eps
      if (y/eps.lt.0.d0) then
        amyf = eps*Log(1 + Exp(2*y/eps))/2
      else
        amyf = y + eps*Log(Exp(-2*y/eps)+1)/2
      endif
      return
      end
c
************************************************************************
      function f1(y,eps)
      double precision f1,y,eps
        f1 = tanh(y/eps)/2.d0 +.5d0
      return
      end
c
************************************************************************
      function f2(y,eps)
      double precision f2,y,eps
      if (dabs(y/eps).lt.10.d0) then
        f2 = 2*Exp(2*y/eps)/(eps*(1 + Exp(2*y/eps))**2)
      else
        f2 = 0.d0
      endif
      return
      end
c
************************************************************************
      function f4(y,eps)
      double precision f4,y,eps
      if (dabs(y/eps).lt.10.d0) then
        f4 = 8*Exp(2*y/eps)*
     &      (1 - 4*Exp(2*y/eps) + Exp(4*y/eps))/
     &   (eps**3*(1 + Exp(2*y/eps))**4)
      else
        f4 = 0.0d0
      endif
      return
      end
c
************************************************************************
       function f6(y,eps)
       double precision f6,y,eps
      if (dabs(y/eps).lt.10.d0) then
         f6 = 32*Exp(2*y/eps)*
     &        (1 - 26*Exp(2*y/eps) +
     &            66*Exp(4*y/eps) -
     &            26*Exp(6*y/eps) + Exp(8*y/eps))/
     &   (eps**5*(1 + Exp(2*y/eps))**6)
      else
        f6 = 0.d0
      endif
      return
      end
c
************************************************************************
      subroutine gendip(p,fa,fm,reftmp)
 
c This is a subroutine for computing the map for a soft edged dipole
c magnet ( F. Neri 5/16/89 ).
c The routine is based on Rob Ryne original gendip, but all the code has been
c rewritten.
c
ctm 1 mar 2012: This version is a modified gendip5 from ML/I made operational
c               in Nov 2011 by ryne, with parameters changed a bit by ctm for ML3.
      include 'impli.inc'
      include 'parm.inc'
      include 'param.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'files.inc'
      include 'dip.inc'
      include 'pie.inc'
      include 'bfield.inc' !used for diagnostic output of field profile
c
c  calling arrays
      dimension p(*)   !cryne 3/22/11 changed from p(6) to p(*)
      dimension fa(monoms), fm(6,6)
      dimension reftmp(6)
c
c  local arrays
      dimension pb(6)
      dimension y(monoms+15)
      dimension yjunk(monoms+15)
      common/fevaldiag/jprintfeval
      common/hmltn3type/nerihmltn3 ! ryne Nov 4, 2011
      logical ltty, ldsk
c
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
c#    za = p(1)
c#    zb = p(2)
c#    ns = 10000
c#    By = p(3)
c#    phi = p(4) * pi180
c#    s = sin(phi)
c#    gap = p(5)
c#    zi  = 0.0d0
c#    zf  = p(6)
      idproc = 0
      ltty = .false.
      ldsk = .false.
      jprintfeval=lunf ! used for printing reftraj in routine feval
      write(6,211) (p(i),i=1,9)
 211  format('P:',5f15.6)
c
      phi1=p(2)*pi180
      phi2=-p(3)*pi180 !cryne 9 April 2011 flipped sign !!!!!!!!!!!!!!!!!!!!!
      gap=p(5)
      ngap=nint(p(6))
      msg=nint(p(7))
      if(msg.eq.1) ltty = .true.
      if(msg.eq.2) ldsk = .true.
      if(msg.gt.2) then
         ltty = .true.
         ldsk = .true.
         lunf = msg
      endif
      lupro=nint(p(8))
      ns=nint(p(9))
c
      if(p(1).eq.-99999.d0 .and. p(4).eq.-99999.d0)then
        write(6,*)'error in gendip: user entered neither L nor B'
        call myexit
      endif
      if(p(1).ne.-99999.d0)then
        By=(sin(phi1)+sin(-phi2))*brho/p(1)
          write(6,*)'(gendip) user entered L'
          write(6,*)'computed value of By=',By
      else
        By=p(4)
        if(msg.eq.1)then
          write(6,*)'(gendip) user entered B'
        endif
        if(By.eq.0.d0)then
          write(6,*)'error: By=0 in routine gendip'
          call myexit
        endif
      endif
c
      zbminusza=-(sin(phi2)-sin(phi1))*brho/By
      zi=0.d0
      efoldlen=gap/sl     !because Filippo passes epsz2=gap/sl to routines
      za=zi+ngap*efoldlen
      zb=za+zbminusza
      zf=zb+ngap*efoldlen
c
      ltty = .true.
      if(ltty)then
        write(6,*) 'gendip parameters:'
        write(6,*) 'phi1,phi2=',phi1,phi2
        write(6,*) 'By=',By
        write(6,*) 'gap,ngap=',gap,ngap
        write(6,*) 'ns=',ns
        write(6,2223) za,zb,zb-za
        write(6,2224) zi,zf,zf-zi
 2223   format('za,zb,zb-za=',3(1pd22.15,1x))
 2224   format('zi,zf,zf-zi=',3(1pd22.15,1x))
      endif
c
      if(lupro.gt.0)then
      do n=0,ns
        zval=zi+(1.d0*n)*(zf-zi)/ns
        call bderivs(zval,yjunk,b(0,0))  !yjunk not used
        write(lupro,*)zval,b(0,0)
      enddo
      endif
c
      jtty = 1
      jdsk = 1
c
c      if((isend.eq.1).or.(isend.eq.3)) jtty = 1
c      if((isend.eq.2).or.(isend.eq.3)) jdsk = 1
c
      h=(zf-zi)/float(ns)
c
c     call VAX system routine for timing report
c
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=0.d0
cryne 3/22/11 does Filippo use "s" anywhere else? It is in dip.inc      y(2)= s
      s=sin(phi1)
      y(2)=s
      y(3)=0.d0
      y(4)=0.d0
      y(5)=0.d0
      y(6)=-1.d0/beta
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c
c  initialize map to the identity map:
crynechannellwalstrom      ne=224
      ne=monoms+15
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c
      if(jprintfeval.ne.0)then
        write(jprintfeval,'(99(1pe12.5,1x))')zi,y(1:6)
      endif
c  do the computation:
      t=zi
      iflag = 3
      nerihmltn3=1  ! Nov 4 2011, ryne this flag is used by feval
      call adam11(h,ns,'start',t,y)
      if(ltty)then
        write(6,*)'y(1),y(2)=',y(1),y(2)
        write(6,*)'y(3),y(4)=',y(3),y(4)
        write(6,*)'y(5),y(6)=',y(5),y(6)
        write(6,*)'y(7),y(8)=',y(7),y(8)
      endif
      call chkdet(y)
      call putmap(y,fa,fm)
      reftmp(1:6)=y(1:6)
cryne 3/22/11 why did Filippo set this to -y(2)?      s1 = -y(2)
      s1 =y(2)
      phi1=asin(s1)
      phi1deg=phi1/pi180
      write(jof,991)phi1deg
      write(jodf,991)phi1deg
991   format('  Final angle is ',1pd22.15,' degrees')
c
c     call VAX system routine for timing report
c
c567  format(' GENDIP integration time = ',f12.2,' sec.')
      nerihmltn3=0  ! Nov 4 2011, ryne reset to default
      return ! Nov 4 2011, ryne added return statement
      end
c
**********************************************************************
c
      subroutine hmltn3neri(t,y,h)
c  this routine is used to specify h(z) for a dipole magnet.
c  Written by F. Neri, 5/16/89.
c  Modified 6/3/89 to handle different gauges.
c  The design orbit is assumed to be in the Y = 0 plane.
c  B Field derivatives on the design orbit are provided by the routine
c  BDERIVS as a function of z,and x = y(1), bderivs is called as
c      call bderivs(z,y,b)
c  The derivatives are stored in the array b(0:*,0:*), defined in the
c  include file bfield.inc, and passed to other parts of the program
c  in the common/bfield/b(0:*,0:*)
c  The array b is arranged so that B(0,0) is By, B(1,0) is d(By)/(dX),
c  b(1,1) is d2(By)/(dX dz), etc.
c
      include 'impli.inc'
      include 'parm.inc'
      include 'param.inc'
crynechannellwalstrom      parameter (itop=209,iplus=224)
      parameter (itop=monoms,iplus=itop+15)
      include 'dip.inc'
      include 'bfield.inc'
c
      dimension h(monoms),y(*)
c
      dimension X(0:itop),YY(0:itop),P1(0:itop),P2(0:itop)
      dimension A(0:12)
c  begin calculation
c
c  compute gradients
      call bderivs(t,y,b)
c
c  initialization
      do 10 i=1,monoms
   10 h(i)=0.d0
c
c cccc
c   Slow method to produce hamiltonian: use polynomial
c   expansion of square root: when this method works we will use it to
c   to check the "hardwired" version.
c cccc
c   Coefficients of expansion of -SQRT(1+X)+1
      A(0) = 0.d0
      A(1) = -1.d0/2.d0
      do 50 i=2, 6
        A(i) =  A(i-1) * (1.d0/2.d0 - (i-1.d0))/i
   50 continue
c
crynechannellwalstrom      do 60 i=0,209
      do 60 i=0,monoms
        X(i) = 0.d0
   60 continue
      X(6) = -2.d0 / beta
c      X(13) = -1.d0
      X(22) = -1.d0
      X(27) =  1.d0
c
      maxord = 4
      nn = 4
cryne Nov 4      maxord = 6
cryne Nov 4      nn = 6
      cos2 = 1 - y(2)**2
c P1 = px - q Ax
      do 101 i=0,itop
  101 P1(i) = 0.d0
c Note constant term in px - q Ax, coming from design orbit
c y(2) = Px = sin(phi).
      P1(0) = y(2)
c
      P1(2) = 1.d0
      P1(18) = b(0,1)/(2.d0*brho)
      P1(39) = b(1,1)/(2.d0*brho)
      P1(95) = b(2,1)/(4.d0*brho)
      P1(175) = -(b(2,1)+b(0,3))/(24.d0*brho)
c P2 = P1**2
      call mypmult(P1,P1,P2,maxord)
c Zero order term subtracted ( Sum starts at 1 ):
c Divide by cos2
      do 102 i=1,itop
  102 X(i) = (X(i) - P2(i))/cos2
c X = pt**2 - 2/beta*pt - (px -q Ax)**2 - (py)**2
c YY = -Sqrt(1+X)
      call mypoly1(nn,A,X,YY,maxord)
c h = -Az
      h(7) = b(1,0)/(2.d0*brho)
      h(18) = -b(1,0)/(2.d0*brho)
      h(28) = b(2,0)/(6.d0*brho)
      h(39) = -b(2,0)/(6.d0*brho)
      h(84) = b(3,0)/(24.d0*brho)
      h(95) = -b(3,0)/(4.d0*brho)
      h(175) = (b(3,0)+b(1,2))/(24.d0*brho)
crynechannellwalstrom there must be some addition terms h(???) here that are missing
c h = -Sqrt(1+X) - Az
c Zero and first order terms subtracted ( Sum starts at 7 ):
c Scale by cos
crynechannellwalstrom      do 70 i=7, 209
      do 70 i=7, monoms
        h(i) = h(i) + dsqrt(cos2)*YY(i)/sl
   70 continue
c
      return
      end
c
c ******************************************************************
c
c  Aux polynomial routines (really a poor man's DA package).
c  They don't really belong here.
c
c ******************************************************************
c
      subroutine mypoly1(N,A,X,Y,maxord)
      include 'param.inc'
      parameter (MN=6)
      include 'lims.inc'
crynechannellwalstrom      double precision X(0:209),Y(0:209)
      double precision X(0:monoms),Y(0:monoms)
      double precision A(0:N)
c
crynechannellwalstrom      double precision Vect(0:209,MN)
      double precision Vect(0:monoms,MN)
c
c   NO COMMENT
      do 100 i=0,top(maxord)
 100    Y(i) = 0.0d0
      do 200 i=0,top(maxord)
 200    Vect(i,1) = X(i)
      do 300 iord=2,N
        call mypmult(X,Vect(0,iord-1),Vect(0,iord),maxord)
 300  continue
      Y(0) = A(0)
      do 400 iord=1,N
        do 500 i=0,top(maxord)
          Y(i) = Y(i) + Vect(i,iord)*A(iord)
 500    continue
 400  continue
      return
      end
c
************************************************************************
      subroutine mypmult(p1,p2,p3,maxord)
      include 'param.inc'
crynechannellwalstrom      double precision p1(0:209),p2(0:209),p3(0:209)
      double precision p1(0:monoms),p2(0:monoms),p3(0:monoms)
c
      include 'lims.inc'
      if (maxord.lt.0) return
c
      do 10 i=1,top(maxord)
   10 p3(i) = 0.0d0
      p3(0) = p1(0) * p2(0)
      do 100 mord=1,maxord
        call mypmadd(p1(1),mord,p2(0),p3(1))
        call mypmadd(p2(1),mord,p1(0),p3(1))
        do 200 nord1 = 1,mord-1
          nord2 = mord - nord1
          call myproduct(p1(1),nord1,p2(1),nord2,p3(1))
  200   continue
  100 continue
      return
      end
c
************************************************************************
      subroutine mypmadd(f,n,coeff,h)
      implicit double precision (a-h,o-z)
      include 'param.inc'
crynechannellwalstrom      dimension f(209),h(209)
      dimension f(monoms),h(monoms)
      include 'len.inc'
      include 'lims.inc'
      if(coeff.eq.1.d0) goto 20
      do 10 i=len(n-1)+1,len(n)
        h(i) = h(i) + f(i)*coeff
 10   continue
      return
 20   continue
      do 30 i = len(n-1)+1, len(n)
        h(i) = h(i) + f(i)
 30   continue
      return
      end
c
************************************************************************
      subroutine myproduct(a,na,b,nb,c)
      include 'impli.inc'
      include 'param.inc'
      include 'len.inc'
      include 'expon.inc'
      include 'vblist.inc'
crynechannellwalstrom      dimension a(209),b(209),c(209),l(6)
      dimension a(monoms),b(monoms),c(monoms),l(6)
      if(na.eq.1) then
        ia1 = 1
      else
        ia1 = len(na-1)+1
      endif
      if(nb.eq.1) then
        ib1 = 1
      else
        ib1 = len(nb-1)+1
      endif
       do 200 ia=ia1,len(na)
           if(a(ia).eq.0.d0) goto 200
           do 20 ib = ib1,len(nb)
               if(b(ib).eq.0.d0) goto 20
               do 2 m=1,6
                   l(m) = expon(m,ia) +  expon(m,ib)
   2            continue
                n = ndex(l)
                c(n) = c(n) + a(ia)*b(ib)
  20        continue
 200   continue
       return
       end
