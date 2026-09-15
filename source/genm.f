************************************************************************
* header              GENM (General GENMAP Routines)                   *
*  All routines common to the various GENMAP programs                  *
************************************************************************
c
      subroutine adam11(h,ns,nf,t,y)
c  Written by Rob Ryne, Spring 1986, based on a routine of Alex Dragt.
c  This integration routine makes local truncation errors at each 
c  step of order h**11.  That is, it is locally correct through
c  order h**10.  Due to round off errors, its true precision is
c  realized only when more than 64 bits are used.
c
      include 'impli.inc'
      character*6 nf
      parameter (ne=224)
c
      dimension y(ne),yp(ne),yc(ne),f1(ne),f2(ne),f3(ne),f4(ne),
     # f5(ne),f6(ne),f7(ne),f8(ne),f9(ne),f10(ne),f11(ne)
c
      dimension a(10),am(10),b(10),bm(10)
c
      data (a(i),i=1,10)/57281.d0,-583435.d0,2687864.d0,
     # -7394032.d0,13510082.d0,-17283646.d0,16002320.d0,
     # -11271304.d0,9449717.d0,2082753.d0/
      data (b(i),i=1,10)/-2082753.d0,20884811.d0,-94307320.d0,
     # 252618224.d0,-444772162.d0,538363838.d0,-454661776.d0,
     # 265932680.d0,-104995189.d0,30277247.d0/
c
      nsa=ns
      if (nf.eq.'cont') go to 20
c
c  rk start
      iqt=5
      qt=float(iqt)
      hqt=h/qt
      call feval(t,y,f1)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f2)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f3)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f4)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f5)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f6)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f7)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f8)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f9)
      call rk78ii(hqt,iqt,t,y)
      call feval(t,y,f10)
      nsa=ns-9
      hdiv=h/7257600.0d+00
      do 10 i=1,10
      am(i)=hdiv*a(i)
  10  bm(i)=hdiv*b(i)
c
c  Adams routine
c
  20  tint=t
      do 100 i=1,nsa
      do 30 j=1,ne
      yp(j)=y(j)+bm(1)*f1(j)+bm(2)*f2(j)+bm(3)*f3(j)
     # +bm(4)*f4(j)+bm(5)*f5(j)+bm(6)*f6(j)+bm(7)*f7(j)
     # +bm(8)*f8(j)+bm(9)*f9(j)+bm(10)*f10(j)
   30 continue
      call feval(t+h,yp,f11)
      do 40 j=1,ne
      yp(j)=y(j)+am(1)*f2(j)+am(2)*f3(j)+am(3)*f4(j)+am(4)*f5(j)
     x +am(5)*f6(j)+am(6)*f7(j)+am(7)*f8(j)+am(8)*f9(j)+am(9)*f10(j)
  40  yc(j)=yp(j)+am(10)*f11(j)
  41  call feval(t+h,yc,f11)
      do 50 j=1,ne
   50 y(j)=yp(j)+am(10)*f11(j)
      do 60 j=1,ne
      f1(j)=f2(j)
      f2(j)=f3(j)
      f3(j)=f4(j)
      f4(j)=f5(j)
      f5(j)=f6(j)
      f6(j)=f7(j)
      f7(j)=f8(j)
      f8(j)=f9(j)
      f9(j)=f10(j)
   60 f10(j)=f11(j)
      t=tint+i*h
  100 continue
      return
      end
c
************************************************************************
c
      subroutine feval(t,yy,f)
c  Written by Rob Ryne, Spring 1986
c  Modified by F. Neri and A. Dragt 9/29/89
c  Generic Vector Potential Version by F. Neri, 3/28/90
      include 'impli.inc'
      include 'param.inc'
      include 'hmflag.inc'
      include 'parm.inc'
c
      include 'bfield.inc'
      include 'gronax.inc'
      include 'vecpot.inc'
c
      dimension h(monoms),hint(monoms),temp(monoms)
      dimension yy(224),f(224)
      dimension ajsm(6,6)
      common/hmltn3type/nerihmltn3 ! ryne Nov 4, 2011
      save ajsm,h ! ryne Nov 4, 2011
c
c  y(1-6) = given (design) trajectory
c  y(7-42) = matrix
c  y(43-98) = f3
c  y(99-224) = f4
c
c  design trajectory:
c  NB We have assumed the coordinates have been selected so that
c  the design trajectory is yy(i)=0.d0 for i=1 to 6.
c  Not for dipole, see code after ---- call hamltn3() ----
c  F. Neri June 3 1989.
      f(1)=0.d0
      f(2)=0.d0
      f(3)=0.d0
      f(4)=0.d0
      f(5)=0.d0
      f(6)=0.d0
c
c  matrix and polynomials:
c  the expressions for f below are from p. 2742 of
c  J. Math.Phys.,Vol 24,No.12,Dec 1983 (A.J.Dragt and E.Forest)
c
c  select and compute hamiltonian, s matrix, and ajsm=-j*s
c  skip these steps if iflag=0
      if (iflag.eq.0) go to 100
c  otherwise carry out the selection and required computations
      go to (10,20,30,40),iflag
      write(6,*) 'trouble with iflag in subroutine feval'
      call myexit
   10 call hmltn1(h)
      call matify(ajsm,h)
      goto 100
   20 call hmltn2(t,yy,h)
      call matify(ajsm,h)
      go to 100
   30 continue
c-------------------------- c ryne Nov 4, 2011
crynechannelwaltrom: select hmltn3 or hmltn3neri based on variable nerihmltn3
c (which gets set in routine lmnt in afro.f when integ or gendip are called)
      if(nerihmltn3.eq.1)then
        call hmltn3neri(t,yy,h)
      elseif(nerihmltn3.eq.0)then
        call hmltn3(t,yy,h)
      else
        write(6,*)'error in feval: do not know which hmltn3 to use'
        call myexit
      endif
c--------------------------
      call matify(ajsm,h)
c Change in design trajectory in the general case: no midplane symmetry
      q = 1./brho
      x = yy(1)
      y = yy(3)
c      t = yy(5)  Oh boy!
      Px = yy(2)
      Py = yy(4)
      Pt = yy(6)
c
      k1 = 1
      k2 = 2
      k3 = 3
      k4 = 4
      k5 = 5
      k6 = 6
c
c Hamilton equations , given vector potentials ( Ax, Ay, Az ) :
c
      Pix = yy(2) - Ax(0)
      Piy = yy(4) - Ay(0)
      xm2 = 1.d0/(gamma*beta)**2
c
      root = Dsqrt(-xm2 + yy(6)**2 - Pix**2 - Piy**2)
c
      f(k1) = Pix/root
      f(k2) = (Ax(1)*Pix+Ay(1)*Piy)/root + Az(1)
      f(k3) = Piy/root
      f(k4) = (Ax(3)*Pix+Ay(3)*Piy)/root + Az(3)
      f(k5) = (-yy(6))/root
      f(k6)=0
c      write(6,*) t, Pix, yy(1)
c      write(36,*) t, Pix, yy(1)
      goto 100
c
c  Insert your code here!
c
   40 call hmltn4(t,yy,h)
      call matify(ajsm,h)
      go to 100
  100 continue
c
c  continue to evaluate f's
c
c  matrix part:  dm/dt = j*s*m = -ajsm*m
      do 110 i=7,42
  110 f(i)=0.d0
      do 120 i=1,6
      do 120 j=1,6
      do 120 k=1,6
c compute  xm2dot(i,j)=-ajsm(i,k)*xm2(k,j)
      ij=i+6*j
  120 f(ij)=f(ij) - ajsm(i,k)*yy(k+6*j)
c
c  compute f3dot and f4dot
      call xform(h,3,yy(7),0,hint)
      call xform(h,4,yy(7),1,hint)
      do 130 i=43,98
  130 f(i)=-hint(i-15)
      call pbkt(yy(16),3,hint,3,temp)
      do 140 i=99,224
  140 f(i)=-hint(i-15)-0.5*temp(i-15)
      return
      end
c
****************************************************************************
c
      subroutine putmap(y,fa,fm)
c  Written by Rob Ryne, Spring 1986
      include 'impli.inc'
      include 'param.inc'
      dimension y(224)
      dimension fa(monoms),fm(6,6)
      call clear(fa,fm)
      do 10 i=28,monoms
   10 fa(i)=y(i+15)
      do 20 i=1,6
      do 20 j=1,6
   20 fm(i,j)=y(i+6*j)
      call revf(1,fa,fm)
      return
      end
c
***********************************************************************
c
      subroutine rk78ii(h,ns,t,y)
c  Written by Rob Ryne, Spring 1986, based on a routine of
c  J. Milutinovic.
c  For a reference, see page 76 of F. Ceschino and J Kuntzmann,
c  Numerical Solution of Initial Value Problems, Prentice Hall 1966.
c  This integration routine makes local truncation errors at each
c  step of order h**7.
c  That is, it is locally correct through terms of order h**6.
c  Each step requires 8 function evaluations.
c
      include 'impli.inc'
c
      parameter (ne=224)   
      dimension y(ne),yt(ne),f(ne),a(ne),b(ne),c(ne),d(ne),
     # e(ne),g(ne),o(ne),p(ne)
c
      tint=t
      do 200 i=1,ns
      call feval(t,y,f)
      do 10 j=1,ne
   10 a(j)=h*f(j)
      do 20 j=1,ne
   20 yt(j)=y(j)+a(j)/9.d+0
      tt=t+h/9.d+0
      call feval(tt,yt,f)
      do 30 j=1,ne
   30 b(j)=h*f(j)
      do 40 j=1,ne
   40 yt(j)=y(j) + (a(j) + 3.d+0*b(j))/24.d+0
      tt=t+h/6.d+0
      call feval(tt,yt,f)
      do 50 j=1,ne
   50 c(j)=h*f(j)
      do 60 j=1,ne
   60 yt(j)=y(j)+(a(j)-3.d+0*b(j)+4.d+0*c(j))/6.d+0
      tt=t+h/3.d+0
      call feval(tt,yt,f)
      do 70 j=1,ne
   70 d(j)=h*f(j)
      do 80 j=1,ne
   80 yt(j)=y(j) + (-5.d+0*a(j) + 27.d+0*b(j) -
     # 24.d+0*c(j) + 6.d+0*d(j))/8.d+0
      tt=t+.5d+0*h
      call feval(tt,yt,f)
      do 90 j=1,ne
   90 e(j)=h*f(j)
      do 100 j=1,ne
  100 yt(j)=y(j) + (221.d+0*a(j) - 981.d+0*b(j) +
     # 867.d+0*c(j)- 102.d+0*d(j) + e(j))/9.d+0
      tt = t+2.d+0*h/3.d+0
      call feval(tt,yt,f)
      do 110 j=1,ne
  110 g(j)=h*f(j)
      do 120 j=1,ne
  120 yt(j) = y(j)+(-183.d+0*a(j)+678.d+0*b(j)-472.d+0*c(j)-
     #  66.d+0*d(j)+80.d+0*e(j) + 3.d+0*g(j))/48.d+0
      tt = t + 5.d+0*h/6.d+0
      call feval(tt,yt,f)
      do 130 j=1,ne
  130 o(j)=h*f(j)
      do 140 j=1,ne
  140 yt(j) = y(j)+(716.d+0*a(j)-2079.d+0*b(j)+1002.d+0*c(j)+
     # 834.d+0*d(j)-454.d+0*e(j)-9.d+0*g(j)+72.d+0*o(j))/82.d+0
      tt = t + h
      call feval(tt,yt,f)
      do 150 j=1,ne
  150 p(j)=h*f(j)
      do 160 j=1,ne
  160 y(j) = y(j)+(41.d+0*a(j)+216.d+0*c(j)+27.d+0*d(j)+
     #  272.d+0*e(j)+27.d+0*g(j)+216.d+0*o(j)+41.d+0*p(j))/840.d+0
      t=tint+i*h
  200 continue
      return
      end
c
c end of file
