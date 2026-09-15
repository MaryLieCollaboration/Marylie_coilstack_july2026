************************************************************************
* header                    GENSOL                                     *
*         (GENMAP for a solenoid magnet with soft fringe fields)       *
*  All routines needed for this special GENMAP                         *
************************************************************************
c
      subroutine gensol(p,fa,fm)
c
c This routine computes the map for a solenoid, by numerical integration.
c F. Neri, 8/18/89; A. Dragt, 10/4/89
c The routine is based on Rob Ryne original solnsc, but all the code
c has been rewritten.
c
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
      include 'hmflag.inc'
      include 'combs.inc'
      include 'parm.inc'
      include 'files.inc'
      include 'sol.inc'
c
c  calling arrays
      dimension p(6)
      dimension fa(monoms), fm(6,6)
c
c  local arrays
      dimension y(monoms+15)
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
      ips = nint(p(5))
      mpole = nint(p(6))
c
c  get other parameters from pset
c
      di =  pst(1,ips)
      tl =  pst(2,ips)
      cl =  pst(3,ips)
      bz0 = pst(4,ips)
      iecho = nint(pst(5,ips))
      ioptr = nint(pst(6,ips))
c
c
c  write out parameters if desired
c
      if (iecho .eq. 1 .or. iecho .eq. 3) then
      write(jof,*)
      write(jof,*) ' zi=',zi,' zf=',zf
      write(jof,*) ' ns=',ns
      write(jof,*) ' di=',di,' length=',tl
      write(jof,*) ' cl=',cl,' B=',bz0
      write(jof,*) ' iopt=',ioptr
      write(jof,*)
      endif
      if (iecho .eq. 2 .or. iecho .eq. 3) then
      write(jodf,*)
      write(jodf,*) ' zi=',zi,' zf=',zf
      write(jodf,*) ' ns=',ns
      write(jodf,*) ' di=',di,' length=',tl
      write(jodf,*) ' cl=',cl,' B=',bz0
      write(jodf,*) ' iopt=',ioptr
      write(jodf,*)
      endif
c
      h=(zf-zi)/float(ns)
c
c  write out gradient and derivatives on file ifile:
      ipflag=0
      if (ifile .ne. 0) then
      if (ifile .lt.0 ) then
      ipflag=1
      ifile=-ifile
      endif
      zz = zi
      do 999 ii = 1, ns+1
        call bz02(zz,b0,b2)
        write(ifile,137) zz, b0, b2, 0., 0., 0.
  137 format(6(1x,1pg12.5))
        zz = zz + h
 999  continue
      write(jof,*) ' profile written on file ',ifile
      endif
c
c  return identity map if ifile was < 0
c
      if (ipflag .eq. 1) then
      call ident(fa,fm)
      return
      endif
c
c  initial values for design orbit (in dimensionless units) :
c
      y(1)=0.d0
      y(2)= 0.
      y(3)=0.d0
      y(4)=0.d0
      y(5)=0.d0
      y(6)=-1.d0/beta
c  set constants
      qbyp=1.d0/brho
      ptg=-1.d0/beta
c
c  initialize map to the identity map:
      ne=224
      do 40 i=7,ne
   40 y(i)=0.d0
      do 50 i=1,6
      j=7*i
   50 y(j)=1.d0
c
c  do the computation:
      t=zi
      iflag = 4
      call adam11(h,ns,'start',t,y)
      call putmap(y,fa,fm)
      call csym(1,fm,ans)
c
      return
      end
c
*************************************************************************
c
      subroutine bz02(z,b0,b2)
c  This routine computes b0(z) = Bz(z)
c  and the second derivative b2(z) on axis
c  Alex Dragt 10/4/89
c
      include 'impli.inc'
      include 'sol.inc'
c
      zz = z - di
      call bump0(zz,cl,tl,ans0)
      call bump2(zz,cl,tl,ans2)
      b0 = bz0*ans0
      b2 = bz0*ans2
c
      return
      end
c
**********************************************************************
c
      subroutine bump0(z,cl,tl,ans0)
c
c This routine computes the soft-edge bump function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
c-----------------------------------------------------------
      sgn0(z,cl)=tanh(z/cl)
c-----------------------------------------------------------
      ans0=( sgn0(z,cl) - sgn0(z-tl,cl) )/2.
c
      return
      end
c
************************************************************************
c
      subroutine bump2(z,cl,tl,ans2)
c
c This routine computes the second derivative of the soft-edge bump function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
      zl=z
      zr=z-tl
      call sgn2(zl,cl,ansl)
      call sgn2(zr,cl,ansr)
      ans2=(ansl - ansr)/2.
c
      return
      end
c
*********************************************************************
c
      subroutine sgn2(z,cl,ans)
c
c This subroutine computes the second derivative of the approximating
c signum function
c Alex Dragt 10/4/89
c
      include 'impli.inc'
c
      ans=0.
      if( abs(z/cl) .lt. 30.) then
      ans=-(2./(cl**2))*(tanh(z/cl))/((cosh(z/cl))**2)
      endif
c
      return
      end
c
***********************************************************************
c
      subroutine hmltn4(t,y,h)
c
c  This routine is used to specify the Hamiltonian h for a solenoid.
c  Written by A. Dragt 9/28/89
c
      include 'impli.inc'
      include 'param.inc'
      include 'parm.inc'
      include 'sol.inc'
c
c calling arrays
      dimension h(monoms)
      dimension y(224)
c
c  begin calculation
c
c  compute gradients
      call bz02(t,b0,b2)
c
c scale gradients
c
      b0=sl*b0/brho
      b2=(sl**3)*b2/brho
c
c compute terms in hamiltonian
c
c terms of degree 2
c
      h(7)= b0**2/(8.*sl)
      h(13)= 1./(2.*sl)
      if (ioptr .eq. 0) then
        h(10)= -b0/(2.*sl)
        h(14)= b0/(2.*sl)
      else
        h(10)= 0.d0
        h(14)= 0.d0
      endif
      h(18)= b0**2/(8.*sl)
      h(22)= 1/(2.*sl)
      h(27)= 1/(2.*beta**2*gamma**2*sl)
c
c terms of degree3
c
      h(33)= b0**2/(8.*beta*sl)
      h(53)= 1./(2.*beta*sl)
      h(57)= b0/(2.*beta*sl)
      h(67)= b0**2/(8.*beta*sl)
      h(45)= -b0/(2.*beta*sl)
      h(76)= 1./(2.*beta*sl)
      h(83)= 1./(2.*beta**3*gamma**2*sl)
c
c terms of degree 4
c
      h(84)= b0*(b0**3 - 4.*b2)/(128.*sl)
      h(90)= b0**2/(16.*sl)
      h(140)= 1./(8.*sl)
      h(91)= (b0**3 - b2)/(16.*sl)
      h(141)= b0/(4.*sl)
      h(95)= b0*(b0**3 - 4.*b2)/(64.*sl)
      h(145)= 3.*b0**2/(16.*sl)
      h(155)= (b0**3 - b2)/(16.*sl)
      h(175)= b0*(b0**3 - 4.*b2)/(128.*sl)
      h(87)= -(b0**3 - b2)/(16.*sl)
      h(107)= -b0/(4.*sl)
      h(111)= -(b0**2)/(4.*sl)
      h(121)= -(b0**3 - b2)/(16.*sl)
      h(99)= 3.*b0**2/(16.*sl)
      h(149)= 1./(4.*sl)
      h(159)= b0/(4.*sl)
      h(179)= b0**2/(16.*sl)
      h(130)= -b0/(4.*sl)
      h(195)= 1./(8.*sl)
      h(104)= -(b0**2*(-3. + beta**2))/(16.*beta**2*sl)
      h(154)= -(-3. + beta**2)/(4.*beta**2*sl)
      h(164)= -(b0*(-3. + beta**2))/(4.*beta**2*sl)
      h(184)= -(b0**2*(-3. + beta**2))/(16.*beta**2*sl)
      h(135)= b0*(-3. + beta**2)/(4.*beta**2*sl)
      h(200)= -(-3. + beta**2)/(4.*beta**2*sl)
      h(209)= -(-5. + beta**2)/(8.*beta**4*gamma**2*sl)
c
c add sextupoles
c      h(28)=fsxnr*sl2
c      h(30)=-3.d0*fsxsk*sl2
c      h(39)=-3.d0*fsxnr*sl2
c      h(64)=fsxsk*sl2
c
c  add octupoles
c
c      h(84)=h(84)+focnr*sl3
c      h(86)=-4.0d0*focsk*sl3
c      h(95)=-6.d0*focnr*sl3
c      h(120)=4.d0*focsk*sl3
c      h(175)=h(175)+focnr*sl3
c
      return
      end
c
c end of file
