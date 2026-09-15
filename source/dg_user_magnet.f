c
c****************************************************************************
c Used in computing the gradient from type 18 (user magnets) with Br and B_phi data
c from a data file
c  May 14, 2024 dgm is the "raw" mth gradient- no dimensionless scaling factor for fitting
c     Version with ndrivs, not broken version with ndriv
c
c  Jan 18 2024 changed name of subroutine bincof to bincoeff to avoid conflict with
c  subroutine bincof in file cylmag_user.
c
c  subroutine bincoeff is in this file, is self-initializing, and computes the same binomial
c  coefficients as bincof
c
c  Don't need  xmu0 (permeability of free space) as with current shhet magnets since the user data
c  are in terms of tesla
c  Jan 2024
c  file dg_user.f: used for type18 only
c  Contains dg_user_magnet and its subroutines getcm_charge,
c  helmgrad, etc. Has two self-initializing subroutines: sdcoeffs, bincof
c
c  Single m
c  Single-m user-magnet gradient Sept. 21, 2023
c     Input single-m cubic coefficient arrays, not mutiple-m arrays

c      dg_user_magnet(zn,abrn,bbrn,cbrn,dbrn,
c     &abtn,bbtn,cbtn,dbtn,rref,m,nlines,z,dgm,ndrivs,iwrite)  
      
      subroutine dg_user_magnet(zn,abrn,bbrn,cbrn,dbrn,
     &abtn,bbtn,cbtn,dbtn,rref,m,npts,z,dgm,ndrivs,iwrite)      
      implicit double precision(a-h,o-z)
c-------------------------------------
c  Inputs:
c  zn=array of z points, ends of intervals for cubic fit to Fourier coeffs.
c     to Br_m and Bphi_m.   Note: there are npts-1 z intervals; zn(1) is the start of
c     the first interval, zn(npts) is the end of the last interval. abrn(npts) is set equal
c     to abrn(npts-1), etc. and the coefficients abrn(npts), etc. are not used
c   
c   abrn,bbrn,cbrn,dbrn,abtn,bbtn,cbtn,dbtn=cubic coeffs. to Br_m and Bphi_m. They are
c   extracted from saved 3-index arrays in store_user_coil_data and retrieved by gtype18
c   
c  rref=reference radius for field data
c  m=multipole number for desired gradient
c  npts=no.of z points 
c  z=field point
c  ndrivs=no. of derivatives of gm (return gm itself + ndrivs derivatives of gm)
c     iwrite=integer code for information typeout (iwrite=0, no typeout; iwrite=1 typeout)
c-------------------------------------
c  Outputs:
c  dgm=vector of gradient and its z derivatives (zero indexed)
c------------
c  maxlines=max. no. of coefficient lines the current field data block
      parameter(maxlines=10000)
      parameter(maxdriv=19)
c  Single-m input arrays
c  z values of ends of intervals
      dimension zn(maxlines)
c  Cubic coeffs for br (single m value)
      dimension abrn(maxlines)
      dimension bbrn(maxlines)
      dimension cbrn(maxlines)
      dimension dbrn(maxlines)
c Cubic coeffs for bt (single m value)
      dimension abtn(maxlines)
      dimension bbtn(maxlines)
      dimension cbtn(maxlines)
      dimension dbtn(maxlines)
c-------------
c  Output array:
c  Single-m gradient and its z derivatives
c  dgm=g and its derivatives at z.
c  dgm(0) is the gradient itself
c  dgm(1) is its first z derivative, etc.
c     dimension dgm(0:21)
      dimension dgm(0:maxdriv)
c-----------------
      if(iwrite.eq.1) write(6,*) 'dg_user_magnet m=',m,' ndrivs=',ndrivs
c  Zero out dgm: i=0 for 0th derivative, etc.
      do i=0,maxdriv
      dgm(i)=0.d0
      enddo
c----------------
c
      if(iwrite.eq.1) write(6,*) 'start dg_user_magnet'
c     write(6,*) 'm=',m
      call getcm_charge(m,cmc)
      if(iwrite.eq.1) write(6,*) 'm=',m,' cmc=',cmc
      if(iwrite.eq.1) write(6,*) 'calling helmgrad'
c  helmgrad computes dgm for a single m using cubic coeffs. for npts-1 intervals
      call helmgrad(zn,abrn,bbrn,cbrn,dbrn,abtn,bbtn,cbtn,dbtn,
     &cmc,rref,z,dgm,m,npts,ndrivs,iwrite)
      if(iwrite.eq.1) then
      write(6,*) 'At end of dg_user_magnet return from helmgrad'
      endif
c
      return
      end
c
c   End subroutine dg_user_mag which computes dgm a single  m
c
c*************************************************************************************************
c
      subroutine helmgrad(zn,abrn,bbrn,cbrn,dbrn,abtn,bbtn,cbtn,dbtn,
     &cmc,rref,z,dgm,m,npts,ndrivs,iwrite)
      implicit double precision(a-h,o-z)
c  Computes dgm=vector of on-axis gradient and its z derivatives for a single m using Helmholtz
c  theorem sources represented by cubic coefficients abrm,bbrm,cbrm,dbrm,abtm,bbtm,cbtm,dbtm
c  at many z' values.
c  z' is source point, z (input) is field point
c  The contributions of the npts-1 intervals to the gradient and its derivatives are summed
c  to get the output dgm
c  The cubic coefficients are in terms of the variable t=z'-zn(n). i.e.
c  brm(z')=abrm(n)+bbrm(n)*t+cbrm(n)*t**2+dbrm(n)*t**3 where z' is a source point
c  between zn(n) and zn(n+1)
c
c  Note: rref replaces a in subroutines called by this routine
c
c  Calls getdgmn_t which computes dgmn_t(0)= integral from 0 to x2-x1 dt of
c  f(t) dt  / [a**2+(t+x1-z)**2]**(m+1/2)
c  with f(t)=(alfi(0) + alfi(1)*t + alfi(2)*t**2 +... alfi(na)*t**na)
c     + ndrivs z derivatives of dgmn_t(0)
c  To get the on-axis gradient for a charge sheet alone, the output of getdgmn_t, dgmn_t, is multiplied by cmc*rref**(m+1)
c  but here it is done first with a factor of rref and later by a factor of rref**m after the charge sheet and dipole sheet
c  contributions are combined
c To get the on-axis gradient for a dipole sheet, getdgmn_t is called 2 times, first with m'=m+1, the second with m
c Then the two outputs are combined (see below)
c 
      parameter(maxlines=10000)
      parameter(maxdriv=19)
c  z values of ends of intervals
      dimension zn(maxlines)
c  Cubic coeffs for brm
      dimension abrn(maxlines),bbrn(maxlines),cbrn(maxlines),
     &dbrn(maxlines)
c  Cubic coeffs for btm
      dimension abtn(maxlines),bbtn(maxlines),cbtn(maxlines),
     &dbtn(maxlines)
c  Note: abrm(npts),bbrm(npts) etc. are not used since the number of z values is one larger
c  than the number of coefficient values. abrm(npts),bbrm(npts) etc. are dummied with their
c  values for the npts-1 th interval 
c  dgm=vector of output g and its derivatives at z. dgm(0) is the gradient itself,
c  dgm(1) is its first z derivative, etc.
      dimension dgm(0:maxdriv)
c  Cubic coefficient vector scratch array. 4th component is always zero since we have cubics, not quartics
      dimension alfi(0:4)
c  Scratch arrays for the gradient and its derivatives from a single interval
      dimension dgmn_t(0:maxdriv),dgmn_mt(0:maxdriv),dgmn_mpt(0:maxdriv)
c--------------------------------------------------------
c  Inputs:
c  zn=array of z points that delineate the intervals. The nth interval is [zn(n),zn(n+1)]
c  abrm,bbrm,cbrm,dbrm=cubic coefficients for brm (mth Fourier component of B_r at (rref,z)
c  abtm,bbtm,cbtm,dbtm=cubic coefficients for btm (mth Fourier component of B_theta at (rref,z)
c  cmc=(2m-1)!!/[ 2^(m+1) (m-1)! ]   cmc is an overall factor in dgm
c  rref=radius of field data cylinder (a=rref in the subroutines called by this subroutine)
c  z=field point z
c  m=multipole number (m=1 for dipole, 2 for quadrupole, etc.)
c  npts=no. of z values in array zn. The number of intervals is npts-1
c  ndrivs=no. of z derivatives of gradient to be computed (ndrivs must not be > 21)
c  iwrite=integer flag for diagnostic typeout: 0, no typeout, 1 typeout
c  Output:
c  dgm=vector of gradient and its z derivatives at z
c
c------------------------------------------
      if(iwrite.eq.1) write(6,*) 'start helmgrad'
c  Zero out dgm
      do nd=0,maxdriv
      dgm(nd)=0.d0
      enddo
c-------------------------------------------
 300  format(1x,i4,1x,1pe20.13)
      iwrite=0
c  Step through source z intervals. There are npts-1 intervals
      do n=1,npts-1
      x1=zn(n)
      x2=zn(n+1)
c  Charge-sheet contribution to dgm (from brm)
      alfi(0)=abrn(n)
      alfi(1)=bbrn(n)
      alfi(2)=cbrn(n)
      alfi(3)=dbrn(n)
      alfi(4)=0.d0
      na=3
c
c      if(iwrite.eq.1) write(6,*) 'in helmgrad 1st call getdgmn_t'
      call getdgmn_t(rref,x1,x2,z,alfi,dgmn_t,m,na,ndrivs,iwrite)
      if(iwrite.eq.1) write(6,*) 'in helmgrad after 1st call getdgmn_t'
c----------------------------------------------------
c  Add charge-sheet contribution for nth interval to dgm. Put in factor of rref (i.e. a)
c  for charge sheet here so that overall power of rref factor is correct when it mutiplies
c  the summed the monopole and dipole contributions outside of the do loop
      do nd=0,ndrivs
      dgm(nd)=dgm(nd)+rref*dgmn_t(nd)
      enddo
c-----------------------------------------------------
c  Dipole-sheet contribution of nth interval
c  Note that the btm (B-theta) cubic coefficients are multiplied by rref/m to convert them
c  into dipole-sheet (stream function) coefficients
c  Note also that there is some repeat calculation of the same quantities. This can be fixed (TBD)
c  Should reduce compute time by ~x3
      alfi(0)=rref*abtn(n)/dfloat(m)
      alfi(1)=rref*bbtn(n)/dfloat(m)
      alfi(2)=rref*cbtn(n)/dfloat(m)
      alfi(3)=rref*dbtn(n)/dfloat(m)
      alfi(4)=0.d0
      na=3
      mp=m+1
c      if(iwrite.eq.1) write(6,*) 'in helmgrad 2nd call getdgmn_t'
      call getdgmn_t(rref,x1,x2,z,alfi,dgmn_mpt,mp,na,ndrivs,iwrite)
c      if(iwrite.eq.1) write(6,*) 'in helmgrad 3rd call getdgmn_t'
      call getdgmn_t(rref,x1,x2,z,alfi,dgmn_mt,m,na,ndrivs,iwrite)
c------------------------------------------------
c  Add dipole-sheet contribution of nth interval to dgm. Note sign of terms
      do nd=0,ndrivs
      dgm(nd)=dgm(nd)+dfloat(2*m+1)*rref**2*dgmn_mpt(nd)-
     &dfloat(m)*dgmn_mt(nd)
      enddo
c-----------------------------------------
c  End loop over intervals
      enddo
c--------------------------------------
c  Multiply dgm vector by cmc*rref**m (recall there is already a factor of rref in
c  the charge-sheet terms)
c      call getcm_charge(m,cmc)
      gmac=cmc*rref**m
c      write(6,*) 'cmc in helmgrad=',cmc
      do nd=0,ndrivs
      dgm(nd)=gmac*dgm(nd)
      enddo
c
      if(iwrite.eq.1) write(6,*) 'At end of helmgrad returning'
      return
      end
c
c  End subroutine helmgrad
c
c*************************************************************************************************
c
c
c  Begin subroutine package grad.f now subsumed in file dg_user_magnet.f
c  Gradient and its z derivatives of a polynomial charge distribution
c  (up to quartic) on the interval [x1,x2]
c  The polynomial is in the variable t=z'-x1, where z' is the source point
c  Contains subroutines:
c  getdgmn_t
c  s0tos4int
c  getnhalf
c  get_dnm
c  getrri
c  sdcoeffs (self-initializing)
c  bincoeff formerly bincof) (self-initializing)
c  getcm_charge
c    + num. check routines
c   gccubic_num.f
c   gdcubic_num.f
c******************************************************************
      subroutine getdgmn_t(a,x1,x2,z,alfi,dgmn_t,m,na,ndrivs,iwrite)
      implicit double precision(a-h,o-z)
      double precision ii0,ii1,ii2,ii3,ii4
      double precision ii0_1,ii1_1,ii2_1,ii3_1,ii4_1
      double precision ii0_2,ii1_2,ii2_2,ii3_2,ii4_2
      parameter(maxdriv=19)
      dimension alfi(0:4)
      dimension ddrm1(0:maxdriv)
      dimension ddrm2(0:maxdriv)
      dimension dgmn_t(0:maxdriv)
c  Charge-sheet gradient for a segment (needs to be multiplied by the appropriate
c  constant)
c !! Nov. 30 2021 uses new subroutine get_dnm to compute the derivatives of
c !!  1/(a**2+s**2)**(m+1/2), now up to the 21 st derivative
c  Computes dgmn_t(0)= integral from 0 to x2-x1 dt of
c  f(t) dt  / [a**2+(t+x1-z)**2]**(m+1/2)
c  with f(t)=(alfi(0) + alfi(1)*t + alfi(2)*t**2 +... alfi(na)*t**na)
c  + ndrivs z derivatives of dgmn_t(0)
c      
c  When summed over all of the intervals of the fit, and when the sum is multiplied
c  by the appropriate constant,the result is the on-axis gradient g_m(z) and its
c  z derivatives
c
c  Max. allowable na is 4 (quartic fit)
c  
c  Change variable from t to s = t+x1-z = t-z1 with z1=z-x1
c  Then the numerator f(t) becomes for na=4 (max. no. of terms)
c  A0(z)+A1(z)*s+A2(z)*s**2+A3(z)*s**3+A4(z)*s**4
c  with
c  A0 = alfi(0) + alfi(1)*z1 + alfi(2)*z1**2 + alfi(3)*z1**3 + alfi(4)*z1**4
c  A1 = alfi(1) + 2*alfi(2)*z1 + 3*alfi(3)*z1**2+ 4*alfi(4)*z1**3
c  A2 = alfi(2) + 3*alfi(3)*z1 + 6*alfi(4)*z1**2
c  A3 = alfi(3) + 4*alfi(4)*z1 
c  A4 = alfi(4)
c
c  Indexing of bcoeff (binomial coefficients) does NOT use zero indexing:
c  (n)                (n)
c  | | = bcoeff(1,n)  | | = bcoeff(2,n)  etc. 
c  (0)                (1)
c  Calls subroutine s0123int to get integrals of s**n ds / (a**2+s**2)**(m+1/2)
c--------------
c  ndrivs should be odd and must be no greater than 21
c  If ndrivs is even, this subroutine will compute derivatives up to
c  and including ndrivs+1
c----------------------
c  No. of derivatives computed is sometimes more than ndrivs:
c  na=0,1,2,3: will always compute up to and including 3rd derivative
c  na=4: will always compute up to and including 4th derivative
c----------------------
c  Nov. 30 2021 uses new subroutine get_dnm to compute the derivatives of
c  1/(a**2+s**2)**(m+1/2), now up to the 21 st derivative
c  P. L. Walstrom  Nov. 20, 2021 checked by numerical integration/differentiation
c  Nov. 17 alfi now has 5 numbers- source function can be quartic
c  Used in computing the on-axis generalized gradient of a charge sheet
c
c  Inputs:
c  a=source cylinder radius
c  x1,x2= x limits of integration
c  z=z of field point
c  alfi(i),i=0,..na vector of coefficients in piecewise polynomial f_m(t) between 0 and x2-x1
c  m=multipole index (=1 for dipole, =2 for quadrupole, etc.)
c  na=degree of polynomial in t=x-x1 that is the fit to the source function on [x1,x2]
c  na=0: polynomial is a constant; na=1: polynomial is linear in x, etc. (max. na=4)
c  ndrivs=no. of z derivatives computed (but not always for smaller ndrivs see comments above)
c
c
c  Outputs:
c  dgmn_t(0)=integral of sum_(i=0)^na alfi(i)t**(i-1) dt /[a**2+(t+x1-z)**2]**(m+1/2)
c          from 0 to x2-x1
c  dgmn_t(1)=d/dz dgmn_t(0)
c  dgmn_t(2)=d**2/dz**2 dgmn_t(0)
c     etc.
c
      if(na.gt.4) then
      write(6,*) 'na > 4 in getdgmn_t- stopped'
      write(6,*) 'na=',na
      stop
      endif
c----------------------
      if(ndrivs.gt.maxdriv) then
      write(6,*) 'ndrivs > maxdriv=19) in getdgmn_t- stopped'
      write(6,*) 'ndrivs=',ndrivs
      stop
      endif
c----------------------
c  Zero out dgmn_t
      do n=0,maxdriv
      dgmn_t(n)=0.d0
      enddo
c----------------------
      s1=x1-z
      s2=x2-z
      z1=z-x1
      dd=x2-x1
c
      if(iwrite.eq.1) write(6,*) 'a,s1,s2,z1,dd=',a,s1,s2,z1,dd
      if(iwrite.eq.1) write(6,*) 'm=',m
c Compute nhalf = (ndrivs-1)/2 if ndrivs is odd, ndrivs/2 if ndrivs is even
c     write(6,*) 'Call getnhalf'
      call getnhalf(ndrivs,nhalf)
      if(iwrite.eq.1) write(6,*) 'ndrivs=',ndrivs,'nhalf=',nhalf
      nhaff=nhalf
      if(nhaff.lt.2) nhaff=2
c  Get repeated derivatives of 1/(a**2+s**2)**(m+1/2), evaluated at
c  both s1 and s2
c       call getsderivs(a,s1,m,rmn1,ddrm1,nsdrivs)
c       call getsderivs(a,s2,m,rmn2,ddrm2,nsdrivs)
      call get_dnm(a,s1,m,ddrm1,nhaff)
      call get_dnm(a,s2,m,ddrm2,nhaff) 
c
c  Compute s integrals int s**n ds/(a**2+s**2)**(m+1/2), n=0,...na
c  s0123int computes only ii0, returns ii1=ii2=ii3=0 if na=0,
c     "     computes only ii0 and ii1, returns ii2=ii3=0 if na=1, etc.
c      call s0123int(s1,a,m,ii0_1,ii1_1,ii2_1,ii3_1,na)
c     call s0123int(s2,a,m,ii0_2,ii1_2,ii2_2,ii3_2,na)
      call s0tos4int(s1,a,m,ii0_1,ii1_1,ii2_1,ii3_1,ii4_1,na)
      call s0tos4int(s2,a,m,ii0_2,ii1_2,ii2_2,ii3_2,ii4_2,na)
      ii0=ii0_2-ii0_1
      ii1=ii1_2-ii1_1
      ii2=ii2_2-ii2_1
      ii3=ii3_2-ii3_1
      ii4=ii4_2-ii4_1
c
c  Compute dgmn_t(n),n=0,1,...ndrivs
c  
      select case(na)
c
c  na=0
      case(0)
c  f(t) is constant alfi(0)
c  gmint=alfi(0)*ii0
      fx1=alfi(0)
      fx2=alfi(0)
      dgmn_t(0)=alfi(0)*ii0
c     
c  The remaining derivatives do not contain ii0 or higher s**n integrals
      dgmn_t(1)=alfi(0)*(ddrm1(0)-ddrm2(0))
      dgmn_t(2)=-alfi(0)*(ddrm1(1)-ddrm2(1))
      dgmn_t(3)=alfi(0)*(ddrm1(2)-ddrm2(2))
      if(ndrivs.lt.4) return
      sign=1.d0
      do i=4,ndrivs
      sign=-sign   
      dgmn_t(i)=sign*(fx1*ddrm1(i-1)-fx2*ddrm2(i-1))
      enddo
c      write(6,*) 'd4gm=',d4gm,' d5gm=',d5gm
c      write(6,*) 'digm(4)=',digm(4),'digm(5)=',digm(5)
      return
c
      case(1)
c  na=1
c  f(t) is linear; alfi(2)=0,alfi(3)=0,alfi(4)=0
c     write(6,*) 'na=',na
c  Nov. 11 2021 new: f(t) version
c  dd defined above: dd = x2-x1 = upper limit of t=Delta in memo
      fx1=alfi(0)
      fx2=alfi(0)+alfi(1)*dd
      fpx1=alfi(1)
      fpx2=alfi(1)
c  Expand powers of t=s-z1 in powers of s, define coefficients for this
      aa0=alfi(0)+alfi(1)*z1
      aa1=alfi(1)
c
c     dgmn_t(0)=(alfi(0)+alfi(1)*z1)*ii0+alfi(1)*ii1
      dgmn_t(0)=aa0*ii0+aa1*ii1
      dgmn_t(1)=fx1*ddrm1(0)-fx2*ddrm2(0)+alfi(1)*ii0
c     
c  The remaining derivatives do not contain ii0 or higher s**n integrals
      dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)
      dgmn_t(3)=fx1*ddrm1(2)-fx2*ddrm2(2)-fpx1*ddrm1(1)+fpx2*ddrm2(1)
      if(ndrivs.lt.4) return
      sign=1.d0
      do i=4,ndrivs
      sign=-sign   
      dgmn_t(i)=sign*(fx1*ddrm1(i-1)-fx2*ddrm2(i-1)-fpx1*ddrm1(i-2)+
     &fpx2*ddrm2(i-2))
      enddo
      return
c
c  na=2
      case(2)
c  f(t) is a quadratic
c     write(6,*) 'na=',na
c  Nov. 17 use coefficients aa0,aa1,aa2
c  dd defined above: dd = x2-x1 = upper limit of t
      fx1=alfi(0)
      fx2=alfi(0)+alfi(1)*dd+alfi(2)*dd**2
      fpx1=alfi(1)
      fpx2=alfi(1)+2.d0*alfi(2)*dd
      fpp=2.d0*alfi(2)
c  Expand powers of t=s-z1 in powers of s, define coefficients for this
      aa0=alfi(0)+alfi(1)*z1+alfi(2)*z1**2
      aa1=alfi(1)+2.d0*alfi(2)*z1
      aa2=alfi(2)
      
c      dgmn_t(0)=(alfi(0)+alfi(1)*z1+alfi(2)*z1**2)*ii0+
c     &(alfi(1)+2.d0*alfi(2)*z1)*ii1+alfi(2)*ii2
       dgmn_t(0)=aa0*ii0+aa1*ii1+aa2*ii2
c
      dgmn_t(1)=fx1*ddrm1(0)-fx2*ddrm2(0)+
     &(alfi(1)+2.d0*alfi(2)*z1)*ii0+2.d0*alfi(2)*ii1
      dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)+
     &2.d0*alfi(2)*ii0
c
c The remaining derivatives do not contain ii0 or higher s**n integrals
      dgmn_t(3)=fx1*ddrm1(2)-fx2*ddrm2(2)-fpx1*ddrm1(1)+fpx2*ddrm2(1)+
     &fpp*(ddrm1(0)-ddrm2(0))
      if(ndrivs.lt.4) return
c
      sign=1.d0
      do i=4,ndrivs
      sign=-sign   
      dgmn_t(i)=sign*(fx1*ddrm1(i-1)-fx2*ddrm2(i-1)-fpx1*ddrm1(i-2)+
     &fpx2*ddrm2(i-2)+fpp*(ddrm1(i-3)-ddrm2(i-3)))
      enddo
      return
c  na=3
      case(3)
c      write(6,*) 'na=',na
c f(t) is a cubic
      fx1=alfi(0)
      fx2=alfi(0)+alfi(1)*dd+alfi(2)*dd**2+alfi(3)*dd**3
      fpx1=alfi(1)
      fpx2=alfi(1)+2.d0*alfi(2)*dd+3.d0*alfi(3)*dd**2
      fppx1=2.d0*alfi(2)
      fppx2=2.d0*alfi(2)+6.d0*alfi(3)*dd
      fppp=6.d0*alfi(3)
c  Expand powers of t=s-z1 in powers of s, define coefficients for this
      aa0=alfi(0)+alfi(1)*z1+alfi(2)*z1**2+alfi(3)*z1**3
      aa1=alfi(1)+2.d0*alfi(2)*z1+3.d0*alfi(3)*z1**2
      aa2=alfi(2)+3.d0*alfi(3)*z1
      aa3=alfi(3)
c
      dgmn_t(0)=aa0*ii0+aa1*ii1+aa2*ii2+aa3*ii3
      if(iwrite.eq.1) write(6,*) 'dgmn_t(0)=',dgmn_t(0)
c
      dgmn_t(1)=fx1*ddrm1(0)-fx2*ddrm2(0)+
     &(alfi(1)+2.d0*alfi(2)*z1+3.d0*alfi(3)*z1**2)*ii0+
     &(2.d0*alfi(2)+6.d0*alfi(3)*z1)*ii1+3.d0*alfi(3)*ii2
c
      daa0=alfi(1)+2.d0*alfi(2)*z1+3.d0*alfi(3)*z1**2
      daa1=2.d0*alfi(2)+6.d0*alfi(3)*z1
      daa2=3.d0*alfi(3)
c      daa3=0.d0
c
      dgmn_t(1)=fx1*ddrm1(0)-fx2*ddrm2(0)+daa0*ii0+daa1*ii1+daa2*ii2
c
      ddaa0=2.d0*alfi(2)+6.d0*alfi(3)*z1
      ddaa1=6.d0*alfi(3)
c
c      dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)+
c     &(2.d0*alfi(2)+6.d0*alfi(3)*z1)*ii0+6.d0*alfi(3)*ii1
      dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)+
     &ddaa0*ii0+ddaa1*ii1
c
      dgmn_t(3)=fx1*ddrm1(2)-fx2*ddrm2(2)-fpx1*ddrm1(1)+fpx2*ddrm2(1)+
     &fppx1*ddrm1(0)-fppx2*ddrm2(0)+6.d0*alfi(3)*ii0
c
c  The remaining derivatives do not contain ii0 or higher s**n integrals
      sign=1.d0
      do i=4,ndrivs
      sign=-sign   
      dgmn_t(i)=sign*(fx1*ddrm1(i-1)-fx2*ddrm2(i-1)-fpx1*ddrm1(i-2)+
     &fpx2*ddrm2(i-2)+fppx1*ddrm1(i-3)-fppx2*ddrm2(i-3)-
     &6.d0*alfi(3)*(ddrm1(i-4)-ddrm2(i-4)))
      enddo
      return
c 
      case(4)
c  na=4
c     f(t) is a quartic
      fx1=alfi(0)
      fx2=alfi(0)+alfi(1)*dd+alfi(2)*dd**2+alfi(3)*dd**3+alfi(4)*dd**4
      fpx1=alfi(1)
      fpx2=alfi(1)+2.d0*alfi(2)*dd+3.d0*alfi(3)*dd**2+4.d0*alfi(4)*dd**3
      fppx1=2.d0*alfi(2)
      fppx2=2.d0*alfi(2)+6.d0*alfi(3)*dd+12.d0*alfi(4)*dd**2
      fpppx1=6.d0*alfi(3)
      fpppx2=6.d0*alfi(3)+24.d0*alfi(4)*dd
      fpppp=24.d0*alfi(4)
c  Expand powers of t=s-z1 in powers of s, define coefficients for this
c  aa0 is coeff. of s**0, aa1 is coeff. of s**1, etc.
      aa0=alfi(0)+alfi(1)*z1+alfi(2)*z1**2+alfi(3)*z1**3+alfi(4)*z1**4
      aa1=alfi(1)+2.d0*alfi(2)*z1+3.d0*alfi(3)*z1**2+4.d0*alfi(4)*z1**3
      aa2=alfi(2)+3.d0*alfi(3)*z1+6.d0*alfi(4)*z1**2
      aa3=alfi(3)+4.d0*alfi(4)*z1
      aa4=alfi(4)
c
      dgmn_t(0)=aa0*ii0+aa1*ii1+aa2*ii2+aa3*ii3+aa4*ii4
c  daa0 is coeff. of s**0, daa1 is coeff. of s**1, etc. in s the integral part of
c  dgmn_t(01)
      daa0=alfi(1)+2.d0*alfi(2)*z1+3.d0*alfi(3)*z1**2+4.d0*alfi(4)*z1**3
      daa1=2.d0*alfi(2)+6.d0*alfi(3)*z1+12.d0*alfi(4)*z1**2
      daa2=3.d0*alfi(3)+12.d0*alfi(4)*z1
      daa3=4.d0*alfi(4)
c Note: daa4=0
      dgmn_t(1)=fx1*ddrm1(0)-fx2*ddrm2(0)+daa0*ii0+daa1*ii1+daa2*ii2+
     &daa3*ii3
c
      ddaa0=2.d0*alfi(2)+6.d0*alfi(3)*z1+12.d0*alfi(4)*z1**2
      ddaa1=6.d0*alfi(3)+24.d0*alfi(4)*z1
      ddaa2=12.d0*alfi(4)
c     ddaa0=2.d0*alfi(2)+6.d0*alfi(3)*z1
c     ddaa1=6.d0*alfi(3)
c     dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)+
c    &ddaa0*ii0+ddaa1*ii1
      dgmn_t(2)=-fx1*ddrm1(1)+fx2*ddrm2(1)+fpx1*ddrm1(0)-fpx2*ddrm2(0)+
     &ddaa0*ii0+ddaa1*ii1+ddaa2*ii2
c
      dddaa0=6.d0*alfi(3)+24.d0*alfi(4)*z1
      dddaa1=24.d0*alfi(4)
      dgmn_t(3)=fx1*ddrm1(2)-fx2*ddrm2(2)-fpx1*ddrm1(1)+fpx2*ddrm2(1)+
     &fppx1*ddrm1(0)-fppx2*ddrm2(0)+dddaa0*ii0+dddaa1*ii1
c 
      d4aa0=24.d0*alfi(4)
      dgmn_t(4)=-fx1*ddrm1(3)+fx2*ddrm2(3)+fpx1*ddrm1(2)-fpx2*ddrm2(2)-
     &fppx1*ddrm1(1)+fppx2*ddrm2(1)+fpppx1*ddrm1(0)-fpppx2*ddrm2(0)+
     &d4aa0*ii0
c  The remaining derivatives do not contain ii0 or higher s**n integrals
      sign=-1.d0
      do i=5,ndrivs
      sign=-sign   
      dgmn_t(i)=sign*(fx1*ddrm1(i-1)-fx2*ddrm2(i-1)-fpx1*ddrm1(i-2)+
     &fpx2*ddrm2(i-2)+fppx1*ddrm1(i-3)-fppx2*ddrm2(i-3)-
     &fpppx1*ddrm1(i-4)+fpppx2*ddrm2(i-4)+d4aa0*(ddrm1(i-5)-ddrm2(i-5)))
      enddo
      case default
      write(6,*) 'na not 0,1,2,3, or 4 in getgmint- stopped'
      stop
      end select
      end
c  End subroutine getdgmn     
c      
c**************************************************
       subroutine s0tos4int(s,a,m,ii0,ii1,ii2,ii3,ii4,na)
      implicit double precision(a-h,o-z)
      double precision ii0,ii1,ii2,ii3,ii4
      double precision numer
c  P. L. Walstrom Checked by numerical integration 11/16/2021
c  Revised version of s0123int, with s**4 integral added
c  Set maxcof=35 in order to agree with self-initializing bincoeff (formerly bincof)
      parameter(maxcof=35)
      dimension bcoff(maxcof,maxcof)
c     Obsolete common- use self-latching bincoeff
c  common /bicof/ bcoeff(maxcof,maxcof)
c  P. L. Walstrom Checked by numerical integration 11/16/2021
c  ii0,ii1,ii2,ii3,ii4 are INDEFINITE integrals
c  Allowable m values are 1,2,...35. This is overkill
c  Allowable na values are 0,1,2,3,4
c  ii0= integral of  ds  / [a**2+s**2]**(m+1/2)
c  ii1= integral of s ds  / [a**2+s**2]**(m+1/2)
c  ii2= integral of s**2 ds  / [a**2+s**2]**(m+1/2)
c  ii3= integral of s**3 ds  / [a**2+s**2]**(m+1/2)
c  ii4= integral of s**4 ds  / [a**2+s**2]**(m+1/2)
c  Note: there is some repeated calculation of some quantities, and this routine
c  could be made a little more efficient.
c  Used in computing the on-axis generalized gradient of a charge sheet
c  Max. allowable na is 4 
c  Use formulas from Gradstein & Ryzhik, p. 86 and p. 87
c
c  Indexing of bcoeff (binomial coefficients) does NOT use zero indexing:
c  (n)                (n)
c  | | = bcoeff(1,n)  | | = bcoeff(2,n)  etc. 
c  (0)                (1)
c
c  Call self-initializing version of bincof
c changed 1-18-24    call bincof(bcoff)
      call bincoeff(bcoff)
c  
      if(na.gt.4) then
      write(6,*) 'na > 4 in s0123int- stopped'
      write(6,*) 'na=',na
      stop
      endif
c
      ii1=0.d0
      ii2=0.d0
      ii3=0.d0
      ii4=0.d0
c
      a2=a**2
      usq=a2+s**2
      u=sqrt(usq)
     
c-----------------------
c  Compute I0
      if(m.eq.1) then
c  s**0 term m=1 requires special treatment since bcoeff(k,0) is undefined
      ii0=s/(u*a2)
c
      else
c
c  m=2 or greater
c  s**0 term Use G&R 2.271 6.
      sign=-1.d0
      ii0=0.d0
      denom=u
      numer=s
      do k=0,m-1
      sign=-sign
      ii0=ii0+sign*bcoff(k+1,m-1)*numer/(denom*dfloat(2*k+1))
      numer=numer*s**2
      denom=denom*usq
      enddo
c ii0=s**0 integral
      ii0=ii0/a2**m
      endif
c      write(6,*) 'ii0=',ii0
      if(na.lt.1) return
c------------
c Compute ii1=s**1 term
      if(m.eq.1) then
      ii1=-1.d0/u
      else
      ii1=-1.d0/(usq**(m-1)*u*dfloat(2*m-1))
      endif
c      write(6,*) 'in line ii1=',ii1
      if(na.lt.2) return
c----------------
c  Compute ii2=s**2 term
c  s**2 term 
      if(m.lt.2) then
c  m=1: can't use G&R 2.272 6.
c  Instead use G&R 2.272 4.
      ii2=-s/u+asinh(s/a)
c      write(6,*) 'm=1 ii2=',ii2
c
      else
c
      if(m.lt.3) then
      ii2=s**3/(usq*u*3.d0*a2)
c      write(6,*) 'm=2 ii2=',ii2
      else
c  m > 2 Use G&R 2.272 6.
      sign=-1.d0
      ii2=0.d0
      denom=usq*u
      numer=s**3
      do k=0,m-2
      sign=-sign
      ii2=ii2+sign*bcoff(k+1,m-2)*numer/(denom*dfloat(2*k+3))
      numer=numer*s**2
      denom=denom*usq
      enddo
c ii2=s**2 integral, m>2
      ii2=ii2/a2**(m-1)
c      write(6,*) 'm > 2 ii2=',ii2
      endif
      endif
      if(na.lt.3) return
c------------
c   Compute ii3=s**3 term G&R 2.273 7.
c
      if(m.eq.1) then
      ii3=(s**2+2.d0*a2)/u
c      write(6,*) 'm=1 ii3=',ii3
      go to 1
      endif
c
      if(m.eq.2) then
      ii3=-1.d0/u+a2/(usq*u*3.d0)
c      write(6,*) 'm=2 ii3=',ii3
      go to 1
      endif
c
c      write(6,*) 'na=3,m>2'
      ii3=-1.d0/(usq**(m-2)*u*dfloat(2*m-3))+
     &a2/(usq**(m-1)*u*dfloat(2*m-1))
c      write(6,*) 'm=3 ii3=',ii3
 1    continue
      if(na.lt.4) return
c
c  Compute s**4 integral ii4
      if(m.eq.1) then
      ii4=0.5d0*s*u+a2*s/u-1.5d0*a2*asinh(s/a)
      return
      endif
c
      if(m.eq.2) then
      ii4=-s/u-s**3/(3.d0*usq*u)+asinh(s/a)
      return
      endif
c
      if(m.eq.3) then
      ii4=0.2d0*s**5/(a2*usq**2*u)
      return
      endif
c
c s**4 integral with m=4 and greater
      sign=-1.d0
      ii4=0.d0
      denom=usq**2*u
      numer=s**5
      do k=0,m-3
      sign=-sign
      ii4=ii4+sign*bcoff(k+1,m-3)*numer/(denom*dfloat(2*k+5))
      numer=numer*s**2
      denom=denom*usq
      enddo
      ii4=ii4/a2**(m-2)
c   
      return
      end
c  End subroutine s0tos4int  
c
c**************************************************
c
      subroutine getnhalf(ndrivs,nhalf)
      integer ievod,ndrivs,nhalf,i
      ievod=ndrivs/2
      i=2*ievod
      if(ndrivs.gt.i) then
      nhalf=(ndrivs-1)/2
      else
      nhalf=ievod
      endif
      return
      end
c
c*************************************************      
c
      subroutine get_dnm(a,s,m,dnm,nhalf)
      implicit double precision(a-h,o-z)
      parameter(maxdriv=19)
      dimension rri(0:maxdriv)
      dimension aij(maxdriv,maxdriv)
      dimension dnm(0:maxdriv)
c  Nov. 30 2021 
c  Computes s derivatives of 1/(a**2+s**2)**(m+1/2) up to 21st derivative
      ndrivs=2*nhalf+1
      dd=a**2+s**2
      call getrri(rri,dd,m,ndrivs)
      call sdcoeffs(aij)
c  Each value of k corresponds to one even derivative and one odd derivative
c  The first value in dnm is the zeroth derivative
      dnm(0)=rri(0)
      dnm(1)=-2.d0*s*rri(1)
      sign=1.d0
      do k=1,nhalf
      sign=-sign
c
c  jth derivative is even: j=2*k, imin=k,imax=2*k, no. of terms=k+1
      j=2*k
      imin=k
      imax=2*k
      dnm(j)=0.d0
      sygn=-1.d0
      nterms=k+1
      imin=k
      do n=1,nterms
      i=imin+n-1 
      sygn=-sygn
      ns=2*n-2
      sgn=sign*sygn
c      write(6,*) 'k,i,j,aij,ns,sgn=',k,i,j,aij(n,j),ns,sgn
      dnm(j)=dnm(j)+sygn*aij(n,j)*rri(i)*s**ns   
      enddo
      dnm(j)=sign*dnm(j)
c
c  jth derivative is odd: j=2*k+1, imin=k+1,imax=2*k+1, no. of terms=k+1
      j=2*k+1
      imin=k+1
      dnm(j)=0.d0
      sygn=-1.d0
      do n=1,nterms
      i=imin+n-1
      sygn=-sygn
      sgn=sign*sygn
c      write(6,*) 'k,i,j,aij,ns,sgn=',k,i,j,aij(n,j),ns,sgn
      ns=2*n-1
      dnm(j)=dnm(j)+sygn*aij(n,j)*rri(i)*s**ns     
      enddo
      dnm(j)=-sign*dnm(j)
c      write(6,*) ' '
      enddo
      return
      end
c  End subroutine get_dnm
c
c*************************************************************************
      subroutine getrri(rri,dd,m,imax)
      implicit double precision(a-h,o-z)
      dimension rri(0:imax)
c  Computes vector array rri(i) from i=0 to i=imax
c  rri(i)=p*(p+1)...(p+i-1)/dd**(p+i)
c  where
c  dd=(a**2+s**2)
c  p=m+1/2
c
      p=dfloat(m)+0.5d0
      rri(0)=1.d0/(dd**m*sqrt(dd))
      do i=1,imax
      rri(i)=(p+dfloat(i-1))*rri(i-1)/dd
      enddo
      return
      end
c
c************************************************************************
c
      subroutine sdcoeffs(aij)
      implicit double precision(a-h,o-z)
      parameter(maxdriv=19)
      dimension aij(maxdriv,maxdriv),aaij(maxdriv,maxdriv)
      data iter/0/
      save aaij
      save iter
c  This is a self-initializing routine with built-in nhalf=10
c  (always computes max. aij array (21x21) in the first call; also returns
c  max. arrray in subsequent calls)
c  Computes aij on first call, after first call returns precomputed
c  and stored aij
c  Computes the array aij by iteration
c  The aij are used in computing numerical coefficients for the
c  repeated s derivatives of 1/(a**2+s**2)**p
c  The aij are independent of p and are all positive integers
c  The alternating signs in the terms in the derivatives are entered
c  later in the subroutine get_dnm
c  This calculation needs to be done only once
c     P. L. Walstrom Nov. 23, 2021 checked against Maxima
c
c  Skip calculation of aij if it has already been done
      if(iter.gt.0) go to 1
c      write(6,*) 'sdcoeffs with iter=0'
c      if(nhalf.gt.10) then
c      write(6,*) 'nhalf > 10'
c      write(6,*) 'nhalf in sdcoeffs =',nhalf
c      write(6,*) 'is too large for array aij'
c      write(6,*) 'Stopped in subroutine sdcoeffs'
c      stop
c      endif
c
      nhalf=10
c
      do j=1,maxdriv
      do i=1,maxdriv
      aaij(i,j)=0.d0
      enddo
      enddo
c  aaij indexing: aaij(i,j): j is order of derivative, i is index of term
c  in jth derivative.
c  1st derivative
      aaij(1,1)=2.0
c  2nd derivative
      aaij(1,2)=2.d0
      aaij(2,2)=4.d0
c  3rd derivative
      aaij(1,3)=12.d0
      aaij(2,3)=8.d0
c  4th derivative
      aaij(1,4)=12.d0
      aaij(2,4)=48.d0
      aaij(3,4)=16.d0
c  5th derivative
      aaij(1,5)=120.d0
      aaij(2,5)=160.d0
      aaij(3,5)=32.d0
c  6th through 2*nhalf+1 derivatives
      do k=3,nhalf
c  2*k th derivative
      j=2*k
      aaij(1,j)=aaij(1,j-1)
      do i=2,k
c      t1=2.d0*aaij(i-1,j)
c      t2=dfloat(2*i-1)*aaij(i,j-1)
c      sum=t1+t2
      aaij(i,j)=2.d0*aaij(i-1,j-1)+dfloat(2*i-1)*aaij(i,j-1)
c      write(6,*)  'i,j,aaij(i,j)=',i,j,aaij(i,j)
      enddo
      aaij(k+1,j)=2.d0*aaij(k,j-1)
c  2*k+1 th derivative
      j=2*k+1
      sign=-1.d0
      do i=1,k
      aaij(i,j)=2.d0*aaij(i,j-1)+dfloat(2*i)*aaij(i+1,j-1)
      enddo
      aaij(k+1,j)=2.d0*aaij(k+1,j-1)
      enddo
c  Transfer saved aaij to aij for return
      do j=1,maxdriv
      do i=1,maxdriv
      aij(i,j)=aaij(i,j)
      enddo
      enddo
c
      iter=1
      return
 1    continue
c      write(6,*) 'sdcoeffs with iter=1'
c  Transfer saved aaij to aij
      do j=1,maxdriv
      do i=1,maxdriv
      aij(i,j)=aaij(i,j)
      enddo
      enddo
      return
      end
c
c     End subroutine sdcoeffs
c********************************************************
c
c     old pre Jan18, 2024      subroutine bincof(bcoff)
      subroutine bincoeff(bcoff)
c  Calculates the same numbers as bincof but eliminates numbered do loops
c  calculates and saves binomial coefficients up to n=35
c  Self initializing: skips calculation if already done
c
c  Does not have zero indexing, so
c  ( n  )
c  ( 0  ) = bcoeff(1,n)
c
c  ( n  )
c  ( 1  ) = bcoeff(2,n)
c
c  ( n )
c  ( k )  = bcoeff(k+1,n)
c
c  ( n  )
c  ( n  ) = bcoeff(n+1,n)
c
c  etc.
c      
      implicit double precision(a-h,o-z)
      parameter(maxcof=35)
      dimension bcoeff(maxcof,maxcof) 
      dimension bcoff(maxcof,maxcof)
      dimension b(maxcof),bold(maxcof)
      data (bcoeff(k,1),k=1,2) /1.d0,1.d0/
      data (bcoeff(k,2),k=1,3) /1.d0,2.d0,1.d0/
      data iter/0/
      save bcoeff
      save iter
      if(iter.eq.0) then
      bold(1)=1.d0
      bold(2)=2.d0
      bold(3)=1.d0
      maxn=maxcof-1
      do n=3,maxn
      b(1)=1.d0
      np1=n+1
      b(np1)=1.d0
      do  k=2,n
      km=k-1
      b(k)=bold(km)+bold(k)
      enddo
      do k=1,np1
      bcoeff(k,n)=b(k)
      bold(k)=b(k)
      enddo
      enddo
      iter=1
      endif
c
      do n=1,maxcof
      do m=1,maxcof
      bcoff(n,m)=bcoeff(n,m)
      enddo
      enddo
      return
      end
c
c  End subroutine bincof(bcoff)
c**********************************************************************
c
      subroutine getcm_charge(m,cmc)
      implicit double precision(a-h,o-z)
c  Analogous to subroutine coeffs in cylmag_usr_singm.f but directly called
c  for a single m, not an array w/ various m values
c  Leading constant for gradient of charge sheet
c  Computes cmc=(2m-1)!!/[ 2^(m+1) (m-1)! ]
c     m must be 1 or greater
      if(m.lt.1) then
      write(6,*) 'm<1 in getcm_charge'
      write(6,*) 'm=',m
      write(6,*) 'stopped'
      stop
      endif
c  m=1 is a special case:
      cmc=0.25d0
      if(m.lt.2) return
      do i=2,m
      xi=dfloat(i)
      cmc=cmc*0.5d0*(2.d0*xi-1.d0)/(xi-1.d0)
      enddo
      return
      end
c     End subroutine getcm_charge
c*****************************************************************************************
c
      
