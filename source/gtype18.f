c
c*****************************************************************************
c
      subroutine gtype18(init,z,a,xl,nr,glprod,mval,ndriv,dg)
      implicit double precision(a-h,o-z)
c  Feb. 20, 2025 multiply dg_user_magnet gradient by scaling factor -glprod/xl
c  Input variables:
c  init=initialization code
c       init=0 initialize (read in user data)
c       init=1 compute gradient
c  z=z coordinate in local coil frame (z relative to coil's nominal center)
c  a=coil winding radius (not used except to check validity of reference radius of coil data)
c  xl=nominal coil length
c  glprod=gradient-length product. glprod=xl if the gradients are unscaled
c  mval=m value of coil (1 for dipole, 2 for quad, etc.) negative mval means gradient of skew field.
c  ndriv=number of times gradient is differentiated
c     iwrite=code for typeout of comments. iwrite=0, no typeout; iwrite >0 , type out comments
c  Output:
c  dg=array of on-axis generalized gradient and its z derivatives. dg(1)=gradient, dg(2)=dg/dz, etc.
c
c     parameter(maxdrv=20,maxdriv=19)
c     set maxdrv=14 to agree with gradients_plw
      parameter(maxdrv=14)
      parameter(maxdriv=19)
      parameter(maxlines=10000)
c  dg index starts with 1. dg(1)=gradient, dg(2)=dg/dz etc.
c  ?TBD i_nr varies from 1 to no_nr
      dimension dg(maxdrv)
      dimension dgm(0:maxdriv)
c  Compatible with cylmag_usr_singm.f in /Users/peterwalstrom/work/Marylie_user_singlem                                                             
      dimension zn(maxlines)
c
      dimension abrn(maxlines)
      dimension bbrn(maxlines)
      dimension cbrn(maxlines)
      dimension dbrn(maxlines)
c
      dimension abtn(maxlines)
      dimension bbtn(maxlines)
      dimension cbtn(maxlines)
      dimension dbtn(maxlines)
c
c  Saved arrays for multiple user coils
      dimension rrefi(10),nlinesi(10)
c
      dimension zni(maxlines,10)
c
      dimension abrni(maxlines,10)
      dimension bbrni(maxlines,10)
      dimension cbrni(maxlines,10)
      dimension dbrni(maxlines,10)
c
      dimension abtni(maxlines,10)
      dimension bbtni(maxlines,10)
      dimension cbtni(maxlines,10)
      dimension dbtni(maxlines,10)
c
c  The integer array nrvals contains the unique nr values that have been entered
c  The subroutine returns with no data file read if an input nr value has
c  already been used and appears in the array nrvals. Max. no. of stored nr values
c  is 10 (for now) Note larger maxnr leads to very large stored user data arrays.
c  Note that if stored user coil data is indexed with icoil instead of the nr index
c  there is possible redundancy and also very large saved arrays since max. icoil is 100
      dimension nrvals(10)
c  Initial value of no_nr is 0
c  initial nrvals is all zeroes
c      data no_nr/0/
c      data nrvals/10*0/
c  ini=0 for first call to this routine; thereafter ini=1
      data ini/0/
c  no_nr=number of distinct nr values=no. of user coil data files; has max. allowed value of 10
c  Each user coil data file has a unique nr, but multiple icoil values can use the same nr
c  and associated user-coil data files
      save nrvals,no_nr
      save rrefi,nlinesi
c
      save zni
c
      save abrni
      save bbrni
      save cbrni
      save dbrni
c
      save abtni
      save bbtni
      save cbtni
      save dbtni
c-----------------------------------------------------------------
c     Zero out saved arrays in first call to this routine
      iwrite=0
      if(ini.eq.0) then
      ini=1
      do i=1,10
      rrefi(i)=0.d0
      nlinesi(i)=0
      do n=1,maxlines
      zni(n,i)=0.d0
c
      abrni(n,i)=0.d0
      bbrni(n,i)=0.d0
      cbrni(n,i)=0.d0
      dbrni(n,i)=0.d0
c
      abtni(n,i)=0.d0
      bbtni(n,i)=0.d0
      cbtni(n,i)=0.d0
      dbtni(n,i)=0.d0
      enddo
      enddo
c  Zero out nrvals: array of distinct nr values
c  Zero out no_nr: number of distinct nr values=no. of user magnet data files
      no_nr=0
      do i=1,10
      nrvals(i)=0
      enddo
c
      endif      
c ---------------------------        
      if(init.gt.0) go to 1      
c  First check to see if the current input nr has been used in a previous call
c  and ensure that a user coil data set labeled by nr is read in only once. If
c  current nr has already been "used", return without reading user data file
      if(iwrite.gt.0) write(6,*) 'In gtype18 initialization nr=',nr
c
      do i=1,10
      if(nr.eq.nrvals(i)) then
c  Current nr has already been used in an initialization call and the data already stored
c  Return with no data file read in
      write(6,*) 'nr value already used'
      write(6,*) 'Return without reading in a user data file'
      return
      endif
c----
      enddo
c ------------- 
c  New nr: increment nr count
      no_nr=no_nr+1
      if(iwrite.gt.0) write(6,*) 'In gtype18 no_nr=',no_nr
c  Stop program if no_nr > 10
      if(no_nr.gt.10) then
      write(6,*) 'no_nr > 10 in gtype18- stopped'
      stop
      endif
c Save current nr value
      nrvals(no_nr)=nr
c Read in data from a user data file- single m value  
      call read_magnet_data(zn,abrn,bbrn,cbrn,
     &dbrn,abtn,bbtn,cbtn,dbtn,rref,nlines,mval,nr,iwrite)
c
      if(iwrite.gt.0) then
      write(6,*) 'no_nr=',no_nr,' rref,nlines=',rref,nlines
      endif
c Load user magnet data into saved arrays
      rrefi(no_nr)=rref
      nlinesi(no_nr)=nlines
c
      do n=1,nlines
      zni(n,no_nr)=zn(n)
c
      abrni(n,no_nr)=abrn(n)
      bbrni(n,no_nr)=bbrn(n)
      cbrni(n,no_nr)=cbrn(n)
      dbrni(n,no_nr)=dbrn(n)
c
      abtni(n,no_nr)=abtn(n)
      bbtni(n,no_nr)=bbtn(n)
      cbtni(n,no_nr)=cbtn(n)
      dbtni(n,no_nr)=dbtn(n)
      enddo
c
      return
c
c  End initialization
c--------------------------------------------------------------------------------
c  Compute gradient
    1 continue
      if(iwrite.gt.0) then
      write(6,*) 'saved no_nr=',no_nr
      write(6,*) 'saved rreffi,nlinesi'
      do i=1,no_nr
      write(6,*) 'i,rrefi(i),nlinesi(i)=',i,rrefi(i),nlinesi(i)
      enddo
      endif
c      if(init.eq.3) go to 4
c  Gradient computation, init > 0
c  Zero out dg and dgm arrrays
      do i=1,maxdrv
      dg(i)=0.d0
      enddo
c
      do i=0,maxdriv
      dgm(i)=0.d0
      enddo
c
c     Check to see if nr is one of the values of data files read in and stored
      do i=1,10
      if(nr.eq.nrvals(i)) go to 3
      enddo
      write(6,*) 'data file for current nr not stored in gtype18'
      write(6,*) 'stopped'
      stop
 3    continue
c  Have found nr value in list of stored nr values of data files read in.
c  Retrieve coil data from arrays with index i
      rref=rrefi(i)
      nlines=nlinesi(i)
c
      do n=1,nlines
      zn(n)=zni(n,i)
c
      abrn(n)=abrni(n,i)
      bbrn(n)=bbrni(n,i)
      cbrn(n)=cbrni(n,i)
      dbrn(n)=dbrni(n,i)
c
      abtn(n)=abtni(n,i)
      bbtn(n)=bbtni(n,i)
      cbtn(n)=cbtni(n,i)
      dbtn(n)=dbtni(n,i)
      enddo    
c  ndrivs is for zero index, ndriv for one index
c  dg has 1,2,...ndriv, dgm has 0,1,2,...ndrivs
c  dg(1)=dgm(0) dg(ndriv)=dgm(ndriv-1) etc.
      ndrivs=ndriv-1
      call dg_user_magnet(zn,abrn,bbrn,cbrn,dbrn,
     &abtn,bbtn,cbtn,dbtn,rref,mval,nlines,z,dgm,ndrivs,iwrite)
c     Index shift and sign flip. Also multiply dg by scaling factor g.
c  For gtypes 1-17 glprod has units of tesla/meter**(m-2) but for type 17 glprod has units of
c  meters
c  Note gg as defined below is dimensionless since dgm is derived from field data
c  (units of teslas) and has dimensions of tesla/meter**(m-2): T/m for quads, etc.
c  Negative sign in in gg is added here to make the sign of the type 18 gradient
c  consistent with that of types 1-17.
      gg=-glprod/xl
      do n=1,ndriv
      dg(n)=gg*dgm(n-1)
      enddo
      return
c
      end
c      
c End subroutine gtype18
c*********************************************************************************
c
c
          
