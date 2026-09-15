c***********************************************************************************
c     August 9 2024 stand-alone user data file read: read_magnet_data+getfortname
c      July 24 2024 single m files change 2 to 3
c      subroutine read_magnet_data_2(zn_nm,abrn_nm,bbrn_nm,cbrn_nm,
c     &dbrn_nm,abtn_nm,bbtn_nm,cbtn_nm,dbtn_nm,rref_nm,mvals_nm,
c     &     nlines_nm,nmvals,nr,iwrite)
c
      subroutine read_magnet_data(zn,abrn,bbrn,cbrn,
     &dbrn,abtn,bbtn,cbtn,dbtn,rref,nlines,mval,nr,iwrite)
      implicit double precision(a-h,o-z)
c  July 24 2024 single m files with m check single data block. Used only in initialization calls
c  7/24/2024 single m, modified header
c  Reads in data file for a single user magnet
c  Data are not stored in this subroutine
c  Default Fortran read assumed. Read in file named fort.nr
c  (not literally fort.nr but fort.99 for nr=99, fort.100 for nr=100, etc.)
c  !! nr must be different for different user-type coils
c  Lines after the header are in pairs with z,a1,b1,c1,d1 followed by line with a2,b2,c2,d2
c  (such a pair of lines is called a double line)
c
c  Reads all lines in the file fort.nr as character strings of length 150- later convert
c  to numbers and shorter character strings as needed by means of Fortran reads from a
c  character string
c  Blank lines separate data blocks with various m values
c  The first line is header, which is followed by a single m data block
c---------------------------------------------------------------
c  July 24, 2024 single m
c  A single call to this subroutine reads in a field-data set for a single m user magnet
c
c  Cubic coefficients on the interval [z_n,z_(n+1)] are in the variable t=z-z_n such that
c  br_m(z)=abrn(n)+bbrn(n)*t+cbrn(n)*t**2+dbrn(n)*t**3
c  bphi_m(z)=abtn(n)+bbtn(n)*t+cbtn(n)*t**2+dbtn(n)*t**3
c  with m=mval (input mval is checked against data file m)
c
c  Input variables:
c  mval=coil m (can be negative to indicate skew coil). Checked against m in header (should be the same)
c  nr=integer that serves as both Fortran read number and unique user coil data set identifier
c  iwrite=integer code for diagnostic typeout: iwrite=1, typeout; iwrite.ne.1, no typeout
c
c  Output variables:
c
c  rref=reference radius for coil data
c  nlines=no. of double lines after header
c  zn=array of z values of data
c  abrn
c  bbrn
c  cbrn
c  dbrn
c  abtn
c  bbtn
c  cbtn
c  dbtn
c---------------------------------------------------------------------------------
     
      logical lexist
      character*150 line
      character*50 ifile
      character*150 blank
      character*150 header
      character*20 filename
c  input file name not used- read in default Fortran file fort.nr
c  fort.99, or fort.111, etc.
c
c  Scratch array of reference radii for a particular user coil
c  Note: Various m values in a particular coil can have different reference radii,
c     although usually the same reference radius is used for all m
c
c maxlines=max. no. of coefficient lines in a field data block
      parameter(maxlines=10000)
c Removed 7/24/2024 single m maxm is maximum number of m values in data arrays.
c      parameter(maxnmvals=10)
c------------------------------------------------------------
c etc. 7/24/2014 return single index arrays below (n steps through z values)
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
c---------------------------------------------------------
c  Check for existance of default Fortran file fort.nr (not literally fort.nr but
c  fort.99 for nr=99, fort.100 for
c  nr=100, etc.). nr must be an integer with either 2 or 3 digits and must be positive
c  !!Note: must not conflict with numbers in Marylie run shellscript for Fortran I/O
      if(iwrite.gt.0) then
      write(6,*) '!!!!!iwrite in read_user_field_data=',iwrite
      endif
c
      write(6,*) 'Enter read_magnet_data'
      call getfortname(nr,filename)
c
      write(6,*) 'In read_magnet_data after getfortname'
      if(iwrite.gt.0) then
      write(6,*) 'Input default Fortran file name=',filename
      endif
c
      inquire(file=filename,exist=lexist)
      if(lexist.eqv.(.false.)) then
      write(6,*) 'File ',filename,' not found'
      write(6,*) 'Stopped in read_magnet_data'
      stop
      endif
      blank=' '
c--------------------------------------------------------------
c     iwrite=1
c  Default Fortran read- e.g. expects to find file with name fort.99 if nr=99, etc.
c  Read header- first line of user data file
      read(nr,100,err=1000) header
      if(iwrite.gt.0) write(6,*) 'header=',header
c  Read from character string header
c  m=multipole index in user file
c     ntheta=nom. no. of integration steps in getting Fourier coeffs.
c     For info., not used here
c  rref=reference radius=radius of user coil data cylinder
c  ifile=typically the name of the raw user data file processed to get cubic. coeffs.
c     For info., not used here
      read(header,*,err=1001) m,rref,ntheta,ifile
      if(iwrite.gt.0) write(6,*) m,rref,ntheta,ifile
      write(6,*) 'm,rref,ntheta,ifile',m,rref,ntheta,ifile
c  Consistency check for input mval and m read from user data file
      if(mval.ne.m) then
      write(6,*) 'In read_magnet_data input mval.ne.data file m'
      write(6,*) 'Stopped'
      stop
      endif
      nlines=0
c  Read character string line
 1    read(nr,100,err=1002,end=2) line
 100  format(200a)
      nlines=nlines+1
c  Read from character string line. Double line 5/4 for a single z
      read(line,*,err=1003) z,a1,b1,c1,d1
c  Read character string line
      read(nr,100,err=1002,end=2) line
c  Read from character string line. 
      read(line,*,err=1003) a2,b2,c2,d2
      zn(nlines)=z
      abrn(nlines)=a1
      bbrn(nlines)=b1
      cbrn(nlines)=c1
      dbrn(nlines)=d1
c
      abtn(nlines)=a2
      bbtn(nlines)=b2
      cbtn(nlines)=c2
      dbtn(nlines)=d2     
      go to 1
c     
 2    continue
      if(iwrite.eq.1) write(6,*) 'Reached end of file'
      if(iwrite.eq.1) write(6,*) 'No. of data lines read in=',nlines
      if(iwrite.gt.0) then
      if(nlines.ge.400) then
      write(6,*) 'typ. data line pair(e.g. n=400)'
      write(6,*) zn(400),abrn(400),bbrn(400),cbrn(400),dbrn(400)
      write(6,*) abtn(400),bbtn(400),cbtn(400),dbtn(400)
      endif
      endif
c
c------------------------------------------------------------------
c     Debug write user data to file not all coeffs
c      if(iwrite.gt.0) then
c      open(unit=21,file='userdat.dat')
c      write(6,*) 'Open file userdat.dat'
c      do nm=1,nmvals
c      nlines=nlines_nm(nm)
c      do n=1,nlines
c     write(21,*) zn_nm(n,nm),abrn_nm(n,nm)
c      enddo
c      enddo
c      close(unit=21)
c      endif
      
c-----------------------------------------------------------------      
      if(iwrite.gt.0) then
      write(6,*) 'Inside read_magnet_data, returning'
      write(6,*) ' '
      endif
      return
c
c  Error messages
 1000 write(6,*) 'In read_magnet_data'
      write(6,*) 'Error reading 1st line (header) of the file- stopped'
      stop
 1001 write(6,*) 'In read_magnet_data'
      write(6,*) 'Error reading m,rref from character string- stopped'
      stop
 1002 write(6,*) 'In read_magnet_data'
      write(6,*) 'Error reading field-data line- stopped'
      stop
 1003 write(6,*) 'In read_magnet_data'
      write(6,*) 'Error reading coefficients from character string'
      write(6,*) 'stopped'
      stop
      end
c
c End subroutine read_magnet_data   
c***********************************************************************************
c
      subroutine getfortname(nr,filename)
      integer nr
      character*5 sname
      character*2 nrname
      character*3 nrrname
      character*7 name
      character*8 nname
      character*20 filename
c  Generates  file name filename=fort.nr 
c  This file name is used to check for existance input file with the
c  default Fortran name
c  nr can be 2 or 3 digits long and is positive
      sname='fort.'
      if(nr.lt.10) then
      write(6,*) 'nr=',nr,' is < 10- not allowed'
      write(6,*) 'stopped in getfortname'
      stop
      endif
      if(nr.gt.999) then
      write(6,*) 'nr=',nr,' is > 999- not allowed'
      write(6,*) 'stopped in getfortname'
      stop
      endif
c----------------------
      if(nr.lt.100) then
c      write(6,*) 'nr < 100'
      write(nrname,101) nr
 101  format(i2)
      name=sname//nrname
c      write(6,*) name
      filename=name
      go to 1
      endif
c------------------------------------------
      if(nr.lt.1000) then
      write(nrrname,102) nr
 102  format(i3)
      nname=sname//nrrname
c      write(6,*) nname
      filename=nname
      endif
c------------------------------------
 1    return
      end
c     End subroutine getfortname
c****************************************************************************
      
