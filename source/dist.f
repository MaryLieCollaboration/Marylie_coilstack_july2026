*******************************************************************************
* header:    DIST                                                             *
*  A collection of programs for generating, modifying, and                    *
*  evaluating particle distributions                                          *
*******************************************************************************
      subroutine bgen(p)
c
c subroutine corresponding to the type code bgen
c written by Alex Dragt 6/14/91
c
      include 'impli.inc'
      include 'param.inc'
      include 'buffer.inc'
      include 'files.inc'
      include 'parset.inc'
c
c calling arrays
      dimension p(6)
      character*3 kynd
c
c working arrays
c use buf1 and buf2: put numerical moments in buf1a and analytic
c moments in buf2a
c use buf3 as working space to store desired eigen emittances
c
c set up control indices
c
      job=nint(p(1))
      iopt=nint(p(2))
      nray=nint(p(3))
      iseed=nint(p(4))
      isend=nint(p(5))
      ipset=nint(p(6))
c
c read contents of pset ipset
c
       ninf=nint(pst(1,ipset))
       nopt=nint(pst(2,ipset))
       nskip=nint(pst(3,ipset))
       nof1=nint(pst(4,ipset))
       nof2=nint(pst(5,ipset))
       sigmax=pst(6,ipset)
c
c clear buf1m,buf2m
c
      call mclear(buf1m)
      call mclear(buf2m)
c
c select job
c
      if (job .eq. 0) then
	  if (iopt .ne. 6) then
c compute moments of particle distribution
      call cmom(buf1a)
c write moments on file nof1
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
	  else
c compute moments of particle distribution including <5> and <6>
      call cmom5(buf1a)
c write moments on file nof1
      mpt=mpo
      mpo=nof1
      call mapout5(0,buf1a,buf1m)
      mpo=mpt
	  endif
	  endif
c
	  if (job .eq. 0) then
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
c
      if (job .ne. 0) then
c generate various distributions and/or moments
c
c get eigen emittances
c test for file read or internal map fetch
      if(ninf .lt. 0) then
      nmap=-ninf
      kynd='gtm'
      call strget(kynd,nmap,buf3a,buf3m)
      else
      mpit=mpi
      mpi=ninf
      call mapin(nopt,nskip,buf3a,buf3m)
      mpi=mpit
      endif
      x2mom=buf3a(7)
      y2mom=buf3a(18)
      t2mom=buf3a(25)
      endif
c use this information to generate distributions and/or compute moments
c
c uniformly filled ellipse or ellipsoids
      if (job .eq. 1) then
c set up scaling factors
      sx=sqrt(4.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call re2d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmre2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call re2d(nray,iseed,sx,sy,st)
      call cmre2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call re2d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call re2d(nray,iseed,sx,sy,st)
      call cmre2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 2) then
c set up scaling factors
      sx=sqrt(6.d0*x2mom)
      sy=sqrt(6.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call re4d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmre4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call re4d(nray,iseed,sx,sy,st)
      call cmre4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call re4d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call re4d(nray,iseed,sx,sy,st)
      call cmre4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 3) then
c set up scaling factors
      sx=sqrt(8.d0*x2mom)
      sy=sqrt(8.d0*y2mom)
      st=sqrt(8.d0*t2mom)
      if(iopt .eq. 1) call re6d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmre6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call re6d(nray,iseed,sx,sy,st)
      call cmre6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call re6d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call re6d(nray,iseed,sx,sy,st)
      call cmre6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
c
c gaussian distributions
      if (job .eq. 4) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call rg2d(nray,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) then
      call cmrg2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rg2d(nray,iseed,sigmax,sx,sy,st)
      call cmrg2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rg2d(nray,iseed,sigmax,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rg2d(nray,iseed,sigmax,sx,sy,st)
      call cmrg2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 5) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=sqrt(y2mom)
      st=0.d0
      if(iopt .eq. 1) call rg4d(nray,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) then
      call cmrg4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rg4d(nray,iseed,sigmax,sx,sy,st)
      call cmrg4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rg4d(nray,iseed,sigmax,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rg4d(nray,iseed,sigmax,sx,sy,st)
      call cmrg4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 6) then
c set up scaling factors
      sx=sqrt(x2mom)
      sy=sqrt(y2mom)
      st=sqrt(t2mom)
      if(iopt .eq. 1) call rg6d(nray,iseed,sigmax,sx,sy,st)
      if(iopt .eq. 2) then
      call cmrg6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rg6d(nray,iseed,sigmax,sx,sy,st)
      call cmrg6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rg6d(nray,iseed,sigmax,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rg6d(nray,iseed,sigmax,sx,sy,st)
      call cmrg6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
c
c systematic uniform tori
      if (job .eq. 7) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call st2d(nray,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call st2d(nray,sx,sy,st)
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call st2d(nray,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call st2d(nray,sx,sy,st)
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 8) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call st4d(nray,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call st4d(nray,sx,sy,st)
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call st4d(nray,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call st4d(nray,sx,sy,st)
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 9) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=sqrt(2.d0*t2mom)
      if(iopt .eq. 1) call st6d(nray,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call st6d(nray,sx,sy,st)
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call st6d(nray,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call st6d(nray,sx,sy,st)
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
c
c random uniform tori
      if (job .eq. 10) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=0.d0
      st=0.d0
      if(iopt .eq. 1) call rt2d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rt2d(nray,iseed,sx,sy,st)
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rt2d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rt2d(nray,iseed,sx,sy,st)
      call cmst2d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 11) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call rt4d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rt4d(nray,iseed,sx,sy,st)
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rt4d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rt4d(nray,iseed,sx,sy,st)
      call cmst4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
      if (job .eq. 12) then
c set up scaling factors
      sx=sqrt(2.d0*x2mom)
      sy=sqrt(2.d0*y2mom)
      st=sqrt(2.d0*t2mom)
      if(iopt .eq. 1) call rt6d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call rt6d(nray,iseed,sx,sy,st)
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call rt6d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call rt6d(nray,iseed,sx,sy,st)
      call cmst6d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
c
c KV distribution in 4-D phase space
      if (job .eq. 13) then
c set up scaling factors
      sx=sqrt(4.d0*x2mom)
      sy=sqrt(4.d0*y2mom)
      st=0.d0
      if(iopt .eq. 1) call kv4d(nray,iseed,sx,sy,st)
      if(iopt .eq. 2) then
      call cmkv4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 12) then
      call kv4d(nray,iseed,sx,sy,st)
      call cmkv4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      mpo=mpt
      endif
      if(iopt .eq. 13) then
      call kv4d(nray,iseed,sx,sy,st)
      call cmom(buf1a)
      mpt=mpo
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      if(iopt .eq. 123) then
      call kv4d(nray,iseed,sx,sy,st)
      call cmkv4d(buf2a,sx,sy,st)
      mpt=mpo
      mpo=nof2
      call mapout(0,buf2a,buf1m)
      call cmom(buf1a)
      mpo=nof1
      call mapout(0,buf1a,buf1m)
      mpo=mpt
c print selected numerical moments at terminal and/or write on drop file
      if(isend .eq. 1 .or. isend .eq. 3) then
      write(jof,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jof,*) buf1a(7),buf1a(8),buf1a(13)
      write(jof,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jof,*) buf1a(18),buf1a(19),buf1a(22)
      write(jof,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jof,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      if(isend .eq. 2 .or. isend .eq. 3) then
      write(jodf,*) ' values of <x*x>, <x*px>, <px*px>:'
      write(jodf,*) buf1a(7),buf1a(8),buf1a(13)
      write(jodf,*) ' values of <y*y>, <y*py>, <py*py>:'
      write(jodf,*) buf1a(18),buf1a(19),buf1a(22)
      write(jodf,*) ' values of <t*t>, <t*pt>, <pt*pt>:'
      write(jodf,*) buf1a(25),buf1a(26),buf1a(27)
      endif
      endif
      endif
c
      return
      end
c
********************************************************************************
c
      subroutine cmre2d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 2D random
c  uniform distribution that fills a 2D ellipse in phase space
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/4.d0
      fa(13)=sx2/4.d0
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/24.d0)
      fa(140)=sx2*sx2*(1.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmre4d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 4D random
c  uniform distribution that fills a 4D ellipsoid in phase space
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/6.d0
      fa(13)=sx2/6.d0
      fa(18)=sy2/6.d0
      fa(22)=sy2/6.d0
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/16.d0)
      fa(90)=sx2*sx2*(1.d0/48.d0)
      fa(95)=sx2*sy2*(1.d0/48.d0)
      fa(99)=sx2*sy2*(1.d0/48.d0)
      fa(140)=sx2*sx2*(1.d0/16.d0)
      fa(145)=sx2*sy2*(1.d0/48.d0)
      fa(149)=sx2*sy2*(1.d0/48.d0)
      fa(175)=sy2*sy2*(1.d0/16.d0)
      fa(179)=sy2*sy2*(1.d0/48.d0)
      fa(195)=sy2*sy2*(1.d0/16.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmre6d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 6D random
c  uniform distribution that fills a 6D ellipsoid in phase space
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/8.d0
      fa(13)=sx2/8.d0
      fa(18)=sy2/8.d0
      fa(22)=sy2/8.d0
      fa(25)=st2/8.d0
      fa(27)=st2/8.d0
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/80.d0)
      fa(90)=sx2*sx2*(1.d0/80.d0)
      fa(95)=sx2*sy2*(1.d0/80.d0)
      fa(99)=sx2*sy2*(1.d0/80.d0)
      fa(102)=sx2*st2*(1.d0/80.d0)
      fa(104)=sx2*st2*(1.d0/80.d0)
      fa(140)=sx2*sx2*(3.d0/80.d0)
      fa(145)=sx2*sy2*(1.d0/80.d0)
      fa(149)=sx2*sy2*(1.d0/80.d0)
      fa(152)=sx2*st2*(1.d0/80.d0)
      fa(154)=sx2*st2*(1.d0/80.d0)
      fa(175)=sy2*sy2*(3.d0/80.d0)
      fa(179)=sy2*sy2*(1.d0/80.d0)
      fa(182)=sy2*st2*(1.d0/80.d0)
      fa(184)=sy2*st2*(1.d0/80.d0)
      fa(195)=sy2*sy2*(3.d0/80.d0)
      fa(198)=sy2*st2*(1.d0/80.d0)
      fa(200)=sy2*st2*(1.d0/80.d0)
      fa(205)=st2*st2*(3.d0/80.d0)
      fa(207)=st2*st2*(1.d0/80.d0)
      fa(209)=st2*st2*(3.d0/80.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg2d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 2D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(140)=sx2*sx2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg4d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 4D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
      fa(18)=sy2
      fa(22)=sy2
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(95)=sx2*sy2
      fa(99)=sx2*sy2
      fa(140)=sx2*sx2*(3.d0)
      fa(145)=sx2*sy2
      fa(149)=sx2*sy2
      fa(175)=sy2*sy2*(3.d0)
      fa(179)=sy2*sy2
      fa(195)=sy2*sy2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmrg6d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 6D random
c  gaussian distribution
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2
      fa(13)=sx2
      fa(18)=sy2
      fa(22)=sy2
      fa(25)=st2
      fa(27)=st2
c  quartic moments
      fa(84)=sx2*sx2*(3.d0)
      fa(90)=sx2*sx2
      fa(95)=sx2*sy2
      fa(99)=sx2*sy2
      fa(102)=sx2*st2
      fa(104)=sx2*st2
      fa(140)=sx2*sx2*(3.d0)
      fa(145)=sx2*sy2
      fa(149)=sx2*sy2
      fa(152)=sx2*st2
      fa(154)=sx2*st2
      fa(175)=sy2*sy2*(3.d0)
      fa(179)=sy2*sy2
      fa(182)=sy2*st2
      fa(184)=sy2*st2
      fa(195)=sy2*sy2*(3.d0)
      fa(198)=sy2*st2
      fa(200)=sy2*st2
      fa(205)=st2*st2*(3.d0)
      fa(207)=st2*st2
      fa(209)=st2*st2*(3.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst2d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 2D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst4d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 4D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
      fa(18)=sy2/2.d0
      fa(22)=sy2/2.d0
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(95)=sx2*sy2*(1.d0/4.d0)
      fa(99)=sx2*sy2*(1.d0/4.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/4.d0)
      fa(149)=sx2*sy2*(1.d0/4.d0)
      fa(175)=sy2*sy2*(3.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/8.d0)
      fa(195)=sy2*sy2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmst6d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 6D systematic
c  distribution on a torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
      st2=st*st
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/2.d0
      fa(13)=sx2/2.d0
      fa(18)=sy2/2.d0
      fa(22)=sy2/2.d0
      fa(25)=st2/2.d0
      fa(27)=st2/2.d0
c  quartic moments
      fa(84)=sx2*sx2*(3.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/8.d0)
      fa(95)=sx2*sy2*(1.d0/4.d0)
      fa(99)=sx2*sy2*(1.d0/4.d0)
      fa(102)=sx2*st2*(1.d0/4.d0)
      fa(104)=sx2*st2*(1.d0/4.d0)
      fa(140)=sx2*sx2*(3.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/4.d0)
      fa(149)=sx2*sy2*(1.d0/4.d0)
      fa(152)=sx2*st2*(1.d0/4.d0)
      fa(154)=sx2*st2*(1.d0/4.d0)
      fa(175)=sy2*sy2*(3.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/8.d0)
      fa(182)=sy2*st2*(1.d0/4.d0)
      fa(184)=sy2*st2*(1.d0/4.d0)
      fa(195)=sy2*sy2*(3.d0/8.d0)
      fa(198)=sy2*st2*(1.d0/4.d0)
      fa(200)=sy2*st2*(1.d0/4.d0)
      fa(205)=st2*st2*(3.d0/8.d0)
      fa(207)=st2*st2*(1.d0/8.d0)
      fa(209)=st2*st2*(3.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine cmkv4d(fa,sx,sy,st)
c  subroutine to compute analytic moments of a 4-variable
c  KV distribution
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'param.inc'
c
c  calling arrays
      dimension fa(monoms)
c
c  procedure
c
c  clear array
      do 10 i=1,monoms
   10 fa(i)=0.d0
c  compute squares
      sx2=sx*sx
      sy2=sy*sy
c  compute nonzero moments
c  quadratic moments
      fa(7)=sx2/4.d0
      fa(13)=sx2/4.d0
      fa(18)=sy2/4.d0
      fa(22)=sy2/4.d0
c  quartic moments
      fa(84)=sx2*sx2*(1.d0/8.d0)
      fa(90)=sx2*sx2*(1.d0/24.d0)
      fa(95)=sx2*sy2*(1.d0/24.d0)
      fa(99)=sx2*sy2*(1.d0/24.d0)
      fa(140)=sx2*sx2*(1.d0/8.d0)
      fa(145)=sx2*sy2*(1.d0/24.d0)
      fa(149)=sx2*sy2*(1.d0/24.d0)
      fa(175)=sy2*sy2*(1.d0/8.d0)
      fa(179)=sy2*sy2*(1.d0/24.d0)
      fa(195)=sy2*sy2*(1.d0/8.d0)
c
      return
      end
c
********************************************************************************
c
      subroutine re2d(nray,jseed,sx,sy,st)
c
c  subroutine to generate a 2D random uniform distribution that fills
c  a 2D ellipse in phase space
c  written by Alex Dragt 7/6/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray
    5 sumsq=0.d0
      do 10 i=1,2
      call myrand(jseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=0.d0
      zblock(j,4)=0.d0
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine re4d(nray,jseed,sx,sy,st)
c
c  subroutine to generate a 4D random uniform distribution that fills
c  a 4D ellipsoid in phase space
c  written by Alex Dragt 7/6/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray
    5 sumsq=0.d0
      do 10 i=1,4
      call myrand(jseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine re6d(nray,jseed,sx,sy,st)
c
c  subroutine to generate a 6D random uniform distribution that fills
c  a 6D ellipsoid in phase space
c  written by Alex Dragt 7/6/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(6)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray
    5 sumsq=0.d0
      do 10 i=1,6
      call myrand(jseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=st*z(5)
      zblock(j,6)=st*z(6)
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine rg2d(nray,jseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 2D gaussian distribution
c  written by Alex Dragt 7/10/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      sigm2=sigmax**2
      do 20 j=1,nray
c generate a ray
   25 call normdv(z,1,jseed)
      ssq = z(1)**2 +z(2)**2
      if(ssq .gt. sigm2) go to 25
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=0.d0
      zblock(j,4)=0.d0
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine rg4d(nray,jseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 4D gaussian distribution
c  written by Alex Dragt 7/10/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      sigm2=sigmax**2
      do 20 j=1,nray
c generate a ray
   25 call normdv(z,2,jseed)
      ssq = z(1)**2 + z(2)**2 + z(3)**2 + z(4)**2
      if(ssq .gt. sigm2) go to 25
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine rg6d(nray,jseed,sigmax,sx,sy,st)
c
c  subroutine to generate a 6D gaussian distribution
c  written by Alex Dragt 7/10/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(6)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      sigm2=sigmax**2
      do 20 j=1,nray
c generate a ray
   25 call normdv(z,3,jseed)
      ssq=z(1)**2+z(2)**2+z(3)**2+z(4)**2+z(5)**2+z(6)**2
      if(ssq .gt. sigm2) go to 25
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=st*z(5)
      zblock(j,6)=st*z(6)
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
******************************************************************************
c
      subroutine rt2d(nray,iseed,sx,sy,st)
c  subroutine to generate a 2D random uniform distribution on a 2-torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray on the unit torus
      call myrand(iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=0.d0
      zblock(j,4)=0.d0
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine rt4d(nray,iseed,sx,sy,st)
c  subroutine to generate a 4D random uniform distribution on a 4-torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray on the unit torus
      call myrand(iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
      call myrand(iseed,ans)
      arg = twopi*ans
      z(3) = cos(arg)
      z(4) = sin(arg)
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine rt6d(nray,iseed,sx,sy,st)
c  subroutine to generate a 6D random uniform distribution on a 6-torus
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(6)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
      twopi = 8.d0*atan(1.d0)
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray on the unit torus
      call myrand(iseed,ans)
      arg = twopi*ans
      z(1) = cos(arg)
      z(2) = sin(arg)
      call myrand(iseed,ans)
      arg = twopi*ans
      z(3) = cos(arg)
      z(4) = sin(arg)
      call myrand(iseed,ans)
      arg = twopi*ans
      z(5) = cos(arg)
      z(6) = sin(arg)
c scale and store ray in zblock
      zblock(j,1)=sx*z(1)
      zblock(j,2)=sx*z(2)
      zblock(j,3)=sy*z(3)
      zblock(j,4)=sy*z(4)
      zblock(j,5)=st*z(5)
      zblock(j,6)=st*z(6)
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
********************************************************************************
c
      subroutine st2d(nray,sx,sy,st)
c  subroutine to generate a 2D systematic uniform distribution on a 2-torus
c written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(2)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on circle
      npts=nray
c
      pidiv = twopi/float(npts)
      ntot = 0
      do 20 i=1,npts
c generate ray on unit 2-torus
      ai = float(i-1)
      arg= ai*pidiv
      z(1) = cos(arg)
      z(2) = sin(arg)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(ntot,1)=sx*z(1)
      zblock(ntot,2)=sx*z(2)
      zblock(ntot,3)=0.d0
      zblock(ntot,4)=0.d0
      zblock(ntot,5)=0.d0
      zblock(ntot,6)=0.d0
   20 continue
      nrays=ntot
      write(6,*) ntot, ' rays generated'
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
      return
      end
c
********************************************************************************
c
      subroutine st4d(nray,sx,sy,st)
c  subroutine to generate a 4D systematic uniform distribution on a 4-torus
c written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(4)
      dimension c(100),s(100)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      if(maxray .gt. 10000) write(6,*)
     # 'error: maxray exceeds 10000, subroutines st* need rewriting'
      return
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on each circle
      pow=1.d0/2.d0
      aray=float(nray)
      npts=int(aray**pow)
      if (npts.gt.100) then
      write(6,*)
     # 'error: maxray exceeds 10000, subroutines st* need rewriting'
      return
      endif
c
c  compute needed sines and cosines
c
      pidiv = twopi/float(npts)
      do 10 i=1,npts
      ai = float(i-1)
      arg= ai*pidiv
      c(i) = cos(arg)
      s(i) = sin(arg)
   10 continue
c
c generate rays on unit 4-torus
c
      ntot = 0
      do 20 i=1,npts
      z(1) = c(i)
      z(2) = s(i)
      do 30 j=1,npts
      z(3) = c(j)
      z(4) = s(j)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(ntot,1)=sx*z(1)
      zblock(ntot,2)=sx*z(2)
      zblock(ntot,3)=sy*z(3)
      zblock(ntot,4)=sy*z(4)
      zblock(ntot,5)=0.d0
      zblock(ntot,6)=0.d0
   30 continue
   20 continue
c
      nrays=ntot
      write(6,*) ntot, ' rays generated'
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
      return
      end
c
********************************************************************************
c
      subroutine st6d(nray,sx,sy,st)
c  subroutine to generate a 6D systematic uniform distribution on a 6-torus
c written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(6)
      dimension c(21),s(21)
c
c  proceedure
c
c  initialize constants, indices, and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      if(maxray .gt. 10000) write(6,*)
     # 'error: maxray exceeds 10000, subroutines st* need rewriting'
      return
      endif
      twopi = 8.d0*atan(1.d0)
c
c fill zblock
c
c compute number of points on each circle
      pow=1.d0/3.d0
      aray=float(nray)
      npts=int(aray**pow)
      if (npts.gt.21) then
      write(6,*)
     # 'error: maxray exceeds 10000, subroutines st* need rewriting'
      return
      endif
c
c  compute needed sines and cosines
c
      pidiv = twopi/float(npts)
      do 10 i=1,npts
      ai = float(i-1)
      arg= ai*pidiv
      c(i) = cos(arg)
      s(i) = sin(arg)
   10 continue
c
c generate rays on unit 6-torus
c
      ntot = 0
      do 20 i=1,npts
      z(1) = c(i)
      z(2) = s(i)
      do 30 j=1,npts
      z(3) = c(j)
      z(4) = s(j)
      do 40 k=1,npts
      z(5) = c(k)
      z(6) = s(k)
      ntot = ntot + 1
c scale and store ray in zblock
      zblock(ntot,1)=sx*z(1)
      zblock(ntot,2)=sx*z(2)
      zblock(ntot,3)=sy*z(3)
      zblock(ntot,4)=sy*z(4)
      zblock(ntot,5)=st*z(5)
      zblock(ntot,6)=st*z(6)
   40 continue
   30 continue
   20 continue
c
      nrays=ntot
      write(6,*) ntot, ' rays generated'
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
      return
      end
c
***************************************************************************
c
      subroutine kv4d(nray,jseed,sx,sy,st)
c subroutine to compute a KV distribution in 4-D phase space
c
c  written by Alex Dragt 7/15/91
c
      include 'impli.inc'
      include 'rays.inc'
c
c  working arrays
      dimension z(4)
c
c  proceedure
c
c  initialize indices and arrays
      if(nray .gt. maxray) then
      write(6,*) 'error: nray exceeds maxray'
      return
      endif
      nrays=nray
c
c  initialize arrrays and set up counters
c
      do 115 k=1,nrays
      istat(k)=0
      ihist(k,1)=0
      ihist(k,2)=0
  115 continue
      iturn=0
      nlost=0
c
c fill zblock
c
      do 20 j=1,nray
c generate a ray
    5 sumsq=0.d0
      do 10 i=1,4
      call myrand(jseed,ans)
      z(i)=2.d0*ans-1.d0
   10 sumsq=sumsq+z(i)*z(i)
      if(sumsq .gt. 1.d0) goto 5
      if(sumsq .eq. 0.d0) goto 5
c scale and store ray in zblock
      rad=sqrt(sumsq)
      zblock(j,1)=sx*z(1)/rad
      zblock(j,2)=sx*z(2)/rad
      zblock(j,3)=sy*z(3)/rad
      zblock(j,4)=sy*z(4)/rad
      zblock(j,5)=0.d0
      zblock(j,6)=0.d0
   20 continue
      write(6,*) nray, ' rays generated'
c
      return
      end
c
***********************************************************************
c
      subroutine myrand(iseed,ans)
c
c this subroutine generates random numbers using either myran1
c or myran3
c written by Alex Dragt, 7/28/91
c
      double precision ans
c
c      write(6,*) iseed
c
      if(iseed .ge. 1) call myran3(iseed,ans)
      if(iseed .eq. 0) then
      write(6,*) 'Error: iseed=0'
      call myexit
      endif
      if(iseed .le. -1) call myran1(iseed,ans)
c
      return
      end
c
***********************************************************************
c
      subroutine myran1(IDUM,ans)
c This subroutine is a minor modification of the FUNCTION RAN1(IDUM)
c given in the book "Numerical Recipes" by W. Press, B. Flannery,
c S. Teukolsky, and W. Vetterling.
c It requires seeds IDUM that are lt -1 to be reset.
c Written by Alex Dragt 7/16/91.
      DIMENSION R(97)
      save IFF,IX1,IX2,IX3,R
      double precision ans
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
      DATA IFF /0/
c
c      write(6,*) IDUM
c
      IF (IDUM.LT.-1.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=-1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1) then
      write (6,*) 'J in myran1 is out of range'
      call myexit
      endif
      ans=dble(R(J))
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
c
      RETURN
      END
c
*************************************************************************
c
      subroutine myran3(IDUM,ans)
c This subroutine is a minor modification of the FUNCTION RAN3(IDUM)
c given in the book "Numerical Recipes" by W. Press, B. Flannery,
c S. Teukolsky, and W. Vetterling.
c It requires seeds IDUM that are gt +1 to be reset.
c written by Alex Dragt 7/16/91.
C         IMPLICIT REAL*4(M)
C         PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=2.5E-7)
      DIMENSION MA(55)
      save IFF,INEXT,INEXTP,MA
      double precision ans
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
      DATA IFF /0/
c
c      write(6,*) IDUM
c
      IF(IDUM.gt.1.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED+IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
          II=MOD(21*I,55)
          MA(II)=MK
          MK=MJ-MK
          IF(MK.LT.MZ)MK=MK+MBIG
          MJ=MA(II)
11      CONTINUE
        DO 13 K=1,4
          DO 12 I=1,55
            MA(I)=MA(I)-MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
12        CONTINUE
13      CONTINUE
        INEXT=0
        INEXTP=31
        IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      ans=dble(MJ*FAC)
c
      RETURN
      END
c
********************************************************************************
c
      subroutine normdv(x,n,jseed)
c  routine to generate 2n independent normal deviates
c  Polar method for normal deviates
c  written by Alex Dragt ca 1986, revised by Alex Dragt 7/15/91
c
c  References:
c  1)Knuth, The Art of Computer Programming (Vol. 2, page 117)
c  2)Irving Haber, NRL Memo Report #3705
c
      include 'impli.inc'
c
c calling arrays
      dimension x(*)
c
c procedure
c
      do 100 k=1,n
   50 call myrand(jseed,ans)
      u1=ans
      call myrand(jseed,ans)
      u2=ans
      v1=2.d0*u1-1.d0
      v2=2.d0*u2-1.d0
      s=v1**2 + v2**2
      if(s .ge. 1.d0)goto 50
      arg=sqrt(-2.d0*log(s)/s)
      ans1=v1*arg
      ans2=v2*arg
      x(k)=ans1
      x(k+n)=ans2
  100 continue
c
      return
      end
c
***********************************************************************
c
      subroutine tic(p)
c
c Translation of initial conditions.
c This subroutine produces a translation in 6-dimensional phase space.
c The parameters p(j) are used to specify translations deltaz(j)
c according to the relations deltaz(j)=p(j).
c The suffixes 'i' and 'f' refer to 'initial' and 'final' respectively.
c Written by Alex Dragt, Fall 1986
c
      include 'impli.inc'
      include 'rays.inc'
      dimension p(6)
c
      do 100 i=1,nrays
c check to see if i'th ray has been lost
      if (istat(i).ne.0) goto 100
c if not, copy the i'th ray out of zblock
      do 110 j=1,6
  110 zi(j)=zblock(i,j)
c
c transform this ray
c
      do 10 j=1,6
   10 zf(j)=zi(j)+p(j)
c
c put the transformed ray back into zblock
      do 120 j=1,6
  120 zblock(i,j)=zf(j)
c
  100 continue
      return
      end
c
c end of file
