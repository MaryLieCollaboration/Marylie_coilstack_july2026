************************************************************************
*  header                    XTRA CODE                                 *
*   MARYLIE code that is not self contained.                           *
************************************************************************
*
      subroutine cputim(t)
c-----------------------------------------------------------------------
c      implicit real(a-h,o-z),logical(l),integer(i,j,k,m,n)
c-----------------------------------------------------------------------
c t is the cpu time used so far
c note, that this routine will work only for VAX VMS-Systems
c written by Petra Schuett, 12 April 1988
c-----------------------------------------------------------------------
c      include '($jpidef)'
c      integer*2 ilen1,ilen2,ilen3,icod1,icod2,icod3
c      integer*4 iadd1,iadd11,iadd2,iadd21
c      common/tcbl/ilen1,icod1,iadd1,iadd11,
c     *            ilen2,icod2,iadd2,iadd21,
c     *            ilen3,icod3
c      data ilen1,ilen2,ilen3/3*4/
c     *     iadd11,iadd21,icod3/3*0/
c      data icod1/jpi$_cpulim/
c     *     icod2/jpi$_cputim/
c
c      iadd1 = %loc(icpulm)
c      iadd2 = %loc(icputm)
c
c      call sys$getjpi(,,,ilen1,,,)
c      t = icputm/100.
c 100  return
      end
c
***********************************************************************
c
      subroutine myexit
c
c  Written by Rob Ryne ca 1984
c
c Exit the program (may be replaced by STOP, if exit is not available).
c This is the only place in MaryLie where exit is called.
c
      call exit
      end
c
**************************************************************
c
      subroutine mytime(p)
c
c This subroutine computes and prints out the run time.
c It calls the VAX specific routine cputim.  This is the
c only place in MaryLie where cputim is called.
c Written by A. Dragt on 7 April 1988.  Modified
c 28 September 1991.
c
      include 'impli.inc'
      include 'time.inc'
      include 'files.inc'
c
      dimension p(6)
      real runtim,tnow
      tnow = 0.0
c
c Compute and write out run time.
c
      iopt = nint(p(1))
      ifile = nint(p(2))
      isend = nint(p(3))
c
      if (isend .ne. 0) then
      call cputim(tnow)
      runtim = tnow-tstart
      if ((isend .eq. 1) .or. (isend .eq. 3))
     # write (jof,100) runtim
      if ((isend .eq. 2) .or. (isend .eq. 3))
     # write (ifile,100) runtim
  100 format(/,1x,'Execution time in seconds = ',f8.2)
      endif
c
c Reset clock if iopt .ne. 0
c
      if (iopt .ne. 0) then
      call cputim(tstart)
      if ((isend .eq. 1) .or. (isend .eq. 3))
     # write (jof,*) ' time reset to zero'
      if ((isend .eq. 2) .or. (isend .eq. 3))
     # write (ifile,*) ' time reset to zero'
      endif
c
      return
      end
c
c end of file
 
