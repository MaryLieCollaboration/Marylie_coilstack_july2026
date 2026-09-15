************************************************************************
* header:                 BOOKKEEP
*  Index manipulation, table creation.
************************************************************************
c
      subroutine binom5
c
c     computes the table of the binomial coefficients
c
      implicit double precision (a-h,o-z)
      integer bin5(24,20)
      common /bin5/ bin5
      save/bin5/
      do 1 i=1,20
      bin5(i,1)=i
      do 1 k=2,20
      if (i-k) 2,3,4
    2 bin5(i,k)=0
      go to 1
    3 bin5(i,k)=1
      go to 1
    4 ip=i-1
      kp=k-1
      bin5(i,k)=bin5(ip,kp)+bin5(ip,k)
    1 continue
      return
      end
c
      integer function ndexn5(num,j)
c  calculates the index of the long vector (dimension 209)
c  for a given set of short indices that indicate a power of
c  each variable, stored in the array j.
c
c  'bin5(i,k)' is the binomial coefficient i select k (i>k)
      integer bin5(24,20),j(6)
      common/bin5/ bin5
c  calculate the index
      n=j(num)
      m=j(num)-1
      do 100 i=2,num
      ib=num+1-i
      m=m+j(ib)
      ib=m+i
      n=n+bin5(ib,i)
  100 continue
      ndexn5=n
      return
      end
c
      integer function ndex5(j)
      integer j(6)
      ndex5=ndexn5(6,j)
      return
      end
c
c
***********************************************************************
c
      subroutine tables5
c  This subroutine creates two tables used in bookkeeping:
c    1) expon5(ind,psv) is the exponent of phase space variable
c       'psv' (1 to 6)  for monomial index 'ind' (1 to top)
c    2) vblist5(ind,vnum) is the variable list for each
c       monomial index number 'ind' (1 to top).  For nth
c       order terms, there will be n non-zero variable
c       numbers.
c  For example, monomial number 109 is X.PX.PX.Pt.
c  Thus, expon5(109,1to6)=1,2,0,0,0,1 and vblist5(109,1to4)=1,2,2,6.
c
c  ------Variables produced by subroutine-------
      include 'expon5.inc'
      include 'vblist5.inc'
      include 'prodex5.inc'
c  ------Variables used internally--------
      integer carry,ind,lnzj,psv,vnum
c  ord, pwr = order, exponent
      integer ord,pwr
c  indpr= index for product
      integer indpr
      include 'lims.inc'
c  j = array of exponents
      integer j(6)
      save j
      data j/6*0/
c
      call binom5
c
c-----------------------------------------------------
c  Sequentially create exponent table & calculate order & rearrangements
      do 150 ind=1, 923
        carry=j(6)
        j(6)=0
        lnzj=0
        do 100 psv=1,5
          if (j(psv).gt.0) lnzj=psv
  100   continue
        if (lnzj.gt.0) j(lnzj)=j(lnzj)-1
        j(lnzj+1)=j(lnzj+1)+1+carry
        ord=0
        do 120 psv=1,6
          pwr=j(psv)
          expon5(psv,ind)=pwr
          ord=ord+pwr
  120   continue
c
c-------------------------------------------------------
c  Create variable list table, using exponent table
        vnum=1
        do 220 psv=1,6
          do 200 k=1,j(psv)
            vblist5(vnum,ind)=psv
            vnum=vnum+1
  200     continue
  220   continue
  150 continue
c
c---------------------------------------------------------
c  Create product table, based on the idea of F. Neri and the
c   method of C. Iselin.
      do 380 psv=1,6
  380   prodex5(psv,0)=psv
c      do 300 ord=1,ordcat-1
      do 300 ord=1,5
        do 320 psv=1,6
          indpr=bottom(ord+1)
          do 340 ind=bottom(ord),top(ord)
  360       if (expon5(psv,indpr).eq.0) then
              indpr=indpr+1
              if (indpr.le.top(ord+1)) goto 360
            endif
            prodex5(psv,ind)=indpr
            indpr=indpr+1
  340     continue
  320   continue
  300 continue
c
      return
      end
c
***********************************************************************
c
      integer function iprod(i1,i2)
	  implicit none
c
	  integer i1, i2, j(6), i, ndex5
	  external ndex5
      include 'expon5.inc'
	  do 100 i = 1, 6
	    j(i) = expon5(i,i1) + expon5(i,i2)
 100  continue
	  iprod = ndex5(j)
	  return
	  end
c
c End of file
c

