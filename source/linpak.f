      subroutine dpofa(a,lda,n,info)
c***begin prologue  dpofa
c***date written   780814   (yymmdd)
c***revision date  861211   (yymmdd)
c***category no.  d2b1b
c***keywords  library=slatec(linpack),
c             type=double precision(spofa-s dpofa-d cpofa-c),
c             linear algebra,matrix,matrix factorization,
c             positive definite
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a double precision symmetric positive definite
c            matrix.
c***description
c
c     dpofa factors a double precision symmetric positive definite
c     matrix.
c
c     dpofa is usually called by dpoco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dpoco) = (1 + 18/n)*(time for dpofa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the symmetric matrix to be factored.  only the
c                diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix  r  so that  a = trans(r)*r
c                where  trans(r)  is the transpose.
c                the strict lower triangle is unaltered.
c                if  info .ne. 0 , the factorization is not complete.
c
c        info    integer
c                = 0  for normal return.
c                = k  signals an error condition.  the leading minor
c                     of order  k  is not positive definite.
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas ddot
c     fortran dsqrt
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  ddot
c***end prologue  dpofa
      integer lda,n,info
      double precision a(lda,1)
c
      double precision ddot,t
      double precision s
      integer j,jm1,k
c     begin block with ...exits to 40
c
c
c***first executable statement  dpofa
         do 30 j = 1, n
            info = j
            s = 0.0d0
            jm1 = j - 1
            if (jm1 .lt. 1) go to 20
            do 10 k = 1, jm1
               t = a(k,j) - ddot(k-1,a(1,k),1,a(1,j),1)
               t = t/a(k,k)
               a(k,j) = t
               s = s + t*t
   10       continue
   20       continue
            s = a(j,j) - s
c     ......exit
            if (s .le. 0.0d0) go to 40
            a(j,j) = dsqrt(s)
   30    continue
         info = 0
   40 continue
      return
      end
      subroutine dposl(a,lda,n,b)
c***begin prologue  dposl
c***date written   780814   (yymmdd)
c***revision date  861211   (yymmdd)
c***category no.  d2b1b
c***keywords  library=slatec(linpack),
c             type=double precision(sposl-s dposl-d cposl-c),
c             linear algebra,matrix,positive definite,solve
c***author  moler, c. b., (u. of new mexico)
c***purpose  solves the double precision symmetric positive definite
c            system a*x=b using the factors computed by dpoco or dpofa.
c***description
c
c     dposl solves the double precision symmetric positive definite
c     system a * x = b
c     using the factors computed by dpoco or dpofa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dpoco or dpofa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal.  technically this indicates
c        singularity, but it is usually caused by improper subroutine
c        arguments.  it will not occur if the subroutines are called
c        correctly and  info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dpoco(a,lda,n,rcond,z,info)
c           if (rcond is too small .or. info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dposl(a,lda,n,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot
c***end prologue  dposl
      integer lda,n
      double precision a(lda,1),b(1)
c
      double precision ddot,t
      integer k,kb
c
c     solve trans(r)*y = b
c
c***first executable statement  dposl
      do 10 k = 1, n
         t = ddot(k-1,a(1,k),1,b(1),1)
         b(k) = (b(k) - t)/a(k,k)
   10 continue
c
c     solve r*x = y
c
      do 20 kb = 1, n
         k = n + 1 - kb
         b(k) = b(k)/a(k,k)
         t = -b(k)
         call daxpy(k-1,t,a(1,k),1,b(1),1)
   20 continue
      return
      end
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
c***begin prologue  dqrdc
c***date written   780814   (yymmdd)
c***revision date  861211   (yymmdd)
c***category no.  d5
c***keywords  library=slatec(linpack),
c             type=double precision(sqrdc-s dqrdc-d cqrdc-c),
c             linear algebra,matrix,orthogonal triangular,
c             qr decomposition
c***author  stewart, g. w., (u. of maryland)
c***purpose  uses householder transformations to compute the qr factori-
c            zation of n by p matrix x.  column pivoting is optional.
c***description
c
c     dqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the user's option.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    double precision(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack.  this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dswap,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot,dnrm2,dscal,dswap
c***end prologue  dqrdc
      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
c
c***first executable statement  dqrdc
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end
      subroutine dqrsl(x,ldx,n,k,qraux,y,qy,qty,b,rsd,xb,job,info)
c***begin prologue  dqrsl
c***date written   780814   (yymmdd)
c***revision date  861211   (yymmdd)
c***category no.  d9,d2a1
c***keywords  library=slatec(linpack),
c             type=double precision(sqrsl-s dqrsl-d cqrsl-c),
c             linear algebra,matrix,orthogonal triangular,solve
c***author  stewart, g. w., (u. of maryland)
c***purpose  applies the output of dqrdc to compute coordinate
c            transformations, projections, and least squares solutions.
c***description
c
c     dqrsl applies the output of dqrdc to compute coordinate
c     transformations, projections, and least squares solutions.
c     for k .le. min(n,p), let xk be the matrix
c
c            xk = (x(jpvt(1)),x(jpvt(2)), ... ,x(jpvt(k)))
c
c     formed from columnns jpvt(1), ... ,jpvt(k) of the original
c     n x p matrix x that was input to dqrdc (if no pivoting was
c     done, xk consists of the first k columns of x in their
c     original order).  dqrdc produces a factored orthogonal matrix q
c     and an upper triangular matrix r such that
c
c              xk = q * (r)
c                       (0)
c
c     this information is contained in coded form in the arrays
c     x and qraux.
c
c     on entry
c
c        x      double precision(ldx,p).
c               x contains the output of dqrdc.
c
c        ldx    integer.
c               ldx is the leading dimension of the array x.
c
c        n      integer.
c               n is the number of rows of the matrix xk.  it must
c               have the same value as n in dqrdc.
c
c        k      integer.
c               k is the number of columns of the matrix xk.  k
c               must not be greater than min(n,p), where p is the
c               same as in the calling sequence to dqrdc.
c
c        qraux  double precision(p).
c               qraux contains the auxiliary output from dqrdc.
c
c        y      double precision(n)
c               y contains an n-vector that is to be manipulated
c               by dqrsl.
c
c        job    integer.
c               job specifies what is to be computed.  job has
c               the decimal expansion abcde, with the following
c               meaning.
c
c                    if a .ne. 0, compute qy.
c                    if b,c,d, or e .ne. 0, compute qty.
c                    if c .ne. 0, compute b.
c                    if d .ne. 0, compute rsd.
c                    if e .ne. 0, compute xb.
c
c               note that a request to compute b, rsd, or xb
c               automatically triggers the computation of qty, for
c               which an array must be provided in the calling
c               sequence.
c
c     on return
c
c        qy     double precision(n).
c               qy contains q*y, if its computation has been
c               requested.
c
c        qty    double precision(n).
c               qty contains trans(q)*y, if its computation has
c               been requested.  here trans(q) is the
c               transpose of the matrix q.
c
c        b      double precision(k)
c               b contains the solution of the least squares problem
c
c                    minimize norm2(y - xk*b),
c
c               if its computation has been requested.  (note that
c               if pivoting was requested in dqrdc, the j-th
c               component of b will be associated with column jpvt(j)
c               of the original matrix x that was input into dqrdc.)
c
c        rsd    double precision(n).
c               rsd contains the least squares residual y - xk*b,
c               if its computation has been requested.  rsd is
c               also the orthogonal projection of y onto the
c               orthogonal complement of the column space of xk.
c
c        xb     double precision(n).
c               xb contains the least squares approximation xk*b,
c               if its computation has been requested.  xb is also
c               the orthogonal projection of y onto the column space
c               of x.
c
c        info   integer.
c               info is zero unless the computation of b has
c               been requested and r is exactly singular.  in
c               this case, info is the index of the first zero
c               diagonal element of r and b is left unaltered.
c
c     the parameters qy, qty, b, rsd, and xb are not referenced
c     if their computation is not requested and in this case
c     can be replaced by dummy variables in the calling program.
c     to save storage, the user may in some cases use the same
c     array for different parameters in the calling sequence.  a
c     frequently occuring example is when one wishes to compute
c     any of b, rsd, or xb and does not need y or qty.  in this
c     case one may identify y, qty, and one of b, rsd, or xb, while
c     providing separate arrays for anything else that is to be
c     computed.  thus the calling sequence
c
c          call dqrsl(x,ldx,n,k,qraux,y,dum,y,b,y,dum,110,info)
c
c     will result in the computation of b and rsd, with rsd
c     overwriting y.  more generally, each item in the following
c     list contains groups of permissible identifications for
c     a single calling sequence.
c
c          1. (y,qty,b) (rsd) (xb) (qy)
c
c          2. (y,qty,rsd) (b) (xb) (qy)
c
c          3. (y,qty,xb) (b) (rsd) (qy)
c
c          4. (y,qy) (qty,b) (rsd) (xb)
c
c          5. (y,qy) (qty,rsd) (b) (xb)
c
c          6. (y,qy) (qty,xb) (b) (rsd)
c
c     in any group the value returned in the array allocated to
c     the group corresponds to the last member of the group.
c
c     linpack.  this version dated 08/14/78 .
c     g. w. stewart, university of maryland, argonne national lab.
c
c     dqrsl uses the following functions and subprograms.
c
c     blas daxpy,dcopy,ddot
c     fortran dabs,min0,mod
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,dcopy,ddot
c***end prologue  dqrsl
      integer ldx,n,k,job,info
      double precision x(ldx,1),qraux(1),y(1),qy(1),qty(1),b(1),rsd(1),
     1                 xb(1)
c
      integer i,j,jj,ju,kp1
      double precision ddot,t,temp
      logical cb,cqy,cqty,cr,cxb
c
c     set info flag.
c
c***first executable statement  dqrsl
      info = 0
c
c     determine what is to be computed.
c
      cqy = job/10000 .ne. 0
      cqty = mod(job,10000) .ne. 0
      cb = mod(job,1000)/100 .ne. 0
      cr = mod(job,100)/10 .ne. 0
      cxb = mod(job,10) .ne. 0
      ju = min0(k,n-1)
c
c     special action when n=1.
c
      if (ju .ne. 0) go to 40
         if (cqy) qy(1) = y(1)
         if (cqty) qty(1) = y(1)
         if (cxb) xb(1) = y(1)
         if (.not.cb) go to 30
            if (x(1,1) .ne. 0.0d0) go to 10
               info = 1
            go to 20
   10       continue
               b(1) = y(1)/x(1,1)
   20       continue
   30    continue
         if (cr) rsd(1) = 0.0d0
      go to 250
   40 continue
c
c        set up to compute qy or qty.
c
         if (cqy) call dcopy(n,y,1,qy,1)
         if (cqty) call dcopy(n,y,1,qty,1)
         if (.not.cqy) go to 70
c
c           compute qy.
c
            do 60 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 50
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qy(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qy(j),1)
                  x(j,j) = temp
   50          continue
   60       continue
   70    continue
         if (.not.cqty) go to 100
c
c           compute trans(q)*y.
c
            do 90 j = 1, ju
               if (qraux(j) .eq. 0.0d0) go to 80
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  t = -ddot(n-j+1,x(j,j),1,qty(j),1)/x(j,j)
                  call daxpy(n-j+1,t,x(j,j),1,qty(j),1)
                  x(j,j) = temp
   80          continue
   90       continue
  100    continue
c
c        set up to compute b, rsd, or xb.
c
         if (cb) call dcopy(k,qty,1,b,1)
         kp1 = k + 1
         if (cxb) call dcopy(k,qty,1,xb,1)
         if (cr .and. k .lt. n) call dcopy(n-k,qty(kp1),1,rsd(kp1),1)
         if (.not.cxb .or. kp1 .gt. n) go to 120
            do 110 i = kp1, n
               xb(i) = 0.0d0
  110       continue
  120    continue
         if (.not.cr) go to 140
            do 130 i = 1, k
               rsd(i) = 0.0d0
  130       continue
  140    continue
         if (.not.cb) go to 190
c
c           compute b.
c
            do 170 jj = 1, k
               j = k - jj + 1
               if (x(j,j) .ne. 0.0d0) go to 150
                  info = j
c           ......exit
                  go to 180
  150          continue
               b(j) = b(j)/x(j,j)
               if (j .eq. 1) go to 160
                  t = -b(j)
                  call daxpy(j-1,t,x(1,j),1,b,1)
  160          continue
  170       continue
  180       continue
  190    continue
         if (.not.cr .and. .not.cxb) go to 240
c
c           compute rsd or xb as required.
c
            do 230 jj = 1, ju
               j = ju - jj + 1
               if (qraux(j) .eq. 0.0d0) go to 220
                  temp = x(j,j)
                  x(j,j) = qraux(j)
                  if (.not.cr) go to 200
                     t = -ddot(n-j+1,x(j,j),1,rsd(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,rsd(j),1)
  200             continue
                  if (.not.cxb) go to 210
                     t = -ddot(n-j+1,x(j,j),1,xb(j),1)/x(j,j)
                     call daxpy(n-j+1,t,x(j,j),1,xb(j),1)
  210             continue
                  x(j,j) = temp
  220          continue
  230       continue
  240    continue
  250 continue
      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
c***begin prologue  dgeco
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  condition,double precision,factor,linear algebra,linpack,
c             matrix
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a double precision matrix by gaussian elimination
c            and estimates the condition of the matrix.
c***description
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  dasum,daxpy,ddot,dgefa,dscal
c***end prologue  dgeco
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c     compute 1-norm of a
c
c***first executable statement  dgeco
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
c***begin prologue  dgesl
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2a1
c***keywords  double precision,linear algebra,linpack,matrix,solve
c***author  moler, c. b., (u. of new mexico)
c***purpose  solves the double precision system  a*x=b or  trans(a)*x=b
c            using the factors computed by dgeco or dgefa.
c***description
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot
c***end prologue  dgesl
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
      double precision ddot,t
      integer k,kb,l,nm1
c***first executable statement  dgesl
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine dsico(a,lda,n,kpvt,rcond,z)
c***begin prologue  dsico
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2b1a
c***keywords  condition,double precision,factor,linear algebra,linpack,
c             matrix,symmetric
c***author  moler, c. b., (u. of new mexico)
c***purpose  factors a d.p. symmetric matrix by elimination with symmet-
c            ric pivoting and estimates the condition of the matrix.
c***description
c
c     dsico factors a double precision symmetric matrix by elimination
c     with symmetric pivoting and estimates the condition of the
c     matrix.
c
c     if  rcond  is not needed, dsifa is slightly faster.
c     to solve  a*x = b , follow dsico by dsisl.
c     to compute  inverse(a)*c , follow dsico by dsisl.
c     to compute  inverse(a) , follow dsico by dsidi.
c     to compute  determinant(a) , follow dsico by dsidi.
c     to compute  inertia(a), follow dsico by dsidi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the symmetric matrix to be factored.
c                only the diagonal and upper triangle are used.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     output
c
c        a       a block diagonal matrix and the multipliers which
c                were used to obtain it.
c                the factorization can be written  a = u*d*trans(u)
c                where  u  is a product of permutation and unit
c                upper triangular matrices, trans(u) is the
c                transpose of  u , and  d  is block diagonal
c                with 1 by 1 and 2 by 2 blocks.
c
c        kpvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack.  this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dsifa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,iabs,dsign
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  dasum,daxpy,ddot,dscal,dsifa
c***end prologue  dsico
      integer lda,n,kpvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
      double precision ak,akm1,bk,bkm1,ddot,denom,ek,t
      double precision anorm,s,dasum,ynorm
      integer i,info,j,jm1,k,kp,kps,ks
c
c     find norm of a using only upper half
c
c***first executable statement  dsico
      do 30 j = 1, n
         z(j) = dasum(j,a(1,j),1)
         jm1 = j - 1
         if (jm1 .lt. 1) go to 20
         do 10 i = 1, jm1
            z(i) = z(i) + dabs(a(i,j))
   10    continue
   20    continue
   30 continue
      anorm = 0.0d0
      do 40 j = 1, n
         anorm = dmax1(anorm,z(j))
   40 continue
c
c     factor
c
      call dsifa(a,lda,n,kpvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  u*d*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve u*d*w = e
c
      ek = 1.0d0
      do 50 j = 1, n
         z(j) = 0.0d0
   50 continue
      k = n
   60 if (k .eq. 0) go to 120
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         kp = iabs(kpvt(k))
         kps = k + 1 - ks
         if (kp .eq. kps) go to 70
            t = z(kps)
            z(kps) = z(kp)
            z(kp) = t
   70    continue
         if (z(k) .ne. 0.0d0) ek = dsign(ek,z(k))
         z(k) = z(k) + ek
         call daxpy(k-ks,z(k),a(1,k),1,z(1),1)
         if (ks .eq. 1) go to 80
            if (z(k-1) .ne. 0.0d0) ek = dsign(ek,z(k-1))
            z(k-1) = z(k-1) + ek
            call daxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
   80    continue
         if (ks .eq. 2) go to 100
            if (dabs(z(k)) .le. dabs(a(k,k))) go to 90
               s = dabs(a(k,k))/dabs(z(k))
               call dscal(n,s,z,1)
               ek = s*ek
   90       continue
            if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
            if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         go to 110
  100    continue
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/a(k-1,k)
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  110    continue
         k = k - ks
      go to 60
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(u)*y = w
c
      k = 1
  130 if (k .gt. n) go to 160
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 150
            z(k) = z(k) + ddot(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     1         z(k+1) = z(k+1) + ddot(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 140
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  140       continue
  150    continue
         k = k + ks
      go to 130
  160 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve u*d*v = y
c
      k = n
  170 if (k .eq. 0) go to 230
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. ks) go to 190
            kp = iabs(kpvt(k))
            kps = k + 1 - ks
            if (kp .eq. kps) go to 180
               t = z(kps)
               z(kps) = z(kp)
               z(kp) = t
  180       continue
            call daxpy(k-ks,z(k),a(1,k),1,z(1),1)
            if (ks .eq. 2) call daxpy(k-ks,z(k-1),a(1,k-1),1,z(1),1)
  190    continue
         if (ks .eq. 2) go to 210
            if (dabs(z(k)) .le. dabs(a(k,k))) go to 200
               s = dabs(a(k,k))/dabs(z(k))
               call dscal(n,s,z,1)
               ynorm = s*ynorm
  200       continue
            if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
            if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         go to 220
  210    continue
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = z(k)/a(k-1,k)
            bkm1 = z(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            z(k) = (akm1*bk - bkm1)/denom
            z(k-1) = (ak*bkm1 - bk)/denom
  220    continue
         k = k - ks
      go to 170
  230 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve trans(u)*z = v
c
      k = 1
  240 if (k .gt. n) go to 270
         ks = 1
         if (kpvt(k) .lt. 0) ks = 2
         if (k .eq. 1) go to 260
            z(k) = z(k) + ddot(k-1,a(1,k),1,z(1),1)
            if (ks .eq. 2)
     1         z(k+1) = z(k+1) + ddot(k-1,a(1,k+1),1,z(1),1)
            kp = iabs(kpvt(k))
            if (kp .eq. k) go to 250
               t = z(k)
               z(k) = z(kp)
               z(kp) = t
  250       continue
  260    continue
         k = k + ks
      go to 240
  270 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dsisl(a,lda,n,kpvt,b)
c***begin prologue  dsisl
c***date written   780814   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d2b1a
c***keywords  double precision,linear algebra,linpack,matrix,solve,
c             symmetric
c***author  bunch, j., (ucsd)
c***purpose  solves the double precision symmetric system
c            a*x=b  using the factors computed by dsifa.
c***description
c
c     dsisl solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dsifa.
c
c     on entry
c
c        a       double precision(lda,n)
c                the output from dsifa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        kpvt    integer(n)
c                the pivot vector from dsifa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero may occur if  dsico  has set rcond .eq. 0.0
c        or  dsifa  has set info .ne. 0  .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dsifa(a,lda,n,kpvt,info)
c           if (info .ne. 0) go to ...
c           do 10 j = 1, p
c              call dsisl(a,lda,n,kpvt,c(1,j))
c        10 continue
c
c     linpack.  this version dated 08/14/78 .
c     james bunch, univ. calif. san diego, argonne nat. lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran iabs
c***references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
c                 *linpack users  guide*, siam, 1979.
c***routines called  daxpy,ddot
c***end prologue  dsisl
      integer lda,n,kpvt(1)
      double precision a(lda,1),b(1)
c
      double precision ak,akm1,bk,bkm1,ddot,denom,temp
      integer k,kp
c
c     loop backward applying the transformations and
c     d inverse to b.
c
c***first executable statement  dsisl
      k = n
   10 if (k .eq. 0) go to 80
         if (kpvt(k) .lt. 0) go to 40
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 30
               kp = kpvt(k)
               if (kp .eq. k) go to 20
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
   20          continue
c
c              apply the transformation.
c
               call daxpy(k-1,b(k),a(1,k),1,b(1),1)
   30       continue
c
c           apply d inverse.
c
            b(k) = b(k)/a(k,k)
            k = k - 1
         go to 70
   40    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 2) go to 60
               kp = iabs(kpvt(k))
               if (kp .eq. k - 1) go to 50
c
c                 interchange.
c
                  temp = b(k-1)
                  b(k-1) = b(kp)
                  b(kp) = temp
   50          continue
c
c              apply the transformation.
c
               call daxpy(k-2,b(k),a(1,k),1,b(1),1)
               call daxpy(k-2,b(k-1),a(1,k-1),1,b(1),1)
   60       continue
c
c           apply d inverse.
c
            ak = a(k,k)/a(k-1,k)
            akm1 = a(k-1,k-1)/a(k-1,k)
            bk = b(k)/a(k-1,k)
            bkm1 = b(k-1)/a(k-1,k)
            denom = ak*akm1 - 1.0d0
            b(k) = (akm1*bk - bkm1)/denom
            b(k-1) = (ak*bkm1 - bk)/denom
            k = k - 2
   70    continue
      go to 10
   80 continue
c
c     loop forward applying the transformations.
c
      k = 1
   90 if (k .gt. n) go to 160
         if (kpvt(k) .lt. 0) go to 120
c
c           1 x 1 pivot block.
c
            if (k .eq. 1) go to 110
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               kp = kpvt(k)
               if (kp .eq. k) go to 100
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  100          continue
  110       continue
            k = k + 1
         go to 150
  120    continue
c
c           2 x 2 pivot block.
c
            if (k .eq. 1) go to 140
c
c              apply the transformation.
c
               b(k) = b(k) + ddot(k-1,a(1,k),1,b(1),1)
               b(k+1) = b(k+1) + ddot(k-1,a(1,k+1),1,b(1),1)
               kp = iabs(kpvt(k))
               if (kp .eq. k) go to 130
c
c                 interchange.
c
                  temp = b(k)
                  b(k) = b(kp)
                  b(kp) = temp
  130          continue
  140       continue
            k = k + 2
  150    continue
      go to 90
  160 continue
      return
      end
      subroutine dswap(n,dx,incx,dy,incy)
****begin prologue  dswap
****date written   791001   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d1a5
****keywords  blas,double precision,interchange,linear algebra,vector
****author  lawson, c. l., (jpl)
*           hanson, r. j., (snla)
*           kincaid, d. r., (u. of texas)
*           krogh, f. t., (jpl)
****purpose  interchange d.p. vectors
****description
*
*                b l a s  subprogram
*    description of parameters
*
*     --input--
*        n  number of elements in input vector(s)
*       dx  double precision vector with n elements
*     incx  storage spacing between elements of dx
*       dy  double precision vector with n elements
*     incy  storage spacing between elements of dy
*
*     --output--
*       dx  input vector dy (unchanged if n .le. 0)
*       dy  input vector dx (unchanged if n .le. 0)
*
*     interchange double precision dx and double precision dy.
*     for i = 0 to n-1, interchange  dx(lx+i*incx) and dy(ly+i*incy),
*     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
*     defined in a similar way using incy.
****references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
*                 *basic linear algebra subprograms for fortran usage*,
*                 algorithm no. 539, transactions on mathematical
*                 software, volume 5, number 3, september 1979, 308-323
****routines called  (none)
****end prologue  dswap
*
      double precision dx(1),dy(1),dtemp1,dtemp2,dtemp3
****first executable statement  dswap
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
*
*       code for unequal or nonpositive increments.
*
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp1 = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp1
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*       code for both increments equal to 1
*
*
*       clean-up loop so remaining vector length is a multiple of 3.
*
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp1 = dx(i)
        dtemp2 = dx(i+1)
        dtemp3 = dx(i+2)
        dx(i) = dy(i)
        dx(i+1) = dy(i+1)
        dx(i+2) = dy(i+2)
        dy(i) = dtemp1
        dy(i+1) = dtemp2
        dy(i+2) = dtemp3
   50 continue
      return
   60 continue
*
*     code for equal, positive, nonunit increments.
*
      ns = n*incx
        do 70 i=1,ns,incx
        dtemp1 = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp1
   70   continue
      return
      end
      integer function idamax(n,dx,incx)
****begin prologue  idamax
****date written   791001   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d1a2
****keywords  blas,double precision,linear algebra,maximum component,
*             vector
****author  lawson, c. l., (jpl)
*           hanson, r. j., (snla)
*           kincaid, d. r., (u. of texas)
*           krogh, f. t., (jpl)
****purpose  find largest component of d.p. vector
****description
*
*                b l a s  subprogram
*    description of parameters
*
*     --input--
*        n  number of elements in input vector(s)
*       dx  double precision vector with n elements
*     incx  storage spacing between elements of dx
*
*     --output--
*   idamax  smallest index (zero if n .le. 0)
*
*     find smallest index of maximum magnitude of double precision dx.
*     idamax =  first i, i = 1 to n, to minimize  abs(dx(1-incx+i*incx)
****references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
*                 *basic linear algebra subprograms for fortran usage*,
*                 algorithm no. 539, transactions on mathematical
*                 software, volume 5, number 3, september 1979, 308-323
****routines called  (none)
****end prologue  idamax
*
      double precision dx(1),dmax,xmag
****first executable statement  idamax
      idamax = 0
      if(n.le.0) return
      idamax = 1
      if(n.le.1)return
      if(incx.eq.1)goto 20
*
*        code for increments not equal to 1.
*
      dmax = dabs(dx(1))
      ns = n*incx
      ii = 1
          do 10 i = 1,ns,incx
          xmag = dabs(dx(i))
          if(xmag.le.dmax) go to 5
          idamax = ii
          dmax = xmag
    5     ii = ii + 1
   10     continue
      return
*
*        code for increments equal to 1.
*
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
          xmag = dabs(dx(i))
          if(xmag.le.dmax) go to 30
          idamax = i
          dmax = xmag
   30 continue
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
****begin prologue  dgefa
****date written   780814   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d2a1
****keywords  double precision,factor,linear algebra,linpack,matrix
****author  moler, c. b., (u. of new mexico)
****purpose  factors a double precision matrix by gaussian elimination.
****description
*
*     dgefa factors a double precision matrix by gaussian elimination.
*
*     dgefa is usually called by dgeco, but it can be called
*     directly with a saving in time if  rcond  is not needed.
*     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
*
*     on entry
*
*        a       double precision(lda, n)
*                the matrix to be factored.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*     on return
*
*        a       an upper triangular matrix and the multipliers
*                which were used to obtain it.
*                the factorization can be written  a = l*u  where
*                l  is a product of permutation and unit lower
*                triangular matrices and  u  is upper triangular.
*
*        ipvt    integer(n)
*                an integer vector of pivot indices.
*
*        info    integer
*                = 0  normal value.
*                = k  if  u(k,k) .eq. 0.0 .  this is not an error
*                     condition for this subroutine, but it does
*                     indicate that dgesl or dgedi will divide by zero
*                     if called.  use  rcond  in dgeco for a reliable
*                     indication of singularity.
*
*     linpack.  this version dated 08/14/78 .
*     cleve moler, university of new mexico, argonne national lab.
*
*     subroutines and functions
*
*     blas daxpy,dscal,idamax
****references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
*                 *linpack users  guide*, siam, 1979.
****routines called  daxpy,dscal,idamax
****end prologue  dgefa
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
*
      double precision t
      integer idamax,j,k,kp1,l,nm1
*
*     gaussian elimination with partial pivoting
*
****first executable statement  dgefa
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
*
*        find l = pivot index
*
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
*
*        zero pivot implies this column already triangularized
*
         if (a(l,k) .eq. 0.0d0) go to 40
*
*           interchange if necessary
*
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
*
*           compute multipliers
*
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
*
*           row elimination with column indexing
*
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dsifa(a,lda,n,kpvt,info)
****begin prologue  dsifa
****date written   780814   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d2b1a
****keywords  double precision,factor,linear algebra,linpack,matrix,
*             symmetric
****author  bunch, j., (ucsd)
****purpose  factors a d.p. symmetric matrix by elimination with
*            symmetric pivoting
****description
*
*     dsifa factors a double precision symmetric matrix by elimination
*     with symmetric pivoting.
*
*     to solve  a*x = b , follow dsifa by dsisl.
*     to compute  inverse(a)*c , follow dsifa by dsisl.
*     to compute  determinant(a) , follow dsifa by dsidi.
*     to compute  inertia(a) , follow dsifa by dsidi.
*     to compute  inverse(a) , follow dsifa by dsidi.
*
*     on entry
*
*        a       double precision(lda,n)
*                the symmetric matrix to be factored.
*                only the diagonal and upper triangle are used.
*
*        lda     integer
*                the leading dimension of the array  a .
*
*        n       integer
*                the order of the matrix  a .
*
*     on return
*
*        a       a block diagonal matrix and the multipliers which
*                were used to obtain it.
*                the factorization can be written  a = u*d*trans(u)
*                where  u  is a product of permutation and unit
*                upper triangular matrices, trans(u) is the
*                transpose of  u , and  d  is block diagonal
*                with 1 by 1 and 2 by 2 blocks.
*
*        kpvt    integer(n)
*                an integer vector of pivot indices.
*
*        info    integer
*                = 0  normal value.
*                = k  if the k-th pivot block is singular.  this is
*                     not an error condition for this subroutine,
*                     but it does indicate that dsisl or dsidi may
*                     divide by zero if called.
*
*     linpack.  this version dated 08/14/78 .
*     james bunch, univ. calif. san diego, argonne nat. lab.
*
*     subroutines and functions
*
*     blas daxpy,dswap,idamax
*     fortran dabs,dmax1,dsqrt
****references  dongarra j.j., bunch j.r., moler c.b., stewart g.w.,
*                 *linpack users  guide*, siam, 1979.
****routines called  daxpy,dswap,idamax
****end prologue  dsifa
      integer lda,n,kpvt(1),info
      double precision a(lda,1)
*
      double precision ak,akm1,bk,bkm1,denom,mulk,mulkm1,t
      double precision absakk,alpha,colmax,rowmax
      integer imax,imaxp1,j,jj,jmax,k,km1,km2,kstep,idamax
      logical swap
*
*     initialize
*
*     alpha is used in choosing pivot block size.
****first executable statement  dsifa
      alpha = (1.0d0 + dsqrt(17.0d0))/8.0d0
*
      info = 0
*
*     main loop on k, which goes from n to 1.
*
      k = n
   10 continue
*
*        leave the loop if k=0 or k=1.
*
*     ...exit
         if (k .eq. 0) go to 200
         if (k .gt. 1) go to 20
            kpvt(1) = 1
            if (a(1,1) .eq. 0.0d0) info = 1
*     ......exit
            go to 200
   20    continue
*
*        this section of code determines the kind of
*        elimination to be performed.  when it is completed,
*        kstep will be set to the size of the pivot block, and
*        swap will be set to .true. if an interchange is
*        required.
*
         km1 = k - 1
         absakk = dabs(a(k,k))
*
*        determine the largest off-diagonal element in
*        column k.
*
         imax = idamax(k-1,a(1,k),1)
         colmax = dabs(a(imax,k))
         if (absakk .lt. alpha*colmax) go to 30
            kstep = 1
            swap = .false.
         go to 90
   30    continue
*
*           determine the largest off-diagonal element in
*           row imax.
*
            rowmax = 0.0d0
            imaxp1 = imax + 1
            do 40 j = imaxp1, k
               rowmax = dmax1(rowmax,dabs(a(imax,j)))
   40       continue
            if (imax .eq. 1) go to 50
               jmax = idamax(imax-1,a(1,imax),1)
               rowmax = dmax1(rowmax,dabs(a(jmax,imax)))
   50       continue
            if (dabs(a(imax,imax)) .lt. alpha*rowmax) go to 60
               kstep = 1
               swap = .true.
            go to 80
   60       continue
            if (absakk .lt. alpha*colmax*(colmax/rowmax)) go to 70
               kstep = 1
               swap = .false.
            go to 80
   70       continue
               kstep = 2
               swap = imax .ne. km1
   80       continue
   90    continue
         if (dmax1(absakk,colmax) .ne. 0.0d0) go to 100
*
*           column k is zero.  set info and iterate the loop.
*
            kpvt(k) = k
            info = k
         go to 190
  100    continue
         if (kstep .eq. 2) go to 140
*
*           1 x 1 pivot block.
*
            if (.not.swap) go to 120
*
*              perform an interchange.
*
               call dswap(imax,a(1,imax),1,a(1,k),1)
               do 110 jj = imax, k
                  j = k + imax - jj
                  t = a(j,k)
                  a(j,k) = a(imax,j)
                  a(imax,j) = t
  110          continue
  120       continue
*
*           perform the elimination.
*
            do 130 jj = 1, km1
               j = k - jj
               mulk = -a(j,k)/a(k,k)
               t = mulk
               call daxpy(j,t,a(1,k),1,a(1,j),1)
               a(j,k) = mulk
  130       continue
*
*           set the pivot array.
*
            kpvt(k) = k
            if (swap) kpvt(k) = imax
         go to 190
  140    continue
*
*           2 x 2 pivot block.
*
            if (.not.swap) go to 160
*
*              perform an interchange.
*
               call dswap(imax,a(1,imax),1,a(1,k-1),1)
               do 150 jj = imax, km1
                  j = km1 + imax - jj
                  t = a(j,k-1)
                  a(j,k-1) = a(imax,j)
                  a(imax,j) = t
  150          continue
               t = a(k-1,k)
               a(k-1,k) = a(imax,k)
               a(imax,k) = t
  160       continue
*
*           perform the elimination.
*
            km2 = k - 2
            if (km2 .eq. 0) go to 180
               ak = a(k,k)/a(k-1,k)
               akm1 = a(k-1,k-1)/a(k-1,k)
               denom = 1.0d0 - ak*akm1
               do 170 jj = 1, km2
                  j = km1 - jj
                  bk = a(j,k)/a(k-1,k)
                  bkm1 = a(j,k-1)/a(k-1,k)
                  mulk = (akm1*bk - bkm1)/denom
                  mulkm1 = (ak*bkm1 - bk)/denom
                  t = mulk
                  call daxpy(j,t,a(1,k),1,a(1,j),1)
                  t = mulkm1
                  call daxpy(j,t,a(1,k-1),1,a(1,j),1)
                  a(j,k) = mulk
                  a(j,k-1) = mulkm1
  170          continue
  180       continue
*
*           set the pivot array.
*
            kpvt(k) = 1 - k
            if (swap) kpvt(k) = -imax
            kpvt(k-1) = kpvt(k)
  190    continue
         k = k - kstep
      go to 10
  200 continue
      return
      end
      double precision function dasum(n,dx,incx)
****begin prologue  dasum
****date written   791001   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d1a3a
****keywords  add,blas,double precision,linear algebra,magnitude,sum,
*             vector
****author  lawson, c. l., (jpl)
*           hanson, r. j., (snla)
*           kincaid, d. r., (u. of texas)
*           krogh, f. t., (jpl)
****purpose  sum of magnitudes of d.p. vector components
****description
*
*                b l a s  subprogram
*    description of parameters
*
*     --input--
*        n  number of elements in input vector(s)
*       dx  double precision vector with n elements
*     incx  storage spacing between elements of dx
*
*     --output--
*    dasum  double precision result (zero if n .le. 0)
*
*     returns sum of magnitudes of double precision dx.
*     dasum = sum from 0 to n-1 of dabs(dx(1+i*incx))
****references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
*                 *basic linear algebra subprograms for fortran usage*,
*                 algorithm no. 539, transactions on mathematical
*                 software, volume 5, number 3, september 1979, 308-323
****routines called  (none)
****end prologue  dasum
*
      double precision dx(1)
****first executable statement  dasum
      dasum = 0.d0
      if(n.le.0)return
      if(incx.eq.1)goto 20
*
*        code for increments not equal to 1.
*
      ns = n*incx
          do 10 i=1,ns,incx
          dasum = dasum + dabs(dx(i))
   10     continue
      return
*
*        code for increments equal to 1.
*
*
*        clean-up loop so remaining vector length is a multiple of 6.
*
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         dasum = dasum + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
         dasum = dasum + dabs(dx(i)) + dabs(dx(i+1)) + dabs(dx(i+2))
     1   + dabs(dx(i+3)) + dabs(dx(i+4)) + dabs(dx(i+5))
   50 continue
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
****begin prologue  daxpy
****date written   791001   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d1a7
****keywords  blas,double precision,linear algebra,triad,vector
****author  lawson, c. l., (jpl)
*           hanson, r. j., (snla)
*           kincaid, d. r., (u. of texas)
*           krogh, f. t., (jpl)
****purpose  d.p computation y = a*x + y
****description
*
*                b l a s  subprogram
*    description of parameters
*
*     --input--
*        n  number of elements in input vector(s)
*       da  double precision scalar multiplier
*       dx  double precision vector with n elements
*     incx  storage spacing between elements of dx
*       dy  double precision vector with n elements
*     incy  storage spacing between elements of dy
*
*     --output--
*       dy  double precision result (unchanged if n .le. 0)
*
*     overwrite double precision dy with double precision da*dx + dy.
*     for i = 0 to n-1, replace  dy(ly+i*incy) with da*dx(lx+i*incx) +
*       dy(ly+i*incy), where lx = 1 if incx .ge. 0, else lx = (-incx)*n
*       and ly is defined in a similar way using incy.
****references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
*                 *basic linear algebra subprograms for fortran usage*,
*                 algorithm no. 539, transactions on mathematical
*                 software, volume 5, number 3, september 1979, 308-323
****routines called  (none)
****end prologue  daxpy
*
      double precision dx(1),dy(1),da
****first executable statement  daxpy
      if(n.le.0.or.da.eq.0.d0) return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
*
*        code for nonequal or nonpositive increments.
*
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
*
*        code for both increments equal to 1
*
*
*        clean-up loop so remaining vector length is a multiple of 4.
*
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
*
*        code for equal, positive, nonunit increments.
*
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          dy(i) = da*dx(i) + dy(i)
   70     continue
      return
      end
      subroutine dscal(n,da,dx,incx)
****begin prologue  dscal
****date written   791001   (yymmdd)
****revision date  820801   (yymmdd)
****category no.  d1a6
****keywords  blas,linear algebra,scale,vector
****author  lawson, c. l., (jpl)
*           hanson, r. j., (snla)
*           kincaid, d. r., (u. of texas)
*           krogh, f. t., (jpl)
****purpose  d.p. vector scale x = a*x
****description
*
*                b l a s  subprogram
*    description of parameters
*
*     --input--
*        n  number of elements in input vector(s)
*       da  double precision scale factor
*       dx  double precision vector with n elements
*     incx  storage spacing between elements of dx
*
*     --output--
*       dx  double precision result (unchanged if n.le.0)
*
*     replace double precision dx by double precision da*dx.
*     for i = 0 to n-1, replace dx(1+i*incx) with  da * dx(1+i*incx)
****references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
*                 *basic linear algebra subprograms for fortran usage*,
*                 algorithm no. 539, transactions on mathematical
*                 software, volume 5, number 3, september 1979, 308-323
****routines called  (none)
****end prologue  dscal
*
      double precision da,dx(1)
****first executable statement  dscal
      if(n.le.0)return
      if(incx.eq.1)goto 20
*
*        code for increments not equal to 1.
*
      ns = n*incx
          do 10 i = 1,ns,incx
          dx(i) = da*dx(i)
   10     continue
      return
*
*        code for increments equal to 1.
*
*
*        clean-up loop so remaining vector length is a multiple of 5.
*
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
      implicit double precision (a-h, o-z)	
c     returns the dot product of double precision dx and dy.
c     ddot = sum for i = 0 to n-1 of  dx(lx+i*incx) * dy(ly+i*incy)
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
      double precision dx(1),dy(1)
      ddot = 0.d0
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c         code for unequal or nonpositive increments.
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
         ddot = ddot + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c        code for both increments equal to 1.
c        clean-up loop so remaining vector length is a multiple of 5.
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
         ddot = ddot + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
         ddot = ddot + dx(i)*dy(i) + dx(i+1)*dy(i+1) +
     $   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
      return
c         code for positive equal increments .ne.1.
   60 continue
      ns = n*incx
          do 70 i=1,ns,incx
          ddot = ddot + dx(i)*dy(i)
   70     continue
      return
      end
      double precision function dnrm2 ( n, dx, incx)
      implicit double precision (a-h, o-z)	
      integer          next
      double precision   dx(1), cutlo, cuthi, hitest, sum, xmax,zero,one
      data   zero, one /0.0d0, 1.0d0/
c     euclidean norm of the n-vector stored in dx() with storage
c     increment incx .
c           c.l.lawson, 1978 jan 08
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
      if(n .gt. 0) go to 10
         dnrm2  = zero
         go to 300
CTM10 assign 30 to next
   10  next = 1
      sum = zero
      nn = n * incx
c                                                 begin main loop
      i = 1
ctm July 2011 20    go to next,(30, 50, 70, 110)  replace assigned go to
   20    go to (30, 50, 70, 110), next
   30 if( dabs(dx(i)) .gt. cutlo) go to 85
CTM   assign 50 to next
      next = 2
      xmax = zero
c                        phase 1.  sum is zero
   50 if( dx(i) .eq. zero) go to 200
      if( dabs(dx(i)) .gt. cutlo) go to 85
c                                prepare for phase 2.
CTM   assign 70 to next
      next = 3
      go to 105
c                                prepare for phase 4.
  100 i = j
CTM   assign 110 to next
       next = 4
      sum = (sum / dx(i)) / dx(i)
  105   xmax = dabs(dx(i))
      go to 115
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
   70 if( dabs(dx(i)) .gt. cutlo ) go to 75
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
  110 if( dabs(dx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / dx(i))**2
         xmax = dabs(dx(i))
         go to 200
  115 sum = sum + (dx(i)/xmax)**2
      go to 200
c                  prepare for phase 3.
   75 sum = (sum * xmax) * xmax
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
   85 hitest = cuthi/float( n )
c                   phase 3.  sum is mid-range.  no scaling.
      do 95 j =i,nn,incx
      if(dabs(dx(j)) .ge. hitest) go to 100
   95    sum = sum + dx(j)**2
      dnrm2 = dsqrt( sum )
      go to 300
  200 continue
      i = i + incx
      if ( i .le. nn ) go to 20
c              end of main loop.
c              compute square root and adjust for scaling.
      dnrm2 = xmax * dsqrt(sum)
  300 continue
      return
      end
      subroutine dcopy(n,dx,incx,dy,incy)
c***begin prologue  dcopy
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a5
c***keywords  blas,copy,double precision,linear algebra,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  d.p. vector copy y = x
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c       dy  double precision vector with n elements
c     incy  storage spacing between elements of dy
c
c     --output--
c       dy  copy of vector dx (unchanged if n .le. 0)
c
c     copy double precision dx to double precision dy.
c     for i = 0 to n-1, copy dx(lx+i*incx) to dy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  dcopy
c
      double precision dx(1),dy(1)
c***first executable statement  dcopy
      if(n.le.0)return
      if(incx.eq.incy) if(incx-1) 5,20,60
    5 continue
c
c        code for unequal or nonpositive increments.
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop so remaining vector length is a multiple of 7.
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
c
c        code for equal, positive, nonunit increments.
c
   60 continue
      ns=n*incx
          do 70 i=1,ns,incx
          dy(i) = dx(i)
   70     continue
      return
      end
 
