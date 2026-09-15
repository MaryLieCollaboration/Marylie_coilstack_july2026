************************************************************************
* header:                 TRACK                                        *
*  Particle tracking, both symplectic and non-symplectic.              *
************************************************************************
c
c  subroutine canx
c  the following line is commented out due to reordering;
c  see the note below.
c     subroutine canx(mh,h)
c  the following line is inserted due to reordering.
      subroutine canx(mh,g)
c  subroutine to establish standard rep of transfer map for
c  use in evaluating ray traces.
c  this subroutine modified on 8/13/86 to change the order
c  in which exp(:f2:) acts.  evalsr was also modified
c  accordingly. eventually this subroutine needs careful
c  rewriting. AJD
c
      include 'impli.inc'
      include 'param.inc'
      double precision mh
      dimension mh(6,6)
      dimension dg(6,monoms)
      dimension h(monoms),g(monoms),gtemp(monoms)
      dimension dgtmp(monoms),crstrm(monoms),f(monoms)
      dimension temp(monoms),dum(monoms),tmp1(monoms),tmp2(monoms)
      include 'deriv.inc'
c  colme = gives column numbers of the non-zero matrix elements, listed
c          by row. Entry row,0 has the total number of non-zero elements
c          in that row.
      integer colme(6,0:6)
c
c  initialize arrays
c
      do 1000 n=1,monoms
c  the following two lines have been commented out
c  due to reordering.
c       g(n)=0.d0
c       gtemp(n)=0.d0
        dgtmp(n)=0.d0
        crstrm(n)=0.d0
        f(n)=0.d0
        temp(n)=0.d0
        tmp1(n)=0.d0
        tmp2(n)=0.d0
        dum(n)=0.d0
        do 1 m=1,6
          dg(m,n)=0.d0
          df(m,n)=0.d0
    1   continue
        do 100 l=1,3
          do 10 m=1,3
            rjac(l,m,n)=0.d0
   10     continue
  100   continue
 1000 continue
c  the following section has been commented out
c  due to reordering.
c
c  transfor nonlinear arrays by linear piece
c
c     do 200 k=3,4
c        do 2 l=1,209
c           gtemp(l)=0.d0
c   2    continue
c        call xform(h,k,mh,k-3,gtemp)
c        do 20 l=1,209
c           g(l)=g(l)+gtemp(l)
c  20    continue
c 200 continue
c
c  determine derivs of array g
c
      do 300 i=1,6
         do 3 j=1,monoms
            dgtmp(j)=0.d0
    3    continue
         ivalue=3
         call svpbkt(g,ivalue,i,1, dgtmp)
         do 30 j=1,monoms
            dg(i,j)=dgtmp(j)
   30    continue
  300 continue
c
c  call prod to compute the crossterm
c
      call prod(dg,crstrm)
c
c  sum all the contributions to the 3rd and 4th orders
c  of the generating fxn
c
      do 4 n=1,monoms
         f(n)=f(n)-g(n)+crstrm(n)
    4 continue
c
c  add in 2nd order piece
c
      f(8)=f(8)+1.d0
      f(19)=f(19)+1.d0
      f(26)=f(26)+1.d0
c
c  now compute the derivs of the generating fxn f (which define
c  the standard representation of the canonical transformation)
c
      do 5000 i=1,5,2
         do 500 k=2,4
            do 5 l=1,monoms
               tmp1(l)=0.d0
               tmp2(l)=0.d0
    5       continue
            ivalue=i+1
            call svpbkt(f,k,ivalue,1, tmp1)
            call svpbkt(f,k,i,1, tmp2)
            do 50 l=1,monoms
               df(i,l)=df(i,l)+tmp1(l)
               df(ivalue,l)=df(ivalue,l)-tmp2(l)
   50       continue
  500    continue
 5000 continue
c
c  now compute the jacobian of the momentum part of the
c  standard rep (which will be used in sub. newton to determine
c  the new momentum)
c
      ivalue=0
      do 6000 i=1,5,2
         do 6 n=1,monoms
            dum(n)=df(i,n)
    6    continue
         ivalue=ivalue+1
         jvalue=0
         do 600 j=1,5,2
            jvalue=jvalue+1
            do 60 k=2,3
               do 58 l=1,monoms
                  temp(l)=0.d0
   58          continue
               call svpbkt(dum,k,j,1, temp)
               do 59 l=1,monoms
                  rjac(ivalue,jvalue,l)=rjac(ivalue,jvalue,l)+temp(l)
   59          continue
   60       continue
  600    continue
 6000 continue
      return
      end
c
********************************************************************
c
      subroutine eval(mh,zi,zf)
c
c     evaluates the image zf of initial
c     data zi under the lie transformation
c     with linear part having matrix representation
c     mh and nonlinear part having generators
c     with coefficients stored in the array h
c     Written by D. Douglas, ca 1982
c
      include 'impli.inc'
      include 'param.inc'
      double precision mh
      dimension zi(6),zf(6),tempz(6)
      dimension mh(6,6)
      dimension vect(monoms)
      include 'pbkh.inc'
      include 'files.inc'
      include 'vblist.inc'
      include 'lims.inc'
      include 'trord.inc'
c
c     initialise arrays
c
      do 10 i=1,6
      tempz(i)=0.0d0
      zf(i)=0.0d0
   10 continue
c
c     generate tempz=mh*zi (linear contribution)
c
      do 30 i=1,6
      do 20 k=1,6
      tempz(i)=tempz(i)+mh(i,k)*zi(k)
   20 continue
   30 continue
c
c     loop to compute final values zf(i), including nonlinearities
c
      do 40 i=1,monoms
      vect(i)=0.0d0
   40 continue
      if(trord.ge.2) then
      do 50 i=7,27
      vect(i)= vect(i)+tempz(vblist(1,i))*tempz(vblist(2,i))
   50 continue
      endif
      if(trord.ge.3) then
      do 55 i=28,83
      vect(i)=vect(i)+tempz(vblist(1,i))*tempz(vblist(2,i))
     $        *tempz(vblist(3,i))
   55 continue
      endif
      if(trord.ge.4) then
      do 57 i=84,209
      vect(i)=vect(i)+tempz(vblist(1,i))*tempz(vblist(2,i))
     #*tempz(vblist(3,i))*tempz(vblist(4,i))
   57 continue
      endif
      if(ibrief.eq.2) write(jof,56)
   56 format(/1h ,'lengthy ray trace data')
      do 1000 i=1,6
      if(ibrief.ne.2)goto 59
      write(jof,58)i
   58 format(/1h ,'i =',i2)
   59 j=i+6
c
c     initialise values of polynomials giving nonlinearities
c
      vph3=0.0d0
      vph4=0.0d0
      vphh=0.0d0
      vphhh=0.0d0
      vph34=0.0d0
c
c     evaluate the polynomials
c
      do 60 k=7,27
      vph3=vph3+pbh(k,i)*vect(k)
   60 continue
      do 70 k=28,83
      vph4=vph4+pbh(k,i)*vect(k)
      vphh=vphh+((pbh(k,j)*vect(k))/2.d0)
   70 continue
      do 80 k=84,209
      vphhh=vphhh+pbh(k,j)*vect(k)/6.d0
      vph34=vph34+pbh(k,i)*vect(k)
   80 continue
c
c     evaluate final conditions zf(i)
c
      zf(i)=tempz(i)+vph3+vph4+vphh+vphhh+vph34
      if(ibrief.ne.2)goto 1000
      write(jof,101) tempz(i)
  101 format(1h ,'rmat = ',d21.15)
      write(jof,102) vph3,vph4
  102 format(1h ,'rf3 = ',d21.15/1x,'rf4 = ',d21.15)
      write(jof,103) vphh,vphhh,vph34
  103 format(1h ,'rf3^2/2 = ',d21.15/1x,'rf3^3/6 = ',d21.15/
     *1x,'rf3f4 = ',d21.15)
 1000 continue
      return
      end
c
***********************************************************************
c
      subroutine evalsr(mh,zi,zf,df,rdf,rrjac)
c  subroutine to evaluate the image zf of initial data zi under the
c  lie transformation with linear part mh and nonlinear part
c  represented by a standard representation stored in common/deriv/df
c  this subroutine modified on 8/13/86 to change the
c  order in which exp(:f2:) acts.
c  canx was also modified accordingly.  this
c  whole subroutine needs rewriting badly AJD
c
      include 'impli.inc'
      include 'param.inc'
      double precision mh
      dimension mh(6,6)
      dimension df(6,83),rdf(3,84),rrjac(3,3,28)
      dimension zi(6),zf(6),ztmp1(6),ztmp2(6),ztmp3(6)
      dimension vect(monoms)
c zlm added for ordering modification 8/13/86
      dimension zlm(6)
      include 'files.inc'
      include 'ind.inc'
      include 'len.inc'
c
c  initialize arrays
c
      do 1 i=1,6
         zf(i)=0.d0
         ztmp1(i)=0.d0
         ztmp2(i)=0.d0
         ztmp3(i)=0.d0
c also initialize zlm
         zlm(i)=0.d0
    1 continue
c  the following section is added due to reordering:
c
c  compute effect of linear part of the map on zi
c
      do 40 i=1,6
         do 4 j=1,6
            zlm(i)=zlm(i) + mh(i,j)*zi(j)
    4    continue
   40 continue
c
c  call newton search routine to find value of new momentum
c
c  the following statement is commented out due to reordering
c     call newt(zi,ztmp1,rdf,rrjac)
c  the folliwing statement is added due to reordering
      call newt(zlm,ztmp1,rdf,rrjac)
c
c  compute vector containing values of basis monomials
c
      do 2 i=1,6
          vect(i) = ztmp1(i)
    2 continue
c
c  (only terms through order imaxi-1 are required)
c
      do 20 i = 7,len(imaxi-1)
      vect(i) = vect(index1(i))*vect(index2(i))
 20   continue
c  compute values of the new coords and old momenta using the standard
c  rep of the transfer map which is stored in df
c
      do 3000 i=1,5,2
         ivalue=i+1
         do 300 j=1,len(imaxi-1)
c
c  new coordinate values
            ztmp2(i)=ztmp2(i)+df(ivalue,j)*vect(j)
c      ztmp2(i)= ddot(len(imaxi-1),vect,1,df(ivalue,1),6)
c  old momentum values (done as a check)
            ztmp3(ivalue)=ztmp3(ivalue)+df(i,j)*vect(j)
  300 continue
c      ztmp3(ivalue)=ddot(len(imaxi-1),vect,1,df(i,1),6)
c
c  transfer new momentum values to ztmp2 (these were returned from
c  newt in ztmp1)
c
         ztmp2(ivalue)=ztmp2(ivalue)+ztmp1(ivalue)
 3000 continue
c
c  at this point, the nonlinearities are fixed.  the image of
c  the initial data under the nonlinear portion of the
c  transformation is stored in ztmp2.  before applying the linear
c  part of the map to this, check to see if the old momentum values
c  which were read in match those computed using the standard
c  representation (done immediately above.)
c
c  nonlinearities check
c
      delpx=zi(2)-ztmp3(2)
      delpy=zi(4)-ztmp3(4)
      delpz=zi(6)-ztmp3(6)
      if(ibrief.ne.2)goto 302
      write(jof,9250)
 9250 format(1h ,' **momentum deviations** ')
      write(jof,9300)delpx,delpy,delpz
 9300 format(1h ,3(d22.15,2x))
  302 continue
c  the following statements are commented out due to reordering
c
c  compute effect of linear part of the map on ztmp2
c
c     do 40 i=1,6
c        do 4 j=1,6
c           zf(i)=zf(i) + mh(i,j)*ztmp2(j)
c   4    continue
c  40 continue
c  add the following lines due to reordering:
      do 137 i=1,6
  137 zf(i)=ztmp2(i)
      return
      end
c
***********************************************************************
c
      subroutine invrs(a,b,d)
c  subroutine to compute the inverse b and determinant d
c  of a 3x3 matrix a.  used by newton search subroutine
c  when searching for the new momentum.
c     Written by D. Douglas, ca 1983
      include 'impli.inc'
      include 'files.inc'
      dimension a(3,3),b(3,3),c(3,3)
      do 10 i=1,3
         do 1 j=1,3
            b(i,j)=0.d0
            c(i,j)=0.d0
    1    continue
   10 continue
c
c  compute determinant of a
c
      d=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)
     # -a(1,3)*a(2,2)*a(3,1)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)
c     write(jof,9100)d
c9100 format(1h ,'determinant of a is ',d21.15)
      abd=dabs(d)
      if(abd.lt.1.d-30)write(6,9000)
 9000 format(1h ,'determinant underflow in subroutine invrs')
      if(abd.lt.1.d-30)return
c
c  compute b=inverse of a
c
      b(1,1)=(a(2,2)*a(3,3)-a(3,2)*a(2,3))/d
      b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))/d
      b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))/d
      b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))/d
      b(1,3)=(a(1,2)*a(2,3)-a(2,2)*a(1,3))/d
      b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))/d
      b(2,3)=(a(2,1)*a(1,3)-a(1,1)*a(2,3))/d
      b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))/d
      b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))/d
c
c  check if this is a reasonable inverse
c
c     do 200 i=1,3
c        do 20 j=1,3
c           do 2 k=1,3
c              c(i,j)=c(i,j) + a(i,k)*b(k,j)
c   2       continue
c  20    continue
c 200 continue
c     write(6,9150)
c9150 format(1h ,'a * ainverse = '/)
c     do 3 i=1,3
c        write(6,9200)(c(i,j),j=1,3)
c9200    format(1h ,3(d21.15,2x))
c   3 continue
      return
      end
c
***********************************************************************
c
      subroutine newt(zi,zo,rdf,rrjac)
c  subroutine to invert map p(old)=df(q(old),p(new))/dq(old)
c  for p(new) using a newton's search procedure
c
c  this is an exemple of a newton search subroutine that will
c  work with polynomials of order up to 10
c  if the arrays are properly dimensioned.
c     Written by D. Douglas, ca 1983 and substantially modified
c  by F. Neri
c
      include 'impli.inc'
      include 'param.inc'
      dimension rdf(3,84),rrjac(3,3,28)
      dimension ztmp(6),zi(6),zo(6)
      dimension rlin(3,3),rmat(3,3),ri(3,3),rinv(3,3)
      dimension pi(3),crtrm(3),pimag(3),p(3),delp(3),cp(3)
      dimension pjac(3,3,20),pdf(3,20)
      dimension qvec(20),pvec(20)
      include 'files.inc'
      include 'ind.inc'
      include 'len.inc'
      include 'len3.inc'
      include 'ind3.inc'
      data ri/1.,3*0.,1.,3*0.,1./
      include 'talk.inc'
c
c  initialize arrays
c
      l31=len3(imaxi-1)
      l32=len3(imaxi-2)
      ivalue=0
      do 1 i=1,3
         ivalue=ivalue+2
         pi(i)=zi(ivalue)
    1 continue
      do 2 i=1,6
         ztmp(i)=zi(i)
    2 continue
c
c generate monomials in the q's
c
      qvec(1) = 1.0d0
      ivalue=-1
      do 10 i=1,3
      ivalue=ivalue+2
      qvec(i+1)=ztmp(ivalue)
 10   continue
      do 11 i=5,len3(imaxi-1)+1
      qvec(i)=qvec(ind31(i))*qvec(ind32(i))
 11   continue
c
c  generate coefficient of polynomials in the momentum
c  variables
c
      do 1001 i1=1,3
      do 1001 i2=1,3
      pjac(i1,i2,1)=0.d0
 1001 continue
c
      do 1002 n=1,4
      do 1003 i=1,3
      pdf(i,n) = 0.d0
 1003 continue
 1002 continue
c
c generate coefficients of reduced df (pdf),function only of p:
c***********************
      do 190 i = 1,3
      do 190 ip = 5,len3(imaxi-1)+1
      pdf(i,ip) = rdf(i,ip)
  190  continue
c***** previous loops replaced by: *************
c     call dcopy(3*(len3(imaxi-1)-3),rdf(1,5),1,pdf(1,5),1)
c**********************
      n0 = len3(imaxi-1)+1
      iq2=1
      do 9 ipoq=1,imaxi-2
      iq1=iq2+1
      iq2=len3(ipoq)+1
      if(imaxi.eq.ipoq+1) l = 1
      if(imaxi.gt.ipoq+1) l = len3(imaxi-1-ipoq)+1
      do 99 iq=iq1,iq2
      qval = qvec(iq)
c********************
      do 191 i = 1,3
      do 191 ip = 1,l
      pdf(i,ip) = pdf(i,ip) + rdf(i,n0+ip)*qval
  191  continue
c**** previous loops replaced by: **********
c     call daxpy(3*l,qval,rdf(1,n0+1),1,pdf(1,1),1)
c********************
          n0 = n0 + l
 99      continue
 9      continue
c********************
      do 999 i = 1,3
      do 999 ip = 1,l31-l32
      pdf(i,1) = pdf(i,1) + rdf(i,n0+ip)*qvec(l32+1+ip)
  999  continue
c**** previous loops replaced by: **********
c     pdf(1,1)=pdf(1,1)+ddot(l31-l32,qvec(l32+2),1,rdf(1,n0+1),3)
c     pdf(2,1)=pdf(2,1)+ddot(l31-l32,qvec(l32+2),1,rdf(2,n0+1),3)
c     pdf(3,1)=pdf(3,1)+ddot(l31-l32,qvec(l32+2),1,rdf(3,n0+1),3)
c********************
c
c   genarate coefficients of reduce  jacobian (rrjac):
c********************
      do 1190 i1 = 1,3
      do 1190 i2 = 1,3
      do 1190 ip = 2,len3(imaxi-2)+1
      pjac(i1,i2,ip) = rrjac(i1,i2,ip)
 1190 continue
c**** previous loops replaced by: **********
c     call dcopy(9*len3(imaxi-2),rrjac(1,1,2),1,pjac(1,1,2),1)
c********************
      n0 =  len3(imaxi-2)+1
      iq2=1
      do 19 ipoq = 1,imaxi-2
      iq1 = iq2 + 1
      iq2 = len3(ipoq)+1
      if(imaxi.eq.ipoq+2) l = 1
      if(imaxi.gt.ipoq+2) l = len3(imaxi-2-ipoq)+1
      do 199 iq = iq1,iq2
      qval = qvec(iq)
c*********************
      do 1999 i1 = 1,3
      do 1999 i2 = 1,3
      do 1999 ip = 1,l
      pjac(i1,i2,ip) = pjac(i1,i2,ip) + rrjac(i1,i2,n0+ip)*qval
 1999 continue
c**** previous loops replaced by: **********
c     call daxpy(9*l,qval,rrjac(1,1,n0+1),1,pjac(1,1,1),1)
c********************
          n0 = n0 + l
 199  continue
 19   continue
c
c
c  loop to apply contraction mapping
c
      do 9999 index=1,12
c
         do 30 i=1,3
            crtrm(i)=0.d0
            pimag(i)=0.d0
            p(i)=0.d0
            delp(i)=0.d0
            cp(i)=0.d0
            do 3 j=1,3
               rlin(i,j)=0.d0
    3       continue
   30    continue
c   compute monomials in the momenta
c
      pvec(1) = 1.d0
c  first order monomials are equal to the p's;
      ivalue=0
      do 4 i=1,3
      ivalue=ivalue+2
      pvec(i+1)=ztmp(ivalue)
 4    continue
c
c  compute the remaining monomials as products of
c  previous ones:
c  need only term up to order imaxi-1
      do 40 i=5,len3(imaxi-1)+1
      pvec(i)=pvec(ind31(i))*pvec(ind32(i))
 40   continue
c
c     compute jacobian
c
         do 500 i=1,3
            do 50 j=1,3
c only terms of order up to imaxi-2 are present in the jacobian,
c the jacobian is produced by taking the second derivatives of a
c imaxi order polynomial.
c********************
           do 5 n =1,len3(imaxi-2)+1
                rlin(i,j) = rlin(i,j) + pjac(i,j,n)*pvec(n)
    5       continue
c**** previous loop replaced by: **********
c        rlin(i,j)=ddot(l32+1,pvec,1,pjac(i,j,1),9)
c********************
c
c
               rmat(i,j)=ri(i,j)-rlin(i,j)
   50       continue
  500    continue
         call invrs(rmat,rinv,det)
         adet=dabs(det)
         if(adet.lt.1.d-30)write(6,501)
  501    format(1h ,'determinant underflow in newt')
         if(adet.lt.1.d-30)return
         ivalue=0
         do 6 i=1,3
            ivalue=ivalue+2
            p(i)=ztmp(ivalue)
    6    continue
         do 600 i=1,3
            pimag(i)=pimag(i)+pi(i)
c only term of order 2 to imaxi-1 are present
c  so the sum only goes up to len(imaxi-1) .   not
c********************
             do 60 n=1,l31+1
                   pimag(i) = pimag(i) - pdf(i,n)*pvec(n)
   60     continue
c**** previous loop repplaced by: **********
c       pimag(i)=pimag(i)-ddot(l31+1,pvec,1,pdf(i,1),3)
c********************
c
c
            delp(i)=p(i)-pimag(i)
  600    continue
         do 70 i=1,3
            do 7 j=1,3
               crtrm(i)=crtrm(i)+rinv(i,j)*delp(j)
    7       continue
   70    continue
c         square=crtrm(1)**2 + crtrm(2)**2 + crtrm(3)**2
c         root=dsqrt(square)
         root=abs(crtrm(1)) + abs(crtrm(2)) + abs(crtrm(3))
         if(ibrief.ne.2)goto 705
         write(jof,71)index,(crtrm(i),i=1,3)
   71    format(1h ,'corrections at iteration ',i2/
     #   1h ,'crtrm(1) = ',d21.15/1h ,'crtrm(2) = ',d21.15/
     #   1h ,'crtrm(3) = ',d21.15)
c         if (1.+root.eq.1.) write(jof,72)index
         if (root .le. 1.d-12) write(jof,72)index
   72    format(1h ,'iteration has converged; i= ',i1)
  705    continue
         do 8 i=1,3
            cp(i)=p(i)-crtrm(i)
    8    continue
         ivalue=0
         do 80 i=1,3
            ivalue=ivalue+2
            ztmp(ivalue)=cp(i)
   80    continue
         do 800 i=1,6
            zo(i)=ztmp(i)
  800    continue
c         if(1.+root.eq.1.) goto 1234
         if(root .le. 1.d-12) goto 1234
 9999 continue
      write(6,906) root
  906 format(1h,' Search did not converge; root= ',e12.5)
c      if (root.gt..1) jwarn=1
      jwarn=1
 1234 continue
      return
      end
c
***********************************************************************
c
      subroutine prod(dg,crstrm)
c  subroutine to compute crossterm in generating fxn of standard
c  rep of transfer map
c     Written by D. Douglas, ca 1982, and modified by Liam Healy
c
      include 'impli.inc'
      include 'param.inc'
      dimension dg(6,*),crstrm(*),l(6)
      include 'expon.inc'
      include 'lims.inc'
c
      do 1 m=1,top(4)
         crstrm(m)=0.d0
    1 continue
      do 2000 i=1,5,2
         ivalue=i+1
         do 200 j=bottom(2),top(2)
            if(dg(i,j).eq.0.d0)goto 200
            do 20 k=bottom(2),top(2)
               if(dg(ivalue,k).eq.0.d0)goto 20
               do 2 m=1,6
                  l(m)=expon(m,j)+expon(m,k)
    2          continue
               n=ndex(l)
               crstrm(n)=crstrm(n) - dg(i,j)*dg(ivalue,k)/2.d0
   20       continue
  200    continue
 2000 continue
      return
      end
c
***********************************************************************
c
      subroutine rearr
c Written by F. Neri, ca 1984
      include 'impli.inc'
      include 'param.inc'
      dimension j(6)
      include 'ja3.inc'
      include 'ind.inc'
      include 'deriv.inc'
      include 'len.inc'
      include 'len3.inc'
      include 'ind3.inc'
c
      do 22 i=1,6
      j(i)=0
 22   continue
      do 1 ip = 2,len3(imaxi-1)+1
      ivalue = 0
      do 11 i = 1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,ip)
 11   continue
      n = ndex(j)
      if(ip.gt.(len3(imaxi-2)+1)) goto 112
      do 111 i1 = 1,3
      do 111 i2 = 1,3
      rrjac(i1,i2,ip) = rjac(i1,i2,n)
 111  continue
 112   continue
      ivalue = -1
      do 10 i = 1,3
      ivalue = ivalue + 2
      rdf(i,ip) = df(ivalue,n)
 10   continue
 1    continue
      n1 = len3(imaxi-1)+1
      n2 = len3(imaxi-2)+1
      iq2 = 1
      do 3 ipoq = 1,imaxi-1
      iq1 = iq2+1
      iq2 = len3(ipoq)+1
      if(imaxi-2.gt.ipoq) then
      l2 = len3(imaxi-2-ipoq)+1
      l1 = len3(imaxi-1-ipoq)+1
      else if(imaxi-2.eq.ipoq) then
      l2 = 1
      l1 = len3(imaxi-1-ipoq)+1
      else
      l1 = 1
      end if
      do 33 iq = iq1,iq2
      if(imaxi-2.ge.ipoq) then
      do 333 ip = 1,l2
      ivalue = -1
      do 3333 i=1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,iq)
      j(ivalue+1) = ja3(i,ip)
 3333 continue
      n = ndex(j)
      do 4 i1=1,3
      do 4 i2 = 1,3
      rrjac(i1,i2,n2+ip) = rjac(i1,i2,n)
 4    continue
 333  continue
      n2 = n2 + l2
      endif
      do 334 ip = 1,l1
      ivalue = -1
      do 3334 i = 1,3
      ivalue = ivalue + 2
      j(ivalue) = ja3(i,iq)
      j(ivalue+1) = ja3(i,ip)
 3334 continue
      n = ndex(j)
      ivalue = -1
      do 44 i = 1,3
      ivalue = ivalue+2
      rdf(i,n1+ip)  = df(ivalue,n)
 44   continue
 334  continue
      n1 =  n1 + l1
 33   continue
 3    continue
      return
      end
c
***********************************************************************
c
      subroutine trace(icfile,norder,ntrace,nwrite,th,tmh)
c     Written by D. Douglas, ca 1982, and modified since by
c     nearly everyone at Maryland.  Could still stand improvement.
c
      include 'impli.inc'
      include 'param.inc'
      include 'deriv.inc'
      include 'rays.inc'
      include 'files.inc'
      include 'talk.inc'
      include 'infin.inc'
      include 'trord.inc'
c calling arrays
      dimension th(monoms),tmh(6,6)
c
c  read initial conditions, if requested:
      if(icfile.ne.0)call raysin(icfile)
c
c  examine various cases for norder
c  norder=-1
      if(norder.eq.-1.and.icfile.ne.0) return
      if(norder.eq.-1.and.icfile.eq.0) then
        write(6,*) 'icfile and norder have inconsistent values'
        write(12,*) 'icfile and norder have inconsistent values'
        call myexit
      endif
c
c  norder=0
c  write final conditions in full precision
      if(norder.eq.0.and.jfcf.lt.0)then
        do 120 i=1,nrays
        do 110 j=1,6
        write(-jfcf,*)zblock(i,j)
  110   continue
  120   continue
        write(jof,130) -jfcf
  130   format(1x,'final conditions written in full precision on',
     #  ' file ',i4)
      endif
c
c  write final conditions in standard format
      if(norder.eq.0.and.jfcf.gt.0)then
        do 135,i=1,nrays
        write(jfcf,136)(zblock(i,j),j=1,6)
  136   format(6(1x,1pe12.5))
  135   continue
      endif
      if(norder.eq.0)return
c
c  cases where norder > 0
c
c  call preliminary routines depending on type of ray trace:
      if(norder.eq.5)then
        call canx(tmh,th)
        call rearr
      else
         trord = norder
         call brkts(th)
      endif
c
c  setup before entering main loop:
      ktrace=ntrace
      if(ktrace.eq.0)ktrace=1
c
c  the two statements below seem to cause a problem
c  and have been commented out
c      kwrite=nwrite
c      if(kwrite.eq.0)kwrite=1
c
c main loop; do 'ktrace' ray traces (of the 'nrays' initial rays)
      do 300 idum=1,ktrace
      do 200 k=1,nrays
      if (nlost.ge.nrays) then
        write (jof,*) 'all particles lost'
        return
      endif
c check to see if particle has already been lost
         if (istat(k).ne.0) goto 200
c
         do 150 l=1,6
  150    zi(l)=zblock(k,l)
c
c see if particle should be 'lost' for having coordinates too large
         jwarn=0
         if ( abs(zi(1)) .gt. xinf ) then
             write(6,*) ' particle lost by xinf overflow'
          jwarn=1
          goto 151
         endif
         if( abs(zi(3)) .gt. yinf ) then
             write(6,*) ' particle lost by yinf overflow'
          jwarn=1
          goto 151
         endif
         if( abs(zi(5)) .gt. tinf ) then
             write(6,*) ' particle lost by tinf overflow'
          jwarn=1
          goto 151
         endif
c
c procedure for a nonsymplectic ray trace
c
      if(norder.le.4) then
      call eval(tmh,zi,zf)
      endif
c
c procedure for a symplectic ray trace
c
         if(norder.gt.4) then
         jwarn=0
         call evalsr(tmh,zi,zf,df,rdf,rrjac)
         endif
c
c procedure if particle was 'lost' either by being outside
c the infinity bounds or by the symplectic ray tracer.
c mark so 'lost' particles by a a distinctive negative number.
c
  151  continue
             if(jwarn .ne. 0) then
c             write(6,*) 'here I am',jwarn
             istat(k)=-idum
             nlost=nlost+1
             ihist(nlost,1)=-idum
             ihist(nlost,2)=k
             endif
c
c procedure for particles that were not lost
c
c 1511  continue
c
c  the line below has been replaced by the following two lines
c         if(mod(idum,kwrite).ne.0)goto 176
          if(nwrite.eq.0)goto 176
          if(mod(idum,nwrite).ne.0)goto 176
c  end of modifications
c  write out final conditions
c  full precision case
         if(jfcf.lt.0) then
           do 152 l=1,6
           write(-jfcf,*) zf(l)
  152      continue
         endif
c  standard format case
         if(jfcf.gt.0) then
           write(jfcf,171)(zf(l),l=1,6)
  171      format(6(1x,1pe12.5))
         endif
         if(ibrief.eq.0)goto 176
c  print lengthy info:
         write(jof,155)
  155    format(/1h ,'initial conditions are: (dimensionless form)')
         write(jof,*)zi(1),zi(2)
         write(jof,*)zi(3),zi(4)
         write(jof,*)zi(5),zi(6)
         write(jof,160)
  160    format(/1h ,'final conditions are: (dimensionless form)')
         write(jof,*)zf(1),zf(2)
         write(jof,*)zf(3),zf(4)
         write(jof,*)zf(5),zf(6)
c  (end of lengthy info)
c
  176    if(ntrace.eq.0)goto 200
         do 180 mm=1,6
  180    zblock(k,mm)=zf(mm)
c
  200 continue
  300 continue
      return
      end
c
c end of file
