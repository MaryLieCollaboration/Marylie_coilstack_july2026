      subroutine mldump(lun)
c-----------------------------------------------------------------------
c  Write Marylie commons to next file on unit lun.
c  User routine CTM 14 May 2012, from original dump by Petra Schuett 26 Oct 1987
c-----------------------------------------------------------------------
      include 'impli.inc'
      include 'param.inc'
      include 'bmline.inc'
      include 'elmnts.inc'
      include 'items.inc'
      include 'coment.inc'
      include 'mldex.inc'
      include 'codes.inc'
      include 'parm.inc'
      include 'files.inc'
c-----------------------------------------------------------------------
c comments
      if(np.gt.0) then
        write(lun,500) ling(1)
        write(lun,510) (mline(i),i=1,np)
  510   format(1x,a)
      endif
c--------------------
c  beam
      write(lun,500) ling(2)
  500 format(a8)
      write(lun,*) brho
      write(lun,*) gamm1
      write(lun,*) achg
      write(lun,*) sl
c--------------------
c  menu
      write(lun,500) ling(3)
      do 10 k=1,na
       write(lun,520) lmnlbl(k),ltc(nt1(k),nt2(k))
       maxp=nrp(nt1(k),nt2(k))
  520 format(1x,a8,1x,a8)
       if(maxp.gt.0) then
        write(lun,522)(pmenu(i+mpp(k)),i=1,maxp)
       endif
  522 format(6(1pg22.15))
   10 continue
c--------------------
c  lines,lumps,loops
      if(nb.ne.0) then
       do 40 ii=2,4
       write(lun,500) ling(ii+2)
        do 40 k=1,nb
        if(ityp(k).eq.ii) then
          write(lun,530) ilbl(k)
  530     format(1x,a8)
          write(lun,532)(irep(l,k),icon(l,k),l=1,ilen(k))
  532     format((1x,5(i5,'*',a8),1x,:'&'))
        endif
   40  continue
   50  continue
      endif
c--------------------
c labor
      if(noble.ne.0) then
      write(lun,500) ling(7)
      do 100 j=1,noble
         write(lun,540) num(j),latt(j)
  540    format(1x,i4,'*',a8)
  100 continue
      endif
      return
      end
