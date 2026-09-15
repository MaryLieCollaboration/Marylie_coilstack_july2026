      subroutine user10(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      include 'files.inc'
      include 'usrdat.inc'
      include 'taylor.inc'
      include 'mldex.inc'
      include 'cxdata.inc'
      dimension p(6), fa(monoms), fm(6,6)
c
      lfile = nint(p(1))
      nread = nint(p(2))
      if(nread.ne.1) go to 20
      write(6,*) 'User10 called to load cxdata block from lun',lfile
      call xdload(lfile,nxd)
c
c         dump cxdata common
c
  20  continue
      write(6,22) nxd
  22  format('USR10:',i4,' records in #xdata buffers:')
      if(nxd.le.0) return
      do 30 j = 1,nxd
         write(6,27) j,lunex(j),mltext(j)
  27     format(i3,': id=',i3,2x,a)
  30  continue
      return
      end
c--------------------------------------------
      subroutine xdload(lfile,nxd)
c
c    read external file into internal text common block.
c    C.T.Mottershead Jan 2014
c----------------------------------------------------------
      character*128 card
      include 'cxdata.inc'
      write(6,*) 'XLOAD call to read file ',lfile
      jrec = 0
      nxd = 0
      jj = 0
  10  read(lfile,117,end=190) card
 117  format(a)
      jj = jj + 1
      write(6,119) jj, card
 119  format('card',i3,a)
      kk = index(card,'ixd>:')
      if(kk.gt.0) then
         read(card(kk+6:),*) kurid
         write(6,*) kk,'=kk',kurid,'=kurid'
         go to 10
      endif
      jrec = jrec + 1
      if(jrec.gt.maxrec) go to 190
      if(nxd.lt.jrec) nxd = jrec
      mltext(jrec) = card
      lunex(jrec) = kurid
      go to 10
      return
 190  continue
      write(6,*) nxd,' lines read into mltext common.'
      return
      end
c------------------------------------------
      subroutine user11(p,fa,fm)
c  P. L. Walstrom 6/24/2011
c  This user routine computes u=R16**2-R12*R56
c  no parameters in call, uses R matrix only
      include 'impli.inc'
      include 'param.inc'
      include 'map.inc'
      include 'usrdat.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Walstrom user11 called'
c
c  map.inc listing:
c c total transfer map
c      real*8     th        ,tmh     ,reftraj   , arclen
c      common/map/th(monoms),tmh(6,6),reftraj(6), arclen
c
c  usrdat.inc listing:
c!cryne 7 April 2011 added integer and real declarations
c      integer nuvar,kusr
c      real*8 ucalc,wa
c      common/usrdat/nuvar,kusr,ucalc(250),wa(15)
c
      u=tmh(1,6)**2-tmh(1,2)*tmh(5,6)
      ucalc(111)=u
      write(6,*) 'user11 ucalc(111) R16**2-R12*R56 value=',u
      return
      end
c------------------------------------------
      subroutine user12(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      lun = nint(p(1))
      call mldump(lun)
      write(6,*) 'MLDUMP called from user12.'
      return
      end
c------------------------------------------
      subroutine user13(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user13 called'
      return
      end
c------------------------------------------
      subroutine user14(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user14 called'
      return
      end
c------------------------------------------
      subroutine user15(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user15 called'
      return
      end
c------------------------------------------
      subroutine user16(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user16 called'
      return
      end
c------------------------------------------
      subroutine user17(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user17 called'
      return
      end
c------------------------------------------
      subroutine user18(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user18 called'
      return
      end
c------------------------------------------
      subroutine user19(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user19 called'
      return
      end
c------------------------------------------
      subroutine user20(p,fa,fm)
      include 'impli.inc'
      include 'param.inc'
      dimension p(6), fa(monoms), fm(6,6)
      write(6,*) 'Dummy user20 called'
      return
      end
c------------------------------------------
