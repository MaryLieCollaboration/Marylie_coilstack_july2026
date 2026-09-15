ctm2014  internal extra text virtual files (used for fitting etc.)
c    mltext(k) = kth extra line of text
c     lunex(k) = virtual file id (logical unit number) of kth record
c        kurid = current id (lun)
c        kurec = current record relative to top of lun
c        kzero = base offset of current file (lun)
c        jrec = actual index of current record
c-----------------------------------------------
      parameter (maxrec=512)
      character*80  mltext
      common/indata/kurid,kurec,kzero,jrec,lunex(maxrec),mltext(maxrec)
ctm   save /indata/
