c
c  keyword data for aimset
c
      parameter (nkeys=56)
      common/keys/key(nkeys)
      character*2 key
      common/kasks/kask(nkeys),maxj(13)
      data maxj/2*6,2*monoms,3*6,83,monoms,250,5,3,4/
      data key/'sm','bm','sf','bf','r(','s(','z(','tm','f(','u(','um',
     *'ls','dz','tx','ty','ts','cx','cy','qx','qy','hh','vv','tt','hv',
     *'ht','vt','ax','bx','gx','ay','by','gy','at','bt','gt','ex','ey',
     *'et','wx','wy','wt','fx','xb','xa','xu','xd','fy','yb','ya','yu',
     *'yd','ft','tb','ta','tu','td'/
      data kask/2*4,6*3,5*2,43*1/
c-----!----------------------------------------------------------------!
ctm 4 Feb 2012: new keyword added for Taylor map: tm(i,j) = tumat(i,j)
