c mnumax = maximum number of elements in #menu
      parameter (mnumax=1024)
c lmnlbl(i) = name of ith element (numbered as they occur in #menu)
c p(k,i)    = kth parameter of ith element
c nt1(i)    = group index of the elements type (see /monics/)
c nt2(i)    = type index  .....
c mpp(i)    = menu parameter pointer for ith element
      common/elmntc/lmnlbl(mnumax)
      character*8 lmnlbl
      common/elmnts/pmenu(6*mnumax),mpp(mnumax),nt1(mnumax),nt2(mnumax)
c
c Global Menu pointer 
      common/imenu/inmenu
