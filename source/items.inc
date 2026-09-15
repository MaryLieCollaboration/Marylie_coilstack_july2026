c itmno  = maximum number of items (lines, lumps, loops)
c itmnln = maximum number of elements per item
      parameter (itmno=300,itmln=1024)
c ilbl(i)   = name of ith item
c icon(k,i) = name of kth element/item in ith item
c irep(k,i) = repetition factor of icon(k,i)
c ityp(k,i) = type code = 2,3,4 if item=line,lump,loop resp.
c ilen(i)   = number of elements/items in ith item (=max k for ith item)
c lmade(i)  = flag, telling whether map of lump is in buffer and where.
c             preset to zero 
      common/itemsc/ilbl(itmno),icon(itmln,itmno)
      character*8 ilbl,icon
      common/items/ irep(itmln,itmno),ityp(itmno),
     &              ilen(itmno),lmade(itmno)
