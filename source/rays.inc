      parameter (maxray=30000)
c zblock(i,k) kth coordinate of ray index i
c zi (k)      initial coordinates
c zf (k)      final coordinates
c nrays       number of rays
c istat(i)    iturn, in which ray index i was lost, 0 otherwise
c ihist(j,1)  iturn, in which the jth lost particle was lost
c ihist(j,2)  index (i) of jth lost particle
c iturn       turn index (not used in cqlate)
c nlost       number of rays lost so far
      common/rays/zblock(maxray,6),zi(6),zf(6),nrays,
     #            istat(maxray),ihist(maxray,2),iturn,nlost
