************************************************************************
* header:                 CONSTRAINTS (CONS)                           *
* Constraint subroutines to be used in conjunction with fitting        *
* routines.                                                            *
************************************************************************
      subroutine con1(p)
c impose constraints amoung parameters in the various parameter sets
      include 'impli.inc'
      include 'files.inc'
      include 'parset.inc'
c
      dimension p(6)
c
c This constraint routine sets the third parameter equal to the
c first parameter in any given parameter set.
c It is useful in fitting a REC triplet for which the two outside
c quads are to be tied together.
c
c select parameter set
      ipset=nint(p(1))
      if((ipset.lt.1) .or. (ipset.gt.maxpst)) then
      write (jof,*) 'WARNING: ipset out of range in command with',
     # ' type code con1'
      return
      endif
c      
c set third parameter equal to the first
      pst(3,ipset)=pst(1,ipset)
      return
      end
c
************************************************************************
c
      subroutine con2(p)
c impose constraints amoung parameters in the various parameter sets
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
c
      dimension p(6)
c
c this constraint relates parameter set p2 to p1 by
c putting parameter p4 of parameter set p2 equal to minus 
c parameter p3 of parameter set p1.
c
c set up control indices
      ipset1=nint(p(1))
      ipset2=nint(p(2))
      ipar1=nint(p(3))
      ipar2=nint(p(4))
c
c relate parameter values
      pst(ipar2,ipset2)=-pst(ipar1,ipset1)
c
      return
      end
c
************************************************************************
c
      subroutine con3(p)
c impose constraints amoung parameters in the various parameter sets
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
c      write(6,*) 'con3 not yet installed'
      ipar=1
      pst(ipar,4)=(1.5-pst(ipar,3))/2.
      return
      end
c
************************************************************************
c
      subroutine con4(p)
c Impose constraints amoung parameters in the various parameter sets.
c This routine imposes constraints among the strengths of 4 octupoles
c in such a way that the sum of the squares of their strengths
c is minimized.
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
      include 'map.inc'
c
c calling array
      dimension p(6)
c
c local arrays
      dimension s(4)
      dimension v1(4),v2(4),v3(4)
      save v1, v2, v3
c
c set up parameters
      ipass=nint(p(1))
      scale=p(2)
      islot=nint(p(3))
      beta1=p(4)
      beta2=p(5)
      beta3=p(6)
c
c assemble the vectors v1,v2,v3
      if (ipass.le.4) then
      v1(ipass)=th(84)/scale
      v2(ipass)=th(95)/scale
      v3(ipass)=th(175)/scale
      endif
c
c debugging code
      if (ipass.eq.4) then
      write(12,*) v1(1),v2(1),v3(1)
      write(12,*) v1(2),v2(2),v3(2)
      write(12,*) v1(3),v2(3),v3(3)
      write(12,*) v1(4),v2(4),v3(4)
      endif
c
c calculate and set octupole strengths
      if (ipass.gt.4) then
      do 10 i=1,4
      s(i)=beta1*v1(i)+beta2*v2(i)+beta3*v3(i)
      pst(islot,i)=s(i)
   10 continue
      endif
c
      return
      end
c
************************************************************************
c
      subroutine con5(p)
c Impose constraints amoung parameters in the various parameter sets.
c This routine imposes constraints among the strengths of 5
c x3 kicks in such a way that the sum of the squares of their
c strengths is minimized.
c
      include 'impli.inc'
      include 'param.inc'
      include 'parset.inc'
      include 'map.inc'
c
c calling array
      dimension p(6)
c
c local arrays
      dimension s(5)
      dimension v1(5),v2(5),v3(5),v4(5)
      save v1, v2, v3, v4
c
c set up parameters
      ipass=nint(p(1))
      scale=p(2)
      ibeta=nint(p(3))
      ians=nint(p(4))
c
c assemble the sensitivity vectors v1,v2,v3,v4
      if (ipass.le.5) then
      v1(ipass)=th(28)/scale
      v2(ipass)=th(29)/scale
      v3(ipass)=th(34)/scale
      v4(ipass)=th(49)/scale
      endif
c
c debugging code
      if (ipass.eq.5) then
      write(12,*) v1(1),v2(1),v3(1),v4(1)
      write(12,*) v1(2),v2(2),v3(2),v4(2)
      write(12,*) v1(3),v2(3),v3(3),v4(3)
      write(12,*) v1(4),v2(4),v3(4),v4(4)
      write(12,*) v1(5),v2(5),v3(5),v4(5)
      endif
c
c calculate and set kick strengths
      if (ipass.gt.5) then
      beta1=pst(1,ibeta)
      beta2=pst(2,ibeta)
      beta3=pst(3,ibeta)
      beta4=pst(4,ibeta)
      do 10 i=1,5
      s(i)=beta1*v1(i)+beta2*v2(i)+beta3*v3(i)+beta4*v4(i)
      pst(i,ians)=s(i)
   10 continue
      endif
c
      return
      end
c
c end of file
