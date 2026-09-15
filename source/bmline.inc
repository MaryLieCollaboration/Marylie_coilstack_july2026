c mlabor = maximum number of tasks in #labor
c latt(i) = the name of the ith task 
c num(i) = the repetition factor of this task
      parameter (mlabor=1000)
      common/bmlinc/latt(mlabor)
      character*8 latt
      common/bmline/num(mlabor)
