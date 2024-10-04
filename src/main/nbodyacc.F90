
program nbodyacc 

 use initial, only:init
 use evolve,  only:evol
 implicit none
 
 ! Initial checks
 ! infile and startdump exist?  

 call init()
 call evol()


end program 