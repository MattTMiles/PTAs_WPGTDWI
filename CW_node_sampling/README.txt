In this directory I'm going to attempt to find an efficient 
sampling algorithm for the phase connected case of CW inference. 

This is difficult as the degeneracy of the phase and distance 
properties result in N nodes (1e50) across the likelihood space,
all of which are potential solutions. 

Justin Ellis created an algorithm to tackle this years ago, Neil
Cornish asserts that this can't be tackled well because of the
modes. So let's see who's correct.

I have a data creation script already, now it's just a process 
of finding an algorithm that works. I think we should start with
Justin's theory and go from there. 