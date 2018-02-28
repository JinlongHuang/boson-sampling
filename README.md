# boson_sampling

In this simulation software we simulated boson sampling from 6 photons up to 10 photons. There are total 5 files from test6.m to test10.m. The number in the end of the file name stands for the number of photons in that simulation. In all files, the transition amplitude from initial state to final state is calculated using both Glynn's algorithm and Gurvits's algorithm. 

For Glynn's algorithms, it has four parameters that can be tuned to a specific simulation instance:  initial state, final state, # of photons n and # of modes m. In the following, specific instructions to adjust these parameters are given. 

The initial state: In the definition of variable "mGenGly", the exponents of Z(i)' specify how the initial state are arranged. For example, if the initial state is set to be [1,2,0,3], i.e. total 6 photons in 4 modes, mGenGly should be defined as mGenGly=Z(i0)'*Z(i1)'^2*Z(i3)'^3*... . The prime in Matlab means complex conjugate. Also, according to the formula for transition amplitude, mGenGly should be further multiplied by a factor.
  
The final state: There is only one place to be modified for final state: the exponents of dot(randU(i,:),ZVec), similarly as the initial state. If the final state is also [1,2,0,3], mGenGly should be defined as mGenGly=\cdots * dot(randU(1,:),ZVec)^1*dot(randU(2,:),ZVec)^2*dot(randU(3,:),ZVec)^0*dot(randU(4,:),ZVec)^3.

\# of photons n: Three places should be modified: 
1) The set of Z that Z(i) can choose from, which is defined before While-loop 
2) The end point for each While-loop 
3) Denominator of TransAmp.

\# of modes m: There are four places should be modified: 
1) \# of While-loop 
2) length of ZVec   
3) \# of dot(randU(i,:),ZVec)
4) Denominator of TransAmp.

For Gurvits's algorithm, besides above four parameters, there is one more parameter one can adjust: \# of samples T. This variable is defined in the beginning of Gurvits's algorithm, and this is the only place should be changed.

To calculate the success probability of Gurvits's algorithm, two additional parameters should be specified: tolerance for error "TOL" and # of experiments "maxExp". In the case of 6 photons, TOL is set from 0.05 to 0.015, and maxExp = 1000. The figure for success probability vs Fidelity is plotted in the end of this report.

The executions of Glynn's algorithm in the case of 8 to 10 photons was not complete, for the time it requires in a standard laptop are 45.8 minutes, 20.6 hours and 27.7 days.

There are two ways to implement a random unitary according to Haar measure. The first one is import from file "RandomUnitary", which is based on file "opt\_args.m". The other method is to firstly implement a random Hermitian matrix from rand(n) method provided in Matlab, and then exponentiate it to get a random unitary matrix.

The number of modes are supposed to be O(n^2) to illustrate quantum supremacy, however in here
for the purpose of demonstrating the efficiency of Gurvits's algorithm, we simplify the simulation by constraining m = n. 

Further suggestions are welcome.

Jin-Long Huang
