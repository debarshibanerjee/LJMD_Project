#Report
In this report we describe the steps followed for benchmarking the LJMD code, and the results obtained.

The benchmark process has been divided in 3 tasks:
 - single-core optimization
 - MPI parallelization
 - OMP parallelization


#How to compile and run the code

A cmake version higher than 3.15 is required to run the code.


#Single-core optimization

This part focuses on the optimization of the serial code. 

We have advanced step by step:

Starting from the baseline code provided by A. Kohlmeier, we first improve the performances by exploiting the compiler optimization flags.

Then, we have considered different optimizations such as:
  - substituting costly mathematical operations with faster ones (for instance reducing the number of calls to sqrt and division operations as possible)
  - improving the data structure (to improve the memory access and avoid memory aliasing)
  - loop unrolling (in azzero function)
  - avoid repetitive operations, with the introduction of prefactors (in ekin, verlet and force function)
  


Then, we have applied our physical knowledge, and we optimize further the force calculation by means of the exploitations of the Newton's third law.


Clearly, we have implemented these improvements having in mind that these shall be mergeable with MPI and OpenMP. 


#MPI benchmark

A benchmark for MPI has been performed measuring computation times for [1, 32] cores for both argon_108.inp and argon_2916.inp. The benchmark took place on 2 computational nodes on the Ulysses cluster.


The results are presented below, both in graphical and tabular form. 

#Hybrid benchmark (with OMP)


