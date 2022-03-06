# Report of the LJMD Project
In this report we describe the steps followed for benchmarking the LJMD code, and the results obtained.

The benchmark process has been divided in 3 tasks:
 - refactoring/single-core optimization (assigned to S. Di Gioia)
 - MPI parallelization/cell-list (assigned to A. Trujillo )
 - OMP parallelization/tests (assigned to D. Banerjee)


## Compilation support for the LJMD code

A cmake version higher than 3.15 is required to run the code.


## Single-core optimization

This part focuses on the optimization of the serial code. 

We have advanced step by step:

Starting from the baseline code provided by A. Kohlmeier, we first improve the performances by exploiting the compiler optimization flags, moving from -O2 to -O3 and adding the flags: -fomit-frame-pointer -msse3

Then, we have considered different optimizations such as:
  - substituting costly mathematical operations with faster ones (for instance reducing the number of calls to sqrt and division operations as possible)
  - improving the data structure (to improve the memory access and avoid memory aliasing)
  - loop unrolling (in azzero function)
  - avoid repetitive operations, with the introduction of prefactors (in ekin, verlet and force function)
  


As last step, we have applied our physical knowledge, and we optimize further the force calculation by means of the exploitations of the Newton's third law.


Below you can see how the runtime for Argon-108 changes with the step-by-step optimization:

![SCO_scaling.png](../plots/SCO-scaling.png)



Once we have reached the point where is difficult to improve further our SCO optimization, we considered this done and we moved to MPI and OpenMP parallelization.
It's important to notice that our efforts have been focused on optimizing the force function and the Verlet functions, more than other utilities functions, since these two functions represent the core of the code.


## MPI benchmark

A benchmark for MPI has been performed measuring computation times for [1, 32] cores for both argon_108.inp and argon_2916.inp. The benchmark took place on 2 computational nodes of the regular partiton 1 on the Ulysses cluster [at SISSA].
We performed 30 runs per each MPI configuration and select 2,4,8,16 cores on single node, and 24 and 32 corse divided by 2 nodes.

The results are presented below, both in graphical and tabular form. 

In the first plot we show the average time scaling of our code, with respect to the number of MPI tasks, where each MPI task is assigned to one distinct physical core.


![time_scaling_MPI.png](../plots/scaling_MPI-1.png)

We can observe that the code shows a better scaling for the simulation with a larger number of atoms,so with bigger size, and fewer steps  (Argon-2916).
This is even more evident in the plot with the average estimated SPEEDUP (S= T(1) /T(N) , where T(1): serial run time, or run time measured with 1 MPI task, and T(N): run time measured with N MPI tasks).

![speedup_plot.png](../plots/SPEEDUP_MPI_LJMD.png)



## Hybrid benchmark (with OMP)

The code presents two different OMP implementation of the Force function, one called the 'default' one and one the 'fancy' one.

The default OMP implementation did not seem to have a very good scaling with increasing the number of threads, like shown in the plot below:

![OMP_scaling.png](../plots/OMP_scaling_LJMD.png)

This has motivated us to implement an alternative OMP approach.
And with this approach we can gain a quite better scaling, but within the 10 % of difference respect to the previous approach.

We think that we need a larger size of the simulation to study properly the effect of the OMP parallelization but this is left for future studies.




