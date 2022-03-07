# Some details on Morse potential:

## Author: Debarshi Banerjee

### General Approach:
The approach taken here was to use the Morse potential to replace Lennard-Jones Potential.

Some practical assumptions:
- The depth of the potential well, &epsilon;, in Lennard-Jones, is the same as D in Morse potential.
- r <sub>e</sub> is the equilibrium bond distance in Morse potential, and that can be taken to be the same as 2<sup>1/6</sup>&sigma; (&sigma; is the same as in LJ)

Using these assumptions, as done in the script `fit_lj_morse.py`, we get a good fitting value of &alpha;, the width of the well in Morse potential. This is done using `curve_fit` from `scipy.optimize`.

This turns out to be &alpha; = 1.7473 and the curve we get is a very good fit:
![fit_morse_lj](https://user-images.githubusercontent.com/45908114/156945969-bc537ddc-0d9c-4fde-bb76-754795c1a0aa.png)


However, when we use this in our `force_morse()` function, we get unexpected results with positive potential energies.

We manually adjust the parameters so that the results given by `force_morse()` are atleast somewhat similar to those given for the LJ potential.

In this process, we discover a good value of &alpha; = 1.05

Now, to see if this is a computational artifact or if it makes physical sense, we plot the LJ and Morse potentials simultaneously, using `plot_morse_lj.py`. 

As we see here, it is still a pretty good fit and the overall nature of the potential is preserved as expected.
![morse_lj_105](https://user-images.githubusercontent.com/45908114/156945927-c960c729-9696-4b24-b6ce-7fb448af6afb.png)

### Further Comments:
Probably the reason that the value of &alpha; derived from the exact comparison is not a good fit is because it places too much emphasis during the fitting process on the highly-positive left-hand side of the potential well, where there are very high values, and also in the asymptotic tail, and while the algorithm tries to minimize the error in these regions, it loses sight of the main aread of interest - in and around the potential well.

Since such problems are highly sensitive to initial value conditions, this is not surprising as far as I understand.

It must be mentioned, that the results for both 108 and 2916 size clusters are very close to LJ values (within 2-3% error past the few starting values).

### Running the code with Morse potential:
The function is implemented in a separate file in `src/force_morse.c`. 

This can be enabled with the `-DENABLE_MORSE=1` flag as shown in the `compile.sh` script. 

Of course, it is turned off by default in `compile.sh` (set to `0`). 

To see the results generated once the code is compiled with this flag turned on, we go to `examples`, and there we do `make -f makefile_morse.mk`, to run the executable for both 108 and 2916 size Argon systems. 

The byte-wise comparison with reference results (using `cmp`) is turned off for obvious reasons, as results are not at all the same. At best the error is very small.




