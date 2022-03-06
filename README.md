# LJMD_Project - Group 2:

![example workflow](https://github.com/debarshibanerjee/LJMD_Project/actions/workflows/linux-build-and-test.yml/badge.svg)

## Participants

The commits in the repository correlate to usernames and people as follows:

- debarshibanerjee - Debarshi Banerjee
- Sera91 - Serafina Di Gioia
- AlexTru96  - Alexander Trujillo


The report, with the description of the procedure followed for the code optimization, and the results of the bencmark,
can be found in directory Report.

## Running the program:

The `compile.sh` script runs the program by compiling it from scratch in a `build` directory. It can be edited to run according to user expectations as explained below.

It takes an argument - `1` or `0`, which respectively enables or disabled OpenMP support.

```
USAGE: ./compile.sh 1
```

The code is always compiled with full MPI support.

The script sets a default number of OpenMP threads using `export OMP_NUM_THREADS = 4`. This can also be modified locally in one's interactive shell of choice.

This can be modified as needed, of course.

The `cmake` command in the script can be modified as follows:
- `ENABLE_TESTING` (`yes` or `no`. Default: `no`)
- `ENABLE_SIMPLE_FORCE` (`yes` or `no`. Default: `no`)
- `ENABLE_OPENMP` (`1` or `0`. Taken from script parameters.)

Testing is disabled in the script by default since it is already automatically carried out by GitHub actions.

The `SIMPLE_FORCE` directive is for using a force function, where the OpenMP parallelization (with an MPI background) is purely a `#pragma omp parallel for` and a `reduction` directive for the main `for` loop.

The default force function uses the OpenMP custom-reduction approach, in addition to MPI, to create an efficient hybrid function.

The `test_verification.sh` script runs `ctest` to verify the tests if `ENABLE_TESTING=yes` is set.

The `result_verification.sh` script runs `make check` in the `examples` directory for 108 and 2916 size Argon systems and compares it with the reference results.
