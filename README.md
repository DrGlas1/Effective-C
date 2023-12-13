# Effective-C
Contains the projects for the course Effective C at LTH. Each folder contains a different working project for the mandatory parts of the course, where each pass the requirements on foresete.cs.lth.se.

## Intopt
The folder intopt contains a readable version of the simplex algorithm in the file `simplex.c` and the integer solver in the file `intopt.c`. Run either if these from `main.c`, an example of inputs can be found in the file i.
When compiled, run `./main < i` to run the algorithm with the inputs from `i`.

The folder also contains a file `forsete.c` which is a single-threaded optimized version of the integer solver, read the pdf file in the folder for the specific performance changes and their impacts on the performance of the program.
