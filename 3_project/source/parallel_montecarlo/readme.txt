The serial program "montecarlo_serial" was compiled using the gcc compiler. 

It can be run by supplying two command-line arguments: name of text file containing the resulting data, and the number of
Monte Carlo iterations to run for.

An example of how to compile and run the program:
-------------------------------------------------------
$ gcc montecarlo_serial.c -lm -o mc_serial.o

$ ./mc_serial.o output.txt 100000
--------------------------------------------------------


The parallel program "montecarlo_parallel" was compiled using the mpicc-compiler and can be run with the same command line arguments 
as the serial program, while also including the number of processors to run the program with.

An example of how to compile and run the program:
-------------------------------------------------------
$ mpicc montecarlo_parallel.c -lm -o mc_parallel.o

$ mpirun ./mc_parallel.o -np 4 output.txt 100000
--------------------------------------------------------
