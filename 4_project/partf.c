#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>


int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

void init(int L, int **lattice)
{
	double E, m, s;
  int i, j;
	double r;


  for (i = 0; i < L; i++){
    for (j = 0; j < L; j++){
      r = (double)rand()/RAND_MAX;
			printf("%lf\n", r);

      s = r > 0.5 ? 1 : -1;

      lattice[i][j] = s;
      m += (double)s;
    } // end i
  } // end j


  for (i = 0; i < L; i++){
    for (j = 0; j < L; j++){
      E -=  (double) spin_matrix[y][x]*
			(spin_matrix[periodic(y, L, -1)][x] +
	 		spin_matrix[y][periodic(x, L, -1)]);
    }
  }

}

void montecarlo()
{
  /* Parallelized Monte Carlo simulation with Metropolis sampling */


}

int main(int argc, char *argv[])
{
  srand((unsigned)time(NULL)) // seed

  int my_rank, num_procs;
  int i;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &num_procs);

  N = atoi(argv[1]); // number of Monte Carlo cycles

  my_N = ceil(N/num_procs); // number of MC cycles per worker

  /* temperature values to study */
  double *T = (double)malloc(6*sizeof(double*));
  for (i = 0; i < 6; i++) T[i] = 2.0 + 0.05*i;


	lattice = (int**)malloc(L*sizeof(int*));
  for (int i = 0; i < L; i++){
    lattice[i] = (int*)malloc(L*sizeof(int));
  }


}
