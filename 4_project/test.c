#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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

}

int main()
{
  int **lattice;
	int L = 20;

	lattice = (int**)malloc(L*sizeof(int*));
  for (int i = 0; i < L; i++){
    lattice[i] = (int*)malloc(L*sizeof(int));
  }

  init(L, lattice);

	for (int i = 0; i < L; i++){
		printf("%d\n", lattice[i][0]);
	}

}
