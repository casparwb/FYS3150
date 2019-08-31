#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

double exact(double x)
{
  /* Exact analytical solution */
  double value = 1 - (1 - exp(-10))*x - exp(-10*x);
  return value;
}

void matvecmulti_gen(int n, double *a, double *b, double *c, double*v,
                   double *d)
{
  /* Solves the general algorithm */
  double w;
  int i;
  double *d_new = (double*)malloc((n)*sizeof(double));

  // copy d
  for (i = 0; i < n; i++) d_new[i] = d[i];

  // forward substitution
  for (i = 1; i < n; i++){
    w = a[i]/b[i-1];
    b[i] -= w*c[i-1];
    d_new[i] -= w*d_new[i-1];
  }

  for (i = 0; i<n; i++) printf("%lf\n", d_new[i]);


  // backward substitution
  v[n-1] = d_new[n-1]/b[n-1];
  for (i = n-2; i >= 0; i--) {
    v[i] = (d_new[i] - c[i]*v[i+1])/b[i];
  }

  free(d_new);
} // end of matvecmulti_gen

void matvecmulti_spec(int n, double *v, double *d)
{
  /* Solves the specialized algorithm */

  double *b;
  double one_over_b;
  int i;

  // precalculating 1/b[i]
  b = (double*)malloc(n*sizeof(double));
  for (i = 0; i < n; i++) b[i] = (i+2.)/(i+1.);
  //for (i = 1; i < n; i++) b[i] -= 1.0/b[i-1];

  // forward substitution
  for (i = 1; i < n; i++){
    d[i] += d[i-1]/b[i-1];
  }

  for (i = 0; i<n; i++) printf("%lf\n", d[i]);
  // backward substitution
  v[n-1] = d[n-1]/b[n-1];

  for (i = n-2; i >= 0; i--){
    v[i] = (d[i] + v[i+1])/b[i];
  }


} // end of matvecmulti_spec

int main(int argc, char *argv[]){

  if (argc < 3){
    printf("Please read output file name\n");
    exit(1);
  }

  int i, n = atoi(argv[1]);
  double t1, t2, h;       // step size
  clock_t start, stop;
  FILE *outfile = fopen(argv[2], "w+");


  /* allocating arrays */
  double *a = (double*)malloc((n)*sizeof(double));
  double *b = (double*)malloc((n)*sizeof(double));
  double *c = (double*)malloc((n)*sizeof(double));
  double *d = (double*)malloc((n)*sizeof(double));
  double *v1 = (double*)malloc((n)*sizeof(double));
  double *v2 = (double*)malloc((n)*sizeof(double));
  double *u = (double*)malloc((n)*sizeof(double));

  h = 1./(n + 1); // set step size

  // initializing arrays
  double hh = 100*h*h;
  for (i = 0; i <= n; i++){

    d[i] = hh*exp(-10*(i+1)*h);
    u[i] = exact((i+1)*h);

    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
    v1[i] = 0;
    v2[i] = 0;
  }

  a[0] = 0;
  c[n-1] = 0;

  // calling the general algorithm
  start = clock();
  matvecmulti_gen(n, a, b, c, v1, d); // perform the multiplication
  stop = clock();

  t1 = (double)((stop - start)*1000)/CLOCKS_PER_SEC;
  printf("\n");
  // calling the specialized algorithm
  start = clock();
  matvecmulti_spec(n, v2, d);
  stop = clock();

  t2 = (double)((stop - start)*1000)/CLOCKS_PER_SEC;


  // writing results to file
  fprintf(outfile, "n = %d\n", n+2);
  fprintf(outfile, "i | Exact value | General algorithm | Specialized algorithm\n");
  fprintf(outfile, "-------------------------------------------------------\n");
  fprintf(outfile, "%d | %lf      %lf            %lf     \n", 0, 0.0, 0.0, 0.0);
  for (i = 0; i < n; i++){
      fprintf(outfile, "%d | %lf      %lf            %lf     \n", i+1, u[i], v1[i], v2[i]);
  }
  fprintf(outfile, "%d | %lf      %lf            %lf     \n", n+1, 0.0, 0.0, 0.0);
  fprintf(outfile, "General time: %lf ms\nSpecialized time: %lf ms\n", t1, t2);

  fclose(outfile);

  // free allocated memory
  free(d); free(v1); free(v2);
  free(a); free(b); free(c); free(u);
  return 0;



}
