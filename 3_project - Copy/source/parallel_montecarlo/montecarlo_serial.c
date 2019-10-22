#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define pi 3.14159
double foo_spherical(double r1, double r2,
                     double theta1, double theta2,
                     double phi1, double phi2){

  // function to integrate (in spherical coords)

  double cosbeta = cos(theta1)*cos(theta2)
                 + sin(theta1)*sin(theta2)
                 * cos(phi1 - phi2);

  double result;
  double r1_sq = r1*r1;
  double r2_sq = r2*r2;
  double dist = r1_sq + r2_sq - 2*r1*r2*cosbeta;

  if (dist < 0 || dist < 1e-8)
  {
    result = 0.;
  }
  else
  {
    result = r1_sq*r2_sq*sin(theta1)*sin(theta2)
                * exp(-4*(r1 + r2))/sqrt(dist);
  }

  return result;
}

double PDF(double x){
  // Probability distribution function

  return 4*exp(-4*x);
}

double montecarlo(double N, double *x, double *stddev){

  double summ1 = 0.;
  double summ2 = 0.;
  double y;
  double jacobian = 4.0 * pi * pi * pi * pi;
  double r1, r2, theta1, theta2, phi1, phi2;

  for (int i = 0; i < N; i++){
    for (int k = 0; k < 6; k++) x[k] = (double)(rand())/RAND_MAX;

    r1 = -0.25*log(1 - x[0]);
    r2 = -0.25*log(1 - x[1]);

    theta1 = pi*x[2];
    theta2 = pi*x[3];

    phi1 = 2.*pi*x[4];
    phi2 = 2.*pi*x[5];

    y = foo_spherical(r1, r2, theta1, theta2, phi1, phi2)/(PDF(r1)*PDF(r2));
    summ1 += y;
    summ2 += y*y;
  } // end for

  summ1 /= N;
  summ2 /= N;


  *stddev = jacobian*sqrt((summ2 - summ2*summ2)/N);
  return summ1*jacobian;


}

int main(int argc, char *argv[]){
  /*
  Serial implementation of the Monte Carlo integral approximation
  with importance sampling
  */
  srand(clock());

  double tot_sum, stddev;
  clock_t t, t_tot = 0;
  double exact = 5.*pi*pi/(16*16);
  double error;

  int N = atoi(argv[1]);
  int m = 100;

  char *filename = argv[2];
  FILE *f = fopen(filename, "w");
  fprintf(f, "log10 N | Total Time | log10 error | standard dev.\n");

  double *x = (double*)malloc(N*sizeof(double));

  // perform m experiments
  for (int j = 0; j < m; j++){

    t = clock();
    tot_sum = montecarlo(N, x, &stddev);
    t_tot = clock() - t;

    error = fabs(tot_sum - exact)/tot_sum;
    // printf("log10 N = %.0lf\n", log10(N));
    // printf("Total time: %lf s\n", t/(double)CLOCKS_PER_SEC);
    // printf("log error: %lf \n", log10(error));
    // printf("Standard dev.: %lf\n", stddev);


    fprintf(f, "%lf | %lf | %lf | %lf\n", log10(N), t_tot/(double)CLOCKS_PER_SEC, log10(error), stddev);

  }

  fclose(f);
  free(x);


  return 0;




}
