#include <stdio.h>
#include <stdlib.h>

void somefunc(double *a){

  for (int i = 0; i < 5; i++) a[i] = i;
}

int main(){


  double *a = (double*)malloc(10*sizeof(double));

  for (int i = 0; i < 10; i++) a[i] = 0;

  somefunc(a);

  for (int i = 0; i < 10; i++) printf("%lf\n", a[i]);

  return 0;
}
