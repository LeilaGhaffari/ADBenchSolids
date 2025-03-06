
#include <stdio.h>

extern float func_(float *x);

int main() {
  float x = 1.5;
  float f = func_(&x);
  printf("here is f = %f\n", f);
  return 0;
}

// gfortran -c test001.f90 -o test001_f.o
// gcc -c test001.c -o test001_c.o
// gcc test001_f.o test001_c.o -lgfortran -o test001_c.out
// ./test001_c.out