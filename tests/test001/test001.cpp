#include <cstdio>
#include <iostream>

extern "C" {
float func_(float *x);
}

int main() {
  float x = 1.5;
  float result = func_(&x);
  std::printf("Here is f = %f\n", result);
  return 0;
}

// gfortran -c test001.f90 -o test001_f.o
// g++ -c test001.cpp -o test001_cpp.o
// g++ test001_f.o test001_cpp.o -lgfortran -o test001_cpp.out
// ./test001_cpp.out
