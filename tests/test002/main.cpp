#include <cstdio>
#include <iostream>

float helping(float *x);

int main() {
  float x = 1.5;
  float help = helping(&x);
  std::printf("Here is help = %f \n", help);
  return 0;
}
