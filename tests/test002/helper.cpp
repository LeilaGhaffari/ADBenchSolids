extern "C" {
float func_(float *x);
}

float helping(float *x) {
  float result = func_(x);
  return result;
}