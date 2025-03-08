// Macros copied over from
// https://github.com/CEED/libCEED/blob/main/include/ceed/types.h

#ifndef TYPES_H
#define TYPES_H

/**
  This macro defines compiler attributes to the BENCH_QFUNCTION to force
inlining for called functions. The `inline` declaration does not necessarily
enforce a compiler to inline a function. This can be detrimental to performance,
so here we force inlining to occur unless inlining has been forced off (like
during debugging).
**/
#ifndef BENCH_QFUNCTION_ATTR
#ifndef __NO_INLINE__
#if defined(__GNUC__) || defined(__clang__)
#define BENCH_QFUNCTION_ATTR __attribute__((flatten))
#elif defined(__INTEL_COMPILER)
#define BENCH_QFUNCTION_ATTR _Pragma("forceinline")
#else
#define BENCH_QFUNCTION_ATTR
#endif
#else
#define BENCH_QFUNCTION_ATTR
#endif
#if defined(__GNUC__) || defined(__clang__)
#define BENCH_QFUNCTION_HELPER_ATTR                                            \
  BENCH_QFUNCTION_ATTR __attribute__((always_inline))
#else
#define BENCH_QFUNCTION_HELPER_ATTR BENCH_QFUNCTION_ATTR
#endif
#endif

/**
  This macro populates the correct function annotations for User QFunction
source for code generation backends or populates default values for CPU
backends. It also creates a variable `name_loc` populated with the correct
source path for creating the respective User QFunction.
**/
#ifndef BENCH_QFUNCTION
#define BENCH_QFUNCTION(name)                                                  \
  static const char name##_loc[] = __FILE__ ":" #name;                         \
  BENCH_QFUNCTION_ATTR static int name
#endif

/**
  This macro populates the correct function annotations for User QFunction
helper function source for code generation backends or populates default values
for CPU backends.
**/
#ifndef BENCH_QFUNCTION_HELPER
#define BENCH_QFUNCTION_HELPER BENCH_QFUNCTION_HELPER_ATTR static inline
#endif

// This macro provides the appropriate SIMD Pragma for the compilation
// environment.
#ifndef BenchPragmaSIMD
#if defined(__INTEL_COMPILER)
#define BenchPragmaSIMD _Pragma("vector")
/// Cannot use Intel pragma ivdep because it miscompiles unpacking symmetric
/// tensors, as in Poisson2DApply, where the SIMD loop body contains temporaries
/// such as the following.
///
///     const double dXdxdXdxT[2][2] = {{qd[i+0*Q], qd[i+2*Q]},
///                                         {qd[i+2*Q], qd[i+1*Q]}};
///     for (int j=0; j<2; j++)
///        vg[i+j*Q] = (du[0] * dXdxdXdxT[0][j] + du[1] * dXdxdXdxT[1][j]);
///
/// Miscompilation with pragma ivdep observed with icc (ICC) 19.0.5.281 20190815
/// at -O2 and above.
#elif defined(__GNUC__) && __GNUC__ >= 5
#define BenchPragmaSIMD _Pragma("GCC ivdep")
#elif defined(__clang__)
#define BenchPragmaSIMD _Pragma("clang loop vectorize(enable)")
#elif defined(_OPENMP) && _OPENMP >= 201307 // OpenMP-4.0 (July, 2013)
#define BenchPragmaSIMD _Pragma("omp simd")
#else
#define BenchPragmaSIMD
#endif
#endif

#endif // TYPES_H
