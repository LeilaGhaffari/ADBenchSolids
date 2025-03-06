#ifndef TYPES_H
#define TYPES_H

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
