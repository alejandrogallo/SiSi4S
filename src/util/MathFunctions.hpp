#ifndef MATH_FUNCTIONS_DEFINED
#define MATH_FUNCTIONS_DEFINED

#include <util/Complex.hpp>
#include <cmath>
#include <ctf.hpp>

/**
 * \brief Common inline math functions.
 * As opposed to the variants in the std namesapce, these functions
 * are type closed, i.e. return the same type as the arguments,
 * which is required by Tensor::sum for univariate functions and
 * by Tensor::contract for bivariate functions.
 */
namespace cc4s {
  // univariate functions
  template <typename F=double>
  inline F sqrt(F const x) {
    return std::sqrt(x);
  }

  template <typename F=double>
  inline F abs(F const x) {
    return std::abs(x);
  }

  template <typename F=double>
  inline F conj(F const x) {
    return std::conj(x);
  }

  inline double conj(double const x) {
    return x;
  }

  // bivariate functions
  template <typename F=double>
  inline F dot(F const x, F const y) {
    return x * conj(y);
  }

  /**
   * \brief Calculates only the real part of x*conj(y).
   */
  template <typename F=double>
  inline F realDot(F const x, F const y) {
    return std::real(x*conj(y));
  }

  template <typename F=double>
  inline F divide(F const x, F const y) {
    return x / y;
  }

  template <typename F>
  inline double frobeniusNorm(CTF::Tensor<F> &t) {
    char *indices(new char[t.order+1]);
    for (int index(0); index < t.order; ++index) indices[index] = 'a' + index;
    indices[t.order] = 0;
    CTF::Bivar_Function<F> fRealDot(&cc4s::realDot<F>);
    CTF::Scalar<F> s(*t.wrld);
    s.contract(1.0, t,indices, t,indices, 0.0,"", fRealDot);
    return std::sqrt(std::real(s.get_val()));
  }
}

#endif
