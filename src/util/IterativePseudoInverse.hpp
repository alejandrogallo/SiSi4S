#ifndef ITERATIVE_PSEUDO_INVERSE_DEFINED
#define ITERATIVE_PSEUDO_INVERSE_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>
#include <random>

template <typename F>
class IterativePseudoInverse {
public:
  IterativePseudoInverse(CTF::Matrix<F> const &matrix, double accuracy = 1e-10);
  CTF::Matrix<F> &get();

  static void test(CTF::World *world);
protected:
  void iterate(double accuracy = 1e-10);
  void iterateQuadratically(double accuracy = 1e-10);
  static void setRandom(
    F &value,
    std::mt19937 &random, std::normal_distribution<double> &normalDistribution
  );
  static void generateHilbertMatrix(CTF::Matrix<F> &matrix);

  CTF::Matrix<F> matrix, square, inverse;
  double alpha;
};

#endif