#include <algorithms/TensorNorm.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;


ALGORITHM_REGISTRAR_DEFINITION(TensorNorm);

TensorNorm::TensorNorm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorNorm::~TensorNorm() {
}

/**
 * \brief Testing environement
 */
void TensorNorm::run() {
  Tensor<> *A(getTensorArgument("Data"));
  double norm(frobeniusNorm(*A));
  LOG(0, "TensorNorm") << "norm = " << norm << std::endl;
  if (isArgumentGiven("Norm")) setRealArgument("Norm", norm);
}
