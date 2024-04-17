#include <algorithms/ClusterSinglesDoublesTriplesQuadruplesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/Tensor.hpp>
#include <Options.hpp>
#include <Sisi4s.hpp>

#include <initializer_list>

using namespace sisi4s;

void ClusterSinglesDoublesTriplesQuadruplesAlgorithm::run() {
  double e(0.0);
  if (in.is_of_type<Tensor<double> *>("PPHHCoulombIntegrals")) {
    run<double>();
  } else {
    run<complex>();
  }
}

template <typename F>
F ClusterSinglesDoublesTriplesQuadruplesAlgorithm::run() {
  int Nv(in.get<Tensor<double> *>("ParticleEigenEnergies")->lens[0]);
  int No(in.get<Tensor<double> *>("HoleEigenEnergies")->lens[0]);

  PTR(const FockVector<F>) amplitudes(
      createAmplitudes<F>({"Singles", "Doubles", "Triples", "Quadruples"},
                          {{Nv, No},
                           {Nv, Nv, No, No},
                           {Nv, Nv, Nv, No, No, No},
                           {Nv, Nv, Nv, Nv, No, No, No, No}},
                          {"ai", "abij", "abcijk", "abcdijkl"}));

  // create a mixer, by default use the linear one
  std::string mixerName(in.get<std::string>("mixer"));
  PTR(Mixer<F>) mixer(MixerFactory<F>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  // number of iterations for determining the amplitudes
  int maxIterationsCount(in.get<int64_t>("maxIterations"));

  F e(0);
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i + 1 << std::endl;
    // call the getResiduum of the actual algorithm,
    // which will be specified by inheriting classes
    auto estimatedAmplitudes(getResiduum(i, amplitudes));
    estimateAmplitudesFromResiduum(estimatedAmplitudes, amplitudes);
    auto amplitudesChange(NEW(FockVector<F>, *estimatedAmplitudes));
    *amplitudesChange -= *amplitudes;
    mixer->append(estimatedAmplitudes, amplitudesChange);
    // get mixer's best guess for amplitudes
    amplitudes = mixer->get();
    e = getEnergy(amplitudes);
  }

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation())
        << "computing energy from given amplitudes" << std::endl;
    e = getEnergy(amplitudes);
  }

  storeAmplitudes(amplitudes, {"Singles", "Doubles", "Triples", "Quadruples"});
  return e;
}
