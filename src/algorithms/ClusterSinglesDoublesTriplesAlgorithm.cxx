#include <algorithms/ClusterSinglesDoublesTriplesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>
#include <util/Emitter.hpp>
#include <util/Exception.hpp>
#include <util/Tensor.hpp>
#include <Options.hpp>
#include <Sisi4s.hpp>

#include <initializer_list>

using namespace sisi4s;

void ClusterSinglesDoublesTriplesAlgorithm::run() {
  if (in.is_of_type<Tensor<double> *>("PPHHCoulombIntegrals")) {
    run<double>();
  } else {
    run<complex>();
  }
}

template <typename F>
F ClusterSinglesDoublesTriplesAlgorithm::run() {
  int Nv(in.get<Tensor<double> *>("ParticleEigenEnergies")->lens[0]);
  int No(in.get<Tensor<double> *>("HoleEigenEnergies")->lens[0]);

  PTR(const FockVector<F>) amplitudes(createAmplitudes<F>(
      {"Singles", "Doubles", "Triples"},
      {{Nv, No}, {Nv, Nv, No, No}, {Nv, Nv, Nv, No, No, No}},
      {"ai", "abij", "abcijk"}));

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

  F amplitudesConvergence(in.get<double>("amplitudesConvergence"));
  F energyConvergence(in.get<double>("energyConvergence"));

  EMIT() << YAML::Key << "maxIterations" << YAML::Value << maxIterationsCount
         << YAML::Key << "amplitudesConvergence" << YAML::Value
         << std::abs(amplitudesConvergence) << YAML::Key << "energyConvergence"
         << YAML::Value << std::abs(energyConvergence);

  EMIT() << YAML::Key << "iterations" << YAML::Value;
  EMIT() << YAML::BeginSeq;

  F e(0), previousE(0);
  int i(0);
  for (; i < maxIterationsCount; ++i) {
    EMIT() << YAML::BeginMap;
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i + 1 << std::endl;
    EMIT() << YAML::Key << "iteration" << YAML::Value << i + 1;
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
    if (std::abs((e - previousE) / e) < std::abs(energyConvergence)
        && std::abs(amplitudesChange->dot(*amplitudesChange)
                    / amplitudes->dot(*amplitudes))
               < std::abs(amplitudesConvergence * amplitudesConvergence)) {
      EMIT() << YAML::EndMap;
      break;
    }
    previousE = e;
    EMIT() << YAML::EndMap;
  }
  EMIT() << YAML::EndSeq;

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation())
        << "computing energy from given amplitudes" << std::endl;
    e = getEnergy(amplitudes);
  } else if (i == maxIterationsCount) {
    LOG(0, getCapitalizedAbbreviation())
        << "WARNING: energy or amplitudes convergence not reached."
        << std::endl;
  }

  storeAmplitudes(amplitudes, {"Singles", "Doubles", "Triples"});
  return e;
}
