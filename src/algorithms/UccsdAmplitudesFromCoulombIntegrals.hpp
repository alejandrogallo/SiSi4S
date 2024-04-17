#ifndef ___CLUSTERSINGLESDOUBLESALGORITHM_DEFINED_PEACE
#define ___CLUSTERSINGLESDOUBLESALGORITHM_DEFINED_PEACE

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#define UCCSD_SPEC_IN                                                          \
  CLUSTER_SINGLES_DOUBLES_INSPEC,                                              \
      {"intermediates", SPEC_VALUE_DEF("TODO: DOC", bool, true)},              \
      {"HoleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},        \
      {"ParticleEigenEnergies", SPEC_VARIN("TODO: DOC", Tensor<double> *)},    \
      {"HHFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},              \
      {"HPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},              \
      {"PPFockMatrix", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},              \
      {"HHHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HHPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HPHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HPPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"HPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PHPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PPHHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PPHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)},      \
      {"PPPHCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)}, {    \
    "PPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<FSPEC> *)           \
  }
#define UCCSD_SPEC_OUT CLUSTER_SINGLES_DOUBLES_OUTSPEC

namespace sisi4s {
class UccsdAmplitudesFromCoulombIntegrals
    : public ClusterSinglesDoublesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UccsdAmplitudesFromCoulombIntegrals);
  using ClusterSinglesDoublesAlgorithm::ClusterSinglesDoublesAlgorithm;
  static AlgorithmInputSpec spec;

  virtual void run();
  virtual std::string getAbbreviation() { return "Uccsd"; }

protected:
  /**
   * \brief Implements the iterate method with the DRCCD iteration.
   * \param[in] i Iteration number
   */
  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes);

  bool usingIntermediates;
  bool onlyPpl;

  virtual PTR(FockVector<complex>)
  getResiduum(const int iteration,
              const PTR(const FockVector<complex>) &amplitudes);

  template <typename F>
  PTR(FockVector<F>)
  getResiduumTemplate(const int iteration,
                      const PTR(const FockVector<F>) &amplitudes);
};

} // namespace sisi4s

#endif
