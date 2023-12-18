#ifndef ___CLUSTERSINGLESDOUBLESALGORITHM_DEFINED_PEACE
#define ___CLUSTERSINGLESDOUBLESALGORITHM_DEFINED_PEACE

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

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
