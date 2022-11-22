/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all
 * rights reserved.*/
#ifndef UCCSDTQ_FROM_COULOMB_INTEGRALS_DEFINED
#define UCCSDTQ_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterSinglesDoublesTriplesQuadruplesAlgorithm.hpp>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>

namespace sisi4s {
class UccsdtqAmplitudesFromCoulombIntegrals
    : public ClusterSinglesDoublesTriplesQuadruplesAlgorithm {
public:
  ALGORITHM_REGISTRAR_DECLARATION(UccsdtqAmplitudesFromCoulombIntegrals);
  UccsdtqAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList);
  virtual ~UccsdtqAmplitudesFromCoulombIntegrals();

  virtual void run();
  virtual std::string getAbbreviation() { return "Uccsdtq"; }

protected:
  /**
   * \brief Implements the iterate method with the DRCCD iteration.
   * \param[in] i Iteration number
   */
  virtual PTR(FockVector<double>)
  getResiduum(const int iteration,
              const PTR(const FockVector<double>) &amplitudes);

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
