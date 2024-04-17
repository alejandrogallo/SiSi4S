#ifndef CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_ALGORITHM_DEFINED
#define CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/FockVector.hpp>
#include <DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <util/Tensor.hpp>

#include <string>
#include <initializer_list>

namespace sisi4s {

#define CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_INSPEC                      \
  CLUSTER_SINGLES_DOUBLES_TRIPLES_INSPEC, {                                    \
    "initialTriplesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<F> *)           \
  }

#define CLUSTER_SINGLES_DOUBLES_TRIPLES_QUADRUPLES_OUTSPEC                     \
  CLUSTER_SINGLES_DOUBLES_TRIPLES_OUTSPEC, {                                   \
    "QuadruplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)              \
  }

class ClusterSinglesDoublesTriplesQuadruplesAlgorithm
    : public ClusterSinglesDoublesAlgorithm {

public:
  using ClusterSinglesDoublesAlgorithm::ClusterSinglesDoublesAlgorithm;
  virtual ~ClusterSinglesDoublesTriplesQuadruplesAlgorithm() {}
  virtual void run();
  virtual std::string getAbbreviation() = 0;

protected:
  template <typename F>
  F run();
};
} // namespace sisi4s

#endif
