#ifndef FLOPS_COUNTER_DEFINED
#define FLOPS_COUNTER_DEFINED

#include <util/Tensor.hpp>

namespace sisi4s {
  class FlopsCounter {
  public:
    FlopsCounter(int64_t *flops, MPI_Comm comm=MPI_COMM_SELF);
    ~FlopsCounter();
  protected:
    int64_t *flops;
    MPI_Comm comm;
    CTF::Flop_counter counter;
  };
}

#endif

