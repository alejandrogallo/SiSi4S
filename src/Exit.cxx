#include <cstdlib>
#include <mpi.h>

namespace sisi4s {

void exit(int const ret) {
  MPI_Finalize();
  std::exit(ret);
}

} // namespace sisi4s
