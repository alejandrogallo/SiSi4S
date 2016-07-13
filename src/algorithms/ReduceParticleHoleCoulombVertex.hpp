/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REDUCE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED
#define REDUCE_PARTICLE_HOLE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ReduceParticleHoleCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ReduceParticleHoleCoulombVertex);
    ReduceParticleHoleCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~ReduceParticleHoleCoulombVertex();
    /**
     * \brief calculates the eigenvalue decomposition of the given
     * energy matrix \f$E_G^{G'}=U_G^H \lambda_H {U^\ast}_H^{G'}\f$,
     * truncates \f$H\f$ to a minimal set of indices to reproduce the
     * the energy \f${\rm Tr}\{E\}\approx\sum_H\lambda_H\f$ within
     * the required accuracy specified.
     * Reduces the particle hole Coulomb vertex by applying the
     * energy matrix reduction transform:
     * \f$\Gamma^{qg}_r = \Gamma^{qG}_r U_G^g\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of reducing the Coulomb vertex.
     */
    virtual void dryRun();
  };
}

#endif

