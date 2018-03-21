#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/MultiCombinations.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalMp2EnergyFromCoulombIntegrals);

ThermalMp2EnergyFromCoulombIntegrals::ThermalMp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalMp2EnergyFromCoulombIntegrals::~ThermalMp2EnergyFromCoulombIntegrals() {
}

void ThermalMp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  beta = 1/getRealArgument("Temperature");

  // compute \Delta^{ab}_{ij} = eps_a + eps_b - eps_i - eps_j
  Dabij = NEW(Tensor<>, false, *Vabij);
  (*Dabij)["abij"] =  (*epsa)["a"];
  (*Dabij)["abij"] += (*epsa)["b"];
  (*Dabij)["abij"] -= (*epsi)["i"];
  (*Dabij)["abij"] -= (*epsi)["j"];

/* tests:
  LOG(1, "FT-MP2") <<
    MultiCombinations(5,3).multinomial(std::vector<unsigned int>({2,0,1,0,0}))
    << std::endl;
  testDerivativeLogZMp2(0);
  testDerivativeLogZMp2(1);
  testDerivativeLogZMp2(2);
  testDerivativeLogZMp2(3);
  testDerivativeLogZHf(0);
  testDerivativeLogZHf(1);
  testDerivativeLogZHf(2);
  testDerivativeLogZHf(3);
  testDerivativeLogZH0(0);
  testDerivativeLogZH0(1);
  testDerivativeLogZH0(2);
*/

  computeFreeEnergy();
  computeEnergyMoments();
}

void ThermalMp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("ThermalPPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}


void ThermalMp2EnergyFromCoulombIntegrals::computeFreeEnergy() {
  real F( -getRealArgument("Temperature") * getDerivativeLogZ(0) );
  setRealArgument("ThermalFreeEnergy", F);
}

void ThermalMp2EnergyFromCoulombIntegrals::computeEnergyMoments() {
  unsigned int n( std::max(getIntegerArgument("maxEnergyMoment", 2),1l) );
  std::vector<real> energyMoments(n);

  for (unsigned int k(1); k <= n; ++k) {
    energyMoments[k-1] = getDerivativeLogZ(k);
  }

  std::vector<int64_t> indices(Dabij->wrld->rank == 0 ? n : 0);
  for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; };
  auto ctfEnergyMoments(new CTF::Tensor<>(1, std::vector<int>({int(n)}).data()));
  ctfEnergyMoments->write(indices.size(), indices.data(), energyMoments.data());
  allocatedTensorArgument("ThermalEnergyMoments", ctfEnergyMoments);
}

real ThermalMp2EnergyFromCoulombIntegrals::getDerivativeLogZ(
  const unsigned int n
) {
  real derivativeLogZ;
  derivativeLogZ = getDerivativeLogZMp2(n);
  derivativeLogZ += getDerivativeLogZHf(n);
  derivativeLogZ += getDerivativeLogZH0(n);
/*
  real mu(getRealArgument("ChemicalPotential"));
  int N(getIntegerArgument("Electrons"));
  switch (n) {
  case 0:
    writeContribution("mu*N", n, -beta*mu*N);
    derivativeLogZ += -beta*mu*N;
    break;
  case 1:
    writeContribution("mu*N", n, mu*N);
    derivativeLogZ += mu*N;
    break;
  }
*/
  return derivativeLogZ;
}

/**
 * \brief Computes the nth derivative of the MP2 contribution to log(Z)
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * w.r.t. (-beta)
 **/
real ThermalMp2EnergyFromCoulombIntegrals::getDerivativeLogZMp2(
  const unsigned int n
) {
  Tensor<> Tabij(false, *getTensorArgument("ThermalPPHHCoulombIntegrals"));
  if (n == 0) {
    // all derivative degrees are zero
    addLogZMp2Amplitudes(Tabij, std::vector<unsigned int>(5));
  } else {
    LOG(2, "FT-MP2") << "computing MP2 contribution of d^"
      << n << " log(Z(beta)) / d(-beta)^" << n << " ..." << std::endl;
    // go through all possiblities to distribute n derivatives among 5 terms
    for (auto &derivedTerms: MultiCombinations(5, n)) {
      // convert to degree of derivative for each term
      std::vector<unsigned int> degrees(5);
      for (auto term: derivedTerms) {
        ++degrees[term];
      };
      std::stringstream contribution;
      real multiplicity( MultiCombinations(5, n).multinomial(degrees) );
      addLogZMp2Amplitudes(Tabij, degrees, multiplicity);
    }
  }

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  energy[""] = (+0.5) * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real direct( energy.get_val() );
  energy[""] = (-0.5) * spins * Tabij["abij"] * (*Vabij)["abji"];
  real exchange( energy.get_val() );
  writeContribution("MP2 direct", n, direct);
  writeContribution("MP2 exchange", n, exchange);
  return direct+exchange;
}

/**
 * \brief Adds the contribution of the thermal MP2 amplitudes
 * \f$ \int_{\tau_1}^\beta{\rm d}\tau_2\int_0^\beta{\rm d}\tau_1\,V^{ab}_{ij}\,
 * f^a\, F^b\, f_i\, f_j\,{\rm e}^{-\beta\Delta^{ab}_{ij}}} \f$
 * to the tensor Tabij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 5 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addLogZMp2Amplitudes(
  Tensor<> &Tabij,
  const std::vector<unsigned int> &degrees, const real multiplicity
) {
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  // start with Vabij
  Tensor<> tabij(*Vabij);
  // apply the derivative of the 5 terms of the respective degrees:
  // Tabij *= f^a = 1/(1+exp(-eps_a*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[0])
    )
  ) (
    (*epsa)["a"], tabij["abij"]
  );
  // Tabij *= f^b
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, true, degrees[1])
    )
  ) (
    (*epsa)["b"], tabij["abij"]
  );
  // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[2])
    )
  ) (
    (*epsi)["i"], tabij["abij"]
  );
  // Tabij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[3])
    )
  ) (
    (*epsi)["j"], tabij["abij"]
  );
  // Tabij *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalMp2Propagation<>(beta, degrees[4])
    )
  ) (
    (*Dabij)["abij"], tabij["abij"]
  );
  // add this contribution
  Tabij["abij"] += multiplicity * tabij["abij"];
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ (-1) \int_0^\beta{\rm d}\tau_1\,f_i\, f_j \f$
 * and returns the result in the tensor Tij.
 **/
real ThermalMp2EnergyFromCoulombIntegrals::getDerivativeLogZHf(
  const unsigned int n
) {
  Tensor<> Tij(2, getTensorArgument("ThermalHHHHCoulombIntegrals")->lens);
  if (n == 0) {
    // all derivative degrees are zero
    addLogZHfAmplitudes(Tij, std::vector<unsigned int>(3));
  } else {
    LOG(2, "FT-MP2") << "computing HF contribution of d^"
      << n << " log(Z(beta)) / d(-beta)^" << n << " ..." << std::endl;
    // go through all possiblities to distribute n derivatives among 3 terms
    for (auto &derivedTerms: MultiCombinations(3, n)) {
      // convert to degree of derivative for each term
      std::vector<unsigned int> degrees(3);
      for (auto term: derivedTerms) {
        ++degrees[term];
      };
      real multiplicity( MultiCombinations(3, n).multinomial(degrees) );
      addLogZHfAmplitudes(Tij, degrees, multiplicity);
    }
  }

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  // we compute -Veff + HF = -2*HF + HF = -HF, so the sign is negative
  energy[""] = (-1.0) * (+0.5) * spins * spins * Tij["ij"] * (*Vijkl)["ijij"];
  real direct( energy.get_val() );
  energy[""] = (-1.0) * (-0.5) * spins * Tij["ij"] * (*Vijkl)["ijji"];
  real exchange( energy.get_val() );

  writeContribution("Hartree", n, direct);
  writeContribution("Exchange", n, exchange);
  return direct+exchange;
}

/**
 * \brief Computes the nth derivative of the thermal HF amplitudes
 * \f$ \int_0^\beta{\rm d}\tau\, f_i\, f_j \f$
 * to the tensor Tij where the derivative degrees w.r.t. \f$(-\beta)$\f of
 * each of the 3 \f$\beta$\f dependent terms is specified by the
 * argument vector.
 **/
void ThermalMp2EnergyFromCoulombIntegrals::addLogZHfAmplitudes(
  Tensor<> &Tij,
  const std::vector<unsigned int> &degrees, const real multiplicity
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> tij(false, Tij);
  // start with 1

  // the first term is (-1) int_0^beta dtau = -beta
  switch (degrees[0]) {
  case 0:
    tij["ij"] = -beta;
    break;
  case 1:
    tij["ij"] = 1;
    break;
  default:
    // all higher derivatives of -beta w.r.t. (-beta) are zero
    // add nothing to Tabij and return
    return;
  }
  // apply the derivative of the other 2 terms of the respective degrees:
  // Tij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[1])
    )
  ) (
    (*epsi)["i"], tij["ij"]
  );
  // Tij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, false, degrees[2])
    )
  ) (
    (*epsi)["j"], tij["ij"]
  );
  // add this contribution
  Tij["ij"] += multiplicity * tij["ij"];
}

real ThermalMp2EnergyFromCoulombIntegrals::getDerivativeLogZH0(
  const unsigned int n
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> Ti(false, *epsi);
  // TODO: H_0 contribution: log(Z_0(beta)) = sum_1^Np log(1+exp(-beta*eps_p))
  switch (n) {
  case 0:
    Transform<real, real>(
      std::function<void(real, real &)>(
        [this](real eps, real &T) {
          T = std::log(1 + std::exp(-beta*eps));
        }
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    break;
  default:
    // the first derivative is eps_i*f_i,
    // use the ThermalContraction clas providing arbitrary derivatives of f_i
    Ti["i"] = (*epsi)["i"];
    Transform<real, real>(
      std::function<void(real, real &)>(
        ThermalContraction<>(beta, false, n-1)
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    break;
  }
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> derivativeLogZH0;
  // sum over all relevant (thermal hole) states
  derivativeLogZH0[""] = spins * Ti["i"];
  writeContribution("single body reference", n, derivativeLogZH0.get_val());
  return derivativeLogZH0.get_val();
}

void ThermalMp2EnergyFromCoulombIntegrals::writeContribution(
  const std::string &contribution,
  const unsigned int n, const real derivativeLogZ
) {
  std::stringstream term;
  real m(derivativeLogZ);
  switch (n) {
  case 0:
    term << "free energy F from " << contribution;
    // convert log(Z) to free energy for writing into log
    m *= -getRealArgument("Temperature");
    break;
  default:
    term << n << ".central moment of H from " << contribution;
    break;
  }
  LOG(0, "FT-MP2") << term.str() << "=" << m << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDerivativeLogZMp2(
  const unsigned int n
) {
  real b(beta), db(0.00001);
  beta = b+db;
  real dLogZ( getDerivativeLogZMp2(n) );
  beta = b-db;
  dLogZ -= getDerivativeLogZMp2(n);
  beta = b;
  dLogZ /= -2*db;
  real exactDLogZ( getDerivativeLogZMp2(n+1) );
  LOG(1, "FT-MP2") << "MP2 numerical(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "MP2     exact(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << exactDLogZ << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDerivativeLogZHf(
  const unsigned int n
) {
  real b(beta), db(0.00001);
  beta = b+db;
  real dLogZ( getDerivativeLogZHf(n) );
  beta = b-db;
  dLogZ -= getDerivativeLogZHf(n);
  beta = b;
  dLogZ /= -2*db;
  real exactDLogZ( getDerivativeLogZHf(n+1) );
  LOG(1, "FT-MP2") << "HF numerical(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "HF     exact(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << exactDLogZ << std::endl;
}

void ThermalMp2EnergyFromCoulombIntegrals::testDerivativeLogZH0(
  const unsigned int n
) {
  real b(beta), db(0.00001);
  beta = b+db;
  real dLogZ( getDerivativeLogZH0(n) );
  beta = b-db;
  dLogZ -= getDerivativeLogZH0(n);
  beta = b;
  dLogZ /= -2*db;
  real exactDLogZ( getDerivativeLogZH0(n+1) );
  LOG(1, "FT-MP2") << "H0 numerical(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << dLogZ << std::endl;
  LOG(1, "FT-MP2") << "H0     exact(d^" << (n+1) << " log(Z(beta)) / d^"
    << (n+1) << ") = " << exactDLogZ << std::endl;
}

