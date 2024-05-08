#include <unistd.h>

#include <Sisi4s.hpp>

#include <equations/SimilarityTransformedHamiltonian.hpp>

#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <algorithms/UccsdtqAmplitudesFromCoulombIntegrals.hpp>

#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>

#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

using F = double;
using FSPEC = double;
DEFSPEC(UccsdtqAmplitudesFromCoulombIntegrals,
        SPEC_IN(UCCSD_SPEC_IN,
                {"InitialTriplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)},
                {"InitialQuadruplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)},
                {"hirataEquations", SPEC_VALUE_DEF("TODO: DOC", bool, false)}),
        SPEC_OUT(UCCSD_SPEC_OUT,
                 {"TriplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)},
                 {"QuadruplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)}));

IMPLEMENT_ALGORITHM(UccsdtqAmplitudesFromCoulombIntegrals) {
  in.set<bool>("unrestricted", true);
  in.set<bool>("antisymmetrize", true);
  ClusterSinglesDoublesTriplesQuadruplesAlgorithm::run();
}

PTR(FockVector<double>) UccsdtqAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

PTR(FockVector<sisi4s::complex>)
UccsdtqAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<sisi4s::complex>) &amplitudes) {
  return getResiduumTemplate<sisi4s::complex>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdtqAmplitudesFromCoulombIntegrals::getResiduumTemplate(
    const int iterationStep,
    const PTR(const FockVector<F>) &amplitudes) {
  // Equations from: hirata group
  // https://github.com/alejandrogallo/hirata

  Tensor<double> *epsi(in.get<Tensor<double> *>("HoleEigenEnergies"));

  Tensor<double> *epsa(in.get<Tensor<double> *>("ParticleEigenEnergies"));

  // Get couloumb integrals
  auto Vijkl(in.get<Tensor<F> *>("HHHHCoulombIntegrals"));
  auto Vabcd(in.get<Tensor<F> *>("PPPPCoulombIntegrals"));
  auto Vijka(in.get<Tensor<F> *>("HHHPCoulombIntegrals"));
  auto Vijab(in.get<Tensor<F> *>("HHPPCoulombIntegrals"));
  auto Viajk(in.get<Tensor<F> *>("HPHHCoulombIntegrals"));
  auto Viajb(in.get<Tensor<F> *>("HPHPCoulombIntegrals"));
  auto Viabc(in.get<Tensor<F> *>("HPPPCoulombIntegrals"));
  auto Vabij(in.get<Tensor<F> *>("PPHHCoulombIntegrals"));
  auto Vabic(in.get<Tensor<F> *>("PPHPCoulombIntegrals"));
  // auto Viabj(in.get<Tensor<F>*>("HPPHCoulombIntegrals"));
  // auto Vaibc(in.get<Tensor<F>*>("PHPPCoulombIntegrals"));
  // auto Vijak(in.get<Tensor<F>*>("HHPHCoulombIntegrals"));
  // auto Vabci(in.get<Tensor<F>*>("PPPHCoulombIntegrals"));

  int Nv = (epsa->lens[0]), //
      No = (epsi->lens[0]), //
      vv[] = {Nv, Nv},      //
      oo[] = {No, No},      //
      syms[] = {NS, NS};
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab")),
      *Fij(new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij")), //
      *Fia;

  if (in.present("HPFockMatrix") && in.present("HHFockMatrix")
      && in.present("PPFockMatrix")) {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    Fia = in.get<Tensor<F> *>("HPFockMatrix");
    Fab = in.get<Tensor<F> *>("PPFockMatrix");
    Fij = in.get<Tensor<F> *>("HHFockMatrix");
  } else {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using hartree fock orbitals" << std::endl;
      LOG(0, getAbbreviation()) << "whatever" << std::endl;
    }
    Fia = NULL;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

      LOG(0, getAbbreviation()) << "whatever2" << std::endl;

  // Create T and R and intermediates
  //
  // Read the amplitudes Tai, Tabij and Tabcijk
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");
  auto Tabcijk(amplitudes->get(2));
  Tabcijk->set_name("Tabcijk");
  auto Tabcdijkl(amplitudes->get(3));
  Tabcdijkl->set_name("Tabcdijkl");

  LOG(0, getAbbreviation()) << "whatever3" << std::endl;

  LOG(0, getAbbreviation()) << "residduum" << std::endl;
  auto residuum(NEW(FockVector<F>, *amplitudes));
  LOG(0, getAbbreviation()) << "before zeroing" << std::endl;
  *residuum *= 0.0;
  LOG(0, getAbbreviation()) << "zeroing" << std::endl;
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  LOG(0, getAbbreviation()) << "Rai" << std::endl;
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");
  LOG(0, getAbbreviation()) << "Rabij" << std::endl;
  auto Rabcijk(residuum->get(2));
  Rabcijk->set_name("Rabcijk");
  LOG(0, getAbbreviation()) << "Rabcijk" << std::endl;
  auto Rabcdijkl(residuum->get(3));
  Rabcdijkl->set_name("Rabcdijkl");
  LOG(0, getAbbreviation()) << "Rabcijk" << std::endl;

  LOG(0, getAbbreviation()) << "residuum" << std::endl;

  if (iterationStep == 0) {
    LOG(0, getAbbreviation()) << "Starting first iteration" << std::endl;
    LOG(0, getAbbreviation()) << "Hirata equations " << in.get<bool>("hirataEquations") << std::endl;
  }


  if (in.get<bool>("hirataEquations")) {
    LOG(1, getAbbreviation())
      << "Using Hirata equations without any kind of intermediates" << std::endl;
#include <equations/ccsdtq_hirata>
  } else {

    LOG(1, getAbbreviation()) << "Using intermediates" << std::endl;
    SimilarityTransformedHamiltonian<F> H(Fij->lens[0], Fab->lens[0]);

    H
        // fock matrices
        .setFij(Fij)
        .setFab(Fab)
        .setFia(Fia)
        // Coulomb Integrals
        .setVabcd(Vabcd)
        .setViajb(Viajb)
        .setVijab(Vijab)
        .setVijkl(Vijkl)
        .setVijka(Vijka)
        .setViabc(Viabc)
        .setViajk(Viajk)
        .setVabic(Vabic)
        // .setVaibc(Vaibc)
        // .setVaibj(Vaibj)
        // .setViabj(Viabj)
        // .setVijak(Vijak)
        // .setVaijb(Vaijb)
        // .setVabci(Vabci)
        .setVabij(Vabij)
        // CC amplitudes
        .setTai(Tai.get())           // singles
        .setTabij(Tabij.get())       // doubles
        .setTabcijk(Tabcijk.get())   // triples
        .setTabcdijkl(Tabcijk.get()) // quadruples
        .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL)
        .withStantonIntermediatesUCCSD(true)
        // end
        ;

    // T1 equations:
    LOG(1, getAbbreviation()) << "<AI|H|0>" << std::endl;
    auto Wai = H.getAI();
    (*Rai)["ai"] += (*Wai)["ai"];
    // These are the residum equations, we have to substract them from Wai
    (*Rai)["bi"] += (-1.0) * (*Fab)["bc"] * (*Tai)["ci"];
    (*Rai)["bi"] += (+1.0) * (*Fij)["ki"] * (*Tai)["bk"];

    // T2 equations:
    LOG(1, getAbbreviation()) << "<ABIJ|H|0>" << std::endl;
    auto Wabij = H.getABIJ();
    (*Rabij)["abij"] += (*Wabij)["abij"];
    // These are the residum equations, substract them from Wabij
    (*Rabij)["cdij"] += (+1.0) * (*Fij)["ii"] * (*Tabij)["cdij"];
    (*Rabij)["cdij"] += (-1.0) * (*Fij)["jj"] * (*Tabij)["cdji"];
    (*Rabij)["cdij"] += (+1.0) * (*Fab)["dd"] * (*Tabij)["dcij"];
    (*Rabij)["cdij"] += (-1.0) * (*Fab)["cc"] * (*Tabij)["cdij"];

    // T3 equations:
    LOG(1, getAbbreviation()) << "<ABCIJK|H|0>" << std::endl;
    auto Wabcijk = H.getABCIJK();
    (*Rabcijk)["abcijk"] += (*Wabcijk)["abcijk"];
    // These are the residum equations, substract them from Wabij
    (*Rabcijk)["defijk"] += (+1.0) * (*Fij)["ii"] * (*Tabcijk)["defijk"];
    (*Rabcijk)["defijk"] += (-1.0) * (*Fij)["jj"] * (*Tabcijk)["defjik"];
    (*Rabcijk)["defijk"] += (-1.0) * (*Fij)["kk"] * (*Tabcijk)["defkji"];
    (*Rabcijk)["defijk"] += (-1.0) * (*Fab)["ff"] * (*Tabcijk)["fdeijk"];
    (*Rabcijk)["defijk"] += (+1.0) * (*Fab)["ee"] * (*Tabcijk)["edfijk"];
    (*Rabcijk)["defijk"] += (+1.0) * (*Fab)["dd"] * (*Tabcijk)["dfeijk"];

    LOG(1, getAbbreviation()) << "<ABCDIJKL|H|0>" << std::endl;
    auto Wabcdijkl = H.getABCDIJKL();
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fij)["Ii"] * (*Tabcdijkl)["efghIjkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fij)["Ij"] * (*Tabcdijkl)["efghIikl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fij)["Ik"] * (*Tabcdijkl)["efghIjil"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fij)["Il"] * (*Tabcdijkl)["efghIjki"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fab)["hA"] * (*Tabcdijkl)["Aefgijkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fab)["gA"] * (*Tabcdijkl)["Aefhijkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fab)["fA"] * (*Tabcdijkl)["Aehgijkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fab)["eA"] * (*Tabcdijkl)["Ahfgijkl"];
  }

  return residuum;
}
