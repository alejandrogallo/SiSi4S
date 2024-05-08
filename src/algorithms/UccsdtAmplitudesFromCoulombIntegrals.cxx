#include <unistd.h>

#include <Sisi4s.hpp>

#include <algorithms/UccsdtAmplitudesFromCoulombIntegrals.hpp>
#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>

#include <equations/SimilarityTransformedHamiltonian.hpp>

#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>

#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <util/Tensor.hpp>

using namespace sisi4s;

#ifdef DEBUG
#  define LDEBUG(msg)                                                          \
    LOG(1, getAbbreviation()) << __LINE__ << ":" << msg << std::endl;
#else
#  define LDEBUG(msg)
#endif

using F = double;
using FSPEC = double;
DEFSPEC(
    UccsdtAmplitudesFromCoulombIntegrals,
    SPEC_IN(
        UCCSD_SPEC_IN,
        {"NonAntiHHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
        {"NonAntiHHPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
        {"NonAntiHPPPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
        {"NonAntiPHHPCoulombIntegrals", SPEC_VARIN("TODO: DOC", Tensor<F> *)},
        {"initialTriplesAmplitudes", SPEC_VARIN("TODO: DOC", Tensor<F> *)}, \
        {"hirataEquations", SPEC_VALUE_DEF("TODO: DOC", bool, false)},
        {"withRingCCSDT", SPEC_VALUE_DEF("TODO: DOC", bool, false)}),
    SPEC_OUT(UCCSD_SPEC_OUT,
             {"TriplesAmplitudes", SPEC_VAROUT("TODO: DOC", Tensor<F> *)}));

IMPLEMENT_ALGORITHM(UccsdtAmplitudesFromCoulombIntegrals) {
  in.set<bool>("unrestricted", true);
  in.set<bool>("antisymmetrize", true);
  ClusterSinglesDoublesTriplesAlgorithm::run();
}

PTR(FockVector<double>) UccsdtAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return in.get<bool>("hirataEquations")
           ? getResiduumTemplate<double>(iterationStep, amplitudes)
           : getResiduumSth<double>(iterationStep, amplitudes);
}

PTR(FockVector<sisi4s::complex>)
UccsdtAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<sisi4s::complex>) &amplitudes) {
  return in.get<bool>("hirataEquations")
           ? getResiduumTemplate<sisi4s::complex>(iterationStep, amplitudes)
           : getResiduumSth<sisi4s::complex>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdtAmplitudesFromCoulombIntegrals::getResiduumSth(
    const int iterationStep,
    const PTR(const FockVector<F>) &amplitudes) {

  auto *epsi = in.get<Tensor<double> *>("HoleEigenEnergies"),
       *epsa = in.get<Tensor<double> *>("ParticleEigenEnergies");

  const bool usingIntermediates = in.get<bool>("intermediates"),
             withRingCCSDT = in.get<bool>("withRingCCSDT");

  if (!usingIntermediates) {
    LOG(0, getAbbreviation())
        << "Not using CCSD stanton intermediates, the code will be much slower."
        << std::endl;
  }

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
  auto Viabj(in.get<Tensor<F> *>("HPPHCoulombIntegrals"));
  auto Vaibc(in.get<Tensor<F> *>("PHPPCoulombIntegrals"));
  auto Vijak(in.get<Tensor<F> *>("HHPHCoulombIntegrals"));
  auto Vabci(in.get<Tensor<F> *>("PPPHCoulombIntegrals"));
  auto Vaibj(in.get<Tensor<F> *>("PHPHCoulombIntegrals"));
  auto Vaijb(in.get<Tensor<F> *>("PHHPCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  Tensor<F> *Fab, *Fij, *Fia;

  if (in.present("HPFockMatrix") && in.present("HHFockMatrix")
      && in.present("PPFockMatrix")) {
    if (iterationStep == 0) {
      LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    Fia = in.get<Tensor<F> *>("HPFockMatrix");
    Fab = in.get<Tensor<F> *>("PPFockMatrix");
    Fij = in.get<Tensor<F> *>("HHFockMatrix");
  } else {
    Fia = nullptr;
    Fab = new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab");
    Fij = new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij");
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");
  auto Tabcijk(amplitudes->get(2));
  Tabcijk->set_name("Tabcijk");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");
  auto Rabcijk(residuum->get(2));
  Rabcijk->set_name("Rabcijk");

  if ((iterationStep == 0) && !in.present("initialDoublesAmplitudes")) {
    LOG(1, getAbbreviation())
        << "Set initial Rabij amplitudes to Vijab" << std::endl;
    (*Rabij)["abij"] = (*Vijab)["ijab"];
    return residuum;
  }

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
      .setVaibc(Vaibc)
      .setVaibj(Vaibj)
      .setViabj(Viabj)
      .setVijak(Vijak)
      .setVaijb(Vaijb)
      .setVabci(Vabci)
      .setVabij(Vabij)
      // CC amplitudes
      .setTai(Tai.get())
      .setTabij(Tabij.get())
      .setTabcijk(Tabcijk.get())
      .setDressing(SimilarityTransformedHamiltonian<F>::Dressing::GENERAL)
      .withStantonIntermediatesUCCSD(usingIntermediates)
      .withRingCCSDT(withRingCCSDT)
      // end
      ;

  if (H.withRingCCSDT()) {
    LOG(1, getAbbreviation()) << "computing RING-CCSDT" << std::endl;
    H
        // set non antisymmetric integrals
        .setVVaijb(in.get<Tensor<F> *>("NonAntiPHHPCoulombIntegrals"))
        .setVVijka(in.get<Tensor<F> *>("NonAntiHHHPCoulombIntegrals"))
        .setVViabc(in.get<Tensor<F> *>("NonAntiHPPPCoulombIntegrals"))
        .setVVijab(in.get<Tensor<F> *>("NonAntiHHPPCoulombIntegrals"));
  }

  // T1 equations:
  auto Wai = H.getAI();
  (*Rai)["ai"] += (*Wai)["ai"];
  // These are the residum equations, we have to substract them from Wai
  (*Rai)["bi"] += (-1.0) * (*Fab)["bc"] * (*Tai)["ci"];
  (*Rai)["bi"] += (+1.0) * (*Fij)["ki"] * (*Tai)["bk"];

  // T2 equations:
  auto Wabij = H.getABIJ();
  (*Rabij)["abij"] += (*Wabij)["abij"];
  // These are the residum equations, substract them from Wabij
  (*Rabij)["cdij"] += (+1.0) * (*Fij)["ii"] * (*Tabij)["cdij"];
  (*Rabij)["cdij"] += (-1.0) * (*Fij)["jj"] * (*Tabij)["cdji"];
  (*Rabij)["cdij"] += (+1.0) * (*Fab)["dd"] * (*Tabij)["dcij"];
  (*Rabij)["cdij"] += (-1.0) * (*Fab)["cc"] * (*Tabij)["cdij"];

  // T3 equations:
  auto Wabcijk = H.getABCIJK();
  (*Rabcijk)["abcijk"] += (*Wabcijk)["abcijk"];
  // These are the residum equations, substract them from Wabij
  (*Rabcijk)["defijk"] += (+1.0) * (*Fij)["ii"] * (*Tabcijk)["defijk"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Fij)["jj"] * (*Tabcijk)["defjik"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Fij)["kk"] * (*Tabcijk)["defkji"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Fab)["ff"] * (*Tabcijk)["fdeijk"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Fab)["ee"] * (*Tabcijk)["edfijk"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Fab)["dd"] * (*Tabcijk)["dfeijk"];

  // TODO: delete Fij and so on
  return residuum;
}

// TODO: remove this when tested against STH
template <typename F>
PTR(FockVector<F>) UccsdtAmplitudesFromCoulombIntegrals::getResiduumTemplate(
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

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  Tensor<F> *Fab(new Tensor<F>(2, vv, syms, *Sisi4s::world, "Fab"));
  Tensor<F> *Fij(new Tensor<F>(2, oo, syms, *Sisi4s::world, "Fij"));
  Tensor<F> *Fia;

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
    }
    Fia = nullptr;
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsi)["i"], (*Fij)["ii"]);
    CTF::Transform<double, F>(std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }))((*epsa)["a"], (*Fab)["aa"]);
  }

  // Create T and R and intermediates
  //
  // Read the amplitudes Tai, Tabij and Tabcijk
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");
  auto Tabcijk(amplitudes->get(2));
  Tabcijk->set_name("Tabcijk");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");
  auto Rabcijk(residuum->get(2));
  Rabcijk->set_name("Rabcijk");

  // Singles
  (*Rai)["bi"] = 0.0;

  if (Fia) {
    (*Rai)["bi"] += (+1.0) * (*Fia)["ib"];
    (*Rai)["bi"] += (+1.0) * (*Fia)["kd"] * (*Tabij)["dbki"];
    (*Rai)["bi"] += (-1.0) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
  }

  // Residum equations
  //(*Rai)["bi"] += ( - 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];
  //(*Rai)["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];

  (*Rai)["bi"] += (-1.0) * (*Tai)["cl"] * (*Viajb)["lbic"];
  (*Rai)["bi"] += (+0.5) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
  (*Rai)["bi"] += (+0.5) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
  (*Rai)["bi"] += (+0.25) * (*Tabcijk)["cdbmni"] * (*Vijab)["mncd"];
  (*Rai)["bi"] += (-1.0) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
  (*Rai)["bi"] += (-1.0) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
  (*Rai)["bi"] += (-0.5) * (*Tabij)["cblm"] * (*Tai)["fi"] * (*Vijab)["lmcf"];
  (*Rai)["bi"] += (-0.5) * (*Tabij)["cdmi"] * (*Tai)["bn"] * (*Vijab)["mncd"];
  (*Rai)["bi"] += (+1.0) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
  (*Rai)["bi"] +=
      (-1.0) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Singles done" << std::endl;
#endif

  // Doubles
  (*Rabij)["cdij"] = 0.0;

  if (Fia) {
    LDEBUG("Fia * T2 * T1")
    (*Rabij)["cdij"] += (-1.0) * (*Fia)["mf"] * (*Tabij)["cdmj"] * (*Tai)["fi"];
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabij)["cdmi"] * (*Tai)["fj"];
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabij)["fcij"] * (*Tai)["dm"];
    (*Rabij)["cdij"] += (-1.0) * (*Fia)["mf"] * (*Tabij)["fdij"] * (*Tai)["cm"];
    LDEBUG("Fia * T3")
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabcijk)["fcdmij"];
  }

  // Residum equations
  //(*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  //(*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  //(*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  //(*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];

  (*Rabij)["cdij"] += (+1.0) * (*Vabij)["cdij"];
  LDEBUG("T1 * V")
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["cm"] * (*Viajk)["mdij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["dm"] * (*Viajk)["mcij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ej"] * (*Vabic)["cdie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ei"] * (*Vabic)["cdje"];
  LDEBUG("T2 * V")
  (*Rabij)["cdij"] += (+0.5) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["edni"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
  LDEBUG("T3 * V")
  (*Rabij)["cdij"] += (+0.5) * (*Tabcijk)["ecdnoj"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += (-0.5) * (*Tabcijk)["ecdnoi"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += (-0.5) * (*Tabcijk)["efcoij"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabcijk)["efdoij"] * (*Viabc)["ocef"];
  LDEBUG("T1 * T1 * V")
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
  LDEBUG("T2 * T1 * V")
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];
  LDEBUG("T3 * T1 * V")
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabcijk)["ecdnoj"] * (*Tai)["hi"] * (*Vijab)["noeh"];
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabcijk)["ecdnoi"] * (*Tai)["hj"] * (*Vijab)["noeh"];
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabcijk)["efcoij"] * (*Tai)["dp"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabcijk)["efdoij"] * (*Tai)["cp"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabcijk)["ecdnij"] * (*Tai)["gp"] * (*Vijab)["npeg"];
  LDEBUG("T2 * T2 * V")
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] +=
      (+0.25) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] +=
      (-0.5) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
  (*Rabij)["cdij"] +=
      (+0.5) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
  LDEBUG("T1 * T1 * T1 * V")
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];
  LDEBUG("T1 * T1 * T2 * V")
  (*Rabij)["cdij"] += (+0.5) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"]
                    * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["dn"] * (*Tabij)["gcpj"]
                    * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tabij)["gdpj"]
                    * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ej"] * (*Tai)["dn"] * (*Tabij)["gcpi"]
                    * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tabij)["gdpi"]
                    * (*Vijab)["npeg"];
  (*Rabij)["cdij"] += (+0.5) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"]
                    * (*Vijab)["mngh"];

  LDEBUG("T2 * T1 * T1 * V")
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["cdpj"] * (*Tai)["ei"] * (*Tai)["fo"]
                    * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["cdpi"] * (*Tai)["ej"] * (*Tai)["fo"]
                    * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["hcij"] * (*Tai)["dm"] * (*Tai)["fo"]
                    * (*Vijab)["mofh"];
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["hdij"] * (*Tai)["cm"] * (*Tai)["fo"]
                    * (*Vijab)["mofh"];

  LDEBUG("T1 * T1 * T1 * T1 * V")
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"]
                    * (*Tai)["dp"] * (*Vijab)["opef"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Doubles done" << std::endl;
#endif

  // Triples
  (*Rabcijk)["edfijk"] = 0.0;

  if (Fia) {
    LDEBUG("Fia * T3 * T1")
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabcijk)["defojk"] * (*Fia)["oh"] * (*Tai)["hi"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabcijk)["defoki"] * (*Fia)["oh"] * (*Tai)["hj"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabcijk)["defoij"] * (*Fia)["oh"] * (*Tai)["hk"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabcijk)["hdeijk"] * (*Fia)["oh"] * (*Tai)["fo"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabcijk)["hdfijk"] * (*Fia)["oh"] * (*Tai)["eo"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabcijk)["hfeijk"] * (*Fia)["oh"] * (*Tai)["do"];
    LDEBUG("Fia * T2 * T2")
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabij)["gfij"] * (*Tabij)["depk"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabij)["geij"] * (*Tabij)["dfpk"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabij)["gdij"] * (*Tabij)["fepk"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabij)["gfik"] * (*Tabij)["depj"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabij)["geik"] * (*Tabij)["dfpj"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabij)["gdik"] * (*Tabij)["fepj"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabij)["gdkj"] * (*Tabij)["fepi"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Tabij)["gekj"] * (*Tabij)["dfpi"] * (*Fia)["pg"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Tabij)["gfkj"] * (*Tabij)["depi"] * (*Fia)["pg"];
  }

  LDEBUG("T2 * V")
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deok"] * (*Viajk)["ofij"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfok"] * (*Viajk)["oeij"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feok"] * (*Viajk)["odij"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deoj"] * (*Viajk)["ofik"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfoj"] * (*Viajk)["oeik"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feoj"] * (*Viajk)["odik"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deoi"] * (*Viajk)["ofkj"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfoi"] * (*Viajk)["oekj"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feoi"] * (*Viajk)["odkj"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdjk"] * (*Vabic)["efig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gejk"] * (*Vabic)["dfig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfjk"] * (*Vabic)["edig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdki"] * (*Vabic)["efjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["geki"] * (*Vabic)["dfjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfki"] * (*Vabic)["edjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdij"] * (*Vabic)["efkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["geij"] * (*Vabic)["dfkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfij"] * (*Vabic)["edkg"];

  // Residum equations
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fij)["oi"] * (*Tabcijk)["defojk"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fij)["oj"] * (*Tabcijk)["defoik"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fij)["ok"] * (*Tabcijk)["defoji"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fab)["fg"] * (*Tabcijk)["gdeijk"];
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fab)["eg"] * (*Tabcijk)["gdfijk"];
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fab)["dg"] * (*Tabcijk)["gfeijk"];

  LDEBUG("T3 * V")
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["defopk"] * (*Vijkl)["opij"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["defopj"] * (*Vijkl)["opik"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["defopi"] * (*Vijkl)["opkj"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["gdepjk"] * (*Viajb)["pfig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gdfpjk"] * (*Viajb)["peig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gfepjk"] * (*Viajb)["pdig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["gdepki"] * (*Viajb)["pfjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gdfpki"] * (*Viajb)["pejg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gfepki"] * (*Viajb)["pdjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["gdepij"] * (*Viajb)["pfkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gdfpij"] * (*Viajb)["pekg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["gfepij"] * (*Viajb)["pdkg"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["ghdijk"] * (*Vabcd)["efgh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["gheijk"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["ghfijk"] * (*Vabcd)["edgh"];

  LDEBUG("T2 * T1 * V")
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deok"] * (*Tai)["fp"] * (*Vijkl)["opij"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["dfok"] * (*Tai)["ep"] * (*Vijkl)["opij"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["feok"] * (*Tai)["dp"] * (*Vijkl)["opij"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoj"] * (*Tai)["fp"] * (*Vijkl)["opik"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["dfoj"] * (*Tai)["ep"] * (*Vijkl)["opik"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["feoj"] * (*Tai)["dp"] * (*Vijkl)["opik"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoi"] * (*Tai)["fp"] * (*Vijkl)["opkj"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["dfoi"] * (*Tai)["ep"] * (*Vijkl)["opkj"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["feoi"] * (*Tai)["dp"] * (*Vijkl)["opkj"];

  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deok"] * (*Tai)["hj"] * (*Viajb)["ofih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["dfok"] * (*Tai)["hj"] * (*Viajb)["oeih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["feok"] * (*Tai)["hj"] * (*Viajb)["odih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deoj"] * (*Tai)["hk"] * (*Viajb)["ofih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["dfoj"] * (*Tai)["hk"] * (*Viajb)["oeih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["feoj"] * (*Tai)["hk"] * (*Viajb)["odih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deok"] * (*Tai)["hi"] * (*Viajb)["ofjh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["dfok"] * (*Tai)["hi"] * (*Viajb)["oejh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["feok"] * (*Tai)["hi"] * (*Viajb)["odjh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoj"] * (*Tai)["hi"] * (*Viajb)["ofkh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["dfoj"] * (*Tai)["hi"] * (*Viajb)["oekh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["feoj"] * (*Tai)["hi"] * (*Viajb)["odkh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoi"] * (*Tai)["hk"] * (*Viajb)["ofjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["dfoi"] * (*Tai)["hk"] * (*Viajb)["oejh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["feoi"] * (*Tai)["hk"] * (*Viajb)["odjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deoi"] * (*Tai)["hj"] * (*Viajb)["ofkh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["dfoi"] * (*Tai)["hj"] * (*Viajb)["oekh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["feoi"] * (*Tai)["hj"] * (*Viajb)["odkh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdjk"] * (*Tai)["ep"] * (*Viajb)["pfig"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gejk"] * (*Tai)["dp"] * (*Viajb)["pfig"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdjk"] * (*Tai)["fp"] * (*Viajb)["peig"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gejk"] * (*Tai)["fp"] * (*Viajb)["pdig"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfjk"] * (*Tai)["dp"] * (*Viajb)["peig"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfjk"] * (*Tai)["ep"] * (*Viajb)["pdig"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdki"] * (*Tai)["ep"] * (*Viajb)["pfjg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geki"] * (*Tai)["dp"] * (*Viajb)["pfjg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdki"] * (*Tai)["fp"] * (*Viajb)["pejg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["geki"] * (*Tai)["fp"] * (*Viajb)["pdjg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfki"] * (*Tai)["dp"] * (*Viajb)["pejg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfki"] * (*Tai)["ep"] * (*Viajb)["pdjg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdij"] * (*Tai)["ep"] * (*Viajb)["pfkg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geij"] * (*Tai)["dp"] * (*Viajb)["pfkg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdij"] * (*Tai)["fp"] * (*Viajb)["pekg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["geij"] * (*Tai)["fp"] * (*Viajb)["pdkg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfij"] * (*Tai)["dp"] * (*Viajb)["pekg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfij"] * (*Tai)["ep"] * (*Viajb)["pdkg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdjk"] * (*Tai)["hi"] * (*Vabcd)["efgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gejk"] * (*Tai)["hi"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfjk"] * (*Tai)["hi"] * (*Vabcd)["edgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdki"] * (*Tai)["hj"] * (*Vabcd)["efgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geki"] * (*Tai)["hj"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfki"] * (*Tai)["hj"] * (*Vabcd)["edgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdij"] * (*Tai)["hk"] * (*Vabcd)["efgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geij"] * (*Tai)["hk"] * (*Vabcd)["dfgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfij"] * (*Tai)["hk"] * (*Vabcd)["edgh"];

  LDEBUG("T3 * T1 * V")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defopk"] * (*Tai)["Aj"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["defopj"] * (*Tai)["Ak"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["defopk"] * (*Tai)["Ai"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defopj"] * (*Tai)["Ai"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defopi"] * (*Tai)["Ak"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["defopi"] * (*Tai)["Aj"] * (*Vijka)["opkA"];

  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepjk"] * (*Tai)["fI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpjk"] * (*Tai)["eI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepjk"] * (*Tai)["dI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepki"] * (*Tai)["fI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpki"] * (*Tai)["eI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepki"] * (*Tai)["dI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepij"] * (*Tai)["fI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpij"] * (*Tai)["eI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepij"] * (*Tai)["dI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["defojk"] * (*Tai)["hI"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["defoki"] * (*Tai)["hI"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["defoij"] * (*Tai)["hI"] * (*Vijka)["oIkh"];

  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepjk"] * (*Tai)["Ai"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpjk"] * (*Tai)["Ai"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepjk"] * (*Tai)["Ai"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepki"] * (*Tai)["Aj"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpki"] * (*Tai)["Aj"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepki"] * (*Tai)["Aj"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepij"] * (*Tai)["Ak"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpij"] * (*Tai)["Ak"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepij"] * (*Tai)["Ak"] * (*Viabc)["pdgA"];

  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["ghdijk"] * (*Tai)["eI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gheijk"] * (*Tai)["dI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["ghdijk"] * (*Tai)["fI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gheijk"] * (*Tai)["fI"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["ghfijk"] * (*Tai)["dI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["ghfijk"] * (*Tai)["eI"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdeijk"] * (*Tai)["hI"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdfijk"] * (*Tai)["hI"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gfeijk"] * (*Tai)["hI"] * (*Viabc)["Idgh"];

  LDEBUG("T2 * T2 * V")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["gfjk"] * (*Tabij)["depI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["gejk"] * (*Tabij)["dfpI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["gdjk"] * (*Tabij)["fepI"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["gfki"] * (*Tabij)["depI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["geki"] * (*Tabij)["dfpI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["gdki"] * (*Tabij)["fepI"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["gfij"] * (*Tabij)["depI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["geij"] * (*Tabij)["dfpI"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["gdij"] * (*Tabij)["fepI"] * (*Vijka)["pIkg"];

  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["efoj"] * (*Tabij)["hdIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["fdoj"] * (*Tabij)["heIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deok"] * (*Tabij)["hfIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deoj"] * (*Tabij)["hfIk"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["fdok"] * (*Tabij)["heIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["efok"] * (*Tabij)["hdIj"] * (*Vijka)["oIih"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["efoi"] * (*Tabij)["hdIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["fdoi"] * (*Tabij)["heIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deok"] * (*Tabij)["hfIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoi"] * (*Tabij)["hfIk"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["fdok"] * (*Tabij)["heIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["efok"] * (*Tabij)["hdIi"] * (*Vijka)["oIjh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["efoi"] * (*Tabij)["hdIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["fdoi"] * (*Tabij)["heIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["deoj"] * (*Tabij)["hfIi"] * (*Vijka)["oIkh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["deoi"] * (*Tabij)["hfIj"] * (*Vijka)["oIkh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["fdoj"] * (*Tabij)["heIi"] * (*Vijka)["oIkh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["efoj"] * (*Tabij)["hdIi"] * (*Vijka)["oIkh"];

  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geij"] * (*Tabij)["hdIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdij"] * (*Tabij)["heIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfij"] * (*Tabij)["hdIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfij"] * (*Tabij)["heIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdij"] * (*Tabij)["hfIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["geij"] * (*Tabij)["hfIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["geik"] * (*Tabij)["hdIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdik"] * (*Tabij)["heIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfik"] * (*Tabij)["hdIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfik"] * (*Tabij)["heIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdik"] * (*Tabij)["hfIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["geik"] * (*Tabij)["hfIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gdkj"] * (*Tabij)["heIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gekj"] * (*Tabij)["hdIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gekj"] * (*Tabij)["hfIi"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabij)["gfkj"] * (*Tabij)["heIi"] * (*Viabc)["Idgh"];

  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["deok"] * (*Tabij)["hAij"] * (*Viabc)["ofhA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["dfok"] * (*Tabij)["hAij"] * (*Viabc)["oehA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["feok"] * (*Tabij)["hAij"] * (*Viabc)["odhA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["deoj"] * (*Tabij)["hAik"] * (*Viabc)["ofhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["dfoj"] * (*Tabij)["hAik"] * (*Viabc)["oehA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["feoj"] * (*Tabij)["hAik"] * (*Viabc)["odhA"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["deoi"] * (*Tabij)["hAkj"] * (*Viabc)["ofhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"] * (*Viabc)["oehA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["feoi"] * (*Tabij)["hAkj"] * (*Viabc)["odhA"];
  LDEBUG("T3 * T2 * V (1)")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["Bfij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIk"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["Beij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIk"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["Bdij"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIk"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["Bfik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIj"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["Beik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIj"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["Bdik"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIj"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["Bfkj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdepIi"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["Bekj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gdfpIi"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["Bdkj"] * (*Vijab)["pIgB"] * (*Tabcijk)["gfepIi"];
  LDEBUG("T3 * T2 * V (2)")
  (*Rabcijk)["defijk"] +=
      (+0.25) * (*Tabcijk)["defopk"] * (*Tabij)["ABij"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["defopj"] * (*Tabij)["ABik"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["defopi"] * (*Tabij)["ABkj"] * (*Vijab)["opAB"];
  LDEBUG("T3 * T2 * V (3)")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["efJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIjk"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["dfJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIjk"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["edJi"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIjk"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["efJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIki"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["dfJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIki"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["edJj"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIki"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabij)["efJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghdIij"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["dfJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["gheIij"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabij)["edJk"] * (*Vijab)["IJgh"] * (*Tabcijk)["ghfIij"];
  LDEBUG("T3 * T2 * V (4)")
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepjk"] * (*Tabij)["AfJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpjk"] * (*Tabij)["AeJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepjk"] * (*Tabij)["AdJi"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepki"] * (*Tabij)["AfJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpki"] * (*Tabij)["AeJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepki"] * (*Tabij)["AdJj"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcijk)["gdepij"] * (*Tabij)["AfJk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gdfpij"] * (*Tabij)["AeJk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] +=
      (-1.0) * (*Tabcijk)["gfepij"] * (*Tabij)["AdJk"] * (*Vijab)["pJgA"];
  LDEBUG("T3 * T2 * V (5)")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defojk"] * (*Tabij)["hAJi"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defoki"] * (*Tabij)["hAJj"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defoij"] * (*Tabij)["hAJk"] * (*Vijab)["oJhA"];
  LDEBUG("T3 * T2 * V (6)")
  (*Rabcijk)["defijk"] +=
      (+0.25) * (*Tabcijk)["ghdijk"] * (*Tabij)["efIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["gheijk"] * (*Tabij)["dfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["ghfijk"] * (*Tabij)["edIJ"] * (*Vijab)["IJgh"];
  LDEBUG("T3 * T2 * V (7)")
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gdeijk"] * (*Tabij)["hfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gdfijk"] * (*Tabij)["heIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gfeijk"] * (*Tabij)["hdIJ"] * (*Vijab)["IJgh"];
  LDEBUG("T1 * T1 * T2 * Vhhhp")
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deIk"] * (*Tai)["gj"]
                        * (*Tai)["fp"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfIk"] * (*Tai)["gj"]
                        * (*Tai)["ep"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feIk"] * (*Tai)["gj"]
                        * (*Tai)["dp"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deIj"] * (*Tai)["gk"]
                        * (*Tai)["fp"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfIj"] * (*Tai)["gk"]
                        * (*Tai)["ep"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feIj"] * (*Tai)["gk"]
                        * (*Tai)["dp"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deIk"] * (*Tai)["gi"]
                        * (*Tai)["fp"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfIk"] * (*Tai)["gi"]
                        * (*Tai)["ep"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feIk"] * (*Tai)["gi"]
                        * (*Tai)["dp"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deIj"] * (*Tai)["gi"]
                        * (*Tai)["fp"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfIj"] * (*Tai)["gi"]
                        * (*Tai)["ep"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feIj"] * (*Tai)["gi"]
                        * (*Tai)["dp"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deIi"] * (*Tai)["gk"]
                        * (*Tai)["fp"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfIi"] * (*Tai)["gk"]
                        * (*Tai)["ep"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feIi"] * (*Tai)["gk"]
                        * (*Tai)["dp"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deIi"] * (*Tai)["gj"]
                        * (*Tai)["fp"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfIi"] * (*Tai)["gj"]
                        * (*Tai)["ep"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feIi"] * (*Tai)["gj"]
                        * (*Tai)["dp"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adjk"] * (*Tai)["eo"]
                        * (*Tai)["fp"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aejk"] * (*Tai)["do"]
                        * (*Tai)["fp"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afjk"] * (*Tai)["do"]
                        * (*Tai)["ep"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adki"] * (*Tai)["eo"]
                        * (*Tai)["fp"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aeki"] * (*Tai)["do"]
                        * (*Tai)["fp"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afki"] * (*Tai)["do"]
                        * (*Tai)["ep"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adij"] * (*Tai)["eo"]
                        * (*Tai)["fp"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aeij"] * (*Tai)["do"]
                        * (*Tai)["fp"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afij"] * (*Tai)["do"]
                        * (*Tai)["ep"] * (*Vijka)["opkA"];
  LDEBUG("T1 * T1 * T2 * Vhppp")
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deIk"] * (*Tai)["gi"]
                        * (*Tai)["hj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfIk"] * (*Tai)["gi"]
                        * (*Tai)["hj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feIk"] * (*Tai)["gi"]
                        * (*Tai)["hj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deIj"] * (*Tai)["gi"]
                        * (*Tai)["hk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfIj"] * (*Tai)["gi"]
                        * (*Tai)["hk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feIj"] * (*Tai)["gi"]
                        * (*Tai)["hk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deIi"] * (*Tai)["gj"]
                        * (*Tai)["hk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfIi"] * (*Tai)["gj"]
                        * (*Tai)["hk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feIi"] * (*Tai)["gj"]
                        * (*Tai)["hk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Adjk"] * (*Tai)["gi"]
                        * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Aejk"] * (*Tai)["gi"]
                        * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adjk"] * (*Tai)["gi"]
                        * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aejk"] * (*Tai)["gi"]
                        * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Afjk"] * (*Tai)["gi"]
                        * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afjk"] * (*Tai)["gi"]
                        * (*Tai)["ep"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adik"] * (*Tai)["gj"]
                        * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aeik"] * (*Tai)["gj"]
                        * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Adik"] * (*Tai)["gj"]
                        * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Aeik"] * (*Tai)["gj"]
                        * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afik"] * (*Tai)["gj"]
                        * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Afik"] * (*Tai)["gj"]
                        * (*Tai)["ep"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Adji"] * (*Tai)["gk"]
                        * (*Tai)["ep"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Aeji"] * (*Tai)["gk"]
                        * (*Tai)["dp"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Adji"] * (*Tai)["gk"]
                        * (*Tai)["fp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Aeji"] * (*Tai)["gk"]
                        * (*Tai)["fp"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Afji"] * (*Tai)["gk"]
                        * (*Tai)["dp"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Afji"] * (*Tai)["gk"]
                        * (*Tai)["ep"] * (*Viabc)["pdgA"];
  LDEBUG("T1 * T1 * T3 * V")
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["defIJk"] * (*Tai)["gi"]
                        * (*Tai)["hj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["defIJj"] * (*Tai)["gi"]
                        * (*Tai)["hk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["defIJi"] * (*Tai)["gj"]
                        * (*Tai)["hk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["AdeJjk"] * (*Tai)["gi"]
                        * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["AdfJjk"] * (*Tai)["gi"]
                        * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["AfeJjk"] * (*Tai)["gi"]
                        * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["AdeJik"] * (*Tai)["gj"]
                        * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["AdfJik"] * (*Tai)["gj"]
                        * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["AfeJik"] * (*Tai)["gj"]
                        * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["AdeJji"] * (*Tai)["gk"]
                        * (*Tai)["fp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["AdfJji"] * (*Tai)["gk"]
                        * (*Tai)["ep"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["AfeJji"] * (*Tai)["gk"]
                        * (*Tai)["dp"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["defJjk"] * (*Tai)["gi"]
                        * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["defJik"] * (*Tai)["gj"]
                        * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["defJji"] * (*Tai)["gk"]
                        * (*Tai)["hI"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["ABdijk"] * (*Tai)["eo"]
                        * (*Tai)["fp"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcijk)["ABeijk"] * (*Tai)["do"]
                        * (*Tai)["fp"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcijk)["ABfijk"] * (*Tai)["do"]
                        * (*Tai)["ep"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabcijk)["Bdeijk"] * (*Tai)["fo"]
                        * (*Tai)["hI"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["Bdfijk"] * (*Tai)["eo"]
                        * (*Tai)["hI"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabcijk)["Bfeijk"] * (*Tai)["do"]
                        * (*Tai)["hI"] * (*Vijab)["oIhB"];
  LDEBUG("T2 * T2 * T1 * V")
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["gfjk"] * (*Tabij)["depI"]
                        * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["gejk"] * (*Tabij)["dfpI"]
                        * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["gdjk"] * (*Tabij)["fepI"]
                        * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["gfki"] * (*Tabij)["depI"]
                        * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["geki"] * (*Tabij)["dfpI"]
                        * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["gdki"] * (*Tabij)["fepI"]
                        * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["gfij"] * (*Tabij)["depI"]
                        * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["geij"] * (*Tabij)["dfpI"]
                        * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["gdij"] * (*Tabij)["fepI"]
                        * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["efoj"] * (*Tabij)["hdIk"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["fdoj"] * (*Tabij)["heIk"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deok"] * (*Tabij)["hfIj"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deoj"] * (*Tabij)["hfIk"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["fdok"] * (*Tabij)["heIj"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["efok"] * (*Tabij)["hdIj"]
                        * (*Tai)["Bi"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["efoi"] * (*Tabij)["hdIk"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["fdoi"] * (*Tabij)["heIk"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deok"] * (*Tabij)["hfIi"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deoi"] * (*Tabij)["hfIk"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["fdok"] * (*Tabij)["heIi"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["efok"] * (*Tabij)["hdIi"]
                        * (*Tai)["Bj"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["efoi"] * (*Tabij)["hdIj"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["fdoi"] * (*Tabij)["heIj"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deoj"] * (*Tabij)["hfIi"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deoi"] * (*Tabij)["hfIj"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["fdoj"] * (*Tabij)["heIi"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["efoj"] * (*Tabij)["hdIi"]
                        * (*Tai)["Bk"] * (*Vijab)["oIhB"];

  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["geij"] * (*Tabij)["hdIk"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdij"] * (*Tabij)["heIk"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gfij"] * (*Tabij)["hdIk"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfij"] * (*Tabij)["heIk"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gdij"] * (*Tabij)["hfIk"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["geij"] * (*Tabij)["hfIk"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["geik"] * (*Tabij)["hdIj"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gdik"] * (*Tabij)["heIj"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfik"] * (*Tabij)["hdIj"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gfik"] * (*Tabij)["heIj"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdik"] * (*Tabij)["hfIj"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["geik"] * (*Tabij)["hfIj"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gdkj"] * (*Tabij)["heIi"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gekj"] * (*Tabij)["hdIi"]
                        * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdkj"] * (*Tabij)["hfIi"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gekj"] * (*Tabij)["hfIi"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfkj"] * (*Tabij)["hdIi"]
                        * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gfkj"] * (*Tabij)["heIi"]
                        * (*Tai)["dJ"] * (*Vijab)["IJgh"];

  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["deok"] * (*Tabij)["hAij"]
                        * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["dfok"] * (*Tabij)["hAij"]
                        * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["feok"] * (*Tabij)["hAij"]
                        * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["deoj"] * (*Tabij)["hAik"]
                        * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["dfoj"] * (*Tabij)["hAik"]
                        * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["feoj"] * (*Tabij)["hAik"]
                        * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabij)["deoi"] * (*Tabij)["hAkj"]
                        * (*Tai)["fJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["dfoi"] * (*Tabij)["hAkj"]
                        * (*Tai)["eJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabij)["feoi"] * (*Tabij)["hAkj"]
                        * (*Tai)["dJ"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gfij"] * (*Tabij)["depk"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["geij"] * (*Tabij)["dfpk"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gdij"] * (*Tabij)["fepk"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfik"] * (*Tabij)["depj"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["geik"] * (*Tabij)["dfpj"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdik"] * (*Tabij)["fepj"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gdkj"] * (*Tabij)["fepi"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["gekj"] * (*Tabij)["dfpi"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["gfkj"] * (*Tabij)["depi"]
                        * (*Tai)["AJ"] * (*Vijab)["pJgA"];

  LDEBUG("T1 * T1 * T1 * V")
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deJk"] * (*Tai)["hj"]
                        * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfJk"] * (*Tai)["hj"]
                        * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feJk"] * (*Tai)["hj"]
                        * (*Tai)["dI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["deJj"] * (*Tai)["hk"]
                        * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["dfJj"] * (*Tai)["hk"]
                        * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["feJj"] * (*Tai)["hk"]
                        * (*Tai)["dI"] * (*Tai)["gi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["deJi"] * (*Tai)["hk"]
                        * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["dfJi"] * (*Tai)["hk"]
                        * (*Tai)["eI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["feJi"] * (*Tai)["hk"]
                        * (*Tai)["dI"] * (*Tai)["gj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Bdjk"] * (*Tai)["ep"]
                        * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Bejk"] * (*Tai)["dp"]
                        * (*Tai)["fI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Bfjk"] * (*Tai)["dp"]
                        * (*Tai)["eI"] * (*Tai)["gi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Bdik"] * (*Tai)["ep"]
                        * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Beik"] * (*Tai)["dp"]
                        * (*Tai)["fI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Bfik"] * (*Tai)["dp"]
                        * (*Tai)["eI"] * (*Tai)["gj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Bdji"] * (*Tai)["ep"]
                        * (*Tai)["fI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tabij)["Beji"] * (*Tai)["dp"]
                        * (*Tai)["fI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tabij)["Bfji"] * (*Tai)["dp"]
                        * (*Tai)["eI"] * (*Tai)["gk"] * (*Vijab)["pIgB"];
#ifdef DEBUG
  LOG(1, getAbbreviation()) << "Triples done" << std::endl;
#endif

  return residuum;
}
