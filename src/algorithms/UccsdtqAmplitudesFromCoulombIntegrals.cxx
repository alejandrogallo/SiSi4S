#include <unistd.h>

#include <Sisi4s.hpp>

#include <algorithms/UccsdtqAmplitudesFromCoulombIntegrals.hpp>
#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>

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
    UccsdtqAmplitudesFromCoulombIntegrals,
    SPEC_IN(UCCSD_SPEC_IN),
    SPEC_OUT(UCCSD_SPEC_OUT));


IMPLEMENT_ALGORITHM(UccsdtqAmplitudesFromCoulombIntegrals) {
  in.set<bool>("unrestricted", true);
  in.set<bool>("antisymmetrize", true);
  ClusterSinglesDoublesTriplesQuadruplesAlgorithm::run();
}

PTR(FockVector<sisi4s::complex>) UccsdtqAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<sisi4s::complex>) &amplitudes) {
  if (iterationStep == 0) {
    LOG(1, getAbbreviation())
        << "WARNING: Using complex version of Uccsdtq" << std::endl;
    LOG(1, getAbbreviation())
        << "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UccsdtqAmplitudesFromCoulombIntegrals::getResiduum(
    const int iterationStep,
    const PTR(const FockVector<double>) &amplitudes) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
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
    Fia = NULL;
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
  auto Tabcdijkl(amplitudes->get(3));
  Tabcdijkl->set_name("Tabcdijkl");

  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");
  auto Rabcijk(residuum->get(2));
  Rabcijk->set_name("Rabcijk");
  auto Rabcdijkl(residuum->get(2));
  Rabcdijkl->set_name("Rabcdijkl");

  if (Fia) {
    (*Rai)["bi"] += (-1.0) * (*Tai)["ci"] * (*Tai)["bl"] * (*Fia)["lc"];
    (*Rai)["bi"] += (+1.0) * (*Fia)["kd"] * (*Tabij)["dbki"];
    (*Rai)["bi"] += (+1.0) * (*Fia)["ib"];
  }
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

  if (Fia) {
    (*Rabij)["cdij"] += (-1.0) * (*Fia)["mf"] * (*Tabij)["cdmj"] * (*Tai)["fi"];
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabij)["cdmi"] * (*Tai)["fj"];
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabij)["fcij"] * (*Tai)["dm"];
    (*Rabij)["cdij"] += (-1.0) * (*Fia)["mf"] * (*Tabij)["fdij"] * (*Tai)["cm"];
    (*Rabij)["cdij"] += (+1.0) * (*Fia)["mf"] * (*Tabcijk)["fcdmij"];
  }
  (*Rabij)["cdij"] += (+1.0) * (*Vabij)["cdij"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["cm"] * (*Viajk)["mdij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["dm"] * (*Viajk)["mcij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ej"] * (*Vabic)["cdie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ei"] * (*Vabic)["cdje"];
  //(*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
  //(*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
  //(*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
  //(*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += (-1.0) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += (+1.0) * (*Tabij)["edni"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabcijk)["ecdnoj"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] += (-0.5) * (*Tabcijk)["ecdnoi"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] += (-0.5) * (*Tabcijk)["efcoij"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] += (+0.5) * (*Tabcijk)["efdoij"] * (*Viabc)["ocef"];
  (*Rabij)["cdij"] += (+0.25) * (*Tabcdijkl)["efcdopij"] * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
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
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];
  (*Rabij)["cdij"] +=
      (-1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
  (*Rabij)["cdij"] +=
      (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];
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
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["fo"] * (*Tabij)["cdpj"]
                    * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["ej"] * (*Tai)["fo"] * (*Tabij)["cdpi"]
                    * (*Vijab)["opef"];
  (*Rabij)["cdij"] += (+0.5) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"]
                    * (*Vijab)["mngh"];
  (*Rabij)["cdij"] += (-1.0) * (*Tai)["dm"] * (*Tai)["fo"] * (*Tabij)["hcij"]
                    * (*Vijab)["mofh"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["cm"] * (*Tai)["fo"] * (*Tabij)["hdij"]
                    * (*Vijab)["mofh"];
  (*Rabij)["cdij"] += (+1.0) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"]
                    * (*Tai)["dp"] * (*Vijab)["opef"];

  if (Fia) {
    (*Rabcijk)["defijk"] += (+1.0) * (*Fia)["oh"] * (*Tabcdijkl)["hdefoijk"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Fia)["oh"] * (*Tabcijk)["defojk"] * (*Tai)["hi"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Fia)["oh"] * (*Tabcijk)["defoki"] * (*Tai)["hj"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Fia)["oh"] * (*Tabcijk)["defoij"] * (*Tai)["hk"];
    (*Rabcijk)["defijk"] +=
        (-1.0) * (*Fia)["oh"] * (*Tabcijk)["hdeijk"] * (*Tai)["fo"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Fia)["oh"] * (*Tabcijk)["hdfijk"] * (*Tai)["eo"];
    (*Rabcijk)["defijk"] +=
        (+1.0) * (*Fia)["oh"] * (*Tabcijk)["hfeijk"] * (*Tai)["do"];
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
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fij)["oi"] * (*Tabcijk)["defojk"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fij)["oj"] * (*Tabcijk)["defoik"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fij)["ok"] * (*Tabcijk)["defoji"];
  //(*Rabcijk)["defijk"] += ( + 1.0  ) * (*Fab)["fg"] * (*Tabcijk)["gdeijk"];
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fab)["eg"] * (*Tabcijk)["gdfijk"];
  //(*Rabcijk)["defijk"] += ( - 1.0  ) * (*Fab)["dg"] * (*Tabcijk)["gfeijk"];
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
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcdijkl)["gdefpIjk"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcdijkl)["gdefpIki"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcdijkl)["gdefpIij"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tabcdijkl)["ghdeIijk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcdijkl)["ghdfIijk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tabcdijkl)["ghfeIijk"] * (*Viabc)["Idgh"];
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
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcdijkl)["gdefpIjk"] * (*Tai)["Bi"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcdijkl)["gdefpIki"] * (*Tai)["Bj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcdijkl)["gdefpIij"] * (*Tai)["Bk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcdijkl)["ghdeIijk"] * (*Tai)["fJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcdijkl)["ghdfIijk"] * (*Tai)["eJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcdijkl)["ghfeIijk"] * (*Tai)["dJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+1.0) * (*Tabcdijkl)["gdefpijk"] * (*Tai)["AJ"] * (*Vijab)["pJgA"];
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
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gdepIk"] * (*Tabij)["Bfij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gdfpIk"] * (*Tabij)["Beij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gfepIk"] * (*Tabij)["Bdij"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gdepIj"] * (*Tabij)["Bfik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gdfpIj"] * (*Tabij)["Beik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gfepIj"] * (*Tabij)["Bdik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gdepIi"] * (*Tabij)["Bfkj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gdfpIi"] * (*Tabij)["Bekj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gfepIi"] * (*Tabij)["Bdkj"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] +=
      (+0.25) * (*Tabcijk)["defopk"] * (*Tabij)["ABij"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["defopj"] * (*Tabij)["ABik"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["defopi"] * (*Tabij)["ABkj"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["ghdIjk"] * (*Tabij)["efJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gheIjk"] * (*Tabij)["dfJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["ghfIjk"] * (*Tabij)["edJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["ghdIki"] * (*Tabij)["efJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gheIki"] * (*Tabij)["dfJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["ghfIki"] * (*Tabij)["edJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["ghdIij"] * (*Tabij)["efJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gheIij"] * (*Tabij)["dfJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["ghfIij"] * (*Tabij)["edJk"] * (*Vijab)["IJgh"];
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
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defojk"] * (*Tabij)["hAJi"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defoki"] * (*Tabij)["hAJj"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["defoij"] * (*Tabij)["hAJk"] * (*Vijab)["oJhA"];
  (*Rabcijk)["defijk"] +=
      (+0.25) * (*Tabcijk)["ghdijk"] * (*Tabij)["efIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["gheijk"] * (*Tabij)["dfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.25) * (*Tabcijk)["ghfijk"] * (*Tabij)["edIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (+0.5) * (*Tabcijk)["gdeijk"] * (*Tabij)["hfIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gdfijk"] * (*Tabij)["heIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] +=
      (-0.5) * (*Tabcijk)["gfeijk"] * (*Tabij)["hdIJ"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["fp"]
                        * (*Tabij)["deIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["ep"]
                        * (*Tabij)["dfIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["dp"]
                        * (*Tabij)["feIk"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["fp"]
                        * (*Tabij)["deIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["ep"]
                        * (*Tabij)["dfIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["dp"]
                        * (*Tabij)["feIj"] * (*Vijka)["pIig"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["fp"]
                        * (*Tabij)["deIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["ep"]
                        * (*Tabij)["dfIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["dp"]
                        * (*Tabij)["feIk"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["fp"]
                        * (*Tabij)["deIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["ep"]
                        * (*Tabij)["dfIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["dp"]
                        * (*Tabij)["feIj"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["fp"]
                        * (*Tabij)["deIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["ep"]
                        * (*Tabij)["dfIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["dp"]
                        * (*Tabij)["feIi"] * (*Vijka)["pIjg"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["fp"]
                        * (*Tabij)["deIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["ep"]
                        * (*Tabij)["dfIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["dp"]
                        * (*Tabij)["feIi"] * (*Vijka)["pIkg"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["eo"] * (*Tai)["fp"]
                        * (*Tabij)["Adjk"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["do"] * (*Tai)["fp"]
                        * (*Tabij)["Aejk"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["do"] * (*Tai)["ep"]
                        * (*Tabij)["Afjk"] * (*Vijka)["opiA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["eo"] * (*Tai)["fp"]
                        * (*Tabij)["Adki"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["do"] * (*Tai)["fp"]
                        * (*Tabij)["Aeki"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["do"] * (*Tai)["ep"]
                        * (*Tabij)["Afki"] * (*Vijka)["opjA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["eo"] * (*Tai)["fp"]
                        * (*Tabij)["Adij"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["do"] * (*Tai)["fp"]
                        * (*Tabij)["Aeij"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["do"] * (*Tai)["ep"]
                        * (*Tabij)["Afij"] * (*Vijka)["opkA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hj"]
                        * (*Tabij)["deIk"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hj"]
                        * (*Tabij)["dfIk"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hj"]
                        * (*Tabij)["feIk"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hk"]
                        * (*Tabij)["deIj"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hk"]
                        * (*Tabij)["dfIj"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hk"]
                        * (*Tabij)["feIj"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["hk"]
                        * (*Tabij)["deIi"] * (*Viabc)["Ifgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["hk"]
                        * (*Tabij)["dfIi"] * (*Viabc)["Iegh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["hk"]
                        * (*Tabij)["feIi"] * (*Viabc)["Idgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["ep"]
                        * (*Tabij)["Adjk"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["dp"]
                        * (*Tabij)["Aejk"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["fp"]
                        * (*Tabij)["Adjk"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["fp"]
                        * (*Tabij)["Aejk"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["dp"]
                        * (*Tabij)["Afjk"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["ep"]
                        * (*Tabij)["Afjk"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["ep"]
                        * (*Tabij)["Adik"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["dp"]
                        * (*Tabij)["Aeik"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["fp"]
                        * (*Tabij)["Adik"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["fp"]
                        * (*Tabij)["Aeik"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["dp"]
                        * (*Tabij)["Afik"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["ep"]
                        * (*Tabij)["Afik"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["ep"]
                        * (*Tabij)["Adji"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["dp"]
                        * (*Tabij)["Aeji"] * (*Viabc)["pfgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["fp"]
                        * (*Tabij)["Adji"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["fp"]
                        * (*Tabij)["Aeji"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["dp"]
                        * (*Tabij)["Afji"] * (*Viabc)["pegA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["ep"]
                        * (*Tabij)["Afji"] * (*Viabc)["pdgA"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tai)["gi"] * (*Tai)["hj"]
                        * (*Tabcijk)["defIJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tai)["gi"] * (*Tai)["hk"]
                        * (*Tabcijk)["defIJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tai)["gj"] * (*Tai)["hk"]
                        * (*Tabcijk)["defIJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["fp"]
                        * (*Tabcijk)["AdeJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["ep"]
                        * (*Tabcijk)["AdfJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["dp"]
                        * (*Tabcijk)["AfeJjk"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["fp"]
                        * (*Tabcijk)["AdeJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["ep"]
                        * (*Tabcijk)["AdfJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["dp"]
                        * (*Tabcijk)["AfeJik"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["fp"]
                        * (*Tabcijk)["AdeJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["ep"]
                        * (*Tabcijk)["AdfJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["dp"]
                        * (*Tabcijk)["AfeJji"] * (*Vijab)["pJgA"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hI"]
                        * (*Tabcijk)["defJjk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["hI"]
                        * (*Tabcijk)["defJik"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["hI"]
                        * (*Tabcijk)["defJji"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tai)["eo"] * (*Tai)["fp"]
                        * (*Tabcijk)["ABdijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (-0.5) * (*Tai)["do"] * (*Tai)["fp"]
                        * (*Tabcijk)["ABeijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (+0.5) * (*Tai)["do"] * (*Tai)["ep"]
                        * (*Tabcijk)["ABfijk"] * (*Vijab)["opAB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["fo"] * (*Tai)["hI"]
                        * (*Tabcijk)["Bdeijk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["eo"] * (*Tai)["hI"]
                        * (*Tabcijk)["Bdfijk"] * (*Vijab)["oIhB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["do"] * (*Tai)["hI"]
                        * (*Tabcijk)["Bfeijk"] * (*Vijab)["oIhB"];
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
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["fI"]
                        * (*Tabij)["deJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["eI"]
                        * (*Tabij)["dfJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hj"] * (*Tai)["dI"]
                        * (*Tabij)["feJk"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["fI"]
                        * (*Tabij)["deJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["eI"]
                        * (*Tabij)["dfJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["hk"] * (*Tai)["dI"]
                        * (*Tabij)["feJj"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["fI"]
                        * (*Tabij)["deJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["eI"]
                        * (*Tabij)["dfJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["hk"] * (*Tai)["dI"]
                        * (*Tabij)["feJi"] * (*Vijab)["IJgh"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["ep"] * (*Tai)["fI"]
                        * (*Tabij)["Bdjk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tai)["fI"]
                        * (*Tabij)["Bejk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gi"] * (*Tai)["dp"] * (*Tai)["eI"]
                        * (*Tabij)["Bfjk"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["ep"] * (*Tai)["fI"]
                        * (*Tabij)["Bdik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tai)["fI"]
                        * (*Tabij)["Beik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gj"] * (*Tai)["dp"] * (*Tai)["eI"]
                        * (*Tabij)["Bfik"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["ep"] * (*Tai)["fI"]
                        * (*Tabij)["Bdji"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (+1.0) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tai)["fI"]
                        * (*Tabij)["Beji"] * (*Vijab)["pIgB"];
  (*Rabcijk)["defijk"] += (-1.0) * (*Tai)["gk"] * (*Tai)["dp"] * (*Tai)["eI"]
                        * (*Tabij)["Bfji"] * (*Vijab)["pIgB"];

  if (Fia) {
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIjkl"] * (*Tai)["Bi"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIlik"] * (*Tai)["Bj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIilj"] * (*Tai)["Bk"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIijk"] * (*Tai)["Bl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Befgijkl"] * (*Tai)["hI"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Befhijkl"] * (*Tai)["gI"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Behgijkl"] * (*Tai)["fI"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Bhfgijkl"] * (*Tai)["eI"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIkl"] * (*Tabij)["Bhij"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIkl"] * (*Tabij)["Bgij"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIkl"] * (*Tabij)["Bfij"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIkl"] * (*Tabij)["Beij"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIlj"] * (*Tabij)["Bhik"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIlj"] * (*Tabij)["Bgik"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIlj"] * (*Tabij)["Bfik"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIlj"] * (*Tabij)["Beik"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIjk"] * (*Tabij)["Bhil"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIjk"] * (*Tabij)["Bgil"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIjk"] * (*Tabij)["Bfil"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIjk"] * (*Tabij)["Beil"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIli"] * (*Tabij)["Bhkj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIli"] * (*Tabij)["Bgkj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIli"] * (*Tabij)["Bfkj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIli"] * (*Tabij)["Bekj"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIik"] * (*Tabij)["Bhlj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIik"] * (*Tabij)["Bglj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIik"] * (*Tabij)["Bflj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIik"] * (*Tabij)["Belj"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIij"] * (*Tabij)["Bhkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIij"] * (*Tabij)["Bgkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIij"] * (*Tabij)["Bfkl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIij"] * (*Tabij)["Bekl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Befjkl"] * (*Tabij)["ghIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Begjkl"] * (*Tabij)["fhIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfjkl"] * (*Tabij)["ehIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Behjkl"] * (*Tabij)["gfIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfjkl"] * (*Tabij)["geIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghjkl"] * (*Tabij)["efIi"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Beflik"] * (*Tabij)["ghIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Beglik"] * (*Tabij)["fhIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgflik"] * (*Tabij)["ehIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behlik"] * (*Tabij)["gfIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhflik"] * (*Tabij)["geIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghlik"] * (*Tabij)["efIj"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Befilj"] * (*Tabij)["ghIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Begilj"] * (*Tabij)["fhIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfilj"] * (*Tabij)["ehIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behilj"] * (*Tabij)["gfIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfilj"] * (*Tabij)["geIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghilj"] * (*Tabij)["efIk"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Befijk"] * (*Tabij)["ghIl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Begijk"] * (*Tabij)["fhIl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfijk"] * (*Tabij)["ehIl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behijk"] * (*Tabij)["gfIl"];
    (*Rabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfijk"] * (*Tabij)["geIl"];
    (*Rabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghijk"] * (*Tabij)["efIl"];
  }

  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIkl"] * (*Viajk)["Ihij"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIkl"] * (*Viajk)["Igij"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIkl"] * (*Viajk)["Ifij"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIkl"] * (*Viajk)["Ieij"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIlj"] * (*Viajk)["Ihik"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIlj"] * (*Viajk)["Igik"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIlj"] * (*Viajk)["Ifik"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIlj"] * (*Viajk)["Ieik"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIjk"] * (*Viajk)["Ihil"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIjk"] * (*Viajk)["Igil"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIjk"] * (*Viajk)["Ifil"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIjk"] * (*Viajk)["Ieil"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIli"] * (*Viajk)["Ihkj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIli"] * (*Viajk)["Igkj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIli"] * (*Viajk)["Ifkj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIli"] * (*Viajk)["Iekj"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIik"] * (*Viajk)["Ihlj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIik"] * (*Viajk)["Iglj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIik"] * (*Viajk)["Iflj"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIik"] * (*Viajk)["Ielj"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIij"] * (*Viajk)["Ihkl"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIij"] * (*Viajk)["Igkl"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIij"] * (*Viajk)["Ifkl"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIij"] * (*Viajk)["Iekl"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefjkl"] * (*Vabic)["ghiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegjkl"] * (*Vabic)["fhiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfjkl"] * (*Vabic)["ehiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehjkl"] * (*Vabic)["gfiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Vabic)["geiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aghjkl"] * (*Vabic)["efiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeflik"] * (*Vabic)["ghjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeglik"] * (*Vabic)["fhjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agflik"] * (*Vabic)["ehjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehlik"] * (*Vabic)["gfjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahflik"] * (*Vabic)["gejA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghlik"] * (*Vabic)["efjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefilj"] * (*Vabic)["ghkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegilj"] * (*Vabic)["fhkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfilj"] * (*Vabic)["ehkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehilj"] * (*Vabic)["gfkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfilj"] * (*Vabic)["gekA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghilj"] * (*Vabic)["efkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefijk"] * (*Vabic)["ghlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegijk"] * (*Vabic)["fhlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfijk"] * (*Vabic)["ehlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehijk"] * (*Vabic)["gflA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfijk"] * (*Vabic)["gelA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghijk"] * (*Vabic)["eflA"];
  //(*Rabcdijkl)["efghijkl"] += ( - 1.0  ) * (*Fij)["Ii"] *
  //(*Tabcdijkl)["efghIjkl"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fij)["Ij"] *
  //(*Tabcdijkl)["efghIikl"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fij)["Ik"] *
  //(*Tabcdijkl)["efghIjil"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fij)["Il"] *
  //(*Tabcdijkl)["efghIjki"];
  //(*Rabcdijkl)["efghijkl"] += ( - 1.0  ) * (*Fab)["hA"] *
  //(*Tabcdijkl)["Aefgijkl"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fab)["gA"] *
  //(*Tabcdijkl)["Aefhijkl"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fab)["fA"] *
  //(*Tabcdijkl)["Aehgijkl"];
  //(*Rabcdijkl)["efghijkl"] += ( + 1.0  ) * (*Fab)["eA"] *
  //(*Tabcdijkl)["Ahfgijkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJkl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJlj"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJjk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJli"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJik"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJij"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJlik"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJlik"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJlik"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJilj"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJilj"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJilj"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJijk"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJijk"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJijk"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABefijkl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABegijkl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABehijkl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABghijkl"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["hJ"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["gJ"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["fJ"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["eJ"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["hJ"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["gJ"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["fJ"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["eJ"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["hJ"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["gJ"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["fJ"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["eJ"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIli"] * (*Tai)["hJ"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIli"] * (*Tai)["gJ"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["fJ"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["eJ"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIik"] * (*Tai)["hJ"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIik"] * (*Tai)["gJ"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["fJ"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["eJ"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIij"] * (*Tai)["hJ"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIij"] * (*Tai)["gJ"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIij"] * (*Tai)["fJ"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIij"] * (*Tai)["eJ"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["Bj"] * (*Viajb)["IhiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["Bj"] * (*Viajb)["IgiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["Bj"] * (*Viajb)["IfiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["Bj"] * (*Viajb)["IeiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["Bk"] * (*Viajb)["IhiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["Bk"] * (*Viajb)["IgiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["Bk"] * (*Viajb)["IfiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["Bk"] * (*Viajb)["IeiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["Bl"] * (*Viajb)["IhiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["Bl"] * (*Viajb)["IgiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["Bl"] * (*Viajb)["IfiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["Bl"] * (*Viajb)["IeiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["Bi"] * (*Viajb)["IhjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["Bi"] * (*Viajb)["IgjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["Bi"] * (*Viajb)["IfjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["Bi"] * (*Viajb)["IejB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["Bi"] * (*Viajb)["IhkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["Bi"] * (*Viajb)["IgkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["Bi"] * (*Viajb)["IfkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["Bi"] * (*Viajb)["IekB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["Bi"] * (*Viajb)["IhlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["Bi"] * (*Viajb)["IglB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["Bi"] * (*Viajb)["IflB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["Bi"] * (*Viajb)["IelB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIli"] * (*Tai)["Bk"] * (*Viajb)["IhjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIli"] * (*Tai)["Bk"] * (*Viajb)["IgjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["Bk"] * (*Viajb)["IfjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["Bk"] * (*Viajb)["IejB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIik"] * (*Tai)["Bl"] * (*Viajb)["IhjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIik"] * (*Tai)["Bl"] * (*Viajb)["IgjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["Bl"] * (*Viajb)["IfjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["Bl"] * (*Viajb)["IejB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIli"] * (*Tai)["Bj"] * (*Viajb)["IhkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIli"] * (*Tai)["Bj"] * (*Viajb)["IgkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["Bj"] * (*Viajb)["IfkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["Bj"] * (*Viajb)["IekB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIik"] * (*Tai)["Bj"] * (*Viajb)["IhlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIik"] * (*Tai)["Bj"] * (*Viajb)["IglB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["Bj"] * (*Viajb)["IflB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["Bj"] * (*Viajb)["IelB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIji"] * (*Tai)["Bl"] * (*Viajb)["IhkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIji"] * (*Tai)["Bl"] * (*Viajb)["IgkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIji"] * (*Tai)["Bl"] * (*Viajb)["IfkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIji"] * (*Tai)["Bl"] * (*Viajb)["IekB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIij"] * (*Tai)["Bk"] * (*Viajb)["IhlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIij"] * (*Tai)["Bk"] * (*Viajb)["IglB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIij"] * (*Tai)["Bk"] * (*Viajb)["IflB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIij"] * (*Tai)["Bk"] * (*Viajb)["IelB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["gJ"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["fJ"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["eJ"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["hJ"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["hJ"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["hJ"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["fJ"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["eJ"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["gJ"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["gJ"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tai)["eJ"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahgjkl"] * (*Tai)["fJ"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["gJ"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["fJ"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agflik"] * (*Tai)["eJ"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["hJ"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["hJ"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agflik"] * (*Tai)["hJ"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["fJ"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["eJ"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["gJ"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["gJ"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghlik"] * (*Tai)["eJ"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahglik"] * (*Tai)["fJ"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["gJ"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["fJ"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["eJ"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["hJ"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["hJ"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["hJ"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["fJ"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["eJ"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["gJ"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["gJ"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghilj"] * (*Tai)["eJ"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahgilj"] * (*Tai)["fJ"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["gJ"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["fJ"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["eJ"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["hJ"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["hJ"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["hJ"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["fJ"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["eJ"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["gJ"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["gJ"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghijk"] * (*Tai)["eJ"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahgijk"] * (*Tai)["fJ"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["Bi"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["Bi"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["Bi"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["Bi"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["Bi"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tai)["Bi"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["Bj"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["Bj"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agflik"] * (*Tai)["Bj"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["Bj"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["Bj"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghlik"] * (*Tai)["Bj"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["Bk"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["Bk"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["Bk"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["Bk"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["Bk"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghilj"] * (*Tai)["Bk"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["Bl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["Bl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["Bl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["Bl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["Bl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghijk"] * (*Tai)["Bl"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJkl"] * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJlj"] * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJjk"] * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJkl"] * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJlj"] * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJjk"] * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJli"] * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJik"] * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJli"] * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIJik"] * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJji"] * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIJij"] * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJlik"] * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJlik"] * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJlik"] * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJilj"] * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJilj"] * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJilj"] * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJijk"] * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJijk"] * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJijk"] * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["efghIjkl"] * (*Tai)["BK"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["efghIlik"] * (*Tai)["BK"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["efghIilj"] * (*Tai)["BK"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["efghIijk"] * (*Tai)["BK"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJlik"] * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJlik"] * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJlik"] * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJilj"] * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJilj"] * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJilj"] * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJijk"] * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJijk"] * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJijk"] * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABefijkl"] * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABegijkl"] * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABefijkl"] * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABegijkl"] * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABehijkl"] * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABehijkl"] * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABghijkl"] * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABhgijkl"] * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["Aefgijkl"] * (*Tai)["BK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["Aefhijkl"] * (*Tai)["BK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["Aehgijkl"] * (*Tai)["BK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["Ahfgijkl"] * (*Tai)["BK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["ghIk"] * (*Tabij)["efJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["hfIk"] * (*Tabij)["egJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["heIk"] * (*Tabij)["gfJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["egIk"] * (*Tabij)["hfJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["efIk"] * (*Tabij)["ghJl"] * (*Vijkl)["IJij"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJl"] * (*Vijkl)["IJik"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJk"] * (*Vijkl)["IJil"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJl"] * (*Vijkl)["IJkj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJk"] * (*Vijkl)["IJlj"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJj"] * (*Vijkl)["IJkl"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJl"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJl"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJl"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJl"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJl"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJl"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJl"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJl"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJl"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJl"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJl"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJl"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJk"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJk"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJk"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJk"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJk"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJk"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJk"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJk"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJk"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJk"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJk"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJk"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJj"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJj"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"] * (*Viajb)["JhiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJj"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJj"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJj"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJj"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"] * (*Viajb)["JgiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJj"] * (*Viajb)["JfiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJj"] * (*Viajb)["JeiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agki"] * (*Tabij)["efJl"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afki"] * (*Tabij)["egJl"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJl"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJl"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJl"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJl"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJl"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJl"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJl"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJl"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJl"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJl"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agli"] * (*Tabij)["efJk"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afli"] * (*Tabij)["egJk"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["gfJk"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["efJk"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["egJk"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["gfJk"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afli"] * (*Tabij)["ehJk"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["hfJk"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agli"] * (*Tabij)["ehJk"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agli"] * (*Tabij)["hfJk"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["ghJk"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afli"] * (*Tabij)["hgJk"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJi"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJi"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJi"] * (*Viajb)["JhjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJi"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJi"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJi"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJi"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJi"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJi"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJi"] * (*Viajb)["JgjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJi"] * (*Viajb)["JfjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJi"] * (*Viajb)["JejA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afij"] * (*Tabij)["egJl"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJl"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJl"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJl"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJl"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJl"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJl"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJl"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afij"] * (*Tabij)["egJk"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJk"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJk"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJk"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJk"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJk"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJk"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJk"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afil"] * (*Tabij)["egJj"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["gfJj"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["egJj"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["gfJj"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["hfJj"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agil"] * (*Tabij)["ehJj"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agil"] * (*Tabij)["hfJj"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["ghJj"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJi"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJi"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJi"] * (*Viajb)["JhkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJi"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJi"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJi"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJi"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJi"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJi"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJi"] * (*Viajb)["JgkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJi"] * (*Viajb)["JfkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJi"] * (*Viajb)["JekA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agki"] * (*Tabij)["efJj"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afki"] * (*Tabij)["egJj"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJj"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJj"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJj"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJj"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJj"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJj"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJj"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJj"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJj"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJj"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJi"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJi"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJi"] * (*Viajb)["JhlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJi"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJi"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJi"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJi"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJi"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJi"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJi"] * (*Viajb)["JglA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJi"] * (*Viajb)["JflA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJi"] * (*Viajb)["JelA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bekl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bekl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"] * (*Vabcd)["feAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bejl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bejl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"] * (*Vabcd)["feAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bekj"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"] * (*Vabcd)["ghAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bekj"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"] * (*Vabcd)["fhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"] * (*Vabcd)["ehAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"] * (*Vabcd)["efAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"] * (*Vabcd)["gfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"] * (*Vabcd)["geAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"] * (*Vabcd)["feAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chjk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgjk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfjk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Cejk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chjl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgjl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfjl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Cejl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chlk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cglk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cflk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Celk"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chki"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgki"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfki"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceki"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chli"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgli"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfli"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Celi"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chij"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgij"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfij"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceij"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chij"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgij"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfij"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Ceij"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chil"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgil"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfil"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceil"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chki"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgki"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfki"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceki"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chlk"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cglk"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cflk"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Celk"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKk"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKi"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKi"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKi"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKk"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKj"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJji"] * (*Tabij)["ghKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJji"] * (*Tabij)["fhKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJji"] * (*Tabij)["ehKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJji"] * (*Tabij)["gfKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJji"] * (*Tabij)["geKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJji"] * (*Tabij)["efKl"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["ghKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["fhKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["ehKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["gfKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["geKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["efKk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKj"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKj"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKj"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKj"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKk"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKk"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKk"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKk"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKl"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKl"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKl"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKl"] * (*Vijka)["IKiB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKi"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKi"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKi"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKi"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKi"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKi"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKi"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKi"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKi"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKi"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKi"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKi"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKk"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKk"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKk"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKk"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKl"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKl"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKl"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKl"] * (*Vijka)["IKjB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKj"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKj"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKj"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKj"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKj"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKj"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKj"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKj"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIji"] * (*Tabij)["BhKl"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIji"] * (*Tabij)["BgKl"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIji"] * (*Tabij)["BfKl"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIji"] * (*Tabij)["BeKl"] * (*Vijka)["IKkB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["efgIij"] * (*Tabij)["BhKk"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["efhIij"] * (*Tabij)["BgKk"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["ehgIij"] * (*Tabij)["BfKk"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["hfgIij"] * (*Tabij)["BeKk"] * (*Vijka)["IKlB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aefjkl"] * (*Tabij)["ghJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aegjkl"] * (*Tabij)["fhJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Agfjkl"] * (*Tabij)["ehJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aehjkl"] * (*Tabij)["gfJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["geJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aghjkl"] * (*Tabij)["efJK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aeflik"] * (*Tabij)["ghJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aeglik"] * (*Tabij)["fhJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Agflik"] * (*Tabij)["ehJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aehlik"] * (*Tabij)["gfJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Ahflik"] * (*Tabij)["geJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghlik"] * (*Tabij)["efJK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aefilj"] * (*Tabij)["ghJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aegilj"] * (*Tabij)["fhJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Agfilj"] * (*Tabij)["ehJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aehilj"] * (*Tabij)["gfJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Ahfilj"] * (*Tabij)["geJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghilj"] * (*Tabij)["efJK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aefijk"] * (*Tabij)["ghJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aegijk"] * (*Tabij)["fhJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Agfijk"] * (*Tabij)["ehJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aehijk"] * (*Tabij)["gfJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Ahfijk"] * (*Tabij)["geJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghijk"] * (*Tabij)["efJK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Chij"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Chij"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Chij"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Chik"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Chik"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Chik"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Chil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Chil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Chil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Cekj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Chkj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Chkj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Chkj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cekj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["Cekj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Cglj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Cflj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Celj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Chlj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Chlj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Chlj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cflj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Celj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cglj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Cglj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["Celj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJik"] * (*Tabij)["Cflj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Cekl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Chkl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Chkl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Chkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cekl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["Cekl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhgJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIkl"] * (*Tabij)["BCij"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIlj"] * (*Tabij)["BCik"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIjk"] * (*Tabij)["BCil"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIli"] * (*Tabij)["BCkj"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIik"] * (*Tabij)["BClj"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIik"] * (*Tabij)["BClj"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIik"] * (*Tabij)["BClj"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIik"] * (*Tabij)["BClj"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IhBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efhIij"] * (*Tabij)["BCkl"] * (*Viabc)["IgBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IfBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hfgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IeBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fgKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["egKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["feKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fhKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["ehKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["hgKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["hgKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["heKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["fhKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["efKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["geKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["fgKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fgKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["egKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["feKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fhKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["ehKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["hgKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["hgKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["heKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["fhKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["efKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["geKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["fgKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fgKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["egKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["feKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fhKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["ehKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["hgKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["hgKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["heKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["fhKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["efKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["geKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["fgKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fgKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["egKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["feKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fhKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["ehKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["hgKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["hgKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["heKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["fhKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["efKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["geKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["fgKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahgjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BgKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BfKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BeKj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BhKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BhKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BhKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BfKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BeKj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BgKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BgKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghlik"] * (*Tabij)["BeKj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahglik"] * (*Tabij)["BfKj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BgKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BfKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BeKk"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BhKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BhKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BhKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BfKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BeKk"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BgKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BgKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghilj"] * (*Tabij)["BeKk"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahgilj"] * (*Tabij)["BfKk"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BgKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BfKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BeKl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BhKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BhKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BhKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BfKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BeKl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BgKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BgKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Aghijk"] * (*Tabij)["BeKl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["Ahgijk"] * (*Tabij)["BfKl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKkl"] * (*Tabij)["Dhij"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKkl"] * (*Tabij)["Dgij"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKkl"] * (*Tabij)["Dfij"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKkl"] * (*Tabij)["Deij"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKlj"] * (*Tabij)["Dhik"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKlj"] * (*Tabij)["Dgik"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKlj"] * (*Tabij)["Dfik"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKlj"] * (*Tabij)["Deik"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKjk"] * (*Tabij)["Dhil"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKjk"] * (*Tabij)["Dgil"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKjk"] * (*Tabij)["Dfil"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKjk"] * (*Tabij)["Deil"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKli"] * (*Tabij)["Dhkj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKli"] * (*Tabij)["Dgkj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKli"] * (*Tabij)["Dfkj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKli"] * (*Tabij)["Dekj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKik"] * (*Tabij)["Dhlj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKik"] * (*Tabij)["Dglj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKik"] * (*Tabij)["Dflj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKik"] * (*Tabij)["Delj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["AefgJKij"] * (*Tabij)["Dhkl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AefhJKij"] * (*Tabij)["Dgkl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AehgJKij"] * (*Tabij)["Dfkl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["AhfgJKij"] * (*Tabij)["Dekl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJkl"] * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJlj"] * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJjk"] * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJli"] * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJik"] * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["efghIJij"] * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABefKjkl"] * (*Tabij)["ghLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABegKjkl"] * (*Tabij)["fhLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABgfKjkl"] * (*Tabij)["ehLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABehKjkl"] * (*Tabij)["gfLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABhfKjkl"] * (*Tabij)["geLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABghKjkl"] * (*Tabij)["efLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABefKlik"] * (*Tabij)["ghLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABegKlik"] * (*Tabij)["fhLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABgfKlik"] * (*Tabij)["ehLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABehKlik"] * (*Tabij)["gfLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABhfKlik"] * (*Tabij)["geLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABghKlik"] * (*Tabij)["efLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABefKilj"] * (*Tabij)["ghLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABegKilj"] * (*Tabij)["fhLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABgfKilj"] * (*Tabij)["ehLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABehKilj"] * (*Tabij)["gfLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABhfKilj"] * (*Tabij)["geLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABghKilj"] * (*Tabij)["efLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABefKijk"] * (*Tabij)["ghLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABegKijk"] * (*Tabij)["fhLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABgfKijk"] * (*Tabij)["ehLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABehKijk"] * (*Tabij)["gfLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["ABhfKijk"] * (*Tabij)["geLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["ABghKijk"] * (*Tabij)["efLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJlik"] * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJlik"] * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJlik"] * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJilj"] * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJilj"] * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJilj"] * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcdijkl)["AefgJijk"] * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AefhJijk"] * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AehgJijk"] * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["efghIjkl"] * (*Tabij)["BCLi"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIlik"] * (*Tabij)["BCLj"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIilj"] * (*Tabij)["BCLk"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["efghIijk"] * (*Tabij)["BCLl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["ABefijkl"] * (*Tabij)["ghKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcdijkl)["ABegijkl"] * (*Tabij)["fhKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcdijkl)["ABgfijkl"] * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcdijkl)["ABehijkl"] * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcdijkl)["ABhfijkl"] * (*Tabij)["geKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcdijkl)["ABghijkl"] * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcdijkl)["Aefgijkl"] * (*Tabij)["BhKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["Aefhijkl"] * (*Tabij)["BgKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["Aehgijkl"] * (*Tabij)["BfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcdijkl)["Ahfgijkl"] * (*Tabij)["BeKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aghijk"] * (*Tabcijk)["BefKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Ahfijk"] * (*Tabcijk)["BegKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aheijk"] * (*Tabcijk)["BgfKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Afgijk"] * (*Tabcijk)["BehKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aegijk"] * (*Tabcijk)["BhfKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["Aefijk"] * (*Tabcijk)["BghKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghijl"] * (*Tabcijk)["BefKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Ahfijl"] * (*Tabcijk)["BegKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aheijl"] * (*Tabcijk)["BgfKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Afgijl"] * (*Tabcijk)["BehKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aegijl"] * (*Tabcijk)["BhfKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aefijl"] * (*Tabcijk)["BghKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghilk"] * (*Tabcijk)["BefKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Ahfilk"] * (*Tabcijk)["BegKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aheilk"] * (*Tabcijk)["BgfKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Afgilk"] * (*Tabcijk)["BehKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aegilk"] * (*Tabcijk)["BhfKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aefilk"] * (*Tabcijk)["BghKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aefljk"] * (*Tabcijk)["BghKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aegljk"] * (*Tabcijk)["BhfKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Afgljk"] * (*Tabcijk)["BehKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aheljk"] * (*Tabcijk)["BgfKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Ahfljk"] * (*Tabcijk)["BegKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["Aghljk"] * (*Tabcijk)["BefKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABhijk"] * (*Tabcijk)["efgKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABgijk"] * (*Tabcijk)["efhKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABfijk"] * (*Tabcijk)["ehgKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABeijk"] * (*Tabcijk)["hfgKLl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABhijl"] * (*Tabcijk)["efgKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABgijl"] * (*Tabcijk)["efhKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABfijl"] * (*Tabcijk)["ehgKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABeijl"] * (*Tabcijk)["hfgKLk"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABhilk"] * (*Tabcijk)["efgKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABgilk"] * (*Tabcijk)["efhKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABfilk"] * (*Tabcijk)["ehgKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABeilk"] * (*Tabcijk)["hfgKLj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABeljk"] * (*Tabcijk)["hfgKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABfljk"] * (*Tabcijk)["ehgKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.25) * (*Tabcijk)["ABgljk"] * (*Tabcijk)["efhKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.25) * (*Tabcijk)["ABhljk"] * (*Tabcijk)["efgKLi"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["fghIij"] * (*Tabcijk)["BCeLkl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hegIij"] * (*Tabcijk)["BCfLkl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehfIij"] * (*Tabcijk)["BCgLkl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIkl"] * (*Tabcijk)["BChLij"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["efgIij"] * (*Tabcijk)["BChLkl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["ehfIkl"] * (*Tabcijk)["BCgLij"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["hegIkl"] * (*Tabcijk)["BCfLij"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["fghIkl"] * (*Tabcijk)["BCeLij"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["fghIik"] * (*Tabcijk)["BCeLjl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hegIik"] * (*Tabcijk)["BCfLjl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehfIik"] * (*Tabcijk)["BCgLjl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIlj"] * (*Tabcijk)["BChLki"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIik"] * (*Tabcijk)["BChLjl"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehfIlj"] * (*Tabcijk)["BCgLki"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hegIlj"] * (*Tabcijk)["BCfLki"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["fghIlj"] * (*Tabcijk)["BCeLki"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["fghIil"] * (*Tabcijk)["BCeLkj"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hegIil"] * (*Tabcijk)["BCfLkj"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehfIil"] * (*Tabcijk)["BCgLkj"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIkj"] * (*Tabcijk)["BChLil"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["efgIil"] * (*Tabcijk)["BChLkj"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["ehfIkj"] * (*Tabcijk)["BCgLil"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-0.5) * (*Tabcijk)["hegIkj"] * (*Tabcijk)["BCfLil"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+0.5) * (*Tabcijk)["fghIkj"] * (*Tabcijk)["BCeLil"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AghJij"] * (*Tabcijk)["CefLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabcijk)["CegLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AheJij"] * (*Tabcijk)["CgfLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AfgJij"] * (*Tabcijk)["CehLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AegJij"] * (*Tabcijk)["ChfLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (+1.0) * (*Tabcijk)["AefJij"] * (*Tabcijk)["CghLkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJik"] * (*Tabcijk)["CefLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabcijk)["CegLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AheJik"] * (*Tabcijk)["CgfLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AfgJik"] * (*Tabcijk)["CehLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJik"] * (*Tabcijk)["ChfLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJik"] * (*Tabcijk)["CghLjl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AghJil"] * (*Tabcijk)["CefLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AhfJil"] * (*Tabcijk)["CegLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AheJil"] * (*Tabcijk)["CgfLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AfgJil"] * (*Tabcijk)["CehLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AegJil"] * (*Tabcijk)["ChfLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] +=
      (-1.0) * (*Tabcijk)["AefJil"] * (*Tabcijk)["CghLkj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKkl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKkl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKkl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKkl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKjl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKjl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKjl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKjl"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKkj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKkj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKkj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKkj"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKkl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKkl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKkl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKkl"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKlj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKlj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKlj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKlj"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKjk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKjk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKjk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKjk"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKli"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKli"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKli"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKli"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKik"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKik"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKik"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKik"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKli"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKli"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKli"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKli"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKik"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKik"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKik"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKik"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKji"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKji"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKji"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKji"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["efgKij"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["efhKij"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["ehgKij"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["hfgKij"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cghjkl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Ceflik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Ceglik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgflik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehlik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chflik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cghlik"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cghilj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cghijk"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tabcijk)["efgKkl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tabcijk)["efhKkl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tabcijk)["ehgKkl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tabcijk)["hfgKkl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tabcijk)["efgKjl"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tabcijk)["efhKjl"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tabcijk)["ehgKjl"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tabcijk)["hfgKjl"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efgKkj"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efhKkj"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tabcijk)["ehgKkj"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tabcijk)["hfgKkj"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tabcijk)["efgKil"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tabcijk)["efhKil"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tabcijk)["ehgKil"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tabcijk)["hfgKil"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efgKki"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efhKki"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tabcijk)["ehgKki"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tabcijk)["hfgKki"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efgKij"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tabcijk)["efhKij"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tabcijk)["ehgKij"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tabcijk)["hfgKij"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cefjkl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cegjkl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cgfjkl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefjkl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegjkl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfjkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cehjkl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Chfjkl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehjkl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfjkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cghjkl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Chgjkl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cefikl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cegikl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cgfikl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefikl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegikl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfikl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cehikl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Chfikl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehikl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfikl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cghikl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Chgikl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cefjil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cegjil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cgfjil"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefjil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegjil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfjil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cehjil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Chfjil"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehjil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfjil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cghjil"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Chgjil"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cefjki"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cegjki"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cgfjki"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cefjki"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cegjki"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcijk)["Cgfjki"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Cehjki"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Chfjki"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Cehjki"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcijk)["Chfjki"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcijk)["Cghjki"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcijk)["Chgjki"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tabcdijkl)["efghKLkl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tabcdijkl)["efghKLjl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tabcdijkl)["efghKLkj"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tabcdijkl)["efghKLil"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tabcdijkl)["efghKLki"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tabcdijkl)["efghKLij"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CefgLjkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CefhLjkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tabcdijkl)["CehgLjkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tabcdijkl)["ChfgLjkl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CefgLikl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CefhLikl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tabcdijkl)["CehgLikl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tabcdijkl)["ChfgLikl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CefgLjil"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CefhLjil"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tabcdijkl)["CehgLjil"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tabcdijkl)["ChfgLjil"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CefgLjki"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CefhLjki"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tabcdijkl)["CehgLjki"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tabcdijkl)["ChfgLjki"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["efghLjkl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["efghLikl"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["efghLjil"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["efghLjki"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["gI"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CDefijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["fI"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CDegijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["eI"] * (*Tai)["hJ"]
                            * (*Tabcdijkl)["CDgfijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["fI"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CDehijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["eI"] * (*Tai)["gJ"]
                            * (*Tabcdijkl)["CDhfijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["eI"] * (*Tai)["fJ"]
                            * (*Tabcdijkl)["CDghijkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["hI"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["Defgijkl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["gI"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["Defhijkl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["Dehgijkl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["BK"]
                            * (*Tabcdijkl)["Dhfgijkl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIk"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIk"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIk"] * (*Tabij)["ghJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJiC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIk"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJk"]
                            * (*Tai)["Ci"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ck"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cl"] * (*Vijka)["IJjC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJl"]
                            * (*Tai)["Cj"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cj"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJj"]
                            * (*Tai)["Cl"] * (*Vijka)["IJkC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ck"] * (*Vijka)["IJlC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKiA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["efJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["egJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"] * (*Tabij)["efJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"] * (*Tabij)["egJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["gfJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["efJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["egJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["gfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"] * (*Tabij)["ehJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["hfJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"] * (*Tabij)["ehJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"] * (*Tabij)["hfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["ghJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"] * (*Tabij)["hgJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKjA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["egJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJl"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJl"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJl"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["egJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJk"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJk"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJk"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["egJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["gfJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["egJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["gfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["hfJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["ehJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["hfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["ghJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKkA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["efJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["egJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJj"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJj"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJj"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJj"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJi"]
                            * (*Tai)["hK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJi"]
                            * (*Tai)["gK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJi"]
                            * (*Tai)["fK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJi"]
                            * (*Tai)["eK"] * (*Vijka)["JKlA"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJl"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJk"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ci"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJl"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"] * (*Tabij)["hgJk"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cj"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJi"]
                            * (*Tai)["Ck"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["efJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["egJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJj"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JhAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JgAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JfAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cl"] * (*Viabc)["JeAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bekl"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bekl"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bekl"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bekl"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bejl"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bejl"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bejl"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bejl"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bekj"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["gK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bekj"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["fK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["eK"] * (*Viabc)["KhAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bekj"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["hK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bekj"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["hK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["hK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["fK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["eK"] * (*Viabc)["KgAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["gK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["gK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["eK"] * (*Viabc)["KfAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["fK"] * (*Viabc)["KeAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chjk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgjk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfjk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Cejk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chjl"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgjl"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfjl"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Cejl"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chlk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cglk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cflk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Celk"]
                            * (*Tai)["Di"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chki"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgki"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfki"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceki"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chli"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgli"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfli"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Celi"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chij"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgij"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfij"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceij"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chij"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgij"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfij"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Ceij"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chil"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgil"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfil"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceil"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chki"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgki"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfki"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceki"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chlk"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cglk"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cflk"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Celk"]
                            * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjl"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjl"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjl"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejl"]
                            * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjk"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjk"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjk"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejk"]
                            * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKj"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKk"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKl"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKi"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKi"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKi"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKk"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKl"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKj"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKj"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJji"] * (*Tabij)["ghKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJji"] * (*Tabij)["fhKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJji"] * (*Tabij)["ehKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJji"] * (*Tabij)["gfKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJji"] * (*Tabij)["geKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJji"] * (*Tabij)["efKl"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["ghKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["fhKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["ehKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["gfKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["geKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["efKk"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKj"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKj"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKj"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKj"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKk"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKk"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKk"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKk"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKl"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKl"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKl"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKl"]
                            * (*Tai)["Di"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKi"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKi"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKi"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKi"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKi"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKi"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKi"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKi"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKi"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKi"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKi"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKi"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKk"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKk"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKk"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKk"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKl"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKl"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKl"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKl"]
                            * (*Tai)["Dj"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKj"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKj"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKj"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKj"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKj"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKj"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKj"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKj"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIji"] * (*Tabij)["BhKl"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIji"] * (*Tabij)["BgKl"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIji"] * (*Tabij)["BfKl"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIji"] * (*Tabij)["BeKl"]
                            * (*Tai)["Dk"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIij"] * (*Tabij)["BhKk"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIij"] * (*Tabij)["BgKk"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIij"] * (*Tabij)["BfKk"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIij"] * (*Tabij)["BeKk"]
                            * (*Tai)["Dl"] * (*Vijab)["IKBD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aefjkl"] * (*Tabij)["ghJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aegjkl"] * (*Tabij)["fhJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Agfjkl"] * (*Tabij)["ehJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aehjkl"] * (*Tabij)["gfJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["geJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aghjkl"] * (*Tabij)["efJK"]
                            * (*Tai)["Di"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aeflik"] * (*Tabij)["ghJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aeglik"] * (*Tabij)["fhJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agflik"] * (*Tabij)["ehJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehlik"] * (*Tabij)["gfJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahflik"] * (*Tabij)["geJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghlik"] * (*Tabij)["efJK"]
                            * (*Tai)["Dj"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aefilj"] * (*Tabij)["ghJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegilj"] * (*Tabij)["fhJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agfilj"] * (*Tabij)["ehJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehilj"] * (*Tabij)["gfJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfilj"] * (*Tabij)["geJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghilj"] * (*Tabij)["efJK"]
                            * (*Tai)["Dk"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aefijk"] * (*Tabij)["ghJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegijk"] * (*Tabij)["fhJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agfijk"] * (*Tabij)["ehJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehijk"] * (*Tabij)["gfJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfijk"] * (*Tabij)["geJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghijk"] * (*Tabij)["efJK"]
                            * (*Tai)["Dl"] * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Cgij"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Cfij"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Ceij"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Chij"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Chij"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Chij"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cfij"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Ceij"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cgij"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Cgij"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["Ceij"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJkl"] * (*Tabij)["Cfij"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Cgik"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Cfik"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Ceik"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Chik"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Chik"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Chik"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cfik"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Ceik"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cgik"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Cgik"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["Ceik"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJlj"] * (*Tabij)["Cfik"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Cgil"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Cfil"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Ceil"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Chil"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Chil"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Chil"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cfil"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Ceil"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cgil"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Cgil"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["Ceil"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJjk"] * (*Tabij)["Cfil"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Cgkj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Cfkj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Cekj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Chkj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Chkj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Chkj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cfkj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cekj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cgkj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cgkj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["Cekj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJli"] * (*Tabij)["Cfkj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Cglj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Cflj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Celj"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Chlj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Chlj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Chlj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cflj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Celj"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cglj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Cglj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["Celj"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJik"] * (*Tabij)["Cflj"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Cgkl"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Cfkl"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Cekl"]
                            * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Chkl"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Chkl"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Chkl"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cfkl"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cekl"]
                            * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cgkl"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cgkl"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["Cekl"]
                            * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJij"] * (*Tabij)["Cfkl"]
                            * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIkl"] * (*Tabij)["BCij"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIkl"] * (*Tabij)["BCij"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BCij"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BCij"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIlj"] * (*Tabij)["BCik"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIlj"] * (*Tabij)["BCik"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BCik"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BCik"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIjk"] * (*Tabij)["BCil"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIjk"] * (*Tabij)["BCil"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BCil"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BCil"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIli"] * (*Tabij)["BCkj"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIli"] * (*Tabij)["BCkj"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIli"] * (*Tabij)["BCkj"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIli"] * (*Tabij)["BCkj"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIik"] * (*Tabij)["BClj"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIik"] * (*Tabij)["BClj"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIik"] * (*Tabij)["BClj"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIik"] * (*Tabij)["BClj"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIij"] * (*Tabij)["BCkl"]
                            * (*Tai)["hL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIij"] * (*Tabij)["BCkl"]
                            * (*Tai)["gL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIij"] * (*Tabij)["BCkl"]
                            * (*Tai)["fL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIij"] * (*Tabij)["BCkl"]
                            * (*Tai)["eL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["Bhij"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["Bgij"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["Bfij"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["Beij"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["Bhik"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["Bgik"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["Bfik"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["Beik"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["Bhil"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["Bgil"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["Bfil"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["Beil"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["Bhkj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["Bgkj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["Bfkj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["Bekj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["Bhlj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["Bglj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["Bflj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["Belj"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIij"] * (*Tabij)["Bekl"]
                            * (*Tai)["CL"] * (*Vijab)["ILBC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fgKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["egKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["feKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fhKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["ehKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["hgKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["hgKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["heKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["fhKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["efKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["geKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["fgKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fgKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["egKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["feKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fhKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["ehKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["hgKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["hgKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["heKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["fhKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["efKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["geKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["fgKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fgKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["egKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["feKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fhKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["ehKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["hgKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["hgKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["heKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["fhKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["efKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["geKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["fgKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fgKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["egKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["feKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fhKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["ehKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["hgKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["hgKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["heKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["fhKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["efKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["geKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["fgKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BgKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BfKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BeKi"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BhKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BhKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BhKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BfKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BeKi"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BgKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BgKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aghjkl"] * (*Tabij)["BeKi"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahgjkl"] * (*Tabij)["BfKi"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BgKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BfKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BeKj"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BhKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BhKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BhKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BfKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BeKj"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BgKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BgKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghlik"] * (*Tabij)["BeKj"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahglik"] * (*Tabij)["BfKj"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BgKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BfKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BeKk"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BhKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BhKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BhKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BfKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BeKk"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BgKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BgKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghilj"] * (*Tabij)["BeKk"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahgilj"] * (*Tabij)["BfKk"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BgKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BfKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BeKl"]
                            * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BhKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BhKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BhKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BfKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BeKl"]
                            * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BgKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BgKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghijk"] * (*Tabij)["BeKl"]
                            * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahgijk"] * (*Tabij)["BfKl"]
                            * (*Tai)["eL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["ghJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["fhJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["ehJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["gfJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["geJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aghjkl"] * (*Tabij)["efJi"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["ghJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["fhJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["ehJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["gfJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["geJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghlik"] * (*Tabij)["efJj"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["ghJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["fhJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["ehJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["gfJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["geJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghilj"] * (*Tabij)["efJk"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["ghJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["fhJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["ehJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["gfJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["geJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghijk"] * (*Tabij)["efJl"]
                            * (*Tai)["CL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agij"] * (*Tabij)["Bekl"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afij"] * (*Tabij)["Bekl"]
                            * (*Tabij)["ghKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"]
                            * (*Tabij)["hgKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agik"] * (*Tabij)["Bejl"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afik"] * (*Tabij)["Bejl"]
                            * (*Tabij)["ghKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"]
                            * (*Tabij)["hgKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"]
                            * (*Tabij)["efKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"]
                            * (*Tabij)["egKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"]
                            * (*Tabij)["gfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agil"] * (*Tabij)["Bekj"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"]
                            * (*Tabij)["ehKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"]
                            * (*Tabij)["hfKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afil"] * (*Tabij)["Bekj"]
                            * (*Tabij)["ghKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"]
                            * (*Tabij)["hgKL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["fgJk"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["geJk"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["geJl"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["fgJl"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["heJk"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["geJl"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["fgJl"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["heJl"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["geJk"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["fgJk"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["fgJj"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["geJj"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["efJl"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["efJj"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["geJl"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["fgJl"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["heJj"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["efJl"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["ehJj"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["geJl"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["fgJl"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["efJj"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["heJl"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["geJj"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["fgJj"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"] * (*Tabij)["heJi"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"] * (*Tabij)["efJl"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"] * (*Tabij)["ehJi"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"] * (*Tabij)["geJl"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"] * (*Tabij)["fgJl"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"] * (*Tabij)["fgJi"]
                            * (*Tabij)["CeLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"] * (*Tabij)["geJi"]
                            * (*Tabij)["CfLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"] * (*Tabij)["efJl"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"] * (*Tabij)["efJi"]
                            * (*Tabij)["CgLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"] * (*Tabij)["geJl"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"] * (*Tabij)["fgJl"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"] * (*Tabij)["fgJi"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"] * (*Tabij)["geJi"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"] * (*Tabij)["efJi"]
                            * (*Tabij)["ChLl"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"] * (*Tabij)["heJl"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"] * (*Tabij)["fhJl"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"] * (*Tabij)["hgJl"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["fgJj"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["geJj"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJk"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["geJk"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["fgJk"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["heJj"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["efJk"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["geJk"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["fgJk"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["heJk"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["geJj"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJk"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["fgJj"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"] * (*Tabij)["heJi"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"] * (*Tabij)["efJk"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"] * (*Tabij)["ehJi"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"] * (*Tabij)["geJk"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"] * (*Tabij)["fgJk"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"] * (*Tabij)["fgJi"]
                            * (*Tabij)["CeLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"] * (*Tabij)["geJi"]
                            * (*Tabij)["CfLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"] * (*Tabij)["efJk"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"] * (*Tabij)["efJi"]
                            * (*Tabij)["CgLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"] * (*Tabij)["geJk"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"] * (*Tabij)["fgJk"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"] * (*Tabij)["fgJi"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"] * (*Tabij)["geJi"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"] * (*Tabij)["ehJk"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"] * (*Tabij)["efJi"]
                            * (*Tabij)["ChLk"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"] * (*Tabij)["heJk"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"] * (*Tabij)["fhJk"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"] * (*Tabij)["hgJk"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekl"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"] * (*Tabij)["hgJi"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekl"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["fgJj"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"] * (*Tabij)["ehJi"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["geJj"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"] * (*Tabij)["fhJi"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"] * (*Tabij)["heJi"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"]
                            * (*Tabij)["ChLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekl"] * (*Tabij)["fgJi"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"] * (*Tabij)["geJi"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"] * (*Tabij)["efJi"]
                            * (*Tabij)["ChLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["heJj"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["fhJj"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"] * (*Tabij)["fgJi"]
                            * (*Tabij)["CeLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"] * (*Tabij)["geJi"]
                            * (*Tabij)["CfLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"]
                            * (*Tabij)["CgLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"] * (*Tabij)["efJi"]
                            * (*Tabij)["CgLj"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["geJj"]
                            * (*Tabij)["CfLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["fgJj"]
                            * (*Tabij)["CeLi"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIk"] * (*Tabij)["efJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIk"] * (*Tabij)["egJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIk"] * (*Tabij)["gfJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIk"] * (*Tabij)["hfJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIk"] * (*Tabij)["ghJl"]
                            * (*Tabij)["CDij"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["ghIj"] * (*Tabij)["efJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["hfIj"] * (*Tabij)["egJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["heIj"] * (*Tabij)["gfJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["egIj"] * (*Tabij)["hfJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["efIj"] * (*Tabij)["ghJl"]
                            * (*Tabij)["CDik"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["ghIi"] * (*Tabij)["efJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["hfIi"] * (*Tabij)["egJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["heIi"] * (*Tabij)["gfJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["egIi"] * (*Tabij)["hfJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["efIi"] * (*Tabij)["ghJl"]
                            * (*Tabij)["CDkj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIj"] * (*Tabij)["efJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIj"] * (*Tabij)["egJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIj"] * (*Tabij)["gfJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIj"] * (*Tabij)["hfJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIj"] * (*Tabij)["ghJk"]
                            * (*Tabij)["CDil"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIi"] * (*Tabij)["efJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIi"] * (*Tabij)["egJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIi"] * (*Tabij)["gfJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIi"] * (*Tabij)["hfJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIi"] * (*Tabij)["ghJk"]
                            * (*Tabij)["CDlj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIi"] * (*Tabij)["efJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIi"] * (*Tabij)["egJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIi"] * (*Tabij)["gfJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIi"] * (*Tabij)["hfJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIi"] * (*Tabij)["ghJj"]
                            * (*Tabij)["CDkl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLkl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLkl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLkl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLkl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLjl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLjl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLjl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLjl"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLkj"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLkj"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLkj"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLkj"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLil"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLil"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLil"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLil"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLki"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLki"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLki"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLki"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tai)["hK"] * (*Tabcijk)["efgLij"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tai)["gK"] * (*Tabcijk)["efhLij"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tai)["fK"] * (*Tabcijk)["ehgLij"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                            * (*Tai)["eK"] * (*Tabcijk)["hfgLij"]
                            * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Defjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Degjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Dgfjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dehjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dhfjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                            * (*Tai)["fK"] * (*Tabcijk)["Dghjkl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Defikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Degikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Dgfikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dehikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dhfikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                            * (*Tai)["fK"] * (*Tabcijk)["Dghikl"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Defjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Degjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Dgfjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dehjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dhfjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                            * (*Tai)["fK"] * (*Tabcijk)["Dghjil"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Defjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Degjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tai)["hK"] * (*Tabcijk)["Dgfjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dehjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tai)["gK"] * (*Tabcijk)["Dhfjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                            * (*Tai)["fK"] * (*Tabcijk)["Dghjki"]
                            * (*Vijab)["JKAD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIk"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dj"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJk"]
                            * (*Tai)["Ci"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJl"]
                            * (*Tai)["Cj"] * (*Tai)["Dk"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cj"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ck"] * (*Tai)["Dl"] * (*Vijab)["IJCD"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJl"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJl"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJl"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJk"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJk"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJk"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJk"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJk"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJk"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJk"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJj"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"]
                            * (*Tai)["Ci"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"]
                            * (*Tai)["Ci"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJj"]
                            * (*Tai)["Ci"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ci"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJl"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJl"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJl"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJl"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJl"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJl"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJl"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["efJk"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["egJk"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"] * (*Tabij)["hgJk"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cj"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJi"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJi"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cj"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cj"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cj"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJl"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJl"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJl"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJl"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJl"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJk"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJk"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJk"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJk"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJk"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["egJj"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["gfJj"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["ehJj"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["hfJj"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["ghJj"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJi"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJi"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJi"]
                            * (*Tai)["Ck"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJi"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJi"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJi"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJi"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJi"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJi"]
                            * (*Tai)["Ck"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJi"]
                            * (*Tai)["Ck"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJi"]
                            * (*Tai)["Ck"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"] * (*Tabij)["efJj"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"] * (*Tabij)["egJj"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJj"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJj"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJj"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJj"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJj"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJj"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJj"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJj"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJj"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJj"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cl"] * (*Tai)["hL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJi"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJi"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJi"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJi"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJi"]
                            * (*Tai)["Cl"] * (*Tai)["gL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJi"]
                            * (*Tai)["Cl"] * (*Tai)["fL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJi"]
                            * (*Tai)["Cl"] * (*Tai)["eL"] * (*Vijab)["JLAC"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bekl"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bekl"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bejl"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bejl"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bekj"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["gK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bekj"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["fK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["eK"] * (*Tai)["hL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["fK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["eK"] * (*Tai)["gL"] * (*Vijab)["KLAB"];
  (*Rabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"]
                            * (*Tai)["eK"] * (*Tai)["fL"] * (*Vijab)["KLAB"];

  return residuum;
}
