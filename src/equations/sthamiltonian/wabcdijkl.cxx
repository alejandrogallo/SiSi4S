#include <equations/SimilarityTransformedHamiltonian.hpp>
#include <math/Complex.hpp>

using namespace sisi4s;

#define LDEBUG(eq) eq

template <typename F>
PTR(Tensor<F>) SimilarityTransformedHamiltonian<F>::getABCDIJKL() {
  if (Wabcdijkl) return Wabcdijkl;
  LOG(1, getAbbreviation()) << "Building Wabcdijkl" << std::endl;
  const int syms[] = {NS, NS, NS, NS, NS, NS, NS, NS};
  const int lens[] = {Nv, Nv, Nv, Nv, No, No, No, No};

  Wabcdijkl = NEW(Tensor<F>, 8, lens, syms, *Sisi4s::world, "Wabcdijkl");

  (*Wabcdijkl)["fuckshit"] = 0.0;

  if (Fia) {
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIjkl"] * (*Tai)["Bi"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIlik"] * (*Tai)["Bj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIilj"] * (*Tai)["Bk"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["efghIijk"] * (*Tai)["Bl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Befgijkl"] * (*Tai)["hI"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Befhijkl"] * (*Tai)["gI"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Behgijkl"] * (*Tai)["fI"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcdijkl)["Bhfgijkl"] * (*Tai)["eI"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIkl"] * (*Tabij)["Bhij"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIkl"] * (*Tabij)["Bgij"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIkl"] * (*Tabij)["Bfij"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIkl"] * (*Tabij)["Beij"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIlj"] * (*Tabij)["Bhik"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIlj"] * (*Tabij)["Bgik"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIlj"] * (*Tabij)["Bfik"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIlj"] * (*Tabij)["Beik"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIjk"] * (*Tabij)["Bhil"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIjk"] * (*Tabij)["Bgil"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIjk"] * (*Tabij)["Bfil"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIjk"] * (*Tabij)["Beil"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIli"] * (*Tabij)["Bhkj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIli"] * (*Tabij)["Bgkj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIli"] * (*Tabij)["Bfkj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIli"] * (*Tabij)["Bekj"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIik"] * (*Tabij)["Bhlj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIik"] * (*Tabij)["Bglj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIik"] * (*Tabij)["Bflj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIik"] * (*Tabij)["Belj"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["efgIij"] * (*Tabij)["Bhkl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["efhIij"] * (*Tabij)["Bgkl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["ehgIij"] * (*Tabij)["Bfkl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["hfgIij"] * (*Tabij)["Bekl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Befjkl"] * (*Tabij)["ghIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Begjkl"] * (*Tabij)["fhIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfjkl"] * (*Tabij)["ehIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Behjkl"] * (*Tabij)["gfIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfjkl"] * (*Tabij)["geIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghjkl"] * (*Tabij)["efIi"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Beflik"] * (*Tabij)["ghIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Beglik"] * (*Tabij)["fhIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgflik"] * (*Tabij)["ehIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behlik"] * (*Tabij)["gfIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhflik"] * (*Tabij)["geIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghlik"] * (*Tabij)["efIj"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Befilj"] * (*Tabij)["ghIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Begilj"] * (*Tabij)["fhIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfilj"] * (*Tabij)["ehIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behilj"] * (*Tabij)["gfIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfilj"] * (*Tabij)["geIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghilj"] * (*Tabij)["efIk"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Befijk"] * (*Tabij)["ghIl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Begijk"] * (*Tabij)["fhIl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bgfijk"] * (*Tabij)["ehIl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Behijk"] * (*Tabij)["gfIl"];
    (*Wabcdijkl)["efghijkl"] +=
        (+1.0) * (*Fia)["IB"] * (*Tabcijk)["Bhfijk"] * (*Tabij)["geIl"];
    (*Wabcdijkl)["efghijkl"] +=
        (-1.0) * (*Fia)["IB"] * (*Tabcijk)["Bghijk"] * (*Tabij)["efIl"];
  }

  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIkl"] * (*Viajk)["Ihij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIkl"] * (*Viajk)["Igij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIkl"] * (*Viajk)["Ifij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIkl"] * (*Viajk)["Ieij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIlj"] * (*Viajk)["Ihik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIlj"] * (*Viajk)["Igik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIlj"] * (*Viajk)["Ifik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIlj"] * (*Viajk)["Ieik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIjk"] * (*Viajk)["Ihil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIjk"] * (*Viajk)["Igil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIjk"] * (*Viajk)["Ifil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIjk"] * (*Viajk)["Ieil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIli"] * (*Viajk)["Ihkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIli"] * (*Viajk)["Igkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIli"] * (*Viajk)["Ifkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIli"] * (*Viajk)["Iekj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIik"] * (*Viajk)["Ihlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIik"] * (*Viajk)["Iglj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIik"] * (*Viajk)["Iflj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIik"] * (*Viajk)["Ielj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIij"] * (*Viajk)["Ihkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIij"] * (*Viajk)["Igkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIij"] * (*Viajk)["Ifkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIij"] * (*Viajk)["Iekl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefjkl"] * (*Vabic)["ghiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegjkl"] * (*Vabic)["fhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfjkl"] * (*Vabic)["ehiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehjkl"] * (*Vabic)["gfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Vabic)["geiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghjkl"] * (*Vabic)["efiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeflik"] * (*Vabic)["ghjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeglik"] * (*Vabic)["fhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agflik"] * (*Vabic)["ehjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehlik"] * (*Vabic)["gfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahflik"] * (*Vabic)["gejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghlik"] * (*Vabic)["efjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefilj"] * (*Vabic)["ghkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegilj"] * (*Vabic)["fhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfilj"] * (*Vabic)["ehkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehilj"] * (*Vabic)["gfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfilj"] * (*Vabic)["gekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghilj"] * (*Vabic)["efkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefijk"] * (*Vabic)["ghlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegijk"] * (*Vabic)["fhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfijk"] * (*Vabic)["ehlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehijk"] * (*Vabic)["gflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfijk"] * (*Vabic)["gelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghijk"] * (*Vabic)["eflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJkl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJlj"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJjk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJli"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJik"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJij"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefgJlik"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefhJlik"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AehgJlik"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefgJilj"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefhJilj"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AehgJilj"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefgJijk"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefhJijk"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AehgJijk"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABefijkl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABegijkl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABehijkl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABghijkl"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["hJ"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["gJ"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["fJ"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["eJ"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["hJ"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["gJ"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["fJ"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["eJ"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["hJ"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["gJ"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["fJ"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["eJ"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIli"] * (*Tai)["hJ"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIli"] * (*Tai)["gJ"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["fJ"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["eJ"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIik"] * (*Tai)["hJ"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIik"] * (*Tai)["gJ"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["fJ"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["eJ"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIij"] * (*Tai)["hJ"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIij"] * (*Tai)["gJ"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIij"] * (*Tai)["fJ"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIij"] * (*Tai)["eJ"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["Bj"] * (*Viajb)["IhiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["Bj"] * (*Viajb)["IgiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["Bj"] * (*Viajb)["IfiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["Bj"] * (*Viajb)["IeiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["Bk"] * (*Viajb)["IhiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["Bk"] * (*Viajb)["IgiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["Bk"] * (*Viajb)["IfiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["Bk"] * (*Viajb)["IeiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["Bl"] * (*Viajb)["IhiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["Bl"] * (*Viajb)["IgiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["Bl"] * (*Viajb)["IfiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["Bl"] * (*Viajb)["IeiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIkl"] * (*Tai)["Bi"] * (*Viajb)["IhjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIkl"] * (*Tai)["Bi"] * (*Viajb)["IgjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tai)["Bi"] * (*Viajb)["IfjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tai)["Bi"] * (*Viajb)["IejB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIlj"] * (*Tai)["Bi"] * (*Viajb)["IhkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIlj"] * (*Tai)["Bi"] * (*Viajb)["IgkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tai)["Bi"] * (*Viajb)["IfkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tai)["Bi"] * (*Viajb)["IekB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIjk"] * (*Tai)["Bi"] * (*Viajb)["IhlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIjk"] * (*Tai)["Bi"] * (*Viajb)["IglB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tai)["Bi"] * (*Viajb)["IflB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tai)["Bi"] * (*Viajb)["IelB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIli"] * (*Tai)["Bk"] * (*Viajb)["IhjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIli"] * (*Tai)["Bk"] * (*Viajb)["IgjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["Bk"] * (*Viajb)["IfjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["Bk"] * (*Viajb)["IejB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIik"] * (*Tai)["Bl"] * (*Viajb)["IhjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIik"] * (*Tai)["Bl"] * (*Viajb)["IgjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["Bl"] * (*Viajb)["IfjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["Bl"] * (*Viajb)["IejB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIli"] * (*Tai)["Bj"] * (*Viajb)["IhkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIli"] * (*Tai)["Bj"] * (*Viajb)["IgkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIli"] * (*Tai)["Bj"] * (*Viajb)["IfkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIli"] * (*Tai)["Bj"] * (*Viajb)["IekB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIik"] * (*Tai)["Bj"] * (*Viajb)["IhlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIik"] * (*Tai)["Bj"] * (*Viajb)["IglB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIik"] * (*Tai)["Bj"] * (*Viajb)["IflB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIik"] * (*Tai)["Bj"] * (*Viajb)["IelB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIji"] * (*Tai)["Bl"] * (*Viajb)["IhkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIji"] * (*Tai)["Bl"] * (*Viajb)["IgkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIji"] * (*Tai)["Bl"] * (*Viajb)["IfkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIji"] * (*Tai)["Bl"] * (*Viajb)["IekB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIij"] * (*Tai)["Bk"] * (*Viajb)["IhlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIij"] * (*Tai)["Bk"] * (*Viajb)["IglB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIij"] * (*Tai)["Bk"] * (*Viajb)["IflB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIij"] * (*Tai)["Bk"] * (*Viajb)["IelB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["gJ"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["fJ"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["eJ"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["hJ"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["hJ"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["hJ"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["fJ"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["eJ"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["gJ"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["gJ"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tai)["eJ"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahgjkl"] * (*Tai)["fJ"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["gJ"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["fJ"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agflik"] * (*Tai)["eJ"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["hJ"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["hJ"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agflik"] * (*Tai)["hJ"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["fJ"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["eJ"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["gJ"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["gJ"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghlik"] * (*Tai)["eJ"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahglik"] * (*Tai)["fJ"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["gJ"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["fJ"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["eJ"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["hJ"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["hJ"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["hJ"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["fJ"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["eJ"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["gJ"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["gJ"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghilj"] * (*Tai)["eJ"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahgilj"] * (*Tai)["fJ"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["gJ"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["fJ"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["eJ"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["hJ"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["hJ"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["hJ"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["fJ"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["eJ"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["gJ"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["gJ"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghijk"] * (*Tai)["eJ"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahgijk"] * (*Tai)["fJ"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tai)["Bi"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tai)["Bi"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tai)["Bi"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tai)["Bi"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tai)["Bi"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tai)["Bi"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeflik"] * (*Tai)["Bj"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeglik"] * (*Tai)["Bj"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agflik"] * (*Tai)["Bj"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehlik"] * (*Tai)["Bj"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahflik"] * (*Tai)["Bj"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghlik"] * (*Tai)["Bj"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefilj"] * (*Tai)["Bk"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegilj"] * (*Tai)["Bk"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfilj"] * (*Tai)["Bk"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehilj"] * (*Tai)["Bk"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tai)["Bk"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghilj"] * (*Tai)["Bk"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefijk"] * (*Tai)["Bl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegijk"] * (*Tai)["Bl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfijk"] * (*Tai)["Bl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehijk"] * (*Tai)["Bl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tai)["Bl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghijk"] * (*Tai)["Bl"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJkl"] * (*Tai)["Cj"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJlj"] * (*Tai)["Ck"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJjk"] * (*Tai)["Cl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJkl"] * (*Tai)["Ci"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJlj"] * (*Tai)["Ci"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJjk"] * (*Tai)["Ci"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJli"] * (*Tai)["Ck"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJik"] * (*Tai)["Cl"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJli"] * (*Tai)["Cj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["efghIJik"] * (*Tai)["Cj"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJji"] * (*Tai)["Cl"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["efghIJij"] * (*Tai)["Ck"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Tai)["hK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Tai)["gK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Tai)["fK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Tai)["eK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJlik"] * (*Tai)["hK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJlik"] * (*Tai)["gK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJlik"] * (*Tai)["fK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Tai)["eK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJilj"] * (*Tai)["hK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJilj"] * (*Tai)["gK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJilj"] * (*Tai)["fK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Tai)["eK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJijk"] * (*Tai)["hK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJijk"] * (*Tai)["gK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJijk"] * (*Tai)["fK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Tai)["eK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["efghIjkl"] * (*Tai)["BK"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["efghIlik"] * (*Tai)["BK"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["efghIilj"] * (*Tai)["BK"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["efghIijk"] * (*Tai)["BK"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefhJjkl"] * (*Tai)["Ci"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AehgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AhfgJjkl"] * (*Tai)["Ci"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJlik"] * (*Tai)["Cj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJlik"] * (*Tai)["Cj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJlik"] * (*Tai)["Cj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJlik"] * (*Tai)["Cj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJilj"] * (*Tai)["Ck"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJilj"] * (*Tai)["Ck"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJilj"] * (*Tai)["Ck"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJilj"] * (*Tai)["Ck"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["AefgJijk"] * (*Tai)["Cl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AefhJijk"] * (*Tai)["Cl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AehgJijk"] * (*Tai)["Cl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["AhfgJijk"] * (*Tai)["Cl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABefijkl"] * (*Tai)["gK"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABegijkl"] * (*Tai)["fK"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Tai)["eK"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABefijkl"] * (*Tai)["hK"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABegijkl"] * (*Tai)["hK"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABgfijkl"] * (*Tai)["hK"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABehijkl"] * (*Tai)["fK"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Tai)["eK"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABehijkl"] * (*Tai)["gK"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcdijkl)["ABhfijkl"] * (*Tai)["gK"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABghijkl"] * (*Tai)["eK"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcdijkl)["ABhgijkl"] * (*Tai)["fK"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcdijkl)["Aefgijkl"] * (*Tai)["BK"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["Aefhijkl"] * (*Tai)["BK"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["Aehgijkl"] * (*Tai)["BK"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcdijkl)["Ahfgijkl"] * (*Tai)["BK"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["ghIk"] * (*Tabij)["efJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["hfIk"] * (*Tabij)["egJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["heIk"] * (*Tabij)["gfJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["fgIk"] * (*Tabij)["ehJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["egIk"] * (*Tabij)["hfJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["efIk"] * (*Tabij)["ghJl"] * (*Vijkl)["IJij"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJl"] * (*Vijkl)["IJik"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["ghIj"] * (*Tabij)["efJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["hfIj"] * (*Tabij)["egJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["heIj"] * (*Tabij)["gfJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["fgIj"] * (*Tabij)["ehJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["egIj"] * (*Tabij)["hfJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["efIj"] * (*Tabij)["ghJk"] * (*Vijkl)["IJil"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJl"] * (*Vijkl)["IJkj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJk"] * (*Vijkl)["IJlj"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["ghIi"] * (*Tabij)["efJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["hfIi"] * (*Tabij)["egJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["heIi"] * (*Tabij)["gfJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["fgIi"] * (*Tabij)["ehJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["egIi"] * (*Tabij)["hfJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["efIi"] * (*Tabij)["ghJj"] * (*Vijkl)["IJkl"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJl"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJl"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJl"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJl"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJl"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJl"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJl"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJl"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJl"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJl"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJl"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJl"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJk"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJk"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJk"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJk"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJk"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJk"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJk"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJk"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJk"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJk"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJk"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJk"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJj"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJj"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJj"] * (*Viajb)["JhiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJj"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJj"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJj"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJj"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJj"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJj"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJj"] * (*Viajb)["JgiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJj"] * (*Viajb)["JfiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJj"] * (*Viajb)["JeiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agki"] * (*Tabij)["efJl"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afki"] * (*Tabij)["egJl"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJl"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJl"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJl"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJl"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJl"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJl"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJl"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJl"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJl"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJl"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agli"] * (*Tabij)["efJk"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afli"] * (*Tabij)["egJk"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeli"] * (*Tabij)["gfJk"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahli"] * (*Tabij)["efJk"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["egJk"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahli"] * (*Tabij)["gfJk"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afli"] * (*Tabij)["ehJk"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["hfJk"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agli"] * (*Tabij)["ehJk"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agli"] * (*Tabij)["hfJk"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeli"] * (*Tabij)["ghJk"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afli"] * (*Tabij)["hgJk"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aelk"] * (*Tabij)["gfJi"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aflk"] * (*Tabij)["egJi"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aglk"] * (*Tabij)["efJi"] * (*Viajb)["JhjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["hfJi"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["ehJi"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aelk"] * (*Tabij)["ghJi"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aflk"] * (*Tabij)["hgJi"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["ehJi"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aglk"] * (*Tabij)["hfJi"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahlk"] * (*Tabij)["efJi"] * (*Viajb)["JgjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["egJi"] * (*Viajb)["JfjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahlk"] * (*Tabij)["gfJi"] * (*Viajb)["JejA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agij"] * (*Tabij)["efJl"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afij"] * (*Tabij)["egJl"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJl"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJl"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJl"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJl"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJl"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJl"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJl"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJl"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJl"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJl"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agij"] * (*Tabij)["efJk"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afij"] * (*Tabij)["egJk"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["gfJk"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["efJk"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["egJk"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["gfJk"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afij"] * (*Tabij)["ehJk"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["hfJk"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agij"] * (*Tabij)["ehJk"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agij"] * (*Tabij)["hfJk"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["ghJk"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afij"] * (*Tabij)["hgJk"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agil"] * (*Tabij)["efJj"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afil"] * (*Tabij)["egJj"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["gfJj"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["efJj"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["egJj"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["gfJj"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afil"] * (*Tabij)["ehJj"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["hfJj"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agil"] * (*Tabij)["ehJj"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agil"] * (*Tabij)["hfJj"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["ghJj"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afil"] * (*Tabij)["hgJj"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejl"] * (*Tabij)["gfJi"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjl"] * (*Tabij)["egJi"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjl"] * (*Tabij)["efJi"] * (*Viajb)["JhkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["hfJi"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["ehJi"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejl"] * (*Tabij)["ghJi"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjl"] * (*Tabij)["hgJi"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["ehJi"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjl"] * (*Tabij)["hfJi"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjl"] * (*Tabij)["efJi"] * (*Viajb)["JgkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["egJi"] * (*Viajb)["JfkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjl"] * (*Tabij)["gfJi"] * (*Viajb)["JekA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agki"] * (*Tabij)["efJj"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afki"] * (*Tabij)["egJj"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeki"] * (*Tabij)["gfJj"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahki"] * (*Tabij)["efJj"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["egJj"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahki"] * (*Tabij)["gfJj"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afki"] * (*Tabij)["ehJj"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["hfJj"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agki"] * (*Tabij)["ehJj"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agki"] * (*Tabij)["hfJj"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeki"] * (*Tabij)["ghJj"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afki"] * (*Tabij)["hgJj"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aejk"] * (*Tabij)["gfJi"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afjk"] * (*Tabij)["egJi"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agjk"] * (*Tabij)["efJi"] * (*Viajb)["JhlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["hfJi"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["ehJi"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aejk"] * (*Tabij)["ghJi"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afjk"] * (*Tabij)["hgJi"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["ehJi"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agjk"] * (*Tabij)["hfJi"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahjk"] * (*Tabij)["efJi"] * (*Viajb)["JglA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["egJi"] * (*Viajb)["JflA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahjk"] * (*Tabij)["gfJi"] * (*Viajb)["JelA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afij"] * (*Tabij)["Bekl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bfkl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agij"] * (*Tabij)["Bekl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bfkl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bgkl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bgkl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bekl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bfkl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahij"] * (*Tabij)["Bgkl"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeij"] * (*Tabij)["Bhkl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afij"] * (*Tabij)["Bhkl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agij"] * (*Tabij)["Bhkl"] * (*Vabcd)["feAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afik"] * (*Tabij)["Bejl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bfjl"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agik"] * (*Tabij)["Bejl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bfjl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bgjl"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bgjl"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bejl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bfjl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahik"] * (*Tabij)["Bgjl"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeik"] * (*Tabij)["Bhjl"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afik"] * (*Tabij)["Bhjl"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agik"] * (*Tabij)["Bhjl"] * (*Vabcd)["feAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Afil"] * (*Tabij)["Bekj"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bfkj"] * (*Vabcd)["ghAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Agil"] * (*Tabij)["Bekj"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bfkj"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bgkj"] * (*Vabcd)["fhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bgkj"] * (*Vabcd)["ehAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bekj"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bfkj"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Ahil"] * (*Tabij)["Bgkj"] * (*Vabcd)["efAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Aeil"] * (*Tabij)["Bhkj"] * (*Vabcd)["gfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabij)["Afil"] * (*Tabij)["Bhkj"] * (*Vabcd)["geAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabij)["Agil"] * (*Tabij)["Bhkj"] * (*Vabcd)["feAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chjk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgjk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfjk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Cejk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chjl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgjl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfjl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Cejl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chlk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cglk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cflk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Celk"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chki"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgki"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfki"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceki"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chli"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgli"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfli"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Celi"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIJl"] * (*Tabij)["Chij"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIJl"] * (*Tabij)["Cgij"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIJl"] * (*Tabij)["Cfij"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIJl"] * (*Tabij)["Ceij"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJk"] * (*Tabij)["Chij"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJk"] * (*Tabij)["Cgij"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJk"] * (*Tabij)["Cfij"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJk"] * (*Tabij)["Ceij"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chil"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgil"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfil"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceil"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJj"] * (*Tabij)["Chki"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJj"] * (*Tabij)["Cgki"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJj"] * (*Tabij)["Cfki"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJj"] * (*Tabij)["Ceki"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chlk"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cglk"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cflk"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Celk"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjl"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjl"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjl"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejl"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efgIJi"] * (*Tabij)["Chjk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efhIJi"] * (*Tabij)["Cgjk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ehgIJi"] * (*Tabij)["Cfjk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["hfgIJi"] * (*Tabij)["Cejk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKk"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["ghKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["fhKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["ehKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["gfKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["geKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["efKi"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["ghKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["fhKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["ehKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["gfKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["geKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["efKi"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["ghKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["fhKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["ehKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["gfKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["geKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["efKi"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKk"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["ghKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["fhKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["ehKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["gfKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["geKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["efKj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["ghKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["fhKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["ehKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["gfKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["geKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["efKj"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJji"] * (*Tabij)["ghKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJji"] * (*Tabij)["fhKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJji"] * (*Tabij)["ehKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJji"] * (*Tabij)["gfKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJji"] * (*Tabij)["geKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJji"] * (*Tabij)["efKl"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["ghKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["fhKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["ehKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["gfKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["geKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["efKk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKj"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKj"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKj"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKj"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKk"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKk"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKk"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKk"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKl"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKl"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKl"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKl"] * (*Vijka)["IKiB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIkl"] * (*Tabij)["BhKi"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIkl"] * (*Tabij)["BgKi"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BfKi"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BeKi"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIlj"] * (*Tabij)["BhKi"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIlj"] * (*Tabij)["BgKi"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BfKi"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BeKi"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIjk"] * (*Tabij)["BhKi"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIjk"] * (*Tabij)["BgKi"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BfKi"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BeKi"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKk"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKk"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKk"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKk"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKl"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKl"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKl"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKl"] * (*Vijka)["IKjB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIli"] * (*Tabij)["BhKj"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIli"] * (*Tabij)["BgKj"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIli"] * (*Tabij)["BfKj"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIli"] * (*Tabij)["BeKj"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efgIik"] * (*Tabij)["BhKj"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efhIik"] * (*Tabij)["BgKj"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["ehgIik"] * (*Tabij)["BfKj"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["hfgIik"] * (*Tabij)["BeKj"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIji"] * (*Tabij)["BhKl"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIji"] * (*Tabij)["BgKl"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIji"] * (*Tabij)["BfKl"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIji"] * (*Tabij)["BeKl"] * (*Vijka)["IKkB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["efgIij"] * (*Tabij)["BhKk"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["efhIij"] * (*Tabij)["BgKk"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["ehgIij"] * (*Tabij)["BfKk"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["hfgIij"] * (*Tabij)["BeKk"] * (*Vijka)["IKlB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aefjkl"] * (*Tabij)["ghJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aegjkl"] * (*Tabij)["fhJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Agfjkl"] * (*Tabij)["ehJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aehjkl"] * (*Tabij)["gfJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["geJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aghjkl"] * (*Tabij)["efJK"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aeflik"] * (*Tabij)["ghJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aeglik"] * (*Tabij)["fhJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Agflik"] * (*Tabij)["ehJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aehlik"] * (*Tabij)["gfJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Ahflik"] * (*Tabij)["geJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aghlik"] * (*Tabij)["efJK"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aefilj"] * (*Tabij)["ghJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aegilj"] * (*Tabij)["fhJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Agfilj"] * (*Tabij)["ehJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aehilj"] * (*Tabij)["gfJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Ahfilj"] * (*Tabij)["geJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aghilj"] * (*Tabij)["efJK"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aefijk"] * (*Tabij)["ghJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aegijk"] * (*Tabij)["fhJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Agfijk"] * (*Tabij)["ehJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Aehijk"] * (*Tabij)["gfJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["Ahfijk"] * (*Tabij)["geJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["Aghijk"] * (*Tabij)["efJK"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJkl"] * (*Tabij)["Chij"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJkl"] * (*Tabij)["Chij"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJkl"] * (*Tabij)["Chij"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJkl"] * (*Tabij)["Cgij"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJkl"] * (*Tabij)["Ceij"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJkl"] * (*Tabij)["Cfij"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJlj"] * (*Tabij)["Chik"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJlj"] * (*Tabij)["Chik"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJlj"] * (*Tabij)["Chik"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJlj"] * (*Tabij)["Cgik"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJlj"] * (*Tabij)["Ceik"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJlj"] * (*Tabij)["Cfik"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJjk"] * (*Tabij)["Chil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJjk"] * (*Tabij)["Chil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJjk"] * (*Tabij)["Chil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJjk"] * (*Tabij)["Cgil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJjk"] * (*Tabij)["Ceil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJjk"] * (*Tabij)["Cfil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Cekj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJli"] * (*Tabij)["Chkj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJli"] * (*Tabij)["Chkj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJli"] * (*Tabij)["Chkj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cekj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJli"] * (*Tabij)["Cgkj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJli"] * (*Tabij)["Cekj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJli"] * (*Tabij)["Cfkj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Cglj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Cflj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Celj"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJik"] * (*Tabij)["Chlj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJik"] * (*Tabij)["Chlj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJik"] * (*Tabij)["Chlj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cflj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Celj"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJik"] * (*Tabij)["Cglj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJik"] * (*Tabij)["Cglj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJik"] * (*Tabij)["Celj"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJik"] * (*Tabij)["Cflj"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Cekl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AefJij"] * (*Tabij)["Chkl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AegJij"] * (*Tabij)["Chkl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AgfJij"] * (*Tabij)["Chkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cekl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AehJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["AhfJij"] * (*Tabij)["Cgkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AghJij"] * (*Tabij)["Cekl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["AhgJij"] * (*Tabij)["Cfkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIkl"] * (*Tabij)["BCij"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIkl"] * (*Tabij)["BCij"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIlj"] * (*Tabij)["BCik"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIlj"] * (*Tabij)["BCik"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIjk"] * (*Tabij)["BCil"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIjk"] * (*Tabij)["BCil"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIli"] * (*Tabij)["BCkj"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIli"] * (*Tabij)["BCkj"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIik"] * (*Tabij)["BClj"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIik"] * (*Tabij)["BClj"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIik"] * (*Tabij)["BClj"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIik"] * (*Tabij)["BClj"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["efgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IhBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["efhIij"] * (*Tabij)["BCkl"] * (*Viabc)["IgBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ehgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IfBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["hfgIij"] * (*Tabij)["BCkl"] * (*Viabc)["IeBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fgKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["egKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["feKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["fhKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["ehKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABejkl"] * (*Tabij)["hgKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABfjkl"] * (*Tabij)["hgKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["heKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABgjkl"] * (*Tabij)["fhKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["efKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["geKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABhjkl"] * (*Tabij)["fgKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fgKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["egKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["feKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["fhKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["ehKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABelik"] * (*Tabij)["hgKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABflik"] * (*Tabij)["hgKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["heKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABglik"] * (*Tabij)["fhKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["efKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["geKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhlik"] * (*Tabij)["fgKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fgKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["egKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["feKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["fhKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["ehKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABeilj"] * (*Tabij)["hgKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABfilj"] * (*Tabij)["hgKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["heKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABgilj"] * (*Tabij)["fhKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["efKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["geKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhilj"] * (*Tabij)["fgKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fgKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["egKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["feKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["fhKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["ehKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABeijk"] * (*Tabij)["hgKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABfijk"] * (*Tabij)["hgKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["heKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+0.5) * (*Tabcijk)["ABgijk"] * (*Tabij)["fhKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["efKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["geKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-0.5) * (*Tabcijk)["ABhijk"] * (*Tabij)["fgKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfjkl"] * (*Tabij)["BhKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfjkl"] * (*Tabij)["BgKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aghjkl"] * (*Tabij)["BeKi"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahgjkl"] * (*Tabij)["BfKi"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BgKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BfKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BeKj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aeflik"] * (*Tabij)["BhKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aeglik"] * (*Tabij)["BhKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agflik"] * (*Tabij)["BhKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BfKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BeKj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehlik"] * (*Tabij)["BgKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahflik"] * (*Tabij)["BgKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghlik"] * (*Tabij)["BeKj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahglik"] * (*Tabij)["BfKj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BgKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BfKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BeKk"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefilj"] * (*Tabij)["BhKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegilj"] * (*Tabij)["BhKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfilj"] * (*Tabij)["BhKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BfKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BeKk"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehilj"] * (*Tabij)["BgKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfilj"] * (*Tabij)["BgKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghilj"] * (*Tabij)["BeKk"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahgilj"] * (*Tabij)["BfKk"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BgKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BfKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BeKl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aefijk"] * (*Tabij)["BhKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aegijk"] * (*Tabij)["BhKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Agfijk"] * (*Tabij)["BhKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BfKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BeKl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Aehijk"] * (*Tabij)["BgKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (-1.0) * (*Tabcijk)["Ahfijk"] * (*Tabij)["BgKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Aghijk"] * (*Tabij)["BeKl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] +=
         (+1.0) * (*Tabcijk)["Ahgijk"] * (*Tabij)["BfKl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKkl"]
                                   * (*Tabij)["Dhij"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKkl"]
                                   * (*Tabij)["Dgij"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKkl"]
                                   * (*Tabij)["Dfij"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKkl"]
                                   * (*Tabij)["Deij"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKlj"]
                                   * (*Tabij)["Dhik"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKlj"]
                                   * (*Tabij)["Dgik"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKlj"]
                                   * (*Tabij)["Dfik"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKlj"]
                                   * (*Tabij)["Deik"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKjk"]
                                   * (*Tabij)["Dhil"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKjk"]
                                   * (*Tabij)["Dgil"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKjk"]
                                   * (*Tabij)["Dfil"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKjk"]
                                   * (*Tabij)["Deil"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKli"]
                                   * (*Tabij)["Dhkj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKli"]
                                   * (*Tabij)["Dgkj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKli"]
                                   * (*Tabij)["Dfkj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKli"]
                                   * (*Tabij)["Dekj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKik"]
                                   * (*Tabij)["Dhlj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKik"]
                                   * (*Tabij)["Dglj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKik"]
                                   * (*Tabij)["Dflj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKik"]
                                   * (*Tabij)["Delj"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["AefgJKij"]
                                   * (*Tabij)["Dhkl"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AefhJKij"]
                                   * (*Tabij)["Dgkl"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AehgJKij"]
                                   * (*Tabij)["Dfkl"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["AhfgJKij"]
                                   * (*Tabij)["Dekl"] * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJkl"]
                                   * (*Tabij)["CDij"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJlj"]
                                   * (*Tabij)["CDik"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJjk"]
                                   * (*Tabij)["CDil"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJli"]
                                   * (*Tabij)["CDkj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJik"]
                                   * (*Tabij)["CDlj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["efghIJij"]
                                   * (*Tabij)["CDkl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABefKjkl"]
                                   * (*Tabij)["ghLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABegKjkl"]
                                   * (*Tabij)["fhLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABgfKjkl"]
                                   * (*Tabij)["ehLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABehKjkl"]
                                   * (*Tabij)["gfLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABhfKjkl"]
                                   * (*Tabij)["geLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABghKjkl"]
                                   * (*Tabij)["efLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABefKlik"]
                                   * (*Tabij)["ghLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABegKlik"]
                                   * (*Tabij)["fhLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABgfKlik"]
                                   * (*Tabij)["ehLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABehKlik"]
                                   * (*Tabij)["gfLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABhfKlik"]
                                   * (*Tabij)["geLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABghKlik"]
                                   * (*Tabij)["efLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABefKilj"]
                                   * (*Tabij)["ghLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABegKilj"]
                                   * (*Tabij)["fhLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABgfKilj"]
                                   * (*Tabij)["ehLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABehKilj"]
                                   * (*Tabij)["gfLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABhfKilj"]
                                   * (*Tabij)["geLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABghKilj"]
                                   * (*Tabij)["efLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABefKijk"]
                                   * (*Tabij)["ghLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABegKijk"]
                                   * (*Tabij)["fhLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABgfKijk"]
                                   * (*Tabij)["ehLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABehKijk"]
                                   * (*Tabij)["gfLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["ABhfKijk"]
                                   * (*Tabij)["geLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["ABghKijk"]
                                   * (*Tabij)["efLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AefgJjkl"]
                                   * (*Tabij)["ChLi"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AefhJjkl"]
                                   * (*Tabij)["CgLi"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AehgJjkl"]
                                   * (*Tabij)["CfLi"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AhfgJjkl"]
                                   * (*Tabij)["CeLi"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AefgJlik"]
                                   * (*Tabij)["ChLj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AefhJlik"]
                                   * (*Tabij)["CgLj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AehgJlik"]
                                   * (*Tabij)["CfLj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AhfgJlik"]
                                   * (*Tabij)["CeLj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AefgJilj"]
                                   * (*Tabij)["ChLk"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AefhJilj"]
                                   * (*Tabij)["CgLk"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AehgJilj"]
                                   * (*Tabij)["CfLk"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AhfgJilj"]
                                   * (*Tabij)["CeLk"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcdijkl)["AefgJijk"]
                                   * (*Tabij)["ChLl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AefhJijk"]
                                   * (*Tabij)["CgLl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AehgJijk"]
                                   * (*Tabij)["CfLl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcdijkl)["AhfgJijk"]
                                   * (*Tabij)["CeLl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["efghIjkl"]
                                   * (*Tabij)["BCLi"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["efghIlik"]
                                   * (*Tabij)["BCLj"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["efghIilj"]
                                   * (*Tabij)["BCLk"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["efghIijk"]
                                   * (*Tabij)["BCLl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["ABefijkl"]
                                   * (*Tabij)["ghKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcdijkl)["ABegijkl"]
                                   * (*Tabij)["fhKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcdijkl)["ABgfijkl"]
                                   * (*Tabij)["ehKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcdijkl)["ABehijkl"]
                                   * (*Tabij)["gfKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcdijkl)["ABhfijkl"]
                                   * (*Tabij)["geKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcdijkl)["ABghijkl"]
                                   * (*Tabij)["efKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcdijkl)["Aefgijkl"]
                                   * (*Tabij)["BhKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["Aefhijkl"]
                                   * (*Tabij)["BgKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["Aehgijkl"]
                                   * (*Tabij)["BfKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcdijkl)["Ahfgijkl"]
                                   * (*Tabij)["BeKL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghijk"]
                                   * (*Tabcijk)["BefKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Ahfijk"]
                                   * (*Tabcijk)["BegKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aheijk"]
                                   * (*Tabcijk)["BgfKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Afgijk"]
                                   * (*Tabcijk)["BehKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aegijk"]
                                   * (*Tabcijk)["BhfKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aefijk"]
                                   * (*Tabcijk)["BghKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aghijl"]
                                   * (*Tabcijk)["BefKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfijl"]
                                   * (*Tabcijk)["BegKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aheijl"]
                                   * (*Tabcijk)["BgfKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Afgijl"]
                                   * (*Tabcijk)["BehKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegijl"]
                                   * (*Tabcijk)["BhfKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aefijl"]
                                   * (*Tabcijk)["BghKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aghilk"]
                                   * (*Tabcijk)["BefKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfilk"]
                                   * (*Tabcijk)["BegKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aheilk"]
                                   * (*Tabcijk)["BgfKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Afgilk"]
                                   * (*Tabcijk)["BehKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegilk"]
                                   * (*Tabcijk)["BhfKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aefilk"]
                                   * (*Tabcijk)["BghKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aefljk"]
                                   * (*Tabcijk)["BghKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegljk"]
                                   * (*Tabcijk)["BhfKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Afgljk"]
                                   * (*Tabcijk)["BehKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aheljk"]
                                   * (*Tabcijk)["BgfKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfljk"]
                                   * (*Tabcijk)["BegKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aghljk"]
                                   * (*Tabcijk)["BefKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABhijk"]
                                   * (*Tabcijk)["efgKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABgijk"]
                                   * (*Tabcijk)["efhKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABfijk"]
                                   * (*Tabcijk)["ehgKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABeijk"]
                                   * (*Tabcijk)["hfgKLl"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABhijl"]
                                   * (*Tabcijk)["efgKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABgijl"]
                                   * (*Tabcijk)["efhKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABfijl"]
                                   * (*Tabcijk)["ehgKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABeijl"]
                                   * (*Tabcijk)["hfgKLk"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABhilk"]
                                   * (*Tabcijk)["efgKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABgilk"]
                                   * (*Tabcijk)["efhKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABfilk"]
                                   * (*Tabcijk)["ehgKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABeilk"]
                                   * (*Tabcijk)["hfgKLj"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABeljk"]
                                   * (*Tabcijk)["hfgKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABfljk"]
                                   * (*Tabcijk)["ehgKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.25) * (*Tabcijk)["ABgljk"]
                                   * (*Tabcijk)["efhKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.25) * (*Tabcijk)["ABhljk"]
                                   * (*Tabcijk)["efgKLi"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["fghIij"]
                                   * (*Tabcijk)["BCeLkl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hegIij"]
                                   * (*Tabcijk)["BCfLkl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehfIij"]
                                   * (*Tabcijk)["BCgLkl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIkl"]
                                   * (*Tabcijk)["BChLij"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIij"]
                                   * (*Tabcijk)["BChLkl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehfIkl"]
                                   * (*Tabcijk)["BCgLij"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hegIkl"]
                                   * (*Tabcijk)["BCfLij"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["fghIkl"]
                                   * (*Tabcijk)["BCeLij"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["fghIik"]
                                   * (*Tabcijk)["BCeLjl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hegIik"]
                                   * (*Tabcijk)["BCfLjl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehfIik"]
                                   * (*Tabcijk)["BCgLjl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIlj"]
                                   * (*Tabcijk)["BChLki"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIik"]
                                   * (*Tabcijk)["BChLjl"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehfIlj"]
                                   * (*Tabcijk)["BCgLki"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hegIlj"]
                                   * (*Tabcijk)["BCfLki"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["fghIlj"]
                                   * (*Tabcijk)["BCeLki"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["fghIil"]
                                   * (*Tabcijk)["BCeLkj"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hegIil"]
                                   * (*Tabcijk)["BCfLkj"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehfIil"]
                                   * (*Tabcijk)["BCgLkj"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIkj"]
                                   * (*Tabcijk)["BChLil"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIil"]
                                   * (*Tabcijk)["BChLkj"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehfIkj"]
                                   * (*Tabcijk)["BCgLil"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hegIkj"]
                                   * (*Tabcijk)["BCfLil"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["fghIkj"]
                                   * (*Tabcijk)["BCeLil"] * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJij"]
                                   * (*Tabcijk)["CefLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJij"]
                                   * (*Tabcijk)["CegLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AheJij"]
                                   * (*Tabcijk)["CgfLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AfgJij"]
                                   * (*Tabcijk)["CehLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJij"]
                                   * (*Tabcijk)["ChfLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJij"]
                                   * (*Tabcijk)["CghLkl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJik"]
                                   * (*Tabcijk)["CefLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJik"]
                                   * (*Tabcijk)["CegLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AheJik"]
                                   * (*Tabcijk)["CgfLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AfgJik"]
                                   * (*Tabcijk)["CehLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJik"]
                                   * (*Tabcijk)["ChfLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJik"]
                                   * (*Tabcijk)["CghLjl"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJil"]
                                   * (*Tabcijk)["CefLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJil"]
                                   * (*Tabcijk)["CegLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AheJil"]
                                   * (*Tabcijk)["CgfLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AfgJil"]
                                   * (*Tabcijk)["CehLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJil"]
                                   * (*Tabcijk)["ChfLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJil"]
                                   * (*Tabcijk)["CghLkj"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKkl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKkl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKkl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKkl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKjl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKjl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKjl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKjl"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKkj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKkj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKkj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKkj"] * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKkl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKkl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKkl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKkl"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKlj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKlj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKlj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKlj"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKjk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKjk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKjk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKjk"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKli"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKli"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKli"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKli"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKik"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKik"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKik"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKik"] * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKli"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKli"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKli"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKli"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKik"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKik"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKik"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKik"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKji"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKji"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKji"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKji"] * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["efgKij"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["efhKij"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["ehgKij"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["hfgKij"] * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cghjkl"] * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Ceflik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Ceglik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgflik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehlik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chflik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cghlik"] * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cghilj"] * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["gI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["fI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["eI"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cghijk"] * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tabcijk)["efgKkl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tabcijk)["efhKkl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tabcijk)["ehgKkl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tabcijk)["hfgKkl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["efgKjl"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["efhKjl"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["ehgKjl"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["hfgKjl"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efgKkj"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efhKkj"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["ehgKkj"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["hfgKkj"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["efgKil"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["efhKil"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["ehgKil"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tabcijk)["hfgKil"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efgKki"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efhKki"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["ehgKki"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["hfgKki"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efgKij"] * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["efhKij"] * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["ehgKij"] * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tabcijk)["hfgKij"] * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cefjkl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cegjkl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cgfjkl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefjkl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegjkl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfjkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cehjkl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Chfjkl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehjkl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfjkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cghjkl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Chgjkl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cefikl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cegikl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cgfikl"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefikl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegikl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfikl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cehikl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Chfikl"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehikl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfikl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cghikl"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Chgikl"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cefjil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cegjil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cgfjil"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefjil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegjil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfjil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cehjil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Chfjil"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehjil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfjil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cghjil"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Chgjil"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cefjki"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cegjki"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cgfjki"] * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cefjki"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cegjki"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcijk)["Cgfjki"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Cehjki"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Chfjki"] * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Cehjki"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcijk)["Chfjki"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcijk)["Cghjki"] * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcijk)["Chgjki"] * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tabcdijkl)["efghKLkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tabcdijkl)["efghKLjl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tabcdijkl)["efghKLkj"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tabcdijkl)["efghKLil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tabcdijkl)["efghKLki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tabcdijkl)["efghKLij"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CefgLjkl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CefhLjkl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tabcdijkl)["CehgLjkl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tabcdijkl)["ChfgLjkl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CefgLikl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CefhLikl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tabcdijkl)["CehgLikl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tabcdijkl)["ChfgLikl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CefgLjil"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CefhLjil"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tabcdijkl)["CehgLjil"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tabcdijkl)["ChfgLjil"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CefgLjki"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CefhLjki"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tabcdijkl)["CehgLjki"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tabcdijkl)["ChfgLjki"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["efghLjkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["efghLikl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["efghLjil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["efghLjki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["gI"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CDefijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["fI"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CDegijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tai)["eI"] * (*Tai)["hJ"]
                                   * (*Tabcdijkl)["CDgfijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["fI"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CDehijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["eI"] * (*Tai)["gJ"]
                                   * (*Tabcdijkl)["CDhfijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tai)["eI"] * (*Tai)["fJ"]
                                   * (*Tabcdijkl)["CDghijkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["hI"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["Defgijkl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["gI"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["Defhijkl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["fI"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["Dehgijkl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["eI"] * (*Tai)["BK"]
                                   * (*Tabcdijkl)["Dhfgijkl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJiC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Ci"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJjC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cj"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Cl"]
                                   * (*Vijka)["IJkC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ck"]
                                   * (*Vijka)["IJlC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKiA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["efJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["egJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["efJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["egJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["ghJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["hgJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKjA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJl"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJl"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJk"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJk"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["efJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["egJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["efJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["egJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["ghJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["hgJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKkA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJj"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJj"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["hK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["gK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["fK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["eK"]
                                   * (*Vijka)["JKlA"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ci"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cj"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Ck"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JhAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JgAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JfAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cl"]
                                   * (*Viabc)["JeAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["gK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["fK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KhAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["hK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["hK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["hK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["fK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KgAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["gK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["gK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["eK"]
                                   * (*Viabc)["KfAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["fK"]
                                   * (*Viabc)["KeAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"]
                                   * (*Tabij)["Chjk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"]
                                   * (*Tabij)["Cgjk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"]
                                   * (*Tabij)["Cfjk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"]
                                   * (*Tabij)["Cejk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"]
                                   * (*Tabij)["Chjl"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"]
                                   * (*Tabij)["Cgjl"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"]
                                   * (*Tabij)["Cfjl"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"]
                                   * (*Tabij)["Cejl"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"]
                                   * (*Tabij)["Chlk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"]
                                   * (*Tabij)["Cglk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"]
                                   * (*Tabij)["Cflk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"]
                                   * (*Tabij)["Celk"] * (*Tai)["Di"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"]
                                   * (*Tabij)["Chki"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"]
                                   * (*Tabij)["Cgki"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"]
                                   * (*Tabij)["Cfki"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"]
                                   * (*Tabij)["Ceki"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"]
                                   * (*Tabij)["Chli"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"]
                                   * (*Tabij)["Cgli"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"]
                                   * (*Tabij)["Cfli"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"]
                                   * (*Tabij)["Celi"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJl"]
                                   * (*Tabij)["Chij"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJl"]
                                   * (*Tabij)["Cgij"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJl"]
                                   * (*Tabij)["Cfij"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJl"]
                                   * (*Tabij)["Ceij"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJk"]
                                   * (*Tabij)["Chij"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJk"]
                                   * (*Tabij)["Cgij"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJk"]
                                   * (*Tabij)["Cfij"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJk"]
                                   * (*Tabij)["Ceij"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"]
                                   * (*Tabij)["Chil"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"]
                                   * (*Tabij)["Cgil"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"]
                                   * (*Tabij)["Cfil"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"]
                                   * (*Tabij)["Ceil"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJj"]
                                   * (*Tabij)["Chki"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJj"]
                                   * (*Tabij)["Cgki"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJj"]
                                   * (*Tabij)["Cfki"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJj"]
                                   * (*Tabij)["Ceki"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJi"]
                                   * (*Tabij)["Chlk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJi"]
                                   * (*Tabij)["Cglk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJi"]
                                   * (*Tabij)["Cflk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJi"]
                                   * (*Tabij)["Celk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIJi"]
                                   * (*Tabij)["Chjl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIJi"]
                                   * (*Tabij)["Cgjl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIJi"]
                                   * (*Tabij)["Cfjl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIJi"]
                                   * (*Tabij)["Cejl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efgIJi"]
                                   * (*Tabij)["Chjk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efhIJi"]
                                   * (*Tabij)["Cgjk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ehgIJi"]
                                   * (*Tabij)["Cfjk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["hfgIJi"]
                                   * (*Tabij)["Cejk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJkl"]
                                   * (*Tabij)["ghKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJkl"]
                                   * (*Tabij)["fhKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJkl"]
                                   * (*Tabij)["ehKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJkl"]
                                   * (*Tabij)["gfKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJkl"]
                                   * (*Tabij)["geKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJkl"]
                                   * (*Tabij)["efKj"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJlj"]
                                   * (*Tabij)["ghKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJlj"]
                                   * (*Tabij)["fhKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJlj"]
                                   * (*Tabij)["ehKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJlj"]
                                   * (*Tabij)["gfKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJlj"]
                                   * (*Tabij)["geKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJlj"]
                                   * (*Tabij)["efKk"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJjk"]
                                   * (*Tabij)["ghKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJjk"]
                                   * (*Tabij)["fhKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJjk"]
                                   * (*Tabij)["ehKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJjk"]
                                   * (*Tabij)["gfKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJjk"]
                                   * (*Tabij)["geKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJjk"]
                                   * (*Tabij)["efKl"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJkl"]
                                   * (*Tabij)["ghKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJkl"]
                                   * (*Tabij)["fhKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJkl"]
                                   * (*Tabij)["ehKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJkl"]
                                   * (*Tabij)["gfKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJkl"]
                                   * (*Tabij)["geKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJkl"]
                                   * (*Tabij)["efKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJlj"]
                                   * (*Tabij)["ghKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJlj"]
                                   * (*Tabij)["fhKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJlj"]
                                   * (*Tabij)["ehKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJlj"]
                                   * (*Tabij)["gfKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJlj"]
                                   * (*Tabij)["geKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJlj"]
                                   * (*Tabij)["efKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJjk"]
                                   * (*Tabij)["ghKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJjk"]
                                   * (*Tabij)["fhKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJjk"]
                                   * (*Tabij)["ehKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJjk"]
                                   * (*Tabij)["gfKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJjk"]
                                   * (*Tabij)["geKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJjk"]
                                   * (*Tabij)["efKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJli"]
                                   * (*Tabij)["ghKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJli"]
                                   * (*Tabij)["fhKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJli"]
                                   * (*Tabij)["ehKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJli"]
                                   * (*Tabij)["gfKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJli"]
                                   * (*Tabij)["geKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJli"]
                                   * (*Tabij)["efKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJik"]
                                   * (*Tabij)["ghKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJik"]
                                   * (*Tabij)["fhKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJik"]
                                   * (*Tabij)["ehKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJik"]
                                   * (*Tabij)["gfKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJik"]
                                   * (*Tabij)["geKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJik"]
                                   * (*Tabij)["efKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJli"]
                                   * (*Tabij)["ghKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJli"]
                                   * (*Tabij)["fhKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJli"]
                                   * (*Tabij)["ehKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJli"]
                                   * (*Tabij)["gfKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJli"]
                                   * (*Tabij)["geKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJli"]
                                   * (*Tabij)["efKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJik"]
                                   * (*Tabij)["ghKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJik"]
                                   * (*Tabij)["fhKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJik"]
                                   * (*Tabij)["ehKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJik"]
                                   * (*Tabij)["gfKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJik"]
                                   * (*Tabij)["geKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AghJik"]
                                   * (*Tabij)["efKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJji"]
                                   * (*Tabij)["ghKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJji"]
                                   * (*Tabij)["fhKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJji"]
                                   * (*Tabij)["ehKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJji"]
                                   * (*Tabij)["gfKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJji"]
                                   * (*Tabij)["geKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJji"]
                                   * (*Tabij)["efKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJij"]
                                   * (*Tabij)["ghKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJij"]
                                   * (*Tabij)["fhKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJij"]
                                   * (*Tabij)["ehKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJij"]
                                   * (*Tabij)["gfKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJij"]
                                   * (*Tabij)["geKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJij"]
                                   * (*Tabij)["efKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIkl"]
                                   * (*Tabij)["BhKj"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIkl"]
                                   * (*Tabij)["BgKj"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIkl"]
                                   * (*Tabij)["BfKj"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIkl"]
                                   * (*Tabij)["BeKj"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIlj"]
                                   * (*Tabij)["BhKk"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIlj"]
                                   * (*Tabij)["BgKk"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIlj"]
                                   * (*Tabij)["BfKk"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIlj"]
                                   * (*Tabij)["BeKk"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIjk"]
                                   * (*Tabij)["BhKl"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIjk"]
                                   * (*Tabij)["BgKl"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIjk"]
                                   * (*Tabij)["BfKl"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIjk"]
                                   * (*Tabij)["BeKl"] * (*Tai)["Di"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIkl"]
                                   * (*Tabij)["BhKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIkl"]
                                   * (*Tabij)["BgKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIkl"]
                                   * (*Tabij)["BfKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIkl"]
                                   * (*Tabij)["BeKi"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIlj"]
                                   * (*Tabij)["BhKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIlj"]
                                   * (*Tabij)["BgKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIlj"]
                                   * (*Tabij)["BfKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIlj"]
                                   * (*Tabij)["BeKi"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIjk"]
                                   * (*Tabij)["BhKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIjk"]
                                   * (*Tabij)["BgKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIjk"]
                                   * (*Tabij)["BfKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIjk"]
                                   * (*Tabij)["BeKi"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIli"]
                                   * (*Tabij)["BhKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIli"]
                                   * (*Tabij)["BgKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIli"]
                                   * (*Tabij)["BfKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIli"]
                                   * (*Tabij)["BeKk"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIik"]
                                   * (*Tabij)["BhKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIik"]
                                   * (*Tabij)["BgKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIik"]
                                   * (*Tabij)["BfKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIik"]
                                   * (*Tabij)["BeKl"] * (*Tai)["Dj"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIli"]
                                   * (*Tabij)["BhKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIli"]
                                   * (*Tabij)["BgKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIli"]
                                   * (*Tabij)["BfKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIli"]
                                   * (*Tabij)["BeKj"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efgIik"]
                                   * (*Tabij)["BhKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efhIik"]
                                   * (*Tabij)["BgKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["ehgIik"]
                                   * (*Tabij)["BfKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["hfgIik"]
                                   * (*Tabij)["BeKj"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIji"]
                                   * (*Tabij)["BhKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIji"]
                                   * (*Tabij)["BgKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIji"]
                                   * (*Tabij)["BfKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIji"]
                                   * (*Tabij)["BeKl"] * (*Tai)["Dk"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIij"]
                                   * (*Tabij)["BhKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIij"]
                                   * (*Tabij)["BgKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIij"]
                                   * (*Tabij)["BfKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIij"]
                                   * (*Tabij)["BeKk"] * (*Tai)["Dl"]
                                   * (*Vijab)["IKBD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aefjkl"]
                                   * (*Tabij)["ghJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aegjkl"]
                                   * (*Tabij)["fhJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Agfjkl"]
                                   * (*Tabij)["ehJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aehjkl"]
                                   * (*Tabij)["gfJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Ahfjkl"]
                                   * (*Tabij)["geJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aghjkl"]
                                   * (*Tabij)["efJK"] * (*Tai)["Di"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aeflik"]
                                   * (*Tabij)["ghJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aeglik"]
                                   * (*Tabij)["fhJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agflik"]
                                   * (*Tabij)["ehJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehlik"]
                                   * (*Tabij)["gfJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahflik"]
                                   * (*Tabij)["geJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghlik"]
                                   * (*Tabij)["efJK"] * (*Tai)["Dj"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aefilj"]
                                   * (*Tabij)["ghJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegilj"]
                                   * (*Tabij)["fhJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agfilj"]
                                   * (*Tabij)["ehJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehilj"]
                                   * (*Tabij)["gfJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfilj"]
                                   * (*Tabij)["geJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghilj"]
                                   * (*Tabij)["efJK"] * (*Tai)["Dk"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aefijk"]
                                   * (*Tabij)["ghJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aegijk"]
                                   * (*Tabij)["fhJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Agfijk"]
                                   * (*Tabij)["ehJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Aehijk"]
                                   * (*Tabij)["gfJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["Ahfijk"]
                                   * (*Tabij)["geJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["Aghijk"]
                                   * (*Tabij)["efJK"] * (*Tai)["Dl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJkl"]
                                   * (*Tabij)["Cgij"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJkl"]
                                   * (*Tabij)["Cfij"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJkl"]
                                   * (*Tabij)["Ceij"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJkl"]
                                   * (*Tabij)["Chij"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJkl"]
                                   * (*Tabij)["Chij"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJkl"]
                                   * (*Tabij)["Chij"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJkl"]
                                   * (*Tabij)["Cfij"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJkl"]
                                   * (*Tabij)["Ceij"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJkl"]
                                   * (*Tabij)["Cgij"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJkl"]
                                   * (*Tabij)["Cgij"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJkl"]
                                   * (*Tabij)["Ceij"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJkl"]
                                   * (*Tabij)["Cfij"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJlj"]
                                   * (*Tabij)["Cgik"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJlj"]
                                   * (*Tabij)["Cfik"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJlj"]
                                   * (*Tabij)["Ceik"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJlj"]
                                   * (*Tabij)["Chik"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJlj"]
                                   * (*Tabij)["Chik"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJlj"]
                                   * (*Tabij)["Chik"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJlj"]
                                   * (*Tabij)["Cfik"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJlj"]
                                   * (*Tabij)["Ceik"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJlj"]
                                   * (*Tabij)["Cgik"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJlj"]
                                   * (*Tabij)["Cgik"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJlj"]
                                   * (*Tabij)["Ceik"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJlj"]
                                   * (*Tabij)["Cfik"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJjk"]
                                   * (*Tabij)["Cgil"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJjk"]
                                   * (*Tabij)["Cfil"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJjk"]
                                   * (*Tabij)["Ceil"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJjk"]
                                   * (*Tabij)["Chil"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJjk"]
                                   * (*Tabij)["Chil"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJjk"]
                                   * (*Tabij)["Chil"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJjk"]
                                   * (*Tabij)["Cfil"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJjk"]
                                   * (*Tabij)["Ceil"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJjk"]
                                   * (*Tabij)["Cgil"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJjk"]
                                   * (*Tabij)["Cgil"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJjk"]
                                   * (*Tabij)["Ceil"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJjk"]
                                   * (*Tabij)["Cfil"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJli"]
                                   * (*Tabij)["Cgkj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJli"]
                                   * (*Tabij)["Cfkj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJli"]
                                   * (*Tabij)["Cekj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJli"]
                                   * (*Tabij)["Chkj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJli"]
                                   * (*Tabij)["Chkj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJli"]
                                   * (*Tabij)["Chkj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJli"]
                                   * (*Tabij)["Cfkj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJli"]
                                   * (*Tabij)["Cekj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJli"]
                                   * (*Tabij)["Cgkj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJli"]
                                   * (*Tabij)["Cgkj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJli"]
                                   * (*Tabij)["Cekj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJli"]
                                   * (*Tabij)["Cfkj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJik"]
                                   * (*Tabij)["Cglj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJik"]
                                   * (*Tabij)["Cflj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJik"]
                                   * (*Tabij)["Celj"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJik"]
                                   * (*Tabij)["Chlj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJik"]
                                   * (*Tabij)["Chlj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJik"]
                                   * (*Tabij)["Chlj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJik"]
                                   * (*Tabij)["Cflj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJik"]
                                   * (*Tabij)["Celj"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJik"]
                                   * (*Tabij)["Cglj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJik"]
                                   * (*Tabij)["Cglj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJik"]
                                   * (*Tabij)["Celj"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJik"]
                                   * (*Tabij)["Cflj"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AefJij"]
                                   * (*Tabij)["Cgkl"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AegJij"]
                                   * (*Tabij)["Cfkl"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AgfJij"]
                                   * (*Tabij)["Cekl"] * (*Tai)["hL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AefJij"]
                                   * (*Tabij)["Chkl"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AegJij"]
                                   * (*Tabij)["Chkl"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AgfJij"]
                                   * (*Tabij)["Chkl"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AehJij"]
                                   * (*Tabij)["Cfkl"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhfJij"]
                                   * (*Tabij)["Cekl"] * (*Tai)["gL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AehJij"]
                                   * (*Tabij)["Cgkl"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["AhfJij"]
                                   * (*Tabij)["Cgkl"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AghJij"]
                                   * (*Tabij)["Cekl"] * (*Tai)["fL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["AhgJij"]
                                   * (*Tabij)["Cfkl"] * (*Tai)["eL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIkl"]
                                   * (*Tabij)["BCij"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIkl"]
                                   * (*Tabij)["BCij"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIkl"]
                                   * (*Tabij)["BCij"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIkl"]
                                   * (*Tabij)["BCij"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIlj"]
                                   * (*Tabij)["BCik"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIlj"]
                                   * (*Tabij)["BCik"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIlj"]
                                   * (*Tabij)["BCik"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIlj"]
                                   * (*Tabij)["BCik"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIjk"]
                                   * (*Tabij)["BCil"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIjk"]
                                   * (*Tabij)["BCil"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIjk"]
                                   * (*Tabij)["BCil"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIjk"]
                                   * (*Tabij)["BCil"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIli"]
                                   * (*Tabij)["BCkj"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIli"]
                                   * (*Tabij)["BCkj"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIli"]
                                   * (*Tabij)["BCkj"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIli"]
                                   * (*Tabij)["BCkj"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIik"]
                                   * (*Tabij)["BClj"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIik"]
                                   * (*Tabij)["BClj"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIik"]
                                   * (*Tabij)["BClj"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIik"]
                                   * (*Tabij)["BClj"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["efgIij"]
                                   * (*Tabij)["BCkl"] * (*Tai)["hL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["efhIij"]
                                   * (*Tabij)["BCkl"] * (*Tai)["gL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ehgIij"]
                                   * (*Tabij)["BCkl"] * (*Tai)["fL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["hfgIij"]
                                   * (*Tabij)["BCkl"] * (*Tai)["eL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIkl"]
                                   * (*Tabij)["Bhij"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIkl"]
                                   * (*Tabij)["Bgij"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIkl"]
                                   * (*Tabij)["Bfij"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIkl"]
                                   * (*Tabij)["Beij"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIlj"]
                                   * (*Tabij)["Bhik"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIlj"]
                                   * (*Tabij)["Bgik"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIlj"]
                                   * (*Tabij)["Bfik"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIlj"]
                                   * (*Tabij)["Beik"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIjk"]
                                   * (*Tabij)["Bhil"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIjk"]
                                   * (*Tabij)["Bgil"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIjk"]
                                   * (*Tabij)["Bfil"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIjk"]
                                   * (*Tabij)["Beil"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIli"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIli"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIli"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIli"]
                                   * (*Tabij)["Bekj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIik"]
                                   * (*Tabij)["Bhlj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIik"]
                                   * (*Tabij)["Bglj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIik"]
                                   * (*Tabij)["Bflj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIik"]
                                   * (*Tabij)["Belj"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["efgIij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["efhIij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["ehgIij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["hfgIij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["CL"]
                                   * (*Vijab)["ILBC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABejkl"]
                                   * (*Tabij)["fgKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfjkl"]
                                   * (*Tabij)["egKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgjkl"]
                                   * (*Tabij)["feKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABejkl"]
                                   * (*Tabij)["fhKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfjkl"]
                                   * (*Tabij)["ehKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABejkl"]
                                   * (*Tabij)["hgKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfjkl"]
                                   * (*Tabij)["hgKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgjkl"]
                                   * (*Tabij)["heKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgjkl"]
                                   * (*Tabij)["fhKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"]
                                   * (*Tabij)["efKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"]
                                   * (*Tabij)["geKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABhjkl"]
                                   * (*Tabij)["fgKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABelik"]
                                   * (*Tabij)["fgKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABflik"]
                                   * (*Tabij)["egKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABglik"]
                                   * (*Tabij)["feKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABelik"]
                                   * (*Tabij)["fhKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABflik"]
                                   * (*Tabij)["ehKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABelik"]
                                   * (*Tabij)["hgKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABflik"]
                                   * (*Tabij)["hgKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABglik"]
                                   * (*Tabij)["heKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABglik"]
                                   * (*Tabij)["fhKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"]
                                   * (*Tabij)["efKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"]
                                   * (*Tabij)["geKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhlik"]
                                   * (*Tabij)["fgKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABeilj"]
                                   * (*Tabij)["fgKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfilj"]
                                   * (*Tabij)["egKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgilj"]
                                   * (*Tabij)["feKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeilj"]
                                   * (*Tabij)["fhKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfilj"]
                                   * (*Tabij)["ehKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeilj"]
                                   * (*Tabij)["hgKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfilj"]
                                   * (*Tabij)["hgKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgilj"]
                                   * (*Tabij)["heKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgilj"]
                                   * (*Tabij)["fhKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"]
                                   * (*Tabij)["efKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"]
                                   * (*Tabij)["geKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhilj"]
                                   * (*Tabij)["fgKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABeijk"]
                                   * (*Tabij)["fgKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABfijk"]
                                   * (*Tabij)["egKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABgijk"]
                                   * (*Tabij)["feKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeijk"]
                                   * (*Tabij)["fhKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfijk"]
                                   * (*Tabij)["ehKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABeijk"]
                                   * (*Tabij)["hgKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABfijk"]
                                   * (*Tabij)["hgKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgijk"]
                                   * (*Tabij)["heKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabcijk)["ABgijk"]
                                   * (*Tabij)["fhKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"]
                                   * (*Tabij)["efKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"]
                                   * (*Tabij)["geKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabcijk)["ABhijk"]
                                   * (*Tabij)["fgKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefjkl"]
                                   * (*Tabij)["BgKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegjkl"]
                                   * (*Tabij)["BfKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfjkl"]
                                   * (*Tabij)["BeKi"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefjkl"]
                                   * (*Tabij)["BhKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegjkl"]
                                   * (*Tabij)["BhKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfjkl"]
                                   * (*Tabij)["BhKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehjkl"]
                                   * (*Tabij)["BfKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfjkl"]
                                   * (*Tabij)["BeKi"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehjkl"]
                                   * (*Tabij)["BgKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfjkl"]
                                   * (*Tabij)["BgKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aghjkl"]
                                   * (*Tabij)["BeKi"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahgjkl"]
                                   * (*Tabij)["BfKi"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeflik"]
                                   * (*Tabij)["BgKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeglik"]
                                   * (*Tabij)["BfKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agflik"]
                                   * (*Tabij)["BeKj"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeflik"]
                                   * (*Tabij)["BhKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeglik"]
                                   * (*Tabij)["BhKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agflik"]
                                   * (*Tabij)["BhKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehlik"]
                                   * (*Tabij)["BfKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahflik"]
                                   * (*Tabij)["BeKj"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehlik"]
                                   * (*Tabij)["BgKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahflik"]
                                   * (*Tabij)["BgKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghlik"]
                                   * (*Tabij)["BeKj"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahglik"]
                                   * (*Tabij)["BfKj"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefilj"]
                                   * (*Tabij)["BgKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegilj"]
                                   * (*Tabij)["BfKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfilj"]
                                   * (*Tabij)["BeKk"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefilj"]
                                   * (*Tabij)["BhKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegilj"]
                                   * (*Tabij)["BhKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfilj"]
                                   * (*Tabij)["BhKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehilj"]
                                   * (*Tabij)["BfKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfilj"]
                                   * (*Tabij)["BeKk"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehilj"]
                                   * (*Tabij)["BgKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfilj"]
                                   * (*Tabij)["BgKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghilj"]
                                   * (*Tabij)["BeKk"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahgilj"]
                                   * (*Tabij)["BfKk"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefijk"]
                                   * (*Tabij)["BgKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegijk"]
                                   * (*Tabij)["BfKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfijk"]
                                   * (*Tabij)["BeKl"] * (*Tai)["hL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefijk"]
                                   * (*Tabij)["BhKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegijk"]
                                   * (*Tabij)["BhKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfijk"]
                                   * (*Tabij)["BhKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehijk"]
                                   * (*Tabij)["BfKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfijk"]
                                   * (*Tabij)["BeKl"] * (*Tai)["gL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehijk"]
                                   * (*Tabij)["BgKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfijk"]
                                   * (*Tabij)["BgKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghijk"]
                                   * (*Tabij)["BeKl"] * (*Tai)["fL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahgijk"]
                                   * (*Tabij)["BfKl"] * (*Tai)["eL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aefjkl"]
                                   * (*Tabij)["ghJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aegjkl"]
                                   * (*Tabij)["fhJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Agfjkl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aehjkl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Ahfjkl"]
                                   * (*Tabij)["geJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aghjkl"]
                                   * (*Tabij)["efJi"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aeflik"]
                                   * (*Tabij)["ghJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aeglik"]
                                   * (*Tabij)["fhJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agflik"]
                                   * (*Tabij)["ehJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehlik"]
                                   * (*Tabij)["gfJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahflik"]
                                   * (*Tabij)["geJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghlik"]
                                   * (*Tabij)["efJj"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefilj"]
                                   * (*Tabij)["ghJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegilj"]
                                   * (*Tabij)["fhJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfilj"]
                                   * (*Tabij)["ehJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehilj"]
                                   * (*Tabij)["gfJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfilj"]
                                   * (*Tabij)["geJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghilj"]
                                   * (*Tabij)["efJk"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aefijk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aegijk"]
                                   * (*Tabij)["fhJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Agfijk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Aehijk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabcijk)["Ahfijk"]
                                   * (*Tabij)["geJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabcijk)["Aghijk"]
                                   * (*Tabij)["efJl"] * (*Tai)["CL"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bgkl"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bfkl"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bekl"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bhkl"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bhkl"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bhkl"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bfkl"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bekl"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bgkl"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bgkl"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bekl"] * (*Tabij)["ghKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bfkl"] * (*Tabij)["hgKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bgjl"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bfjl"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bejl"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bhjl"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bhjl"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bhjl"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bfjl"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bejl"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bgjl"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bgjl"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bejl"] * (*Tabij)["ghKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bfjl"] * (*Tabij)["hgKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bgkj"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bfkj"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bekj"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bhkj"] * (*Tabij)["efKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bhkj"] * (*Tabij)["egKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bhkj"] * (*Tabij)["gfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bfkj"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bekj"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bgkj"] * (*Tabij)["ehKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bgkj"] * (*Tabij)["hfKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bekj"] * (*Tabij)["ghKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bfkj"] * (*Tabij)["hgKL"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["geJk"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJk"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["geJl"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["heJk"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJl"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJk"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["geJl"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJk"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["heJl"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["geJk"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["geJj"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["efJj"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["geJl"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["heJj"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["efJl"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["ehJj"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["geJl"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["efJj"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["heJl"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["geJj"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["heJi"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["efJl"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["ehJi"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["geJl"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["CeLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["geJi"] * (*Tabij)["CfLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["efJi"] * (*Tabij)["CgLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["geJl"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkj"]
                                   * (*Tabij)["fgJl"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["geJi"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["efJi"] * (*Tabij)["ChLl"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["heJl"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkj"]
                                   * (*Tabij)["fhJl"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekj"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkj"]
                                   * (*Tabij)["hgJl"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["geJj"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["efJk"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["efJj"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["geJk"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["heJj"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["efJk"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["ehJj"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["geJk"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["efJj"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["heJk"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["geJj"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["ehJk"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["heJi"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["efJk"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["ehJi"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["geJk"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["CeLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["geJi"] * (*Tabij)["CfLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["efJk"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["efJi"] * (*Tabij)["CgLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["geJk"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlj"]
                                   * (*Tabij)["fgJk"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["geJi"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["ehJk"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["efJi"] * (*Tabij)["ChLk"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["heJk"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglj"]
                                   * (*Tabij)["fhJk"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelj"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflj"]
                                   * (*Tabij)["hgJk"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekl"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"]
                                   * (*Tabij)["hgJi"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aekl"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"]
                                   * (*Tabij)["ehJi"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["geJj"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"]
                                   * (*Tabij)["fhJi"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"]
                                   * (*Tabij)["heJi"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJj"] * (*Tabij)["ChLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aekl"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afkl"]
                                   * (*Tabij)["geJi"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJj"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJj"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agkl"]
                                   * (*Tabij)["efJi"] * (*Tabij)["ChLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["heJj"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["fhJj"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"]
                                   * (*Tabij)["fgJi"] * (*Tabij)["CeLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"]
                                   * (*Tabij)["geJi"] * (*Tabij)["CfLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJj"] * (*Tabij)["CgLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahkl"]
                                   * (*Tabij)["efJi"] * (*Tabij)["CgLj"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["geJj"] * (*Tabij)["CfLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["fgJj"] * (*Tabij)["CeLi"]
                                   * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIk"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIk"]
                                   * (*Tabij)["egJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIk"]
                                   * (*Tabij)["gfJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIk"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIk"]
                                   * (*Tabij)["hfJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIk"]
                                   * (*Tabij)["ghJl"] * (*Tabij)["CDij"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJl"] * (*Tabij)["CDik"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+0.5) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJl"] * (*Tabij)["CDkj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJk"] * (*Tabij)["CDil"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJk"] * (*Tabij)["CDlj"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-0.5) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJj"] * (*Tabij)["CDkl"]
                                   * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bj"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLkl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLjl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLjl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLjl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bk"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLjl"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLkj"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLkj"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLkj"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["Bl"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLkj"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bk"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLil"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["Bl"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLki"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tai)["hK"] * (*Tabcijk)["efgLij"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tai)["gK"] * (*Tabcijk)["efhLij"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tai)["fK"] * (*Tabcijk)["ehgLij"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["Bl"]
                                   * (*Tai)["eK"] * (*Tabcijk)["hfgLij"]
                                   * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["gJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Defjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Degjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Dgfjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["fJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dehjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dhfjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ai"] * (*Tai)["eJ"]
                                   * (*Tai)["fK"] * (*Tabcijk)["Dghjkl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["gJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Defikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Degikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Dgfikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["fJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dehikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dhfikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Aj"] * (*Tai)["eJ"]
                                   * (*Tai)["fK"] * (*Tabcijk)["Dghikl"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["gJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Defjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Degjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Dgfjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["fJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dehjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dhfjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Ak"] * (*Tai)["eJ"]
                                   * (*Tai)["fK"] * (*Tabcijk)["Dghjil"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["gJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Defjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Degjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tai)["hK"] * (*Tabcijk)["Dgfjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["fJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dehjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tai)["gK"] * (*Tabcijk)["Dhfjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tai)["Al"] * (*Tai)["eJ"]
                                   * (*Tai)["fK"] * (*Tabcijk)["Dghjki"]
                                   * (*Vijab)["JKAD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dj"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIj"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIj"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIj"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIj"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIj"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIj"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Ci"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Cj"]
                                   * (*Tai)["Dk"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cj"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["ghIi"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["hfIi"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["heIi"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["fgIi"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["egIi"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["efIi"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ck"]
                                   * (*Tai)["Dl"] * (*Vijab)["IJCD"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJk"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJk"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ci"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ci"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ci"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ci"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJl"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJl"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahli"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agli"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeli"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afli"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cj"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aelk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aflk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aglk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cj"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cj"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahlk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cj"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJl"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJl"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJl"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJl"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJl"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJl"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJl"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["efJk"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["egJk"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["gfJk"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["ehJk"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["hfJk"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["ghJk"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["hgJk"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["efJj"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["egJj"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["Ck"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejl"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjl"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjl"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["efJi"] * (*Tai)["Ck"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["egJi"] * (*Tai)["Ck"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjl"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Ck"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["efJj"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["egJj"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["efJj"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["egJj"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahki"]
                                   * (*Tabij)["gfJj"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["ehJj"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agki"]
                                   * (*Tabij)["hfJj"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeki"]
                                   * (*Tabij)["ghJj"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afki"]
                                   * (*Tabij)["hgJj"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cl"]
                                   * (*Tai)["hL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aejk"]
                                   * (*Tabij)["ghJi"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afjk"]
                                   * (*Tabij)["hgJi"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["ehJi"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agjk"]
                                   * (*Tabij)["hfJi"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["efJi"] * (*Tai)["Cl"]
                                   * (*Tai)["gL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["egJi"] * (*Tai)["Cl"]
                                   * (*Tai)["fL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahjk"]
                                   * (*Tabij)["gfJi"] * (*Tai)["Cl"]
                                   * (*Tai)["eL"] * (*Vijab)["JLAC"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bekl"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bfkl"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahij"]
                                   * (*Tabij)["Bgkl"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agij"]
                                   * (*Tabij)["Bhkl"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bejl"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bfjl"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahik"]
                                   * (*Tabij)["Bgjl"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agik"]
                                   * (*Tabij)["Bhjl"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["gK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["fK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["eK"]
                                   * (*Tai)["hL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bekj"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bfkj"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Ahil"]
                                   * (*Tabij)["Bgkj"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Aeil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["fK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (-1.0) * (*Tabij)["Afil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["eK"]
                                   * (*Tai)["gL"] * (*Vijab)["KLAB"]);
  LDEBUG((*Wabcdijkl)["efghijkl"] += (+1.0) * (*Tabij)["Agil"]
                                   * (*Tabij)["Bhkj"] * (*Tai)["eK"]
                                   * (*Tai)["fL"] * (*Vijab)["KLAB"]);

  return Wabcdijkl;
}

// instantiate
template class sisi4s::SimilarityTransformedHamiltonian<sisi4s::complex>;
template class sisi4s::SimilarityTransformedHamiltonian<double>;
