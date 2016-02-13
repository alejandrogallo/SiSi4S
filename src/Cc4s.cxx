/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Cc4s.hpp>
#include <Algorithm.hpp>
#include <Parser.hpp>
#include <util/Timer.hpp>
#include <util/FlopsCounter.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <fstream>

// TODO: to be removed from the main class
#include <util/MathFunctions.hpp>

using namespace cc4s;
using namespace CTF;

Cc4s::Cc4s() {
}

Cc4s::~Cc4s() {
}

void Cc4s::run() {
  printBanner();
  Parser parser(options->file);
  std::vector<Algorithm *> algorithms(parser.parse());
  LOG(0, "root") <<
    "execution plan read, steps=" << algorithms.size() << std::endl;

  int64_t rootFlops, totalFlops;
  double totalTime;
  {
    FlopsCounter rootCounter(&rootFlops);
    FlopsCounter totalCounter(&totalFlops, world->comm);
    Timer totalTimer(&totalTime);

    for (unsigned int i(0); i < algorithms.size(); ++i) {
      LOG(0, "root") << "step=" << (i+1) << ", " << algorithms[i]->getName() << std::endl;

      int64_t flops;
      double time;
      {
        FlopsCounter flopsCounter(&flops);
        Timer timer(&time);
        algorithms[i]->run();
      }

      LOG(1, "root") << "step=" << (i+1) << ", realtime=" << time << " s"
        << ", operations=" << flops / 1e9 << " GFLOPS/core"
        << ", speed=" << flops / 1e9 / time << " GFLOPS/s/core" << std::endl;
    }
  }

  printStatistics(rootFlops, totalFlops, totalTime);
}


void Cc4s::printBanner() {
  OUT() << "                __ __      " << std::endl
        << "     __________/ // / _____" << std::endl
        << "    / ___/ ___/ // /_/ ___/" << std::endl
        << "   / /__/ /__/__  __(__  ) " << std::endl
        << "   \\___/\\___/  /_/ /____/  " << std::endl
        << "  Coupled Cluster for Solids" << std::endl << std::endl;
  LOG(0, "root") << "version=" << CC4S_VERSION <<
    ", date=" << CC4S_DATE << std::endl;
  LOG(0, "root") << "build date=" << __DATE__ << " " << __TIME__ << std::endl;
  LOG(0, "root") << "compiler=" << COMPILER_VERSION << std::endl << std::endl;
}

void Cc4s::printStatistics(
  int64_t rootFlops, int64_t totalFlops, double totalTime
) {
  std::string pid, comm, state, ppid, pgrp, session, ttyNr,
    tpgid, flags, minflt, cminflt, majflt, cmajflt,
    utime, stime, cutime, cstime, priority, nice,
    O, itrealvalue, starttime;
  int64_t vsize, rss;
  // assuming LINUX 
  std::ifstream statStream("/proc/self/stat", std::ios_base::in);
  statStream >> pid >> comm >> state >> ppid >> pgrp >> session >> ttyNr
    >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
    >> utime >> stime >> cutime >> cstime >> priority >> nice
    >> O >> itrealvalue >> starttime >> vsize >> rss;
  statStream.close();
  // in case x86-64 is configured to use 2MB pages
  int64_t pageSize = sysconf(_SC_PAGE_SIZE);
  LOG(0, "root") << "total realtime=" << totalTime << " s" << std::endl;
  LOG(0, "root") << "total operations=" << rootFlops / 1e9 << " GFLOPS/core"
    << " speed=" << rootFlops/1e9 / totalTime << " GFLOPS/s/core" << std::endl;
  LOG(0, "root") << "physical memory=" << rss * pageSize / 1e9 << " GB/core"
    << ", virtual memory: " << vsize / 1e9 << " GB/core" << std::endl;

  int64_t globalVSize, globalRss;
  MPI_Reduce(&vsize, &globalVSize, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  MPI_Reduce(&rss, &globalRss, 1, MPI_LONG_LONG, MPI_SUM, 0, world->comm);
  LOG(0, "root") << "overall operations=" << totalFlops / 1.e9 << " GFLOPS"
    << std::endl;
  LOG(0, "root") << "overall physical memory="
    << globalRss * pageSize / 1e9 << " GB"
    << ", overall virtual memory=" << globalVSize / 1e9 << " GB" << std::endl;
}


void Cc4s::iterateRpa() {
  {
    int syms[] = {NS, NS, NS, NS};
    // Define Tensors
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    //Allocate Tensors for RPA
    Tensor<> Rabij = Tensor<>(V->abij);
    Tensor<> Cabij = Tensor<>(V->abij);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);

    LOG(0) << "Solving RPA Amplitude Equations:" << std::endl;


    Rabij["abij"] = (*V)["abij"];
    Rabij["abij"] += 2.0 * (*V)["acik"] * (*T)["cbkj"];
    Cabij["abij"] =  2.0 * (*V)["cbkj"] * (*T)["acik"];
    Rabij["abij"] += Cabij["abij"];
    //(*V)["cbkj"]*(*T)["acjk"];
    Rabij["abij"] += 2.0 * Cabij["acik"] * (*T)["cbkj"];


    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fDivide(&divide<double>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
}

void Cc4s::iterateRccd() {
  {
    int syms[] = {NS, NS, NS, NS};
    // Define Tensors
    //Allocate Tensors for T1 amplitude equations
    Tensor<> Dai(V->ai);
    Tensor<> Rai(T->ai);
    //Allocate Tensors for T2 amplitudes
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Tensor<> Rabij(V->abij);
    Tensor<> Fba(V->ab);
    Tensor<> Fji(V->ij);
    Tensor<> Fai(V->ai);
    //intermediates
    Tensor<> Cai(T->ai);
    Tensor<> Lac(V->ab);
    Tensor<> Kac(V->ab);
    Tensor<> Lki(V->ij);
    Tensor<> Kki(V->ij);
    Tensor<> Cklij(V->ijkl);
    Tensor<> Cakic(V->aijb);
    Tensor<> Cakci(V->aibj);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);


    LOG(0) << "Solving restricted T2 CCD Amplitude Equations:" << std::endl;



//Build Kac
//    Kac["ac"] = Fba["ac"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["adkl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["adkl"];

//Build Lac
    Lac["ac"] = Kac["ac"];

//Build Kki
//    Kki["ki"] = Fji["ki"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["cdil"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["cdil"];

//Build Lki
    Lki["ki"] = Kki["ki"];

//  Contract L_ac with T2 Amplitudes
    Rabij["abij"] = Lac["ac"] * (*T)["cbij"];

//  Contract L_ki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*T)["abkj"];

//Build C_akic
    Cakic["akic"] = (*V)["akic"];
    Cakic["akic"] -= 0.5 * (*V)["lkdc"] * (*T)["dail"];
    Cakic["akic"] += (*V)["lkdc"] * (*T)["adil"];
    Cakic["akic"] -= 0.5 * (*V)["lkcd"] * (*T)["adil"];

//Build C_akci
    Cakci["akci"] = (*V)["akci"];
    Cakci["akci"] -= 0.5 * (*V)["lkcd"] * (*T)["dail"];

//  Contract C_akic and C_akci intermediates with T2 amplitudes

    Rabij["abij"] += 2.0 * Cakic["akic"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakic["akic"] * (*T)["bckj"];

    Rabij["abij"] -= Cakci["akci"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakci["bkci"] * (*T)["ackj"];

// Symmetrize Rabij by applying permutation operator
    // to save memory we use Cakci as intermediate for the permutation operator 
    Cakci["aibj"] = Rabij["abij"];
    Rabij["abij"] += Cakci["bjai"]; 

//////////////////////////////////////////////////////////////////////
// Now add all terms to Rabij that do not need to be symmetrized with
// the permutation operator
//////////////////////////////////////////////////////////////////////

//  Rabij are the Tabij amplitudes for the next iteration and need to be build
    Rabij["abij"] += (*V)["abij"];

//  Build Chi_klij intermediate
    Cklij["klij"] = (*V)["klij"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["cdij"];

//  Contract Chi_klij with T2 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["abkl"];

    if (V->abcd) {
      Tensor<> Cabcd(V->abcd);
  //  Build Chi_abcd intermediate
      Cabcd["abcd"] = (*V)["abcd"];
  //  Contract Chi_abcd with T2 Amplitudes
      Rabij["abij"] += Cabcd["abcd"] * (*T)["cdij"];
    } else {
  // Slicing:
      for (int b(0); b < chiReal->nv; b += options->nw) {
        for (int a(b); a < chiReal->nv; a += options->nw) {
          LOG(0) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          Tensor<> Vxycd(V->getSlice(a, b));
          int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);
          int origin[] = {0, 0, 0, 0};
          int lens[] = {na, nb, chiReal->no, chiReal->no};
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
          Rxyij["xyij"] = Vxycd["xycd"] * (*T)["cdij"];

          int rBegin[] = {a, b, 0, 0};
          int rEnd[] = {a+na, b+nb, chiReal->no, chiReal->no};
          // R["abij"] += R["xyij"] at current x,y
          Rabij.slice(rBegin,rEnd,1.0, Rxyij,origin,lens,1.0);
          if (a>b) {
            // add the same slice at (b,a,j,i):
            rBegin[0] = b; rBegin[1] = a;
            rEnd[0] = b+nb; rEnd[1] = a+na;
            // note that na may be != nb
            lens[0] = nb; lens[1] = na;
            Tensor<> Ryxji(4,lens,syms,*world, "Ryxij", Vxycd.profile);
            Ryxji["yxji"] = Rxyij["xyij"];
            Rabij.slice(rBegin,rEnd,1.0, Ryxji,origin,lens,1.0);
          }
        }
      }
    }

    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fDivide(&divide<double>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
}


void Cc4s::iterateRccsd() {
  {
    int syms[] = {NS, NS, NS, NS};
    LOG(0) << "Allocating stuff:" << std::endl;
    // Define Tensors
    //Allocate Tensors for T1 amplitude equations
    Tensor<> Rai(T->ai);
    //Allocate Tensors for T2 amplitudes
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Tensor<> Rabij(V->abij);
    Tensor<> Fba(V->ab);
    Tensor<> Fji(V->ij);
    Tensor<> Fai(V->ai);
    LOG(0) << "Allocating stuff 2:" << std::endl;
    //intermediates
    LOG(0) << "Allocating stuff 2a:" << std::endl;
    Tensor<> Dai(T->ai);
    Tensor<> Lac(V->ab);
    Tensor<> Kac(V->ab);
    Tensor<> Lki(V->ij);
    Tensor<> Kki(V->ij);
    Tensor<> Kck(T->ai);
    LOG(0) << "Allocating stuff 2b:" << std::endl;
    Tensor<> Cklij(V->ijkl);
    Tensor<> Cabcd(V->abcd);
    Tensor<> Cakic(V->aijb);
    Tensor<> Cakci(V->aibj);
    //Tensor<> Chi(4, V->abij->lens, syms, *world, "Cabij");
    //Chi = new Amplitudes(V, world, options);
    LOG(0) << "Allocating stuff 3:" << std::endl;


//********************************************************************************
//  T2 amplitude equations
//********************************************************************************

    LOG(0) << "Solving restricted T2 CCSD Amplitude Equations:" << std::endl;

//Build Kac
    Kac["ac"] = Fba["ac"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["adkl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["adkl"];
    Kac["ac"] -= 2.0 * (*V)["klcd"] * (*T)["ak"] * (*T)["dl"];
    Kac["ac"] += (*V)["kldc"] * (*T)["ak"] * (*T)["dl"];

    LOG(0) << "Building Lac:" << std::endl;
//Build Lac
    Lac["ac"] = Kac["ac"];
    Lac["ac"] -= Fai["ck"] * (*T)["ak"];
    Lac["ac"] += 2.0 * (*V)["akcd"] * (*T)["dk"];
    Lac["ac"] -= (*V)["akdc"] * (*T)["dk"];


    LOG(0) << "Building Kki:" << std::endl;
//Build Kki
    Kki["ki"] = Fji["ki"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["cdil"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["cdil"];
    Kki["ki"] += 2.0 * (*V)["klcd"] * (*T)["ci"] * (*T)["dl"];
    Kki["ki"] -= (*V)["kldc"] * (*T)["ci"] * (*T)["dl"];

    LOG(0) << "Building Lki:" << std::endl;
//Build Lki
    Lki["ki"] = Kki["ki"];
    Lki["ki"] += Fai["ck"] * (*T)["ci"];
    Lki["ki"] += 2.0 * (*V)["lkci"] * (*T)["cl"];
    Lki["ki"] -= (*V)["klci"] * (*T)["cl"];

    LOG(0) << "Building Rabij:" << std::endl;
//  Contract L_ac with T2 Amplitudes
    Rabij["abij"] = Lac["ac"] * (*T)["cbij"];

//  Contract L_ki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*T)["abkj"];

//  Contract Coulomb integrals with T2 amplitudes

    Rabij["abij"] += (*V)["baci"] * (*T)["cj"];

    Rabij["abij"] -= (*V)["bkci"] * (*T)["ak"] * (*T)["cj"];

    Rabij["abij"] -= (*V)["akij"] * (*T)["bk"];

    Rabij["abij"] += (*V)["akic"] * (*T)["cj"] * (*T)["bk"];

    LOG(0) << "Building Cakic:" << std::endl;
//Build C_akic
    Cakic["akic"] = (*V)["akic"];
    Cakic["akic"] -= (*V)["klci"] * (*T)["al"];
    Cakic["akic"] += (*V)["akdc"] * (*T)["di"];
    Cakic["akic"] -= 0.5 * (*V)["lkdc"] * (*T)["dail"];
    Cakic["akic"] -= (*V)["lkdc"] * (*T)["di"] * (*T)["al"];
    Cakic["akic"] += (*V)["lkdc"] * (*T)["adil"];
    Cakic["akic"] -= 0.5 * (*V)["lkcd"] * (*T)["adil"];

//Build C_akci
    Cakci["akci"] = (*V)["akci"];
    Cakci["akci"] -= (*V)["lkci"] * (*T)["al"];
    Cakci["akci"] += (*V)["akcd"] * (*T)["di"];
    Cakci["akci"] -= 0.5 * (*V)["lkcd"] * (*T)["dail"];
    Cakci["akci"] -= (*V)["lkcd"] * (*T)["di"] * (*T)["al"];

//  Contract C_akic and C_akci intermediates with T2 amplitudes

    Rabij["abij"] += 2.0 * Cakic["akic"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakic["akic"] * (*T)["bckj"];

    Rabij["abij"] -= Cakci["akci"] * (*T)["cbkj"];
    Rabij["abij"] -= Cakci["bkci"] * (*T)["ackj"];

// Symmetrize Rabij by applying permutation operator
    // to save memory we use Cakci as intermediate for the permutation operator 
    Cakci["aibj"] = Rabij["abij"];
    Rabij["abij"] += Cakci["bjai"]; 



//////////////////////////////////////////////////////////////////////
// Now add all terms to Rabij that do not need to be symmetrized with
// the permutation operator
//////////////////////////////////////////////////////////////////////

//  Rabij are the Tabij amplitudes for the next iteration and need to be build
    Rabij["abij"] += (*V)["abij"];

//  Build Chi_klij intermediate
    Cklij["klij"] = (*V)["klij"];
    Cklij["klij"] += (*V)["lkci"] * (*T)["cj"];
    Cklij["klij"] += (*V)["klcj"] * (*T)["ci"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["cdij"];
    Cklij["klij"] += (*V)["klcd"] * (*T)["ci"] * (*T)["dj"]; 

//  Contract Chi_klij with T2 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["abkl"];

//  Contract Chi_klij with T1 Amplitudes
    Rabij["abij"] += Cklij["klij"] * (*T)["ak"] * (*T)["bl"];

//*****************************
//@Felix start slicing here
//*****************************
    if (V->abcd && V->aibc) {
  // copy from Vabcd including data
      Tensor<> Cabcd(V->abcd);
  //  Build Chi_abcd intermediate
      Cabcd["abcd"] -= (*V)["akcd"] * (*T)["bk"];
      Cabcd["abcd"] -= (*V)["dcbk"] * (*T)["ak"];

  //  Contract Chi_abcd with T2 Amplitudes
      Rabij["abij"] += Cabcd["abcd"] * (*T)["cdij"];
      Rabij["abij"] += Cabcd["abcd"] * (*T)["ci"] * (*T)["dj"];
    } else {
    }
//*****************************
//@Felix stop slicing here
//*****************************

    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    Bivar_Function<> fDivide(&divide<double>);
    T->abij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);


//********************************************************************************
//  T1 amplitude equations
//********************************************************************************

    LOG(0) << "Solving restricted T1 CCSD Amplitude Equations:" << std::endl;

    Rai["ai"] = Fai["ai"];
    
    Rai["ai"] -= 2.0 * Fai["ck"] * (*T)["ak"] * (*T)["ci"];
    
    Rai["ai"] += Kac["ac"] * (*T)["ci"];

    Rai["ai"] -= Kki["ki"] * (*T)["ak"];

//Build Kck
    Kck["ck"] = Fai["ck"];
    Kck["ck"] += 2.0 * (*V)["cdkl"] * (*T)["dl"];
    Kck["ck"] -= (*V)["cdlk"] * (*T)["dl"];

    Rai["ai"] += 2.0 * Kck["ck"] * (*T)["caki"];
    
    Rai["ai"] -= Kck["ck"] * (*T)["caik"];
    
    Rai["ai"] += Kck["ck"] * (*T)["ci"] * (*T)["ak"];
    
    Rai["ai"] += 2.0 * (*V)["akic"] * (*T)["ck"];
    
    Rai["ai"] -= (*V)["akci"] * (*T)["ck"];

    Rai["ai"] += 2.0 * (*V)["cdak"] * (*T)["cdik"];

    Rai["ai"] -= (*V)["dcak"] * (*T)["cdik"];
  
    Rai["ai"] += 2.0 * (*V)["cdak"] * (*T)["ci"] * (*T)["dk"];

    Rai["ai"] -= (*V)["dcak"] * (*T)["ci"] * (*T)["dk"];

    Rai["ai"] -= 2.0 * (*V)["lkci"] * (*T)["ackl"];

    Rai["ai"] += (*V)["klci"] * (*T)["ackl"];

    Rai["ai"] -= 2.0 * (*V)["lkci"] * (*T)["ak"] * (*T)["cl"];

    Rai["ai"] += (*V)["klci"] * (*T)["ak"] * (*T)["cl"];
//  Calculate Kki
  
//    Rai["ai"] += Kki["ac"] * (*T)["ak"];
    
//   (*V)["acik"] * (*T)["cbkj"];


//    Rabij["abij"] += 2.0 * (*V)["acik"] * (*T)["cbkj"];
//    Cabij["abij"] =  2.0 * (*V)["cbkj"] * (*T)["acik"];
//    Rabij["abij"] += Cabij["abij"];
    //(*V)["cbkj"]*(*T)["acjk"];
//    Rabij["abij"] += 2.0 * Cabij["acik"] * (*T)["cbkj"];


    Dai["ai"] += (*V)["i"];
    Dai["ai"] -= (*V)["a"];

    T->ai->contract(1.0, Rai, "ai", Dai, "ai", 0.0, "ai", fDivide);


  }
}


void Cc4s::iterateMp2() {
  {
    int syms[] = {NS, NS, NS, NS};
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = 0.5 * Dabij["abij"];

    Bivar_Function<> fDivide(&divide<double>);
    T->abij->contract(1.0, *V->abij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
}

void Cc4s::iterateCcsd() {
  int no = chiReal->no;
  int nv = chiReal->nv;
  Tensor<> T21 = Tensor<>(T->abij);
  // NOTE: ctf double counts if lhs tensor is AS
  T21["abij"] += 0.5 * (*T)["ai"] * (*T)["bj"];
  Tensor<> tZabij = Tensor<>(V->abij);

  if (!V->abcd) {
    for (int b(0); b < chiReal->nv; b += chiReal->no) {
  //    for (int a(b); a < nv; a += no) {
      for (int a(0); a < nv; a += no) {
        LOG(0) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
        Tensor<> Vxycd(V->getSlice(a, b));
        Vxycd.set_name("Vxycd");
        int na(Vxycd.lens[0]), nb(Vxycd.lens[1]);
        int Tbegin[] = {0, 0, 0, 0};
        int lens[] = {na, nb, no, no};
//        int syms[] = {Vxycd.sym[0], NS, AS, NS};
        int syms[] = {NS, NS, NS, NS};
        Tensor<> Txyij(4, lens, syms, *world, "Txyij", Vxycd.profile);
        Txyij["xyij"] = Vxycd["xyef"] * T21["efij"];

        int tzBegin[] = {a, b, 0, 0};
        int tzEnd[] = {a+na, b+nb, chiReal->no, chiReal->no};
        tZabij.slice(
          tzBegin,tzEnd,1.0, Txyij,Tbegin,lens,0.5
        );
	//TODO: if a!=b remove double counting of ctf, use b>=a loop
      }
    }
  } else {
    tZabij["abij"] += 0.5 * (*V)["abef"] * T21["efij"];
  }

  {
    int syms[] = {SH, NS, SH, NS};
    Tensor<> Dabij(4, V->abij->lens, syms, *world, "Dabij");
    Dabij["abij"] += (*V)["i"];
    Dabij["abij"] += (*V)["j"];
    Dabij["abij"] -= (*V)["a"];
    Dabij["abij"] -= (*V)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = 0.5 * Dabij["abij"];

    Bivar_Function<> fDivide(&divide<double>);
    T->abij->contract(1.0, tZabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
  }
} 


World *Cc4s::world;
Options *Cc4s::options;
Chi *Cc4s::chiReal, *Cc4s::chiImag;
CoulombIntegrals *Cc4s::V;
Amplitudes *Cc4s::T;


int main(int argumentCount, char **arguments) {
  MPI_Init(&argumentCount, &arguments);

  try {
    Cc4s::world = new World(argumentCount, arguments);
    Cc4s::options = new Options(argumentCount, arguments);
    Log::setRank(Cc4s::world->rank);
    Log::setLogStream(
      new LogStream(Cc4s::options->logFile, Cc4s::options->logLevel)
    );
    Cc4s cc4s;
    cc4s.run();
  } catch (DetailedException *cause) {
    LOG(0) << std::endl << cause->getMessage() << std::endl;
  }
  Log::setLogStream(nullptr);
  MPI_Finalize();
  return 0;
}

