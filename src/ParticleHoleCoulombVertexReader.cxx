#include <ParticleHoleCoulombVertexReader.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <fstream>

using namespace cc4s;

char const *BinaryFtodReader::Header::MAGIC = "cc4sFTOD";
char const *BinaryFtodReader::Chunk::REALS_MAGIC = "FTODreal";
char const *BinaryFtodReader::Chunk::IMAGS_MAGIC = "FTODimag";
char const *BinaryFtodReader::Chunk::REALSIA_MAGIC = "FTIAreal";
char const *BinaryFtodReader::Chunk::IMAGSIA_MAGIC = "FTIAimag";
char const *BinaryFtodReader::Chunk::EPSILONS_MAGIC = "FTODepsi";

ParticleHoleCoulombVertexReader::ParticleHoleCoulombVertexReader(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  
}

ParticleHoleCoulombVertexReader::~ParticleHoleCoulombVertexReader() {
}

/**
 * \brief Reads the Fourier transformed overlap densities from disk.
 */
void ParticleHoleCoulombVertexReader::run() {
  std::string fileName(getStringArgument("file"));
  LOG(0) <<
    "Reading particle hole Coulomb vertex from file " << fileName << " ...";
  std::ifstream file(fileName, std::ios::binary|std::ios::in);
  if (!file.is_open()) throw new Exception("Failed to open FTODDUMP file");
  // read header
  Header header;
  file.read(reinterpret_cast<char *>(&header), sizeof(header));
  if (strncmp(header.magic, Header::MAGIC, sizeof(header.magic)) != 0)
    throw new Exception("Invalid file format of FTODDUMP file");
  nG = header.nG;
  no = header.no;
  nv = header.nv;
  np = no+nv;

  // allocate output tensors
  TensorData *aiCoulombVertexRealData(
    getTensorDataArgument("aiCoulombVertexReal")
  );
  TensorData *aiCoulombVertexImagData(
    getTensorDataArgument("aiCoulombVertexImag")
  );
  TensorData *epsData(getTensorDataArgument("eps"));
  int vertexLens[] = { nG, nv, no };
  int vertexSyms[] = { NS, NS, NS };
  aiCoulombVertexRealData->value = new Tensor<>(
    3, vertexLens, vertexSyms, *Cc4s::world, "SvRgai"
  );
  aiCoulombVertexImagData->value = new Tensor<>(
    3, vertexLens, vertexSyms, *Cc4s::world, "SvIgai"
  );
  epsData->value = new Vector<>(np, *Cc4s::world, "epsp");
  // FIXME: continue here ...

  Chunk chunk;
  while (file.read(reinterpret_cast<char *>(&chunk), sizeof(chunk))) {
    if (strncmp(chunk.magic, Chunk::REALS_MAGIC, sizeof(chunk.magic)) == 0) {
      readChiChunk(file, Cc4s::chiReal);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGS_MAGIC, sizeof(chunk.magic)) == 0) {
      readChiChunk(file, Cc4s::chiImag);
    } else
    if (strncmp(chunk.magic, Chunk::REALSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(0) << "Found ia chunk";
      //readChiChunk(file, Cc4s::chiIAReal);
      readChiAiChunkBlocked(file, Cc4s::chiAiReal);
    } else
    if (strncmp(chunk.magic, Chunk::IMAGSIA_MAGIC, sizeof(chunk.magic)) == 0) {
      LOG(0) << "Found ia chunk";
      //readChiChunk(file, Cc4s::chiIAImag);
      readChiAiChunkBlocked(file, Cc4s::chiAiImag);
    } else
    if (strncmp(chunk.magic, Chunk::EPSILONS_MAGIC, sizeof(chunk.magic)) == 0) {
      readEpsChunk(file);
    }
  }
  file.close();
  LOG(0) << " OK" << std::endl;
}


// TODO: use several write calls instead of one big to reduce int64 requirement
void BinaryFtodReader::readChiAiChunkBlocked(
  std::ifstream &file, ChiAi *chiAi
) {
  // TODO: separate distribution from reading
  // allocate local indices and values of the chi tensors
  int64_t nvPerNode(nv / chiAi->nv);
  int64_t nvLocal(
    Cc4s::world->rank+1 < chiAi->nv ?
      nvPerNode : nv - Cc4s::world->rank * nvPerNode
  );
  int64_t nvToSkipBefore(Cc4s::world->rank * nvPerNode);
  int64_t nvToSkipAfter(nv - (nvToSkipBefore + nvLocal));
  double *values(new double[nvLocal*no*nG]);
  int64_t *indices(new int64_t[nvLocal*no*nG]);
  file.seekg(sizeof(double)*nvToSkipBefore*no*nG, file.cur);
  file.read(reinterpret_cast<char *>(values), sizeof(double)*nvLocal*no*nG);
  file.seekg(sizeof(double)*nvToSkipAfter*no*nG, file.cur);
  for (int64_t i(0); i < nvLocal*no*nG; ++i) {
    indices[i] = i + nvToSkipBefore*no*nG;
  }
  chiAi->gai->write(nvLocal*no*nG, indices, values);
  delete[] values; delete[] indices;
}


void BinaryFtodReader::readEpsChunk(std::ifstream &file) {
  // allocate local indices and values of eigenenergies
  double *iValues(new double[no]);
  double *aValues(new double[nv]);
  int64_t *iIndices(new int64_t[no]);
  int64_t *aIndices(new int64_t[nv]);

  if (Cc4s::world->rank == 0) {
    file.read(reinterpret_cast<char *>(iValues), no*sizeof(double));
    for (int i(0); i < no; ++i) iIndices[i] = i;
    file.read(reinterpret_cast<char *>(aValues), nv*sizeof(double));
    for (int a(0); a < nv; ++a) aIndices[a] = a;
  } else {
    // skip the data otherwise
    file.seekg(sizeof(double)*np, file.cur);
  }
  int64_t iValuesCount(Cc4s::world->rank == 0 ? no : 0);
  int64_t aValuesCount(Cc4s::world->rank == 0 ? nv : 0);
  Cc4s::V->i->write(iValuesCount, iIndices, iValues);
  Cc4s::V->a->write(aValuesCount, aIndices, aValues);
  delete[] iValues; delete[] aValues;
}
