AUTOMAKE_OPTIONS = subdir-objects

AM_CXXFLAGS = -fmax-errors=1 \
              -D_POSIX_C_SOURCE=200112L     \
              -D__STDC_LIMIT_MACROS         \
              -DFTN_UNDERSCORE=1            \
              "-DSISI4S_DATE=\"${DATE}\""


bin_PROGRAMS = sisi4s

# CORE
sisi4s_SOURCES =                                                 \
./Arguments.cxx                                                  \
./AlgorithmInputSpec.cxx                                         \
./NewData.cxx                                                    \
./algorithms/Algorithm.cxx                                       \
./Sisi4s.cxx                                                     \
./Exit.cxx                                                       \
./Options.cxx                                                    \
./DryTensor.cxx                                                  \
./Parser.cxx                                                     \
./util/Tensor.cxx                                                \
./util/BlacsWorld.cxx                                            \
./util/Emitter.cxx                                               \
./util/Parsing.cxx                                               \
./util/ScaLapackMatrix.cxx                                       \
./util/Log.cxx                                                   \
./util/TensorIo.cxx                                              \
./util/FlopsCounter.cxx                                          \
./util/Timer.cxx                                                 \
./util/BasisSet.cxx                                              \
./algorithms/Exit.cxx                                            \
./algorithms/Nop.cxx

sisi4s_SOURCES +=                                                 \
./util/AngularMomentum.cxx                                       \
./turbomole/MosParser.cxx                                        \
./nwchem/MovecsParser.cxx                                        \
./nwchem/BasisParser.cxx                                         \
./mixers/LinearMixer.cxx                                         \
./mixers/DiisMixer.cxx                                           \
./mixers/Mixer.cxx                                               \
./math/IterativePseudoInverse.cxx                                \
./math/PseudoInverseSvd.cxx                                      \
./math/ComplexTensor.cxx                                         \
./math/MathFunctions.cxx                                         \
./equations/StantonIntermediatesUCCSD.cxx                        \
./equations/sthamiltonian/wabij.cxx                             \
./equations/sthamiltonian/wabcijk.cxx                           \
./equations/sthamiltonian/wabcdijkl.cxx                           \
./equations/SimilarityTransformedHamiltonian.cxx                \
./algorithms/FromComplexTensor.cxx                               \
./algorithms/TensorNorm.cxx                                      \
./algorithms/TensorSum.cxx                                       \
./algorithms/TensorReader.cxx                                    \
./algorithms/TensorWriter.cxx                                    \
./algorithms/GenerateRandomMatrix.cxx                            \
./algorithms/GenerateRandomTensor.cxx                            \
./algorithms/CoulombIntegralsFromVertex.cxx                      \
./algorithms/OneBodyRotation.cxx                                 \
./algorithms/CoulombIntegralsFromRotatedCoulombIntegrals.cxx     \
./algorithms/CcsdPreconditioner.cxx                              \
./algorithms/HartreeFockFromCoulombIntegrals.cxx                 \
./algorithms/ClusterSinglesDoublesAlgorithm.cxx                  \
./algorithms/ClusterSinglesDoublesTriplesQuadruplesAlgorithm.cxx \
./algorithms/CcsdEnergyFromCoulombIntegralsReference.cxx         \
./algorithms/CcsdEnergyFromCoulombIntegrals.cxx                  \
./algorithms/ClusterSinglesDoublesTriplesAlgorithm.cxx           \
./algorithms/CcsdtEnergyFromCoulombIntegrals.cxx                 \
./algorithms/cc4s/DefineHolesAndParticles.cxx                    \
./algorithms/Read.cxx                                            \
./algorithms/FiniteSizeCorrection.cxx                            \
./algorithms/BasisSetExtrapolation.cxx                           \
./algorithms/TensorUnrestricter.cxx                              \
./algorithms/FcidumpReader.cxx                                   \
./algorithms/FcidumpWriter.cxx                                   \
./algorithms/MoReader.cxx                                        \
./algorithms/Mp2NaturalOrbitals.cxx                              \
./algorithms/TensorAntisymmetrizer.cxx                           \
./algorithms/FockMatrixFromCoulombIntegrals.cxx                  \
./algorithms/UccsdAmplitudesFromCoulombIntegrals.cxx             \
./algorithms/UccsdtAmplitudesFromCoulombIntegrals.cxx            \
./algorithms/UccsdtqAmplitudesFromCoulombIntegrals.cxx            \
./algorithms/MeanCorrelationHoleDepth.cxx                        \
./algorithms/CcsdEquationOfMotionDavidson.cxx                    \
./algorithms/PQRSCoulombIntegralsToVertex.cxx \
./algorithms/UegVertexGenerator.cxx

## sisi4s_SOURCES =                                                 \
## ./algorithms/UnrestrictedCIS.cxx                                 \
## ./algorithms/UCcsdIPEquationOfMotionDavidson.cxx                 \
## ./algorithms/UCcsdEAEquationOfMotionDavidson.cxx                 \
## ./math/CanonicalPolyadicDecomposition.cxx                        \
## ./math/PseudoInverseHermitianSvd.cxx                             \
## ./math/RegularizedAlternatingLeastSquares.cxx                    \
## ./algorithms/PseudoInverse.cxx                                   \
## ./algorithms/ParticleHoleCoulombIntegrals.cxx                    \
## ./algorithms/TensorGetMax.cxx                                    \
## ./algorithms/ParticleHoleCoulombVertexSingularVectors.cxx        \
## ./algorithms/NaturalTransitionOrbitalsFromRhoAI.cxx              \
## ./algorithms/ComplexTensorNorm.cxx                               \
## ./algorithms/TensorSlicer.cxx                                    \
## ./algorithms/ParticleHoleCoulombVertexReader.cxx                 \
## ./algorithms/CoulombVertexReader.cxx                             \
## ./algorithms/SliceCoulombVertex.cxx                              \
## ./algorithms/ComplexTensorWriter.cxx                             \
## ./algorithms/ApproximateCoulombVertex.cxx                        \
## ./algorithms/CcsdPerturbativeTriples.cxx                         \
## ./algorithms/UnrestrictedEquationOfMotionSinglesFromRpa.cxx      \
## ./algorithms/UrpaAmplitudesFromCoulombIntegrals.cxx              \
## ./algorithms/Delete.cxx                                          \
## ./algorithms/UccsdtAmplitudesFromCoulombIntegrals.cxx            \
## ./algorithms/TensorReduceijij.cxx                                \
## ./algorithms/CcsdDiagrammaticDecomposition.cxx                   \

if DISABLE_LIBINT
else
sisi4s_SOURCES +=                                   \
./util/Eigen.cxx                          \
./algorithms/XYZReader.cxx                          \
./algorithms/OneBodyFromGaussian.cxx                \
./algorithms/CoulombIntegralsFromGaussian.cxx       \
./algorithms/HartreeFockFromGaussian.cxx            \
./algorithms/NuclearRepulsionEnergy.cxx
endif

sisi4s_CPPFLAGS = $(EXTERNAL_CPPFLAGS)
sisi4s_LDADD = $(EXTERNAL_LDFLAGS)



## ./algorithms/UPerturbativeTriples.cxx                            \
## ./algorithms/DoublesAmplitudesFromVertex.cxx                     \
## ./algorithms/ComplexTensorSum.cxx                                \
## ./algorithms/TensorContraction.cxx                               \
## ./algorithms/ApproximateParticleHoleCoulombVertex.cxx            \
## ./algorithms/PerturbativeTriples.cxx                             \
## ./algorithms/CoulombVertexDecomposition.cxx                      \
## ./algorithms/ParticleHoleCoulombVertexDecomposition.cxx          \
## ./algorithms/DoublesAmplitudesDecomposition.cxx                  \
## ./algorithms/CcsdPerturbativeTriplesComplex.cxx                  \
## ./algorithms/ComplexTensorContraction.cxx                        \
## ./algorithms/ParenthesisTriples.cxx                              \
## ./algorithms/FromComplexTensor.cxx                               \
## ./algorithms/SingleParticleOccupancies.cxx                       \
## ./algorithms/DrccdEnergyFromCoulombIntegrals.cxx                 \
## ./algorithms/CoulombVertexSingularVectors.cxx                    \
## ./algorithms/ClusterSinglesDoublesTriplesQuadruplesAlgorithm.cxx \

##
## These are the only algorithms using Tcc
##
## ./algorithms/CoulombVertexFromFactors.cxx                        \
## ./algorithms/ParticleHoleCoulombVertexFromFactors.cxx            \
