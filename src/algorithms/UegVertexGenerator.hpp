/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef UEG_VERTEX_GENERATOR_DEFINED
#define UEG_VERTEX_GENERATOR_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/TensorExpression.hpp>
#include <Complex.hpp>
#include <Integer.hpp>
#include <TensorSet.hpp>


namespace cc4s {
  using ivec  = array<Integer<>,3>;
  using dvec  = array<Real<>,4>;

  class UegVertexGenerator: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UegVertexGenerator)
    UegVertexGenerator(std::vector<Argument> const &argumentList);
    void run();

    template<typename F>
    void run();

  protected:
    Real<> evalMadelung(double volume);
    Real<> Vijji(const dvec a, const dvec b, const Real<> v);

    bool halfGrid;
    Natural<> No, Nv, NF;
    Real<> rs, madelung;
  };
}

#endif

