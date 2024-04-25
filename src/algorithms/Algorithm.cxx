#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <regex>

#include <algorithms/Algorithm.hpp>
#include <NewData.hpp>
#include <math/Complex.hpp>
#include <DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>

using namespace sisi4s;

Algorithm::Algorithm(Arguments const &in_, Arguments const &out_)
    : in(in_)
    , out(out_) {}

Algorithm::~Algorithm() {}

/**
 * \brief The dryRun estimates resource consumption, especially
 * memory and processor time.
 */
void Algorithm::dryRun() {
  LOG(0, getName()) << "dry run not implemented" << std::endl;
}

AlgorithmFactory::AlgorithmMap *AlgorithmFactory::algorithmMap;

std::string AlgorithmFactory::normalize_name(std::string const &name) {
  std::vector<std::string> words;
  const std::regex words_regex("[^\\s_-]+");
  const auto words_begin =
                 std::sregex_iterator(name.begin(), name.end(), words_regex),
             words_end = std::sregex_iterator();

  for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
    std::string word = i->str();
    word[0] = std::toupper(word[0]);
    words.push_back(word);
  }
  std::string new_name;
  for (auto const &word : words) new_name += word;
  return new_name;
}
