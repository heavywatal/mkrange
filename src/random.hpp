#pragma once
#ifndef MKRANGE_RANDOM_HPP_
#define MKRANGE_RANDOM_HPP_

#include <random>
#include <limits>

namespace wtl {

template <class URBG>
bool bernoulli(double p, URBG& engine) {
  return std::generate_canonical<double, std::numeric_limits<double>::digits>(engine) < p;
}

template <class URBG>
bool randombit(URBG& engine) {
  return engine() & 1u;
}

}

#endif// MKRANGE_RANDOM_HPP_
