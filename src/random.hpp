#include <random>

namespace wtl {

template <class URBG>
bool bernoulli(double p, URBG& engine) {
  return std::generate_canonical<double, std::numeric_limits<double>::digits>(engine) < p;
}

template <class URBG>
bool randombit(URBG& engine) {
  return bernoulli(0.5, engine);
}

}
