#include <random>

namespace wtl {

constexpr double PI = 3.14159265358979323846;

template <class URBG>
bool bernoulli(double p, URBG& engine) {
  return std::generate_canonical<double, std::numeric_limits<double>::digits>(engine) < p;
}

template <class URBG>
bool randombit(URBG& engine) {
  return bernoulli(0.5, engine);
}

template <class URBG>
std::vector<int> randomove(double r, URBG& engine) {
  std::uniform_real_distribution<double> uniform(0.0, 2.0 * PI);
  double angle = uniform(engine);
  int dx = std::round(r * std::cos(angle));
  int dy = std::round(r * std::sin(angle));
  return {dx, dy};
}

}
