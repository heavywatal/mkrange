// -*- mode: c++; coding: utf-8 -*-
#ifndef READ_ARRAY_HPP_
#define READ_ARRAY_HPP_

#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>

namespace wtl {

inline std::vector<std::string> readlines(std::istream& ist) {
    std::vector<std::string> lines;
    std::string buffer;
    while (std::getline(ist, buffer)) {
        lines.push_back(buffer);
    }
    return lines;
}

template <class T>
inline std::vector<T> split(const std::string& src) {
    std::stringstream sst(src);
    return std::vector<T>(std::istream_iterator<T>(sst), std::istream_iterator<T>());
}

} // namespace wtl

inline std::vector<std::vector<int> > read_int_array(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    ifs.exceptions(std::ios_base::badbit);
    std::vector<std::string> lines(wtl::readlines(ifs));
    const size_t n = lines.size();
    std::vector<std::vector<int> > output;
    output.reserve(n);
    for (size_t i=0; i<n; ++i) {
        if (lines[i].empty()) continue;
        output.push_back(wtl::split<int>(lines[i]));
    }
    return output;
}

#endif /* READ_ARRAY_HPP_ */
