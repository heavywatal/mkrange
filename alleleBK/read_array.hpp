// -*- mode: c++; coding: utf-8 -*-
#ifndef READ_ARRAY_HPP_
#define READ_ARRAY_HPP_

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <cerrno>

/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

class Fios: public std::fstream{
  public:
	explicit Fios(
			const std::string& filename, 
			const std::ios::openmode mode = std::ios::in | std::ios::out | std::ios::binary
            ) : std::fstream(filename.c_str(), mode)
    {
        if (this->fail()) {
            std::cerr << filename << ": " << std::strerror(errno) << std::endl;
            exit(errno);
        }
    };
};

class Fin: public Fios{
  public:
	explicit Fin(
			const std::string& filename,
			const std::ios::openmode mode = std::ios::binary
			)
    : Fios(filename, mode | std::ios::in){}
	
	std::basic_istream<char>& readline(std::string* buffer, const char delimiter='\n'){
		return std::getline(*this, *buffer, delimiter);
	}
	
	template<class T>
	std::basic_istream<char>& readlines(T* lines, const char delimiter='\n'){
		std::string line;
		while (readline(&line, delimiter)) {
			if (!line.empty()) lines->push_back(line);
		}
		return *this;
	}
	
	std::basic_istream<char>& read(std::string* buffer, const char delimiter='\0'){
		return readline(buffer, delimiter);
	}
	
	std::string readline(const char delimiter='\n'){
		std::string buffer;
		readline(&buffer, delimiter);
		return buffer;
	}
	
	std::vector<std::string> readlines(const char delimiter='\n'){
		std::vector<std::string> lines;
		readlines(&lines, delimiter);
		return lines;
	}
	
	std::string read(const char delimiter='\0'){ return readline(delimiter); }
};

inline std::vector<std::string> split(const std::string& src, const std::string& delimiter=" \t\n") {
    std::vector<std::string> dst;
    for (size_t delim_pos=0, start=0; delim_pos!=src.npos;) {
        delim_pos = src.find_first_of(delimiter, start);
        const int substr_size(delim_pos - start);
        if (substr_size) {dst.push_back(src.substr(start, substr_size));}
        ++(start = delim_pos);
    }
    return dst;
}

inline std::vector<int> split_int(const std::string& src, const std::string& delimiter=" \t\n") {
    std::vector<int> dst;
    for (size_t delim_pos=0, start=0; delim_pos!=src.npos;) {
        delim_pos = src.find_first_of(delimiter, start);
        const int substr_size(delim_pos - start);
        if (substr_size) {dst.push_back(std::atoi(src.substr(start, substr_size).c_str()));}
        ++(start = delim_pos);
    }
    return dst;
}

inline std::vector<std::vector<int> > read_int_array(const std::string& filename) {
    std::vector<std::string> vs(Fin(filename).readlines());
    const size_t n(vs.size());
    std::vector<std::vector<int> > vvi(n);
    for (unsigned int i=0; i<n; ++i) {
        vvi[i] = split_int(vs[i]);
    }
    return vvi;
}


/////////1/////////2/////////3/////////4/////////5/////////6/////////7/////////8

#endif /* READ_ARRAY_HPP_ */
