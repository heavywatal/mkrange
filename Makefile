CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2 -g

all: a.out
	./$< steep

a.out: src/main.o src/uball.o
	$(LINK.cpp) $^

src/main.o: src/uball.h src/random.hpp
src/uball.o: src/uball.h src/random.hpp src/read_array.hpp
