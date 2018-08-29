CXX = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wfloat-equal -O2 -march=native -g

all: a.out
	@:

a.out: src/main.o src/uball.o
	$(LINK.cpp) $^

src/main.o: src/uball.h src/global.hpp
src/uball.o: src/uball.h src/global.hpp src/random.hpp src/read_array.hpp

clean:
	$(RM) src/*.o a.out

flat: a.out
	ln -sf Inputfile-$@ Inputfile
	./$< flat TestInput.txt

linear: a.out
	ln -sf Inputfile-$@ Inputfile
	./$< linear

steep: a.out
	ln -sf Inputfile-$@ Inputfile
	./$< steep
