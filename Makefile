PROFILE = -fno-omit-frame-pointer -ggdb
FLAGS = -O3 -march=native -lm -Wall $(PROFILE)
CFLAGS = $(FLAGS) -std=c11
CXXFLAGS = $(FLAGS) -std=c++14
.PHONY: clean test
all: strassen strassen_cpp
strassen: strassen.c strassen.h
	$(CC) $(CFLAGS) $< -o $@
strassen_cpp: strassen.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
test: strassen
	time ./$< 2000 < test/2000.in > result
release: strassen.c strassen.h test/1000.in Makefile README.md
	zip -9 strassen.zip $^ 
clean:
	rm -rf strassen result
