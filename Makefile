CXX = g++
CXXFLAGS = -O2 -Wall -fopenmp

msfinder: msfinder.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

test: msfinder
	for i in 3 4 5 6 7 8 9; do \
		./msfinder -n $$i; \
	done

all: msfinder test

.PHONY: clean
clean:
	rm -f msfinder