CXX = g++
CXXFLAGS = -O2 -Wall -fopenmp

msfinder: msfinder.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

.PHONY: clean
clean:
	rm -f msfinder