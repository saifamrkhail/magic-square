all: clean msfinder run 

msfinder:
	g++ -o2 -Wall -fopenmp msfinder.cpp -o msfinder

run:
	./msfinder

clean:
	rm msfinder