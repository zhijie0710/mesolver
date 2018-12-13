test.exe: timeEvolution.o test.o
	g++ test.o timeEvolution.o -o test.exe -larmadillo

test.o: test.cpp
	g++ -c -std=c++11 test.cpp -larmadillo

timeEvolution.o: timeEvolution.cpp timeEvolution.h
	g++ -c -std=c++11 timeEvolution.cpp -larmadillo

clean:
	rm *.o *.exe

run: test.exe
	./test.exe
