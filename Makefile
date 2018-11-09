main.exe: mesolver.o main.o
	g++ main.o mesolver.o -o main.exe -larmadillo

main.p: main.cpp
	g++ -c -std=c++11 main.cpp -larmadillo

mesolver.o: mesolver.cpp mesolver.h
	g++ -c -std=c++11 mesolver.cpp -larmadillo

clean:
	rm *.o *.exe

run: main.exe
	./main.exe
