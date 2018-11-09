main.exe: mesolver.o test.o
	g++ test.o mesolver.o -o test.exe -larmadillo

test.o: test.cpp
	g++ -c -std=c++11 main.cpp -larmadillo

mesolver.o: mesolver.cpp mesolver.h
	g++ -c -std=c++11 mesolver.cpp -larmadillo

clean:
	rm *.o *.exe

run: test.exe
	./test.exe
