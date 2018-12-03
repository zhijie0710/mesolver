test2.exe: mesolver.o test.o
	g++ test.o mesolver.o -o test2.exe -larmadillo

test.o: test.cpp
	g++ -c -std=c++11 test.cpp -larmadillo

mesolver.o: mesolver.cpp mesolver.h
	g++ -c -std=c++11 mesolver.cpp -larmadillo

clean:
	rm *.o *.exe

run: test2.exe
	./test2.exe
