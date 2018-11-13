INC = -I ~/Eigen/Eigen

test.exe: mesolver.o test.o
	g++ $(INC) test.o mesolver.o -o test.exe -larmadillo

test.o: test.cpp
	g++ $(INC) -c -std=c++11 test.cpp -larmadillo

mesolver.o: mesolver.cpp mesolver.h
	g++ $(INC) -c -std=c++11 mesolver.cpp -larmadillo

clean:
	rm *.o *.exe

run: test.exe
	./test.exe
