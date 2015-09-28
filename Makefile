#CFLAGS = -O2 -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -ansi -pedantic
CFLAGS = -m64 -O3 -flto -march=native -funroll-loops -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -ansi -pedantic
LIBS = -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system -lgfortran

default:
	gfortran -m64 -O3 -flto -march=native -funroll-loops -O3 -c fastell.f
	g++ -c $(CFLAGS) *.cpp Lenses/*.cpp Sources/*.cpp
	g++ -o main *.o $(LIBS)
	rm -rf *.o

