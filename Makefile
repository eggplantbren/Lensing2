CFLAGS = -O3 -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -ansi -pedantic
LIBS = -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system -lgfortran

default:
	gfortran -O3 -c fastell.f
	g++ -c $(CFLAGS) *.cpp Lenses/*.cpp Sources/*.cpp
	g++ -o main *.o $(LIBS)
	rm -rf *.o

