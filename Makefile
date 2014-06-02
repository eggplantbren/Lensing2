CFLAGS = -O2 -DNDEBUG -Wall -Wextra -ansi -pedantic
LIBS = -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system

default:
	g++ -c $(CFLAGS) *.cpp Lenses/*.cpp Sources/*.cpp
	g++ -o main *.o $(LIBS)
	rm -rf *.o

