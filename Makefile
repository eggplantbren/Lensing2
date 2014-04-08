CFLAGS = -O2 -DNDEBUG -Wall -Wextra -ansi -pedantic
LIBS = -lrjobject -ldnest3 -lgsl -lgslcblas -lboost_thread -lboost_system

default:
	g++ -c $(CFLAGS) Lenses/*.cpp Sources/*.cpp
	g++ -o main main.cpp $(LIBS)
	rm -rf *.o

