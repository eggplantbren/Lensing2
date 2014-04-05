CFLAGS = -O2 -DNDEBUG -Wall -Wextra -ansi -pedantic
LIBS = -lgsl -lgslcblas -lboost_thread -lboost_system -ldnest3

default:
	g++ -c $(CFLAGS) Lenses/*.cpp
	g++ -o main main.cpp $(LIBS)
	rm -rf *.o


