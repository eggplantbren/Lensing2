CC = g++
CFLAGS = -std=c++11 -O3 -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -pedantic
LIBS = -ldnest4 -lpthread -lgfortran

default:
	gfortran -O3 -c fastell.f
	$(CC) $(CFLAGS) -I$(DNEST4_PATH) -c *.cpp Sources/Blobby.cpp Lenses/BlobbySPEMD.cpp
	$(CC) -L$(DNEST4_PATH)/DNest4/code -o main *.o $(LIBS)
	rm -rf *.o

