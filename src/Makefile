CC = g++
CFLAGS = -std=c++11 -O3 -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -pedantic
OBJECTS = BasicCircular.o Blobby.o Data.o MyModel.o BasicUniform.o BlobbySPEMD.o fastell.o PSF.o
LIBS = -ldnest4 -lpthread -lgfortran

default:
	gfortran -O3 -c fastell.f
	$(CC) $(CFLAGS) -I$(DNEST4_PATH) -c *.cpp Sources/Blobby.cpp Lenses/BlobbySPEMD.cpp
	$(CC) -L$(DNEST4_PATH)/DNest4/code -o Lensing2 main.o $(OBJECTS) $(LIBS)
	$(CC) -L$(DNEST4_PATH)/DNest4/code -o loadrow loadrow.o $(OBJECTS) $(LIBS)
	rm -rf *.o

