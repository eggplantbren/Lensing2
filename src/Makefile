CC = g++
CFLAGS = -std=c++14 -I /usr/include/eigen3 -O3 -march=native -DEIGEN_NO_DEBUG -DARMA_NO_DEBUG -DNDEBUG -Wall -Wextra -pedantic
OBJECTS = BasicCircular.o Blobby.o Data.o MyModel.o LensBlobConditionalPrior.o BlobbySPEMD.o NoiseModel.o NoiseModel2.o fastell.o PSF.o
LIBS = -ldnest4 -lpthread -lgfortran -lyaml-cpp

default:
	make -C DNest4/code noexamples
	gfortran -O3 -march=native -DNDEBUG -c fastell.f
	$(CC) $(CFLAGS) -I . -c *.cpp Sources/Blobby.cpp Lenses/BlobbySPEMD.cpp
	$(CC) $(CFLAGS) -I . -I CorrelatedNoise -c CorrelatedNoise/NoiseModel*.cpp
	ar rcs liblensing2.a $(OBJECTS)
	$(CC) -L DNest4/code -o Lensing2 main.o $(OBJECTS) $(LIBS)
	rm -rf *.o

