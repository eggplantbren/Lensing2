CC = g++
CFLAGS = -std=c++14 -Wall -Wextra -pedantic
OPTIM = -O0 -g
LIBS = -llensing2 -ldnest4

default:
	$(CC) -I$(DNEST4_PATH) -L$(DNEST4_PATH)/DNest4/code -I .. $(CFLAGS) -L .. $(OPTIM) -o BasicCircular BasicCircular.cpp $(LIBS)

clean:
	rm -rf BasicCircular

