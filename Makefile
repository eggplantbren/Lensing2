CFLAGS = -O2 -Wall -Wextra -ansi -pedantic

default:
	g++ -c $(CFLAGS) Lenses/*.cpp

