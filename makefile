CC = g++
CFLAGS = -O3 -v -Wall -Wextra -Werror -g
OBJECTS = main.o

lattice : $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o benchmark

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o benchmark
