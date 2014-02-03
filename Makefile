TARGET = fpldpc
CC = g++
CFLAGS = -O3 

OBJECTS = ArrayLDPC.o

SOURCE = ArrayLDPC.cpp rngs.cpp rvgs.cpp PerfTest.cpp


  
HEADERS = 

all:

clean:
	rm $(TARGET)

compile:
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) -lm

