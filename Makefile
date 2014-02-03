TARGET = fpldpc
CC = g++
CFLAGS = -O3 

OBJECTS = TurboCode.o Noise.o Decoder.o Encoder.o Interleaver.o Array.o SubBlkInt.o

SOURCE = ArrayLDPC.cpp rvgs.cpp Memory.cpp PerfTest.cpp


  
HEADERS = 

all:

clean:
	rm $(TARGET)

compile:
	$(CC) $(CFLAGS) $(SOURCE) -o $(TARGET) -lm

