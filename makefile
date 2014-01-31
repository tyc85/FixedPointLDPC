rcdir = .
prefix = /usr/local
exec_prefix=${prefix}

CC=g++
CFLAGS=-g -O2 -I. -fPIC -Wall
LIBS=arrayldpc.o perftest.o timetrial.o rvgs.o rngs.o wrapper.o

SLIBS=libldpc.so

all: $(SLIBs)

clean:
        rm -f *.o $(SHARED_LIB)


# for Linux et al
libldpc.so: $(LIBS)
        gcc -shared -Xlinker -soname=$@ -o $@ -Wl,-whole-archive $^ -Wl,-no-whole-archive -lc


arrayldpc.o: ArrayLDPC.cpp ArrayLDPC.h ArrayLDPCMacro.h

#memory.o: Memory.cpp Memory.h
perftest.o: PerfTest.cpp ArrayLDPC.h ArrayLDPCMacro.h

rvgs.o: rvgs.cpp rvgs.h 

rngs.o: rngs.cpp rngs.h

#wrapper.o: Wrapper.cpp ArrayLDPCMacro.h 
