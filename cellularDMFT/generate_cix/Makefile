C++ = g++
CFLAGS = -O3
LIBS = -framework Accelerate 

all : hub

hub : hub.o
	$(C++) $(CFLAGS) -o $@ $? $(LIBS)

clean :
	rm -f hub hub.o


.SUFFIXES : .cc
.cc.o:
	$(C++) $(CFLAGS) -c $<


