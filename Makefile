.SUFFIXES:
.SUFFIXES: .C .o


#CFLAGS = -U_FORTIFY_SOURCE -g -pg -Wall -Wextra -pedantic
CFLAGS = -U_FORTIFY_SOURCE -O3 -mtune=generic


# GNU Scientific Library (GSL) support:
GSL_C := $(shell gsl-config --cflags)
GSL_L := $(shell gsl-config --libs)


# Object files:
OBJS = WLSim.o      


# Target file:
TARGET = WLSim

.C.o:
	$(CXX) $(CFLAGS) $(GSL_C) -c $*.C

$(TARGET): $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(GSL_L)


WLSim.o: WLSim.C 

cleanup:
	rm -f *.o $(TARGET)
