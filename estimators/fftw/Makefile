CC    = gcc #run without -mp-6 on splinter
FLAGS        = #-std=gnu99 -Iinclude -I/usr/local/include
CFLAGS       = -fPIC -g -pedantic -Wall -Wextra -ggdb3 -O3 -ffast-math -fopenmp -I/usr/include/gsl/ -I/usr/lib/x86_64-linux-gnu/ -I/home/ganjalf/anaconda2/include -I/usr/local/include/

LFLAGS       =  -lgomp -lm -fopenmp -L/usr/lib/x86_64-linux-gnu/ -L/usr/include/gsl/ -lgsl -lgslcblas -L/usr/local/lib -lfftw3 -lfftw3_threads -lfftw3_omp
LDFLAGS      = -shared

DEBUGFLAGS   = -O0 -D _DEBUG
RELEASEFLAGS = -O2 -D NDEBUG -combine -fwhole-program

TARGET  = dft.so
SOURCES = $(wildcard *.c)
HEADERS = $(wildcard *.h)
OBJECTS = $(SOURCES:.c=.o)


all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(CFLAGS) -o $(TARGET) $(OBJECTS) $(LFLAGS)
