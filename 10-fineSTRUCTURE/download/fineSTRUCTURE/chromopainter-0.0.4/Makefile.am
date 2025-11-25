ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = chromopainter

OPTIMIZATION = -O3 -lm -lz -Wall -Wno-write-strings -mfpmath=sse -msse -msse2 -funroll-loops -fomit-frame-pointer -ftree-vectorize -funsafe-math-optimizations
chromopainter_CXXFLAGS = $(GSL_CFLAGS) $(OPTIMIZATION) -Wall
chromopainter_LDFLAGS = $(OPTIMIZATION)
chromopainter_LDADD = $(GSL_LIBS) -lz
chromopainter_SOURCES = ChromoPainterMain.c

dist_data_DATA = ChromoPainterMutEM.c ChromoPainterMutEM.h
