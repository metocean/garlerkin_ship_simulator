headers = *.h
#CC = mpicc
#CC = nvcc -ccbin /usr/bin/gcc-4.6
CC = gcc
#CFLAGS = -g -Wall
#CFLAGS = -g - Wall -O
#LFLAGS = `pkg-config --cflags --libs glib-2.0` `pkg-config --cflags --libs gsl`  `pkg-config --cflags --libs gts` -lm
LFLAGS = -D_GNU_SOURCE -I/usr/include/superlu -L/usr/lib -lsuperlu_4.3 -lblas `pkg-config --cflags --libs glib-2.0` `pkg-config --cflags --libs gsl` `pkg-config --cflags --libs plasma` -lgfortran -lnetcdf -linterp2d -lm
#CFLAGS = -g -Wall $(LFLAGS)
CFLAGS = -g $(LFLAGS)
COMPILE = $(CC) $(CFLAGS) -c
#COMPILE = $(CC) -O -c
#OBJFILES := $(patsubst %.c,%.o,$(wildcard *.c))
OBJFILES = ship.o surface.o ncinput.o boundaries.o motion.o hull.o b4.o linearproblem.o patch.o structures.o
HEADERS = ship.h surface.h ncinput.c boundaries.h motion.h hull.h b4.h linearproblem.h patch.h structures.h
SRCS_C = ship.c hull.c surface.c ncinput.h boundaries.c motion.c b4.c linearproblem.c patch.c structures.c
SRCS = $(SRCS_C) $(HEADERS)

all: gss.exe

#gss.exe : ship.o surface.o
#	$(COMPILE) -o gss.exe ship.o surface.o

gss.exe: $(OBJFILES)
	$(CC) $(OBJFILES) $(LFLAGS) -o gss.exe

structures.o : structures.c structures.h
	$(COMPILE) structures.c

patch.o : patch.c patch.h structures.h
	$(COMPILE) patch.c

b4.0 : b4.c b4.h patch.h structures.h
	$(COMPILE) b4.c

linearproblem.o : linearproblem.c linearproblem.h structures.h
	$(COMPILE) linearproblem.c

#hypre.o : hypre.c hypre.h linearproblem.h structures.h
#	$(COMPILE) hypre.c

motion.o : motion.c motion.h hull.h structures.h
	$(COMPILE) motion.c

boundaries.o : boundaries.c boundaries.h structures.h
	$(COMPILE) boundaries.c

surface.o : surface.c surface.h boundaries.h linearproblem.h structures.h
	$(COMPILE) surface.c

hull.o : hull.c hull.h motion.h patch.h structures.h
	$(COMPILE) hull.c

#surface.h boundaries.h linearproblem.h

ncinput.o : ncinput.c ncinput.h patch.h structures.h
	$(COMPILE) ncinput.c

#potential.o : potential.c potential.h surface.h boundaries.h linearproblem.h structures.h
#	$(COMPILE) potential.c

ship.o : ship.c ship.h ncinput.h potential.h surface.h boundaries.h linearproblem.h structures.h
	$(COMPILE) ship.c

clean:
	rm -f *.o

TAGS : $(SRCS)
	etags $(SRCS)

tags : $(SRCS)
	etags $(SRCS)
