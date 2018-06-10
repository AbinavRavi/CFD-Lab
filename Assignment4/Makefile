INCLUDES = \
	-I. \
	-I$(PRECICE_ROOT)/src/precice

LIBS = \
	-L$(PRECICE_ROOT)/build/last -lprecice\
	-lm

CC = gcc
CFLAGS = -pedantic -Werror $(INCLUDES)
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = \
	helper.o\
	init.o\
	boundary_val.o\
	uvp.o\
	main.o\
	visual.o\
	sor.o \
	precice_adapter.o


all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ) $(LIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o $(LIBS)

clean:
	rm $(OBJ)

helper.o      		: helper.h
init.o        		: helper.h init.h
boundary_val.o		: helper.h boundary_val.h
uvp.o         		: helper.h uvp.h
visual.o      		: helper.h
sor.o	      		: sor.h
precice_adapter.o 	: precice_adapter.h

main.o        : helper.h init.h boundary_val.h uvp.h visual.h sor.h precice_adapter.h
