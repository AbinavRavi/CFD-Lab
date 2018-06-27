CC = gcc
CFLAGS = -Wall -pedantic -Werror
.c.o:  ; $(CC) -c $(CFLAGS) $<

OBJ = 	helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
	sor.o\
      	main.o\
      	visual.o

rm_OBJ = helper.o\
      	init.o\
      	boundary_val.o\
      	uvp.o\
	sor.o\
      	main.o\
      	visual.o\
	sim\
	*.log
	

all:  $(OBJ)
	$(CC) $(CFLAGS) -o sim $(OBJ) -lm

%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	rm $(rm_OBJ) 

helper.o      : helper.h 
init.o        : helper.h init.h 
boundary_val.o: helper.h boundary_val.h 
uvp.o	      : uvp.h boundary_val.h
visual.o      : helper.h
sor.o	      : sor.h boundary_val.h
main.o        : helper.h init.h boundary_val.h uvp.h visual.h

