###########################################################
# Makefile
###########################################################
### Variables ###
#Created CSG
# Use GCC compiler


CC =  gcc
CFLAGS = -Wall 
OBJECTS = ODE.o gillespie.o
LDFLAGS =  -lgsl -lgslcblas -lm 
FILE=main




$(FILE): $(OBJECTS) $(FILE).c
	$(CC) $(CFLAGS) $(OBJECTS) -o $(FILE) $(FILE).c $(LDFLAGS) $(DEBUG)

%.o : %.c
	$(CC) $(CFLAGS) -c $< $(DEBUG)

run: $(FILE)
	./$(FILE)

clean:
	rm -fv *.o *.tab *~ $(FILE)

