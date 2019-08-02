VPATH=.
SF_HEADER=funcHO.h
C_MATH_LIB=-lm
CC=gcc
CFLAGS=-Wall -pedantic #-fopenmp 

run: harm_osc


harm_osc : harm_osc.o funcHO.o init.o  Makefile
	${CC} -o harm_osc harm_osc.o funcHO.o init.o  $(C_MATH_LIB) #-fopenmp

harm_osc.o : harm_osc.c Makefile
	${CC} ${CFLAGS} -c harm_osc.c 

init.o : init.c init.h Makefile
	${CC} ${CFLAGS} -c init.c 

funcHO.o : funcHO.c funcHO.h Makefile
	${CC} ${CFLAGS} -c funcHO.c 


clean :
	rm harm_osc.o funcHO.o init.o







