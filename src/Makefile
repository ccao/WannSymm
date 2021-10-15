include ../make.sys

#ifeq ($(shell uname -s), Darwin)
#	include ../make.sys.osx
#else
#	include ../make.sys
#endif


objects = readinput.o rotate_basis.o rotate_ham.o rotate_orbital.o rotate_spinor.o vector.o wanndata.o wannorb.o readsymm.o usefulio.o matrix.o bndstruct.o

all:main

wannsymm: main.c $(objects)
	@echo Making wannsymm.x 
	@$(CC) $(CFLAGS) -o wannsymm.x main.c $(INCLUDES) $(objects) $(LIBS) -lm

main:wannsymm
	@if [ ! -d ../bin ]; then mkdir ../bin; fi
	mv wannsymm.x ../bin/

clean:
	rm -rf *.o *.mod *.x
	rm ../bin/*.x
	
