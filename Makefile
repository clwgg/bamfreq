CC= gcc
CFLAGS= -g -Wall #-O2
INC= 
LIB= -Lhtslib -lhts -lpthread -lz -lm
OBJS= freq.o utils.o

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(INC)

freq: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(INC) $(LIB)

clean:
	rm -f $(OBJS)

submodules:
	cd htslib && make libhts.a ; cd ..

# DO NOT DELETE THIS LINE -- make depend depends on it.

freq.o utils.o: utils.h

