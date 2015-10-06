.SUFFIXES: .c .o

OBJS = parse.o dictionary.o gen_inc.o fortprintf.o ezxml/ezxml.o

all: parse

ezxml/ezxml.o:
	(cd ezxml; $(CC) -c -g ezxml.c)

parse: $(OBJS)
	$(CC) -g -o $@ $(OBJS)

clean:
	$(RM) *.o ezxml/*.o parse

.c.o:
	$(CC) -c -g $<
