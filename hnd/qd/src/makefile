include makefile.in

TARGET     = quark

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)
quark: main.o hierarchy.o util.o bucket.o graph.o cycle.o acyclic.o outp.o cyclep.o inp.o cyclepp.o
	$(LD) $(LDFLAGS) -o ../$@ main.o hierarchy.o util.o bucket.o graph.o cycle.o acyclic.o outp.o cyclep.o inp.o cyclepp.o $(LIBS)
	/bin/rm *.o

clean: 
	/bin/rm -f ../$(TARGET) *.o 
