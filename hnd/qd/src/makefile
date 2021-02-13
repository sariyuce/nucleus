include makefile.in

TARGET     = dinuc

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)
dinuc: main.o hierarchy.o util.o bucket.o graph.o cycle_truss.o acyclic_truss.o outp_truss.o cyclep_truss.o inp_truss.o cyclepp_truss.o
	$(LD) $(LDFLAGS) -o ../$@ main.o hierarchy.o util.o bucket.o graph.o cycle_truss.o acyclic_truss.o outp_truss.o cyclep_truss.o inp_truss.o cyclepp_truss.o $(LIBS)
	/bin/rm *.o

clean: 
	/bin/rm -f ../$(TARGET) *.o 
