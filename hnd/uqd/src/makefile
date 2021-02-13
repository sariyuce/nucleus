include makefile.in

TARGET     = undirected_quark

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

undirected_quark: main.o ktruss.o hierarchy.o util.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o ktruss.o hierarchy.o util.o bucket.o graph.o  $(LIBS)
	/bin/rm *.o
	
clean: 
	/bin/rm -f ../$(TARGET) *.o 
