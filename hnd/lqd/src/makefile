include makefile.in

TARGET     = labeled_quark

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

labeled_quark: main.o ktruss.o k3-4.o hierarchy.o util.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o ktruss.o k3-4.o hierarchy.o util.o bucket.o graph.o  $(LIBS)
	/bin/rm *.o
	
clean: 
	/bin/rm -f ../$(TARGET) *.o 
