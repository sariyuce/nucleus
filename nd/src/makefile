include makefile.in

INCLUDES = -I/global/homes/a/aerdem/sparsehash-master/src # REPLACE THIS WITH YOUR PATH TO SPARSEHASH library

TARGET     = nucleus

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

nucleus: main.o kcore.o k1-3.o k1-4.o ktruss.o k2-4.o k3-4.o hierarchy.o util.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o kcore.o k1-3.o k1-4.o ktruss.o k2-4.o k3-4.o hierarchy.o util.o bucket.o graph.o  $(LIBS)
	/bin/rm *.o
	
clean: 
	/bin/rm -f ../$(TARGET) *.o 
