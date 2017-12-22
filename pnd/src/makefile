include makefile.in

INCLUDES = -I/user/erdem/sparsehash-master/src # REPLACE THIS WITH YOUR PATH TO SPARSEHASH library

TARGET = pnd 

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

all: $(TARGET) 

pnd: main.o local12.o hier.o util.o pitruss.o kcore.o ktruss.o k3-4.o pi34.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o local12.o hier.o util.o pitruss.o kcore.o ktruss.o k3-4.o pi34.o bucket.o graph.o $(LIBS)
	/bin/rm *.o  

clean: 
	/bin/rm -f *.o ../pnd *~ $(TARGET) 

depend: $(SRCS)
	makedepend -Y -- $(INCLUDES) -- $(SRCS) 


#********************************************************************
# Dependencies
#********************************************************************
# DO NOT DELETE THIS LINE -- make depend depends on it.
