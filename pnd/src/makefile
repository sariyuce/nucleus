include makefile.in

TARGET = pnd 

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

all: $(TARGET) 

pnd: main.o local12.o local23.o local34.o bucket.o graph.o
	$(LD) $(LDFLAGS) -o ../$@ main.o local12.o local23.o local34.o bucket.o graph.o $(LIBS)
	/bin/rm *.o  

clean: 
	/bin/rm -f *.o ../pnd *~ $(TARGET) 

depend: $(SRCS)
	makedepend -Y -- $(INCLUDES) -- $(SRCS) 


#********************************************************************
# Dependencies
#********************************************************************
# DO NOT DELETE THIS LINE -- make depend depends on it.
