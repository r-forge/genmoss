TOPDIR=../..
include $(TOPDIR)/make.inc

#######################################################################
#  This is the makefile to create a wrapper for the CBLAS.
#######################################################################

all: libcblaswr.so

libcblaswr.so: cblaswr.o
	cc -shared -o libcblaswr.so cblaswr.o
	
libfblaswr.a: fblaswr.o
	$(ARCH) $(ARCHFLAGS) libfblaswr.a fblaswr.o
	$(RANLIB) libfblaswr.a

clean:
	rm -f *.o *.a *.so

.c.o: 
	$(CC) $(CFLAGS) -c $*.c

