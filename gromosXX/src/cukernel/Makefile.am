AM_CXXFLAGS = @MY_CXXFLAGS@
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS =
SUFFIXES = .cu

.cu.o:
	$(NVCC) -arch sm_30 $(CFLAGS) $ -dc --ptxas-options=-v --compiler-options '-fPIC -O2 -D DNDEBUG -lcuda -lcudart' $< -o $@

##libgrocukernel.o:
	##$(NVCC) -arch sm_30 $(CFLAGS) $ --ptxas-options=-v --compiler-options '-fPIC -O2 -D DNDEBUG -lcuda -lcudart' -dlink lib/utils.o -o libgrocukernel.o

##libgrocukernel.a:
	##$(NVCC) --lib -o libgrocukernel.a cudaKernel.o lib/utils.o libgrocukernel.o

libcukernel.o: cudaKernel.o lib/utils.o
	$(NVCC) $(CFLAGS) $ --ptxas-options=-v --compiler-options '-fPIC -O2 -D DNDEBUG -lcuda -lcudart' -dlink cudaKernel.o lib/utils.o -o libcukernel.o

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src

xxincludedir = $(includedir)/md++/cukernel

xxinclude_HEADERS = cudaKernel.h\
	lib/utils.h

noinst_LIBRARIES = libgrocukernel.a

libgrocukernel_a_SOURCES = cudaKernel.cu cudaKernel.h\
	lib/utils.cu

libgrocukernel_a_LIBADD = libcukernel.o

libgrocukernel_a_CFLAGS = $(CUDA_CFLAGS) -I$(top_srcdir)/src
libgrocukernel_a_AR = $(NVCC) -arch=sm_30 --lib -o
##libgrocukernel_a_LINK = $(NVCC) -arch=sm_30 --lib -o
##$(NVCC) --lib -o libcukernel.a cudaKernel.o lib/utils.o libcukernel.o

bin_PROGRAMS =

check_PROGRAMS =

TESTS =