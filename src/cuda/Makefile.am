AM_CXXFLAGS = @MY_CXXFLAGS@
NVCC = /usr/local/cuda/bin/nvcc
NVCCFLAGS = -arch sm_20
SUFFIXES = .cu

.cu.o:
	$(NVCC) -o $@ -c $< $(NVCCFLAGS)
.cu.lo:
	$(LIBTOOL) --tag=CC --mode=compile $(NVCC) -o $@ -c $< $(NVCCFLAGS)

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)

xxincludedir = $(includedir)/md++/cuda

xxinclude_HEADERS = gromos_cuda.h\
		    cuda_mpi.h\
                    c_pairlist.h


noinst_LTLIBRARIES = libgrocuda.la

libgrocuda_la_SOURCES = gromos_cuda.cc\
			cuda_mpi.cc\
			cuda_pairlist.cu

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
