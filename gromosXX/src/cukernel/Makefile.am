
SUFFIXES = .cu
AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src
NVCC_CFLAGS = @NVCC_CFLAGS@

## Tracking dependencies correctly as with CXX would need a macro for NVCC
## This would be too tedious and I hope that CMake will be handling this better
## For now, we just delete the intermediate files to always enforce recompilation
## 
##	$(NVCC) $(NVCCFLAGS) -std=c++11 -dc --compiler-options '-std=c++11 -I../.. -I../../../src $(CXXFLAGS) $(NVCC_CFLAGS)' $< -o $@
##	$(NVCC) $(NVCCFLAGS) --compiler-options '-std=c++11 -I../.. -I../../../src $(CXXFLAGS) $(NVCC_CFLAGS)' -dlink $^ -o $@
.cu.o:
	$(NVCC) $(NVCCFLAGS) -std=c++11 -dc --compiler-options '$(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(CXXFLAGS) $(NVCC_CFLAGS)' $< -o $@

libgrocukernel.o: $(libgrocukernel_a_OBJECTS)
	$(NVCC) $(NVCCFLAGS) --compiler-options '$(CXXFLAGS) $(NVCC_CFLAGS)' -dlink $^ -o $@

libgrocukernel.a: $(libgrocukernel_a_OBJECTS) libgrocukernel.o
	$(NVCC) $(NVCCFLAGS) --lib -o $@ $^
	rm -f $^ # Remove this line if you dont want to recompile every time

if WITH_CUDA
noinst_LIBRARIES = libgrocukernel.a
libgrocukernel_a_SOURCES =  cuda_kernel.cu \
							io.cu \
							pairlist.cu \
							interaction.cu \
							constraints.cu \
							cudaShake.cu \
							lib/utils.cu
else
xxincludedir = $(includedir)/md++/cukernel
##xxinclude_HEADERS =	cuda_kernel.h
noinst_LTLIBRARIES = libgrocukernel.la
libgrocukernel_la_SOURCES = cuda_kernel.cc \
							cukernel.cc
libgrocukernel_la_CXXFLAGS = $(AM_CXXFLAGS)
endif

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
