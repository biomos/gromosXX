
SUFFIXES = .cu
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src
NVCC_CFLAGS = @NVCC_CFLAGS@

## Tracking dependencies correctly as with CXX would need a macro for NVCC
## This would be too tedious and I hope that CMake will be handling this better
## For now, we just delete the intermediate files to always enforce recompilation
## to avoid confusion after changes or reconfiguration
.cu.o:
	$(MKDIR_P) $(@D)
	$(NVCC) $(NVCCFLAGS) -dc --compiler-options '$(NVCC_CFLAGS)' $< -o $@

libgrocukernel.o: $(libgrocukernel_a_OBJECTS)
	$(NVCC) $(NVCCFLAGS) --compiler-options '$(NVCC_CFLAGS)' -dlink $^ -o $@

libgrocukernel.a: $(libgrocukernel_a_OBJECTS) libgrocukernel.o
	$(NVCC) $(NVCCFLAGS) --lib -o $@ $^
	rm -f $^

noinst_LIBRARIES = libgrocukernel.a

libgrocukernel_a_SOURCES =  cudaKernel.cu\
							lib/utils.cu

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
