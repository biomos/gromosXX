
SUFFIXES = .cu
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src
NVCC_CFLAGS = @NVCC_CFLAGS@

.cu.o:
	$(MKDIR_P) $(@D)
	$(NVCC) $(NVCCFLAGS) -dc --compiler-options '$(NVCC_CFLAGS)' $< -o $@

libgrocukernel.o: $(libgrocukernel_a_OBJECTS)
	$(NVCC) $(NVCCFLAGS) --compiler-options '$(NVCC_CFLAGS)' -dlink $^ -o $@

libgrocukernel.a: $(libgrocukernel_a_OBJECTS) libgrocukernel.o
	$(NVCC) $(NVCCFLAGS) --lib -o $@ $^

xxincludedir = $(includedir)/md++/cukernel
xxinclude_HEADERS = cudaKernel.h \
                    macros.h \
                    parameter.h \
                    pairlist.h \
                    interaction.h \
                    constraints.h \
                    cudaShake.h \
                    gpu_status.h \
                    io.cu \
                    pairlist.cu \
                    interaction.cu \
                    constraints.cu \
                    cudaShake.cu

xxculibdir = $(includedir)/md++/cukernel/lib
xxculib_HEADERS = lib/utils.h \
                  lib/math.h \
                  lib/types.h \
                  lib/float3.h \
                  lib/double3.h

noinst_LIBRARIES = libgrocukernel.a

libgrocukernel_a_SOURCES = cudaKernel.cu\
	lib/utils.cu

bin_PROGRAMS =

check_PROGRAMS =

TESTS =
