SUFFIXES = .cu

AM_CXXFLAGS = @MY_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(top_builddir)
NVCC_CFLAGS = @NVCC_CFLAGS@

## Tracking dependencies correctly as with CXX would need a macro for NVCC
## This would be too tedious and I hope that CMake will be handling this better
## For now, we just delete the intermediate files to always enforce recompilation
## 
##	$(NVCC) $(NVCCFLAGS) -std=c++11 -dc --compiler-options '-std=c++11 -I../.. -I../../../src $(CXXFLAGS) $(NVCC_CFLAGS)' $< -o $@
##	$(NVCC) $(NVCCFLAGS) --compiler-options '-std=c++11 -I../.. -I../../../src $(CXXFLAGS) $(NVCC_CFLAGS)' -dlink $^ -o $@
.cu.o:
	$(NVCC) $(NVCCFLAGS) -dc -Xcompiler '$(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $(CXXFLAGS) $(NVCC_CFLAGS)' $< -o $@

# CUDA object list
libgrocuda_a_SOURCES = \
    gpu.cc \
    cuda/utils.cu \
    cuda/manager/cuda_manager.cc \
    cuda/kernels/reduction.cu \
    cuda/memory/container.cu \
    cuda/manager/cuda_device_worker.cu \
    cuda/manager/cuda_device_manager.cu \
    cuda/manager/cuda_memory_manager.cu \
    cuda/algorithm/constraints.cu

# Manual dlink step if needed for .a library
libgrocuda.o: $(libgrocuda_a_OBJECTS)
	$(NVCC) $(NVCCFLAGS) -Xcompiler '$(CXXFLAGS) $(NVCC_CFLAGS)' -dlink $^ -o $@

libgrocuda.a: $(libgrocuda_a_OBJECTS) libgrocuda.o
	$(NVCC) $(NVCCFLAGS) -Xcompiler '$(CXXFLAGS) $(NVCC_CFLAGS)' --lib -o $@ $^
	rm -f $^  # optional: remove intermediates to force recompile

# Conditional library building
if WITH_CUDA
    noinst_LIBRARIES = libgrocuda.a
else
    noinst_LTLIBRARIES = libgrocuda.la

    libgrocuda_la_SOURCES = \
        gpu.cc \
        cuda/algorithm/constraints.cc \
        cuda/manager/cuda_manager.cc \
        cuda/manager/cuda_device_worker.cc \
        cuda/manager/cuda_device_manager.cc \
        cuda/manager/cuda_memory_manager.cc

    libgrocuda_la_CPPFLAGS = $(AM_CPPFLAGS)
    libgrocuda_la_CXXFLAGS = $(AM_CXXFLAGS)
endif

bin_PROGRAMS =
check_PROGRAMS =
TESTS =
