CC=nvcc   
CFLAGS=-O3 # -g
CUSOURCES=AngleData.cu BondData.cu DihedralData.cu gpu_settings.cu NeighborListBinnedGPU.cu NeighborList.cu NeighborListNsqGPU.cu ParticleData.cu
CUOBJECTS=$(CUSOURCES:%.cu=%.cu.o)
CXX=g++
CXXFLAGS=-I/usr/local/cuda/include -DENABLE_STATIC -DSINGLE_PRECISION -DCUDA_ARCH=10 -DENABLE_CUDA -DLARGE_EXCLUSION_LIST -O3 # -g
# CUDA_ARCH = 10, 12, or 13 depending on compute capability of hardware (one in uqmd12 is 1.0 so CUDA_ARCH=10 is used)
LDFLAGS=-lcudart -L/usr/local/cuda/lib64 -lboost_thread-mt -lboost_filesystem-mt -lboost_signals-mt -lpthread 
SOURCES=AngleData.cc BinnedNeighborList.cc BinnedNeighborListGPU.cc BondData.cc ClockSource.cc Compute.cc DihedralData.cc ExecutionConfiguration.cc GPUWorker.cc Initializers.cc NeighborList.cc NeighborListNsqGPU.cc ParticleData.cc Profiler.cc WallData.cc
OBJECTS=$(SOURCES:.cc=.o)

all: libhoomd.so

.cc.o:
	$(CXX) $(CXXFLAGS) -fPIC $< -c

%.cu.o: %.cu
	$(CC) $(CFLAGS) --compiler-options '-fPIC' -c -o $@ $<

#static
libhoomd.a: $(OBJECTS) $(CUOBJECTS)
	ar -cvq $@ $(OBJECTS) $(CUOBJECTS) 
#shared
libhoomd.so: $(OBJECTS) $(CUOBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -shared -Wl -o $@ $(OBJECTS) $(CUOBJECTS) 

clean:
	rm -f $(OBJECTS) $(CUOBJECTS) libhoomd.so
