# Documentation for the new GROMOS GPU implementation

# Introduction

## namespace gpu
Everything we do on gpu is in namespace gpu. Cuda code is in separate subdirectory for a possiblity of other gpu accelerators (OpenCL,...).

## Classes:

- CudaManager - Top-level interface to all GPU devices and streams (maybe rename to CudaInterface).
- CudaDeviceManager - Per-device metadata, streams, and hooks to memory/worker
- CudaMemoryManager - Manages allocations and sync on a per-device basis
- CudaDeviceWorker - Launches kernels, handles async tasks, manages streams

## Data structure classes:
- cuvector<T> - CUDA managed unified memory version of std::vector, directly accessible from both CPU and GPU. Possible penalty due to non-ideal automatic copy scheduling.
- cuhvector<T> - Pinned version of std::vector, fast access from GPU for async manual copies.
- container<T> - 2D jagged array, a light-weight device-only version of std::vector< std::vector >

- possibly a pair of <device cuhvector,cuarray> encapsulated in a single class could be used for performance-critical arrays


## Constant Memory
CUDA provides a fast, read-only constant memory space for storing constants. However, on modern GPUs, global memory caching often makes this unnecessary.
If you feel we still need it, possible implementation is in `memory/constants.h`.

## Coalesced Reads/Writes
CUDA performs best when adjacent threads access adjacent memory locations. Historically, `float3` required padding due to its 12-byte size versus the 16-byte memory access granularity. On GPUs from the last decade, this is no longer a concern — `float3` is now about 25% faster than float4 in many cases.

Packing additional data (like atomic charges) to "fill" `float4` is no longer beneficial. Similarly, reordering data from an AoS layout (e.g., x1, y1, z1...) to SoA (x1, x2, ..., y1, y2, ...) is generally unnecessary; modern CUDA compilers and schedulers optimize access patterns automatically.

## Floating Point Precision
CUDA defaults to single-precision (32-bit), while CPUs often use double-precision (64-bit). CUDA supports doubles, but:
- They consume more memory and registers,
- Run slower, and
- Can reduce occupancy.

### Mixed precision
A practical approach is to compute bulk quantities (e.g., forces, energies) in single precision, and accumulate results in double precision. This maintains accuracy while improving performance.
In some cases, half precision (16-bit) can be used for non-critical values.
Pure single precision is simple and fast but performs poorly in molecular dynamics (MD) due to accumulated numerical errors.

### Quantization
Quantization involves converting floats to integers, leveraging the fact that integer operations are 2–16× faster on CUDA.
It offers better performance than mixed precision but is harder to implement correctly, requiring precise control over value ranges and scaling. Best suited when data ranges are known and stable.

### Work scheduling
CUDA is able to overlap work with CPU. Traditional linear Algorithm_Sequence is therefore not optimal.
In the beginning of the algorithm sequence, we create CudaManager.
Each algorithm, that needs CUDA will submit the tasks to this CudaManager.
Additionally, we need to know data dependencies to overlap work efficiently.
If alg1 on cpu updates positions, and only alg4 needs them,alg1 have to start copying the positions to GPU already, and only then passes
control to alg2.
The dependency map (Directed Acyclic Graph) can be created at the start.

#### Algorithm class modifications
We can give these members to Algorithm class:
```cpp
virtual bool uses_cuda() const { return false; }

virtual int enqueue_cuda_work(topology::Topology &topo,
                              configuration::Configuration &conf,
                              simulation::Simulation &sim,
                              CudaManager &cuda_mgr) {
  return 0; // Default: not using CUDA
}
```

And data dependencies:
```cpp
enum class ArrayID { Positions, Velocities, Forces, Count };

struct DataRequirement {
  ArrayID id;
  enum class Access { Read, Write } access;
};

virtual std::vector<DataRequirement> get_data_requirements() const {
  return {}; // default: no special dependencies
}
```


## CudaManager
Is the first and the only contact point to the GPU(s). It should wrap all CUDA functionalities into nicely readable and easily usable functions. This includes:
- device control - the cuda manager probes all available GPU resources and sets itself accordingly. Then it exposes all functionalities to the GROMOS code.
- memory access - variable creation should be transparent, users should not operate with pointers and have direct access to the variables. The copying and everything is wrapped and hidden in the background.
- work submission - accepting jobs from the gromos code. Gromos only provides data, dependencies and the function to be called.
- scheduling - based on the data dependencies, submits the functions, kernels and provides input/output
- load balancing - balances the workload across multiple GPUS

There are multiple ways on how to do this. Ideally, the user has to do as least C++ coding as possible and use the wrappers.
E.g. variables can be defined dynamically and map<string,VAR>.

We have to think how to do the scheduling.
For every algorithm, one should define what are the output and input dependencies. Based on that, the jobs are ideally run asynchronously on GPU. We can e.g. create an execution plan based on the algorithm sequence and use directed acyclic graph (DAG).


Possible declaration of string-named arrays:
```cpp
class CudaManager {
public:
    void async_copy_to_device(const std::string& array_name);
    void async_copy_to_host(const std::string& array_name);
    bool is_on_device(const std::string& array_name);
    bool is_on_host(const std::string& array_name);
    void enqueue_kernel(Algorithm* alg);
    void synchronize_array(const std::string& array_name); // only if needed
};
```
So we access arrays by name strings. This allows easier creation of new, custom, user arrays. Alternatively, we can use enums:
```cpp
enum class ArrayID {
    Positions,
    Velocities,
    Forces,
    Temperature,
    Count
};

// Optionally map to/from string for UI/debugging
static const std::unordered_map<std::string, ArrayID> name_to_id = {
    {"positions", ArrayID::Positions},
    {"velocities", ArrayID::Velocities},
    {"forces", ArrayID::Forces},
};
```

During the md sequence preparation, we need to decide, if GPU variants are used. Also, mutual data dependencies tracking
allows us to schedule overlap in CPU and GPU execution.


Then, we add to md_seq Algorithm with appropriate backend:
```cpp
if (sim.cuda_enabled()) {
    // Use GPU backend
    md_seq.emplace_back(std::make_unique<M_Shake<gpuBackend>>();) // we might also want to switch to smart pointers
} else {
    // Use CPU backend
    md_seq.emplace_back(std::make_unique<M_Shake<cpuBackend>>();)
}
```

`M_Shake<GPU>` then calls cuda_shake() declared in gpu/cuda/.../shake.h:

```cpp
// shake.h
namespace gpu {
    void shake();
}
```

```cpp
// shake.cu
#include "cuda_shake.h"
#include <iostream>

void gpu::shake() {
    std::cout << "Running CUDA shake implementation\n";
    // ...
}
```

```cpp
// shake.cc
void gpu::shake() {
    DISABLED_VOID(); // macro of function doing nothing
}
```
## Code splitting for CPU and GPU
Conditional compilation for CPU-only or CPU+GPU requires careful split of the code. The source code for CUDA
have to be invisible in the CPU-only compilation. Traditionally, this is achieved by using the `#ifdef USE_CUDA`
macro, which is very easy to use, but with tighter and broader integration of CUDA-enabled code, the created clutter makes the code hard to read. So we try another strategy.
The single communication point, the `CudaManager` is split into `.cc` and `.cu` compilation. The `.cc` source
is for CPU-only compilation and contains empty method definitions using the `DISABLED()` macro. This macro
is supposed to throw a runtime error any time you try to interact with the `CudaManager` in the CPU-only
compilation. The `.cu` source is compiled in CUDA-enabled build and contains full working method definitions.

The CUDA-capable `Algorithm`s and `Interaction`s are defined as templates with common bases to achieve clean code separation.
Based on the compile time constants and/or runtime settings, we then dispatch either `MyAlgorithm<cpuBackend>`
or `MyAlgorithm<gpuBackend>`. `MyAlgoritm` defaults to `MyAlgorithm<cpuBackend>`.
The common bases are `IAlgorithm`, `IInteraction` and `IForcefield` allowing sequencing their base pointers
in the execution sequence.

## Algoritms can be CPU or GPU
The `Algorithm` class has been reworked such that it is a template `AlgorithmT`. CPU and GPU variants are initialized:
```cpp
// CPU variant
AlgorithmT<util::cpuBackend> alg(os);

// GPU variant
AlgorithmT<util::gpuBackend> alg(os);

// also legacy invocation is possible
Algorithm alg(os); // invokes AlgorithmT<util::cpuBackend> and works as before.
```

The common interface is `IAlgorithm`, so `Algorithm_Sequence` becomes `std::vector<std::unique_ptr<IAlgorithm>>`.

Every algoritm requires a single template declaration:

```cpp
// shake.h
  template<typename Backend>
  class Shake : public AlgorithmT<Backend>
  {
  public:
    /**
     * Constructor is inherited
     */
    using AlgorithmT<Backend>::AlgorithmT;

    /**
     * Destructor.
     */
    virtual ~Shake() {}
```

Then in the source
```cpp
//shake_cpu.cc
Shake<util::cpuBackend>::apply(...) {
    ...
}
//shake_gpu.cc
Shake<util::gpuBackend>::apply(...) {
    ...
}
```

The reason for this is to have some compile-time control. Every algorithm now has to be explicitly defined as CPU or GPU variant.
In a CPU-only compilation, the GPU variants should not get compiled at all, and an error is thrown,
if you try to instantiate them. This is achieved by the `is_valid_algorithm_backend` type trait:
```cpp
// util/backend.h
template <typename Backend>
struct is_valid_algorithm_backend
{
    static constexpr bool value =
    std::is_same_v<Backend, cpuBackend> ||
#ifdef USE_CUDA
      std::is_same_v<Backend, gpuBackend>;
#else
    false;
#endif
};

// algorithm.h
  class AlgorithmT : public IAlgorithm {
    /**
     * @brief compile time check for CUDA enabled
     * 
     */
    static_assert(util::is_valid_algorithm_backend<Backend>::value,
                "This backend is not supported in the current build configuration.");
  }
```
Instantiating the GPU variant is not allowed in CPU-only compilation:
```cpp
// create_md_sequence.cc
if (sim.cuda_enabled()) {
    md_seq.push_back(Shake<util::gpuBackend>()); // this does not compile in USE_CUDA=OFF
} else {
    md_seq.push_back(Shake<util::cpuBackend>());
}
```

Instead, you should use the factory function, that checks the availability of the GPU variant at the compile time, and compiles the 
appropriate factory function:
```cpp
  template <template <typename> class AlgT, typename... Args>
  IAlgorithm* make_algorithm(
                            simulation::Simulation & sim, 
                            Args&&... args) {
    if constexpr (util::has_gpu_backend_v<AlgT>) {
      if (sim.cuda_enabled()) {
        return new AlgT<util::gpuBackend>(std::forward<Args>(args)...);
      }
    }
    return new AlgT<util::cpuBackend>(std::forward<Args>(args)...);
  }
```
Usage example:
```cpp
algorithm::make_algorithm<algorithm::Remove_COM_Motion>(sim, os);
```
The provided arguments are passed to the constructor of `Remove_COM_Motion`.

If somehow a `Algorithm<util::gpuBackend>` is created somewhere by accident, the CPU-only build should not compile any
`AlgorithmT<util::gpuBackend>`.


## How to implement GPU variant of an algorithm?
1. Change the class declaration into a template:
```cpp
// in my_algorithm.h
// Change this into
  class MyAlgorithm : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    MyAlgorithm() : Algorithm("MyAlgorithm") {}
  }
// this
  template <typename Backend = util::cpuBackend>
  class MyAlgorithm : public AlgorithmT<Backend>
  {
  public:
    /**
     * Constructor.
     */
    MyAlgorithm() : AlgorithmT<Backend>("MyAlgorithm") {}
  }
```
2. Change all definitions into template definitions and define the variants implementations.
```cpp
// from this
int algorithm::MyAlgorithm::init()

//into this
// preferably in my_algorithm_cpu.cc
template<>
int algorithm::Remove_COM_Motion<util::cpuBackend>::init() {
    ...
}
// preferably in my_algorithm_gpu.cc
template<>
int algorithm::Remove_COM_Motion<util::gpuBackend>::init() {
    ...
}

// or this if the function is the same or similar for both backends
template<typename Backend>
int algorithm::Remove_COM_Motion<Backend>::init() {
    ... // Use Backend in logic, where you need conditionals
}
```

3. You can then either create the `AlgoritmT<util::cpuBackend>` or `AlgoritmT<util::gpuBackend>` directly, or use:
```cpp
algorithm::make_algorithm<algorithm::MyAlgorithm>(sim/*, optional arguments */);

```
to let the program create the variant based on user input and availability.

4. Update the file lists accordingly in `Makefile.am` and `CMakeLists.txt`
5. Check compilations with and without CUDA.

For examples, see commit `#6d3bb58f8cb037b256e9a7cd70870ba6d653da73`
The template of `Lattice_Shift_Tracker` there specifies indentical CPU and GPU variants.

For `Remove_COM_Motion`, variants are completely separate.

6. Sometimes you might need to change m_timer to this->m_timer. Also, if you do not provide explicit specializations, 
an explicit instantiation at the bottom of the `.cc` file is neccessary for the linker:
```cpp
template class algorithm::Lattice_Shift_Tracker<util::cpuBackend>;
```


## Cuda Manager conditional compilation

There is also a soft restriction on CudaManager. In CPU-only build, all its methods are stubs throwing a runtime error.
If you anyhow use them somewhere, they compile, but throw a runtime error as soon as you try to use them.
Turning them into a hard restriction (compile time) should be considered.

## Lattice_Shift_Tracker

The `Lattice_Shift_Tracker` is a class that tracks the shift of a lattice in a simulation. It has both CPU and GPU variants, which are specified using template specialization. The CPU variant is defined as follows:
