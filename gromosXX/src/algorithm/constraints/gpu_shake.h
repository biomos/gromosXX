/**
 * @file gpu_shake.h
 * m_shake algroithm on the gpu
 */
#ifndef _GPU_SHAKE_H
#define	_GPU_SHAKE_H

#ifdef HAVE_LIBCUKERNEL
#include <cudaKernel.h>
#else
#define gpu_status void
#endif
#include <algorithm/constraints/gpu_shake_thread.h>

namespace interaction {
  struct bond_type_struct;
}

namespace algorithm
{
  /**
   * @class GPU_Shake
   * implements the m_shake algorithm.
   */
  class GPU_Shake : public Algorithm
  {
  public:
    /**
     * Constructor.
     */
    GPU_Shake(double const tolerance = 0.000001,
	  int const max_iterations = 1000,
	  std::string const name = "GPU_Shake");

    /**
     * Destructor.
     */
    virtual ~GPU_Shake();

    /**
     * apply gpu_shake.
     */
    virtual int apply(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim);
    /**
     * set the tolerance.
     */
    void tolerance(double const tol);

    /**
     * tolerance.
     */
    double const & tolerance()const {return m_tolerance;}
    /**
     * max iterations.
     */
    int const & max_iterations()const {return m_max_iterations;}

    /**
     * the const bond type parameter.
     */
    std::vector<interaction::bond_type_struct> const &parameter()const
    {
      return m_parameter;
    }
    /**
     * the bond type parameter.
     */
    std::vector<interaction::bond_type_struct> & parameter()
    {
      return m_parameter;
    }
    /**
     * initialize startup positions and velocities
     * if required.
     */
    virtual int init(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     std::ostream & os = std::cout,
		     bool quiet = false);

  protected:

    /**
     * shake tolerance
     */
    double m_tolerance;
    /**
     * max iterations
     */
    const int m_max_iterations;
    /**
     * bond parameter
     */
    std::vector<interaction::bond_type_struct> m_parameter;
    /**
     * rank and size for parallelization
     */
    int m_rank, m_size;

    void solvent
    (
     topology::Topology const & topo,
     configuration::Configuration & conf,
     simulation::Simulation & sim,
     double dt, int const max_iterations,
     int & error
     );

    /**
     * Pointer to all the GPU_Shake_Threads
     */
    std::vector<GPU_Shake_Thread *> m_shake_set;


  };

} //algorithm

#endif	/* _GPU_SHAKE_H */

