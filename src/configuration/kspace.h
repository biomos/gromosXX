/**
 * @file kspace.h
 * dealing with k space
 */

#ifndef INCLUDED_KSPACE_H
#define	INCLUDED_KSPACE_H

namespace configuration {
  /**
   * @class KSpace_Element
   * an element in k-space.
   * This class is many used for caching values that have to be calculated
   * multiple times in lattice sum methods
   */
  class KSpace_Element {
  public:
    /**
     * the k vector itself: @f$ \mathbf{k} @f$
     */
    math::Vec k;
    /**
     * inverted squared length of the k vector: @f$ k^{-2} @f$
     */
    double k2i;
    /**
     * k squared;
     */
    double k2;
    /**
     * absolute value of k
     */
    double abs_k;
    /**
     * absolute value of k aplified by charge width
     */
    double ak;
    /**
     * the fourier coefficient @f$ \hat{\gamma}(ak) @f$
     */
    double fourier_coefficient;
    /**
     * the fourier coefficient derivative @f$ \hat{\gamma}'(ak) @f$
     */
    double fourier_coefficient_derivative;
    /**
     * @f$ k^{-2}\hat{\gamma}(ak) @f$
     */
    double k2i_gammahat;
    /**
     * @f$ ak\hat{\gamma}'(ak) @f$
     */
    double ak_gammahat_prime;
  };
  
  /**
   * calculates the k space vectors and adds them to a list
   *
   * This functions loops over the three dimensions and creates @f$ \mathbf{l} @f$
   * vectors with @f$ l_{\mu} = -l_{\mu} ... +l_{\mu} @f$ where @f$ \mu @f$ denotes the
   * coordinate and @f$ l_{\mu} @f$ is given in the input file. From these @f$ \mathbf{l} @f$
   * vectors the @f$ \mathbf{k} @f$ vector is calculated and its values and fourier 
   * coefficient is stored if the absolute value of k is lower than a certain cutoff.
   *
   * If size is != 0, stride parallelization will be applied.
   */
  void calculate_k_space(
            const topology::Topology & topo,
            configuration::Configuration & conf,
            const simulation::Simulation & sim,
            int rank, int size);  
  
  /**
   * @class KSpace_Utils
   * static functions to deal with k space.
   */
  class KSpace_Utils {
  public: 
    /**
     * calculates the matrix that maps @f$\mathbf{l}@f$ vectors to k vectors via 
     * @f$\mathbf{k} = 2\pi {}^{t}\mathbf{\underline{L}}^{-1} \mathbf{l}@f$. 
     *
     * @return the matrix @f$2\pi {}^{t}\mathbf{\underline{L}}^{-1}@f$
     */
    static math::Matrix l_to_k_matrix(const math::Box & box, math::boundary_enum b);
  };
  
} // namespace configuration
#endif	/* INCLUDED_KSPACE_H */

