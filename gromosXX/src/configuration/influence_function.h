/**
 * @file influence_function.h
 * the optimal influence function class
 */

#ifndef INCLUDED_INFLUENCE_FUNCTION_H
#define	INCLUDED_INFLUENCE_FUNCTION_H

namespace configuration {
  /**
   * @class Influence_Function
   * the optimal influence function and logic
   */
  class Influence_Function {
  public:
    /**
     * Constructor
     */
    Influence_Function();
    /**
     * initialize the arrays
     */
    void init(const simulation::Parameter & param);
    /**
     * calculates the optimal influence function for a certain
     * configuration and stores it in the configuration.
     * The quality constant Q is also calculated.
     *
     * The optimal influence function is calculated as
     * @f[  \hat{G}_g^0(\mathbf{k}_{\mathbf{l}}) =
             \frac{\hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}}) \cdot
             [\sum_{\mathbf{m} \in Z^3} \mathbf{k}_{\mathbf{l},\mathbf{m}}
             k_{\mathbf{l},\mathbf{m}}^{-2}
             \hat{\gamma}(ak_{\mathbf{l},\mathbf{m}})\hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})]}
             {\mathbf{\hat{D}}_g^2(\mathbf{k}_{\mathbf{l}})
             [\sum_{\mathbf{m} \in Z^3} 
             \hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})]^2} @f]
     *
     * In practice, summation is restricted to @f$ \mathbf{m} @f$ vectors with
     * integer components @f$ [-m_{max}...m_{max}] @f$.
     *
     * If required (virial calculation), it also calculates the derivative
     * by evaluation of
     * @f[ \underline{\hat{\Gamma}}_g^0(\mathbf{k}_{\mathbf{l}}) =
             \frac{1}{\mathbf{\hat{D}}_g^2(\mathbf{k}_{\mathbf{l}})
             [\sum_{\mathbf{m} \in Z^3} 
             \hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})]^2}
             \sum_{\mathbf{m} \in Z^3} 
             \frac{\hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})}{k_{\mathbf{l},\mathbf{m}}^2}
             \left[ \left[ \mathbf{k}_{\mathbf{l},\mathbf{m}} \otimes \hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}})
             + \hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}}) \otimes \mathbf{k}_{\mathbf{l},\mathbf{m}}
             - 2 \mathbf{k}_{\mathbf{l},\mathbf{m}} \cdot \hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}}) \left[
             \frac{\mathbf{k}_{\mathbf{l},\mathbf{m}} \otimes \mathbf{k}_{\mathbf{l},\mathbf{m}}}{k_{\mathbf{l},\mathbf{m}}^2} + 
             \frac{\hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}}) \otimes \hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}})}{\hat{D}^2_g(\mathbf{k}_{\mathbf{l}})}
             \right]\right] \hat{\gamma}(ak_{\mathbf{l},\mathbf{m}})
             + \mathbf{k}_{\mathbf{l},\mathbf{m}} \cdot \hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}})
             \frac{\mathbf{k}_{\mathbf{l},\mathbf{m}} \otimes \mathbf{k}_{\mathbf{l},\mathbf{m}}}{k_{\mathbf{l},\mathbf{m}}^2}
             ak_{\mathbf{l},\mathbf{m}}\hat{\gamma}'(ak_{\mathbf{l},\mathbf{m}}) \right] @f]
     *
     * The quality of the newly calculated optimal influence function can
     * by estimated from the value @f$ Q @f$ which is calculated together with
     * the influence function by evaluation of
     *
     * @f[ Q = 16\pi^2V^{-1} \sum_{\mathbf{l} \in G, l \ne 0} \left[
             \sum_{\mathbf{m} \in Z^3} k_{\mathbf{l},\mathbf{m}}^{-2}
             \hat{\gamma}^2(ak_{\mathbf{l},\mathbf{m}})
             - \frac{\left[\hat{\mathbf{D}}_g(\mathbf{k}_{\mathbf{l}}) \cdot
             [\sum_{\mathbf{m} \in Z^3} \mathbf{k}_{\mathbf{l},\mathbf{m}}
             k_{\mathbf{l},\mathbf{m}}^{-2}
             \hat{\gamma}(ak_{\mathbf{l},\mathbf{m}})\hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})]\right]^2}
             {\mathbf{\hat{D}}_g^2(\mathbf{k}_{\mathbf{l}})
             [\sum_{\mathbf{m} \in Z^3} 
             \hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})]^2}
             \right] @f]
     *
     * Finally the estimated RMS force error is calculated by
     * @f[ \Delta F = \frac{\tilde{S}^2}{4\pi\epsilon_0} 
              \sqrt{\frac{Q}{N_q V}} @f]
     */
    template<class MeshType>
    void calculate(
      const topology::Topology & topo,
      configuration::Configuration & conf,
      const simulation::Simulation & sim);

    /**
     * accessor to the initial influence function
     */
    inline double & ghat_0(unsigned int x, unsigned int y, unsigned int z) {
      return ghat(x, y, z);
    }
    /**
     * const accessor to the initial influence function
     */
    inline const double & ghat_0(unsigned int x, unsigned int y, unsigned int z) const {
      return ghat(x, y, z);
    }
    /**
     * accessor to the initial influence function derivative
     */
    inline math::SymmetricMatrix & gammahat_0(unsigned int x, unsigned int y, unsigned int z) {
      return gammahat(x, y, z);
    }
    /**
     * const accessor to the initial influence function derivative
     */
    inline const math::SymmetricMatrix & gammahat_0(unsigned int x, unsigned int y, unsigned int z) const {
      return gammahat(x, y, z);
    }
    /**
     * accessor to the influence function
     *
     * if the box is scaled to influence function is corrected:
     * @f[\hat{G}_g^\dagger(\mathbf{k}_{\mathbf{l}}) =
             \hat{G}_g^0(\mathbf{k}_{\mathbf{l}}) - \mathrm{Tr}\left[
             \underline{\hat{\Gamma}}_g^0(\mathbf{k}_{\mathbf{l}})
             ({}^t\underline{L}_0)^{-1}({}^t\underline{L} - {}^t\underline{L}_0)
             \right] @f]
     */
    inline double operator()(const unsigned int & x,
            const unsigned int & y,
            const unsigned int & z) const {
      double result = ghat(x,y,z);
      if (do_scale) {
        math::Matrix gamma(gammahat(x,y,z));
        result -= math::trace(math::product(math::product(gamma, tLi_0),tbox_dev));
      }
      return result;
    }
    /**
     * accessor to the influence function, takes peridocity into account.
     */
    inline double operator()(const math::GenericVec<int> & p) const {
      return (*this)((p(0) + ghat.x()) % ghat.x(),
                (p(1) + ghat.y()) % ghat.y(),
                (p(2) + ghat.z()) % ghat.z());
    }
    /**
     * set the current box
     */
    void setBox(const math::Box & box);
    /**
     * accessor to the quality
     */
    inline double quality() const {
      return q;
    }
    /**
     * calculate the RMS force error from Q
     */
    inline double rms_force_error() const {
      return force_error;
    }
    /**
     * evaluate the quality the optimal influence function. 
     * The influence function is taken from the configuration and evaluated against
     * the current box. The quality constant Q is calculated and stored. 
     *
     * The value Q is calculated as
     * @f[ Q = 16\pi^2V^{-1} \sum_{\mathbf{l} \in G, l \ne 0} \left[
             \sum_{\mathbf{m} \in Z^3} k_{\mathbf{l},\mathbf{m}}^{-2}
             \hat{\gamma}^2(ak_{\mathbf{l},\mathbf{m}}) +
             [\hat{G}_g^\dagger]^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})
             \mathbf{\hat{D}}_g^2(\mathbf{k}_{\mathbf{l}}) \left[
             \sum_{\mathbf{m} \in Z^3}\hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})\right]^2
             -2 \hat{G}_g^\dagger(\mathbf{k}_{\mathbf{l},\mathbf{m}})
             \mathbf{\hat{D}}_g(\mathbf{k}_{\mathbf{l}}) 
             \cdot \sum_{\mathbf{m} \in Z^3}
             \mathbf{k}_{\mathbf{l},\mathbf{m}}
             k_{\mathbf{l},\mathbf{m}}^{-2}
             \hat{\gamma}(ak_{\mathbf{l},\mathbf{m}})\hat{P}^2(\mathbf{k}_{\mathbf{l},\mathbf{m}})
             \right] @f]
     *
     * The RMS force error is calculated as in described 
     * @ref configuration::Influence_Function::calculate "here"
     *
     * @sa configuration::Influence_Function::calculate
     */
    template<class MeshType>
    void evaluate_quality(const topology::Topology & topo,
      configuration::Configuration & conf,
      const simulation::Simulation & sim);
  protected:
    /**
     * the influence function generated on update
     * @f$ \hat{G}^{0}_{g} @f$
     */
    GenericMesh<double> ghat;
    /**
     * the influence function initial derivative by the initial box
     * @f$ \hat{\Gamma}^{0}_{g} @f$
     */
    GenericMesh<math::SymmetricMatrix> gammahat;
    /**
     * is the derivative calculated?
     */
    bool do_virial;
    /**
     * is the influence function corrected to the box
     */
    bool do_scale;
    /**
     * initial triclinic lattice @f$ ({}^{t}L^{0})^{-1} @f$
     */
    math::Matrix tLi_0;
    /**
     * transpose of initial box @f$ {}^{t}L^{0} @f$
     */
    math::Matrix tL_0;
    /**
     * transpose of current box
     */
    math::Matrix tL;
    /**
     * deviation from the initial box
     */
    math::Matrix tbox_dev;
    /**
     * the quality measure q
     */
    double q;
    /**
     * the RMS force error
     */
    double force_error;
    
  };
}

#endif	/* INCLUDED_INFLUENCE_FUNCTION_H */

