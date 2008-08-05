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
     */
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
    inline math::Matrix & gammahat_0(unsigned int x, unsigned int y, unsigned int z) {
      return gammahat(x, y, z);
    }
    /**
     * const accessor to the initial influence function derivative
     */
    inline const math::Matrix & gammahat_0(unsigned int x, unsigned int y, unsigned int z) const {
      return gammahat(x, y, z);
    }
    /**
     * accessor to the influence function
     */
    inline double operator()(const unsigned int & x,
            const unsigned int & y,
            const unsigned int & z) const {
      double result = ghat(x,y,z);
      if (hasDerivative) {
        result -= math::trace(math::product(math::product(gammahat(x,y,z), tLi_0),
                (tL - tL_0)));
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
     * evaluate the quality the optimal influence function. 
     * The influence function is taken from the configuration and evaluated against
     * the current box. The quality constant Q is calculated and stored in the
     * configuration.
     */
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
    GenericMesh<math::Matrix> gammahat;
    /**
     * is the derivative calculated?
     */
    bool hasDerivative;
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
     * the quality measure q
     */
    double q;
    
  };
}

#endif	/* INCLUDED_INFLUENCE_FUNCTION_H */

