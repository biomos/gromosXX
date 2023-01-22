/**
 * @file nonbonded_outerloop.h
 * the non bonded outerloops.
 */

#include "latticesum.h"


#ifndef INCLUDED_NONBONDED_OUTERLOOP_H
#define INCLUDED_NONBONDED_OUTERLOOP_H

namespace topology
{
  class Topology;
}
namespace configuration
{
  class Configuration;
}
namespace simulation
{
  class Simulation;
}
namespace util {
  class Algorithm_Timer;
}

namespace interaction
{
  class Nonbonded_Parameter;
  class Storage;
  class Pairlist;
  class Lattice_Sum;
  struct KSpace_Element;
  
  /**
   * @class Nonbonded_Outerloop
   * loops the nonbonded interactions...
   */
  class Nonbonded_Outerloop
  {
  public:    
    /**
     * Constructor.
     */
    Nonbonded_Outerloop(Nonbonded_Parameter & nbp);
    
    /**
     * calculate the lj crf interactions.
     */
    void lj_crf_outerloop(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  Pairlist const & pairlist_solute,
                          Pairlist const & pairlist_solvent,
			  Storage & storage,
                          bool longrange,
                          util::Algorithm_Timer & timer,
                          bool master);
    
    /**
     * calculate only lj interactions.
     */
    void lj_outerloop(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation & sim,
		      Pairlist const & pairlist_solute,
                      Pairlist const & pairlist_solvent,
		      Storage & storage,
                      bool longrange,
                      util::Algorithm_Timer & timer,
                      bool master);

    /**
     * helper function to calculate the forces and energies from the
     * 1,4 interactions.
     */
    void cg_exclusions_outerloop(topology::Topology & topo,
				 configuration::Configuration & conf,
				 simulation::Simulation & sim,
				 Storage & storage);

    /**
     * calculate the sasa interactions.
     */
    void sasa_outerloop(topology::Topology & topo,
                        configuration::Configuration & conf,
                        simulation::Simulation & sim,
                        Storage & storage);

    /**
     * calculate the 1,4-interactions.
     */
    void one_four_outerloop(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim,
			    Storage & storage,
                            int rank, int size);

    /**
     * calculate the LJ exceptions
     */
    void lj_exception_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage,
            int rank, int size);

    /**
     * calculate the RF contributions for excluded atoms.
     */
    void RF_excluded_outerloop(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       Storage & storage,
                               int rank, int size);
    
    /**
     * calculate the self energy (polarisation).
     */
    void self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    /**
     * calculate the electric field (polarisation)
     * and reposition the COS till field is consistent
     */
    void electric_field_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
		            PairlistContainer const & pairlist,
                            Storage & storage,
                            Storage & storage_lr,
                            int rank);
    
    /**
     * calculate the ls interactions. (lattice sum)
     * in real space
     */
    void ls_real_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Pairlist const & pairlist_solute,
            Pairlist const & pairlist_solvent,
            Storage & storage, int rank, int size);
    
    /**
     * calculate the ls interactions kspace interactions
     * using the Ewald method
     *
     * The energy is calculated as
     * @f[ E_\gamma = (2\epsilon_0V)^(-1)
             \sum_{\mathbf{l} \in W, l \ne 0} k^{-2}\hat{\gamma}(ak)\left[
             C^2(\mathbf{k}) + S^2(\mathbf{k})
             \right] @f]
     * with 
     * @f[ C(\mathbf{k}) = \sum_{i=1}^{N_q} q_i \cos(\mathbf{k} \cdot \mathbf{r}_i) @f]
     * @f[ S(\mathbf{k}) = \sum_{i=1}^{N_q} q_i \sin(\mathbf{k} \cdot \mathbf{r}_i) @f]
     *
     * The force acting on atom @f$ i @f$ is given by
     * @f[ \mathbf{F}_{\gamma, i} = (\epsilon_0V)^{-1} q_i
             \sum_{\mathbf{l} \in W, l \ne 0} k^{-2}\hat{\gamma}(ak)\mathbf{k}\left[
             C(\mathbf{k}) \sin(\mathbf{k} \cdot \mathbf{r}_i) -
             S(\mathbf{k}) \cos(\mathbf{k} \cdot \mathbf{r}_i)
             \right] @f]
     *
     * If required the reciprocal space virial is calculated as
     * @f[ \underline{W}_\gamma = -(4\epsilon_0V)^(-1) \sum_{\mathbf{l} \in W, l \ne 0}
             \left[ C^2(\mathbf{k}) + S^2(\mathbf{k}) \right]
             \left[\frac{ak\hat{\gamma}'(ak) - 2\hat{\gamma}(ak)}{k^4} \mathbf{k} \otimes \mathbf{k}
             + \mathrm{diag}(\frac{\hat{\gamma}(ak)}{k^2})\right] @f]
     *
     *
     * The methodology dependent A term @f$ \tilde{A}_2 @f$ is given by
     * @f[ \tilde{A}_2 = \frac{4\pi}{V} \sum_{\mathbf{l} \in W, l \ne 0} k^{-2}\hat{\gamma}(ak) @f]
     *
     * If required the derivative of the @f$ \tilde{A}_2 @f$ term is given by
     * @f[ \underline{\frac{d\tilde{A}_2}{dL}} = \frac{4\pi}{V}
               \sum_{\mathbf{l} \in W, l \ne 0} \left[
               \frac{(ak\hat{\gamma}'(ak) - 2\hat{\gamma}(ak)}{k^4} 
               \mathbf{k} \otimes \mathbf{k} + 
               \mathrm{diag}(\frac{\hat{\gamma}(ak)}{k^2}) \right]
     * @f]
     */
    void ls_ewald_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
    /**
     * calculate the ls interactions kspace interactions
     * using the P3M method
     */
    void ls_p3m_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size,
            util::Algorithm_Timer & timer,
            bool & is_ok);
    
    /**
     * calculate the ls self and constant interactions/energies
     *
     * The lattice sum self term @f$ \Delta G_{\mathrm{\it slf}} @f$ and constant
     * term @f$ E_A @f$ are calculated by the following equations:
     * 
     * @f[\Delta G_{\mathrm{\it slf}} = (8\pi\epsilon_o)^{-1}\,(A_1+A_2+A_3)\,\tilde{S}^2@f]
     * @f[E_A = (8\pi\epsilon_o)^{-1}\,[A_1\,S^2 - (A_1+\tilde{A}_2)\,\tilde{S}^2]@f]
     *
     * where @f$ A_1 @f$, @f$ A_2 @f$ and @f$ A_3 @f$ are quantities dependent on the
     * charge shape. @f$ A_1 @f$ and @f$ A_3 @f$ are calculated from tabulated values:
     *
     * <table border="0">
     * <tr><td>@f$ N_{\gamma} @f$</td><td>@f$ -V \pi^{-1} a^{-2} A_1 @f$</td><td>@f$ -a A_3 @f$</td></tr>
     * <tr><td>-1</td><td>@f$ 1 @f$</td><td>@f$ 2 \pi^{-1/2} @f$</td></tr>
     * <tr><td>0</td><td>@f$ 2/5 @f$</td><td>@f$ 3/2 @f$</td></tr>
     * <tr><td>1</td><td>@f$ 4/15 @f$</td><td>@f$ 2 @f$</td></tr>
     * <tr><td>2</td><td>@f$ 4/21 @f$</td><td>@f$ 5/2 @f$</td></tr>
     * <tr><td>3</td><td>@f$ 3/14 @f$</td><td>@f$ 9/4 @f$</td></tr>
     * <tr><td>4</td><td>@f$ 1/6 @f$</td><td>@f$ 21/8 @f$</td></tr>
     * <tr><td>5</td><td>@f$ 2/15 @f$</td><td>@f$ 3 @f$</td></tr>
     * <tr><td>6</td><td>@f$ 8/55 @f$</td><td>@f$ 45/16 @f$</td></tr>
     * <tr><td>7</td><td>@f$ 4/33 @f$</td><td>@f$ 25/8 @f$</td></tr>
     * <tr><td>8</td><td>@f$ 4/39 @f$</td><td>@f$ 55/16 @f$</td></tr>
     * <tr><td>9</td><td>@f$ 10/91 @f$</td><td>@f$ 105/32 @f$</td></tr>
     * <tr><td>10</td><td>@f$ 2/21 @f$</td><td>@f$ 455/128 @f$</td></tr>
     * </table>
     * where @f$ V @f$ is the box volume and @f$ a @f$ is the charge shape width.
     *
     * The @f$A2@f$ is calculated by the following scheme: There are several
     * options depending on the user input parameters
     * - it is set to 0.0
     * - it is set to @f$\tilde{A}_2@f$
     * - it is calculated by direct summation of k vectors.
     *
     * @f[A_2 = 4\pi\,V^{-1}\,\sum_{\mathbf{l} \in Z^3\,,\,l\ne 0}\,
               k^{-2}\,  \hat{\gamma}(ak)@f]
     *
     * In practice this sum is evaluated over increasing volumes till the user
     * specified precision is reached. For cubic periodic boundary conditions 
     * with box edge length @f$ b @f$ the value of @f$ A_2 @f$ is calculated as
     *
     * @f[A_2 = -2.83729748 / b - A_1 - A_3@f]
     *
     * The @f$ \tilde{A}_2 @f$ is methodology depedent and calculated in the 
     * Ewald and P3M routines.
     *
     * These terms have no force contribution.
     *
     * The virial contribution of these terms is calculated as:
     *
     * @f[\underline{W}_{\mathrm{\it slf}} = \underline{W}^{A_1}_{\mathrm{\it slf}} + \underline{W}^{A_2} @f]
     * @f[\underline{W}_A = \underline{W}^{A_1}_A - \underline{W}^{\tilde{A}_2} @f]
     *
     * where the the virials of the @f$ A_1 @f$ terms is calculated as 
     *
     * @f[\underline{W}^{A_1}_{\mathrm{\it slf}} = -\frac{1}{16\pi\epsilon_0}\mathrm{diag}(A_1)\tilde{S}^2@f]
     * @f[\underline{W}^{A_1}_A = -\frac{1}{16\pi\epsilon_0}\mathrm{diag}(A_1)(S^2-\tilde{S}^2)@f]
     *
     * and the virial from the @f$ A_2 @f$ term is calculated as
     *
     * @f[ \underline{W}^{A_2} = -\frac{\tilde{S}^2}{4V\epsilon_0}
               \sum_{\mathbf{l} \in Z^3\,,\,l\ne 0}\,\left[
               \frac{(ak\hat{\gamma}'(ak) - 2\hat{\gamma}(ak)}{k^4} 
               \mathbf{k} \otimes \mathbf{k} + 
               \mathrm{diag}(\frac{\hat{\gamma}(ak)}{k^2}) \right] @f]
     *
     * The virial from the @f$ \tilde{A}_2 @f$ term is valculated as
     * @f[ \underline{W}^{\tilde{A}_2} = - \frac{\tilde{S}^2}{16\pi\epsilon_0}
               \underline{\frac{d\tilde{A}_2}{dL}} @f]
     * The calculation of the box derivative @f$ \underline{\frac{d\tilde{A}_2}{dL}} @f$
     * is carried out in the Ewald and P3M routines.
     */
    void ls_self_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);

    /** 
     * Calculate the ls surface force, energy and virial.
     * 
     * The lattice sum surface free energy @f$ \Delta G_{\mathrm{\it srf}} @f$
     * is calculated using the following equation:
     * @f[ \Delta G_{\mathrm{\it srf}} =
              \frac{1}{2\epsilon_0(2\epsilon_{LS}+1)V} \mathbf{M}^2 @f]
     * where @f$ \mathbf(M) @f$ is the box dipole moment
     * @f[ \mathbf{M} = \sum_{i=1}^{N_q}\,q_i\,(\mathbf{r}_i - \mathbf{r}_c) @f]
     * and @f$ \mathbf{r}_c @f$ is the box centre.
     *
     * The force is calculated as
     * @f[ \mathbf{F}_{\mathrm{\it srf},i} =
              -\frac{1}{\epsilon_0(2\epsilon_{LS}+1)V} q_i \mathbf{M} @f]
     *
     * The virial is calculated as
     * @f[ \underline{W}_{\mathrm{\it srf}} = 
              \frac{1}{2\epsilon_0(2\epsilon_{LS}+1)V}\left[
               \frac{1}{2}\mathrm{diag}(M^2_x,M^2_y,M^2_z)
               - \frac{1}{4}\mathrm{diag}(\mathbf{M}^2)\right]@f]
     */
    void ls_surface_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
    /**
     * calculate the interaction for a given atom pair.
     * SLOW! as it has to create the periodicity...
     */
    int calculate_interaction(topology::Topology & topo,
			      configuration::Configuration & conf,
			      simulation::Simulation & sim,
			      unsigned int atom_i, unsigned int atom_j,
			      math::Vec & force, 
			      double &e_lj, double &e_crf);

    /**
     * calculate the hessian for a given pair
     */
    int calculate_hessian(topology::Topology & topo,
			  configuration::Configuration & conf,
			  simulation::Simulation & sim,
			  unsigned int atom_i, unsigned int atom_j,
			  math::Matrix & hessian,
			  PairlistContainer const & pairlist);

  protected:
    /**
     * the nonbonded parameter.
     */
    Nonbonded_Parameter & m_param;

#ifdef OMP
    /**
     * OMP shared total electric field
     */
    static math::VArray electric_field;
    /**
     * OMP shared convergence criterion
     */
    static double minfield;
#endif

    template<typename t_interaction_spec>
    void _lj_crf_outerloop(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Pairlist const & pairlist_solute,
                           Pairlist const & pairlist_solvent,
			   Storage & storage,
                           bool longrange, util::Algorithm_Timer & timer,
                           bool master);

    /**
     * In principle, this is the specialization with
     *    interaction::Interaction_Spec<math::rectangular, simulation::lj_crf_func>
     * Unfortunately, I am not able to properly include it as a specialization.
     * Might be a g++ bug ( https://gcc.gnu.org/bugzilla/show_bug.cgi?id=56480 ) or
     * an error in the implementation. For now, this is a workaround...
     */
    template<simulation::charge_type_enum t_charge_type>
    void _lj_crf_outerloop_fast(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   Pairlist const & pairlist_solute,
                           Pairlist const & pairlist_solvent,
			   Storage & storage,
                           bool longrange, util::Algorithm_Timer & timer,
                           bool master);

    template<typename t_interaction_spec>
    void _lj_outerloop(topology::Topology & topo,
                  configuration::Configuration & conf,
                  simulation::Simulation & sim,
                  Pairlist const & pairlist_solute,
                  Pairlist const & pairlist_solvent,
                  Storage & storage,
                  bool longrange, util::Algorithm_Timer & timer, bool master);

    template<typename t_interaction_spec>
    void _one_four_outerloop(topology::Topology & topo,
			     configuration::Configuration & conf,
			     simulation::Simulation & sim,
			     Storage & storage,
                             int rank, int size);

    template<typename t_interaction_spec>
    void _lj_exception_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage,
            int rank, int size);

    template<typename t_interaction_spec>
    void _cg_exclusions_outerloop(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
				  Storage & storage);

    template<typename t_interaction_spec>
    void _sasa_outerloop(topology::Topology & topo,
                         configuration::Configuration & conf,
                         simulation::Simulation & sim,
                         Storage & storage);

    template<typename t_interaction_spec>
    void _RF_excluded_outerloop(topology::Topology & topo,
				configuration::Configuration & conf,
				simulation::Simulation & sim,
				Storage & storage,
                                int rank, int size);
    
    template<typename t_interaction_spec>
    void _self_energy_outerloop(topology::Topology & topo,
		            configuration::Configuration & conf,
		            simulation::Simulation & sim, 
			    Storage & storage);

    template<typename t_interaction_spec>
    void _electric_field_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            PairlistContainer const & pairlist,
            Storage & storage,
            Storage & storage_lr,
            int rank);
    
    template<typename t_interaction_spec>
    void _ls_real_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Pairlist const & pairlist_solute,
            Pairlist const & pairlist_solvent,
            Storage & storage, int rank, int size);
    
    template<typename t_interaction_spec>
    void _ls_ewald_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);

    template<typename t_interaction_spec>
    void _ls_p3m_kspace_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size,
            util::Algorithm_Timer & timer,
            bool & is_ok);

    template<typename t_interaction_spec>
    void _ls_surface_outerloop(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim,
            Storage & storage, int rank, int size);
    
    template<typename t_interaction_spec>
    int _calculate_interaction(topology::Topology & topo,
			       configuration::Configuration & conf,
			       simulation::Simulation & sim,
			       unsigned int atom_i, unsigned int atom_j,
			       math::Vec & force,
			       double & e_lj, double & e_crf);

    template<typename t_interaction_spec>
    int _calculate_hessian(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   unsigned int atom_i, unsigned int atom_j,
			   math::Matrix & hessian,
			   PairlistContainer const & pairlist);
    
  };
  
} // interaction

#endif
