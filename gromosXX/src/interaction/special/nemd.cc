/**
 * @file nemd.cc
 * template NEMD_Interaction
 */


#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"
#include "../../interaction/interaction.h"

#include "../../math/periodicity.h"

// special interactions
#include "../../interaction/interaction_types.h"

#include "../../interaction/special/nemd.h"

#include "../../util/template_split.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special

/**
 * calculate NEMD
 */
int interaction::NEMD_Interaction::
calculate_interactions(topology::Topology& topo,
                       configuration::Configuration& conf,
                       simulation::Simulation& sim) {

  // Perturbation is applied every PERTFRQ timesteps
  // Commented out: now it is performed inside each method
  // if((sim.steps() % sim.param().nemd.pertfrq) == 0) {
    math::Periodicity<math::rectangular> periodicity(conf.current().box);

    
    /* Several methods for non-equilibrium molecular dynamics were implemented
     * in the context of viscosity calculation: (i) weak-coupling to an external
     * bath; (ii) periodic perturbation method; and (iii) internal reservoir (reverse NEMD). 
     * Method (i) has several problems and convergence is very difficult; method (ii) works
     * well at least for LJ fluids, but can be in principle applied to (small)
     * molecules; method (iii) works well for LJ fluid under NVT boundary conditions,
     * and can theoretically be applied for molecules, but some implementation is needed.
     * 
     * In principle, NEMD can be used to calculate other properties such as thermal
     * conductivities, electric conductivity... Therefore, this file can be extended
     * to incorporate more combinations of methods and properties. The input flags 
     * give the user a lot of freedom on the choice of those methods, but care has to be 
     * taken as some of the methods do not work for some properties.
     * 
     * In the input file the user will see the block (that can be extended to incorporate
     * new features):
     * 
     *  
     * NEMD
        # NEMD 0,1 controls the use of non-equilibrium molecular dynamics.
        #    0 : not used [default]
        #    1 : nemd is used
        # PROPERTY 0- select property to calculate
        #    0 : viscosity      
        #    1 : thermal conductivity
        # METHOD 0- select method of NEMD.
        #    0 : internal reservoir method 
        #    1 : weak-coupling method 
        #    2 : periodic perturbation method
        # SLABNUM >=1 number of slabs used in the discretization along z-direction.
        #             the effective number is 2xSLABNUM due to periodicity
        # PERTFRQ >=1 perturbation frequency: apply perturbation every PERTFRQth timestep
        # REFTEMP >=0 reference temperature for thermal conductivity calculation
        # AMPLI >=0 amplitude applied field
        # STDYAFT >=0 first STDYAFTth steps do not contribute for accumulated averages
        # WRITE >=1 write flux and average velocities to special trajectory every WRITEth timestep
        # NEMD     PROPERTY  METHOD
              1         0        0
        # SLABNUM  PERTFRQ  REFTEMP  AMPBATH   STDYAFT   WRITE
         10       20   298.15       10      1000     200
        END
     * 
     * 
     * That will not be exactly like that for now as some methods were not yet tested
     * but the example above serves as a guide for future implementations. For now, the
     * block will be simpler and will only consider viscosity calculations with the two
     * methods that were found to be robust.
     * 
     * NEMD
        # NEMD 0,1 controls the use of non-equilibrium molecular dynamics.
        #    0 : not used [default]
        #    1 : nemd is used
        # PROPERTY 0- select property to calculate
        #    0 : viscosity      
        # METHOD 0- select method of NEMD.
        #    0 : periodic perturbation method (PPM)
        #    1 : internal reservoir method (IRM)
        # SLABNUM >=1 number of slabs used in the discretization along z-direction.
        #             the effective number is 2xSLABNUM due to periodicity
        # PERTFRQ >=1 perturbation frequency: apply perturbation every PERTFRQth timestep
        #             [this flag is ignored by the PPM method, but a value must be provided]
        # AMPLI   >=0 amplitude of applied field
        #             [this flag is ignored by the IRM method, but a value must be provided]
        # STDYAFT >=0 first STDYAFTth steps do not contribute for accumulated averages
        # WRITE >=1 write flux and average velocities to special trajectory every WRITEth timestep
        # NEMD     PROPERTY  METHOD
              1         0        0
        # SLABNUM  PERTFRQ    AMPLI   STDYAFT   WRITE
         10       20          10      1000     200
        END
     *  
     * 
     * 
     * 
     */
    
    
    /*
     * METHOD SELECTION 
     * 
     */
    
    int prop_vs_method = 0;
    if (sim.param().nemd.property == 0){
      if (sim.param().nemd.method == 0)
        prop_vs_method = 0;
      if (sim.param().nemd.method == 1)
        prop_vs_method = 1;
      //enter selection for other methods for this property here
    }
    //enter selection for other properties with respective methods here
    

    switch(prop_vs_method) {
      
      
      case 0:
      {
        /*
         *
         * This method consists of an applied periodic acceleration field
         * as described in JCC 20, 786 (1999), which gives a periodic velocity
         * profile corresponding to the velocity profile obtained directly by
         * scaling the velocities as described in reference "Berendsen in
         * M. Meyer and V. Pontikis, Computer simulation in Material Science
         * 139-155 (1991).
         *
         *
         */
        unsigned int nslab = 2 * sim.param().nemd.slabnum;
        double ampbath = sim.param().nemd.ampbath; // AAF: applied acceleration field
        double k = 2 * math::Pi / conf.current().box(2)(2); // k = 2*pi*L_z
        
        
        

        // In this case ampbath is an acceleration amplitude

        for(unsigned int i = 0; i < topo.num_atoms(); ++i){
            math::Vec pos = conf.current().pos(i);
            periodicity.put_into_positive_box(pos);
            double zbath = pos(2);
            double xforce = ampbath * cos(k*zbath) * topo.mass(i); //F_i = A0 cos(2pi*z_i/L_z) * mass
            math::Vec shearforce(xforce, 0.0, 0.0);
            
            
            conf.current().force(i) += shearforce;
        
        }
        
        /*
         * This grid will just store the velocities in such a way
         * that we can easily avarage over frames (after stdyaft steps)
         * and over the molecules of the grid.
         * That is not really necessary to calculate the amplitudes
         * but can be useful to visualize the effect of the periodic perturbation.
         * 
         */

        std::vector<std::vector<unsigned int> > grid(nslab);
        for(unsigned int i = 0; i < grid.size(); ++i) {
          grid[i].clear();
        }
        
           
        //DEBUG
        /*double vel =0.0;
        for(unsigned int i=0; i<topo.num_atoms(); ++i){
        math::Vec pos = conf.current().pos(i);
        periodicity.put_into_positive_box(pos);
        double zbath = pos(2);
        double a = conf.current().vel(i)(0) * cos(zbath*k);
        std::cout << "CHECK   " << a << "   " <<conf.current().vel(i)(0)<< "   " << zbath << "\n";
        vel += a;
        }
        vel /= topo.num_atoms();
        std::cout << vel << "   "  << k << "   " << "\n";
        */
        //ENDDEBUG


      // Accumulate after stdyaft steps
        if(sim.steps() > sim.param().nemd.stdyaft) {
          //Put atom indexes into grid

          for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
            math::Vec pos = conf.current().pos(i);
            periodicity.put_into_positive_box(pos);
            double z = pos(2);
            int bin = int(z / conf.current().box(2)(2) * nslab);
            grid[bin].push_back(i);
          }

          for(unsigned int slice = 0; slice < nslab; ++slice) {
            double vel_i = 0.0;
            int counter = 0;
            for(unsigned int i = 0; i < grid[slice].size(); ++i) {
              vel_i += conf.current().vel(grid[slice][i])(0);
              counter++;
            }
            double av_vel = vel_i / counter;
            conf.special().nemd_conf.stored_data_per_bin[slice] += av_vel;
            //vx_per_slab.push_back(av_vel);
          }
        }


        break;
      }
      
      case 1:
      {
        /*
         * Internal reservoir method: not applicable to molecules in the 
         * current stage. Viscosity of LJ fluids can be accurately calculated
         * using NVT ensemble. The details of the method can be found in reference
         * Muller-Plathe, F. Phys. Rev. E, volume 59, number 5, 1999.
         * Although not implemented here, an extension of the present method
         * for the treatment of molecular fluids can be found in ref.
         * Bordat, P. and Muller-Plathe, F. J. Chem. Phys. 116, 3362.  
         */
        //number of grid points
        unsigned int nslab = 2 * sim.param().nemd.slabnum;
        //double zlength = conf.current().box(2)(2);
        //double width = zlength/nslab;

        // constructing a grid with nslab number of grid points
        std::vector<std::vector<unsigned int> > grid(nslab);

        for(unsigned int i = 0; i < grid.size(); ++i) {
          grid[i].clear();
        }

        //Put atom indexes into grid
        for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
          math::Vec pos = conf.current().pos(i);
          periodicity.put_into_positive_box(pos);
          double z = pos(2);
          int bin = int(z / conf.current().box(2)(2) * nslab);
          grid[bin].push_back(i);
        }


        /*
        std::cout << "DEBUG" << "\n";
        for(unsigned int i = 0; i < grid.size(); ++i) {
          std::cout << grid[i].size() << "   " ;
        }
        std::cout << "\n";
        */
        if((sim.steps() % sim.param().nemd.pertfrq) == 0) {        
        // loop over member of 1st slice (grid[0])
        // and get the atom/molecule with most negative x-component
        double min_vel = 1E9; // Big positive number
        int min_index = -1; // Sentinel
        for(unsigned int i = 0; i < grid[0].size(); ++i) {
          //double vel_i = math::abs(conf.current().vel(grid[0][i]);
          double vel_i = conf.current().vel(grid[0][i])(0);
          if(min_vel > vel_i) {
            min_vel = vel_i;
            min_index = i;
          }
        }

        // loop over member of Lz/2 slice (grid[nslab/2])
        // and get the atom/molecule with most positive x-component
        double max_vel = -1E9; // Big negative number
        int max_index = -1; // Sentinel
        for(unsigned int i = 0; i < grid[nslab / 2].size(); ++i) {
          double vel_i = conf.current().vel(grid[nslab / 2][i])(0);
          if(max_vel < vel_i) {
            max_vel = vel_i;
            max_index = i;
          }
        }


        if (max_index != -1 && min_index != -1){
        // swap x-velocities of slice 1 and slice Lz/2; // std::swap(a,b);
        double x1 = conf.current().vel(grid[0][min_index])(0);
        double y1 = conf.current().vel(grid[0][min_index])(1);
        double z1 = conf.current().vel(grid[0][min_index])(2);
        double x2 = conf.current().vel(grid[nslab / 2][max_index])(0);
        double y2 = conf.current().vel(grid[nslab / 2][max_index])(1);
        double z2 = conf.current().vel(grid[nslab / 2][max_index])(2);
        math::Vec vel1(x2, y1, z1);
        math::Vec vel2(x1, y2, z2);

        if (x1 < x2){
        conf.current().vel(grid[0][min_index]) = vel1;
        conf.current().vel(grid[nslab / 2][max_index]) = vel2;

        double mass_particle = topo.mass(grid[0][min_index]);
        conf.special().nemd_conf.Px += mass_particle * (x2 - x1);
        }
        
        if(topo.mass(grid[0][min_index]) != topo.mass(grid[nslab / 2][max_index]))
          io::messages.add("NEMD_interaction", "Cannot swap velocities of particles with different mass", io::message::error);

        }
        } 
            
        if(sim.steps() > sim.param().nemd.stdyaft) {
        for(unsigned int slice = 0; slice < nslab; ++slice) {
          double vel_i = 0.0;
          int counter = 0;
          for(unsigned int i = 0; i < grid[slice].size(); ++i) {
            vel_i += conf.current().vel(grid[slice][i])(0);
            counter++;
          }
          double av_vel = vel_i / counter;
          conf.special().nemd_conf.stored_data_per_bin[slice] += av_vel;
          //conf.special().nemd_conf.stored_data_per_bin[slice]=av_vel;
        }
      }


        //math::Vec vel = conf.current().vel(max_index[1]);
        //conf.current().vel(max_index[1]) = conf.current().vel(max_index[10]);
        //conf.current().vel(max_index[10]) = vel;

        // get maximum for every slice
        /*
        std::vector<double> max_vel(num_grid_pints, 0.0);
                std::vector<unsigned int> max_index(num_grid_pints, 0.0);
        for(unsigned int slice = 0; slice < num_grid_points; ++slice) {
          for(unsigned int i = 0; i < grid[slice].size(); ++i) {
            double vel_i = math::abs(conf.current().vel(grid[slice][i]);
            if(max_vel[slice] < vel_i) {
              max_vel[slice] = vel_i;
                      max_index[slice] = i;
            }
          }
        }
         */

        break;
      }
      
      //Enter future algorithms here
      default: break;

    }
  return 0;
}

