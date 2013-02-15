/*
 * GROMOSXX Interface to HOOMD Code
 * By: Matthew Breeze (November 2009)
 * 
 * See documentation lyx/pdf
*/

#ifdef HAVE_HOOMD
#ifndef INCLUDED_HOOMD_GROMOSXX_INTERFACE_H
#define INCLUDED_HOOMD_GROMOSXX_INTERFACE_H

#include <HOOMD_GROMOSXX_processor.h>

#include <ParticleData.h>
#include <NeighborList.h>
#include <NeighborListNsqGPU.h>
#include <BinnedNeighborList.h>
#include <BinnedNeighborListGPU.h>

#include <boost/utility.hpp>
#include <stdheader.h>
#include <set>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <simulation/parameter.h>
#include <configuration/configuration.h>

#include <math/periodicity.h>

#include <interaction/nonbonded/pairlist/pairlist.h>
#include <interaction/nonbonded/pairlist/pairlist_algorithm.h>

namespace interaction
{
  /**
   * @class HOOMD_Pairlist_Algorithm
   * creates a pairlist.
   */
  class HOOMD_Pairlist_Algorithm : public Pairlist_Algorithm
  {
  protected:
    boost::shared_ptr<ParticleData> pdata; // Particle Data
	boost::shared_ptr<NeighborList> nlist; // Neighbour list algorithm
	// The chosen neighbour list implementation
	enum nlist_choice {
		bc, bg, nc, ng
	} choice;
	// Are the box dimensions the same since ParticleData was last built?
	bool box_size_same; 
	// HOOMD has an internal check to see if the pairlist has already 
	// been built for the given time step, so to ensure it always builds
	// it, we simply keep a counter and increment it every neighbour list 
	// update call
	unsigned int timestep; 
		
  public:
    /**
     * Constructor.
     */
    HOOMD_Pairlist_Algorithm(simulation::Simulation const & sim)
		   : Pairlist_Algorithm() {
		if (sim.param().hoomd.processor == simulation::unknown) {
			std::cerr << "HOOMD used without HOOMD block in input file" 
					<< std::endl;
			assert(false);
			return;
		}
		if (sim.param().hoomd.processor == simulation::cpu) {
			if (sim.param().pairlist.grid == 0) {
				choice = nc; // CPU, O(n^2) (not binned)
			} else {
				choice = bc; // CPU, O(n) (binned)
			}
		} else {
			if (sim.param().pairlist.grid == 0) {
				choice = ng; // GPU, O(n^2) (not binned)
			} else {
				choice = bg; // GPU, O(n) (binned)
			}
		}
		box_size_same = false; // Need to create new ParticleData
		timestep = 0;
	}
    
	/**
     * destructor. As shared_pointers are used no explicit destruction
	 * is needed
     */
    virtual ~HOOMD_Pairlist_Algorithm() {}

    /**
     * init
     */
    virtual int init(topology::Topology &topo, 
		     configuration::Configuration &conf,
		     simulation::Simulation &sim,
		     std::ostream &os = std::cout,
		     bool quiet = false) { 
		if (conf.boundary_type != math::rectangular) {
			std::cerr << "HOOMD only supports rectangular boxes" << std::endl;
			assert(false);
			return -1;
		}
		g_gpu_error_checking = false; 
		return 0; 
	};

    /**
     * prepare the pairlist(s). always called just before update
     */
    virtual int prepare(topology::Topology & topo,
			 configuration::Configuration & conf,
			 simulation::Simulation &sim) { 
		//double start_time = util::now();
		// if box size changed or if it could be rescaled by P coupling, need to 
		// recreate ParticleData object and consequently NeighbourList object
		if (!box_size_same || sim.param().pcouple.scale != math::pcouple_off) { 
			pdata = boost::shared_ptr<ParticleData>(
					new ParticleData(topo.num_chargegroups() - topo.num_solute_chargegroups(),
							      	 BoxDim(conf.current().box(0,0), 
								   	 		conf.current().box(1,1), 
											conf.current().box(2,2)
											), 
									 1, 0, 0, 0, 0, sim.proc->exec_conf));
			box_size_same = true;
			// Construct neighbour list algorithm object
			switch (choice) { 
				case bg: nlist = boost::shared_ptr<BinnedNeighborListGPU>(new BinnedNeighborListGPU(pdata, sim.param().pairlist.cutoff_short, 0.0)); break;
				case bc: nlist = boost::shared_ptr<BinnedNeighborList>(new BinnedNeighborList(pdata, sim.param().pairlist.cutoff_short, 0.0)); break;
				case ng: nlist = boost::shared_ptr<NeighborListNsqGPU>(new NeighborListNsqGPU(pdata, sim.param().pairlist.cutoff_short, 0.0)); break;
				case nc: nlist = boost::shared_ptr<NeighborList>(new NeighborList(pdata, sim.param().pairlist.cutoff_short, 0.0)); break;
			}
			if (choice == bg || choice == ng) // gpu-based 
				nlist->setStorageMode(NeighborList::full); // needs to be converted to half for GromosXX
			else // cpu-based (cpu-based can also do NeighbourList::full)
				nlist->setStorageMode(NeighborList::half);
		} 
		// Put molecules into box (copied from Extended_Grid_Pairlist_Algorithm::prepare_grid)
		// and whilst this is done put the first atom of each solvent molecule into HOOMD
		ParticleDataArrays arrays = pdata->acquireReadWrite();
		math::Periodicity<math::rectangular> periodicity(conf.current().box);
  	
		math::VArray &pos = conf.current().pos;
 		math::Vec v, v_box, trans;
  		topology::Chargegroup_Iterator cg_it = topo.chargegroup_it(topo.num_solute_chargegroups()),
    	cg_to = topo.chargegroup_end();
  		unsigned int i = 0;
  		// solvent chargegroups
  		for( ; cg_it != cg_to; ++cg_it, ++i){

    		// cog is first atom
    		v = pos(**cg_it);
    		v_box = v;
    		periodicity.put_into_box(v_box);
    		trans = v_box - v;
			
			// copy atomic coordinates into HOOMD
			arrays.x[i] = v_box(0); // GROMOS uses the first atom of each solvent molecule for pairlist
			arrays.y[i] = v_box(1);
			arrays.z[i] = v_box(2);

    		// loop over the atoms (i.e. translate all atoms in this solvent molecule to the amount needed
			// by the first atom for it to be inside the box)
    		topology::Atom_Iterator at_it = cg_it.begin(),
    		  at_to = cg_it.end();
   			for( ; at_it != at_to; ++at_it){
      			pos(*at_it) += trans;
    		} // atoms
  		} // solvent cg's
		pdata->release(); // release and let HOOMD transfer to the GPU
		//std::cout << "prepare coordinates: " << (util::now() - start_time) << std::endl;
		return 0; 
	}

    /**
     * update the pairlist. solvent-solvent short-range pairs only (single solvent type expected)
     */
    void update(topology::Topology & topo,
			configuration::Configuration & conf,
			simulation::Simulation &sim,
			interaction::PairlistContainer &pairlist,
			unsigned int begin, unsigned int end, 
			unsigned int stride) {
		// Generate neighbour list in HOOMD
		//double start_time = util::now();
		nlist->compute(++timestep);
		//std::cout << "computing pairlist: " << (util::now() - start_time) << std::endl;
		// Output list from HOOMD (ultimately the list would go straight to CUDA force routines on the GPU)
		//start_time = util::now();
		std::vector<std::vector<unsigned int> > const &list = nlist->getList();
		// convert from full to half list (does nothing if list is half originally)
		std::vector<std::set <unsigned int> > done(topo.num_chargegroups() - topo.num_solute_chargegroups());
		bool isfull = (choice == bg || choice == ng);
  		// empty the pairlist
  		pairlist.clear();
		
		// GROMOS XX requires atom-atom pairs and no (j,i) if (i,j) exists
		int first_cg = topo.num_solute_chargegroups();
		for(int i = 0, ito = list.size(); i < ito; i++) { // for each mol i
			int a1 = topo.chargegroups()[first_cg + i], // start atom
				a1to = topo.chargegroups()[first_cg + i + 1]; // one after end
			for(int j = 0, jto = list[i].size(); j < jto; j++) { // each mol j
				if (!isfull || done[list[i][j]].find(i) == done[list[i][j]].end()) { // if reversed pair does not exist, add regular pair 
					int a2 = topo.chargegroups()[first_cg + list[i][j]], // start atom
						a2to = topo.chargegroups()[first_cg + list[i][j] + 1]; // one after end
					for(int a1tmp = a1; a1tmp != a1to; a1tmp++) { // add atom pairs
						for (int a2tmp = a2; a2tmp != a2to; a2tmp++) { 
							pairlist.solvent_short[a1tmp].push_back(a2tmp);
						}
					}
					if (isfull)
						done[i].insert(list[i][j]); // add regular pair
				}
			}
		}
		//std::cout << "receiving and formatting pairlist: " << (util::now() - start_time) << std::endl;
	}

    /**
     * update the pairlist, separating perturbed and non-perturbed interactions
     */
    virtual void update_perturbed(topology::Topology & topo,
				  configuration::Configuration & conf,
				  simulation::Simulation & sim,
                                  interaction::PairlistContainer & pairlist,
				  interaction::PairlistContainer & perturbed_pairlist,
				  unsigned int begin, unsigned int end, 
				  unsigned int stride) {
      std::cerr << "update_perturbed not yet supported" << std::endl;
      assert(false);
    }
  };
} // interaction

#endif // INCLUDED_HOOMD_GROMOSXX_INTERFACE_H
#endif // HAVE_HOOMD
