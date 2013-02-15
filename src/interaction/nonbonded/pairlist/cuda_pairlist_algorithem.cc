/**
 * @file standard_cuda_algorithm.cc
 * standard pairlist algorithm
 */
#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../math/periodicity.h"

#include "../../../interaction/nonbonded/pairlist/pairlist.h"

#include "../../../interaction/nonbonded/pairlist/pairlist_algorithm.h"
#include "../../../interaction/nonbonded/pairlist/cuda_pairlist_algorithem.h"

#include "../../../util/debug.h"
#include "../../../util/template_split.h"

#ifdef XXCUDA
#include "../../../cuda/gromos_cuda.h"
#include "../../../cuda/c_pairlist.h"
#endif

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE pairlist

interaction::cuda_pairlist::
cuda_pairlist()
: interaction::Pairlist_Algorithm(), 
 m_solvent_solvent_timing(0.0)
#ifdef XXCUDA
, gpumemory(20)
#endif
 {}

interaction::cuda_pairlist::~cuda_pairlist() {

#ifdef XXCUDA
   gcuda::freecuda(gpumemory, atomic);
#endif
}
// init
int interaction::cuda_pairlist::init(topology::Topology &topo,
         configuration::Configuration &conf,
	 simulation::Simulation &sim,
         std::ostream &os,
	 bool quiet) {

  #ifdef XXCUDA
    int number_gpus = gcuda::has_gpu();
    if (number_gpus != 0) {
        gcuda::devices available_devices;
        available_devices = gcuda::get_dev_properties(number_gpus);
        if (available_devices.error != true) {
            os << "CUDA Pairlist\n";
            os << "\tAvailable GPUs:\n";
            for (int i = 0; i < number_gpus; i++) { // signed because number_gpus -1 = error
                os << "\t\tName: " << available_devices.property[i].name << "\n";
                os << "\t\tCompute capability: ";
                os << available_devices.property[i].major;
                os << "." <<  available_devices.property[i].minor << "\n";
                os << "\t\tAvailable Memory: " << available_devices.property[i].totalGlobalMem << "\n";
                os << "\t\tMode: " << available_devices.property[i].computeMode << "\n";
                os << "\n";
            }
        } else {
             io::messages.add("CUDA ERROR", "create_nonbonded - cuda init",
             io::message::error);
        }
    } else {
        io::messages.add("NO GPU FOUND", "create_nonbonded - cuda init",
             io::message::error);
    }

    // chargegroups or atomic ?
    if (sim.param().pairlist.atomic_cutoff){

        atomic = true;
        
        // exclusions
        gcuda::TwoDArray<unsigned int> h_exclusions(topo.num_solute_atoms(), 20);
        prepare_exclusions(h_exclusions, topo);

        // set cutoffs
        set_cutoff(sim.param().pairlist.cutoff_short, sim.param().pairlist.cutoff_long);

        // host array pointers
        gpumemory.add("host_exclusions", (void*)*h_exclusions.ptr , h_exclusions.pitch);
        gpumemory.add("host_exclusions_elements", (void*)h_exclusions.elements, 4);

        // calculate results pitch
        size_t host_results_pitch = sim.param().pairlist.pairlist_height * sizeof(unsigned int);

        // GPU function
        gcuda::initcuda_atomic(gpumemory, topo.num_atoms(), topo.num_solute_atoms(), cutoff, host_results_pitch);
        //

            os << "CUDA ATOMIC INIT";

    } else {

        atomic = false;

        // chargegroup array, containing an index and the corresponding first atom number of the chargegroup
            gcuda::TwoDArray<unsigned int> h_chargegroups(topo.num_chargegroups(), 2);

         // fill chargegroupe array
            for (unsigned int i = 0; i < topo.num_chargegroups(); i++) {
                h_chargegroups.ptr[i][0] = topo.chargegroup(i);
                h_chargegroups.ptr[i][1] = topo.chargegroup(i+1);
            }

        // exclusions
        gcuda::TwoDArray<unsigned int> h_exclusions(topo.num_solute_atoms(), 20);
        prepare_exclusions(h_exclusions, topo);

        // set cutoffs
        set_cutoff(sim.param().pairlist.cutoff_short, sim.param().pairlist.cutoff_long);

        // host array pointers
        gpumemory.add("host_chargegroups", (void*)*h_chargegroups.ptr , h_chargegroups.pitch);
        gpumemory.add("host_exclusions", (void*)*h_exclusions.ptr , h_exclusions.pitch);
        gpumemory.add("host_exclusions_elements", (void*)h_exclusions.elements, 4);

        // calculate results pitch
        size_t host_results_pitch = sim.param().pairlist.pairlist_height * sizeof(unsigned int);
        // GPU function
        gcuda::initcuda_cg(gpumemory, topo.num_atoms(), topo.num_solute_atoms(), topo.num_chargegroups(), topo.num_solute_chargegroups(), cutoff, host_results_pitch);
        //
            os << "CUDA CHARGEGROUP INIT";
    }

    return 0;
#else
      io::messages.add("CUDA SUPPORT NOT COMPILED", "create_nonbonded - cuda init",
             io::message::error);
    return 0;
#endif
}

// prepare
int interaction::cuda_pairlist::prepare(topology::Topology & topo,
            configuration::Configuration & conf,
	    simulation::Simulation &sim){

   // CUDA is unprepared :)
 
    return 0;

}

// update
void interaction::cuda_pairlist::update(topology::Topology & topo,
	    configuration::Configuration & conf,
            simulation::Simulation & sim,
            interaction::PairlistContainer & pairlist,
            unsigned int begin,
            unsigned int end,
	    unsigned int stride) {

  if (sim.param().pairlist.atomic_cutoff){
    // atomic or chargegroup based
    update_atomic(topo, conf, sim, pairlist, begin, end, stride);
  }
  else{
    update_cg(topo, conf, sim, pairlist, begin, end, stride);
  }

}

void interaction::cuda_pairlist::update_cg(topology::Topology & topo,
	    configuration::Configuration & conf,
            simulation::Simulation & sim,
            interaction::PairlistContainer & pairlist,
            unsigned int begin,
            unsigned int end,
	    unsigned int stride) {

#ifdef XXCUDA
    //timer
     timer().start("pairlist");
    //update chargegroup
     
     // gpu functionv   // nr atoms nr solute nr solvent atoms do not copy nratoms but copy only solute or only solvent (ajust pointers n pitches)
            gcuda::update_cg(gpumemory, &conf.current().pos(0)(0), conf.current().box, topo.num_atoms(), topo.num_solute_atoms(), topo.num_chargegroups(), &pairlist.solute_short.container(), &pairlist.solute_long.container(), &pairlist.solvent_short.container(), &pairlist.solvent_long.container());
//
  
      std::cout << "############SOLUTE SHORT###############" << std::endl;
  std::cout <<  pairlist.solute_short.size()  << std::endl;
  std::cout << "############SOLUTE LONG###############" << std::endl;
  std::cout <<  pairlist.solute_long.size()  << std::endl;
   std::cout << "############SOLVENT SHORT###############" << std::endl;
  std::cout <<  pairlist.solvent_short.size()  << std::endl;
  std::cout << "############SOLVENT LONG###############" << std::endl;
  std::cout <<  pairlist.solvent_long.size()  << std::endl;
  std::cout << "#######################################" << std::endl;
             //test
             std::cout << "############SOLUTE SHORT###############" << std::endl;
            for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
                std::cout << "# " << i << " #" ;
                for (unsigned int j = 0; j < pairlist.solute_short.size(i) ; j++) {
                        std::cout << pairlist.solute_short.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                  std::cout << "############SOLUTE LONG###############" << std::endl;
            for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
                std::cout << "# " << i << " #" ;
               for (unsigned int j = 0; j < pairlist.solute_long.size(i) ; j++) {
                        std::cout << pairlist.solute_long.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                     std::cout << "############SOLVENT SHORT###############" << std::endl;
            for (unsigned int i = topo.num_solute_atoms(); i < topo.num_solvent_atoms(); i++) {
                std::cout << "# " << i << " #" ;
                for (unsigned int j = 0; j < pairlist.solvent_short.size(i) ; j++) {
                        std::cout << pairlist.solvent_short.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                  std::cout << "############SOLVENT LONG###############" << std::endl;
            for (unsigned int i = topo.num_solute_atoms(); i < topo.num_solvent_atoms(); i++) {
                std::cout << "# " << i << " #" ;
               for (unsigned int j = 0; j < pairlist.solvent_long.size(i) ; j++) {
                        std::cout << pairlist.solvent_long.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
            std::cout << "##############################" << std::endl;
            exit(0);
            
     //timer
    timer().stop("pairlist");
    return;

#else
      io::messages.add("CUDA SUPPORT NOT COMPILED", "create_nonbonded - cuda update chargegroups",
             io::message::error);
    return;
#endif
}


void interaction::cuda_pairlist::update_atomic(topology::Topology & topo,
	    configuration::Configuration & conf,
            simulation::Simulation & sim,
            interaction::PairlistContainer & pairlist,
            unsigned int begin,
            unsigned int end,
	    unsigned int stride) {

#ifdef XXCUDA
    //timer
     timer().start("pairlist");
    //update atomic
    // gpu functionv   // nr atoms nr solute nr solvent atoms do not copy nratoms but copy only solute or only solvent (ajust pointers n pitches)
            gcuda::update_atomic(gpumemory, &conf.current().pos(0)(0), conf.current().box, topo.num_atoms(), topo.num_solute_atoms(), &pairlist.solute_short.container(), &pairlist.solute_long.container(), &pairlist.solvent_short.container(), &pairlist.solvent_long.container());
    //

  /*            std::cout << "############SOLUTE SHORT###############" << std::endl;
  std::cout <<  pairlist.solute_short.size()  << std::endl;
  std::cout << "############SOLUTE LONG###############" << std::endl;
  std::cout <<  pairlist.solute_long.size()  << std::endl;
   std::cout << "############SOLVENT SHORT###############" << std::endl;
  std::cout <<  pairlist.solvent_short.size()  << std::endl;
  std::cout << "############SOLVENT LONG###############" << std::endl;
  std::cout <<  pairlist.solvent_long.size()  << std::endl;
  std::cout << "#######################################" << std::endl;
             //test
             std::cout << "############SOLUTE SHORT###############" << std::endl;
            for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
                std::cout << "# " << i << " #" ;
                for (unsigned int j = 0; j < pairlist.solute_short.size(i) ; j++) {
                        std::cout << pairlist.solute_short.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                  std::cout << "############SOLUTE LONG###############" << std::endl;
            for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
                std::cout << "# " << i << " #" ;
               for (unsigned int j = 0; j < pairlist.solute_long.size(i) ; j++) {
                        std::cout << pairlist.solute_long.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                     std::cout << "############SOLVENT SHORT###############" << std::endl;
            for (unsigned int i = topo.num_solute_atoms(); i < topo.num_solvent_atoms(); i++) {
                std::cout << "# " << i << " #" ;
                for (unsigned int j = 0; j < pairlist.solvent_short.size(i) ; j++) {
                        std::cout << pairlist.solvent_short.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
                  std::cout << "############SOLVENT LONG###############" << std::endl;
            for (unsigned int i = topo.num_solute_atoms(); i < topo.num_solvent_atoms(); i++) {
                std::cout << "# " << i << " #" ;
               for (unsigned int j = 0; j < pairlist.solvent_long.size(i) ; j++) {
                        std::cout << pairlist.solvent_long.container().ptr[i][j] << " " ;
                }
              std::cout << std::endl;
            }
            std::cout << "##############################" << std::endl;

     */
     //timer
    timer().stop("pairlist");
    return;

#else
      io::messages.add("CUDA SUPPORT NOT COMPILED", "create_nonbonded - cuda update atomic",
             io::message::error);
    return;
#endif
    
}

// update perturbed
void interaction::cuda_pairlist::update_perturbed(topology::Topology & topo,
	    configuration::Configuration & conf,
            simulation::Simulation & sim,
            interaction::PairlistContainer & pairlist,
            interaction::PairlistContainer & perturbed_pairlist,
            unsigned int begin,
            unsigned int end,
	    unsigned int stride){

    io::messages.add("work in progress :)", "create_nonbonded - cuda update_perturbed",
             io::message::error);
    return;
}


#ifdef XXCUDA
// std::set to array
void interaction::cuda_pairlist::prepare_exclusions(gcuda::TwoDArray<unsigned int> & earray, topology::Topology & topo) {

    for (unsigned int i = 0; i < topo.num_solute_atoms(); i++) {
        unsigned int j = 0;

        std::set<int>::reverse_iterator e = topo.all_exclusion(i).rbegin(), e_to = topo.all_exclusion(i).rend();
            for ( ; e != e_to; ++e) {
                earray.ptr[i][j] = *e;
                j++;
            }
        earray.elements[i] = j;
    }
}
#endif

//copied funtion will be changed later
bool interaction::cuda_pairlist
::excluded_solute_pair(topology::Topology & topo,
		       unsigned int i, unsigned int j)
{
  assert(i<j);

  std::set<int>::reverse_iterator
    e = topo.all_exclusion(i).rbegin(),
    e_to = topo.all_exclusion(i).rend();

  for( ; e != e_to; ++e){
    if (j > unsigned(*e)) break;
    if (j == unsigned(*e)){
      DEBUG(11, "\texcluded");
      return true;
    }

  }
  DEBUG(12, "\tnot excluded");
  return false;
}