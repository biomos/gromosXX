/**
 * @file lincs.cc
 * contains the template methods for
 * the class Lincs.
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../interaction/interaction.h"
#include "../../interaction/interaction_types.h"

#include "../../math/periodicity.h"

#include "../../algorithm/constraints/lincs.h"

#include "../../util/template_split.h"
#include "../../util/error.h"
#include "../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

struct coupling_struct
{
  std::vector<double> a;
};

/**
 * Constructor.
 */
algorithm::Lincs::Lincs()
  : Algorithm("Lincs")
{}

/**
 * Destructor.
 */
algorithm::Lincs::~Lincs()
{
}

int _solve_lincs(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::vector<topology::two_body_term_struct> const & constr,
		  math::VArray & B,
		  std::vector<coupling_struct> & A,
		  math::SArray rhs[],
		  math::SArray & sol,
		  topology::Compound::lincs_struct const & lincs,
		  int lincs_order,
		  unsigned int offset = 0)
{
  DEBUG(8, "LINCS SOLVE");
  
  unsigned int w = 1;

  const unsigned int num_constr = unsigned(constr.size());
  math::VArray & pos = conf.current().pos;
  
  for(unsigned int rec=0; int(rec) < lincs_order; ++rec){
    for(unsigned int i=0; i < num_constr; ++i){
      DEBUG(9, "rhs[" << i << "] = " << rhs[1-w](i));
      rhs[w](i) = 0.0;

      for(unsigned int n=0; n < lincs.coupled_constr[i].size(); ++n){
	rhs[w](i) += A[i].a[n] * rhs[1-w](lincs.coupled_constr[i][n]);
      }

      sol(i) = sol(i) + rhs[w](i);
    }
    w = 1 - w;
  }

  for(unsigned int i=0; i < num_constr; ++i){
    pos(constr[i].i + offset) -= B(i) * lincs.sdiag[i] / 
      topo.mass()(constr[i].i + offset) * sol(i);
    pos(constr[i].j + offset) += B(i) * lincs.sdiag[i] / 
      topo.mass()(constr[i].j + offset) * sol(i);
  }
    
  return 0;
}

template<math::boundary_enum bound>
int _lincs(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::vector<topology::two_body_term_struct> const & constr,
		  topology::Compound::lincs_struct const & lincs,
		  int lincs_order,
		  util::Algorithm_Timer & m_timer,
		  unsigned int offset = 0)
{
  
  const unsigned int num_constr = unsigned(constr.size());
  math::VArray & old_pos = conf.old().pos;
  math::VArray & pos = conf.current().pos;

  math::Vec r, ref_r;

  math::Periodicity<bound> periodicity(conf.current().box);
  
  math::VArray B(num_constr);

  std::vector<coupling_struct> A(num_constr);

  math::SArray rhs[2];
  rhs[0].resize(num_constr);
  rhs[1].resize(num_constr);

  math::SArray sol(num_constr);

  for(unsigned int i=0; i<num_constr; ++i){
    periodicity.nearest_image(old_pos(constr[i].i + offset), 
			      old_pos(constr[i].j + offset), ref_r);

    DEBUG(12, "i=" << constr[i].i << " j=" << constr[i].j << " offset=" << offset);
    DEBUG(12, "pos i = " << math::v2s(old_pos(constr[i].i + offset)));
    DEBUG(12, "pos j = " << math::v2s(old_pos(constr[i].j + offset)));
    DEBUG(12, "ref_r = " << math::v2s(ref_r));
    
    B(i) = ref_r / sqrt(math::abs2(ref_r));
    DEBUG(12, "B(" << i << ") = " << B(i)(0) << " / " << B(i)(1) << " / " << B(i)(2));
  }
  
  for(unsigned int i=0; i<num_constr; ++i){
    // coupled distance constraints
    const unsigned int ccon = unsigned(lincs.coupled_constr[i].size());

    for(unsigned int n=0; n < ccon; ++n){
      DEBUG(8, "A[" << i << "]: coupled: " << lincs.coupled_constr[i][n]);
      DEBUG(8, "constraint length: " << topo.bond_types_harm()[constr[i].type].r0);
      
      A[i].a.push_back(lincs.coef[i][n] * 
		       math::dot(B(i), B(lincs.coupled_constr[i][n])));
    }
    
    periodicity.nearest_image(pos(constr[i].i + offset), pos(constr[i].j + offset), r);

    rhs[0](i) = lincs.sdiag[i] * 
      (math::dot(B(i), r) - topo.bond_types_harm()[constr[i].type].r0);

    sol(i) = rhs[0](i);

  }

  _solve_lincs(topo, conf, sim, constr, B, A, rhs, sol, lincs, lincs_order, offset);
  
  // correction for rotational lengthening
  double p = 0.0;
  unsigned int count = 0;
  
  for(unsigned int i=0; i<num_constr; ++i){

    periodicity.nearest_image(pos(constr[i].i + offset), pos(constr[i].j + offset), r);

    const double diff = 2 * topo.bond_types_harm()[constr[i].type].r0 * topo.bond_types_harm()[constr[i].type].r0 -
      math::abs2(r);
    if (diff > 0.0)
      p = sqrt(diff);
    else{
      p = 0;
      ++count;
    }
    
    rhs[0](i) = lincs.sdiag[i] * (topo.bond_types_harm()[constr[i].type].r0 - p);
    sol(i) = rhs[0](i);

  }

  _solve_lincs(topo, conf, sim, constr, B, A, rhs, sol, lincs, lincs_order, offset);

  if (count)
    std::cout << "LINCS:\ttoo much rotation in " << count << " cases!\n";

  return 0;
}

template<math::boundary_enum B, math::virial_enum V>
void _solvent(topology::Topology & topo,
		     configuration::Configuration & conf,
		     simulation::Simulation & sim,
		     util::Algorithm_Timer & m_timer, int & error)
{
  m_timer.start("solvent");
  // the first atom of a solvent
  unsigned int first = unsigned(topo.num_solute_atoms());

  // for all solvents
  for(unsigned int i=0; i<topo.num_solvents(); ++i){
    
    // loop over the molecules
    for(unsigned int nm=0; nm<topo.num_solvent_molecules(i);
	++nm, first+=topo.solvent(i).num_atoms()){

      _lincs<B>(topo, conf, sim, topo.solvent(i).distance_constraints(),
		topo.solvent(i).lincs(), 
		sim.param().constraint.solvent.lincs_order,
		m_timer, first);
    }
  }
 
  m_timer.stop("solvent");
  error = 0;
}

/**
 * apply the Lincs algorithm
 */
int algorithm::Lincs::apply(topology::Topology & topo,
			    configuration::Configuration & conf,
			    simulation::Simulation & sim)
{
  DEBUG(7, "applying LINCS");
  m_timer.start();
  
  // check whether we lincs solute
  if (topo.solute().distance_constraints().size() && 
      sim.param().constraint.solute.algorithm == simulation::constr_lincs &&
      sim.param().constraint.ntc > 1){

    DEBUG(8, "\twe need to lincs SOLUTE");
    m_timer.start("solute");
    SPLIT_BOUNDARY(_lincs, topo, conf, sim, topo.solute().distance_constraints(),
		   topo.solute().lincs(),
		   sim.param().constraint.solute.lincs_order, m_timer);
    m_timer.stop("solute");
 
  }

  // SOLVENT
  if (sim.param().system.nsm &&
      sim.param().constraint.solvent.algorithm == simulation::constr_lincs){

    DEBUG(8, "\twe need to lincs SOLVENT");
    int error = 0;
    SPLIT_VIRIAL_BOUNDARY(_solvent,
			  topo, conf, sim, m_timer, error);
    
  }
  
  // inclsed velocity:
  // stochastic dynamics, energy minimisation, analysis needs to shake without
  // velocity correction (once; it shakes twice...)
  if (!sim.param().stochastic.sd && !sim.param().minimise.ntem &&
      !sim.param().analyze.analyze) {
        std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.current().vel(*it) = (conf.current().pos(*it) - conf.old().pos(*it)) / 
	sim.time_step_size();
    }
  }

  // return success!
  m_timer.stop();
  return 0;	   
}

static void _setup_lincs(topology::Topology const & topo,
			 topology::Compound::lincs_struct & lincs,
			 std::vector<topology::two_body_term_struct> const & constr,
			 unsigned int offset = 0)
{
  const unsigned int num_constr = unsigned(constr.size());  
  lincs.coupled_constr.resize(num_constr);
  lincs.coef.resize(num_constr);

  for(unsigned int i=0; i < num_constr; ++i){
    // diagonal matrix
    lincs.sdiag.push_back(1.0 / sqrt(1.0 / topo.mass()(constr[i].i + offset) +
				     1.0 / topo.mass()(constr[i].j + offset)));

    DEBUG(8, "diag[" << i << "] = " << lincs.sdiag[i]);

  }
  
  for(unsigned int i=0; i < num_constr; ++i){

    // connected constraints
    for(unsigned int j=i+1; j < num_constr; ++j){

      int con = -1;
      
      if (constr[i].i == constr[j].i)
	con = constr[i].i;
      else if (constr[i].j == constr[j].i)
	con = constr[i].j;
      else if (constr[i].i == constr[j].j)
	con = constr[i].i;
      else if (constr[i].j == constr[j].j)
	con = constr[i].j;
      
      if (con != -1){
	DEBUG(8, "constraint " << i << ": " << constr[i].i << " - " << constr[i].j << "\tconnected with");
	DEBUG(8, "constraint " << j << ": " << constr[j].i << " - " << constr[j].j);
	
	DEBUG(8, "\tthrough atom " << con);
	
	lincs.coupled_constr[i].push_back(j);
	lincs.coupled_constr[j].push_back(i);

	double c = 1.0 / topo.mass()(con + offset) * lincs.sdiag[i] * lincs.sdiag[j];
	
	if (constr[i].i == constr[j].i ||
	    constr[i].j == constr[j].j)
	  c *= -1;
	
	lincs.coef[i].push_back(c);
	lincs.coef[j].push_back(c);

	DEBUG(8, "coef[" << i << "] = " << c);
	
      }
    } 
  }
}


int algorithm::Lincs::init(topology::Topology & topo,
			   configuration::Configuration & conf,
			   simulation::Simulation & sim,
			   std::ostream & os,
			   bool quiet)
{
  if (!quiet){
    os << "LINCS\n"
	      << "\tsolute\t";
    if (sim.param().constraint.solute.algorithm == simulation::constr_lincs)
      os << "ON\n";
    else os << "OFF\n";
  
    os << "\t\torder = " << sim.param().constraint.solute.lincs_order << "\n";
  
    os << "\tsolvent\t";
    if (sim.param().constraint.solvent.algorithm == simulation::constr_lincs)
      os << "ON\n";
    else os << "OFF\n";
    
    os << "\t\torder = " << sim.param().constraint.solvent.lincs_order << "\n";
  }

  if (sim.param().constraint.solute.algorithm == simulation::constr_lincs) {
    // loop over the constraints to find out which atoms are constrained
    std::vector<topology::two_body_term_struct>::const_iterator
    it = topo.solute().distance_constraints().begin(),
            to = topo.solute().distance_constraints().end();
    for (; it != to; ++it) {
      constrained_atoms().insert(it->i);
      constrained_atoms().insert(it->j);
    }
  }

  if (sim.param().constraint.solvent.algorithm == simulation::constr_lincs) {
    // this means that all solvent atoms are constrained. Add the to the list
    for (unsigned int i = topo.num_solute_atoms(); i < topo.num_atoms(); ++i) {
      constrained_atoms().insert(i);
    }
  }

  // setup lincs
  DEBUG(8, "setting up lincs");
  _setup_lincs(topo, topo.solute().lincs(),
	       topo.solute().distance_constraints());
  
  // the first atom of a solvent
  unsigned int first = topo.num_solute_atoms();

  // for all solvents
  for(unsigned int i=0; i<topo.num_solvents(); ++i){  //Initialize lincs for each solvent
    if (topo.num_solvent_molecules(i) != 0){
      _setup_lincs(topo, topo.solvent(i).lincs(),
		    topo.solvent(i).distance_constraints(),
		    first);
      first+=topo.solvent(i).num_atoms();
    } else {
      io::messages.add("Solvent " + std::to_string(i) + ": No molecules found; does not initialize LINCS for solvent "  + std::to_string(i) + ".", 
    		 	  "algorithm::Lincs::init" , io::message::warning);
    }
  }
  
  //================================================================================
  // initialize positions / velocities
  //================================================================================

  if (sim.param().start.shake_pos){
    if (!quiet)
      os << "\n\tshaking (lincs) initial positions\n";

    // old and current pos and vel are the same for constrained atoms...
    std::set<unsigned int>::const_iterator it = constrained_atoms().begin(),
            to = constrained_atoms().end();
    for (; it != to; ++it) {
      conf.old().pos(*it) = conf.current().pos(*it);
      conf.old().vel(*it) = conf.current().vel(*it);
    }

    // lincs the current ones
    apply(topo, conf, sim);

    it = constrained_atoms().begin();
    for (; it != to; ++it) {
      // restore the velocities
      conf.current().vel(*it) = conf.old().vel(*it);
      // take a step back
      conf.old().pos(*it) = conf.current().pos(*it);
    }
    
    if (sim.param().start.shake_vel){
      if (!quiet)
	os << "\tshaking (lincs) initial velocities\n";

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        conf.current().pos(*it) = conf.old().pos(*it) -
                sim.time_step_size() * conf.old().vel(*it);
      }
    
      // lincs again
      apply(topo, conf, sim);

      it = constrained_atoms().begin();
      for (; it != to; ++it) {
        // restore the positions
        conf.current().pos(*it) = conf.old().pos(*it);
        // velocities are in opposite direction (in time)
        conf.current().vel(*it) = -1.0 * conf.current().vel(*it);
        conf.old().vel(*it) = conf.current().vel(*it);
      }
    } // lincs initial velocities
  } else if (sim.param().start.shake_vel) {
    io::messages.add("shaking velocities without shaking positions illegal.",
            "lincs", io::message::error);
  }

  if (!quiet)
    os << "END\n";

  return 0;
}


