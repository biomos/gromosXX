/**
 * @file lincs.tcc
 * contains the template methods for
 * the class Lincs.
 */

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
template<math::virial_enum do_virial>
algorithm::Lincs<do_virial>
::Lincs()
  : Algorithm("Lincs"),
    m_solvent_timing(0.0)
{
}

/**
 * Destructor.
 */
template<math::virial_enum do_virial>
algorithm::Lincs<do_virial>
::~Lincs()
{
}

static int _solve(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::vector<topology::two_body_term_struct> const & constr,
		  math::VArray & B,
		  std::vector<coupling_struct> & A,
		  math::SArray rhs[],
		  math::SArray & sol,
		  topology::Compound::lincs_struct const & lincs,
		  size_t const offset = 0)
{
  DEBUG(8, "LINCS SOLVE");
  
  size_t w = 1;

  const size_t num_constr = constr.size();
  math::VArray & pos = conf.current().pos;
  
  for(size_t rec=0; int(rec) < sim.param().shake.lincs_order; ++rec){
    for(size_t i=0; i < num_constr; ++i){
      DEBUG(9, "rhs[" << i << "] = " << rhs[1-w](i));
      rhs[w](i) = 0.0;

      for(size_t n=0; n < lincs.coupled_constr[i].size(); ++n){
	rhs[w](i) = rhs[w](i) + A[i].a[n] * rhs[1-w](lincs.coupled_constr[i][n]);
      }

      sol(i) = sol(i) + rhs[w](i);
    }
    w = 1 - w;
  }

  for(size_t i=0; i < num_constr; ++i){
    pos(constr[i].i + offset) -= B(i) * lincs.sdiag[i] / 
      topo.mass()(constr[i].i + offset) * sol(i);
    pos(constr[i].j + offset) += B(i) * lincs.sdiag[i] / 
      topo.mass()(constr[i].j + offset) * sol(i);
  }
    
  return 0;
  
}

template<math::boundary_enum b>
static int _lincs(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim,
		  std::vector<topology::two_body_term_struct> const & constr,
		  topology::Compound::lincs_struct const & lincs,
		  std::vector<interaction::bond_type_struct> const & param,
		  double & timing,
		  size_t const offset = 0)
{
  const double start = util::now();
  
  const size_t num_constr = constr.size();
  math::VArray const & old_pos = conf.old().pos;
  math::VArray & pos = conf.current().pos;

  math::Vec r, ref_r;

  math::Periodicity<b> periodicity(conf.current().box);
  
  math::VArray B(num_constr);

  std::vector<coupling_struct> A(num_constr);

  math::SArray rhs[2];
  rhs[0].resize(num_constr);
  rhs[1].resize(num_constr);

  math::SArray sol(num_constr);
  

  for(size_t i=0; i<num_constr; ++i){
    periodicity.nearest_image(old_pos(constr[i].i + offset), 
			      old_pos(constr[i].j + offset), ref_r);

    B(i) = ref_r / sqrt(math::dot(ref_r, ref_r));
    
    DEBUG(12, "B(" << i << ") = " << B(i)(0) << " / " << B(i)(1) << " / " << B(i)(2));
  }
  
  for(size_t i=0; i<num_constr; ++i){
    // coupled distance constraints
    const size_t ccon = lincs.coupled_constr[i].size();

    for(size_t n=0; n < ccon; ++n){
      DEBUG(8, "A[" << i << "]: coupled: " << lincs.coupled_constr[i][n]);
      DEBUG(8, "constraint length: " << param[constr[i].type].r0);
      
      A[i].a.push_back(lincs.coef[i][n] * 
		       math::dot(B(i), B(lincs.coupled_constr[i][n])));
    }
    
    periodicity.nearest_image(pos(constr[i].i + offset), pos(constr[i].j + offset), r);

    rhs[0](i) = lincs.sdiag[i] * 
      (math::dot(B(i), r) - param[constr[i].type].r0);

    sol(i) = rhs[0](i);

  }

  _solve(topo, conf, sim, constr, B, A, rhs, sol, lincs, offset);
  
  // correction for rotational lengthening
  double p;
  size_t count = 0;
  
  for(size_t i=0; i<num_constr; ++i){

    periodicity.nearest_image(pos(constr[i].i + offset), pos(constr[i].j + offset), r);

    const double diff = 2 * param[constr[i].type].r0 * param[constr[i].type].r0 -
      math::dot(r, r);
    if (diff > 0.0)
      p = sqrt(diff);
    else{
      p = 0;
      ++count;
    }
    
    rhs[0](i) = lincs.sdiag[i] * (param[constr[i].type].r0 - p);
    sol(i) = rhs[0](i);

  }

  _solve(topo, conf, sim, constr, B, A, rhs, sol, lincs, offset);

  if (count)
    std::cout << "LINCS:\ttoo much rotation in " << count << " cases!\n";

  timing += util::now() - start;

  return 0;
}

template<math::virial_enum do_virial, math::boundary_enum b>
static int _solvent(topology::Topology & topo,
		    configuration::Configuration & conf,
		    simulation::Simulation & sim,
		    std::vector<interaction::bond_type_struct> const & param,
		    double & timing)
{
  // the first atom of a solvent
  size_t first = topo.num_solute_atoms();

  // for all solvents
  for(size_t i=0; i<topo.num_solvents(); ++i){
    
    // loop over the molecules
    for(size_t nm=0; nm<topo.num_solvent_molecules(i);
	++nm, first+=topo.solvent(i).num_atoms()){

      _lincs<b>(topo, conf, sim, topo.solvent(i).distance_constraints(),
		topo.solvent(i).lincs(), param, timing, first);
      
    }
  }
 
  return 0;
  
}

/**
 * apply the Lincs algorithm
 */
template<math::virial_enum do_virial>
int algorithm::Lincs<do_virial>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "applying LINCS");

  bool do_vel = false;
  
  // check whether we "shake" solute
  if (topo.solute().distance_constraints().size() && 
      sim.param().shake.ntc > 1){
    DEBUG(8, "\twe need to shake SOLUTE");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	_lincs<math::vacuum>(topo, conf, sim, topo.solute().distance_constraints(),  
			     topo.solute().lincs(), parameter(), m_timing);
	break;
      case math::rectangular:
	_lincs<math::rectangular>(topo, conf, sim, topo.solute().distance_constraints(),  
				  topo.solute().lincs(), parameter(), m_timing);
	break;
      case math::triclinic:
	_lincs<math::triclinic>(topo, conf, sim, topo.solute().distance_constraints(),  
				topo.solute().lincs(), parameter(), m_timing);
	break;
      default:
	throw std::string("wrong boundary type");
    }
  }

  if (sim.param().system.nsm){
    DEBUG(8, "\twe need to shake SOLVENT");
    do_vel = true;
    switch(conf.boundary_type){
      case math::vacuum:
	_solvent<do_virial, math::vacuum>(topo, conf, sim, parameter(), m_solvent_timing);
	break;
      case math::rectangular:
	_solvent<do_virial, math::rectangular>(topo, conf, sim, parameter(), m_solvent_timing);
	break;
      case math::triclinic:
	_solvent<do_virial, math::triclinic>(topo, conf, sim, parameter(), m_solvent_timing);
	break;
      default:
	throw std::string("wrong boundary type");
    }
  }
  
  // "shake" velocities
  if (do_vel)
    conf.current().vel = (conf.current().pos - conf.old().pos) / 
      sim.time_step_size();

  // return success!
  return 0;
		   
}

static void _setup_lincs(topology::Topology const & topo,
			 topology::Compound::lincs_struct & lincs,
			 std::vector<topology::two_body_term_struct> const & constr,
			 size_t const offset = 0)
{
  const size_t num_constr = constr.size();  
  lincs.coupled_constr.resize(num_constr);
  lincs.coef.resize(num_constr);

  for(size_t i=0; i < num_constr; ++i){
    // diagonal matrix
    lincs.sdiag.push_back(1.0 / sqrt(1.0 / topo.mass()(constr[i].i + offset) +
				     1.0 / topo.mass()(constr[i].j + offset)));

    DEBUG(8, "diag[" << i << "] = " << lincs.sdiag[i]);

  }
  
  for(size_t i=0; i < num_constr; ++i){

    // connected constraints
    for(size_t j=i+1; j < num_constr; ++j){

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
	DEBUG(8, "constraint " << i << ": " << constr[i].i << " - " << constr[i].j);
	DEBUG(8, "constraint " << j << ": " << constr[j].i << " - " << constr[j].j);
	
	DEBUG(8, "connected: " << i << " - " << j << " with " << con);
	
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


template<math::virial_enum do_virial>
int algorithm::Lincs<do_virial>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim)
{

  // setup lincs
  _setup_lincs(topo, topo.solute().lincs(),
	       topo.solute().distance_constraints());
  
  // the first atom of a solvent
  size_t first = topo.num_solute_atoms();

  // for all solvents
  for(size_t i=0; i<topo.num_solvents(); ++i){
    
    _setup_lincs(topo, topo.solvent(i).lincs(),
		 topo.solvent(i).distance_constraints(),
		 first);

    first+=topo.solvent(i).num_atoms();
  }
  
  //================================================================================
  // initialize positions / velocities
  //================================================================================

  if (sim.param().start.shake_pos){
    std::cout << "shaking(lincs) initial positions\n";

    // old and current pos and vel are the same...
    // shake the current ones
    apply(topo, conf, sim);

    // restore the velocities
    conf.current().vel = conf.old().vel;
    
    // take a step back
    conf.old().pos = conf.current().pos;
    
    if (sim.param().start.shake_vel){
      std::cout << "shaking(lincs) initial velocities\n";

      conf.current().pos = conf.old().pos - 
	sim.time_step_size() * conf.old().vel;
    
      // shake again
      apply(topo, conf, sim);
    
      // restore the positions
      conf.current().pos = conf.old().pos;
    
      // velocities are in opposite direction (in time)
      conf.current().vel = -1.0 * conf.current().vel;
      conf.old().vel = conf.current().vel;
    }
    
  }
  else if (sim.param().start.shake_vel){
    io::messages.add("shaking velocities without shaking positions illegal.",
		     "lincs", io::message::error);
  }
  
  return 0;
}

template<math::virial_enum do_virial>
void algorithm::Lincs<do_virial>
::print_timing(std::ostream & os)
{
  os << std::setw(40) << std::left << "Lincs::solute"
     << std::setw(20) << m_timing << "\n"
     << std::setw(40) << std::left << "Lincs::solvent"
     << std::setw(20) << m_solvent_timing << "\n";
}
