/**
 * @file rottrans.cc
 * contains the template methods for
 * the class Rottrans_Constraints.
 */

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE constraints

/**
 * Constructor.
 */
template<math::virial_enum do_virial>
algorithm::Rottrans_Constraints<do_virial>
::Rottrans_Constraints(std::string const name)
  : Algorithm(name)
{
}

/**
 * Destructor.
 */
template<math::virial_enum do_virial>
algorithm::Rottrans_Constraints<do_virial>
::~Rottrans_Constraints()
{
}

template<math::boundary_enum b, math::virial_enum do_virial>
static int _apply(topology::Topology & topo,
		  configuration::Configuration & conf,
		  simulation::Simulation & sim)
{
  // calculate c
  DEBUG(8, "Roto-constrational translaints");
  
  blitz::TinyVector<double, 6U> c(0.0), lambda(0.0);
  
  for(unsigned int i=0; i<topo.num_solute_atoms(); ++i){
    const math::Vec diff = conf.current().pos(i) - conf.old().pos(i);
    
    c(0) += topo.mass()(i) * diff(0);
    c(1) += topo.mass()(i) * diff(1);
    c(2) += topo.mass()(i) * diff(2);

    c(3) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(1) * diff(2) -
			      conf.special().rottrans_constr.pos(i)(2) * diff(1));

    c(4) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(2) * diff(0) -
			      conf.special().rottrans_constr.pos(i)(0) * diff(2));

    c(5) += topo.mass()(i) * (conf.special().rottrans_constr.pos(i)(0) * diff(1) -
			      conf.special().rottrans_constr.pos(i)(1) * diff(0));
  }

  // lambda = - theta_inv * c
  for(int d1=0; d1 < 6; ++ d1){
    for(int d2=0; d2 < 6; ++d2){
      lambda(d1) -= conf.special().rottrans_constr.theta_inv(d1, d2) * c(d2);
    }
    DEBUG(10, "\tlambda[" << d1 << "] = " << lambda(d1) << "\n");
  }
  
  // update the positions
  for(unsigned int i=0; i<topo.num_solute_atoms(); ++i){

    conf.current().pos(i)(0) += // sim.time_step_size() * sim.time_step_size() *
      ( lambda(0) +
	lambda(4) * conf.special().rottrans_constr.pos(i)(2) -
	lambda(5) * conf.special().rottrans_constr.pos(i)(1)
	);

    conf.current().pos(i)(1) += // sim.time_step_size() * sim.time_step_size() *
      ( lambda(1) -
	lambda(3) * conf.special().rottrans_constr.pos(i)(2) +
	lambda(5) * conf.special().rottrans_constr.pos(i)(0)
	);
    
    conf.current().pos(i)(2) += // sim.time_step_size() * sim.time_step_size() *
      ( lambda(2) +
	lambda(3) * conf.special().rottrans_constr.pos(i)(1) -
	lambda(4) * conf.special().rottrans_constr.pos(i)(0)
	);
  }

  // update velocities
  conf.current().vel = (conf.current().pos - conf.old().pos) / 
    sim.time_step_size();

  //==================================================
  // CHECK
  //==================================================
  
  math::Vec v = 0.0;
  for(unsigned int i=0; i<topo.num_solute_atoms(); ++i){
    v += topo.mass()(i) * (math::cross(conf.special().rottrans_constr.pos(i), conf.current().pos(i)));
  }
  
  DEBUG(10, "testestest: " << math::v2s(v));

  return 0;
}

/**
 * apply roto-translational constraints
 */
template<math::virial_enum do_virial>
int algorithm::Rottrans_Constraints<do_virial>
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  DEBUG(7, "roto-constrational translaints");

  // the roto-constrational translaints
  if (sim.param().rottrans.rottrans){
    
    DEBUG(8, "\tapplying roto-constrational translaints");
    
    switch(conf.boundary_type){
      case math::vacuum:
	_apply<math::vacuum, do_virial>(topo, conf, sim);
	break;
      case math::rectangular:
	_apply<math::rectangular, do_virial>(topo, conf, sim);
	break;
      case math::triclinic:
	_apply<math::triclinic, do_virial>(topo, conf, sim);
	break;
      default:
	throw std::string("wrong boundary type");
    }
  }
  
  return 0;		   
}


template<math::boundary_enum b, math::virial_enum do_virial>
static int _init(topology::Topology & topo,
		 configuration::Configuration & conf,
		 simulation::Simulation & sim,
		 bool quiet)
{
  if (!quiet)
    std::cout << "Roto-translational constraints\tON\n";

  math::Periodicity<b> periodicity(conf.current().box);
  configuration::State_Properties sp(conf);

  math::Vec com_pos;
  math::Matrix com_e_kin;
  
  periodicity.put_chargegroups_into_box(conf, topo);
  sp.center_of_mass(0, topo.num_solute_atoms(), topo.mass(), com_pos, com_e_kin);

  DEBUG(12, "com: " << math::v2s(com_pos));

  conf.current().pos -= com_pos;
  periodicity.put_chargegroups_into_box(conf, topo);
  conf.exchange_state();
  conf.current().pos -= com_pos;
  periodicity.put_chargegroups_into_box(conf, topo);
  conf.exchange_state();

  // check
  sp.center_of_mass(0, topo.num_solute_atoms(), topo.mass(), com_pos, com_e_kin);
  DEBUG(10, "com (centered): " << math::v2s(com_pos));

  // store initial (reference) positions
  conf.special().rottrans_constr.pos.resize(topo.num_solute_atoms());
  for(unsigned int i=0; i < topo.num_solute_atoms(); ++i){
    conf.special().rottrans_constr.pos(i) = conf.current().pos(i);
  }
  
  // now calculate the thetas
  std::vector<math::Vec> theta(6, math::Vec(0.0));
  const int X = 0;
  const int Y = 1;
  const int Z = 2;
  
  for(unsigned int i=0; i < topo.num_solute_atoms(); ++i){
    theta[0](X) += topo.mass()(i);
    theta[1](Y) += topo.mass()(i);
    theta[2](Z) += topo.mass()(i);
    
    theta[3](X) += topo.mass()(i) * (conf.current().pos(i)(Z) * conf.current().pos(i)(Z) +
				     conf.current().pos(i)(Y) * conf.current().pos(i)(Y));
    
    theta[3](Y) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(X));
    
    theta[3](Z) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Z));
    
    
    theta[4](X) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Y));
    
    theta[4](Y) += topo.mass()(i) * (conf.current().pos(i)(X) * conf.current().pos(i)(X) +
				     conf.current().pos(i)(Z) * conf.current().pos(i)(Z));
    
    theta[4](Z) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(Z));
    

    theta[5](X) += topo.mass()(i) * (-conf.current().pos(i)(X) * conf.current().pos(i)(Z));
    
    
    theta[5](Y) += topo.mass()(i) * (-conf.current().pos(i)(Y) * conf.current().pos(i)(Z));
    
    theta[5](Z) += topo.mass()(i) * (conf.current().pos(i)(X) * conf.current().pos(i)(X) +
				     conf.current().pos(i)(Y) * conf.current().pos(i)(Y));
    
  }

  DEBUG(10, "theta[4] " << math::v2s(theta[3]));
  DEBUG(10, "theta[5] " << math::v2s(theta[4]));
  DEBUG(10, "theta[6] " << math::v2s(theta[5]));

  const double d = math::dot(math::cross(theta[3], theta[4]), theta[5]);

  // inverse that
  blitz::TinyMatrix<double, 6U, 6U> & theta_inv = conf.special().rottrans_constr.theta_inv;

  for(int d1=0; d1<6; ++d1)
    for(int d2=0; d2<6; ++d2)
      theta_inv(d1, d2) = 0.0;
  
  theta_inv(0, 0) = theta_inv(1, 1) = theta_inv(2, 2) = 1.0 / theta[0](0);
  
  math::Vec h = math::cross(theta[4], theta[5]) / d;
  theta_inv(3, 3) = h(0);
  theta_inv(4, 3) = h(1);
  theta_inv(5, 3) = h(2);
  
  h = math::cross(theta[3], theta[5]) / (-d);
  
  theta_inv(3, 4) = h(0);
  theta_inv(4, 4) = h(1);
  theta_inv(5, 4) = h(2);
  
  h = math::cross(theta[3], theta[4]) / d;
  
  theta_inv(3, 5) = h(0);
  theta_inv(4, 5) = h(1);
  theta_inv(5, 5) = h(2);

  return 0;
}

template<math::virial_enum do_virial>
int algorithm::Rottrans_Constraints<do_virial>
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation & sim,
       bool quiet)
{
  switch(conf.boundary_type){
    case math::vacuum:
      _init<math::vacuum, do_virial>(topo, conf, sim, quiet);
      break;
    case math::rectangular:
      _init<math::rectangular, do_virial>(topo, conf, sim, quiet);
      break;
    case math::triclinic:
      _init<math::triclinic, do_virial>(topo, conf, sim, quiet);
      break;
    default:
      throw std::string("wrong boundary type");
  }
  
  return 0;
}
  

