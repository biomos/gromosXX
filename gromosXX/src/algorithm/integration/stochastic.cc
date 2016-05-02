/**
 * @file stochastic.cc
 * contains the implementation
 * of the stochastic dynamics algorithm
 */

#include "../../stdheader.h"

#include "../../algorithm/algorithm.h"
#include "../../topology/topology.h"
#include "../../simulation/simulation.h"
#include "../../configuration/configuration.h"

#include "../../simulation/multibath.h"
#include "../../simulation/parameter.h"

#include "../../math/periodicity.h"
#include "../../util/template_split.h"

#include "stochastic.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

 
 //as the friction coefficients only have to be calculated ones per cycle, let's do it for the velocities (as they are first in the cycle...)
template<math::boundary_enum B>
int algorithm::Stochastic_Dynamics_Vel1
::calc_friction_coeff(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation &sim)
{
  const unsigned int size = topo.num_atoms();
  math::Vec rij;
  std::vector<unsigned int> neigh(size, 0);
  
  math::Periodicity<B> periodicity(conf.current().box);

  double const cutoff2 = sim.param().stochastic.rcutf * sim.param().stochastic.rcutf;
  
  //first we count the number of neighbours if SASA weighting is applied...
 
  if (sim.param().stochastic.ntfr == 3) {
    
    // loop over all atoms and calculate neighbors
    for (unsigned int i = 0; i < size; ++i){
      for (unsigned int j = i+1; j < size; ++j){
        periodicity.nearest_image(conf.current().pos(i), conf.current().pos(j), rij);
	const double r2 = math::abs2(rij);

	if (r2 < cutoff2){
	  neigh[i] += 1;
          neigh[j] += 1;
	}
      }
    }
  }
  //determine the gammas
  
  for (unsigned int i = 0; i < size; ++i) {
    DEBUG(10, "neighbours(" << i << ") = " << neigh[i]);
    double xh = 0.0;
    switch(sim.param().stochastic.ntfr) {
      case 0 : // set gamma to 0.0
        xh = 0.0;
        break;
      case 1 : // set gamma to cfric
        xh = 1.0;
        break;
      case 2 : // set gamma to gamma0 x cfric
        xh = topo.stochastic().gamma(i);
        break;
      case 3 :
        xh = std::max(0.0, 1.0 - neigh[i]/double(sim.param().stochastic.nbref));
        break;
      default :
        io::messages.add("invalid NTFR in STOCHASTIC block",
		         "calc_fricition_coeff", io::message::error);
    }
    topo.stochastic().gamma(i) = sim.param().stochastic.cfric * xh;
    DEBUG(10, "gamma: " << topo.stochastic().gamma(i) << " xh=" << xh);
  }

  //we have the gammas...
  //...now generate SD coefficients
  for (unsigned int i = 0; i < size; ++i){
    DEBUG(12, "stochastic coefficient for atom i = " << i);
    const double gg = fabs(topo.stochastic().gamma(i));
    const double gdt = gg * sim.time_step_size();
    const double gdth = gdt * 0.5;

    const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));

    if (fabs(gdt) > 0.05){
      DEBUG(12, "\tdoing the analytical formulas");
      const double emdth = exp(-gdth);
      const double epdth = exp(+gdth);
      const double emdt  = emdth * emdth;
      const double epdt  = epdth * epdth;
      const double omdt  = 1.0 - emdt;
      const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
      const double ddth  = 2.0 - epdth - emdth;
      const double bpdth = gdt * (epdt - 1.0) - 4.0 * (epdth - 1.0) * (epdth - 1.0);
      const double bmdth = gdt * (1.0 - emdt) - 4.0 * (emdth -1.0) * (emdth -1.0);
      topo.stochastic().c1(i) = emdt;
      topo.stochastic().c2(i) = omdt / gdt;
      topo.stochastic().c3(i) = sqrt(fabs(omdt));
      topo.stochastic().c4(i) = sqrt(fabs(bpdth/cdth));
      topo.stochastic().c5(i) = gg * ddth/cdth;
      topo.stochastic().c6(i) = (epdth - emdth) / gdt;
      topo.stochastic().c7(i) = sqrt(fabs(cdth)) / gg;
      topo.stochastic().c8(i) = sqrt(fabs(bmdth/omdt)) / gg;
      topo.stochastic().c9(i) = -ddth/(gg * omdt);
      
      if (sim.steps() == 0 && sim.param().stochastic.generate_integral){
	DEBUG(10, "initial stochint: gamma(" << i << ") = " << topo.stochastic().gamma(i));
	const double sd = kToverM / topo.stochastic().gamma(i) * sqrt(cdth);
        m_rng->stddev(sd);
	conf.current().stochastic_integral(i) = m_rng->get_gaussian_vec();
      }
      
    } else {
      DEBUG(12, "\tdoing the power series");
      //this is a power series expansion for the coefficients used
      //in the SD algortihm. it should be used when gamdt < 0.05, otherwise
      //the analytical expression from above is good enough
      const double gdth2 = gdth * gdth;
      const double gdth3 = gdth2 * gdth;
      const double gdth4 = gdth2 * gdth2;
      const double gdth5 = gdth2 * gdth3;
      const double gdth6 = gdth3 * gdth3;
      const double gdth7 = gdth4 * gdth3;

      topo.stochastic().c1(i) = exp(-gdt);

      topo.stochastic().c2(i) = 1.0 - gdth + gdth2 * 2.0/3.0 - gdth3/3.0 
	+ gdth4 * 2.0/15.0 - gdth5 * 2.0/45.0 + gdth6 * 4.0/315.0;

      topo.stochastic().c3(i) = sqrt(fabs(topo.stochastic().c2(i)) * 2.0 * gdth);

      topo.stochastic().c4(i) = sqrt(fabs(gdth/2.0 + gdth2 * 7.0/8.0 + gdth3 * 367.0/480.0 + gdth4 * 857.0/1920.0 
					  + gdth5 * 52813.0/268800.0 + gdth6 * 224881.0/3225600.0 +
					  gdth7 * 1341523.0/64512000.0));

      topo.stochastic().c5(i) = -2.0 / sim.time_step_size() *
	(1.5 + gdth * 9.0/8.0 + gdth2 * 71.0/160.0 + gdth3 * 81.0/640.0 + gdth4 * 7807.0/268800.0 
	 + gdth5 * 1971.0/358400.0 + gdth6 * 56417.0/64512000.0);

      topo.stochastic().c6(i) = 1.0 + gdth2/6.0 + gdth4/10.0 + gdth6/5040.0;

      topo.stochastic().c7(i) = sim.time_step_size() * 0.5 * 
	sqrt(fabs(gdth * 2.0/3.0 - gdth2/2.0 + gdth3 * 7.0/30.0 -gdth4/12.0 + gdth5 * 31.0/1260.0 - 
		  gdth6/160.0 + gdth7 * 127.0/90720.0));

      topo.stochastic().c8(i) = sim.time_step_size() * 0.5 * 
	sqrt(fabs(gdth/6.0 - gdth3/60.0 + gdth5 * 17.0/10080.0 - gdth7 * 31.0/181440.0));

      topo.stochastic().c9(i) = sim.time_step_size() * 0.5 *
	(0.5 + gdth/2.0 + gdth2 * 5.0/24.0 + gdth3/24.0 + gdth4/240.0 + gdth5/720.0 + gdth6 * 5.0/8064.0);
 
      if (sim.steps() == 0 && sim.param().stochastic.generate_integral){
	// just an initial guess. no need(?) for high precision??
	DEBUG(10, "gdt = " << gdt << " in initial stochastic integral");
	if (topo.stochastic().gamma(i) < math::epsilon){
	  conf.current().stochastic_integral(i) = 0.0;
	}
	else{
	  const double emdth = exp(-gdth);
	  const double emdt  = emdth * emdth;
	  const double cdth  = gdt - 3.0 + 4.0 * emdth - emdt;
	  
	  const double sd = kToverM / topo.stochastic().gamma(i) * sqrt(cdth);
          m_rng->stddev(sd);
	  conf.current().stochastic_integral(i) = m_rng->get_gaussian_vec();
	}
      }
    }   
  }
  return 0;
}

/*SD Velocities 1
*/

algorithm::Stochastic_Dynamics_Vel1::Stochastic_Dynamics_Vel1(
const simulation::Parameter& param)
  : Algorithm("Stochastic_Dynamics_Vel1")
{
  m_rng = math::RandomGenerator::create(param, "0");
  m_rng->mean(0.0);    
}

/**
 * SD init
 */
int algorithm::Stochastic_Dynamics_Vel1
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream &os,
       bool quiet)
{
  if (!quiet)
    os << "STOCHASTIC DYNAMICS\n";
  
  if (!quiet && sim.param().stochastic.generate_integral)
    os << "\tgenerating initial stochastic integrals\n";
  
  if (sim.param().stochastic.ntfr == 2 && !quiet)
    os << "\tcalculating friction coefficients by GAM0 x CFRIC\n";
  
  if (sim.param().stochastic.ntfr == 3 && !quiet)
    os << "\tcalculating friction coefficients using FRIC\n";
  
  // set the seed
  if (sim.param().stochastic.generate_integral) {
    // use seed from input file
    std::ostringstream seed; seed << sim.param().start.ig;
    m_rng->seed(seed.str());
  } else {
    // use seed from configuration
    try {
      m_rng->seed(conf.current().stochastic_seed);
    } catch (std::runtime_error &exp) {
      io::messages.add(exp.what(), "Stochastic_Dynamics", io::message::error);
      return 1;
    }
  }
  
  // anyway call this, it will initialize the
  // friction coefficients with cfric otherwise...
  // and also initialise the stochastic integrals
  // if necessary
  SPLIT_BOUNDARY(calc_friction_coeff, topo, conf, sim);
  
#ifndef NDEBUG
  for(unsigned int i = 0; i < topo.num_atoms(); ++i) {
    DEBUG(13, "STOCHASTIC STRUCT for atom i = " << i);
    DEBUG(13, "--------------------------------");
    DEBUG(13, "C1: " << topo.stochastic().c1(i));
    DEBUG(13, "C2: " << topo.stochastic().c2(i));
    DEBUG(13, "C3: " << topo.stochastic().c3(i));
    DEBUG(13, "C4: " << topo.stochastic().c4(i));
    DEBUG(13, "C5: " << topo.stochastic().c5(i));
    DEBUG(13, "C6: " << topo.stochastic().c6(i));
    DEBUG(13, "C7: " << topo.stochastic().c7(i));
    DEBUG(13, "C8: " << topo.stochastic().c8(i));
    DEBUG(13, "C9: " << topo.stochastic().c9(i));
    DEBUG(13, "GAMMA: " << topo.stochastic().gamma(i));
  }
#endif
  
  if (!quiet)
    os << "\tseed: " << m_rng->seed() << "\n\n";
  
  if (sim.param().constraint.solute.algorithm != simulation::constr_off &&
      sim.param().constraint.solute.algorithm != simulation::constr_shake){
    io::messages.add("Constraints: only SHAKE or unconstrained supported",
		     "Stochastic_Dynamics",
		     io::message::error);
  }

  // this has been commented out to enable SD with SOLVENT simulation
  //if (topo.num_solvent_atoms() > 0){
  //  io::messages.add("SD only supported with no SOLVENT",
  //		     "Stochastic_Dynamics",
  //		     io::message::error);
  //}
  
  if (!quiet)
    os << "END\n";

  // resize the random vectors
  m_vrand.resize(topo.num_atoms());
  return 0;
}
/**
 * SD: calculate velocities 1
*/ 

int algorithm::Stochastic_Dynamics_Vel1
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
        simulation::Simulation &sim) {
  SPLIT_BOUNDARY(_apply, topo, conf, sim);
  return 0;
}

/**
 * SD: calculate velocities 1
*/ 
template<math::boundary_enum B>
int algorithm::Stochastic_Dynamics_Vel1
::_apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  DEBUG(7, "doing stochastic dynamics velocities")
  
  m_timer.start();

  // calculate atomic friction oefficients?
  if (sim.param().stochastic.ntfr == 3){
    if ((sim.steps() % sim.param().stochastic.nsfr) == 0)
      calc_friction_coeff<B>(topo, conf, sim);
  }

  conf.exchange_state();
  // copy the box
  conf.current().box = conf.old().box;     
  for (unsigned int i=0; i < topo.num_atoms(); ++i){
    if(topo.stochastic().gamma(i) != 0.0) {
	
      const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));

      //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
      double cf = sim.time_step_size() * topo.stochastic().c2(i) / topo.mass(i);
      
      //CC3(NATTOT) = delivered with SQRT(1-EXP(-GDT))
      //this is 2.11.2.8
      double sd1 = topo.stochastic().c3(i) * kToverM;
      
      //CC4(NATTOT) = delivered with SQRT(B(+GDT/2)/C(+GDT/2)) (SEE PAPER)
      // this is 2.11.2.9
      double sd2 = kToverM * topo.stochastic().c4(i);
      
      DEBUG(10, "\ti: " << i << " kT/m(i): " << kToverM <<
	    " sd1: " << sd1 << " sd2: " << sd2);
      
      //we sample the V'i vector from eq. 2.11.2.20 from a Gaussian with 0.0 mean
      //and width sigma2_sq (sd1)
      m_rng->stddev(sd1);
      m_vrand.vrand1(i) = m_rng->get_gaussian_vec();
      //we sample the Vi vector to be used in eq. 2.11.2.2 from a Gaussian with 0.0 mean
      //and width rho1_sq (sd2)
      m_rng->stddev(sd2);
      m_vrand.vrand2(i) = m_rng->get_gaussian_vec();
      
      DEBUG(10, "vrand1=" <<  math::v2s(m_vrand.vrand1(i)));
      DEBUG(10, "vrand2=" << math::v2s(m_vrand.vrand2(i)));
      
      //CC6(NATTOT) = delivered with (EXP(+GDT/2)-EXP(-GDT/2))/GDT
      //double cvdt = topo.stochastic().c6(i) * sim.time_step_size();
      
      DEBUG(9, "old vel(" << i << ") " << math::v2s(conf.old().vel(i)));
      DEBUG(9, "old pos(" << i << ") " << math::v2s(conf.old().pos(i)));
      DEBUG(9, "old frc(" << i << ") " << math::v2s(conf.old().force(i)));
      DEBUG(10, "c5=" << topo.stochastic().c5(i));
      DEBUG(10, "stochastic integral=" << math::v2s(conf.old().stochastic_integral(i)));
      
      math::Vec svh = conf.old().stochastic_integral(i) * topo.stochastic().c5(i) + m_vrand.vrand2(i);
      //assign vrand to the Stochastic Integral array
      conf.current().stochastic_integral(i) = m_vrand.vrand1(i);
      
      //calculate new velocity using eq. 2.11.2.2 from the GROMOS96 book
      //CC1(NATTOT) = delivered with EXP(-GDT)  (GDT=GAM(J)*DT)
      //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
      //slightly rewritten form of first line of 2.11.2.3 minus
      //the first term in the third line of 2.11.2.3
      
      DEBUG(10, "svh=" << math::v2s(svh));
      DEBUG(10, "cf=" << cf);
      
      conf.current().vel(i) = (conf.old().vel(i) - svh) * topo.stochastic().c1(i)
	// force * m-1 * dt * (1-EXP(-GDT))/GDT, i.e.
	+ conf.old().force(i) * cf
	// last term of 2.11.2.2
	+ m_vrand.vrand1(i);
      
      const double sd3 = kToverM * topo.stochastic().c7(i);
      const double sd4 = kToverM * topo.stochastic().c8(i);
      
      DEBUG(10, "sd3: " << sd3 << " sd4: " << sd4);   
    
      //we sample the R'i vector from eq. 2.11.2.25 from a Gaussian with 0.0 mean
      //and width rho2_sq (sd3)
      m_rng->stddev(sd3);
      m_vrand.vrand3(i) = m_rng->get_gaussian_vec();
      //we sample the Ri vector to be used in eq. 2.11.2.26 from a Gaussian with 0.0 mean
      //and width rho1_sq (sd4)
      m_rng->stddev(sd4);
      m_vrand.vrand4(i) = m_rng->get_gaussian_vec();
      
      DEBUG(10, "vrand3=" << math::v2s(m_vrand.vrand3(i)));
      DEBUG(10, "vrand4=" << math::v2s(m_vrand.vrand4(i)));
    } else {
      conf.current().vel(i) = conf.old().vel(i) + conf.old().force(i) * sim.time_step_size() / topo.mass(i);
      m_vrand.vrand1(i)=math::Vec(0.0, 0.0, 0.0);
      m_vrand.vrand2(i)=math::Vec(0.0, 0.0, 0.0);
      m_vrand.vrand3(i)=math::Vec(0.0, 0.0, 0.0);
      m_vrand.vrand4(i)=math::Vec(0.0, 0.0, 0.0);

    }
  } // loop over atoms
 
  // save the seed
  conf.current().stochastic_seed = m_rng->seed();

  m_timer.stop();
  return 0;
}


/**
 * SD: calculate positions 1
*/ 
int algorithm::Stochastic_Dynamics_Pos1
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  DEBUG(7, "doing stochastic dynamics positions")  
  m_timer.start(); 
  for (unsigned int i=0; i < topo.num_atoms(); ++i){
    //calc new positions
    //according to step 7 in leap frog for SD (GROMOS96 book)
    conf.current().pos(i) = conf.old().pos(i) 
    			+ conf.current().vel(i) * sim.time_step_size() * topo.stochastic().c6(i);    
    DEBUG(9, "old pos(" << i << ") " << math::v2s(conf.old().pos(i)));
    DEBUG(9, "new pos(" << i << ") " << math::v2s(conf.current().pos(i)));      
  } // loop over atoms
  m_timer.stop();
  return 0;
}


/*SD Velocities 2
*/

int algorithm::Stochastic_Dynamics_Vel2
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
        simulation::Simulation &sim) {
  SPLIT_BOUNDARY(_apply, topo, conf, sim);
  return 0;
}

template<math::boundary_enum B>
int algorithm::Stochastic_Dynamics_Vel2
::_apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  DEBUG(7, "doing stochastic dynamics velocities (from SHAKE)")
  
  m_timer.start();
  //now comes the second part, after calling shake
  //as given by 2.11.2.24
  
  math::Periodicity<B> periodicity(conf.current().box);

  if (sim.param().constraint.solute.algorithm == simulation::constr_shake){ 
    // --- call SHAKE before ! ---
    // (without velocity correction)

    // do the SD-velocity correction right away...
    // velocity correction, if constraints are turned on
    // this is 2.11.2.24
    for (unsigned int i=0; i < topo.num_atoms(); ++i) {
      double cinv = 1.0 / (topo.stochastic().c6(i) * sim.time_step_size());
      math::Vec r;
      periodicity.nearest_image(conf.current().pos(i), conf.old().pos(i), r);
      conf.current().vel(i) = r * cinv;
      DEBUG(10, "velocity SHAKEN" << math::v2s(conf.current().vel(i)))
    }
  } // constraints
  else{
    if (sim.param().constraint.solute.algorithm != simulation::constr_off){
      io::messages.add("Constraints: only SHAKE or unconstrained supported",
		       "SD integration",
		       io::message::error);
    }
  }

  m_timer.stop();
  return 0;
}

/**
 * SD: calculate positions 2
*/ 

int algorithm::Stochastic_Dynamics_Pos2
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  DEBUG(7, "doing stochastic dynamics positions")
  
  m_timer.start();

  for (unsigned int i=0; i < topo.num_atoms(); ++i){
    //this is 2.11.2.25
    if(topo.stochastic().gamma(i) != 0.0) {
      math::Vec sxh = conf.current().stochastic_integral(i) * topo.stochastic().c9(i)
      + m_vrand->vrand4(i);
      conf.current().stochastic_integral(i) = m_vrand->vrand3(i);  
      //this is 2.11.2.26
      conf.current().pos(i) += m_vrand->vrand3(i) - sxh;
    }     
  } // loop over atoms

  m_timer.stop();
  return 0;
}
