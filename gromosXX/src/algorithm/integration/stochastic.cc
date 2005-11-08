/**
 * @file stochastic.cc
 * contains the implementation
 * of the stochastic dynamics algorithm
 */

#include <stdheader.h>

#include <algorithm/algorithm.h>
#include <topology/topology.h>
#include <simulation/simulation.h>
#include <configuration/configuration.h>

#include "stochastic.h"

#undef MODULE
#undef SUBMODULE
#define MODULE algorithm
#define SUBMODULE integration

algorithm::Stochastic_Dynamics_Pos::Stochastic_Dynamics_Pos()
  : Algorithm("Stochastic_Dynamics_Pos")
{

  // get a random number generator type
  const gsl_rng_type *rng_type;
	
  // enable control via environment variables
  gsl_rng_env_setup();
  rng_type = gsl_rng_default;

  // get the rundom number generator
  m_rng = gsl_rng_alloc(rng_type);

}

/**
 * SD init
 */
int algorithm::Stochastic_Dynamics_Pos
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream &os,
       bool quiet)
{
  os << "Stochastic Dynamics\n";
  os << "\tbloody bloody bla...\n";      
  

  if (sim.param().stochastic.ntfr == 3)
    os << "\tcalculating friction coefficients\n";
  // anyway call this, it will initialize the
  // friction coefficients with cfric otherwise...
  // and also initialise the stochastic integrals
  // if necessary
  calc_friction_coeff(topo, conf, sim);

  os << "\t\trandom number generator: " << gsl_rng_name(m_rng);

  // set seed if not set by environment variable
  if (gsl_rng_default_seed == 0){
    gsl_rng_set(m_rng, sim.param().start.ig);
    os << "\n\t\tseed: " << sim.param().start.ig << "\n\n";
  }
  else{
    os << "\n\t\t(default) seed: " << gsl_rng_default_seed << "\n\n";
  }
  
  if (sim.param().constraint.solute.algorithm != simulation::constr_off &&
      sim.param().constraint.solute.algorithm != simulation::constr_shake){
    io::messages.add("Constraints: only SHAKE or unconstrained supported",
		     "SD integration",
		     io::message::error);
  }

  if (topo.num_solvent_atoms() > 0){
    io::messages.add("SD only supported with no SOLVENT",
		     "SD integration",
		     io::message::error);
  }

  return 0;
}

/**
 * SD: calculate positions
 */
int algorithm::Stochastic_Dynamics_Pos
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation &sim)
{
  const double start = util::now();

  if (sim.param().stochastic.ntfr == 3){
    if ((sim.steps() % sim.param().stochastic.nsfr) == 0)
      calc_friction_coeff(topo, conf, sim);
  }

  conf.exchange_state();
  // copy the box
  conf.current().box = conf.old().box;

  double oldvel = 0, svh = 0, cf = 0, sd1 = 0, sd2 = 0, cvdt = 0;
     
  for (unsigned int i=0; i < topo.num_atoms(); ++i){
    
    const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));

    //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
    cf = sim.time_step_size() * topo.stochastic().c2(i) / topo.mass(i);

    //CC3(NATTOT) = delivered with SQRT(1-EXP(-GDT))
    //this is 2.11.2.8
    sd1 = topo.stochastic().c3(i) * kToverM;

    //CC4(NATTOT) = delivered with SQRT(B(+GDT/2)/C(+GDT/2)) (SEE PAPER)
    // this is 2.11.2.9
    sd2 = kToverM * topo.stochastic().c4(i);

    //CC5 = GAM(J)*D(+GDT/2)/C(+GDT/2)
    //this is 2.11.2.10
    //csx = cc5[i];
    //System.out.println("SDS: " + cc1[i] + " " + cf + " " + sd1 + " " + sd2 + " " + csx + " " + kToverM[i]);
      
    math::Vec vrand, vvrand;
    
    //dont really know what this thing will eventually do

    //we sample the V'i vector from eq. 2.11.2.20 from a Gaussian with 0.0 mean
    //and width sigma2_sq (sd1)
    for(int d=0; d<3; ++d){
      vrand(d) = gsl_ran_gaussian(m_rng, sd1);
    }
    //we sample the Vi vector to be used in eq. 2.11.2.2 from a Gaussian with 0.0 mean
    //and width rho1_sq (sd2)
    for(int d=0; d<3; ++d){
      vvrand(d) = gsl_ran_gaussian(m_rng, sd2);
    }

    //CC6(NATTOT) = delivered with (EXP(+GDT/2)-EXP(-GDT/2))/GDT
    cvdt = topo.stochastic().c6(i) * sim.time_step_size();

    DEBUG(9, "old vel(" << i << ") " << math::v2s(conf.old().vel(i)));
    DEBUG(9, "old pos(" << i << ") " << math::v2s(conf.old().pos(i)));
    DEBUG(9, "old frc(" << i << ") " << math::v2s(conf.old().force(i)));

    for (int j=0; j < 3; ++j) {
      // this corresponds to eq. 2.11.2.20 from the GROMOS96 book
      // all ccX coefficients are calculated in the FrictionCoeff routines
      DEBUG(10, "stochastic integral=" << conf.old().stochastic_integral(i)(j));
      DEBUG(10, "c5=" << topo.stochastic().c5(i));
      DEBUG(10, "vvrand=" <<  vvrand(j) << " vrand=" << vrand(j));
      svh = conf.old().stochastic_integral(i)(j) * topo.stochastic().c5(i) + vvrand(j);
	      
      //assign vrand to the Stochastic Integral array
      conf.current().stochastic_integral(i)(j) = vrand(j);

      oldvel = conf.old().vel(i)(j);

      //calculate new velocity using eq. 2.11.2.2 from the GROMOS96 book
      //CC1(NATTOT) = delivered with EXP(-GDT)  (GDT=GAM(J)*DT)
      //CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
      //slightly rewritten form of first line of 2.11.2.3 minus 
      //the first term in the third line of 2.11.2.3
      //as usual with GROMOS96 documentation standards, there is no comment in the code...          
      
      DEBUG(10, "svh=" << svh << " c1=" << topo.stochastic().c1(i));
      DEBUG(10, "cf=" << cf);

      conf.current().vel(i)(j) = (oldvel-svh) * topo.stochastic().c1(i) +
	// force * m-1 * dt * (1-EXP(-GDT))/GDT, i.e.
	+ conf.old().force(i)(j) * cf
	// last term of 2.11.2.2
	+ vrand(j);
      
      //calc new positions
      //this is 2.11.2.21
      DEBUG(10, "cvdt=" << cvdt);
      conf.current().pos(i)(j) = conf.old().pos(i)(j) + conf.current().vel(i)(j) * cvdt;

    } // loop over dimenstions

    DEBUG(9, "new vel(" << i << ") " << math::v2s(conf.current().vel(i)));
    DEBUG(9, "new pos(" << i << ") " << math::v2s(conf.current().pos(i)));
    
  } // loop over atoms

  m_timing += util::now() - start;
  return 0;
}

/**
 * friction coefficients
 *   calculates the friction coefficients from
 *   GAM(J) = CFRIC * MAX(0;1-NB(J)/NBREF)
 *   where NB(J) denotes the number of neighbour atoms of atom J
 *   within a cut-off distance RCUTF. This values are then
 *   used in calculating the CC# arrays.
 *   some variables that are used:
 *   CC1(NATTOT) = delivered with EXP(-GDT)  (GDT=GAM(J)*DT)
 *   CC2(NATTOT) = delivered with (1-EXP(-GDT))/GDT
 *   CC3(NATTOT) = delivered with SQRT(1-EXP(-GDT))
 *   CC4(NATTOT) = delivered with SQRT(B(+GDT/2)/C(+GDT/2)) (SEE PAPER)
 *   CC5(NATTOT) = delivered with GAM(J)*D(+GDT/2)/C(+GDT/2)
 *   CC6(NATTOT) = delivered with (EXP(+GDT/2)-EXP(-GDT/2))/GDT
 *   CC7(NATTOT) = delivered with SQRT(C(+GDT/2))/GAM(J)
 *   CC8(NATTOT) = delivered with SQRT(B(-GDT/2)/(1-EXP(-GDT)))/GAM(J)
 *   CC9(NATTOT) = delivered with D(+GDT/2)/(EXP(-GDT)-1)/GAM(J)
 *
 */
int algorithm::Stochastic_Dynamics_Pos
::calc_friction_coeff(topology::Topology & topo,
		      configuration::Configuration & conf,
		      simulation::Simulation &sim)
{
  // copied from mika,
  // mistakes are mine (at least some of them)

  const int size = topo.num_atoms();
  math::Vec rij;
  std::vector<int> neigh(size, 0);

  double const cutoff2 = sim.param().stochastic.rcutf * sim.param().stochastic.rcutf;
  
  //first we count the number of neighbours if SASA weighting is applied...
 
  if (conf.boundary_type != math::vacuum){
    io::messages.add("SD only implemented for vacuum!",
		     "Stochastic_Dynamics",
		     io::message::critical);
    return 1;
  }
  
  if (sim.param().stochastic.ntfr == 3){

    for (int i = 0; i < size; ++i){
      for (int j = i+1; j < size; ++j){

	rij = conf.current().pos(i) - conf.current().pos(j);
	const double r2 = math::abs2(rij);

	if (r2 < cutoff2){
	  neigh[i] += 1;
          neigh[j] += 1;
	}
      }
    }
  }

  //dont know what this rbref really means yet
  //determine the gammas
  topo.stochastic().resize(size);
  
  double xh = 1.0;
  for (int i = 0; i < size; ++i) {
    DEBUG(10, "neighbours(" << i << ") = " << neigh[i]);
    if (sim.param().stochastic.ntfr == 3) xh = std::max(0.0, 1.0 - neigh[i]/double(sim.param().stochastic.nbref));
    topo.stochastic().gamma(i) = sim.param().stochastic.cfric * xh;
    DEBUG(10, "gamma: " << topo.stochastic().gamma(i) << " xh=" << xh);
  }

  //we have the gammas...
  //...now generate SD coefficients
  for (int i = 0; i < size; ++i){
    const double gg = fabs(topo.stochastic().gamma(i));
    const double gdt = gg * sim.time_step_size();
    const double gdth = gdt * 0.5;

    const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));

    if (fabs(gdt) > 0.05){
      
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
	for(int d=0; d<3; ++d){
	  conf.current().stochastic_integral(i)(d) = gsl_ran_gaussian(m_rng, sd);
	}
      }
      
    }
    else {
      //this crap is a power series expansion for the coefficients used
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
	  for(int d=0; d<3; ++d){
	    conf.current().stochastic_integral(i)(d) = gsl_ran_gaussian(m_rng, sd);
	  }
	}
      }
    }
      
  }
  return 0;
}

int algorithm::Stochastic_Dynamics_Int
::init(topology::Topology & topo,
       configuration::Configuration & conf,
       simulation::Simulation &sim,
       std::ostream &os,
       bool quiet)
{
  return 0;
}

/**
 * SD: add stochastic integrals
 */
int algorithm::Stochastic_Dynamics_Int
::apply(topology::Topology & topo,
	configuration::Configuration & conf,
	simulation::Simulation & sim)
{
  const double start = util::now();

  double cinv = 0, sd3 = 0, sd4 = 0, csv = 0, sxh = 0;
  //now comes the second part, after calling shake
  //as given by 2.11.2.24

  if (sim.param().constraint.solute.algorithm == simulation::constr_shake){ 
    // --- call SHAKE before ! ---
    // System.out.println("SHAKE!");
    // cons.docorrection();

    // do the SD-velocity correction right away...
    // velocity correction, if constraints are turned on
    // this is 2.11.2.24
    for (unsigned int i=0; i < topo.num_atoms(); ++i) {

      cinv = 1.0 / (topo.stochastic().c6(i) * sim.time_step_size());
      
      for (int j=0; j < 3; ++j) {
	conf.current().vel(i)(j) = (conf.current().pos(i)(j) - conf.old().pos(i)(j)) * cinv;
      }
    }
  } // constraints
  else{
    if (sim.param().constraint.solute.algorithm != simulation::constr_off){
      io::messages.add("Constraints: only SHAKE or unconstrained supported",
		       "SD integration",
		       io::message::error);
    }
  }

  for (unsigned int i=0; i < topo.num_atoms(); ++i){

    const double kToverM = sqrt(math::k_Boltzmann * sim.param().stochastic.temp / topo.mass(i));
    
    //CVINV = 1.0/(CC6(I)*DT)
    cinv = 1.0 / (topo.stochastic().c6(i) * sim.time_step_size());
    sd3 = kToverM * topo.stochastic().c7(i);
    sd4 = kToverM * topo.stochastic().c8(i);

    //CSV   = CC9(I)
    //CC9(NATTOT) = delivered with D(+GDT/2)/(EXP(-GDT)-1)/GAM(J)
    csv = topo.stochastic().c9(i);

    math::Vec rrand, rrrand;

    //we sample the R'i vector from eq. 2.11.2.25 from a Gaussian with 0.0 mean
    //and width rho2_sq (sd3)
    for(int d=0; d<3; ++d){
      rrand(d) = gsl_ran_gaussian(m_rng, sd3);
    }
    //we sample the Ri vector to be used in eq. 2.11.2.26 from a Gaussian with 0.0 mean
    //and width sigma1_sq (sd4)
    for(int d=0; d<3; ++d){
      rrrand(d) = gsl_ran_gaussian(m_rng, sd4);
    }
         
    //loop over dimensions...
    for (int j=0; j < 3; ++j){
      
      //this is 2.11.2.25
      sxh = conf.current().stochastic_integral(i)(j) * csv + rrrand(j);

      //this is 2.11.2.26
      conf.current().pos(i)(j) += (rrand(j) - sxh);
      conf.current().stochastic_integral(i)(j) = rrand(j);
    }
  } // loop over atoms
         
  // --- call SHAKE again ---
  // if (constrain) cons.docorrection();

  m_timing += util::now() - start;
  return 0;
}
