/**
 * @file pressure/berendsen.tcc
 * methods of the berendsen barostat.
 */

inline algorithm::Berendsen_Barostat::Berendsen_Barostat()
{
}

template<typename t_simulation>
inline void algorithm::Berendsen_Barostat
::apply(t_simulation &sim, double const dt)
{
  // calculate the pressure (tensor)
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      sim.system().pressure()(i,j) = 2 * (sim.system().molecular_kinetic_energy()(i,j)
					  - sim.system().virial()(i,j)) /
	sim.system().periodicity().volume();

  if (m_ntp){

    math::Matrix mu;

    if (m_ntp == 1){
      std::cout << "\tisotropic pressure coupling...\n";
      
      double total_pressure =  (sim.system().pressure()(0,0)
				+ sim.system().pressure()(1,1)
				+ sim.system().pressure()(2,2)) / 3.0;

      double mu_iso = pow(1.0 - m_comp*dt/m_tau
			  * (m_pres0 - total_pressure), 1.0/3.0);


      for(int i=0; i<3; ++i){
	for(int j=0; j<3; ++j){
	  if (i != j) mu(i,j) = 0.0;
	  else
	    mu(i,j) = mu_iso;
	}
      }
    }
    else if(m_ntp == 2){
    std::cout << "\tanisotropic pressure coupling...\n";

      for(int i=0; i<3; ++i){
	for(int j=0; j<3; ++j){
	  if (i != j) mu(i,j) = 0.0;
	  else
	    mu(i,j) = pow(1.0 - m_comp*dt/m_tau
			  * (m_pres0 - sim.system().pressure()(i,j)), 1.0/3.0);
	}
      }
    }
      

    // std::cout.precision(20);
    // std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);

    // std::cout << "comp*dt/tau = " << m_comp*dt/m_tau << std::endl;
    // std::cout << "pres0=" << m_pres0 << " presto=" << total_pressure << std::endl;
    // std::cout << "mu: " << mu << std::endl;

    // scale the box
    math::Box box = sim.system().periodicity().box();
    box = math::product(mu, sim.system().periodicity().box());
    sim.system().periodicity().box(box);

    // scale the positions
    math::VArray &pos = sim.system().pos();
    for(int i=0; i<pos.size(); ++i)
      pos(i) = math::product(mu, pos(i));

  }
  
}

inline void algorithm::Berendsen_Barostat
::initialize(int ntp, double pres0, double comp, double tau)
{
  m_ntp = ntp;
  m_pres0 = pres0;
  m_comp = comp;
  m_tau = tau;
}


