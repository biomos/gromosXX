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

    math::VArray &pos = sim.system().pos();
    math::Matrix &pressure = sim.system().pressure();
    math::Box box = sim.system().periodicity().box();

    if (m_ntp == 1){ // isotropic
      // std::cout << "\tisotropic pressure coupling...\n";
      
      double total_pressure =  (pressure(0,0)
				+ pressure(1,1)
				+ pressure(2,2)) / 3.0;

      double mu = pow(1.0 - m_comp*dt/m_tau
		      * (m_pres0(0,0) - total_pressure), 1.0/3.0);


      // scale the box
      box = mu * box;
      sim.system().periodicity().box(box);

      // scale the positions
      for(int i=0; i<pos.size(); ++i)
	pos(i) = mu * pos(i);

    }
    else if(m_ntp == 2){ // anisotropic
      // std::cout << "\tanisotropic pressure coupling...\n";
    
      math::Vec mu;

      for(int i=0; i<3; ++i){
	mu(i) = pow(1.0 - m_comp*dt/m_tau
		    * (m_pres0(i,i) - pressure(i,i)), 1.0/3.0);
      }

      // scale the box
      for(int i=0; i<3; ++i)
	box(i) = box(i) * mu;
      sim.system().periodicity().box(box);

      // scale the positions
      for(int i=0; i<pos.size(); ++i)
	pos(i) = mu * pos(i);

    }
    else if(m_ntp == 3){ // fully anisotropic
      // std::cout << "\tfull anisotropic pressure coupling...\n";
      
      math::Matrix mu;

      for(int i=0; i<3; ++i){
	for(int j=0; j<3; ++i){
	  
	  mu(i, j) = pow(1.0 - m_comp*dt/m_tau
		      * (m_pres0(i,j) - pressure(i,j)), 1.0/3.0);
	}
      }

      // scale the box
      box = math::product(mu, box);
      sim.system().periodicity().box(box);

      // scale the positions
      for(int i=0; i<pos.size(); ++i)
	pos(i) = math::product(mu, pos(i));

    }
    
  }
  
}

inline void algorithm::Berendsen_Barostat
::initialize(int ntp, double pres0, double comp, double tau)
{
  if (ntp < 0 || ntp > 3){
    io::messages.add("Invalid pressure coupling scheme requested",
		     "Berendsen_Barostat",
		     io::message::error);
  }
    
  m_ntp = ntp;
  m_pres0(0,0) = pres0;
  m_pres0(1,1) = pres0;
  m_pres0(2,2) = pres0;
  m_comp = comp;
  m_tau = tau;
}

