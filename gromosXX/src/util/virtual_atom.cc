/**
 * 
 * 
 */

#include <stdheader.h>
#include "virtual_atom.h"
#include <util/debug.h>

#undef MODULE
#undef SUBMODULE
#define MODULE util
#define SUBMODULE util

static const double TETHCO=0.577350;
static const double TETHSI=0.816497;




util::Virtual_Atom::Virtual_Atom()
  :  m_type(0),
     m_atom(),
     m_dish(0.1),
     m_disc(0.153),
     m_orientation(0)
{
}

util::Virtual_Atom::Virtual_Atom(int type, std::vector<int> atom,
				 double dish, double disc,int orientation)
  :  m_type(type),
     m_atom(atom),
     m_dish(dish),
     m_disc(disc),
     m_orientation(orientation)
{
}

math::Vec util::Virtual_Atom::pos(math::VArray const & position)const

{
  math::Vec s,t;

  switch(m_type){
    
    case 0: // explicit atom
    case 7: // rotating ring
      assert(m_atom.size()>0);
      return position(m_atom[0]); 

    case 1: // CH1
    
      assert(m_atom.size()>3);
      s = 3.0 * position(m_atom[0]) - position(m_atom[1])
	- position(m_atom[2]) - position(m_atom[3]);
      return  position(m_atom[0]) + m_dish / math::abs(s) * s;

    case 2: // aromatic H
      assert(m_atom.size()>2);
       
      s = 2.0 * position(m_atom[0]) - position(m_atom[1])
	- position(m_atom[2]);
      return  position(m_atom[0]) + m_dish / math::abs(s) * s;
    
    case 3: // non-stereospecific CH2
      assert(m_atom.size()>2);
      s = 2.0 * position(m_atom[0]) - position(m_atom[1])
	- position(m_atom[2]);
      return  position(m_atom[0]) + m_dish * TETHCO / math::abs(s) * s;
       
    case 4: // stereospecific CH2
      assert(m_atom.size()>2);
      
      s = 2.0 * position(m_atom[0]) - position(m_atom[1]) - position(m_atom[2]);
      DEBUG(10, "\ts = " << math::v2s(s));
      
      t = math::cross(position(m_atom[0]) - position(m_atom[1]), position(m_atom[0]) - position(m_atom[2]));
      DEBUG(10, "\tq = " << math::v2s(t));
      
      DEBUG(10, "\tDISH = " << m_dish << "\tTETHCO = " << TETHCO << "\tTETHSI = " << TETHSI);
      return position(m_atom[0]) + m_dish * TETHCO / math::abs(s) * s + m_dish * TETHSI / math::abs(t)*t;
 
      break;
      
    case 5: // CH3
      assert(m_atom.size()>1);
      
      s =  position(m_atom[0]) - position(m_atom[1]);
      return position(m_atom[0]) + m_dish / (3 * math::abs(s)) * s;
      
    case 6: // non-stereospecific CH3 (Leu, Val)
      assert(m_atom.size()>2);
      
      s = 2.0 * position(m_atom[0]) - position(m_atom[1])
	- position(m_atom[2]);
      return position(m_atom[0]) - TETHCO * (m_disc + m_dish / 3.0) / math::abs(s) * s;
      
    case 8: // NH2-group (one pseudosite)
      assert(m_atom.size()>1);
      
      s = 2.0 * position(m_atom[0]) -position(m_atom[1])
	- position(m_atom[2]);
      return position(m_atom[0]) - (m_dish * 0.5) * s / math::abs(s);
      
    case 9: // (CH3)3-group (one psuedosite)
      assert(m_atom.size()>1);
      
      s = position(m_atom[0]) -position(m_atom[1]);
      return position(m_atom[0]) +  ( m_disc + m_dish/3.0 )/ (3 * math::abs(s)) * s;
      
    default:
      io::messages.add("Virual Atom", "wrong type", io::message::critical);    
  }
  
  return math::Vec(0,0,0);
  
}

void util::Virtual_Atom::force(math::VArray const & position, math::Vec const f, math::VArray & force)const
{
 math::Vec s,t,a,b, calc1, calc2, calc3, calch;
 double abs_s, abs_t;
 double m_c; //multiplication constant
 

 switch(m_type){

   case 0: // explicit atom
   case 7: // rotating ring
     force(m_atom[0])+=f;
     break;
   case 1: // CH1
     s = 3.0 * position(m_atom[0]) - position(m_atom[1])
       - position(m_atom[2]) - position(m_atom[3]);
     abs_s = math::abs(s);
     m_c=3*m_dish;
     
     calc1= math::Vec(m_c*(abs_s*abs_s - s(0)*s(0)),
	     -m_c*s(1)*s(0), 
	     -m_c*s(2)*s(0))/(abs_s*abs_s*abs_s) + math::Vec(1,0,0);
     calc2= math::Vec(-m_c*s(1)*s(0),
	     m_c*(abs_s*abs_s - s(1)*s(1)),
	     -m_c*s(2)*s(1))/(abs_s*abs_s*abs_s) + math::Vec(0,1,0);
     calc3= math::Vec(-m_c*s(2)*s(0),
	     -m_c*s(2)*s(1),
	     m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s) + math::Vec(0,0,1) ;     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));


     calc1=m_c/3*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		  s(1)*s(0), 
		  s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2=m_c/3*math::Vec(s(1)*s(0),
		  -(abs_s*abs_s - s(1)*s(1)),
		  s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3=m_c/3*math::Vec(s(2)*s(0),
		  s(2)*s(1),
		  -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[3])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     break;

   case 2: // aromatic H
     s = 2.0 * position(m_atom[0]) - position(m_atom[1])
       - position(m_atom[2]);
     abs_s = math::abs(s);

     calc1= math::Vec(2*m_dish*(abs_s*abs_s - s(0)*s(0)),
	     -2*m_dish*s(1)*s(0), 
	     -2*m_dish*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(-2*m_dish*s(1)*s(0),
	     2*m_dish*(abs_s*abs_s - s(1)*s(1)),
	     -2*m_dish*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(-2*m_dish*s(2)*s(0),
	     -2*m_dish*s(2)*s(1),
	     2*m_dish*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));


     calc1=m_dish*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		   s(1)*s(0), 
		   s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2=m_dish*math::Vec(s(1)*s(0),
		   -(abs_s*abs_s - s(1)*s(1)),
		   s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3=m_dish*math::Vec(s(2)*s(0),
		   s(2)*s(1),
		   -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     break;
   case 3: // non-stereospecific CH2
     s = 2.0 * position(m_atom[0]) - position(m_atom[1])
       - position(m_atom[2]);
     abs_s = math::abs(s);
     m_c=TETHCO*m_dish;
     
     calc1= math::Vec(2*m_c*(abs_s*abs_s - s(0)*s(0)),
	     -2*m_c*s(1)*s(0), 
	     -2*m_c*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(-2*m_c*s(1)*s(0),
	     2*m_c*(abs_s*abs_s - s(1)*s(1)),
	     -2*m_c*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(-2*m_c*s(2)*s(0),
	     -2*m_c*s(2)*s(1),
	     2*m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));


     calc1= m_c*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		 s(1)*s(0), 
		 s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(s(1)*s(0),
		 -(abs_s*abs_s - s(1)*s(1)),
		 s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(s(2)*s(0),
		 s(2)*s(1),
		 -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     break;

   case 4: // stereospecific CH2
     
     DEBUG(8, "FORCE REDISTRIBUTION: case 4!!!");
     
     double m_c_2, abs_a, abs_b;
     
     s = 2.0 * position(m_atom[0]) - position(m_atom[1])- position(m_atom[2]);
     abs_s = math::abs(s);
     DEBUG(10, "\ts = " << math::v2s(s));
     
     t = math::cross(position(m_atom[0]) - position(m_atom[1]), position(m_atom[0]) - position(m_atom[2]));
     abs_t = math::abs(t);
     DEBUG(10, "\tq = " << math::v2s(t));
     
     b = position(m_atom[2])- position(m_atom[1]);
     abs_b = math::abs(b);
     DEBUG(10, "\tb = " << math::v2s(b));

     a = math::cross(t,b);
     abs_a = math::abs(a) ;
     DEBUG(10, "\ta = " << math::v2s(a));
     
     m_c = TETHCO*m_dish;
     m_c_2 = TETHSI*m_dish;
      
     calc1= math::Vec(2*m_c*(abs_s*abs_s - s(0)*s(0)),
	     -2*m_c*s(1)*s(0), 
	     -2*m_c*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);

     calc2= math::Vec(-2*m_c*s(1)*s(0),
	     2*m_c*(abs_s*abs_s - s(1)*s(1)),
	     -2*m_c*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);

     calc3= math::Vec(-2*m_c*s(2)*s(0),
             -2*m_c*s(2)*s(1),
	     2*m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);
      
     DEBUG(10, "A = " << math::v2s(calc1));
     DEBUG(10, "A = " << math::v2s(calc2));
     DEBUG(10, "A = " << math::v2s(calc3));

     calch= -m_c_2*math::Vec(t(0)*a(0), 
		     t(1)*a(0), 
		     t(2)*a(0))/(abs_t*abs_t*abs_t);
     DEBUG(10, "B = " << math::v2s(calch));

     calc1+=calch;     

     calch= -m_c_2*math::Vec(t(0)*a(1),
		     t(1)*a(1),
		     t(2)*a(1))/(abs_t*abs_t*abs_t);
     DEBUG(10, "B = " << math::v2s(calch));
     calc2+=calch;    

     calch= -m_c_2*math::Vec(t(0)*a(2),
		     t(1)*a(2),
		     t(2)*a(2))/(abs_t*abs_t*abs_t);
     DEBUG(10, "B = " << math::v2s(calch));
     calc3+=calch;
       
     calch = m_c_2*math::Vec(0, 
		    b(2), 
		    -b(1))/(abs_t);
     DEBUG(10, "C = " << math::v2s(calch));
     calc1+=calch;

     calch= m_c_2*math::Vec(-b(2),
		    0,
		    b(0))/(abs_t);
     DEBUG(10, "C = " << math::v2s(calch));
     calc2+=calch; 

     calch= m_c_2*math::Vec(b(1),
		    b(0),
		    0)/(abs_t);
     DEBUG(10, "C = " << math::v2s(calch));
     calc3+=calch;
     
     DEBUG(10, "drn/dri = " << math::v2s(calc1));
     DEBUG(10, "drn/dri = " << math::v2s(calc2));
     DEBUG(10, "drn/dri = " << math::v2s(calc3));
     
     force(m_atom[0]) += math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     DEBUG(10, "f(i) = " << math::v2s(math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f))));
      
     calc1= m_c*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		 s(1)*s(0), 
		 s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(s(1)*s(0),
		 -(abs_s*abs_s - s(1)*s(1)),
		 s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(s(2)*s(0),
 		 s(2)*s(1),
		 -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);

     DEBUG(10, "D = " << math::v2s(calc1));
     DEBUG(10, "D = " << math::v2s(calc2));
     DEBUG(10, "D = " << math::v2s(calc3));

     b = position(m_atom[0])- position(m_atom[2]);
     abs_b =  math::abs(b);
     a = math::cross(t,b);
     abs_a = math::abs(a) ;
      
     calch= -m_c_2*math::Vec(t(0)*a(0), 
		     t(1)*a(0), 
		     t(2)*a(0))/(abs_t*abs_t*abs_t);
     DEBUG(10, "E = " << math::v2s(calch));
     calc1+=calch;
     calch= -m_c_2*math::Vec(t(0)*a(1),
		     t(1)*a(1),
		     t(2)*a(1))/(abs_t*abs_t*abs_t); 
     DEBUG(10, "E = " << math::v2s(calch));
     calc2+=calch; 
     calch= -m_c_2*math::Vec(t(0)*a(2),
		     t(1)*a(2),
		     t(2)*a(2))/(abs_t*abs_t*abs_t);
     DEBUG(10, "E = " << math::v2s(calch));
     calc3+=calch;
 
     calch= m_c_2*math::Vec(0, 
		    b(2), 
		    -b(1))/(abs_t);
     DEBUG(10, "F = " << math::v2s(calch));
     calc1+=calch;
     calch= m_c_2*math::Vec(-b(2),
		    0,
		    b(0))/(abs_t);
     DEBUG(10, "F = " << math::v2s(calch));
     calc2+=calch;
     calch= m_c_2*math::Vec(b(1),
		    -b(0),
		    0)/(abs_t);
     DEBUG(10, "F = " << math::v2s(calch));
     calc3+=calch;
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     DEBUG(10, "f(j) = " << math::v2s(math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f))));

      
     calc1= m_c*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		 s(1)*s(0), 
		 s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(s(1)*s(0),
		 -(abs_s*abs_s - s(1)*s(1)),
		 s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(s(2)*s(0),
 		 s(2)*s(1),
		 -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);

     DEBUG(10, "D = " << math::v2s(calc1));
     DEBUG(10, "D = " << math::v2s(calc2));
     DEBUG(10, "D = " << math::v2s(calc3));
     
     b = position(m_atom[1])- position(m_atom[0]);
     abs_b = math::abs(b) ;
     a = math::cross(t,b );
     abs_a = math::abs(a) ;
     
     calch= -m_c_2*math::Vec(t(0)*a(0), 
		    t(1)*a(0), 
		    t(2)*a(0))/(abs_t*abs_t*abs_t) ;
     DEBUG(10, "G = " << math::v2s(calch));
     calc1+=calch;
     calch = -m_c_2*math::Vec(t(0)*a(1),
		     t(1)*a(1),
		     t(2)*a(1))/(abs_t*abs_t*abs_t);
     DEBUG(10, "G = " << math::v2s(calch));
     calc2+=calch;
     calch= -m_c_2*math::Vec(t(0)*a(2),
		    t(1)*a(2),
		    t(2)*a(2))/(abs_t*abs_t*abs_t);
     DEBUG(10, "G = " << math::v2s(calch));
     calc3+=calch;			

	
     calch= m_c_2*math::Vec(0, 
		    b(2), 
		    -b(1))/(abs_t);
     DEBUG(10, "H = " << math::v2s(calch));
     calc1+=calch;
     calch= m_c_2*math::Vec(-b(2),
		    0,
		    b(0))/(abs_t);
     DEBUG(10, "H = " << math::v2s(calch));
     calc2+=calch;
     calch= m_c_2*math::Vec(b(1),
		    -b(0),
		    0)/(abs_t) ;
     DEBUG(10, "H = " << math::v2s(calch));
     calc3+=calch;			
 
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     DEBUG(10, "f(k) = " << math::v2s(math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f))));

     break;

   case 5: // CH3
     s =  position(m_atom[0]) - position(m_atom[1]);
     abs_s = math::abs(s);
     m_c= m_dish/3;
     
     calc1= math::Vec(m_c*(abs_s*abs_s - s(0)*s(0)),
	     -m_c*s(1)*s(0), 
	     -m_c*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(-m_c*s(1)*s(0),
	     m_c*(abs_s*abs_s - s(1)*s(1)),
	     -m_c*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(-m_c*s(2)*s(0),
	     -m_c*s(2)*s(1),
	     m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     calc1= m_c*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		 s(1)*s(0), 
		 s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(s(1)*s(0),
		 -(abs_s*abs_s - s(1)*s(1)),
		 s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(s(2)*s(0),
		 s(2)*s(1),
		 -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     break;

   case 6: // non-stereospecific CH3 (Leu, Val)
     s = 2.0 * position(m_atom[0]) - position(m_atom[1])
       - position(m_atom[2]);
     abs_s=math::abs(s);    
     m_c=-TETHCO* (m_disc + m_dish / 3.0);
     
     calc1= math::Vec(2*m_c*(abs_s*abs_s - s(0)*s(0)),
	     -2*m_c*s(1)*s(0), 
	     -2*m_c*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(-2*m_c*s(1)*s(0),
	     2*m_c*(abs_s*abs_s - s(1)*s(1)),
	     -2*m_c*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(-2*m_c*s(2)*s(0),
	     -2*m_c*s(2)*s(1),
	     2*m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     calc1= m_c*math::Vec(-m_dish*(abs_s*abs_s - s(0)*s(0)),
		 m_dish*s(1)*s(0), 
		 m_dish*s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(m_dish*s(1)*s(0),
		 -m_dish*(abs_s*abs_s - s(1)*s(1)),
		 m_dish*s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(m_dish*s(2)*s(0),
		 m_dish*s(2)*s(1),
		 -m_dish*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     break;
   case 8: // NH2-group (one pseudosite)
     s = 2.0 * position(m_atom[0]) -position(m_atom[1])
       - position(m_atom[2]);
     abs_s=math::abs(s);
     
     calc1= math::Vec(-m_dish*(abs_s*abs_s - s(0)*s(0)),
	     m_dish*s(1)*s(0), 
	     m_dish*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(m_dish*s(1)*s(0),
	     -m_dish*(abs_s*abs_s - s(1)*s(1)),
	     m_dish*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(m_dish*s(2)*s(0),
	     m_dish*s(2)*s(1),
	     -m_dish*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     calc1=-0.5*m_dish*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
			s(1)*s(0), 
			s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2=-0.5*m_dish*math::Vec(s(1)*s(0),
			-(abs_s*abs_s - s(1)*s(1)),
			s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3=-0.5*m_dish*math::Vec(s(2)*s(0),
			s(2)*s(1),
			-(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     force(m_atom[2])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     break;

   case 9: // (CH3)3-group (one psuedosite)
     s =  position(m_atom[0]) - position(m_atom[1]);
     abs_s = math::abs(s);
     m_c=  ( m_disc + m_dish/3.0 )/ 3;
  
     calc1= math::Vec(m_c*(abs_s*abs_s - s(0)*s(0)),
	     -m_c*s(1)*s(0), 
	     -m_c*s(2)*s(0))/(abs_s*abs_s*abs_s)+math::Vec(1,0,0);
     calc2= math::Vec(-m_c*s(1)*s(0),
	     m_c*(abs_s*abs_s - s(1)*s(1)),
	     -m_c*s(2)*s(1))/(abs_s*abs_s*abs_s)+math::Vec(0,1,0);
     calc3= math::Vec(-m_c*s(2)*s(0),
	     -m_c*s(2)*s(1),
	     m_c*(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s)+math::Vec(0,0,1);     
     force(m_atom[0])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));

     calc1= m_c*math::Vec(-(abs_s*abs_s - s(0)*s(0)),
		 s(1)*s(0), 
		 s(2)*s(0))/(abs_s*abs_s*abs_s);
     calc2= m_c*math::Vec(s(1)*s(0),
		 -(abs_s*abs_s - s(1)*s(1)),
		 s(2)*s(1))/(abs_s*abs_s*abs_s);
     calc3= m_c*math::Vec(s(2)*s(0),
		 s(2)*s(1),
		 -(abs_s*abs_s - s(2)*s(2)))/(abs_s*abs_s*abs_s);
     force(m_atom[1])+=math::Vec(math::dot(calc1,f),math::dot(calc2,f),math::dot(calc3,f));
     break;
     
   default:
     std::cerr <<"Type not implemented";
     assert(false);
     break;
         
     
 }
 
}
