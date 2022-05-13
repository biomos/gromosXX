/**
 * @file in_perturbation.cc
 * implements methods of In_Perturbation.
 */

#include <climits>
#include <stdheader.h>

#include <topology/topology.h>
#include <simulation/multibath.h>
#include <simulation/parameter.h>
#include <io/instream.h>

#include <io/blockinput.h>

#include "in_perturbation.h"
#include <algorithm>


#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

static std::set<std::string> block_read;

/**
 * Constructor.
 */

io::In_Perturbation::In_Perturbation(std::istream &is) 
  : GInStream(is) 
{
  // read the whole file at beginning
  readStream();
};

template<typename T>
bool check_type(std::vector<std::string> const & buffer, std::vector<T> term)
{
  if (buffer.size()){
    std::istringstream is(buffer[1]);
    int num = 0;
    if (!(is >> num) || num < 0)
      return false;
      
    for(typename std::vector<T>::const_iterator
	  it = term.begin(),
	  to = term.end();
	it != to;
	++it){
      
      if (int(it->A_type) >= num && int(it->B_type) >= num)
	return false;
    }
  }
  else return false;
  return true;
}

/**
 * @class bondMatcher
 * helper function object to match two_body_term_structs
 */
class bondMatcher {
  public:
  bondMatcher(const topology::two_body_term_struct& arg) : a(arg) {}
  bool operator()(const topology::two_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j);
  }
  bool operator()(const topology::perturbed_two_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j);
  }
  private:
  topology::two_body_term_struct a;
};

/**
 * @class angleMatcher
 * helper function object to match three_body_term_structs
 */
class angleMatcher {
  public:
  angleMatcher(const topology::three_body_term_struct& arg) : a(arg) {}
  bool operator()(const topology::three_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j && a.k == x.k);
  }
  bool operator()(const topology::perturbed_three_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j && a.k == x.k);
  }
  private:
  topology::three_body_term_struct a;
};

/**
 * @class dihedralMatcher
 * helper function object to match four_body_term_structs
 */
class dihedralMatcher {
  public:
  dihedralMatcher(const topology::four_body_term_struct& arg) : a(arg) {}
  bool operator()(const topology::four_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j && a.k == x.k && a.l == x.l);
  }
  bool operator()(const topology::perturbed_four_body_term_struct &x) {
    return (a.i == x.i && a.j == x.j && a.k == x.k && a.l == x.l);
  }
  private:
  topology::four_body_term_struct a;
};

void
io::In_Perturbation::read(topology::Topology &topo,
	simulation::Parameter &param,
        std::ostream & os)
{

  if (!param.perturbation.perturbation && !param.eds.eds && !param.reeds.reeds){
    io::messages.add("Ignoring perturbation topology because perturbation is not enabled.",
		     "In_Perturbation",
		     io::message::warning);
    return;
  }

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  DEBUG(7, "Reading PERTURBATION");
  if (!quiet){
    os << "PERTURBATION TOPOLOGY\n";
    os << title << "\n";
  }
  
  
  // let's do a warning because state A information is overwritten?
  bool warn = false;
  if(!param.eds.eds){
    // prepare arrays
    topo.is_perturbed().resize(topo.num_solute_atoms(), false);
    
    { // PERTBONDSTRETCH(H)

      // Chris: the next two lines seemed to be missing...
      DEBUG(10, "PERTBOND03 block");
      buffer = m_block["PERTBOND03"];
      if (buffer.size()){
        block_read.insert("PERTBOND03");
        io::messages.add("The PERTBOND03 block was renamed to PERTBONDSTRETCH.",
                "In_Perturbation", io::message::error);
      }
      
      std::vector<std::string> pertbondstretch;
      pertbondstretch.push_back("PERTBONDSTRETCHH");
      pertbondstretch.push_back("PERTBONDSTRETCH");
      for (unsigned int hh=0; hh < pertbondstretch.size(); hh++) {
        buffer = m_block[pertbondstretch.at(hh)];
        if (buffer.size()){
          
          block_read.insert(pertbondstretch.at(hh));
          
          if (!quiet)
            os << "\t" << pertbondstretch.at(hh) << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (param.constraint.ntc == 2){
            io::messages.add("No perturbed distance constraints for "
            "NTC = 2 from perturbed bonds",
                    "in_perturbation",
                    io::message::warning);
          }
          else if (param.constraint.ntc == 3){
            if (!quiet)
              os << "\n\t\t"
              << num
              << " perturbed bonds from " << pertbondstretch.at(hh) << " block added to "
              << "perturbed distance constraints.";
          }
          
          if (!quiet)
            os << "\n\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, t_A = 0, t_B = 0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> t_A >> t_B;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in " + pertbondstretch.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
            }
            
            topology::two_body_term_struct b(i-1, j-1, t_A-1);
            
            if (param.constraint.ntc != 3){
              std::vector<topology::two_body_term_struct>::iterator b_it
                      = std::find_if(topo.solute().bonds().begin(),
                      topo.solute().bonds().end(),
                      bondMatcher(b));
              
              if (b_it == topo.solute().bonds().end()){
              // check if this bond is already in the soft perturbed ones
              std::vector<topology::perturbed_two_body_term_struct>::iterator pb_it
                    = std::find_if(topo.perturbed_solute().softbonds().begin(),
                    topo.perturbed_solute().softbonds().end(),
                    bondMatcher(b));
              if (pb_it != topo.perturbed_solute().softbonds().end()){
                io::messages.add("Bond found both in PERTBONDSTRETCH(H) and PERTBONDSOFT blocks!",
                    "In_Perturbation", io::message::error);
              } else {
                io::messages.add("Perturbation of a non-existing bond "
                "in " + pertbondstretch.at(hh) + " block.",
                        "In_Perturbation", io::message::error);
              }
                return;
              }
              
              if (b_it->type != b.type)
                warn = true;
              
              topo.solute().bonds().erase(b_it);
              topology::perturbed_two_body_term_struct
              pb(i-1, j-1, t_A-1, t_B-1);
              
              if (!quiet)
                os << "\t"
                << std::setw(10) << pb.i+1
                << std::setw(10) << pb.j+1
                << std::setw(10) << pb.A_type+1
                << std::setw(10) << pb.B_type+1
                << "\n";
              
              topo.perturbed_solute().bonds().push_back(pb);
            } else { // we have constraints
              std::vector<topology::two_body_term_struct>::iterator b_it
                      = std::find_if(topo.solute().distance_constraints().begin(),
                      topo.solute().distance_constraints().end(),
                      bondMatcher(b));
              
              if (b_it == topo.solute().distance_constraints().end()){
              // check if this bond is already in the soft perturbed ones
              std::vector<topology::perturbed_two_body_term_struct>::iterator pb_it
                    = std::find_if(topo.perturbed_solute().softbonds().begin(),
                    topo.perturbed_solute().softbonds().end(),
                    bondMatcher(b));
              if (pb_it != topo.perturbed_solute().softbonds().end()){
                io::messages.add("Bond found both in PERTBONDSTRETCH(H) and PERTBONDSOFT blocks!",
                    "In_Perturbation", io::message::error);
              } else {
                io::messages.add("Perturbation of a non-existing distance "
                "constraint in " + pertbondstretch.at(hh) + " block.",
                        "In_Perturbation", io::message::error);
                }
                return;
              }
              
              if (b_it->type != b.type)
                warn = true;
              
              topo.solute().distance_constraints().erase(b_it);
              topology::perturbed_two_body_term_struct
              pb(i-1, j-1, t_A-1, t_B-1);
              
              if (!quiet)
                os << "\t"
                << std::setw(10) << pb.i+1
                << std::setw(10) << pb.j+1
                << std::setw(10) << pb.A_type+1
                << std::setw(10) << pb.B_type+1
                << "\n";
              
              topo.perturbed_solute().distance_constraints().push_back(pb);
            }
          }
          
          if (n != num){
            io::messages.add("Wrong number of bonds in " + pertbondstretch.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in " + pertbondstretch.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet){
            os << "\n\t\tbonds :                          "
               << unsigned(topo.solute().bonds().size())
               << "\n\t\tperturbed bonds :                "
               << unsigned(topo.perturbed_solute().bonds().size())
               << "\n\t\tdistance constraints :           "
               << unsigned(topo.solute().distance_constraints().size())
               << "\n\t\tperturbed distance constraints : "
               << unsigned(topo.perturbed_solute().distance_constraints().size())
               << "\n\n";

            //write a warning if we do perturbed constraints and PRECALCLAM
            if(topo.perturbed_solute().distance_constraints().size()
               && param.precalclam.nr_lambdas) {
              os << "\tWARNING: Perturbed distance constraints in combination with\n\t\t precalculated lambdas will give a contribution to energy\n\t\t derivatives, but not to the energies. Use in BAR or reweighting\n\t\t in extended TI will be inappropriate\n\n";

              io::messages.add("Perturbed SHAKE in combination with precalculated\n                            lambdas only contributes to derivatives not energies", "In_Perturbation", io::message::warning);
            }

            os << "\tEND\n";
          }

          
        } // if block present
      }  // loop over H/non H blocks
    } // PERTBONDSTRETCH(H)
    
    
    { // PERTBONDSOFT
        buffer = m_block["PERTBONDSOFT"];
        if (buffer.size()){
          
          block_read.insert("PERTBONDSOFT");
          
          if (!quiet)
            os << "\t" << "PERTBONDSOFT" << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\n\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << std::setw(10) << "ALB"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, t_A = 0, t_B = 0;
            double alphab = 0.0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> t_A >> t_B >> alphab;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in PERTBONDSOFT block.",
                      "In_Perturbation", io::message::error);
            }
            
            // does a bond or constraint exist between the two atoms?
            // if yes we remove it, it will be handled in the perturbed-soft-bonds
            topology::two_body_term_struct b(i-1, j-1, t_A-1);
            bool bond_exists=false, constraint_exists=false;
            std::vector<topology::two_body_term_struct>::iterator b_it, c_it;
            b_it = std::find_if(topo.solute().bonds().begin(),
                      topo.solute().bonds().end(),
                      bondMatcher(b));
                  if (b_it != topo.solute().bonds().end()){
                     bond_exists=true;
                     if (b_it->type != b.type) warn = true;
                     topo.solute().bonds().erase(b_it);
                  }
            c_it = std::find_if(topo.solute().distance_constraints().begin(),
                      topo.solute().distance_constraints().end(),
                      bondMatcher(b));
              
                  if (c_it != topo.solute().distance_constraints().end()){
                     constraint_exists=true;
                     if (c_it->type != b.type) warn = true;
                     topo.solute().distance_constraints().erase(c_it); 
                  }
              // check if this bond is already in the perturbed ones
              std::vector<topology::perturbed_two_body_term_struct>::iterator pb_it
                    = std::find_if(topo.perturbed_solute().bonds().begin(),
                    topo.perturbed_solute().bonds().end(),
                    bondMatcher(b));
              std::vector<topology::perturbed_two_body_term_struct>::iterator pc_it
                    = std::find_if(topo.perturbed_solute().distance_constraints().begin(),
                    topo.perturbed_solute().distance_constraints().end(),
                    bondMatcher(b));
              if (pb_it != topo.perturbed_solute().bonds().end() || pc_it != topo.perturbed_solute().distance_constraints().end()){
                io::messages.add("Bond found both in PERTBONDSTRETCH(H) and PERTBONDSOFT blocks!",
                    "In_Perturbation", io::message::error);
                    return;
              }

            
            
            if (t_A > 0 && ! bond_exists && ! constraint_exists) {
                    io::messages.add("Perturbation of a non-existing bond "
                    "in PERTBONDSOFT block.",
                        "In_Perturbation", io::message::error);
                    return;
            }
            if (t_A == 0 && t_B != 0 && (bond_exists || constraint_exists)){
                    io::messages.add("PERTBONDSOFT block: Adding a bond that already exists",
                        "In_Perturbation", io::message::warning);
                   // return;
            }

            int t_A_copy=t_A, t_B_copy=t_B;
            // if given bond_type is <= 0 create a dummy bond with the same
            // length as in the other state and K=0
            if (t_A <= 0 && t_B <= 0) {
                io::messages.add("Conversion from bond type 0 to bond type 0 "
                "in PERTBONDSOFT block.",
                        "In_Perturbation", io::message::error);
                return;
            } else if (t_A <= 0) {
              t_A=INT_MAX;
              //store type of stateB, we will retrieve the bond length for the new bond from there
              topo.perturbed_solute().soft_bond_types().push_back(t_B-1);
              io::messages.add("You are making a bond, make sure you adapt the excluded atoms to this change.", io::message::warning);
            } else if (t_B <= 0) {
              t_B=INT_MAX;
              topo.perturbed_solute().soft_bond_types().push_back(t_A-1);
              io::messages.add("You are breaking a bond, make sure you adapt the excluded atoms to this change.", io::message::warning);
            }          
             topo.perturbed_solute().alpha_bond().push_back(alphab);
              
              topology::perturbed_two_body_term_struct
              pb(i-1, j-1, t_A-1, t_B-1);
              
              if (!quiet)
                os << "\t"
                << std::setw(10) << pb.i+1
                << std::setw(10) << pb.j+1
                << std::setw(10) << t_A_copy
                << std::setw(10) << t_B_copy
                << std::setw(10) << alphab
                << "\n";
              
              topo.perturbed_solute().softbonds().push_back(pb);
          }
          
          if (n != num){
            io::messages.add("Wrong number of bonds in PERTBONDSOFT block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in PERTBONDSOFT block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet)
            os << "\n\t\tbonds :                          "
            << unsigned(topo.solute().bonds().size())
            << "\n\t\tperturbed soft bonds :           "
            << unsigned(topo.perturbed_solute().softbonds().size())
            << "\n\t\tdistance constraints :           "
            << unsigned(topo.solute().distance_constraints().size())
            << "\n\t\tperturbed distance constraints : "
            << unsigned(topo.perturbed_solute().distance_constraints().size())
            << "\n\n"
            << "\tEND\n";
          
        } // if block present
    } // PERTBONDSOFT
    
    { // PERTCONSTRAINT03
      DEBUG(10, "PERTCONSTRAINT03 block");
      buffer = m_block["PERTCONSTRAINT03"];
      
      if (buffer.size() && param.constraint.ntc != 1){
        block_read.insert("PERTCONSTRAINT03");
        
        it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;
        
        if (!quiet)
          os << "\tPERTCONSTRAINT03\n\t\t"
          << num
          << " bonds in PERTCONSTRAINT03 block."
          << "\n\t\ttotal of perturbed constraint bonds : "
          << unsigned(num + topo.perturbed_solute().distance_constraints().size())
          << "\n"
          << "\t"
          << std::setw(10) << "atom i"
          << std::setw(10) << "atom j"
          << std::setw(10) << "type A"
          << std::setw(10) << "type B"
          << "\n";
        
        for(n=0; it != buffer.end() - 1; ++it, ++n){
          int i = 0, j = 0, t_A = 0, t_B = 0;
          
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> t_A >> t_B;
          
          if (_lineStream.fail()){
            io::messages.add("Bad line in PERTCONSTRAINT03 block",
                    "In_Topology", io::message::error);
          }
          
          if (i > int(topo.num_solute_atoms()) ||
                  j > int(topo.num_solute_atoms()) ||
                  i < 1 || j < 1){
            io::messages.add("Atom number out of range in PERTCONSTRAINT03 "
            " block", "In_Topology", io::message::error);
          }
          
          topology::two_body_term_struct b(i-1, j-1, t_A-1);
          
          std::vector<topology::two_body_term_struct>::iterator b_it
                  = std::find_if(topo.solute().distance_constraints().begin(),
                  topo.solute().distance_constraints().end(),
                  bondMatcher(b));
          
          if (b_it == topo.solute().distance_constraints().end()){
            io::messages.add("Perturbation of a non-existing distance "
            "constraint in PERTCONSTRAINT03 block.",
                    "In_Perturbation", io::message::error);
            
          }
          
          if (b_it->type != b.type)
            warn = true;
          
          topo.solute().distance_constraints().erase(b_it);
          topology::perturbed_two_body_term_struct
          pb(i-1, j-1, t_A-1, t_B-1);
          
          topo.perturbed_solute().distance_constraints().push_back(pb);
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << pb.i+1
            << std::setw(10) << pb.j+1
            << std::setw(10) << pb.A_type+1
            << std::setw(10) << pb.B_type+1
            << "\n";
          
        }
        
        if(n != num){
          io::messages.add("Wrong number of bonds in PERTCONSTRAINT03 block",
                  "In_Perturbation", io::message::error);
        }
      }
      
    } // PERTCONSTRAINT03
    
    { // PERTBONDANGLE(H)
      DEBUG(10, "PERTBANGLE03 block");
      buffer = m_block["PERTBANGLE03"];
      if (buffer.size()){
        
        block_read.insert("PERTBANGLE03");
        io::messages.add("The PERTBANGLE03 block was renamed to PERTBONDANGLE.",
                "In_Perturbation", io::message::error);
      }
      
      std::vector<std::string> pertbondangle;
      pertbondangle.push_back("PERTBONDANGLEH");
      pertbondangle.push_back("PERTBONDANGLE");
      for (unsigned int hh=0; hh < pertbondangle.size(); hh++) {
        buffer = m_block[pertbondangle.at(hh)];
        if (buffer.size()){
          block_read.insert(pertbondangle.at(hh));
          
          if (!quiet)
            os << "\t" << pertbondangle.at(hh) << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "atom k"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, k = 0, t_A = 0, t_B = 0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> k >> t_A >> t_B;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in " + pertbondangle.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
            }
            
            topology::three_body_term_struct a(i-1, j-1, k-1, t_A-1);
            std::vector<topology::three_body_term_struct>::iterator a_it
                    = std::find_if(topo.solute().angles().begin(),
                    topo.solute().angles().end(),
                    angleMatcher(a));
            
            if (a_it == topo.solute().angles().end()){
              // check if this angle is already in the soft perturbed ones
              std::vector<topology::perturbed_three_body_term_struct>::iterator pa_it
                    = std::find_if(topo.perturbed_solute().softangles().begin(),
                    topo.perturbed_solute().softangles().end(),
                    angleMatcher(a));
              if (pa_it != topo.perturbed_solute().softangles().end()){
                io::messages.add("Angle found both in PERTBONDANGLE(H) and PERTANGLESOFT blocks!",
                    "In_Perturbation", io::message::error);
              } else {
                io::messages.add("Perturbation of a non-existing angle in "
                + pertbondangle.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
              }
              return;
            }
            
            if (a_it->type != a.type)
              warn = true;
            
            topo.solute().angles().erase(a_it);
            topology::perturbed_three_body_term_struct pa(i-1, j-1, k-1, t_A-1, t_B-1);
            topo.perturbed_solute().angles().push_back(pa);
            
            if (!quiet)
              os << "\t"
              << std::setw(10) << pa.i+1
              << std::setw(10) << pa.j+1
              << std::setw(10) << pa.k+1
              << std::setw(10) << pa.A_type+1
              << std::setw(10) << pa.B_type+1
              << "\n";
            
          }
          if (n != num){
            io::messages.add("Wrong number of angles in " + pertbondangle.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in " + pertbondangle.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet)
            os << "\tEND\n";
          
        } // if block present
      } // loop over H/non H blocks
    } // PERTBONDANGLE(H)
    
    
    { // PERTANGLESOFT
        buffer = m_block["PERTANGLESOFT"];
        if (buffer.size()){
          
          block_read.insert("PERTANGLESOFT");
          
          if (!quiet)
            os << "\t" << "PERTANGLESOFT" << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\n\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "atom k"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << std::setw(10) << "ALA"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, k = 0, t_A = 0, t_B = 0;
            double alphaa = 0.0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> k >> t_A >> t_B >> alphaa;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in PERTANGLESOFT block.",
                      "In_Perturbation", io::message::error);
            }
            
            
            // does an angle exist between the atoms?
            // if yes we remove it, it will be handled in the perturbed-soft-angles
            topology::three_body_term_struct a(i-1, j-1, k-1, t_A-1);
            bool angle_exists=false;
            std::vector<topology::three_body_term_struct>::iterator b_it;
            b_it = std::find_if(topo.solute().angles().begin(),
                      topo.solute().angles().end(),
                      angleMatcher(a));
                  if (b_it != topo.solute().angles().end()){
                     angle_exists=true;
                     if (b_it->type != a.type) warn = true;
                     topo.solute().angles().erase(b_it);
                  }
              // check if this angle is already in the perturbed ones
              std::vector<topology::perturbed_three_body_term_struct>::iterator pa_it
                    = std::find_if(topo.perturbed_solute().angles().begin(),
                    topo.perturbed_solute().angles().end(),
                    angleMatcher(a));
              if (pa_it != topo.perturbed_solute().angles().end()){
                io::messages.add("Angle found both in PERTBONDANGLE(H) and PERTANGLESOFT blocks!",
                    "In_Perturbation", io::message::error);
                    return;
              }

            
            
            if (t_A > 0 && ! angle_exists) {
                    io::messages.add("Perturbation of a non-existing angle "
                    "in PERTANGLESOFT block.",
                        "In_Perturbation", io::message::error);
                    return;
            }
            if (t_A == 0 && t_B != 0 && angle_exists) {
                    io::messages.add("PERTANGLESOFT block: Adding an angle that already exists",
                        "In_Perturbation", io::message::warning);
                   // return;
            }
         

            int t_A_copy=t_A, t_B_copy=t_B;   
            // if given angle_type is <= 0 create a dummy angle with the same
            // length as in the other state and K=0
            if (t_A <= 0 && t_B <= 0) {
                io::messages.add("Conversion from angle type 0 to angle type 0 "
                "in PERTANGLESOFT block.",
                        "In_Perturbation", io::message::error);
                return;
            } else if (t_A <= 0) {
              t_A=INT_MAX;
              //store type of stateB, we will retrieve the angle length for the new angle from there
              topo.perturbed_solute().soft_angle_types().push_back(t_B-1);
            } else if (t_B <= 0) {
              t_B=INT_MAX;
              topo.perturbed_solute().soft_angle_types().push_back(t_A-1);
            }          
             topo.perturbed_solute().alpha_angle().push_back(alphaa);
              
              topology::perturbed_three_body_term_struct
              pb(i-1, j-1, k-1, t_A-1, t_B-1);
              
              if (!quiet)
                os << "\t"
                << std::setw(10) << pb.i+1
                << std::setw(10) << pb.j+1
                << std::setw(10) << pb.k+1
                << std::setw(10) << t_A_copy
                << std::setw(10) << t_B_copy
                << std::setw(10) << alphaa
                << "\n";
              
              topo.perturbed_solute().softangles().push_back(pb);
          }
          
          if (n != num){
            io::messages.add("Wrong number of angles in PERTANGLESOFT block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in PERTANGLESOFT block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet)
            os << "\n\t\tangles :                          "
            << unsigned(topo.solute().angles().size())
            << "\n\t\tperturbed soft angles :           "
            << unsigned(topo.perturbed_solute().softangles().size())
            << "\n\n"
            << "\tEND\n";
          
        } // if block present
    } // PERTANGLESOFT
    
    { // PERTIMPROPERDIH(H)
      DEBUG(10, "PERTIMPDIHEDRAL03 block");
      buffer = m_block["PERTIMPDIHEDRAL03"];
      if (buffer.size()){
        
        block_read.insert("PERTIMPDIHEDRAL03");
        io::messages.add("The PERTIMPDIHEDRAL03 block was renamed to PERTIMPROPERDIH.",
                "In_Perturbation", io::message::error);
      }
      
      std::vector<std::string> pertimproperdih;
      pertimproperdih.push_back("PERTIMPROPERDIHH");
      pertimproperdih.push_back("PERTIMPROPERDIH");
      for (unsigned int hh=0; hh < pertimproperdih.size(); hh++) {
        buffer = m_block[pertimproperdih.at(hh)];
        if (buffer.size()){
          block_read.insert(pertimproperdih.at(hh));
          
          if (!quiet)
            os << "\t" << pertimproperdih.at(hh) << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "atom k"
            << std::setw(10) << "atom l"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, k = 0, l = 0, t_A = 0, t_B = 0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> k >> l >> t_A >> t_B;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in " + pertimproperdih.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
            }
            
            topology::four_body_term_struct id(i-1, j-1, k-1, l-1, t_A-1);
            std::vector<topology::four_body_term_struct>::iterator id_it
                    = std::find_if(topo.solute().improper_dihedrals().begin(),
                    topo.solute().improper_dihedrals().end(),
                    dihedralMatcher(id));
            
            if (id_it == topo.solute().improper_dihedrals().end()){
              
              // check if this improper is already in the soft perturbed ones
              std::vector<topology::perturbed_four_body_term_struct>::iterator pid_it 
                    = std::find_if(topo.perturbed_solute().softimpropers().begin(),
                    topo.perturbed_solute().softimpropers().end(),
                    dihedralMatcher(id));
              if (pid_it != topo.perturbed_solute().softimpropers().end()){
                io::messages.add("Improper dihedral found both in PERTIMPROPERDIH and PERTIMPROPERDIHSOFT blocks!",
                    "In_Perturbation", io::message::error);
              } else {
              io::messages.add("Perturbation of a non-existing improper dihedral in "
              + pertimproperdih.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
              }
              return;
            }
            
            if (id_it->type != id.type)
              warn = true;
            
            topo.solute().improper_dihedrals().erase(id_it);
            topology::perturbed_four_body_term_struct pid(i-1, j-1, k-1, l-1,
                    t_A-1, t_B-1);
            
            if (!quiet)
              os << "\t"
              << std::setw(10) << pid.i+1
              << std::setw(10) << pid.j+1
              << std::setw(10) << pid.k+1
              << std::setw(10) << pid.l+1
              << std::setw(10) << pid.A_type+1
              << std::setw(10) << pid.B_type+1
              << "\n";
            
            topo.perturbed_solute().improper_dihedrals().push_back(pid);
          }
          
          if (n != num){
            io::messages.add("Wrong number of bonds in " + pertimproperdih.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in " + pertimproperdih.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          if (!quiet)
            os << "\tEND\n";
          
        } // if block present
      } // loop over H/non H blocks
    } // PERTIMPROPERDIH(H)
    
        
    { // PERTIMPROPERDIHSOFT
        buffer = m_block["PERTIMPROPERDIHSOFT"];
        if (buffer.size()){
          
          block_read.insert("PERTIMPROPERDIHSOFT");
          
          if (!quiet)
            os << "\t" << "PERTIMPROPERDIHSOFT" << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\n\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "atom k"
            << std::setw(10) << "atom l"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << std::setw(10) << "ALI"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, k = 0, l = 0, t_A = 0, t_B = 0;
            double alphaa = 0.0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> k >> l >> t_A >> t_B >> alphaa;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in PERTIMPROPERDIHSOFT block.",
                      "In_Perturbation", io::message::error);
            }
            
            // does the improper exist?
            // if yes we remove it, it will be handled in the perturbed-soft-impropers
            topology::four_body_term_struct b(i-1, j-1, k-1, l-1, t_A-1);
            bool imp_exists=false;
            std::vector<topology::four_body_term_struct>::iterator b_it;
            b_it = std::find_if(topo.solute().improper_dihedrals().begin(),
                      topo.solute().improper_dihedrals().end(),
                      dihedralMatcher(b));
                  if (b_it != topo.solute().improper_dihedrals().end()){
                     imp_exists=true;
                     if (b_it->type != b.type) warn = true;
                     topo.solute().improper_dihedrals().erase(b_it);
                  }
              // check if this improper is already in the perturbed ones
              std::vector<topology::perturbed_four_body_term_struct>::iterator pid_it 
                    = std::find_if(topo.perturbed_solute().improper_dihedrals().begin(),
                    topo.perturbed_solute().improper_dihedrals().end(),
                    dihedralMatcher(b));
              if (pid_it != topo.perturbed_solute().improper_dihedrals().end()){
                io::messages.add("Improper dihedral found both in PERTIMPROPERDIH and PERTIMPROPERDIHSOFT blocks!",
                    "In_Perturbation", io::message::error);
                    return;
              }
          
            if (t_A > 0 && ! imp_exists) {   
                    io::messages.add("Perturbation of a non-existing improper "
                    "in PERTIMPROPERDIHSOFT block.",
                        "In_Perturbation", io::message::error);
                    return;
            }
            if (t_A <= 0 && t_B > 0 && imp_exists) {
                     io::messages.add("PERTIMPROPERDIHSOFT block: Adding an improper that exists in the topology",
                        "In_Perturbation", io::message::warning);
                   // return;
            }

            int t_A_copy=t_A, t_B_copy=t_B;   
            // if given improper_type is <= 0 create a dummy improper with the same
            // length as in the other state and K=0
            if (t_A <= 0 && t_B <= 0) {
                io::messages.add("Conversion from improper type 0 to improper type 0 "
                "in PERTIMPROPERDIHSOFT block.",
                        "In_Perturbation", io::message::error);
                return;
            } else if (t_A <= 0) {
              t_A=INT_MAX;
              //store type of stateB, we will retrieve the value for the new improper from there
              topo.perturbed_solute().soft_improper_types().push_back(t_B-1);
            } else if (t_B <= 0) {
              t_B=INT_MAX;
              topo.perturbed_solute().soft_improper_types().push_back(t_A-1);
            }          
             topo.perturbed_solute().alpha_improper().push_back(alphaa);
              
              topology::perturbed_four_body_term_struct
              pb(i-1, j-1, k-1, l-1, t_A-1, t_B-1);
              
              if (!quiet)
                os << "\t"
                << std::setw(10) << pb.i+1
                << std::setw(10) << pb.j+1
                << std::setw(10) << pb.k+1
                << std::setw(10) << pb.l+1
                << std::setw(10) << t_A_copy
                << std::setw(10) << t_B_copy
                << std::setw(10) << alphaa
                << "\n";
              
              topo.perturbed_solute().softimpropers().push_back(pb);
          }
          
          if (n != num){
            io::messages.add("Wrong number of improper dihedrals in PERTIMPROPERDIHSOFT block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in PERTIMPROPERDIHSOFT block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet)
            os << "\n\t\timproper dihedrals :                          "
            << unsigned(topo.solute().improper_dihedrals().size())
            << "\n\t\tperturbed soft improper dihedrals :           "
            << unsigned(topo.perturbed_solute().softimpropers().size())
            << "\n\n"
            << "\tEND\n";
          
        } // if block present
    } // PERTIMPROPERDIHSOFT
    
    { // PERTPROPERDIH(H)
      
      std::vector<std::string> pertproperdih;
      pertproperdih.push_back("PERTPROPERDIHH");
      pertproperdih.push_back("PERTPROPERDIH");
      for (unsigned int hh=0; hh < pertproperdih.size(); hh++) {
        buffer = m_block[pertproperdih.at(hh)];
        if (buffer.size()){
          block_read.insert(pertproperdih.at(hh));
          
          if (!quiet)
            os << "\t" << pertproperdih.at(hh) << "\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << "atom i"
            << std::setw(10) << "atom j"
            << std::setw(10) << "atom k"
            << std::setw(10) << "atom l"
            << std::setw(10) << "type A"
            << std::setw(10) << "type B"
            << "\n";
          
          for(n=0; it != buffer.end() -1; ++it, ++n){
            int i = 0, j = 0, k = 0, l = 0, t_A = 0, t_B = 0;
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> k >> l >> t_A >> t_B;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in " + pertproperdih.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
            }
            
            topology::four_body_term_struct id(i-1, j-1, k-1, l-1, t_A-1);
            std::vector<topology::four_body_term_struct>::iterator id_it
                    = std::find_if(topo.solute().dihedrals().begin(),
                    topo.solute().dihedrals().end(),
                    dihedralMatcher(id));
            
            if (id_it == topo.solute().dihedrals().end()){
              io::messages.add("Perturbation of a non-existing dihedral in "
              + pertproperdih.at(hh) + " block.",
                      "In_Perturbation", io::message::error);
              return;
            }
            
            if (id_it->type != id.type)
              warn = true;
            
            topo.solute().dihedrals().erase(id_it);
            topology::perturbed_four_body_term_struct pid(i-1, j-1, k-1, l-1,
                    t_A-1, t_B-1);
            
            if (!quiet)
              os << "\t"
              << std::setw(10) << pid.i+1
              << std::setw(10) << pid.j+1
              << std::setw(10) << pid.k+1
              << std::setw(10) << pid.l+1
              << std::setw(10) << pid.A_type+1
              << std::setw(10) << pid.B_type+1
              << "\n";
            
            topo.perturbed_solute().dihedrals().push_back(pid);
          }
          
          if (n != num){
            io::messages.add("Wrong number of bonds in " + pertproperdih.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in " + pertproperdih.at(hh) + " block.",
                    "In_Perturbation", io::message::error);
          }
          
          if (!quiet)
            os << "\tEND\n";
          
        } // if block present
      } // loop over H/non H blocks
    } // PERTPROPERDIH(H)
    
    { // PERTATOMPAIR
      // has to be read in before(!!) PERTATOMPARAM
      // because the exclusions and 1,4 exclusions have to be adapted...
      
      DEBUG(10, "PERTATOMPAIR03 block");
      buffer = m_block["PERTATOMPAIR03"];
      
      if (buffer.size()){
        
        block_read.insert("PERTATOMPAIR03");
        io::messages.add("The PERTATOMPAIR03 block was renamed to PERTATOMPAIR.",
                "In_Perturbation", io::message::error);
      }
      
      buffer = m_block["PERTATOMPAIR"];
      if (buffer.size()){
        block_read.insert("PERTATOMPAIR");
        
        if (!quiet)
          os << "\tPERTATOMPAIR\n";
        
        it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;

        int i = 0, j = 0, A = 0, B = 0;
        
        if (!quiet)
          os << "\t"
          << std::setw(10) << "atom i"
          << std::setw(10) << "atom j"
          << std::setw(10) << "type A"
          << std::setw(10) << "type B"
          << "\n";
        
        for(n = 0; it != buffer.end() - 1; ++it, ++n){
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> i >> j >> A >> B;
          
          if (_lineStream.fail()){
            io::messages.add("Bad line in PERTATOMPAIR block.",
                    "In_Perturbation", io::message::error);
          }
          
          if(j<i) { int t=j; j=i; i=t; }
          
          topology::perturbed_two_body_term_struct ap(i-1, j-1, A, B);
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << ap.i+1
            << std::setw(10) << ap.j+1
            << std::setw(10) << ap.A_type
            << std::setw(10) << ap.B_type
            << std::endl;
          
          topo.perturbed_solute().atompairs().push_back(ap);
          
          // make sure it's excluded
          if (topo.all_exclusion(ap.i).insert(ap.j)){
            DEBUG(7, "excluding perturbed pair " << ap.i << " and " << ap.j);
          }
          else{
            // it was already excluded, let's remove it from the
            // exclusions or 1,4 pairs...
            
            // is it in the exclusions
            if (topo.exclusion(ap.i).find_and_remove(ap.j)){
              DEBUG(7, "removing perturbed pair from exclusion "
                      << ap.i << " and " << ap.j);
            }
            if (topo.one_four_pair(ap.i).find_and_remove(ap.j)){
              DEBUG(7, "removing perturbed pair from one four "
                      << ap.i << " and " << ap.j);
            }
            
          }
        }
        
        if (n != num){
          io::messages.add("Wrong number of pairs in PERTATOMPAIR block.",
                  "In_Perturbation", io::message::error);
        }
        else if (_lineStream.fail()){
          io::messages.add("Bad line in PERTATOMPAIR block.",
                  "In_Perturbation", io::message::error);
        }
        
        if (!quiet)
          os << "\tEND\n";
        
      } // if block present
    } // PERTATOMPAIR
    
    { // PERTATOMPARAM
      
      DEBUG(10, "PERTATOM03 block");
      buffer = m_block["PERTATOM03"];
      
      if (buffer.size()){
        
        block_read.insert("PERTATOM03");
        io::messages.add("The PERTATOM03 block was renamed to PERTATOMPARAM.",
                "In_Perturbation", io::message::error);
      }
      
      buffer = m_block["PERTATOMPARAM"];
      if (buffer.size()){
        block_read.insert("PERTATOMPARAM");
        
        if (!quiet)
          os << "\tPERTATOMPARAM\n";
        DEBUG(7, "PERTATOMPARAM block");
        
        it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;
       
        int seq = 0, a_iac_int = 0, b_iac_int = 0; 
        unsigned int res = 0, a_iac = 0, b_iac = 0;
        double a_mass = 0.0, b_mass = 0.0, a_charge = 0.0, b_charge = 0.0;
        double lj_soft = 0.0, crf_soft = 0.0;
        std::string name;
        
        if (!quiet)
          os << "\t"
          << std::setw(5) << "seq"
          << std::setw(8) << "IAC(A)"
          << std::setw(10) << "mass(A)"
          << std::setw(10) << "charge(A)"
          << std::setw(8) << "IAC(B)"
          << std::setw(10) << "mass(B)"
          << std::setw(10) << "charge(B)"
          << std::setw(10) << "LJ(soft)"
          << std::setw(10) << "CRF(soft)"
          << "\n";
        
        for(n = 0; it != buffer.end() - 1; ++it, ++n){
          DEBUG(10, "\treading a line: " << n);
          
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> seq >> res >> name >> a_iac_int >> a_mass >> a_charge
                  >> b_iac_int >> b_mass >> b_charge
                  >> lj_soft >> crf_soft;
          
          if (_lineStream.fail()){
            io::messages.add("Bad line in PERTATOMPARAM block.",
                    "In_Perturbation", io::message::error);
          }
          
          --seq;
          --a_iac_int;
          --b_iac_int;
          
          // weight by input
          lj_soft *= param.perturbation.soft_vdw;
          crf_soft *= param.perturbation.soft_crf;
          
          if (seq < 0 || seq >= int(topo.num_solute_atoms())){
            io::messages.add("atom sequence number wrong in PERTATOMPARAM block",
                    "In_Perturbation", io::message::critical);
            return;
          }
          
          if (a_iac_int < 0 || b_iac_int < 0){
            os << "Problematic line: n=" << n+1 << " a_iac=" << a_iac_int+1 << " b_iac=" << b_iac_int+1 << std::endl;
            io::messages.add("integer atom code wrong in PERTATOMPARAM block",
                    "In_Perturbation", io::message::critical);
            return;
          }
          else{
            a_iac = (unsigned int)a_iac_int;
            b_iac = (unsigned int)b_iac_int;
          }	
          
          if (a_iac >= topo.atom_names().size() || b_iac >= topo.atom_names().size()){
            os << "Problematic line: n=" << n+1 << " a_iac=" << a_iac+1 << " b_iac=" << b_iac+1 << std::endl;
            io::messages.add("integer atom code wrong in PERTATOMPARAM block",
                    "In_Perturbation", io::message::critical);
            return;
          }

          topology::Perturbed_Atom atom(seq, a_iac, a_mass, a_charge,
                  topo.polarisability(seq), topo.damping_level(seq),
                  b_iac, b_mass, b_charge,
                  topo.polarisability(seq), topo.damping_level(seq),
                  lj_soft, crf_soft);
          
          DEBUG(10, "\tcreated an atom");
          
          if (topo.mass(seq) != atom.A_mass() ||
                  topo.iac(seq) != int(atom.A_IAC()) ||
                  topo.charge(seq) != atom.A_charge())
            warn = true;
          
          if (!quiet)
            os << "\t"
            << std::setw(5) << seq + 1
            << std::setw(5) << a_iac + 1
            << "   "
            << std::setw(10) << a_mass
            << std::setw(10) << a_charge
            << std::setw(5) << b_iac + 1
            << "   "
            << std::setw(10) << b_mass
            << std::setw(10) << b_charge
            << std::setw(10) << lj_soft
            << std::setw(10) << crf_soft
            << "\n";
          
          atom.exclusion() = topo.exclusion(seq);
          topo.exclusion(seq).clear();
          DEBUG(10, "\treplace the exclusions to perturbation");
          
          topology::excl_cont_t & ex = topo.exclusion();
          int seq2=0;
          
          for(topology::excl_cont_t::iterator eit=ex.begin(),
                  eto=ex.end(); eit!=eto; ++eit, ++seq2){
            if(eit->find_and_remove(seq)){
              atom.exclusion().insert(seq2);
            }
          }
          DEBUG(10, "\tadapted perturbed exclusions");
          
          atom.one_four_pair() = topo.one_four_pair(seq);
          topo.one_four_pair(seq).clear();
          DEBUG(10, "\treplaced the 14 interactions");
          
          topology::excl_cont_t & ofp = topo.one_four_pair();
          seq2=0;
          
          for(topology::excl_cont_t::iterator pit=ofp.begin(),
                  pito= ofp.end(); pit!=pito; ++pit, ++seq2){
            if(pit->find_and_remove(seq)){
              atom.one_four_pair().insert(seq2);
            }
          }
          DEBUG(10, "\tadapted 14 interactions");
          
          
          topo.perturbed_solute().atoms()[seq] = atom;
          
          assert(seq<int(topo.is_perturbed().size()));
          topo.is_perturbed()[seq] = true;
          
        }
        if (n != num){
          io::messages.add("Wrong number of bonds in PERTATOMPARAM block.",
                  "In_Perturbation", io::message::error);
        }
        else if (_lineStream.fail()){
          io::messages.add("Bad line in PERTATOMPARAM block.",
                  "In_Perturbation", io::message::error);
        }
        
        
        if (!quiet)
          os << "\tEND\n";
        
      } // if block present
      
    } // PERTATOMPARAM
    
    { // SCALEDINTERACTIONS
      
      buffer = m_block["SCALEDINTERACTIONS"];
      if (buffer.size()){
        block_read.insert("SCALEDINTERACTIONS");
        
        if(!param.perturbation.scaling){
          io::messages.add("Scaled interactions not turned on, ignoring SCALEDINTERACTIONS block.",
                  "InPerturbation", io::message::warning);
        }
        else{
          if (!quiet)
            os << "\tSCALED INTERACTIONS\n";
          
          it = buffer.begin() + 1;
          _lineStream.clear();
          _lineStream.str(*it);
          int num = 0, n = 0;
          _lineStream >> num;
          ++it;
          
          int i = 0, j = 0;
          double A = 0.0, B = 0.0;
          
          if (!quiet)
            os << "\t"
            << std::setw(10) << "group i"
            << std::setw(10) << "group j"
            << std::setw(10) << "scale A"
            << std::setw(10) << "scale B"
            << "\n";
          
          for(n = 0; it != buffer.end() - 1; ++it, ++n){
            _lineStream.clear();
            _lineStream.str(*it);
            _lineStream >> i >> j >> A >> B;
            
            if (_lineStream.fail()){
              io::messages.add("Bad line in SCALEDINTERACTIONS block.",
                      "In_Perturbation", io::message::error);
            }
            
            --i;
            --j;
            
            std::pair<int, int> energy_pair(i, j);
            std::pair<int, int> energy_pair2(j, i);
            
            std::pair<double, double> scale_pair(A, B);
            
            topo.energy_group_scaling()[energy_pair]=scale_pair;
            topo.energy_group_scaling()[energy_pair2]=scale_pair;
            
            if (!quiet)
              os << "\t"
              << std::setw(10) << i+1
              << std::setw(10) << j+1
              << std::setw(10) << A
              << std::setw(10) << B
              << std::endl;
            
          }
          
          if (n != num){
            io::messages.add("Wrong number of pairs in SCALEDINTERACTIONS block.",
                    "In_Perturbation", io::message::error);
          }
          else if (_lineStream.fail()){
            io::messages.add("Bad line in SCALEDINTERACTIONS block.",
                    "In_Perturbation", io::message::error);
          }
          if (!quiet)
            os << "\tEND\n";
        } // if scaling turned on
        
      } // if block present
      else{
        if(param.perturbation.scaling){
          io::messages.add("Scaling turned on but no SCALEDINTERACTIONS block.",
                  "In_Perturbation", io::message::error);
        }
      }
    } // SCALEDINTERACTIONS
    
    { // LAMBDADEP
      
      buffer = m_block["LAMBDADEP"];
      if (buffer.size()){
        block_read.insert("LAMBDADEP");
	io::messages.add("LAMBDADEP block in perturbation topology is deprecated. Use LAMBDAS block in the input file",
			 "InPerturbation", io::message::error);
      }
    } // LAMBDADEP
    
    { // PERTPOLPARAM
      DEBUG(10, "PERTPOLPARAM block");
      
      buffer = m_block["PERTPOLPARAM"];
      if (buffer.size()){
        block_read.insert("PERTPOLPARAM");
        
        if (!quiet)
          os << "\tPERTPOLPARAM\n";
        DEBUG(7, "PERTPOLPARAM block");
        
        it = buffer.begin() + 1;
        _lineStream.clear();
        _lineStream.str(*it);
        int num = 0, n = 0;
        _lineStream >> num;
        ++it;
        
        int seq = 0, res = 0;
        double a_pol = 0.0, b_pol = 0.0, a_lev = 0.0, b_lev = 0.0;
        std::string name;
        
        if (!quiet)
          os << "\t"
          << std::setw(5) << "seq"
          << std::setw(12) << "POL(A)"
          << std::setw(12) << "DAMPLEV(A)"
          << std::setw(12) << "POL(B)"
          << std::setw(12) << "DAMPLEV(B)"
          << "\n";
        
        for(n = 0; it != buffer.end() - 1; ++it, ++n){
          DEBUG(10, "\treading a line: " << n);
          
          _lineStream.clear();
          _lineStream.str(*it);
          _lineStream >> seq >> res >> name
                  >> a_pol >> a_lev
                  >> b_pol >> b_lev;
          
          if (_lineStream.fail()){
            io::messages.add("Bad line in PERTPOLPARAM block.",
                    "In_Perturbation", io::message::error);
          }
          
          --seq;
          if (seq < 0 || seq >= int(topo.num_solute_atoms())){
            io::messages.add("atom sequence number wrong in PERTPOLPARAM block",
                    "In_Perturbation", io::message::error);
            return;
          }
          
          if (!topo.is_perturbed(seq)) {
            std::ostringstream msg;
            msg << "Atom " << seq + 1 << " appears in the PERTPOLPARAM block but"
                    << " is not perturbed i.e. does not appear in the PERTATOMPARAM block.";
            io::messages.add(msg.str(), "In_Perturbation", io::message::error);
            return;
          }
          
          if (a_pol < 0.0 || b_pol < 0.0){
            io::messages.add("PERTPOLPARAM block: polarisability must be >= 0.0",
                    "In_Perturbation", io::message::error);
            return;
          }
          
          if (a_lev < 0.0 || b_lev < 0.0){
            io::messages.add("PERTPOLPARAM block: damping level must be >= 0.0",
                    "In_Perturbation", io::message::error);
            return;
          }
          
          topo.perturbed_solute().atom(seq).A_polarisability(a_pol/ math::four_pi_eps_i);
          topo.perturbed_solute().atom(seq).B_polarisability(b_pol/ math::four_pi_eps_i);
          
          topo.perturbed_solute().atom(seq).A_damping_level(a_lev * sqrt(math::four_pi_eps_i));
          topo.perturbed_solute().atom(seq).B_damping_level(b_lev * sqrt(math::four_pi_eps_i));
          
          topo.is_polarisable()[seq] = true;
          
          if (topo.polarisability(seq) !=
                  topo.perturbed_solute().atom(seq).A_polarisability() ||
                  topo.damping_level(seq) !=
                  topo.perturbed_solute().atom(seq).A_damping_level())
            warn = true;
          
          DEBUG(10, "\tassigned perturbed polarisation parameters to atom");
          
          if (!quiet)
            os << "\t"
            << std::setw(5) << seq + 1
            << std::setw(12) << topo.perturbed_solute().atom(seq).A_polarisability()* math::four_pi_eps_i
            << std::setw(12) << topo.perturbed_solute().atom(seq).A_damping_level()/sqrt(math::four_pi_eps_i) 
            << std::setw(12) << topo.perturbed_solute().atom(seq).B_polarisability()* math::four_pi_eps_i
            << std::setw(12) << topo.perturbed_solute().atom(seq).B_damping_level()/sqrt(math::four_pi_eps_i)
            << "\n";
        }
        if (n != num){
          io::messages.add("Wrong number of bonds in PERTPOLPARAM block.",
                  "In_Perturbation", io::message::error);
        }
        else if (_lineStream.fail()){
          io::messages.add("Bad line in PERTPOLPARAM block.",
                  "In_Perturbation", io::message::error);
        }
        
        
        if (!quiet)
          os << "\tEND\n";
        
      } // if block present
      
    } // PERTPOLPARAM
  } // end of if(!param.eds.eds)
  { // MPERTATOM
    DEBUG(10, "MPERTATOM block");

    buffer = m_block["MPERTATOM"];
    if (buffer.size()) {
      block_read.insert("MPERTATOM");

      if (!quiet)
        os << "MPERTATOM\n";
      DEBUG(7, "MPERTATOM block");

      it = buffer.begin() + 1;
      _lineStream.clear();
      _lineStream.str(*it);
      unsigned int numstates = 0, numatoms = 0, n = 0;
      //numstates = param.eds.numstates;

      _lineStream >> numatoms >> numstates;
      ++it;
      if (_lineStream.fail()){
        io::messages.add("Bad line (numatoms, numstates) in MPERTATOM block.",
                         "In_Perturbation", io::message::error);
        return;
      }

      if (numstates != param.eds.numstates) {
        std::ostringstream msg;
        msg << "Number of perturbed states given in perturbation topology ("
                << numstates << ") and input file (" << param.eds.numstates
                << ") do not match.";
        io::messages.add(msg.str(),
                "In_Perturbation", io::message::error);
        return;
      }
      
      // read in the name to identify the perturbation
      std::vector<std::string> identifier(numstates);
      _lineStream.clear();
      _lineStream.str(*it);
      for (unsigned int i = 0; i < identifier.size(); i++) {
        _lineStream >> identifier[i];
      }
      ++it;
      if (_lineStream.fail()) {
        io::messages.add("Bad line (PTNAME i.e. perturbation identifier name) in MPERTATOM block.",
                "In_Perturbation", io::message::error);
        return;
      }
      if (!quiet) {
        os << "\t" << std::setw(5) << "name";
        for (unsigned int i = 0; i < identifier.size(); i++) {
          os << std::setw(14) << identifier[i];
        }
        os << "\n";
      }
      int seq = 0;
      std::string name;
      std::vector<unsigned int> m_iac(numstates);
      std::vector<double> m_charge(numstates);
      double lj_soft = 0.0, crf_soft = 0.0;
      
      // prepare arrays
      topo.is_eds_perturbed().resize(topo.num_solute_atoms(), false);
      
      if (!quiet) {
        os << "\t"
                << std::setw(5) << "seq";
        for (unsigned int i = 0; i < numstates; i++) {
          os << std::setw(8) << "IAC(" << i << ")"
                  << std::setw(10) << "charge(" << i << ")";
        }
        os << std::setw(10) << "LJ(soft)"
                << std::setw(10) << "CRF(soft)"
                << "\n";
      }

      for (n = 0; it != buffer.end() - 1; ++it, ++n) {
        DEBUG(10, "\treading a line: " << n);
        DEBUG(10, "\tline content: " << *it);

        _lineStream.clear();
        _lineStream.str(*it);
        _lineStream >> seq >> name;
        for (unsigned int i = 0; i < numstates; i++) {
          _lineStream >> m_iac[i] >> m_charge[i];
          DEBUG(12,"\t iac = " << m_iac[i] << ", charge = " << m_charge[i]);
        }
        _lineStream >> lj_soft >> crf_soft;
        
        if (_lineStream.fail()) {
          io::messages.add("Bad line in MPERTATOM block.",
                  "In_Perturbation", io::message::error);
          return;
        }
        
        // weight by input
        lj_soft *= param.eds.soft_vdw;
        crf_soft *= param.eds.soft_crf;

        --seq;
        for (unsigned int i = 0; i < numstates; i++) {
          --m_iac[i];
          if (m_iac[i] < 0.0) {
            io::messages.add("Negative IAC in MPERTATOM block.",
                    "In_Perturbation", io::message::critical);
            return;
          }
        }

        if (seq < 0 || seq >= int(topo.num_solute_atoms())) {
          io::messages.add("atom sequence number wrong in MPERTATOM block",
                  "In_Perturbation", io::message::error);
          return;
        }

        if (!quiet) {
          os << "\t"
                  << std::setw(5) << seq + 1;
          for (unsigned int i = 0; i < numstates; i++) {
            os << std::setw(8) << m_iac[i] + 1
                    << std::setw(10) << m_charge[i];
          }
          os << std::setw(10) << lj_soft
                  << std::setw(10) << crf_soft
                  << "\n";
        }

        // create an eds perturbed atom
        topology::EDS_Perturbed_Atom atom(seq, m_iac, m_charge, lj_soft, crf_soft);
        DEBUG(10, "\tcreated an atom");
              
        
        atom.exclusion() = topo.exclusion(seq);
        topo.exclusion(seq).clear();
        DEBUG(10, "\treplace the exclusions to perturbation");
        
        topology::excl_cont_t & ex = topo.exclusion();
        int seq2=0;
        
        for(topology::excl_cont_t::iterator eit=ex.begin(),
                eto=ex.end(); eit!=eto; ++eit, ++seq2){
          if(eit->find_and_remove(seq)){
            atom.exclusion().insert(seq2);
          }
        }
        DEBUG(10, "\tadapted perturbed exclusions");
        
        atom.one_four_pair() = topo.one_four_pair(seq);
        topo.one_four_pair(seq).clear();
        DEBUG(10, "\treplaced the 1,4 interactions");
        
        topology::excl_cont_t & ofp = topo.one_four_pair();
        seq2=0;
        
        for(topology::excl_cont_t::iterator pit=ofp.begin(),
                pito= ofp.end(); pit!=pito; ++pit, ++seq2){
          if(pit->find_and_remove(seq)){
            atom.one_four_pair().insert(seq2);
          }
        }
        DEBUG(10, "\tadapted 1,4 interactions");
        
        topo.eds_perturbed_solute().atoms()[seq] = atom;
        assert(seq < int(topo.is_eds_perturbed().size()));
        topo.is_eds_perturbed()[seq] = true;
      }
      if (n != numatoms) {
        io::messages.add("Wrong number of perturbed atoms in MPERTATOM block.",
                         "In_Perturbation", io::message::error);
        return;
      }
      else if (_lineStream.fail()){
        io::messages.add("Bad line in MPERTATOM block.",
                         "In_Perturbation", io::message::error);
        return;
      }
           
      if (!quiet)
        os << "\tEND\n";
         
      DEBUG(10,"We have " << unsigned(topo.eds_perturbed_solute().atoms().size()) 
              << " eds perturbed atoms.");
    } // if block present
    
  } // MPERTATOM

  for(std::map<std::string, std::vector<std::string> >::const_iterator
	it = m_block.begin(),
	to = m_block.end();
      it != to;
      ++it){
    
    if (block_read.count(it->first) == 0 && it->second.size()){
      if(param.eds.eds){
        io::messages.add("block " + it->first + " not supported for eds!",
                "In_Perturbation (Topology)",
                io::message::error);
      }
      else{
        io::messages.add("block " + it->first + " not supported!",
                "In_Perturbation (Topology)",
                io::message::error);
      }
      
    }
  }
  
  if (!quiet)
    os << "END\n";

  if (warn) {
    io::messages.add("Some parameters in state A do not match the unperturbed topology.",
                     "In_Perturbation", io::message::warning);
  }  
  
  // and update the properties for lambda
  //topo.update_for_lambda();
  
}


