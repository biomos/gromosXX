/**
 * @file in_ifp.cc
 * read in an Interaction Function Parameter file
 */

#include <stdheader.h>

#include <interaction/interaction_types.h>

#include <io/instream.h>
#include <io/blockinput.h>

#include "in_ifp.h"

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE forcefield


void io::In_IFP
::read_harmonic_bonds(std::vector<interaction::bond_type_struct> &b)
{
  
  DEBUG(10, "BONDTYPECODE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDTYPECODE"];

  if (buffer.size()==0)
    io::messages.add("BONDTYPE block not found",
		     "In_IFP",
		     io::message::error);

  for (it = buffer.begin() + 1; 
       it != buffer.end() - 1; ++it) {
    
    int i;
    double k, r;
    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> i >> k >> r;
      
    if (_lineStream.fail()){
      io::messages.add("Bad line in BONDTYPECODE block: " + *it,
		       "In_IFP", io::message::error);
    }
    /*
    if (! _lineStream.eof()){
      io::messages.add("eof not reached in BONDTYPECODE block: " + *it,
		       "In_IFP", io::message::warning);
    }
    */
    // we are reading into harmonic bond term, so convert k
    k *= 2 * r * r;
    
    // and add...
    b.push_back(interaction::bond_type_struct(k, r));
  }
}

void io::In_IFP
::read_g96_bonds(std::vector<interaction::bond_type_struct> &b)
{
  DEBUG(10, "BONDTYPECODE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDTYPECODE"];

  if (buffer.size()==0)
    io::messages.add("BONDTYPE block not found",
		     "In_IFP",
		     io::message::error);

  for (it = buffer.begin() + 1; 
       it != buffer.end() - 1; ++it) {

    int i;
    double k, r;

    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> i >> k >> r;
      
    if (_lineStream.fail()){
      io::messages.add("Bad line in BONDTYPECODE block: " + *it,
		       "In_IFP", io::message::error);
    }
    /*
    if (! _lineStream.eof()){
      io::messages.add("eof not reached in BONDTYPECODE block: " + *it,
		       "In_IFP", io::message::warning);
    }
    */
    // and add...
    b.push_back(interaction::bond_type_struct(k, r));
  }

}

void io::In_IFP
::read_angles(std::vector<interaction::angle_type_struct> &b)
{
  DEBUG(10, "BONDANGLETYPECOD block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["BONDANGLETYPECOD"];

  if (buffer.size()==0)
    io::messages.add("BONDANGLETYPECOD block not found",
		     "In_IFP",
		     io::message::error);

  for (it = buffer.begin() + 1; 
       it != buffer.end() - 1; ++it) {
    
    int i;
    double k, cos0;

    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> i >> k >> cos0;
      
    if (_lineStream.fail()){
      io::messages.add("Bad line in BONDANGLETYPECOD block: " + *it,
		       "In_IFP", io::message::error);
    }
    /*
    if (! _lineStream.eof()){
      io::messages.add("eof not reached in BONDANGLETYPECOD block: " + *it,
		       "In_IFP", io::message::warning);
    }
    */
    // and add...
    b.push_back(interaction::angle_type_struct(k, cos(cos0 * 2 * math::Pi / 360.0)));
  }
}

void io::In_IFP
::read_improper_dihedrals(std::vector<interaction::improper_dihedral_type_struct> &imp)
{
  DEBUG(10, "IMPDIHEDRALTYPEC block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["IMPDIHEDRALTYPEC"];

  if (buffer.size()==0)
    io::messages.add("IMPDIHEDRALTYPEC block not found",
		     "In_IFP",
		     io::message::error);

  for (it = buffer.begin() + 1; 
       it != buffer.end() - 1; ++it) {

    int i;
    double k, q0;

    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> i >> k >> q0;

    if (_lineStream.fail()){
      io::messages.add("Bad line in IMPDIHEDRALTYPEC block: " + *it,
		       "In_IFP", io::message::error);
    }
    /*
    if (! _lineStream.eof()){
      io::messages.add("eof not reached in IMPDIHEDRALTYPEC block: " + *it,
		       "In_IFP", io::message::warning);
    }
    */
    
    // and add...
    imp.push_back(interaction::improper_dihedral_type_struct(k*180*180/math::Pi/math::Pi,
							   q0 * math::Pi / 180.0));
  }
}

void io::In_IFP
::read_dihedrals(std::vector<interaction::dihedral_type_struct> &d)
{
  DEBUG(10, "DIHEDRALTYPECODE block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  buffer = m_block["DIHEDRALTYPECODE"];

  if (buffer.size()==0)
    io::messages.add("DIHEDRALTYPECODE block not found",
		     "In_IFP",
		     io::message::error);

  for (it = buffer.begin() + 1; 
       it != buffer.end() - 1; ++it) {
    
    int i;
    double k, pd;
    int m;

    _lineStream.clear();
    _lineStream.str(*it);
      
    _lineStream >> i >> k >> pd >> m;
      
    if (_lineStream.fail()){
      io::messages.add("Bad line in DIHEDRALTYPECODE block: " + *it,
		       "In_IFP", io::message::error);
    }
    /*
    if (! _lineStream.eof()){
      io::messages.add("eof not reached in DIHEDRALTYPECODE block: " + *it,
		       "In_IFP", io::message::warning);
    }
    */

    // and add...
    d.push_back(interaction::dihedral_type_struct(k, pd, m));
  }
}

void io::In_IFP
::read_lj_parameter(std::vector<std::vector
		    <interaction::lj_parameter_struct> > 
		    & lj_parameter)
{
  DEBUG(10, "SINGLEATOMLJPAIR block");

  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;

  { // SINGLEATOMLJPAIR

    buffer = m_block["SINGLEATOMLJPAIR"];

    if (buffer.size() < 3)
      io::messages.add("SINGLEATOMLJPAIR block not found",
		       "In_IFP",
		       io::message::error);
    
    it = buffer.begin() + 1;
    _lineStream.clear();
    _lineStream.str(*it);

    int atom_types;
    _lineStream >> atom_types;
    ++it;
    
    DEBUG(11, "resize lj_parameter to " << atom_types);
    lj_parameter.resize(atom_types);
    for(int i=0; i<atom_types; ++i)
      lj_parameter[i].resize(atom_types);
    
    std::vector<double> sc6(atom_types), sc12[3], scs6(atom_types), scs12(atom_types);
    sc12[0].resize(atom_types);
    sc12[1].resize(atom_types);
    sc12[2].resize(atom_types);
    
    std::vector<std::vector<int> >    pl(atom_types);
	for(int ii=0; ii<atom_types; ++ii){
      pl[ii].resize(atom_types);
	}
    
    // get the whole thing in one...
    std::string ljblock;
    io::concatenate(it, buffer.end()-1, ljblock);

    _lineStream.clear();
    _lineStream.str(ljblock);

    int num;
    std::string s;
    
    DEBUG(10, "read in sc6, sc12, scs6, scs12");
    for(int n=0; n<atom_types; n++){

      DEBUG(15, "atom type: " << n);

      _lineStream >> num >> s >> sc6[n]>> sc12[0][n] >> sc12[1][n] >> sc12[2][n]
		  >> scs6[n] >> scs12[n];

      if (_lineStream.fail()){
	io::messages.add("Bad line in SINGLEATOMLJPAIR block",
			 "In_IFP", io::message::error);
      }
      /*
      if (! _lineStream.eof()){
	io::messages.add("eof not reached in SINGLEATOMLJPAIR block",
			 "In_IFP", io::message::warning);
      }
      */

      if(num != n+1){
	io::messages.add("atom types not sequential in SINGLEATOMLJPAIR block",
			 "In_IFP", io::message::error);
      }

      DEBUG(15, "connection matrix");
      for(int k=0; k<atom_types; k++)
	_lineStream >> pl[n][k];

      if(_lineStream.fail()){
	io::messages.add("pair matrix error in SINGLEATOMLJPAIR block (file corrupt)",
			 "In_IFP", io::message::error);
      }

      DEBUG(15, "combine data");
      for(int k=0; k<=n; k++){
	const double c6 = sc6[n]               * sc6[k];
	const double c12 = sc12[pl[n][k]-1][n] * sc12[pl[k][n]-1][k];
	const double cs6 = scs6[n]             * scs6[k];
	const double cs12 = scs12[n]           * scs12[k];

	interaction::lj_parameter_struct s(c6, c12, cs6, cs12);
	DEBUG(20, "add pair " << n << " : " << k);
	lj_parameter[n][k] = s;
	DEBUG(20, "add pair " << k << " : " << n);
	lj_parameter[k][n] = s;
      }
    }
  } // SINGLEATOMLJPAIR
  { // MIXEDATOMLJPAIR block
    // this one does not have to be there
    if(m_block.count("MIXEDATOMLJPAIR")){
      
      DEBUG(8, "MIXEDATOMLJPAIR block");
      
      buffer = m_block["MIXEDATOMLJPAIR"];
      
      for(it = buffer.begin() + 1; it != buffer.end() - 1; ++it){

	_lineStream.clear();
	_lineStream.str(*it);
	
	int i, j;
	double c6, c12, cs6, cs12;

	_lineStream >> i >> j >> c12 >> c6 >> cs12 >> cs6;

	if(_lineStream.fail()){
	}

	interaction::lj_parameter_struct s(c6, c12, cs6, cs12);
	lj_parameter[i-1][j-1] = s;
	lj_parameter[j-1][i-1] = s;
	
      }
    } // MIXEDATOMLJPAIR
  }

}

