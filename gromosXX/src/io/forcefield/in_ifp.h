// gio_InParameter.h

#ifndef INCLUDED_IN_IFP_H
#define INCLUDED_IN_IFP_H

#include <gromosXX/io/ifp.h>

namespace io{

  /**
   * Class In_IFP
   * defines an instream that can read a GROMOS96 ifp-file
   *
   * The GROMOS96 ifp file is read in and stored in a GromosForceField
   * format. This means that vdw-parameters are already calculated to 
   * the individual pairs, taking all special tables in the manual into 
   * account'
   *
   * @class In_IFP
   * @author B.C. Oostenbrink
   */
  class In_IFP : public GInStream, public IFP {
  public:
    In_IFP(std::istream& is) : GInStream(is) { readStream(); }
    virtual ~In_IFP() {}
    
  private:
    
  };
} // io

#endif
