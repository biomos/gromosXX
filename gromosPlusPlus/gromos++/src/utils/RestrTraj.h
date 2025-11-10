/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

// utils_RestrTraj.h

#ifndef INCLUDED_UTILS_RESTRTRAJ
#define INCLUDED_UTILS_RESTRTRAJ

#include <string>
#include <vector>

#include "../gromos/Exception.h"

namespace utils{
  class Time;
}

namespace utils{

  class RestrTraj_i;
  /**
   * Class RestrTraj
   * Defines an instream that can read the GROMOS09 restraint trajectory file.
   *
   * The instream can (currently) handle the TIMESTEP and JVALUERESEPS blocks
   *
   * @class RestrTraj
   * @author J. Allison, N. Schmid
   * @sa 
   * @todo Add more block
   */

  /**
   * Class JValueRestrData
   * A class to store restraint data read in from a restraint trajectory (*.trs)
   *
   * @class JValueRestrData
   * @author J. Allison, N. Schmid
   * @ingroup utils
   */
  class JValueRestrData{
    public:
      // struct for storing J-value LE potential data
      struct JValueEps {
        unsigned int i, j, k, l;
        std::vector<double> epsilon;
        JValueEps() : i(0), j(0), k(0), l(0) {}
        JValueEps(const JValueEps & jeps) : i(jeps.i), j(jeps.j), k(jeps.k),
                l(jeps.l), epsilon(jeps.epsilon) {
        }
        JValueEps & operator=(const JValueEps & jeps) {
          i = jeps.i; j = jeps.j; k = jeps.k; l = jeps.l;
          epsilon = jeps.epsilon;
          return *this;
        }
      };
    /**
     * Constructor
     */
    JValueRestrData() {}
    /**
     *  RestrData copy constructor
     */
    JValueRestrData(const JValueRestrData &jvaluerestrdata) {
      m_data = jvaluerestrdata.m_data;
    }
    /**
     *  RestrData deconstructor
     */
    ~JValueRestrData() {}
    /**
     * const accessor to jvalue eps data
     */
    const std::vector<JValueEps> & data() const {
      return m_data;
    }
    /**
     * accessor to jvalue eps data
     */
    std::vector<JValueEps> & data() {
      return m_data;
    }
  private:
    std::vector<JValueEps> m_data;
  };

    /**
   * Class XrayRestrData
   * A class to store restraint data read in from a restraint trajectory (*.trs)
   *
   * @class XrayRestrData
   * @author N. Schmid
   * @ingroup utils
   */
  class XrayRestrData{
    public:
      /**
       * struct for storing XrayState data like R-values and scaling constants
       */
      struct XrayState {
        double scale_inst, scale_free_inst, scale_avg, scale_free_avg;
        double r_inst, r_free_inst, r_avg, r_free_avg;
        XrayState() : scale_inst(0.0), scale_free_inst(0.0), 
                      scale_avg(0.0), scale_free_avg(0.0),
                      r_inst(0.0), r_free_inst(0.0),
                      r_avg(0.0), r_free_avg(0.0) {}
        XrayState(const XrayState & r) : scale_inst(r.scale_inst),scale_free_inst(r.scale_free_inst),
                scale_avg(r.scale_avg), scale_free_avg(r.scale_free_avg),
                r_inst(r.r_inst), r_free_inst(r.r_free_inst),
                r_avg(r.r_avg), r_free_avg(r.r_free_avg) {
        }
        XrayState & operator=(const XrayState & r) {
          scale_inst = r.scale_inst;
          scale_avg = r.scale_avg;
          scale_free_inst = r.scale_free_inst;
          scale_free_avg = r.scale_free_avg;
          r_inst = r.r_inst;
          r_avg = r.r_avg;
          r_free_inst = r.r_free_inst;
          r_free_avg = r.r_free_avg;
          return *this;
        }
      };
    /**
     * Constructor
     */
    XrayRestrData() {}
    /**
     *  RestrData copy constructor
     */
    XrayRestrData(const XrayRestrData &data) {
      m_state = data.m_state;
    }
    /**
     *  RestrData deconstructor
     */
    ~XrayRestrData() {}
    /**
     * const accessor to state
     */
    const XrayState & state() const {
      return m_state;
    }
    /**
     * accessor to state
     */
    XrayState & state() {
      return m_state;
    }
  private:
    XrayState m_state;
  };
  
  class RestrTraj{
    /**
     * pIMPL
     */
    RestrTraj_i *d_this;
    /**
     * skip
     */
    int d_skip;
    /**
     * stride
     */
    int d_stride;
    /**
     * end of file during striding (or skipping)
     */
    int d_stride_eof;

    /**
     * copy constructor
     * not implemented
     */
    RestrTraj(const RestrTraj&);
    /**
     * operator =
     * not implemented
     */
    RestrTraj &operator=(const RestrTraj&);
    
  public:
    /**
     * Constructor
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    RestrTraj(int skip=0, int stride=1);
    /**
     * Constructor
     * @param name   : trajectory file
     * @param skip   : skip frames
     * @param stride : take only every n-th frame
     */
    RestrTraj(const std::string &name, int skip=0, int stride=1);
    /**
     * Destructor
     */
    ~RestrTraj();

    // Methods
    /**
     * open a trajectory file
     * skip and stride continue
     */
    void open(const std::string &name);

    /**
     * read the trajectory
     */
    void read();
    /**
     * close
     */
    void close();
    /**
     * read a jvalue restraints frame
     */
    RestrTraj &operator>>(utils::JValueRestrData &jvaluerestrdata);
    /**
     * get the time information from a frame
     */
    RestrTraj &operator>>(utils::Time &time);
    /**
     * read a xray restraints frame
     */
    RestrTraj &operator>>(utils::XrayRestrData &xrayrestrdata);

    // Accessors
    /**
     * title (of the trajectory)
     */
    std::string title()const;
    /**
     * name?
     */
    std::string name()const;
    /**
     * end of file
     */
    bool eof()const;

    /**
     * take every n-th frame
     * starts with frame 0, then n, 2n, 3n, ...
     */
    int stride()const;
    /**
     * skip n frames
     * after skipping, skip will be 0.
     * to skip again, just set it again...
     */
    int skip()const;

    /**
     * set stride
     */
    void stride(int stride);
    /**
     * set skip
     */
    void skip(int skip);

    /**
     * encountered end of file
     * during striding
     */
    bool stride_eof()const;

    //Exceptions
    /**
     * Exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("RestrTraj", what_arg){}
    };
  };
}
#endif
