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

/* 
 * File:   JvalFuncs.h
 * Author: jallison
 *
 * Created on November 27, 2009, 12:50 PM
 */

#ifndef INCLUDED_UTILS_JVALFUNCS
#define INCLUDED_UTILS_JVALFUNCS

#include <vector>
#include <string>
#include "../gromos/Exception.h"
#include "../args/Arguments.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

namespace gcore {
    class System;
}
namespace args {
    class Arguments;
}
namespace utils {
    class PropertyContainer;
}

// define class for storing Jvalue data and the corresponding atoms and parameters
namespace utils {

    /**
     * Class JvalFuncs
     * Defines a data structure for storing Jvalues and related information
     * and some functions for SVD-fitting to Jvalues to get Karplus parameters.
     *
     * @class JvalFuncs
     * @author J. Allison
     */

    /**
     * Class JvalData
     * A class to store Jvalue data read in from a Jvalue specification file
     *
     * @class JvalData
     * @author J. Allison
     * @ingroup utils
     */
    class JvalData {
    public:

        // struct for storing Jval data
        struct jvalues {
            // atoms defining the dihedral angle
            int mol;
            int i;
            int j;
            int k;
            int l;
            // parameters
            double w; // weight factor
            double exp; // experimental jvalue

            jvalues() : mol(0), i(0), j(0), k(0), l(0) {
            }

            jvalues(const jvalues & jval) : mol(jval.mol), i(jval.i), j(jval.j),
            k(jval.k), l(jval.l), w(jval.w), exp(jval.exp) {
            }

            jvalues & operator=(const jvalues & jval) {
                mol = jval.mol;
                i = jval.i;
                j = jval.j;
                k = jval.k;
                l = jval.l;
                w = jval.w;
                exp = jval.exp;
                return *this;
            }
        };

        /**
         * Constructor
         */
        JvalData() {
        }

        /**
         *  JvalData copy constructor
         */
        JvalData(const JvalData &jvaldata) {
            m_data = jvaldata.m_data;
        }

        /**
         *  JvalData deconstructor
         */
        ~JvalData() {
        }

        /**
         * const accessor to jvalues data
         */
        const std::vector<jvalues> & data() const {
            return m_data;
        }

        /**
         * accessor to jvalues data
         */
        std::vector<jvalues> & data() {
            return m_data;
        }

    private:
        std::vector<jvalues> m_data;

    };

    /**
     * Class JvalWeights
     * A class to store weights for individual frames of a trajectory
     *
     * @class JvalWeights
     * @author J. Allison
     * @ingroup utils
     */
    class JvalWeights {
    public:


        // struct for storing weights (for frames of a trajectory)

        struct weights {
            int frame;
            double weight;

            weights() : frame(0), weight(0) {
            }

            weights(const weights & w) : frame(w.frame), weight(w.weight) {
            }

            weights & operator=(const weights & w) {
                frame = w.frame;
                weight = w.weight;
                return *this;
            }

        };

        /**
         * Constructor
         */
        JvalWeights() {
        }

        /**
         *  JvalWeights copy constructor
         */
        JvalWeights(const JvalWeights &weights) {
            m_weights = weights.m_weights;
        }

        /**
         *  JvalWeights deconstructor
         */
        ~JvalWeights() {
        }

        /**
         * const accessor to weight data
         */
        const std::vector<weights> & wdata() const {
            return m_weights;
        }

        /**
         * accessor to weights data
         */
        std::vector<weights> & wdata() {
            return m_weights;
        }

    private:

        std::vector<weights> m_weights;

    };

    /**
     * Class JvalFuncs
     * A class of functions for fitting to and calculating Jvalues and Karplus parameters
     *
     * @class JvalFuncs
     * @author J. Allison
     * @ingroup utils
     */
    class JvalFuncs {
        /**
         * copy constructor
         * not implemented
         */
        JvalFuncs(const JvalFuncs&);

        /**
         * operator =
         * not implemented
         */
        JvalFuncs & operator=(const JvalFuncs&);

        public:

        /**
         * JvalFuncs Constructor
         */
        JvalFuncs(gcore::System &sys, args::Arguments &args);

        /**
         * JvalFuncs Destructor
         */
        ~JvalFuncs();

        // Methods

        // function to read in jvalue data and get one global delta
        void read_jval(std::vector<std::string> buffer, const gcore::System &sys,
                std::vector<utils::JvalData::jvalues> &jval,
                utils::PropertyContainer &props, double &delta);

        // function to read weights for individual frames from file
        void read_weights(std::vector<std::string> buffer, std::vector<JvalWeights::weights> &weight_data);

        // function to calculate Jvalue coefficients for the SVD fit
        void calc_coef(const gcore::System &sys, utils::PropertyContainer &fit_props,
        gsl_matrix *coef_mat, int njval, double w, double delta);

        // my version of this function
        void fill_jvalvec(const std::vector<utils::JvalData::jvalues> &J, gsl_vector *v,
        int njval);

        // calculate Q value (goodness of fit)
        double calc_Q(gsl_vector *calc, gsl_vector *expt);

        /**
         * @struct Exception
         * Throws an exception if something is wrong
         */
        struct Exception : public gromos::Exception {

            /**
             * @exception If called says JvalFuncs, followed by the argument
             * @param what The string that is thrown
             */
            Exception(const std::string & what) :
            gromos::Exception("JvalFuncs", what) {
            }
        };

    }; // end class JvalFuncs

}


#endif	/* _JVALFUNCS_H */

