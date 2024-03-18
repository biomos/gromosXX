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

/**
 * @file tf_rdc_restraint_interaction.h
 * tensor-free RDC restraining
 */
#ifndef TF_RDC_RESTRAINT_INTERACTION_H
#define	TF_RDC_RESTRAINT_INTERACTION_H

#include "../../math/random.h"

namespace interaction
{
    /**
    * @class TF_RDC_Restraint_Interaction
    * calculates tensor-free RDC restraining interaction
    */
    // All constant variables are defined in GROMOS units, see book
    // math::eps0_i is defined as 1/\epsilon_0 = 1.7459 *10^3 (kJ nm)/(e^2 mol)
    // math::spd_l is defined as c = 2.9979 *10^5 nm/ps
    // const double m0 = math::eps0_i / (math::spd_l * math::spd_l); //1.9426 *10^-8 (kJ ps^2)/(e^2 mol nm)
    // math::h_bar is defined as  hbar*N_A = 6.351 *10^-2 (kJ ps)/mol
    // const double h_planck = math::h_bar * 2.0 * math::Pi; // 0.399 (kJ ps)/mol
    // returns D_c (rdc_max)  in units of ps^-1
    inline double D_c(std::vector<topology::tf_rdc_restraint_struct>::iterator it){
        return - (math::eps0_i * math::h_bar * it->gyri * it->gyrj) / (
            pow(it->normalisation_distance,3.0)*pow(2.0 * math::Pi * math::spd_l, 2));
    }
    class TF_RDC_Restraint_Interaction :
    public Interaction
    {
    public:
        /**
         * Constructor.
         */
        TF_RDC_Restraint_Interaction() : Interaction("TFRDCRestraint") {}

        /**
        * Destructor.
        */
        virtual ~TF_RDC_Restraint_Interaction() {}

        /**
        * init
        */
        virtual int init(topology::Topology &topo,
            configuration::Configuration &conf,
            simulation::Simulation &sim,
            std::ostream &os = std::cout,
            bool quiet = false);

        /**
        * calculate the interactions.
        */
        virtual int calculate_interactions(topology::Topology & topo,
            configuration::Configuration & conf,
            simulation::Simulation & sim);

        /**
         * random number generator
         */
        math::RandomGenerator * m_rng;

  };
} // interaction


#endif	/* TF_RDC_RESTRAINT_INTERACTION_H */
