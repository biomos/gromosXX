/**
 * @file mopac_worker.cc
 * interface to the MOPAC software package
 */

#include "../../../stdheader.h"

#include "../../../algorithm/algorithm.h"
#include "../../../topology/topology.h"
#include "../../../simulation/simulation.h"
#include "../../../configuration/configuration.h"

#include "../../../interaction/interaction.h"
#include "../../../util/timing.h"
#include "../../../util/system_call.h"
//#include "../../../math/periodicity.h"

// special interactions
#include "qm_storage.h"
#include "mm_atom.h"
#include "qm_worker.h"
#include "mopac_worker.h"
#include "../../../util/system_call.h"
#include "../../../io/blockinput.h"
#include "../../../util/debug.h"

#undef MODULE
#undef SUBMODULE
#define MODULE interaction
#define SUBMODULE special
#define MAXPATH 10240


double interaction::MOPAC_Worker::pointchg_pot(const int j,
                const std::vector<interaction::MM_Atom> & mm_atoms,
                const topology::Topology & topo,
                const configuration::Configuration &conf) {
        double x=0.0 ;
        for (std::vector<interaction::MM_Atom>::const_iterator
                it=mm_atoms.begin(),to=mm_atoms.end();it!=to;++it  ){
            math::Vec dR=conf.current().pos(it->index)-conf.current().pos(j);
            x+= (it->charge /abs(dR));
        }
return x*33.2063707338; // 1E10*e^2*NA/(4*pi*eps0*4184)
}


math::Vec interaction::MOPAC_Worker::pointchg_force(const int j,
                                                    std::vector<double> qm_chg,
                                               const topology::Topology & topo,
                                               const configuration::Configuration &conf) {
    math::Vec rij;
    math::Vec fij;
    fij[0]=0.0;
    fij[1]=0.0;
    fij[2]=0.0;
    unsigned int pi=0;
    for (std::set<topology::qm_atom_struct>::const_iterator
            it=topo.qm_zone().begin(),to=topo.qm_zone().end();it!=to;++it,pi++){
            rij=conf.current().pos(j)-conf.current().pos(it->index);
            fij+=(rij/abs(rij))*topo.charge(j)*qm_chg[pi]/(abs2(rij)) ;

    }
    return fij*41.84; //conversion from kcal/(mol*A) -> kJ/(mol*nm) - thermodynamic calorie
}


int interaction::MOPAC_Worker::run_QM(topology::Topology& topo,
                                     configuration::Configuration& conf,
                                     simulation::Simulation& sim,
                                     const math::VArray & qm_pos,
                                     const std::vector<interaction::MM_Atom> & mm_atoms,
                                     interaction::QM_Storage& storage,
                                     interaction::QM_Storage& LA_storage,
                                     const configuration::Configuration & qmmm_conf){
    {

        this->input_file = sim.param().qmmm.mopac.input_file;
        this->output_file = sim.param().qmmm.mopac.output_file;
        this->output_gradient_file = sim.param().qmmm.mopac.output_gradient_file;
        this->header_file = sim.param().qmmm.mopac.input_header;
        this->molin_file = sim.param().qmmm.mopac.molin_file;
        std::ofstream inp(input_file.c_str());
        if (!inp.is_open()) {
            io::messages.add("Could not create input file for MOPAC at the location "
                             + input_file, "MOPAC_Worker", io::message::critical);
            return 1;
        }
        bool verbose = false;
        // get the number of point charges
        unsigned int num_charge = mm_atoms.size();
        for (unsigned int i = 0; i < mm_atoms.size(); ++i) {
            if (topo.is_polarisable(mm_atoms[i].index)) {
                ++num_charge;
            }
        }
        std::ostringstream oss;
        std::string header(this->header_file);
        oss.str("");
        oss.clear();
 //       oss << topo.qm_mm_pair().size();

        // write header
        inp << header << std::endl;

        // write QM zone
        std::ofstream molin(this->molin_file.c_str());
        molin << std::endl;
        molin << topo.qm_zone().size() << "  " << topo.qm_mm_pair().size() << std::endl;
        double len_to_qm = 1.0 / sim.param().qmmm.unit_factor_length;
        double ptchg = 0.0;
        unsigned int pi = 0;
        for (std::set<topology::qm_atom_struct>::const_iterator
             it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it, ++pi) {
            inp.setf(std::ios::fixed, std::ios::floatfield);
            inp.precision(8);
            ptchg = 0.0;
            inp << std::setw(2) << std::left << it->atomic_number;
            molin.setf(std::ios::fixed, std::ios::floatfield);
            molin.precision(8);
            molin << std::setw(2) << std::left << it->atomic_number;
            for (unsigned int i = 0; i < 3; ++i) {
                inp << std::setw(14) << std::right << qm_pos(pi)(i) * len_to_qm << " +1";
                molin << std::setw(14) << std::right << qm_pos(pi)(i) * len_to_qm;
            }
            ptchg = pointchg_pot(it->index, mm_atoms, topo, conf);
            molin << std::setw(14) << std::right << ptchg;
            inp << std::endl;
            molin << std::endl;
        }

        pi = 0;

// append list of link-atoms to the QM coordinates
//   the order of the capping H atoms follows the same
//  order as the QM atoms in the "qmmm" file.
//
        math::Vec posCap, dR;
        unsigned int m1, q1;
        const double rch = sim.param().qmmm.cap_dist; //Carbon-Hydrogen bond length

        if (verbose) {
            std::cout << "number of link atoms: " << topo.qm_mm_pair().size() << std::endl;
        }
        for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                     it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it, ++pi) {
            q1 = it->first;
            m1 = it->second;
            dR = conf.current().pos(m1) - conf.current().pos(q1);
            posCap = (rch / abs(dR)) * dR + conf.current().pos(q1);
            inp << std::setw(5) << std::left << "1 ";  //link_atoms[i];
            molin << std::setw(5) << std::left << "1 ";  //link_atoms[i];
            for (unsigned int i = 0; i < 3; ++i) {
                inp << std::setw(14) << std::right << posCap(i) * len_to_qm << " +1";
                molin << std::setw(14) << std::right << posCap(i) * len_to_qm;
            }
            ptchg = pointchg_pot(q1, mm_atoms, topo, qmmm_conf);
            molin << std::setw(14) << std::right << ptchg;
            inp << std::endl;
            molin << std::endl;
        }
     molin.close();
    }


    // the termination line. just a couple of zeros
    /*
    inp << std::setw(2) << std::left << 0;
    for (unsigned int i = 0; i < 3; ++i) {
        inp << std::setw(14) << std::right << 0.0 << " 0";
    }
     */
    //write list of link atoms

    std::vector<double > qm_chg;
    double chg_to_qm = 1.0 / sim.param().qmmm.unit_factor_charge;
    // starting MOPAC

    std::ostringstream command_to_launch;
    command_to_launch << sim.param().qmmm.mopac.binary << " " << this->input_file << " >  "
                      << this->output_file << " 2> mopac.err" ;
    system(command_to_launch.str().c_str());
    /*
    if (result != 0) {
        std::ostringstream msg;
        msg << "MOPAC failed with code " << result << ". See output file "
            << sim.param().qmmm.mopac.output_file << " for details.";
        io::messages.add(msg.str(), "MOPAC_Worker", io::message::error);
        return 1;
    }
*/
    // read output gradient
    std::ifstream output(this->output_file.c_str());
    // std::ifstream output(sim.param().qmmm.mopac.output_file.c_str());
    if (!output.is_open()) {
        io::messages.add("Cannot open MOPAC output file", "MOPAC_Worker",
                         io::message::error);
        return 1;
    }
    while(1) {
        std::string dummy;
        double charge = 0.0;
        std::string line;
        std::getline(output, line);
        if (output.eof() || output.fail()) {
            break;
        }
        // get gradients
        math::Vec gradient;
        if (line.find("FINAL HEAT OF FORMATION") != std::string::npos){
            double energy;
            std::istringstream is(line);
            is >> dummy >> dummy >> dummy >> dummy >> dummy >> energy;
            storage.energy = energy * sim.param().qmmm.unit_factor_energy;
        }
        else if (line.find("PARAMETER") != std::string::npos) {
            // read atoms
            for (std::set<topology::qm_atom_struct>::const_iterator
                         it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
                for (unsigned int ixyz = 0; ixyz < 3; ixyz++) {
                    std::getline(output, line);
                    std::istringstream is(line);
                    is >> dummy >> dummy >> dummy >> dummy >> dummy >>
                       dummy >> gradient[ixyz];
                    if (output.fail()) {
                        std::ostringstream msg;
                        msg << "Failed to read charge line of " << (it->index + 1) << "th atom.";
                        io::messages.add(msg.str(), "MOPAC_Worker", io::message::error);
                        return 1;
                    }
                }
                DEBUG(5, "\t grad " << gradient[0] <<
                                    gradient[1] << "  " << gradient[2] << std::endl);
                storage.force(it->index) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                                                        sim.param().qmmm.unit_factor_length);
            }//read partial charges on QM-atoms required for force on mm-pointcharges
            //read now the forces on the capping atom
            for (unsigned int i = 0; i < topo.qm_mm_pair().size(); i++) {
                for (unsigned int ixyz = 0; ixyz < 3; ixyz++) {
                    std::getline(output, line);
                    std::istringstream is(line);
                    is >> dummy >> dummy >> dummy >> dummy >> dummy >>
                       dummy >> gradient[ixyz];
                }
                if (output.fail()) {
                    std::ostringstream msg;
                    msg << "Failed to read charge line of " << (i + 1) << "th linkatom.";
                    io::messages.add(msg.str(), "MOPAC_Worker", io::message::error);
                    return 1;
                }
                LA_storage.force(i) = gradient * (-(sim.param().qmmm.unit_factor_energy) /
                                                  sim.param().qmmm.unit_factor_length);
            } // for atoms
        }
        else if (line.find("ELECS")  != std::string::npos ) {
            for (std::set<topology::qm_atom_struct>::const_iterator
                         it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
                std::getline(output, line);
                std::istringstream is(line);
                is >> dummy >> dummy >> charge;
                qm_chg.push_back(charge);
                if (is.fail()) {
                    std::ostringstream msg;
                    msg << "Failed to parse charge line of " << (it->index + 1) << "th atom.";
                    io::messages.add(msg.str(), "MOPAC_Worker", io::message::error);
                    return 1;
                }
            }

            for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator
                         it = topo.qm_mm_pair().begin(); it != topo.qm_mm_pair().end(); ++it) {
                std::getline(output, line);
                std::istringstream is(line);
                is >> dummy >> dummy >> charge;
                qm_chg.push_back(0.0); //cap atom should not act on MM-atoms
               // qm_chg.push_back(charge);
                if (is.fail()) {
                    std::ostringstream msg;
                    msg << "Failed to parse charge line of " << it->first << " link atom.";
                    io::messages.add(msg.str(), "MOPAC_Worker", io::message::error);
                    return 1;
                }

            }
            }
        }
    output.close();
   //gradient on MM-pointcharges caused by QM-charges
    for (std::vector<interaction::MM_Atom>::const_iterator
        it=mm_atoms.begin(),to=mm_atoms.end();it!=to;++it  ){
        math::Vec gradient;
        gradient=pointchg_force(it->index,qm_chg,topo,qmmm_conf);
        DEBUG(5,"\t grad " << gradient[0] <<
                            gradient[1] << "  " << gradient[2] << std::endl);
        storage.force(it->index) += gradient ;
    }
        //print out the gradient
  /*
  for (std::set<topology::qm_atom_struct>::const_iterator
        it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
    std::cout << "index :   " << it->index << "  " << storage.force(it->index)[0] << "  "
            << storage.force(it->index)[1] << "   " <<
               storage.force(it->index)[2] << std::endl;

  }
  */
    if (false) {
        std::cout << "QM/MM Forces " << std::endl;
        //print out the gradient
        for (std::set<topology::qm_atom_struct>::const_iterator
                     it = topo.qm_zone().begin(), to = topo.qm_zone().end(); it != to; ++it) {
            std::cout << "index :   " << it->index << "  " << storage.force(it->index)[0] << "  "
                      << storage.force(it->index)[1] << "   " <<
                      storage.force(it->index)[2] << std::endl;
        }
    }
    qm_chg.clear();
    return 0;
}

int interaction::MOPAC_Worker::init(topology::Topology& topo,
                                   configuration::Configuration& conf, simulation::Simulation& sim) {
    this->input_file = sim.param().qmmm.mopac.input_file;
    this->output_file = sim.param().qmmm.mopac.output_file;
    this->output_gradient_file = sim.param().qmmm.mopac.output_gradient_file;
    this->header_file=sim.param().qmmm.mopac.input_header;
    if (this->get_qmID() == 2) {
        this->input_file = sim.param().qmmm.mopac.input_file2;
        this->output_file = sim.param().qmmm.mopac.output_file2;
        this->output_gradient_file = sim.param().qmmm.mopac.output_gradient_file2;
        this->header_file=sim.param().qmmm.mopac.input_header2;
    }

    return 0;
}

interaction::MOPAC_Worker::~MOPAC_Worker() {
    this->del_qmID();
}
