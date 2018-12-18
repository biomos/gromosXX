/**
 * @file mopac_worker.h
 * The worker class for the MOPAC software
 */


#ifndef MOPAC_WORKER_H
#define	MOPAC_WORKER_H

#include "qm_worker.h"


namespace interaction {
    class QM_Worker;
    /**
     * @class MNDO_Worker
     * a worker class which calls the MNDO software
     */
    class MOPAC_Worker : public QM_Worker {
    public:
        /**
         * Constructor
         */
        MOPAC_Worker() : QM_Worker("MOPAC Worker") {
          this->get_new_qmID();
        }
        /**
         * Destructor
         */
        virtual ~MOPAC_Worker();
        /**
         * initialise the QM worker
         * @return 0 if successful, non-zero on failure
         */
        virtual int init(topology::Topology & topo,
                         configuration::Configuration & conf,
                         simulation::Simulation & sim);
        /**
         * run a QM job in MOPAC
         * @param qm_pos a vector containing the QM atom positions
         * @param mm_atoms the MM atoms to include
         * @param storage the energies, forces, charges obtained
         * @return 0 if successful, non-zero if not.
         */
        virtual int run_QM(topology::Topology & topo,
                           configuration::Configuration & conf,
                           simulation::Simulation & sim,
                           const math::VArray & qm_pos,
                           const std::vector<MM_Atom> & mm_atoms,
                           interaction::QM_Storage & storage);
    private:
        double pointchg_pot(const int j,
                        const std::vector<MM_Atom>&mm_atoms,
                        const topology::Topology & topo,
                        const configuration::Configuration &conf);
        math::Vec pointchg_force(const int j,
                                 std::vector< double> qm_chg,
                            const topology::Topology & topo,
                            const configuration::Configuration &conf);
        /**
         * file name for MOPAC input file
         */
        std::string input_file;
        /**
         * file name for MOPAC output file
         */
        std::string output_file;
        /**
         * file name for MOPAC gradient output file
         */
        std::string output_gradient_file;
        /**
         * file name for MOPAC header file
         */
        std::string header_file;
        /**
         * file name for MOPAC header file
         */
        std::string molin_file;
        /**
         * file name for MOPAC header file
         */
        std::string molin_file2;
    };
}

#endif	/* MOPAC_WORKER_H */
