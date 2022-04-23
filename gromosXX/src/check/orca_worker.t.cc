#include "../stdheader.h"

#include "check.h"

#include "qm_worker.t.h"


int main(int argc, char** argv) {
  // set up simulation files
  #define ROOT_FOLDER "src/check/data/orca/"
  std::string topo_f, conf_f, input_f, qmmm_f, trc_f, tre_f, fin_f;
  GETFILEPATH(topo_f, "md.top", ROOT_FOLDER);
  GETFILEPATH(conf_f, "md.cnf", ROOT_FOLDER);
  GETFILEPATH(input_f, "md.imd", ROOT_FOLDER);
  GETFILEPATH(qmmm_f, "md.qmmm", ROOT_FOLDER);
  GETFILEPATH(trc_f, "md.trc", ROOT_FOLDER);
  GETFILEPATH(tre_f, "md.tre", ROOT_FOLDER);
  GETFILEPATH(fin_f, "md_final.cnf", ROOT_FOLDER);
  
  // initialize test
  test::QM_Worker_Test qm_worker_test(argv[0], "QM/MM Test", topo_f, input_f, conf_f, qmmm_f, trc_f, tre_f, fin_f);
  qm_worker_test.test_friend();
}