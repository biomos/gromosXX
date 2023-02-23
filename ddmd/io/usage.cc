/**
 * @file usage.cc
 * get usage string
 */

#include <stdheader.h>
#include "usage.h"

void get_usage(Known const &knowns, std::string &usage, std::string name)
{
  usage = "#\n# " + name + "\n\n";

  if (knowns.count("topo")){
    usage += "\t# topology data\n";
    usage += "\t@topo        filename\n\n";
  }
  
  if (knowns.count("pttopo")){
    usage += "\t# perturbation topology data\n";
    usage += "\t# @pttopo    filename\n\n";
  }
  
  if (knowns.count("conf")){
    usage += "\t# coordinates\n";
    usage += "\t@conf        filename\n\n";
  }
  
  if (knowns.count("input")){
    usage += "\t# input parameter\n";
    usage += "\t@input       filename\n\n";
  }

  if (knowns.count("fin")){
    usage += "\t# output finale coordinates\n";
    usage += "\t@fin         filename\n\n";
  }
  
  if (knowns.count("trj")){
    usage += "\t# output coordinates trajectory\n";
    usage += "\t@trj         filename\n\n";
  }
  
  if (knowns.count("trv")){
    usage += "\t# output velocity trajectory\n";
    usage += "\t# @trv       filename\n\n";
  }
  
  if (knowns.count("trf")){
    usage += "\t# output force trajectory\n";
    usage += "\t# @trf       filename\n\n";
  }
  
  if (knowns.count("tre")){
    usage += "\t# output energy trajectory\n";
    usage += "\t# @tre       filename\n\n";
  }
  
  if (knowns.count("bae")){
    usage += "\t# output block averaged energy trajectory\n";
    usage += "\t# @bae       filename\n\n";
  }
  
  if (knowns.count("trg")){
    usage += "\t# output free energy trajectory\n";
    usage += "\t# @trg       filename\n\n";
  }
  
  if (knowns.count("bag")){
    usage += "\t# output block averaged free energy trajectory\n";
    usage += "\t# @bag       filename\n\n";    
  }
  
  if (knowns.count("posres")){
    usage += "\t# position restraints specification\n";
    usage += "\t# @posres    filename\n\n";
  }
  
  if (knowns.count("distrest")){
    usage += "\t# distance restraints specification\n";
    usage += "\t# @distrest  filename\n\n";
  }
    //ORIOL_GAMD
  if (knowns.count("gamd")){
    usage += "\t# GAMD restraints specification\n";
    usage += "\t# @gamd      filename\n\n";
  }
  
  if (knowns.count("jval")){
    usage += "\t# J-value restraints specification\n";
    usage += "\t# @jval      filename\n\n";
  }
  
  if (knowns.count("rep")){
    usage += "\t# deprecated: replica exchange information\n";
    usage += "\t# @rep       filename\n\n";
  }

  if (knowns.count("master")){
    usage += "\t# replica exchange: master process\n";
    usage += "\t# @master    name\n\n";
  }

  if (knowns.count("slave")){
    usage += "\t# replica exchange: slave process\n";
    usage += "\t# @slave     name\n\n";
  }
  
  if (knowns.count("print")){
    usage += "\t# print additional information\n";
    usage += "\t# @print     <pairlist/force>\n\n";
  }
  
  if (knowns.count("trp")){
    usage += "\t# write additional information to file\n";
    usage += "\t# @trp       filename\n\n";
  }
  
  if (knowns.count("anatrj")){
    usage += "\t# re-analyze trajectory\n";
    usage += "\t# @anatrj    filename\n\n";
  }
  
  if (knowns.count("verb")){
    usage += "\t# control verbosity (in debug builds)\n";
    usage += "\t# @verb      <[module:][submodule:]level>\n\n";
  }
  
  if (knowns.count("version")){
    usage += "\t# print version information\n";
    usage += "\t# @version\n\n";
  }
  
  if (knowns.count("print")){
    usage += "\t# print pairlist / force\n";
    usage += "\t# @print pairlist/force\n";
  }

  usage += "#\n\n";

}


