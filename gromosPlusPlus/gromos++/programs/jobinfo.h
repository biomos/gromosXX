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


#ifndef INCLUDED_JOBINFO_H
#define INCLUDED_JOBINFO_H

#include "block.h"

class Jobinfo
{
public:
  int read(std::vector<std::string> & input)
  {
    if (input.size() < 3)
      throw gromos::Exception("jobinfo", "not enough jobinfo data");
    
    // read data format
    std::istringstream is(input[1]);
    
    while (!is.eof()){
      std::string w;
      is >> w;
      
      if (w == "JOBID" || w == "jobid" || w == "job_id"){
	data.push_back(variable("JOBID", int_type));
      }
      else if (w == "dir" || w == "DIR"){
	data.push_back(variable("DIR", word_type));
      }
      else if (w == "runafter" || w == "RUNAFTER"){
	data.push_back(variable("RUNAFTER", int_type));
      }
      else{
	data.push_back(parse_format(w));
	std::cerr << "\tjoblist variable " << data.back().name << std::endl;
      }
      
    }
    std::cerr << "\treading " << input.size() - 3 << " jobs" << std::endl;
    
    // int_data.resize(input.size() - 3);
    // real_data.resize(input.size() - 3);
    // word_data.resize(input.size() - 3);
    
    for(unsigned int i=2; i<input.size()-1; ++i){
      std::istringstream is(input[i]);

      // first HAS to be JOBID !
      int jobid;
      if (!(is >> jobid))
	throw gromos::Exception("Jobinfo", "could not parse JOBID");

      int_data[jobid]["JOBID"] = jobid;
      
      for(unsigned int d = 1; d<data.size(); ++d){
	parse_data(is, data[d], int_data[jobid], real_data[jobid], word_data[jobid]);
      }
    }

    return 0;
  }

  int size() const
  {
    return int_data.size();
  }
  
  void back_substitute(std::map<std::string, double> & parameter,
		       std::map<std::string, data_type> & type)
  {
    std::cerr << "BACKSUBSTITUTE" << std::endl;

    std::map<std::string, data_type>::const_iterator
      it = type.begin(),
      to = type.end();
    for( ; it!=to; ++it){
      
      std::cerr << "\tbacksubstitute: " << it->first << std::endl;
      std::string::size_type dot = it->first.find('.');
      std::string::size_type dot2 = it->first.find('.', dot + 1);
      
      if (dot == std::string::npos || dot2 == std::string::npos ||
	  dot2 <= dot)
	throw gromos::Exception("Jobinfo", "Could not backsubstitute " + it->first);
      
      int j;
      {
	std::istringstream is(it->first.substr(2, dot));
	if (!(is >> j))
	  throw gromos::Exception("Jobinfo", "Could not parse ID ('" + it->first.substr(2, dot) + "')");
      }
      std::string name;
      {
	std::istringstream is(it->first.substr(dot + 1, dot2 - dot - 1));
	if (!(is >> name))
	  throw gromos::Exception("Jobinfo", "Could not parse block name ('"
				  + it->first.substr(dot+1, dot2-dot-1) + "')");
      }
      std::string varname;
      {
	std::istringstream is(it->first.substr(dot2+1, std::string::npos));
	if (!(is >> varname))
	  throw gromos::Exception("Jobinfo", "Could not parse varialbe name ('"
				  + it->first.substr(dot2+1, std::string::npos) + "')");
      }
      
      std::string common_name = name + "." + varname;
      std::cerr << "\t\tsubstituting: " << j << " " << name << " " << varname << std::endl;
      std::cerr << "\t\twith common name : " << common_name << std::endl;
      std::cerr << "\t\tand value : " << parameter[it->first] << std::endl;

      switch(it->second){
	case int_type:
	  int_data[j][common_name] = int(parameter[it->first]);
	  break;
	case real_type:
	  real_data[j][common_name] = parameter[it->first];
	  break;
	case word_type:
	  throw gromos::Exception("Jobinfo", "Can not substitute word...");
	  break;
	case new_line:
	  throw gromos::Exception("Jobinfo", "new_line does not belong here...");
      }
      

    }
    
  }
  
  std::vector<variable> data;

  std::map<int, std::map<std::string, int> > int_data;
  std::map<int, std::map<std::string, double> > real_data;
  std::map<int, std::map<std::string, std::string> > word_data;
  
};

#endif
