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


#ifndef INCLUDED_BLOCK_H
#define INCLUDED_BLOCK_H

#include <string>
#include <vector>

struct variable;
variable parse_format(std::string);

void parse_data(std::istringstream & is, variable & data,
		std::map<std::string, int> & int_data,
		std::map<std::string, double> & real_data,
		std::map<std::string, std::string> & word_data);


enum data_type { int_type, real_type, word_type, new_line };
struct variable
{
  variable(std::string name, data_type type)
    : name(name), type(type)
  {}
  
  std::string name;
  data_type type;
};

class Block
{
public:

  Block()
    : name(""),
      present(false),
      active(false),
      required(false),
      activate("")
  {
  }


  Block & operator<<(std::vector<std::string> & input)
  {
    if (input.size() < 2 || input[0] != name){
      return *this;
    }
    if (read(input))
      throw gromos::Exception(name, "could not read block data");

    present = true;

    // check activate

    if (activate != ""){

      bool eq = true;
      std::string::size_type it = activate.find(':');
      if (it == std::string::npos){
	it = activate.find('!');
	if (it == std::string::npos)
	  throw gromos::Exception(name, "could not parse activate!");
	eq = false;
      }
    
      std::string w1 = activate.substr(0, it);
      std::string w2 = activate.substr(it+1, std::string::npos);
      
      std::istringstream is(w2);
      if (int_data.count(w1) > 0){
	int d;
	if (!(is >> d))
	  throw gromos::Exception(name, "could not parse activate value");
	if ((eq && d == int_data[w1]) || (!eq && d != int_data[w1]))
	  active = true;
      }
      else if (real_data.count(w1) > 0){
	double d;
	if (!(is >> d))
	  throw gromos::Exception(name, "could not parse activate value");
	if ((eq && d == real_data[w1]) || (!eq && d != real_data[w1]))
	  active = true;
      }
      else if (word_data.count(w1) > 0){
	std::string d;
	if (!(is >> d))
	  throw gromos::Exception(name, "could not parse activate value");
	if ((eq && d == word_data[w1]) || (!eq && d != word_data[w1]))
	  active = true;
      }
      else
	throw gromos::Exception(name, "activation variable not found");
    }
    else
      active = true;
    
    return *this;
  }
  
  int read(std::vector<std::string> & input) 
  {
    std::string s;
    gio::concatenate(input.begin() + 1, input.end() - 1, s);
    std::istringstream is(s);
    
    for(unsigned int i=0; i<data.size(); ++i){

      parse_data(is, data[i], int_data, real_data, word_data);

    }
    
    return 0;
  }

  std::ostream & write(std::ostream & os = std::cout)
  {
    os << name << "\n";

    std::ostringstream os1, os2;
    
    for(unsigned int i=0; i<data.size(); ++i){

      switch(data[i].type){
	case int_type:
	  {
	    os1 << std::setw(8) << data[i].name;
	    if (int_data.count(data[i].name) == 0)
	      os2 << std::setw(8) << data[i].name;
	    else
	      os2 << std::setw(8) << int_data[data[i].name];
	    break;
	  }
	case real_type:
	  {
	    os1 << std::setw(8) << data[i].name;
	    if (real_data.count(data[i].name) == 0)
	      os2 << std::setw(8) << data[i].name;
	    else
	      os2 << std::setw(8) << real_data[data[i].name];
	    break;
	  }
	case word_type:
	  {
	    os1 << std::setw(10) << data[i].name;
	    if (word_data.count(data[i].name) == 0)
	      os2 << std::setw(10) << data[i].name;
	    else
	      os2 << std::setw(10) << word_data[data[i].name];
	    break;
	  }
	case new_line:
	  {
	    os << "#  " << os1.str() << "\n"
	       << "   " << os2.str() << "\n";
	    
	    os1.str("");
	    os2.str("");
	    break;
	  }
      }
    }
    
    os << "END" << std::endl;

    return os;
  }
  

  int check(Message & message, std::vector<Block> & blocks)
  {
    std::cerr << "checking " << name << std::endl;
    
    if (required && !present){
      std::cerr << "\trequired and not present" << std::endl;
      message.add(Message::error, "block " + name + " is required!");
    }
    
    if (present && active){
      for(unsigned int i=0; i<depends.size(); ++i){
	for(unsigned int j=0; j<blocks.size(); ++j){
	  if (depends[i] == blocks[j].name && !blocks[j].active)
	    message.add(Message::error, "block " + name + " requires " + blocks[j].name);
	}
      }
      for(unsigned int i=0; i<conflicts.size(); ++i){
	for(unsigned int j=0; j<blocks.size(); ++j){
	  if (conflicts[i] == blocks[j].name && blocks[j].active)
	    message.add(Message::error, "block " + name + " conflicts with " + blocks[j].name);
	}
      }
    }

    return 0;
  }

  int update(std::map<std::string, int> & idat,
	     std::map<std::string, double> & ddat,
	     std::map<std::string, std::string> & wdat)
  {
    {
      std::map<std::string, int>::iterator 
	it = int_data.begin(),
	to = int_data.end();
      for( ; it!=to; ++it){
	std::cerr << "\ttrying to update " << it->first << " (" << it->second << ")" << std::endl;
	if (idat.count(name + "." + it->first) > 0) it->second = idat[name + "." + it->first];
      }
    }
    {
      std::map<std::string, double>::iterator 
	it = real_data.begin(),
	to = real_data.end();
      for( ; it!=to; ++it){
	std::cerr << "\ttrying to update " << it->first << " (" << it->second << ")" << std::endl;	
	if (ddat.count(name + "." + it->first) > 0) it->second = ddat[name + "." + it->first];
      }
    }
    {
      std::map<std::string, std::string>::iterator 
	it = word_data.begin(),
	to = word_data.end();
      for( ; it!=to; ++it){
	std::cerr << "\ttrying to update " << it->first << " (" << it->second << ")" << std::endl;	
	if (wdat.count(name + "." + it->first) > 0) it->second = wdat[name + "." + it->first];
      }
    }

    return 0;
  }

  void extract_var(int runid,
		   std::map<std::string, double> & parameter)
  {
    std::ostringstream os;
    if (runid >= 0)
      os << "ID" << runid << "." << name << ".";
    else
      os << name << ".";
    
    {
      std::map<std::string, int>::iterator 
	it = int_data.begin(),
	to = int_data.end();
      for( ; it!=to; ++it){
	std::cerr << "PARAM(int): " << os.str() << it->first << std::endl;
	parameter[os.str() + it->first] = double(it->second);
      }
    }
    {
      std::map<std::string, double>::iterator 
	it = real_data.begin(),
	to = real_data.end();
      for( ; it!=to; ++it){
	std::cerr << "PARAM(real): " << os.str() << it->first << std::endl;
	parameter[os.str() + it->first] = double(it->second);
      }
    }
    /*
    {
      std::map<std::string, std::string>::iterator 
	it = word_data.begin(),
	to = word_data.end();
      for( ; it!=to; ++it)
	parameter[os.str() + it->first] = double(it->second);
    }
    */
  }

  void init(std::vector<std::string> & format)
  {
    if (format.size() <= 2)
      throw gromos::Exception("block", "could not parse this! (block size < 3)");
    
    name = format[0];

    if (format.back() != "END")
      throw gromos::Exception(name, "no END at the end of block");
    
    for(unsigned int i=1; i<format.size()-1; ++i){

      std::cerr << "\tparsing line " << format[i] << std::endl;
      
      std::istringstream is(format[i]);

      bool data_line = false;

      while(!is.eof()){
	
	std::string w;
	is >> w;
	
	if (w == "DEPENDS" || w == "depends"){
	  while(!is.eof()){
	    if(is >> w){
	      std::cerr << "\t\tdepends on " << w << std::endl;
	      depends.push_back(w);
	    }
	  }
	}
	else if (w == "CONFLICTS" || w == "conflicts"){
	  while(!is.eof()){
	    if (is >> w){
	      std::cerr << "\t\tconflicts with " << w << std::endl;
	      conflicts.push_back(w);
	    }
	  }
	}
	else if (w == "REQUIRED" || w == "required"){
	  std::cerr << "\t\trequired = true" << std::endl;
	  required = true;
	}
	else if (w == "ACTIVATE" || w == "active"){

	  if (is >> w){
	    activate = w;	    
	  }
	  while(!is.eof())
	    is >> w;

	}
	else{
	  std::cerr << "\t\tparsing variable " << w << std::endl;
	  data_line = true;
	  data.push_back(parse_format(w));

	}
	
      }
      if (data_line)
	data.push_back(variable("cr/lf", new_line));

    }
  }
  
  
  std::string name;
  
  bool present;
  bool active;
  bool required;
  std::string activate;
  
  std::vector<std::string> depends;
  std::vector<std::string> conflicts;
  
  std::vector<variable> data;

  std::map<std::string, int>    int_data;
  std::map<std::string, double> real_data;
  std::map<std::string, std::string> word_data;

};


  
#endif
