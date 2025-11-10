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


#ifndef INCLUDED_MESSAGE_H
#define INCLUDED_MESSAGE_H

#include <iostream>
#include <string>
#include <vector>

class Message
{
public:
  enum level { warning = 0, error=1 };
  
  void add(level l, std::string message)
  {
    if (l) m_errors.push_back(message);
    else m_warnings.push_back(message);
  }
  
  void display(std::ostream & os = std::cout)
  {
    for(unsigned int i=0; i<m_errors.size(); ++i)
      os << "ERROR:   " << m_errors[i] << std::endl;
    for(unsigned int i=0; i<m_warnings.size(); ++i)
      os << "WARNING: " << m_warnings[i] << std::endl;
  }
  
private:
  std::vector<std::string> m_errors;
  std::vector<std::string> m_warnings;
    
};

#endif
