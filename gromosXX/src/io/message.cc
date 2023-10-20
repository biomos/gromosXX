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
 * @file message.cc
 * definition of the message class.
 */

#include "../stdheader.h"

#include "message.h"

namespace io
{  
/**
 * global message class.
 */
  message messages;
}

/**
 * add a message.
 */
void io::message::add(std::string msg, std::string source,
		  severity_enum severity)
{
  std::string s = source + " : " + msg;
  m_message.insert(message_type(severity, s)); 
}

void io::message::add(std::string msg, severity_enum severity)
{
  std::string s = msg;
  m_message.insert(message_type(severity, s));
}

/**
 * display the messages.
 */
io::message::severity_enum io::message::display(std::ostream &os)
{
  severity_enum r = none;

  std::multimap<severity_enum, std::string>::const_iterator
    it = m_message.lower_bound(critical),
    it_to = m_message.upper_bound(critical);
  
  if (it != it_to && critical > r) r = critical;
  for( ; it != it_to; ++it){
    os << std::setw(10) << "CRITICAL " << it->second << std::endl;
  }
  
  it = m_message.lower_bound(error);
  it_to = m_message.upper_bound(error);
  
  if(it != it_to && error > r) r = error;
  for( ; it != it_to; ++it)
    os << std::setw(10) << "ERROR " << it->second << std::endl;

  it = m_message.lower_bound(warning);
  it_to = m_message.upper_bound(warning);
  
  if(it != it_to && warning > r) r = warning;
  for( ; it != it_to; ++it)
    os << std::setw(10) << "WARNING " << it->second << std::endl;

  it = m_message.lower_bound(notice);
  it_to = m_message.upper_bound(notice);
  
  if(it != it_to && notice > r) r = notice;
  for( ; it != it_to; ++it)
    os << std::setw(10) << "NOTICE " << it->second << std::endl;
  
  it = m_message.lower_bound(develop); 
  it_to = m_message.upper_bound(develop); 
	   
  if(it != it_to && develop > r) r = develop; 
  for( ; it != it_to; ++it) 
    os << std::setw(10) << "DEVELOP " << it->second << std::endl;     

  return r;
  
}

void io::message::clear()
{
  m_message.clear();
}

bool io::message::contains(severity_enum severity)
{
  return m_message.lower_bound(severity) != m_message.end();
}
