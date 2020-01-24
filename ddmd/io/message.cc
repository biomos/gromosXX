/**
 * @file message.cc
 * definition of the message class.
 */

#include <stdheader.h>

#include "message.h"

/**
 * global Message class.
 */
Message messages;

/**
 * add a message.
 */
void Message::add(std::string msg, std::string source,
		  severity_enum severity)
{
  std::string s = source + " : " + msg;
  m_message.insert(message_type(severity, s));
  
  if (severity == critical) display();
}

/**
 * display the messages.
 */
Message::severity_enum Message::display(std::ostream &os)
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
    
  return r;
  
}

void Message::clear()
{
  m_message.clear();
}
