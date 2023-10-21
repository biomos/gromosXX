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
 * @file message.h
 * container for error messages and warnings.
 */

#ifndef INCLUDED_MESSAGE_H
#define INCLUDED_MESSAGE_H

namespace io
{
  /**
   * @class message
   * collects messages.
   */
  class message
  {
  public:
    /**
     * @enum severity_enum
     * severity level.
     */
    enum severity_enum{ 
      /**
       * no errors and warnings.
       */
      none,
      /**
       * informative.
       */
      notice, 
      /**
       * considered possibly wrong but
       * simulation can go on.
       */
      warning, 
      /**
       * informative message for code under development -> program exits 
       */ 
      develop,  
      /** 
       * an error. stop at the next
       * critical point.
       */
      error, 
      /**
       * a critical error. stop immediately.
       */
      critical 
    };
    /**
     * add a message.
     */
    void add(std::string msg, std::string source, severity_enum severity);
    void add(std::string msg, severity_enum severity); 
    /**
     * display the message, return
     * the highest severity of the messages
     * in the container.
     */
    severity_enum display(std::ostream &os = std::cerr);
    /**
     * clear the messages.
     */
    void clear();
    /**
     * Checks whether the messages contain a message with a certain severity.
     */
    bool contains(severity_enum severity);
  private:
    typedef std::multimap<severity_enum, std::string>
      ::value_type message_type;

    std::multimap<severity_enum, std::string> m_message;
    
  };
  
  extern message messages;

} // io

#endif

  
  
    
