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

  
  
    
