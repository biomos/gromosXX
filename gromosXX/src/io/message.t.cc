/**
 * @file message.t.cc
 * test routines for the message class
 */

#include <iostream>
#include <iomanip>
#include <map>

#include "message.h"

/**
 * testing message
 */
int message_test()
{
  int result = 0;
  
  io::messages.add("a warning.", "message_test",io::message::warning);
  io::messages.add("a notice.", "message_test", io::message::notice);
  
  if (io::messages.display(std::cout) != io::message::warning) ++result;
  
  io::messages.add("a critical.", "message_test", io::message::critical);

  return result;
  
}

int main()
{
  int r1;
  if ((r1 = message_test()))
    std::cout << "message_test failed" << std::endl;
  
  return r1;
}
