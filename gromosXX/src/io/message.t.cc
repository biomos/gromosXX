/**
 * @file message.t.cc
 * test routines for the message class
 */

#include <iostream>
#include <iomanip>
#include <map>

#include "message.h"

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>

using namespace boost::unit_test_framework;

/**
 * testing message
 */
void message_test()
{
  io::messages.add("a warning.", "message_test",io::message::warning);
  io::messages.add("a notice.", "message_test", io::message::notice);
  
  BOOST_CHECK_EQUAL(io::messages.display(std::cout), io::message::warning);
  
  io::messages.add("a critical.", "message_test", io::message::critical);
  
}

test_suite*
init_unit_test_suite(int argc, char*argv[])
{
  test_suite* test=BOOST_TEST_SUITE("message test");
  test->add(BOOST_TEST_CASE( &message_test));
  
  return test;
}
