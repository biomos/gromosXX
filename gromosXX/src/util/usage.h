/**
 * @file usage.h
 */

#ifndef INCLUDED_USAGE_H
#define INCLUDED_USAGE_H

namespace util
{
  /**
   * @class Known
   * Known
   */
  class Known : public std::set<std::string>
  {
  public:
    Known & operator << (std::string s)
    {
      insert(s);
      return *this;
    }
  };

  void get_usage(Known const &knowns, std::string &s, std::string name);

  void print_title(bool mpi, std::ostream & os = std::cout,
          bool repex = false);
  
}

#endif

  
