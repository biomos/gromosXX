/**
 * @file GInStream.h
 * basic input stream class.
 */

#ifndef INCLUDED_GINSTREAM_H
#define INCLUDED_GINSTREAM_H

namespace io {

  class GInStream {

  public:

    /*
     * Default constructor opens the GInStream with stdin.
     */
    GInStream(std::istream& is = std::cin) { stream(is); } 


    /*
     * Accessors to the input stream
     */
    std::istream& stream() { return *_is; }
    void stream(std::istream& is) { _is = &is; readTitle(); }

    /*
     * Read a title block from the input stream,
     * concatenate and store it in title.
     */
    void readTitle();
    std::string title;

  protected:
    std::istringstream _lineStream;

  private:   
    std::istream* _is;
  };

} // io

#endif
