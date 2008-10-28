/**
 * @file inframe.h
 * frame based input / output for trajectories
 */

#ifndef INCLUDED_INFRAME_H
#define INCLUDED_INFRAME_H

#include "../instream.h"

namespace io {

  class In_Frame : public GInStream {

  public:

    /*
     * default constructor
     */
    In_Frame() : GInStream(){} 

    /*
     * Constructor
     */
    In_Frame(std::istream& is) : GInStream(is) {
      io::getline(*_is, _next_block);
      _first_block = _next_block;
    } 
    
    /*
     * Accessors to the input stream
     */
    void stream(std::istream& is) { 
      GInStream::stream(is);io::getline(*_is, _next_block);
      _first_block = _next_block;  
    }

    /**
     * read frame.
     */
    int read_frame();
    

  protected:
    /**
     * next block name in file.
     */
    std::string _next_block;
    /**
     * first block name in file (after title block)
     */
    std::string _first_block;
    
  };

} // io

#endif
