/**
 * @file inframe.cc
 * basic input stream class definition.
 */

#include <stdheader.h>

#include "../blockinput.h"
#include "../instream.h"
#include "inframe.h"
#include <util/error.h>

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE configuration

int io::In_Frame::read_frame()
{
  std::vector<std::string> buffer;

  // clear the map
  m_block.clear();
  
  DEBUG(7, "reading frame until block " << _first_block);
  
  while(!GInStream::stream().eof()){

    // set the blockname
    DEBUG(8, "reading block " << _next_block);
    
    try{
      io::getblock(GInStream::stream(), buffer);
    }
    catch(std::runtime_error e){
      if (buffer.size() && buffer[0] != ""){
	std::string s = "error while reading frame : " + buffer[0];
	io::messages.add(s, "inframe", io::message::error);
	return E_INPUT_ERROR;
      }
      
      break;
    }

    trimblock(buffer);
    
    m_block[_next_block] = buffer;
    
    buffer.clear();

    io::getline(*_is, _next_block);
    // one frame read ?
    if (_next_block == _first_block)
      break;
    
  }
  return 0;
}

