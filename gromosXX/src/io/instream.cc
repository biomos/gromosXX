/**
 * @file instream.cc
 * basic input stream class definition.
 */

#include <stdheader.h>

#include "blockinput.h"
#include "instream.h"

void io::GInStream::readTitle() {

  std::vector<std::string> _b;

  io::getblock(*_is, _b);
  if (_b[0] != "TITLE")
    io::messages.add("title block expected: found " + _b[0],
		     "instream",
		     io::message::error);
  title = io::concatenate(_b.begin() + 1, _b.end() - 1, title);
}

void io::GInStream::readStream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    if (!io::getblock(stream(), buffer)){
      if (buffer.size() && buffer[0] != ""){
	std::cerr << "invalid block " + buffer[0] << " in input file?"
		  << std::endl;
      }
      break;
    }

    trimblock(buffer);
    
    m_block[buffer[0]] = buffer;
    
    buffer.clear();
    
  }
  
}

