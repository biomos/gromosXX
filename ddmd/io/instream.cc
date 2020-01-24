/**
 * @file instream.cc
 * basic input stream class definition.
 */

#include "../stdheader.h"

#include "blockinput.h"
#include "instream.h"

void GInStream::readTitle() {

  std::vector<std::string> _b;

  getblock(*_is, _b);
  if (_b[0] != "TITLE")
    messages.add("title block expected: found " + _b[0],
		     "instream",
		     Message::error);
  title = concatenate(_b.begin() + 1, _b.end() - 1, title);
}

void GInStream::readStream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    if (!getblock(stream(), buffer)){
      if (buffer.size() && buffer[0] != ""){
	std::cerr << "invalid block " + buffer[0] << " in input file?"
		  << std::endl;
      }
      break;
    }

    trimblock(buffer);
    
    std::string n = buffer[0];
    m_block[buffer[0]] = buffer;
    
    buffer.clear();

  }
  
}

