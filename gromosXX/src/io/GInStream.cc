/**
 * @file GInStream.cc
 * basic input stream class definition.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <map>

#include "blockinput.h"
#include "GInStream.h"

void io::GInStream::readTitle() {

  std::vector<std::string> _b;

  io::getblock(*_is, _b);
  if (_b[0] != "TITLE")
    throw std::runtime_error("TITLE block expected. Found: " + _b[0]);
  title = io::concatenate(_b.begin() + 1, _b.end() - 1, title);
}

void io::GInStream::readStream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    try{
      io::getblock(stream(), buffer);
    }
    catch(std::runtime_error e){
      break;
    }
    
    trimblock(buffer);
    // std::cerr << buffer[0] << std::endl;
    
    m_block[buffer[0]] = buffer;
    
    buffer.clear();
    
  }
}

