/**
 * @file InInput.tcc
 * implements methods of InInput.
 */

/**
 * read the stream into blocks.
 */
inline void io::InInput::read_stream()
{
  std::vector<std::string> buffer;
  
  while(!stream().eof()){

    try{
      io::getblock(stream(), buffer);
    }
    catch(std::runtime_error e){
      break;
    }
    
    m_block[buffer[0]] = buffer;    
    buffer.clear();
    
  }
}

/**
 * read the STEP block.
 */
inline void io::InInput::read_STEP(int &num_steps, double &t0, double &dt)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["STEP"];
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> num_steps >> t0 >> dt;
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in BOND block");
  
}

    
