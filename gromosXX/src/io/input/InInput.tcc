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
 * read the SYSTEM block.
 */
inline void io::InInput::read_SYSTEM(int &nsm)
{
  std::vector<std::string> buffer;
  buffer = m_block["SYSTEM"];
  _lineStream.clear();
  _lineStream.str(buffer[1]);
  
  int npm;
  _lineStream >> npm >> nsm;
 
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in SYSTEM block");
  
  if (npm != 1)
    io::messages.add("SYSTEM: only NPM=1 allowed",
		     "io::InInput::read_SYSTEM",
		     io::message::error);
  
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
    throw std::runtime_error("bad line in STEP block");
  
}

/**
 * the SHAKE block.
 */
inline void io::InInput::read_SHAKE(int &ntc, double &tolerance)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  buffer = m_block["SHAKE"];
  
  it = buffer.begin() + 1;
  _lineStream.clear();
  _lineStream.str(*it);
  
  _lineStream >> ntc >> tolerance;
  
  if (_lineStream.fail() || ! _lineStream.eof())
    throw std::runtime_error("bad line in SHAKE block");

}
