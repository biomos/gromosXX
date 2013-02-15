/**
 * @file instream.h
 * basic input stream class.
 */

#ifndef INCLUDED_GINSTREAM_H
#define INCLUDED_GINSTREAM_H

namespace io {

  class GInStream {

  public:

    /*
     * default constructor
     */
    GInStream() : quiet(false), _auto_delete(false) {} 

    /*
     * Constructor
     */
    GInStream(std::istream& is) : quiet(false), _auto_delete(false) { stream(is); } 

    /**
     * Destructor.
     */
    ~GInStream() { if (_auto_delete) delete _is; }
    
    /*
     * Accessors to the input stream
     */
    std::istream& stream() { return *_is; }
    void stream(std::istream& is) { _is = &is; readTitle(); }

    /**
     * Read a title block from the input stream,
     * concatenate and store it in title.
     * @section title TITLE block
 @verbatim
 TITLE
  your title
 END
 @endverbatim
     */
    void readTitle();
    std::string title;
    /**
     * read the entire stream and store the blocks in the map.
     */
    void readStream();
    /**
     * auto delete accessor.
     */
    void auto_delete(bool b) { _auto_delete = b; }

    /**
     * output during reading?
     */
    bool quiet;

  protected:
    std::istringstream _lineStream;
    /**
     * stores the blocks if read_stream is called.
     */
    std::map<std::string, std::vector<std::string> > m_block;

    /**
     * the stream.
     */
    std::istream* _is;
    /**
     * delete the stream in the end?
     */
    bool _auto_delete;
    
  };

} // io

#endif
