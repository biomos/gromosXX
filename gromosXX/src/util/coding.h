/**
 * @file coding.h
 * de- and encoding of binary data
 */

#ifndef INCLUDED_CODING_H
#define INCLUDED_CODING_H

namespace util {
  /**
   * encode a binary data using RFC1113 base64
   */
  std::string base64_encode(unsigned char const* , unsigned int len);
  /**
   * decode a RFC1113 base64 encoded string
   */
  std::string base64_decode(std::string const& s);
}

#endif

