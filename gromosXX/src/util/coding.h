/**
 * @file coding.h
 * de- and encoding of binary data
 */

#ifndef INCLUDED_CODING_H
#define INCLUDED_CODING_H

namespace util {
  /**
   * encode a binary data using RFC1113 base64
   * @return the base64 encoded string with "B64:" at the begin.
   */
  std::string base64_encode(unsigned char const* , unsigned int len);
  /**
   * decode a RFC1113 base64 encoded string
   * @param a string with "B64:" at the begin.
   * @return decoded string or "" if the parameter is not in the right format
   */
  std::string base64_decode(std::string const& s);
  /**
   * puts the text into a frame
   * @param str the string to format
   * @param indent the indent
   * @param boxwidth the width of the box
   * @return the formated text
   */
  std::string frame_text(const std::string & str, unsigned int indent = 8, unsigned int boxwidth = 50);
}

#endif

