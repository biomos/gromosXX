/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

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
   * @param s a string with "B64:" at the begin.
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

