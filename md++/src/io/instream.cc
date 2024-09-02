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
 * @file instream.cc
 * basic input stream class definition.
 */

#include "../stdheader.h"
#include <vector>
#include <sstream>

#include "blockinput.h"
#include "instream.h"
#include "message.h"

#include "../util/coding.h"
#include <fstream>

#undef MODULE
#undef SUBMODULE
#define MODULE io
#define SUBMODULE topology

void io::GInStream::readTitle() {

  std::vector<std::string> _b;

  io::getblock(*_is, _b);

  if (_b.empty()) {
    io::messages.add("GInStream", "An input file was empty!", io::message::critical);
    return;
  }

  if (_b[0] != "TITLE") {
    if (_b[0] == "") {
      std::ostringstream msg;
      msg << "Empty input file detected.";
      io::messages.add(msg.str(), "GInStream", io::message::critical);
    } else {
      io::messages.add("title block expected: found " + _b[0],
              "GInStream", io::message::error);
    }
    title = "no title";
  } else {
    title = io::concatenate(_b.begin() + 1, _b.end() - 1, title);
  }
}

void io::GInStream::readStream() {
  std::vector<std::string> buffer;
  while (!stream().eof()) {

    if (!io::getblock(stream(), buffer)) {
      if (buffer.size() && buffer[0] != "") {
        std::ostringstream oss;
        oss << "Corrupted block " + buffer[0] << " in file\n" << util::frame_text(title);
        io::messages.add(oss.str(), "GInStream", io::message::critical);
      }
      break;
    }

    trimblock(buffer);

    // empty blocks may cause problems
    if (buffer.size() == 2) {
      std::ostringstream oss;
      oss << "Empty block " + buffer[0] << " in file\n" << util::frame_text(title);
      io::messages.add(oss.str(), "GInStream", io::message::error);
    } else {
      std::string n = buffer[0];
      DEBUG(10, "reading block -" << buffer[0] << "- size " << buffer.size());
      if (m_block.find(buffer[0]) != m_block.end())
        io::messages.add("Found more than one "+buffer[0]+" block.", "GInStream", io::message::warning);
      m_block[buffer[0]] = buffer;
    }
    buffer.clear();
  }
}

