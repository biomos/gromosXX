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

#include "CommandLine.h"

#include <cctype>
#include <iostream>
#include <string>
#include <ostream>

#include "../../config.h"

#include <algorithm>

#if defined(HAVE_READLINE) && defined(HAVE_LIBREADLINE)
#include <cstdio>
#include <cstdlib>

#include <readline/readline.h>
#include <readline/history.h>
#endif

std::string utils::CommandLine::getLine(const std::string& prompt, std::ostream & out) {
  std::string line;
#if defined(HAVE_READLINE) && defined(HAVE_LIBREADLINE)
  char * buf;
  //rl_bind_key('\t',rl_abort);//disable auto-complete
  buf = readline(prompt.c_str());
  if (buf == NULL) return line;
  line = buf;
  if (buf[0]!=0) add_history(buf);
  free(buf);
#else
  if (!prompt.empty()) out << prompt;
  std::getline(std::cin, line);
#endif
  return line;
}

char to_lower(char c) { return std::tolower(c); }

bool utils::CommandLine::getYesNo(const std::string& prompt, std::ostream& out, const std::string& error) {
  while(true) {
    std::string a = getLine(prompt, out);
    std::transform(a.begin(), a.end(), a.begin(), to_lower);
    if (a[0] == 'y') return true;
    if (a[0] == 'n') return false;
    out << error;
  }
}
