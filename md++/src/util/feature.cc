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
 * @file feature.cc
 * features and their cross checks
 */

#include "../stdheader.h"
#include "feature.h"

void util::FeatureChecker::lock_all() {
  m_locked.clear();
  std::vector<Feature>::const_iterator i_it = m_features.begin(), j_it,
                                       to = m_features.end();
  
  // loop over feature pairs
  for(i_it = m_features.begin(); i_it != to; ++i_it) {
    for(j_it = i_it + 1; j_it != to; ++j_it) {
      if (i_it->key() < j_it->key()) { // always use the smaller key as the first key
        m_locked[i_it->key()][j_it->key()] = fc_locked;
      } else {
        m_locked[j_it->key()][i_it->key()] = fc_locked;
      }
    }
  }
}

void util::FeatureChecker::unlock(std::string key1, std::string key2, state s) {
  bool found1 = false, found2 = false;
  std::vector<Feature>::const_iterator it = m_features.begin(),
                                       to = m_features.end();
  
  for(; it != to; ++it) {
    if (it->key() == key1) found1 = true;
    if (it->key() == key2) found2 = true;
  }
  
  // this should never happen! the program will crash at every start.
  if (!found1)
    throw std::runtime_error("Attempt to unlock a nonexisting feature: " + key1);
  if (!found2)
    throw std::runtime_error("Attempt to unlock a nonexisting feature: " + key2);
  
  if (key1 < key2) {
    m_locked[key1][key2] = s;
  } else {
    m_locked[key2][key1] = s;
  }
}

bool util::FeatureChecker::check() {
  bool result = true;
  std::vector<Feature>::const_iterator i_it = m_features.begin(), j_it,
                                       to = m_features.end();
  
  // loop over feature pairs
  for(; i_it != to; ++i_it) {
    if (!i_it->is_active()) continue;
    for(j_it = i_it + 1; j_it != to; ++j_it) {
      if (!j_it->is_active()) continue;
      // both are active -> so do the check
      
      // get the state
      state s = fc_locked;
      if (i_it->key() < j_it->key()) {
        s = m_locked[i_it->key()][j_it->key()];
      } else {
        s = m_locked[j_it->key()][i_it->key()];
      }
      
      switch(s) {
        case fc_locked : // raise an error
        {
          std::ostringstream msg;
          msg << i_it->description() << " does not work with " 
              << j_it->description() << ".";
          io::messages.add(msg.str(), "FeatureChecker", io::message::error);  
          result = false;
          break;
        }
        case fc_unlocked_warning : // raise a warning
        {
          std::ostringstream msg;
          msg << i_it->description() << " used with " 
              << j_it->description() << ".";
          io::messages.add(msg.str(), "FeatureChecker", io::message::error);  
          break;
        }
        case fc_unlocked :
        default: ;
          // fine. do nothing.
      }
    }
  }  
  return result;
}

