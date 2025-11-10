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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   StringOps.cc
 * Author: bschroed
 * 
 * Created on December 18, 2017, 4:04 PM
 */

#include "StringOps.h"

#include <sstream>
#include <string>

namespace utils{
    string StringOps::replaceall(string input, string replacepattern, string insertpattern){

        int replacepattternlength= replacepattern.length();
        int startpos = 0;
        int found_pos = input.find( replacepattern, startpos); //initial search
        ostringstream finalstr;

        //get all start positions
        while(found_pos != string::npos){
            finalstr << input.substr(startpos, found_pos-startpos) << insertpattern; //found a position, build up new string in stream final
            startpos = found_pos + replacepattternlength; // set new star position for search
            found_pos = input.find( replacepattern, startpos); // get new position
        }

        finalstr << input.substr(startpos, input.length()-startpos);    //append suffix
        return finalstr.str();
    }
}
