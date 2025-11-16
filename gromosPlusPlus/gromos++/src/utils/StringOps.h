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
 * File:   StringOps.h
 * Author: bschroed
 *
 * Created on December 18, 2017, 3:46 PM
 */



#ifndef STRINGOPS_H
#define STRINGOPS_H

#include <string>
#include <list>
#include <iostream>
#include <sstream>

using namespace std;

namespace utils {
/**
     * @class StrigOps
     * 
     * contains some internal string operations
     * @param
     * @return 
     */
    class StringOps {
    public:
        //StringOps();
        //StringOps(const StringOps& orig);
        //virtual ~StringOps();

        //function replaces all replacepattern occurences with insertpatterns
        static string replaceall(std::string input, std::string replacepattern, std::string insertpattern);

    private:

    };
}
#endif /* STRINGOPS_H */
