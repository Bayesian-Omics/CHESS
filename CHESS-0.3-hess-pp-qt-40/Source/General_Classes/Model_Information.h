/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * CHESS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CHESS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CHESS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MODEL_INFORMATION_H
#define MODEL_INFORMATION_H

#include <stdlib.h>
#include <string>
#include <algorithm>
#include <vector>

using namespace std;

class Model_Information
{
public:
    Model_Information();
    Model_Information(string _Model_Tag,string _Regression_Model);

    string Model_Tag;
    string Regression_Model;

    //for HESS:
    vector < vector < double > > * Omega_PerLine;
    vector < vector < double > > * Rho_PerCol;
    vector < unsigned int > * Active_Yk;



};

#endif // MODEL_INFORMATION_H
