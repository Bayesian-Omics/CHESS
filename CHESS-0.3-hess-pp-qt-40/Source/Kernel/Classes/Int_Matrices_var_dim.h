/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Int_Matrices_var_dim.h in the ESS++ program
 *      Copyright (c) Marc Chadeau (m.chadeau@imperial.ac.uk)
 *                    Leonardo Bottolo (l.bottolo@imperial.ac.uk)
 *                    David Hastie (d.hastie@imperial.ac.uk)
 *      2010
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

#ifndef INT_MATRICES_VAR_DIM_H
#define INT_MATRICES_VAR_DIM_H 1

#include <sstream>
#include <iostream>
#include <fstream>
#include <string.h>
#include <vector>
#include <stdlib.h>

using namespace std;

class Int_Matrices_Var_Dim
{
public:
  Int_Matrices_Var_Dim();
  ~Int_Matrices_Var_Dim(){};

  int nb_rows;
  vector < vector < unsigned int > > matrix;
  
  void Free_matrix();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Write_to_file(char *filename);

};

#endif /* !defined INT_MATRICES_VAR_DIM_H */
