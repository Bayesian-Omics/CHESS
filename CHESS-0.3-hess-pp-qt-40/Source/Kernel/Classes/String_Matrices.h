/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from String_Matrices.h in the ESS++ program
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

#ifndef STRING_MATRICES_H
#define STRING_MATRICES_H 1


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <string.h>

using namespace std;

class String_Matrices
{
public:
  String_Matrices();
  ~String_Matrices(){};

  int nb_rows;
  int nb_columns;
  string **matrix;
  

  void Alloc_string_matrix(int Rows,
			   int Columns);
  void Free_string_matrix();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Display_matrix_header();
  void Write_to_file(char *filename);

};

#endif /* !defined STRING_MATRICES_H */
