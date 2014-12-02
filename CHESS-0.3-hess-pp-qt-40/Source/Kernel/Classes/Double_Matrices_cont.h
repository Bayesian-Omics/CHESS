/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Double_Matrices_cont.h in the ESS++ program
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


#ifndef DOUBLE_MATRICES_CONT_H
#define DOUBLE_MATRICES_CONT_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <algorithm>
#include <string.h>
#include "Double_Matrices.h"

using namespace std;

class Double_Matrices_cont
{
public:
  Double_Matrices_cont();
  ~Double_Matrices_cont(){};

  int nb_rows;
  int nb_columns;
  double *matrix;
  

  void Alloc_double_matrix_cont(int Rows,
				int Columns);
  void Replicate_double_matrix_cont(Double_Matrices_cont Source);
  void Copy_from_double_matrix(Double_Matrices Source);
  void Free_double_matrix_cont();
  void Read_from_file(char *filename);
  void Display_matrix();
  void Display_matrix_header();
  void Write_to_file(char *filename);

};

#endif /* !defined DOUBLE_MATRICES_H */
