/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from struc.h in the ESS++ program
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

#ifndef STRUC_H
#define STRUC_H

struct data
{
  int nb_rows;
  int nb_columns;
  float **matrix;
};
typedef struct data data;

struct data_double
{
  int nb_rows;
  int nb_columns;
  double **matrix;
};
typedef struct data_double data_double;

struct data_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  float ***matrix;
};
typedef struct data_3D data_3D;

struct data_double_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  double ***matrix;
};
typedef struct data_double_3D data_double_3D;

struct data_integer_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  int ***matrix;
};

typedef struct data_integer_3D data_integer_3D;

struct data_integer
{
  int nb_rows;
  int nb_columns;
  int **matrix;
};
typedef struct data_integer data_integer;

struct vect_double
{
  double nb_rows;
  double *matrix;
};
typedef struct vect_double vect_double;

struct data_long
{
  int nb_rows;
  int nb_columns;
  long **matrix;
};
typedef struct data_long data_long;

struct data_long_double
{
  int nb_rows;
  int nb_columns;
  long double **matrix;
};
typedef struct data_long_double data_long_double;

struct data_long_double_3D
{
  int nb_rows;
  int nb_columns;
  int vect_size;
  long double ***matrix;
};
typedef struct data_long_double_3D data_long_double_3D;

#endif /* STRUC_H */
 

