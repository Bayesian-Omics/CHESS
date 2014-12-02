/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from AdMH.h in the ESS++ program
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

#ifndef AdMH_H
#define AdMH_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <math.h>
#include "../Routines/matrix_handling.h"
#include <gsl/gsl_statistics_double.h>
#include <vector>
#include <algorithm>

#include "../Classes/Move_monitor.h"

using namespace std;

class AdMH
{
public:
  AdMH();
  ~AdMH(){};
  unsigned int  n_batch;
  double optimal;
  double ls;
  double delta_n;

  unsigned int tilda_accept;
  unsigned int tilda_accept_ins;
  unsigned int tilda_n_sweep;
  unsigned int tilda_n_sweep_ins;
  vector < double > M;
  vector < double > Ls;
  
  void set_AdMH(int g_sample,
		  unsigned int n_batch_from_read,
          double AdMH_optimal_from_read,
          double AdMH_ls_from_read,
		  unsigned int pX,
		  unsigned int burn_in,
          double M_min_input,
          double M_max_input);
  

  void display_AdMH();

  void Adapt_AdMH_ESS(unsigned int sweep, Move_monitor *My_Move_monitor);

  double Adapt_AdMH(unsigned int sweep);


};

#endif /* !defined AdMH_H */
