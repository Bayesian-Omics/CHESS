/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from CM.cc in the ESS++ program
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

#include "CM.h"
#define DEBUG 0

using namespace std;

CM::CM()
{
  n_max_breakpoint=0;
  n_possible_CM_moves=0;
}

void CM::set_CM(unsigned int k_max_from_read,
		vector < unsigned int  > &list_CM_moves_enabled_from_read)
{
  unsigned int Is_k_point_CM_in=0;
  n_max_breakpoint=k_max_from_read;
  //Getting the number of possible CM moves
  for(unsigned int col=0;col<list_CM_moves_enabled_from_read.size();col++){
    if(list_CM_moves_enabled_from_read[col]==1){
      Is_k_point_CM_in=1;
      if(n_max_breakpoint==0){
	cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	     << "           WARNING" << endl
	     << "  k-point CM enabled, but" << endl
	     << "    n_max_breakpoint is set to 0" << endl
	     << " ****** Run stopped ****** " << endl
	     << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	exit(1);
      }
    }
  }
  //Calculating the pbty of each move: col in 0:k-1: k-point CM; k haplotype
  n_possible_CM_moves=Is_k_point_CM_in*(n_max_breakpoint)+(list_CM_moves_enabled_from_read.size()-1);

  unit_move_pbty_cum.resize(n_possible_CM_moves);
  if(unit_move_pbty_cum.size()>1){
    for(unsigned int col=0;col<unit_move_pbty_cum.size();col++){
      unit_move_pbty_cum[col]=(double)(col+1)/(double)(n_possible_CM_moves);
    }
  }
 
}

void CM::display_CM()
{
  cout << endl << "**********************************************************" << endl
       << "********************** CM parameters *********************" << endl 
       << "\tn_max_breakpoint = " << n_max_breakpoint << endl
       << "\tnb_possible_CM_move = " << n_possible_CM_moves << endl
       << "\tInit CM pbty cum:" << endl
       << "\t";
  for(unsigned int col=0;col< unit_move_pbty_cum.size();col++){
    cout << unit_move_pbty_cum[col] << " ";
  }
  cout << endl
       << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
  
}
