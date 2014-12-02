/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from DR.cc in the ESS++ program
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

#include "DR.h"
#define DEBUG 0

using namespace std;

DR::DR()
{
  nb_calls_adj=0;
  nb_calls=0;

}

void DR::set_DR(unsigned int n_chains)
{
  nb_calls_adj=0;
  nb_calls=0;
  mat_moves_accepted.resize(n_chains);
  for(unsigned int row=0;row<mat_moves_accepted.size();row++){
    mat_moves_accepted[row].resize(n_chains);
  }
  mat_moves_proposed.resize(n_chains);
  for(unsigned int row=0;row<mat_moves_proposed.size();row++){
    mat_moves_proposed[row].resize(n_chains);
  }
  
}

void DR::display_DR()
{
  cout << endl << "**********************************************************" << endl
       << "********************** DR parameters *********************" << endl 
       << "\tnb_calls_adj = " << nb_calls_adj << endl
       << "\tnb_calls = " << nb_calls << endl;
  cout << endl
       << "\tMat_moves_accepted" 
       << " n_rows " << mat_moves_accepted.size() 
       << " -- ncol= " << mat_moves_accepted[0].size() << endl;
  for(unsigned int row=0;row<mat_moves_accepted.size();row++){
    cout <<  "\t";
    for(unsigned int col=0;col<mat_moves_accepted[row].size();col++){
      cout << mat_moves_accepted[row][col] << "\t";
    }
    cout << endl;
  }
 cout << endl
       << "\tMat_moves_proposed" 
       << " n_rows " << mat_moves_proposed.size() 
       << " -- ncol= " << mat_moves_proposed[0].size() << endl;
  for(unsigned int row=0;row<mat_moves_proposed.size();row++){
    cout <<  "\t";
    for(unsigned int col=0;col<mat_moves_proposed[row].size();col++){
      cout << mat_moves_proposed[row][col] << "\t";
    }
    cout << endl;
  }
  cout << endl;
  cout << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
  
}
