/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Int_Matrices_var_dim.cc in the ESS++ program
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

#include "Int_Matrices_var_dim.h"
#define DEBUG 0

using namespace std;

Int_Matrices_Var_Dim::Int_Matrices_Var_Dim()
{
  nb_rows=0;
}

void Int_Matrices_Var_Dim::Free_matrix()
{
  unsigned int n_rows=matrix.size();
  for(unsigned int row=0;row<n_rows;row++){
    matrix[row].clear();
  }

  matrix.clear();
}

void Int_Matrices_Var_Dim::Read_from_file(char *filename)
{
  
  ifstream INFILE;
  INFILE.open(filename, ios::in);
  //Checking the file path
  unsigned int n_rows=0;
  if(INFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  string line;
  while(getline(INFILE,line)){
    n_rows++;
    matrix.resize(n_rows);
    istringstream iss(line);
    unsigned int tmp_value=0;
    while(iss >> tmp_value){
      matrix[n_rows-1].push_back(tmp_value);
    }
  }
}

void Int_Matrices_Var_Dim::Display_matrix()
{
  unsigned int nb_rows=matrix.size();

  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    unsigned int nb_columns=matrix[current_line].size();
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line][current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void Int_Matrices_Var_Dim::Write_to_file(char *filename)
{
  
  ofstream OUTFILE;
  OUTFILE.open(filename, ios::out);
  //Checking the file path
  
  if(OUTFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  
  unsigned int nb_rows=matrix.size();
  //Writing the core of the matrix.
  for(unsigned int current_line=0;current_line<nb_rows;current_line++){
    unsigned int nb_columns=matrix[current_line].size();
    for(unsigned int current_column=0;current_column<nb_columns;current_column++){
      OUTFILE << matrix[current_line][current_column] << " ";
    }
    OUTFILE << endl;
  } 
  OUTFILE.close();
}
