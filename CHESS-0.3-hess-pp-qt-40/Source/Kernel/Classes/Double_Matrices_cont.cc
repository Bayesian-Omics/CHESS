/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Double_Matrices_cont.cc in the ESS++ program
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

#include "Double_Matrices_cont.h"
#include "Double_Matrices.h"
#define DEBUG 0

using namespace std;

Double_Matrices_cont::Double_Matrices_cont()
{
  nb_rows=0;
  nb_columns=0;
}

void Double_Matrices_cont::Alloc_double_matrix_cont(int Rows,
						    int Columns)
{
  
  nb_rows=Rows;
  nb_columns=Columns;
  
  matrix=new double[Rows*Columns];
}

void Double_Matrices_cont::Replicate_double_matrix_cont(Double_Matrices_cont Source)
{
  
  Alloc_double_matrix_cont(Source.nb_rows,Source.nb_columns);
  for(int row=0;row<Source.nb_rows;row++){
    for(int col=0;col<Source.nb_columns;col++){
      matrix[row*nb_columns+col]=Source.matrix[row*nb_columns+col];
    }
  }
}

void Double_Matrices_cont::Copy_from_double_matrix(Double_Matrices Source)
{
  
  Alloc_double_matrix_cont(Source.nb_rows,Source.nb_columns);
  for(int row=0;row<Source.nb_rows;row++){
    for(int col=0;col<Source.nb_columns;col++){
      matrix[row*nb_columns+col]=Source.matrix[row][col];
    }
  }
}

void Double_Matrices_cont::Free_double_matrix_cont()
{
  delete[] matrix;
}

void Double_Matrices_cont::Read_from_file(char *filename)
{
  
  ifstream INFILE;
  INFILE.open(filename, ios::in);
  //Checking the file path
  
  if(INFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  
  //  Reading the first element: nb_rows
  INFILE >> nb_rows;
  //  Reading the second element: nb_columns
  INFILE >> nb_columns;
  cout << "Reading file: " << filename
       << " -- #rows: " << nb_rows
       << " -- #cols: " << nb_columns << endl;
  //Allocating memory for the Matrix
  Alloc_double_matrix_cont(nb_rows,
			   nb_columns);
  
  //Storing the map file in the matrix.
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      if(!INFILE.eof()){
	if(DEBUG){
	  cout << " ligne " << current_line
	       << " -- colonne " << current_column
	       << " --Test " << INFILE.eof() << endl;
	}
	INFILE>>matrix[current_line*nb_columns+current_column];
      }
      else{
	cout << "Missing element at line " << current_line
	     << " and column " << current_column
	     << " in file " << filename << "." << endl
	     << "!!!!Run Stopped!!!!" << endl;
	exit(1);
      }
    }
  } 
}

void Double_Matrices_cont::Display_matrix()
{
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line+nb_columns+current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void Double_Matrices_cont::Display_matrix_header()
{
  
  cout << nb_rows << endl;
  cout << nb_columns << endl;
  
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      cout << matrix[current_line*nb_columns+current_column] << " ";
    }
    cout << endl;
  }
  cout << endl;
  
}

void Double_Matrices_cont::Write_to_file(char *filename)
{
  
  ofstream OUTFILE;
  OUTFILE.open(filename, ios::out);
  //Checking the file path
  
  if(OUTFILE.fail()){
    cout << "Invalid Path and/or permission rights for " << filename << " -- run stopped." << endl;
    exit(1);
  }
  
  //Writing the first element: nb_rows
  OUTFILE << nb_rows << endl;
  //Writing the second element: nb_columns
  OUTFILE << nb_columns << endl;
  
  //Writing the core of the matrix.
  for(int current_line=0;current_line<nb_rows;current_line++){
    for(int current_column=0;current_column<nb_columns;current_column++){
      OUTFILE << matrix[current_line*nb_columns+current_column] << " ";
    }
    OUTFILE << endl;
  } 
  OUTFILE.close();
}
