/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from dyn_name.cc in the ESS++ program
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

#include "dyn_name.h"
#include "struc.h"
#include <string>
#include <iostream>
#include <sstream> 

string Get_dyn_name(string File_name,
		    int Number,
		    string Extension)
{
  string Result;
  ostringstream Int_to_String;

  Int_to_String << Number;

  Result= File_name+Int_to_String.str()+Extension;

  return Result;
}

string Get_stddzed_name(string File_name,
			int Number,
			string Name_number,
			string Output_type,
			string Extension)
{
  string Result;
  ostringstream Int_to_String;
  Int_to_String << Number;
  string separator="_";


  Result= File_name+separator+Int_to_String.str()+separator+Name_number+separator+Output_type+Extension;

  return Result;
}

void Write_dyn_double(data_double * object,
		      char * File_name,
		      string Name_number1 ,
		      int Number1,
		      string Name_number2 ,
		      int Number2,
		      string Extension)
{
  char *temp;
  
  string temp_res1=File_name;
  ostringstream ostr1, ostr2;
  
  ostr1 << Number1;
  ostr2 << Number2;
  
  
  string temp_res2=temp_res1+Name_number1+ostr1.str()+Name_number2+ostr2.str()+Extension;
  (temp)=(char*)(temp_res2.c_str());
  
  
  Write_matrix_double(object,temp);
  
}

void Write_dyn_double_bis(data_double * object,
			  char * File_name,
			  string Name_number1 ,
			  int Number1,
			  string Name_number2 ,
			  int Number2,
			  string Name_number3 ,
			  int Number3,
			  string Extension)
{
  char *temp;
  
  string temp_res1=File_name;
  ostringstream ostr1, ostr2, ostr3;
  
  ostr1 << Number1;
  ostr2 << Number2;
  ostr3 << Number3;
  
  
  string temp_res2=temp_res1+Name_number1+ostr1.str()+Name_number2+ostr2.str()+Name_number3+ostr3.str()+Extension;
  (temp)=(char*)(temp_res2.c_str());
  
  
  Write_matrix_double(object,temp);
  
}

void Write_matrix_double(data_double *M, char *file_name)
{
  int i;
  int j;
  FILE *flot;
  
  /* flot est un pointeur vers un fichier de type FILE*/
  
  flot=fopen(file_name,"w");
  fprintf(flot, "%d \n",(*M).nb_rows);
  fprintf(flot, "%d \n",(*M).nb_columns);
  
  for(i=0;i<(*M).nb_rows;i++) {
    for(j=0;j<(*M).nb_columns;j++) {
      fprintf(flot,"%f ",(double)((*M).matrix[i][j]));
    }
    fprintf(flot,"\n");
  }
  fclose(flot);
}
