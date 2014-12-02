/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Temperatures.cc in the ESS++ program
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

#include "Temperatures.h"
#define DEBUG 0

using namespace std;

Temperatures::Temperatures()
{
  C=0;
  b_t=0.0;
  a_t_den=0.0;
  //a_t.resize(2);
  //t.resize(2);
  nbatch=1;
  M.resize(2);
  delta_n=0.0;
  optimal=0.0;
  c_idx=0;
}

void Temperatures::set_Temp_param(unsigned int nb_chains,
				  unsigned int pX,
				  double b_t_input,
				  double a_t_den_inf_5k,
				  double a_t_den_5_10k,
				  double a_t_den_sup_10k,
				  unsigned int nbatch_input,
				  vector < double > &M_input,
				  unsigned int burn_in,
				  double optimal_input,
				  bool iso_T_Flag)
{
  C=nb_chains;

  if(nb_chains==1){
    t.resize(1);
    t[0]=1.0;
    c_idx=0;
  }
  else{
    b_t=b_t_input;
    if(pX<=5000){
      a_t_den= a_t_den_inf_5k;
    }
    else if(pX>5000 && pX<=10000){
      a_t_den=a_t_den_5_10k;
    }
    else{//pX>10,000
      a_t_den= a_t_den_sup_10k;
    }

    a_t.resize(nb_chains);
    t.resize(nb_chains);

    for(unsigned int chain=0;chain<nb_chains;chain++){
      a_t[chain]=(double)(chain)/a_t_den;
      if(!iso_T_Flag){
	t[chain]=pow(b_t,a_t[chain]);
      }
      else{
	t[chain]=1;
      }
    }
    nbatch=nbatch_input;
    unsigned int tmp_M_input_size=M_input.size();
    if(tmp_M_input_size==2){
      for(unsigned int row=0;row<M_input.size();row++){
	M[row]=M_input[row];
      }
    }
    else{
      cout << "Dimension mismatch M_input has " <<  tmp_M_input_size << " --elements, " << 2 << " expected." << endl
	   << "Run stopped." << endl;
      exit(1);
    }
    double tmp1=fabs(log(M[0]/b_t)/log(2));
    double tmp2=fabs(log(M[1]/b_t)/log(2));
    double numer=max(tmp1,tmp2);
    delta_n=numer/((double)(burn_in)/nbatch);
    optimal=optimal_input;
    for(unsigned int chain=0;chain<nb_chains;chain++){
      if(t[chain]==1.0){
	c_idx=chain;
      }
    }
  }
  
}

void Temperatures::display_Temp_param()
{

  cout << endl << "**********************************************************" << endl
       << "******************** Temp parameters ********************" << endl 
       << "\tC = " << C << endl;
  if(C>1){
    cout << "\tb_t = " << b_t << endl
	 << "\ta_t_den = " << a_t_den << endl
	 << "\tnbatch = " << nbatch << endl
	 << "\tdelta_n = " << delta_n << endl
	 << "\toptimal = " << optimal << endl
	 << "\tc_idx = " << c_idx << endl
	 << "\ta_t: [";
    for(unsigned int chain=0;chain<C;chain++){
      cout << a_t[chain] << " ";
    }
    cout << "]" << endl;
  }
  cout << "\tt: [";
  for(unsigned int chain=0;chain<C;chain++){
    cout << t[chain] << " ";
  }
  cout << "]" << endl;
  if(C>1){   
    cout << "\tM: [";
    for(unsigned int chain=0;chain<M.size();chain++){
      cout << M[chain] << " ";
    }
    cout << "]" << endl;
  }
  cout << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
  
}

