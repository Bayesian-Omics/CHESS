/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Prior_param.cc in the ESS++ program
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

#include "Prior_param.h"
#define DEBUG 0

using namespace std;

Prior_param::Prior_param()
{
  E_p_gam=0.0;
  Sd_p_gam=0.0;
  a=0.0;
  b=0.0;
  a_pi=0.0;
  b_pi=0.0;
  k=0.0;
  delta=0.0;
  alpha=0.0;
  beta=0.0;
  w.resize(2);
}

void Prior_param::set_PR_param(double E_p_gam_from_read,
			       double Sd_p_gam_from_read,
			       unsigned int pX,
			       unsigned int pY,
			       unsigned int nX,
			       double lambda,
			       gsl_vector *vect_RMSE,
			       double P_mutation_from_read,
			       double Prob_sel_from_read,
			       double P_crsv_r_from_read,
			       double P_DR_from_read)
{
  E_p_gam = (double)(E_p_gam_from_read);
  Sd_p_gam = (double)(Sd_p_gam_from_read);
  

  a=(double)(E_p_gam)/(double)(pX);

  b=((pow((double)(Sd_p_gam),2.0)-(double)(pX)*a*(1.0 - (double)(pX)*a))/
        ((double)(pX)*(double)(pX - 1)*a));



  a_pi=((a + b - a*b - 1)/(a - b)*
           a / (1 - a));


  b_pi=((a + b - a*b -1)/
           (a - b));


  //bool stop = (a_pi>0.0 && b_pi>0.0);
  if(DEBUG){
    cout << "*****************************" << endl
         << "   Initializing PR object" << endl
         << "*****************************" << endl;

    cout << "\ta_init=" << a << "\t"
         << "b_init=" << b << "\t"
         << "a_pi_init=" << a_pi << "\t"
         << "b_pi_init=" << b_pi << "\t";

    //cout << "\tSTOP= " << stop << endl;
  }

  unsigned int count_loop=0;

//  while(!stop){
  while(a_pi <= 0 || b_pi <= 0)
  {
    count_loop++;
    if(DEBUG){
      cout << "\t\tLoop # " << count_loop << endl;
    }
    Sd_p_gam+=1;

    if(Sd_p_gam>pX){
      a_pi=1.0;
      b_pi=1.0;
      //stop=true;
      break;
    }

    b=(pow((double)(Sd_p_gam),2.0)-
       (double)(pX)*a*(1.0 - (double)(pX)*a))/
      ((double)(pX)*(double)(pX - 1)*a);
    
    a_pi=((a+b-a*b-1.0)/(a-b)*
          (a)/(1.0-a));
    b_pi=((a + b - a*b -1)/
             (a - b));

    if(DEBUG){
      cout << "\t\ta_tmp=" << a << "\t"
           << "b_tmp=" << b << "\t"
           << "a_pi_tmp=" << a_pi << "\t"
           << "b_pi_tmp=" << b_pi << endl;
    }

//    if(a_pi>0.0 && b_pi>0.0)
//    {
//      stop=true;
//    }
  }
  
  w[0] = a;
  w[1] = a_pi + b_pi;


  if(pY>1)
  {
    delta=3.0;
    gsl_sort_vector(vect_RMSE);
    if(DEBUG){    
      cout << endl << "Vect RMSE SORTED" << endl;
      display_gsl_vector(vect_RMSE);
      cout << endl;
    }
    double med_RMSE=gsl_stats_median_from_sorted_data(vect_RMSE->data,
                                                      1,
                                                      vect_RMSE->size);
    
    if(DEBUG){
      cout << "Median RMSE " << med_RMSE << endl;
    }
    k=pow(med_RMSE,2.0);

  }
  else
  {
    delta=1e-10;
    k=1e-3;
  }

  alpha=0.5;
  beta=pow((double)(nX),-lambda)/2.0;

  Prob_mut=P_mutation_from_read;
  Prob_sel=Prob_sel_from_read;
  Prob_crsv_r=P_crsv_r_from_read;
  Prob_DR=P_DR_from_read;



  if(DEBUG){
    cout << endl << "****************************************************************" << endl
         << "******************** Prior Hyper parameters ********************" << endl
         << "\t\tE_p_gam = " << E_p_gam << endl
         << "\t\tSd_p_gam = " << Sd_p_gam << endl
         << "\t\ta = " << a << endl
         << "\t\tb = " << b << endl
         << "\t\ta_pi = " << a_pi << endl
         << "\t\tb_pi = " << b_pi << endl
         << "\t\tw = [" << w[0] << " ; " << w[1] << "]" << endl
         << "\t\tk = " << k << endl
         << "\t\tdelta = " << delta << endl
         << "\t\talpha = " << alpha << endl
         << "\t\tbeta = " << beta << endl
         << "\t\tPROB = " << beta << endl
         << "\t\t\t-Prob_mut " << Prob_mut << endl
         << "\t\t\t-Prob_sel " << Prob_sel << endl
         << "\t\t\t-Prob_crsv_R " << Prob_crsv_r << endl
         << "\t\t\t-Prob_DR " << Prob_DR << endl;
    cout << endl;
    cout << "****************************************************************" << endl
         << "****************************************************************" << endl << endl;
  }
}

void Prior_param::display_prior_param()
{

  cout << endl << "****************************************************************" << endl
       << "******************** Prior hyper parameters ********************" << endl 
       << "\tE_p_gam = " << E_p_gam << endl
       << "\tSd_p_gam = " << Sd_p_gam << endl
       << "\ta = " << a << endl
       << "\tb = " << b << endl
       << "\ta_pi = " << a_pi << endl
       << "\tb_pi = " << b_pi << endl
       << "\tw = [" << w[0] << " ; " << w[1] << "]" << endl
       << "\tTEST = " << a_pi-w[0]*w[1] << " " << b_pi-(1.0-w[0])*w[1] << endl
       << "\tk = " << k << endl
       << "\tdelta = " << delta << endl
       << "\talpha = " << alpha << endl
       << "\tbeta = " << beta << endl
       << "\t##### Prob #####" << endl
       << "\t\t- Prob_mut " << Prob_mut << endl
       << "\t\t- Prob_sel " << Prob_sel << endl
       << "\t\t- Prob_crsv_r " << Prob_crsv_r << endl
       << "\t\t- Prob_DR " << Prob_DR << endl;
  cout << endl; 
  cout << "\t################" << endl
       << "****************************************************************" << endl
       << "****************************************************************" << endl << endl;
  
}
