/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from AdMH.cc in the ESS++ program
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


#include "AdMH.h"
#define DEBUG 0

using namespace std;

AdMH::AdMH()
{
  n_batch=0;
  optimal=0.0;
  ls=0.0;
  delta_n=0.0;
  tilda_accept=0.0;
  tilda_accept_ins=0.0;
  tilda_n_sweep=0;
  tilda_n_sweep_ins=0;
  M.resize(2);

}
  
void AdMH::set_AdMH(int g_sample,
			unsigned int n_batch_from_read,
            double AdMH_optimal_from_read,
            double AdMH_ls_from_read,
			unsigned int pX,
			unsigned int burn_in,
            double M_min_input,
            double M_max_input)
{

  if(g_sample==1){
    n_batch= n_batch_from_read;
    optimal=AdMH_optimal_from_read;
    ls=AdMH_ls_from_read;

    if(fabs(M_min_input-0)<1e-10){
      M[0]=-0.5*log(pX);
    }
    else{
      M[0]=M_min_input;
    }

    if(fabs(M_max_input-0)<1e-10){
      M[1]=0.5*log(pX);
    }
    else{
      M[1]=M_max_input;
    }
    double temp1=fabs(M[0]-ls);
    double temp2=fabs(M[1]-ls);
    delta_n=max(temp1,temp2);
    delta_n/=(double)(burn_in)/(double)(n_batch);
 
  }
}

void AdMH::display_AdMH()
{

  cout << endl << "**********************************************************" << endl
       << "******************** AdMH parameters ********************" << endl 
       << "\tn_batch = " << n_batch << endl
       << "\toptimal = " << optimal << endl
       << "\tls = " << ls << endl
       << "\tM[0] = " << M[0] << " -- " << "M[1] = " << M[1] << endl
       << "\tdelta_n = " << delta_n << endl
       << "\tG_tilda_accept = " << tilda_accept << endl
       << "\tG_tilda_accept_ins = " << tilda_accept_ins << endl
       << "\tG_tilda_n_sweep = " << tilda_n_sweep << endl
       << "\tG_tilda_n_sweep_ins = " << tilda_n_sweep_ins << endl
       << "\tLs " << endl << "\t";
  for(unsigned int col=0;col<Ls.size();col++){
    cout << Ls[col] << " ";
  }
  cout << endl
       << "**********************************************************" << endl
       << "**********************************************************" << endl << endl;
  
}

void AdMH::Adapt_AdMH_ESS(unsigned int sweep, Move_monitor *My_Move_monitor)
{

    double Old_ls=ls;

    double local_accept_rate=Adapt_AdMH(sweep);

    (*My_Move_monitor).g_adapt_history[0].push_back((double)(sweep));
    (*My_Move_monitor).g_adapt_history[1].push_back(local_accept_rate);
    (*My_Move_monitor).g_adapt_history[2].push_back(Old_ls);


}

double AdMH::Adapt_AdMH(unsigned int sweep)
{
    if(DEBUG){
      cout << "G-adaptation" << endl;
      cout << "sweep=" << sweep
       << " -- My_g_AdMH.n_batch " << n_batch
       << " -- Test " << sweep%n_batch << endl;
    }
    double local_accept_rate=(double)(tilda_accept_ins)/(double)(tilda_n_sweep_ins);

    double batchNumber = (double)sweep/(double)n_batch;

    double delta =  delta_n;
    if(local_accept_rate<optimal)
    {
      //Updating AdMH.ls
      ls-=delta;
      if(ls < M[0])
      {
        ls = M[0];
      }
    }
    else{
      ls+=delta;
      if(ls > M[1])
      {
        ls = M[1];
      }
    }
    tilda_accept_ins=0;
    tilda_n_sweep_ins=0;
    if(DEBUG){
      cout << "Sweep " << sweep << " -- Update the AdMH parameters" << endl;
      cout << "Updated g_ls " << ls << endl;
    }

    return local_accept_rate;

}
