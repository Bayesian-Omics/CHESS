/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from moves.cc in the ESS++ program
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

#include "moves.h"
#include "struc.h"

#include "../../General_Classes/Model_Information.h"
#include "../../General_Classes/Kernel_Single_Gamma.h"

#define DEBUG 0
#define DEBUG_HESS 0

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF


using namespace std;

void Gibbs_move(unsigned int Response,
        Double_Matrices mat_log_marg,
		Double_Matrices mat_log_cond_post,
		gsl_matrix *mat_X,
		gsl_matrix *mat_Y,
		Temperatures *t_tun,
		gsl_permutation *MyPerm,
		vector < vector <unsigned int> > &vect_gam,
		unsigned int sweep,
		bool gPriorFlag,
		bool indepPriorFlag,
		bool gSampleFlag,
		double lambda,
		double g,
		Prior_param PR,
		Move_monitor *My_Move_monitor,
		vector <unsigned int > &chain_idx,
		vector <unsigned int > &n_Models_visited,
        bool cudaFlag,
        unsigned int nConfounders,
        unsigned int maxPX,
        gsl_rng *RandomNumberGenerator,
        Model_Information &Model_Data)
 {
    unsigned int pX=mat_X->size2;
    gsl_vector *prop_log_marg_condPost=gsl_vector_calloc(2);
    vector < unsigned int > list_columns_X_gam;
    unsigned int count_0_1=0;
    unsigned int count_1_0=0;
    unsigned int count_unchanged=0;

    (*My_Move_monitor).Gibbs_nb_sweep++;
    (*My_Move_monitor).Gibbs_move_history[0].push_back(sweep);
    unsigned int pos_current_chain=chain_idx[0];
    if(DEBUG)
    {
        cout << "!!!!!!!!!!!!" << endl;
        cout << "CHAIN " << 0
        << " --pos " << pos_current_chain << endl;
        cout << "!!!!!!!!!!!!" << endl << endl;
    }
    //Step1: get the random order of variable: shuffle the permutation object
    My_Permut_unsigned_int(MyPerm,RandomNumberGenerator);
    //display_gsl_perm(MyPerm);
    if(DEBUG)
    {
        cout << endl;
    }
    double current_log_cond_post=mat_log_cond_post.matrix[pos_current_chain][sweep];
    double current_log_marg=mat_log_marg.matrix[pos_current_chain][sweep];

    //Step 2: for each variable calculating the moving pbty

    for(unsigned int current_variable=0;current_variable<MyPerm->size;current_variable++){

    (*My_Move_monitor).Gibbs_nb_model++;

    unsigned int pos_curr_var=MyPerm->data[current_variable];
    if(pos_curr_var<nConfounders)
    {
        // confounders don't get updated as always in
        continue;
    }
    unsigned int current_value=vect_gam[pos_current_chain][pos_curr_var];
    if(DEBUG){
        cout << "Initial gamma:" << endl;
        for(unsigned int col=0;col<vect_gam[pos_current_chain].size();col++)
        {
            if(vect_gam[pos_current_chain][col]==1){
                cout << col << " ";
            }
        }
        cout << endl;
        cout << "Variable rank " << current_variable
        << " -- position " << pos_curr_var
        << " -- current_value " << current_value
        << " -- current_log_cond_post " << current_log_cond_post;

    }

    //Step 3: Calculate the alternative long_cond_post
    //Step 3.1: change Gamma, and get the new proposed X_gam
    vect_gam[pos_current_chain][pos_curr_var]=1-vect_gam[pos_current_chain][pos_curr_var];
    if(DEBUG){
        cout << " -- modified value " << vect_gam[pos_current_chain][pos_curr_var] << endl;
    }

    //Initial Calculation of the logMarg and log_cond_post;

    //Step 1: Getting X_gamma
    get_list_var_in(list_columns_X_gam,
        vect_gam[pos_current_chain]);

    unsigned int n_vars_in=list_columns_X_gam.size();

    double propLogLik,propLogPost;
    double log_condPost0,log_condPost1;
    double current_t,argument,theta;
    if(n_vars_in<=maxPX){
          computeLogPosterior(vect_gam[pos_current_chain],0,Response,
                              mat_X,propLogLik,propLogPost,mat_Y,PR,gPriorFlag,
                              indepPriorFlag,gSampleFlag,lambda,g,
                              //pX,nX,pY,
                              cudaFlag,
                              Model_Data);
        prop_log_marg_condPost->data[0]=propLogLik;
        prop_log_marg_condPost->data[1]=propLogPost;


        //Step 3: Sample the new status
        // We prevent n_var_in moving to more than maxPX later below
        if(DEBUG){
            cout << "\tCurrent log_marg' " << current_log_marg
                 << " -- cond_post " << current_log_cond_post << endl;
            cout << "\tProposed cond_post " << prop_log_marg_condPost->data[1]
                 << " -- log_marg' " << prop_log_marg_condPost->data[0] << endl;
        }
        log_condPost0=0.0;
        log_condPost1=0.0;

        //The sampling theta = 1/(1+exp(condpost(0)-condpost(1)))
        if(current_value==0){//i.e. 0->1 move
            log_condPost0=current_log_cond_post;
            log_condPost1=prop_log_marg_condPost->data[1];
        }
        else{//i.e. 1->0 move
            log_condPost1=current_log_cond_post;
            log_condPost0=prop_log_marg_condPost->data[1];
        }
        current_t=(*t_tun).t[0];
        argument=(log_condPost0-log_condPost1)/(current_t);
        theta=1.0/(1.0+exp(argument));
    }

    unsigned int sampled_value;
    if(current_value==0&&n_vars_in>maxPX){
        // Stop from adding more variables than permitted (as this would have 0 prob)
        sampled_value=0;
    }else{
        // Otherwise perform the sample
        sampled_value=(unsigned int)(genBernoulli(theta,RandomNumberGenerator));
    }
    if(DEBUG){
        cout << "\tTheta " << theta
        << " -- sampled gamma " << sampled_value
        << " -- current value " << current_value
        << " -- current_T " << current_t
        << endl;
    }
    //Step 4: updating vectors
     if(current_value==sampled_value){
        if(DEBUG){
        cout << "\tNO CHANGE SAMPLED" << endl;
        }
        //Revert vect_gam
        vect_gam[pos_current_chain][pos_curr_var]=1-vect_gam[pos_current_chain][pos_curr_var];
        count_unchanged++;
    }
    else{
        if(DEBUG){
            cout << "\tMOVE SAMPLED: RETAINED GAMMA" << endl;
            for(unsigned int col=0;col<vect_gam[pos_current_chain].size();col++){
                if(vect_gam[pos_current_chain][col]==1){
                cout << col << " ";
                }
            }
            cout << endl;
        }
        //updating log_marg and log_condPost
        current_log_cond_post=prop_log_marg_condPost->data[1];
        current_log_marg=prop_log_marg_condPost->data[0];
        if(current_value==0){
            count_0_1++;
            (*My_Move_monitor).Gibbs_nb_0_1++;
        }
        else{
            count_1_0++;
            (*My_Move_monitor).Gibbs_nb_1_0++;
        }
    }

    list_columns_X_gam.clear();

    }//end of for variable
    mat_log_cond_post.matrix[pos_current_chain][sweep]=current_log_cond_post;
    mat_log_marg.matrix[pos_current_chain][sweep]=current_log_marg;


    if(DEBUG){
        unsigned int final_n_vars_in=sum_line_std_mat(vect_gam,
                          pos_current_chain);
        cout << "**********************************************************************" << endl;
        cout << "End of the variables, Final Results for chain 1 -- pos " << pos_current_chain << endl;
        cout << "\tlog_cond_post " << mat_log_cond_post.matrix[pos_current_chain][sweep]
        << " -- log_marg " << mat_log_marg.matrix[pos_current_chain][sweep]
        << " -- # var in " << final_n_vars_in << endl;
        cout << "\t";

        for(unsigned int col=0;col<vect_gam[pos_current_chain].size();col++){
            if(vect_gam[pos_current_chain][col]==1){
            cout << col << " ";
            }
        }
        cout << endl;
        cout << "**********************************************************************" << endl << endl;
    }
    (*My_Move_monitor).Gibbs_move_history[1].push_back(count_0_1);
    (*My_Move_monitor).Gibbs_move_history[2].push_back(count_1_0);
    (*My_Move_monitor).Gibbs_move_history[3].push_back(count_unchanged);
    n_Models_visited[sweep]+=pX;

    gsl_vector_free(prop_log_marg_condPost);
}

void intialize_chain_idx(vector <unsigned int > &chain_idx,
			 unsigned int nb_chains)
{
  chain_idx.resize(nb_chains);

  for(unsigned int curr_chain=0;curr_chain<nb_chains;curr_chain++){
    chain_idx[curr_chain]=curr_chain;
  }
}

void All_exchange_move(gsl_matrix *description_exchange_move,
		       vector <unsigned int > &chain_idx,
		       Double_Matrices mat_log_cond_post,
		       Temperatures *t_tun,
		       unsigned int sweep,
		       Move_monitor *My_Move_monitor,
               gsl_rng *RandomNumberGenerator,
               Model_Information &Model_Data)
{

  //cout << "Description exchange move init" << endl;
  //display_gsl_matrix(description_exchange_move);
  if(DEBUG){
    cout << "//////////////////////////////////////////" << endl
	 << "//  All exchange move" << endl
	 << "//////////////////////////////////////////" << endl << endl;
  }


  double cum_pbty=1.0;
  unsigned int rank=1;
  unsigned int nb_chains=chain_idx.size();
  unsigned int nb_cols=description_exchange_move->size2;
  description_exchange_move->data[2*description_exchange_move->size2]=1.0;
  (*My_Move_monitor).All_exchange_nb_sweep++;
  //Step 1: Calculating the exchange move pbties
  for(unsigned int c1=0;c1<nb_chains-1;c1++){
    for(unsigned int c2=c1+1;c2<nb_chains;c2++){
      unsigned int pos_c1=chain_idx[c1];
      unsigned int pos_c2=chain_idx[c2];
      // WARNING the last log_cond_post is taken into account
      // i.e. all echange always after local move
      double argument=((mat_log_cond_post.matrix[pos_c2][sweep]-mat_log_cond_post.matrix[pos_c1][sweep])*
		       (1.0/(*t_tun).t[c1] - 1.0/(*t_tun).t[c2]));
      double pbty=exp(argument);
      if(isnan(pbty)==1){
	pbty=0.0;
      }
      description_exchange_move->data[2*nb_cols+rank]=pbty;
      cum_pbty+=pbty;
      if(DEBUG){
	cout << "c1 " << c1 << " -- c2 " << c2
	     << " -- pos_c1 " << pos_c1 << " -- pos_c2 " << pos_c2
	     << " -- pbty " << pbty << " -- cum_pbty " << cum_pbty
	     << " -- Test1 " << description_exchange_move->data[2*description_exchange_move->size2+rank] << endl;
      }
      rank++;    
    }
  }
  //Step 2: Standardizing the pbty vector
  for(unsigned int col=0;col<description_exchange_move->size2;col++){
    description_exchange_move->data[2*nb_cols+col]/=cum_pbty;
    if(DEBUG){
      cout << "col " << col << " -- cum " << cum_pbty
	   << " -- pbty" << description_exchange_move->data[2*description_exchange_move->size2+col] << endl;
    }
  }
  //Step 3: Sampling the exchange move
  int my_sample=SampleFromDiscrete_All_exchange(description_exchange_move,RandomNumberGenerator);

  //step 4: Updating chain_idx

  unsigned int c1=description_exchange_move->data[my_sample];
  unsigned int c2=description_exchange_move->data[nb_cols+my_sample];
  if(DEBUG){
    cout << "Sampled exchange move " << my_sample
	 << " -- C1 " << c1
	 << " -- C2 " << c2 << endl;
  }
  (*My_Move_monitor).All_exchange_move_history[0].push_back(sweep);
  (*My_Move_monitor).All_exchange_move_history[1].push_back(c1);
  (*My_Move_monitor).All_exchange_move_history[2].push_back(c2);
  if(my_sample>0){
    (*My_Move_monitor).All_exchange_n_accept++;
    (*My_Move_monitor).All_exchange_freq[c1]++;
    (*My_Move_monitor).All_exchange_freq[c2]++;
    unsigned int pos_c1_i=chain_idx[c1];
    unsigned int pos_c2_i=chain_idx[c2];
    if(DEBUG){
      cout << " -- pos_c1_i " << pos_c1_i
     	   << " -- pos_c2_i " << pos_c2_i;
    }
    chain_idx[c1]=pos_c2_i;
    chain_idx[c2]=pos_c1_i;
  }
  if(DEBUG){
    cout << endl;
    cout << "vect_idx" << endl;
    for(unsigned int col=0;col<chain_idx.size();col++){
      cout << chain_idx[col] << " ";
    }
    cout << endl;
  }
  //cout << "Description exchange move final" << endl;
  //display_gsl_matrix(description_exchange_move);
  if(DEBUG){
    cout << "//////////////////////////////////////////" << endl
	 << "// end of All Exchange move" << endl
	 << "//////////////////////////////////////////" << endl << endl;
  }

}

void DR_move(DR *My_DR,
	     vector <unsigned int > &chain_idx,
	     Double_Matrices mat_log_cond_post,
	     Temperatures *t_tun,
	     unsigned int sweep,
	     Move_monitor *My_Move_monitor,
         gsl_rng *RandomNumberGenerator,
         Model_Information &Model_Data)
{
  (*My_Move_monitor).DR_nb_sweep++;
  (*My_DR).nb_calls++;
  if(DEBUG){
    cout << "//////////////////////////////////////////" << endl
	 << "// DR move" << endl
	 << "//////////////////////////////////////////" << endl << endl;
  }
  
   //Step 1: sampling uniformly two different chains:
  unsigned int nb_chains=chain_idx.size();
  unsigned int c_1=0;
  unsigned int c_2=0;
  unsigned int c_1_adj=0;
  unsigned int c_2_adj=0;
  int delta=0;

  unsigned int stop_sampling=0;
  while(stop_sampling==0){
    double rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
    c_1=(unsigned int)(rand_test);
    rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
     c_2=(unsigned int)(rand_test);
    int delta=c_1-c_2;
    if(delta!=0){
      stop_sampling=1;
    }
  }
  unsigned int pos_c1=chain_idx[c_1];
  unsigned int pos_c2=chain_idx[c_2];

  //Step2: Calculating the move Probability
  double argument=((mat_log_cond_post.matrix[pos_c2][sweep]-mat_log_cond_post.matrix[pos_c1][sweep])*
		   (1.0/(*t_tun).t[c_1] - 1.0/(*t_tun).t[c_2]));
  double pbty_1=min(1.0,exp(argument));

  //Step3: Sampling Acceptance for the first proposed move.
  double rand_accept_1=myrand(RandomNumberGenerator);

  
  if(DEBUG){
    cout << "First Round; proposed move:" << endl;
    cout << "cl " << c_1 << " -- cr " << c_2
	 << " -- pos_cl " << pos_c1 << " -- pos_cr " << pos_c2
	 << " -- pbty_1 " << pbty_1 
	 << " -- rand_1 " << rand_accept_1
	 << endl;
  }
  
  if(isnan(pbty_1)==1 || pbty_1==0.0){
    
    if(DEBUG){
      cout << "pbty_1 is nan or null -- move rejected" << endl;
    }
  }
  else{//pbty_1 is numeric and >0.0
    if(rand_accept_1<pbty_1){//Accept c_1 <-> c_2 move
      if(DEBUG){
	cout << "First proposed move ACCEPTED" << endl;
      }
      delta=c_1-c_2;
      if(abs(delta)==1){
	(*My_DR).nb_calls_adj++;
      }
      (*My_DR).mat_moves_accepted[c_1][c_2]++;
      (*My_DR).mat_moves_accepted[c_2][c_1]++;
      (*My_DR).mat_moves_proposed[c_1][c_2]++;
      (*My_DR).mat_moves_proposed[c_2][c_1]++;
      (*My_Move_monitor).DR_move_history[0].push_back(sweep);
      (*My_Move_monitor).DR_move_history[1].push_back(c_1);
      (*My_Move_monitor).DR_move_history[2].push_back(c_2);
      (*My_Move_monitor).DR_n_accept++;
      (*My_Move_monitor).DR_freq[c_1]++;
      (*My_Move_monitor).DR_freq[c_2]++;
      unsigned int pos_c1_i=chain_idx[c_1];
      unsigned int pos_c2_i=chain_idx[c_2];  
      if(DEBUG){
	cout << "pos_c1_i " << pos_c1_i
	     << " -- pos_c2_i " << pos_c2_i << endl;
      }
      chain_idx[c_1]=pos_c2_i;
      chain_idx[c_2]=pos_c1_i;
      if(DEBUG){
	cout << "Updated vect_idx" << endl;
	for(unsigned int col=0;col<chain_idx.size();col++){
	  cout << chain_idx[col] << " ";
	}
	cout << endl;
      }
    }
    else{
      if(DEBUG){
	cout << "First proposed move Rejected: " << endl;
      }
      //Re-sampling tow adjascent chains
      double rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
      c_1_adj=(unsigned int)(rand_test);
      rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
      c_2_adj=(unsigned int)(rand_test);
      delta=c_1_adj-c_2_adj;
      if(abs(delta)!=1){
	if(c_1_adj==0){//if first chain sampled, the other one is second
	  c_2_adj=1;
	}
	else if(c_1_adj==nb_chains-1){//if last chain sampled, the other one is next to last
	  c_2_adj=c_1_adj-1;
	}
	else{//otherwise one takes with probability 50% the one below or above
	  double rand_bin=myrand(RandomNumberGenerator);
	  if(rand_bin<0.5){
	    c_2_adj=c_1_adj-1;
	  }
	  else{
	    c_2_adj=c_1_adj+1;
	  }
	}
      }
      (*My_DR).nb_calls_adj++;
      
      if(DEBUG){
	cout << "Second Round; proposed move:" << endl;
	cout << "c1_adj " << c_1_adj << " -- c2_adj " << c_2_adj
	     << endl;
      }
      (*My_DR).mat_moves_proposed[c_1_adj][c_2_adj]++;
      (*My_DR).mat_moves_proposed[c_2_adj][c_1_adj]++;
  
      unsigned int pos_c1_adj=chain_idx[c_1_adj];
      unsigned int pos_c2_adj=chain_idx[c_2_adj];
      double arg_adj=((mat_log_cond_post.matrix[pos_c2_adj][sweep]-mat_log_cond_post.matrix[pos_c1_adj][sweep])*
		      (1.0/(*t_tun).t[c_1_adj] - 1.0/(*t_tun).t[c_2_adj]));
      
      
      //Sampling two other candidates:
      stop_sampling=0;
      while(stop_sampling==0){
	rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
	c_1=(unsigned int)(rand_test);
	rand_test=myrand(RandomNumberGenerator)*(double)(nb_chains);
	c_2=(unsigned int)(rand_test);
	int delta=c_1-c_2;
	if(delta!=0){
	  stop_sampling=1;
	}
	
      }
      pos_c1=chain_idx[c_1];
      pos_c2=chain_idx[c_2];
    
      double arg_num=((mat_log_cond_post.matrix[pos_c2][sweep]-mat_log_cond_post.matrix[pos_c1][sweep])*
		      (1.0/(*t_tun).t[c_1] - 1.0/(*t_tun).t[c_2]));
      
      double pbty_2=min(1.0,(exp(arg_num)/exp(arg_adj)));//alpha_DR_star=alpha_DR_num
      
      double tmp_pbty_final= exp(arg_adj)*((1.0-pbty_2)/(1.0-pbty_1));
      double pbty_final=min(1.0,tmp_pbty_final);
      double rand_accept_2=myrand(RandomNumberGenerator);
      
      if(DEBUG){
	cout << "c_l_star " << c_1
	     << " -- c_r_star " << c_2 
	     << " -- pos_cl_star " << pos_c1
	     << " -- pos_cr_star " << pos_c2
	     << " -- alpha_DR_star " << pbty_2
	     << endl
	     << "c_ll " << c_1_adj
	     << " -- c_rr " << c_2_adj
	     << " -- pos_cll " << pos_c1_adj
	     << " -- pos_crr " << pos_c2_adj
	     << endl
	     << "arg_num " << arg_num
	     << " -- alpha_DR_star " << pbty_2
	     << " -- arg_adj " << arg_adj
	     << " -- tmp_pbty_final " << tmp_pbty_final
	     << " -- alpha_DR_2 " << pbty_final 
	     << " -- rand_2 " << rand_accept_2
	     << endl;
      }
      if(isnan(pbty_final)==1){
	if(DEBUG){
	  cout << "Pbty for second DR move is nan: move rejected" << endl;
	}
      }
      else{
	if(rand_accept_2<pbty_final){
	  if(DEBUG){
	    cout << "Second proposed move ACCEPTED" << endl;
	  }
	  (*My_DR).mat_moves_accepted[c_1_adj][c_2_adj]++;
	  (*My_DR).mat_moves_accepted[c_2_adj][c_1_adj]++;

	  (*My_Move_monitor).DR_move_history[0].push_back(sweep);
	  (*My_Move_monitor).DR_move_history[1].push_back(c_1_adj);
	  (*My_Move_monitor).DR_move_history[2].push_back(c_2_adj);
	  (*My_Move_monitor).DR_n_accept++;
	  (*My_Move_monitor).DR_freq[c_1_adj]++;
	  (*My_Move_monitor).DR_freq[c_2_adj]++;
	  unsigned int pos_c1_i=chain_idx[c_1_adj];
	  unsigned int pos_c2_i=chain_idx[c_2_adj];  
	  if(DEBUG){
	    cout << "pos_c1_i " << pos_c1_i
		 << " -- pos_c2_i " << pos_c2_i << endl;
	  }
	  chain_idx[c_1_adj]=pos_c2_i;
	  chain_idx[c_2_adj]=pos_c1_i;
	  if(DEBUG){
	    cout << "Updated vect_idx" << endl;
	    for(unsigned int col=0;col<chain_idx.size();col++){
	      cout << chain_idx[col] << " ";
	    }
	    cout << endl;
	  }
	}
	else{
	  if(DEBUG){
	    cout << "Second proposed move Rejected: " << endl;
	  }
	}
      }
    }//end of else: second move tried;
  }//pbty_1 is >0.0 and numeric
    
  if(DEBUG){
    cout << "//////////////////////////////////////////" << endl
	 << "// End of DR move" << endl
	 << "//////////////////////////////////////////" << endl << endl;
  }
  
}

gsl_matrix *description_exch_moves(unsigned int nb_chains)
{
  unsigned int nb_columns=(nb_chains*(nb_chains - 1)/2) + 1;
  gsl_matrix *result=gsl_matrix_calloc(3,nb_columns); 
  unsigned int rank=1;
  for(unsigned int c1=0;c1<nb_chains-1;c1++){
    for(unsigned int c2=c1+1;c2<nb_chains;c2++){
      result->data[rank]=c1;
      result->data[nb_columns+rank]=c2;
      rank++;
    }    
  }
  return result;
}

/*
void sample_g_HESS(unsigned int n_chains,
          //AdMH *My_g_AdMH,
          vector <AdMH *> My_g_AdMH,
          vector <Kernel_Single_Gamma > Gammas,
          gsl_matrix *mat_X,
          vector <gsl_matrix *> &Vector_Data_Y,
          bool gPriorFlag,
          bool indepPriorFlag,
          bool gSampleFlag,
          double lambda,
          //double &g,
          vector <double> &g,
          vector <Prior_param> PR_per_Resp,
          unsigned int sweep,
              bool cudaFlag,
              gsl_rng *RandomNumberGenerator,
              Model_Information &Model_Data,
          bool Single_g)
{

    //cout << "Just entered sample_g_HESS" << endl;


  bool DEBUG_sample_g=false;
  /*
  cout << "in Sample_g_HESS:  DEBUG_sample_g = " << DEBUG_sample_g << endl;
  cout << endl;

  cout << "X(1,1) = " << gsl_matrix_get(mat_X,0,0) << endl;
  cout << "Y_1(1,1) = " << gsl_matrix_get(Vector_Data_Y[0],0,0) << endl;

  cout << endl;

  cout << "omega(1,1) = " << (*(Model_Data.Omega_PerLine))[0][0] << endl;
  cout << "rho(1,1) = " << (*(Model_Data.Omega_PerLine))[0][0] << endl;

  cout << endl;
  *

  if(DEBUG_HESS){
    cout << "//////////////////////////////////////////" << endl
     << "//  Adaptive M-H: sampling g" << endl
     << "//////////////////////////////////////////" << endl << endl;
  }

  double norm_mean;
  double norm_sd;

  double sample_tmp;
  double g_prop;

  double alpha_g;

  double logPG;
  double logPG_Prop;

  if (Single_g==true)
  {
      norm_mean=log(g[0]);
      norm_sd=exp((*My_g_AdMH[0]).ls);

      //cout << "mean proposal: " << norm_mean << endl;
      cout << "sd proposal: " << norm_sd << endl;

      sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);
      g_prop=exp(sample_tmp);

      if (DEBUG_sample_g==true)
      {
        g_prop=g[0]+norm_sd;
      }

  }

  vector < vector < double > > store_log_marg;
  vector < vector < double > > store_log_cond_post;

  unsigned int n_Responses=Vector_Data_Y.size();

  store_log_marg.resize(n_Responses);
  store_log_cond_post.resize(n_Responses);

  double cum_diff;
  if (Single_g==true)
      cum_diff=0.0;


  //Loop on responses
  for(unsigned int k=0;k<n_Responses;k++)
  {
      //cout << "k=" << k << endl;

      if(DEBUG_HESS)
      {
          cout << "***************************" << endl;
          cout << "Response Nb " << k << endl;
          cout << "***************************" << endl;
      }


      if (Single_g==false)
      {
          norm_mean=log(g[k]);
          norm_sd=exp((*My_g_AdMH[k]).ls);

          sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);
          g_prop=exp(sample_tmp);

          cum_diff=0.0;
      }

      for(unsigned int current_chain=0;current_chain<n_chains;current_chain++)
      {
        unsigned int pos_current_chain=Gammas[k].chain_idx[current_chain];
        if(DEBUG_HESS){
          cout << "CHAIN " << current_chain
           << " -- pos " << pos_current_chain << endl;
        }

        double current_log_cond_post=(Gammas[k].mat_log_cond_post.matrix)[pos_current_chain][sweep];
        double current_log_marg=(Gammas[k].mat_log_marg.matrix)[pos_current_chain][sweep];

        gsl_matrix *Data_Y_k=Vector_Data_Y[k];

        double propLogMarg,propLogPost;
        if(DEBUG_HESS)
        {
            cout << current_chain << " ";
        }

        //unsigned int pY=Data_Y_k->size2;
        computeLogPosterior(Gammas[k].vect_gam[pos_current_chain],current_chain,k,
                            mat_X,propLogMarg,propLogPost,Data_Y_k,
                            PR_per_Resp[k],
                            gPriorFlag,
                            indepPriorFlag,gSampleFlag,lambda,
                            g_prop,
                            cudaFlag,
                            Model_Data);
        //propLogMarg = logP(Y_k | Gamma_k,g)
        //propLogPost = log[P(Y_k | Gamma_k,g).P(Gamma_k | Omega).P(g)]

        Gammas[k].n_Models_visited[sweep]++;

        cum_diff+=(propLogMarg-current_log_marg) / (*(Gammas[k].t_tun)).t[current_chain];

        /*
          cout << "k =  " <<  k+1 //<< ", Temp = " << (*(Gammas[k].t_tun)).t[current_chain]
                  << ", curr_log_marg = " << current_log_marg
                  << ", propLogMarg = " << propLogMarg
                  << ", Delta = " << propLogMarg-current_log_marg
                  <<", cum_diff = " << cum_diff << endl;
        *


        store_log_marg[k].push_back(propLogMarg);
        store_log_cond_post[k].push_back(propLogPost);

        if(DEBUG_HESS){
          cout << "\tCurrent log_marg' " << current_log_marg
           << " -- cond_post " << current_log_cond_post << endl;
          cout << "\tProposed log_marg' " << propLogMarg << endl;
           cout << " -- cond_post " << propLogPost << endl;
          cout << "\tcum_diff " << cum_diff << endl;
        }

      }//end of for chains

      if (Single_g==false)
      {
          logPG=getPriorG(PR_per_Resp[k],gSampleFlag,g[k]);
          logPG_Prop=getPriorG(PR_per_Resp[k],gSampleFlag,g_prop);

          cum_diff=cum_diff+logPG_Prop-logPG+
                            log(g_prop)-log(g[k]);

          alpha_g=min(1.0,exp(cum_diff));

          if(!isnan(alpha_g))
          {
            double rand_test = myrand(RandomNumberGenerator);
            if(rand_test<alpha_g && g_prop>0.0)
            {
              g[k]=g_prop;

              for(unsigned int chain=0;chain<n_chains;chain++)
              {
                unsigned int pos_chain=Gammas[k].chain_idx[chain];
                Gammas[k].mat_log_cond_post.matrix[pos_chain][sweep]=store_log_cond_post[k][chain];
                Gammas[k].mat_log_marg.matrix[pos_chain][sweep]=store_log_marg[k][chain];
              }
              (*My_g_AdMH[k]).tilda_accept++;
              (*My_g_AdMH[k]).tilda_accept_ins++;
            }
          }else{
              if(DEBUG_HESS){cout << "alpha_g is NaN" << endl;}
          }

          (*My_g_AdMH[k]).tilda_n_sweep++;
          (*My_g_AdMH[k]).tilda_n_sweep_ins++;

      }// End sampling g[k]


      //Removed:
      //gsl_vector_free(prop_log_marg_condPost);


  }//end of for responses

  //cout << "After loop on responses" << endl;
  if (DEBUG_sample_g==true)
  {
    cout << "With single g, nb_chains = " << n_chains << endl;
    cout << "cum_diff = " << cum_diff << endl;

    cout << endl << endl;
  }


  if (Single_g==true)
  {
    // Begin Change from ESS
    logPG=getPriorG(PR_per_Resp[0],gSampleFlag,g[0]);
    logPG_Prop=getPriorG(PR_per_Resp[0],gSampleFlag,g_prop);

    /*
    cout << "g = " << g[0] << endl;
    cout << "logP(g) = " << logPG << endl;
    cout << "g_prop = " << g_prop << endl;
    cout << "logP(g_prop) = " << logPG_Prop << endl;
    *


    cum_diff=cum_diff+logPG_Prop-logPG;
    // End Change from ESS

    cum_diff=cum_diff+log(g_prop)-log(g[0]);

    /*
    cout << "Log Ratio = " << cum_diff << endl;
    cout << endl << endl;
    *

    alpha_g=min(1.0,exp(cum_diff));

    if(DEBUG_HESS){cout << endl;cout << "alpha_g= " << alpha_g << endl;}

    if(!isnan(alpha_g))
    {
        double rand_test = myrand(RandomNumberGenerator);
        if(rand_test<alpha_g && g_prop>0.0)
        {
          g[0]=g_prop;
          for(unsigned int k=0;k<n_Responses;k++)
          {
              for(unsigned int chain=0;chain<n_chains;chain++)
              {
                unsigned int pos_chain=Gammas[k].chain_idx[chain];
                Gammas[k].mat_log_cond_post.matrix[pos_chain][sweep]=store_log_cond_post[k][chain];
                Gammas[k].mat_log_marg.matrix[pos_chain][sweep]=store_log_marg[k][chain];
              }
          }

          (*My_g_AdMH[0]).tilda_accept++;
          (*My_g_AdMH[0]).tilda_accept_ins++;
        }
    }
    else{
        if(DEBUG_HESS){cout << "alpha_g is NaN" << endl;}
    }
    //  (*My_Move_monitor).g_sample_history.push_back(g);

    (*My_g_AdMH[0]).tilda_n_sweep++;
    (*My_g_AdMH[0]).tilda_n_sweep_ins++;

  }

  store_log_marg.clear();
  store_log_cond_post.clear();
//gsl_vector_free(prop_log_marg_condPost);

  ////////////////////////////////////////////////////
  //  Adaptation: every My_g_AdMH.n_batch sweeps
  ////////////////////////////////////////////////////

  for (unsigned int k=0;k<My_g_AdMH.size();k++)
  {

    if(sweep%(*My_g_AdMH[k]).n_batch==0)
    {
      if(DEBUG_HESS){
        cout << "G-adaptation" << endl;
        cout << "sweep=" << sweep
         << " -- My_g_AdMH.n_batch " << (*My_g_AdMH[k]).n_batch
         << " -- Test " << sweep%(*My_g_AdMH[k]).n_batch << endl;
      }
      double local_accept_rate=(double)((*My_g_AdMH[k]).tilda_accept_ins)/(double)((*My_g_AdMH[k]).tilda_n_sweep_ins);

      double batchNumber = (double)sweep/(double)(*My_g_AdMH[k]).n_batch;
      double delta = (pow(batchNumber,-0.5)<(*My_g_AdMH[k]).delta_n ? pow(batchNumber,-0.5) : (*My_g_AdMH[k]).delta_n );
      if(local_accept_rate<(*My_g_AdMH[k]).optimal){
        //Updating AdMH.ls
        (*My_g_AdMH[k]).ls-=delta;
        if((*My_g_AdMH[k]).ls < (*My_g_AdMH[k]).M[0])
        {
          (*My_g_AdMH[k]).ls = (*My_g_AdMH[k]).M[0];
        }
      }
      else{
        (*My_g_AdMH[k]).ls+=delta;
        if((*My_g_AdMH[k]).ls > (*My_g_AdMH[k]).M[1])
        {
          (*My_g_AdMH[k]).ls = (*My_g_AdMH[k]).M[1];
        }
      }
      (*My_g_AdMH[k]).tilda_accept_ins=0;
      (*My_g_AdMH[k]).tilda_n_sweep_ins=0;
      if(DEBUG_HESS){
        cout << "Sweep " << sweep << " -- Update the AdMH parameters" << endl;
        cout << "Updated g_ls " << (*My_g_AdMH[k]).ls << endl;
      }

    }
  }


}
*/


void sample_g_ESS(unsigned int Response,
          Double_Matrices mat_log_marg,
	      Double_Matrices mat_log_cond_post,
	      AdMH *My_g_AdMH,
	      Temperatures *t_tun,
          vector < vector <unsigned int> > &vect_gam,
          gsl_matrix *mat_X,
	      gsl_matrix *mat_Y,
	      bool gPriorFlag,
	      bool indepPriorFlag,
	      bool gSampleFlag,
	      double lambda,
	      double &g,
	      Prior_param PR,
	      unsigned int sweep,
	      Move_monitor *My_Move_monitor,
	      vector <unsigned int > &chain_idx,
	      vector <unsigned int > &n_Models_visited,
              bool cudaFlag,
              gsl_rng *RandomNumberGenerator,
              Model_Information &Model_Data)
{

  if(DEBUG){
    cout << "//////////////////////////////////////////" << endl
	 << "//  Adaptive M-H: sampling g" << endl
	 << "//////////////////////////////////////////" << endl << endl;
  }
  double norm_mean=log(g);
  double norm_sd=exp((*My_g_AdMH).ls);
  
  double sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);
  double g_prop=exp(sample_tmp);
  
  unsigned int n_chains=vect_gam.size();
  gsl_vector *prop_log_marg_condPost=gsl_vector_calloc(2);
  //vector < unsigned int > list_columns_X_gam;
  vector < double > store_log_marg;
  vector < double > store_log_cond_post;
  
  if(DEBUG){
    cout << "norm_mean " << norm_mean
	 << " -- norm sd " << norm_sd
	 << " -- sample tmp " << sample_tmp
	 << " -- g_prop " << g_prop
	 << endl;
  }
  double cum_diff=0.0;

  for(unsigned int current_chain=0;current_chain<n_chains;current_chain++){
    unsigned int pos_current_chain=chain_idx[current_chain];
    if(DEBUG){
      cout << "CHAIN " << current_chain 
	   << " -- pos " << pos_current_chain << endl;
    }
    double current_log_cond_post=mat_log_cond_post.matrix[pos_current_chain][sweep];
    double current_log_marg=mat_log_marg.matrix[pos_current_chain][sweep];
    
    //Step 1: Getting X_gamma  

	
    //Step 2: Calculate log_cond_post if move
    double propLogLik,propLogPost;
    if(DEBUG)
    {
        cout << current_chain << " ";
    }
    computeLogPosterior(vect_gam[pos_current_chain],current_chain,Response,
                        mat_X,propLogLik,propLogPost,mat_Y,PR,gPriorFlag,
                        indepPriorFlag,gSampleFlag,lambda,g_prop,
                        cudaFlag,
                        Model_Data);
    prop_log_marg_condPost->data[0]=propLogLik;
    prop_log_marg_condPost->data[1]=propLogPost;

    n_Models_visited[sweep]++;
    

    // BIG DISCUSSION HERE:

    //Current ESS
    cum_diff+=(prop_log_marg_condPost->data[1]-current_log_cond_post)/(*t_tun).t[current_chain];
    //Proposed change
    //cum_diff+=(prop_log_marg_condPost->data[0]-current_log_marg)/(*t_tun).t[current_chain];
    store_log_marg.push_back(prop_log_marg_condPost->data[0]);
    store_log_cond_post.push_back(prop_log_marg_condPost->data[1]);
    if(DEBUG){
      cout << "\tCurrent log_marg' " << current_log_marg
	   << " -- cond_post " << current_log_cond_post << endl;
      cout << "\tProposed log_marg' " << prop_log_marg_condPost->data[0]
	   << " -- cond_post " << prop_log_marg_condPost->data[1] << endl;
      cout << "\tcum_diff " << cum_diff << endl;
    }
    
  }//end of for chain

  (*My_Move_monitor).g_sample_nb_sweep++;

// DISCUSSION :
// Begin Change
//  double logPG=getPriorG(PR,gSampleFlag,g);
//  double logPG_Prop=getPriorG(PR,gSampleFlag,g_prop);
//  cum_diff=cum_diff+logPG_Prop-logPG;
// End Change
  double alpha_g=min(1.0,exp(cum_diff+log(g_prop)-log(g)));
  if(DEBUG){
    cout << "\talpha_g= " << alpha_g << endl;
  }
  if(!isnan(alpha_g)){
    double rand_test = myrand(RandomNumberGenerator);
    if(rand_test<alpha_g && g_prop>0.0){
      (*My_Move_monitor).g_sample_nb_accept++;
      g=g_prop;
      for(unsigned int chain=0;chain<n_chains;chain++){
	unsigned int pos_chain=chain_idx[chain];
	mat_log_cond_post.matrix[pos_chain][sweep]=store_log_cond_post[chain];
	mat_log_marg.matrix[pos_chain][sweep]=store_log_marg[chain];
      }
    
      (*My_g_AdMH).tilda_accept++;
      (*My_g_AdMH).tilda_accept_ins++;
    }
    
  }
  else{
    if(DEBUG){
      cout << "alpha_g is NaN" << endl;
    }
  }
  (*My_Move_monitor).g_sample_history.push_back(g);
  
  (*My_g_AdMH).tilda_n_sweep++;
  (*My_g_AdMH).tilda_n_sweep_ins++;
  
  if(DEBUG){
    cout << "**********************************************************************" << endl;
    cout << "End of g_sample " << endl;
    cout << "new g= " << g << endl;
    cout << endl << "log_marg:" << endl << "\t";
    for(unsigned int chain=0;chain<n_chains;chain++){
      unsigned int pos_chain=chain_idx[chain];
      cout << mat_log_marg.matrix[pos_chain][sweep] << " ";
    }
    cout << endl;
    cout << "log_cond_post: " << endl << "\t";
    for(unsigned int chain=0;chain<n_chains;chain++){
      unsigned int pos_chain=chain_idx[chain];
      cout << mat_log_cond_post.matrix[pos_chain][sweep] << " ";
    }
    cout << endl;
    cout << "**********************************************************************" << endl << endl;
  }
  store_log_marg.clear();
  store_log_cond_post.clear();
  gsl_vector_free(prop_log_marg_condPost);
  
  ////////////////////////////////////////////////////
  //  Adaptation: every My_g_AdMH.n_batch sweeps
  ////////////////////////////////////////////////////
  
  if(sweep%(*My_g_AdMH).n_batch==0){

      (*My_g_AdMH).Adapt_AdMH_ESS(sweep,My_Move_monitor);
    
    (*My_Move_monitor).g_adapt_history[3].push_back((*My_g_AdMH).ls);
    
  }
  

}

void FSMH_move(unsigned int Response,
           Double_Matrices mat_log_marg,
	       Double_Matrices mat_log_cond_post,
	       gsl_matrix *mat_X,
	       gsl_matrix *mat_Y,
	       Temperatures *t_tun,
	       gsl_permutation *MyPerm,
	       vector < vector <unsigned int> > &vect_gam,
	       unsigned int sweep,
	       bool gPriorFlag,
	       bool indepPriorFlag,
	       bool gSampleFlag,
	       double lambda,
	       double g,
	       Prior_param PR,
	       Move_monitor *My_Move_monitor,
	       vector <unsigned int > &chain_idx,
	       vector <unsigned int > &n_Models_visited,
               bool cudaFlag,
               unsigned int nConfounders,
               unsigned int maxPX,
               gsl_rng *RandomNumberGenerator,
               Model_Information &Model_Data)
{

    double omega_k;
    vector<double> * Rho_V;

  unsigned int n_chains=vect_gam.size();
  unsigned int pX=mat_X->size2;
  gsl_vector *prop_log_marg_condPost=gsl_vector_calloc(2);
  vector < unsigned int > list_columns_X_gam;
  double alpha=PR.w[1]*PR.w[0];
  double beta=PR.w[1]*(1.0-PR.w[0]);

  unsigned int count_nb_model=0;
  unsigned int count_nb_accept=0;

  unsigned int count_nb_model_01=0;
  unsigned int count_nb_accept_01=0;

  unsigned int count_nb_model_10=0;
  unsigned int count_nb_accept_10=0;

  if(DEBUG){
    cout << "alpha=" << alpha << " -- beta=" << beta << endl;
  }
  for(unsigned int current_chain=0;current_chain<n_chains;current_chain++)
  {

      if (Model_Data.Model_Tag=="HESS")
      {
          omega_k=(*(Model_Data.Omega_PerLine))[current_chain][Response];
          Rho_V=&(*(Model_Data.Rho_PerCol))[current_chain];
      }


    unsigned int pos_current_chain=chain_idx[current_chain];
    double tVal = (1.0/(*t_tun).t[current_chain]);
    if(DEBUG){
      cout << "!!!!!!!!!!!!" << endl;
      cout << "CHAIN " << current_chain 
	   << " -- pos " << pos_current_chain << endl;
      cout << "!!!!!!!!!!!!" << endl;
    }
    //Step1: get the random order of variable: shuffle the permutation object
    My_Permut_unsigned_int(MyPerm,RandomNumberGenerator);
    double current_log_cond_post=mat_log_cond_post.matrix[pos_current_chain][sweep];
    double current_log_marg=mat_log_marg.matrix[pos_current_chain][sweep];
    unsigned int sum_gam=sum_line_std_mat(vect_gam,
					  pos_current_chain);

    for(unsigned int current_variable=0;current_variable<MyPerm->size;current_variable++){
      unsigned int pos_curr_var=MyPerm->data[current_variable];
      if(pos_curr_var<nConfounders){
        continue;
      }
      unsigned int current_value=vect_gam[pos_current_chain][pos_curr_var];
      double rand_test=myrand(RandomNumberGenerator);
      double alpha_FSMH=0.0;
      bool doMove =false;
 

      if(DEBUG){
	cout << "Variable rank " << current_variable+1
	     << " -- position " << pos_curr_var
	     << " -- current_value " << current_value
	     << " -- rand_test " << rand_test << endl;
      }
      
      //Step 2: Calculate alpha_FSMH: probability to make a move
      bool autoReject = false;
      if(current_value==0){
        // Note in the next line, because p_gamma (sum gam) is the value with the
        // current covariate set to 0, we need to add 1 to the numerator
          double theta_1=0;
          if (Model_Data.Model_Tag=="ESS")
          {
            theta_1=max(0.0, (double)(sum_gam+alpha) / (double)(pX+alpha+beta-1));
          }
          else if (Model_Data.Model_Tag=="HESS")
          {
            theta_1=max(0.0,omega_k * (*Rho_V)[current_variable]);
          }

        double theta_0=1.0-theta_1;
        //double theta_1_tilda=pow(theta_1,tVal) / (pow(theta_0,tVal)+pow(theta_1,tVal));
        //double theta_0_tilda=1.0-theta_1_tilda;
        double theta_1_tilda=theta_1;
        double theta_0_tilda=theta_0;
	//Step 2.1: in the variable is not in p(alpha_FSMH=0.0)=theta_1_tilda
	if(DEBUG){
	  cout << "\ttheta_1_tilda=" << theta_1_tilda << endl;
	}
	if(rand_test<theta_1_tilda){
	  doMove = true;
	  if(DEBUG){
	    cout << "\tTry to change gamma 0->1";
	  }
	  count_nb_model_01++;
	  
	  //double log_L=current_log_marg*(1.0/(*t_tun).t[current_chain]);
	  double log_P=current_log_cond_post*(1.0/(*t_tun).t[current_chain]);
	  //Step 2.2: if alpha_FSMH!=0.0, alpha_FSMH=f(log_prob,theta_1...)
	  vect_gam[pos_current_chain][pos_curr_var]=1-vect_gam[pos_current_chain][pos_curr_var];
	  if(DEBUG){
	    cout << " -- modified value " << vect_gam[pos_current_chain][pos_curr_var] << endl;
	  }
	  //Getting X_gamma  
      get_list_var_in(list_columns_X_gam,
              vect_gam[pos_current_chain]);

      unsigned int n_vars_in=list_columns_X_gam.size();
	  if(n_vars_in>maxPX){
	    autoReject=true;
	    alpha_FSMH=0.0;
	  }else{

	    //Calculate log_cond_post for the proposed move
	    double propLogLik,propLogPost;
        computeLogPosterior(vect_gam[pos_current_chain],current_chain,Response,
                            mat_X,propLogLik,propLogPost,mat_Y,PR,
                                gPriorFlag,indepPriorFlag,gSampleFlag,lambda,g,
                                //pX,nX,pY,
                                cudaFlag,
                                Model_Data);
	    prop_log_marg_condPost->data[0]=propLogLik;
	    prop_log_marg_condPost->data[1]=propLogPost;

        n_Models_visited[sweep]++;

        if(DEBUG){
          cout << "\tproposed log_marg " << prop_log_marg_condPost->data[0] << endl
              << "\tproposed log_cond_post " << prop_log_marg_condPost->data[1] << endl;
        }
        double temp=log(theta_0_tilda)-log(theta_1_tilda)+prop_log_marg_condPost->data[1]*tVal-log_P;
        alpha_FSMH=min(1.0, exp(temp));
	  }
	}
	else{
	  if(DEBUG){
	    cout << "\tNothing to do" << endl;
	  }
	}
      }
      else{//current_value==1
        //double theta_1=max(0.0, (double)(sum_gam+alpha-1) / (double)(pX+alpha+beta-1));
        double theta_1=0;
        if (Model_Data.Model_Tag=="ESS")
        {
          theta_1=max(0.0, (double)(sum_gam+alpha-1) / (double)(pX+alpha+beta-1));
        }
        else if (Model_Data.Model_Tag=="HESS")
        {
          theta_1=max(0.0,omega_k * (*Rho_V)[current_variable]);
        }

        double theta_0=1.0-theta_1;
        //double theta_1_tilda=pow(theta_1,tVal) / (pow(theta_0,tVal)+pow(theta_1,tVal));
        //double theta_0_tilda=1.0-theta_1_tilda;
        double theta_1_tilda=theta_1;
        double theta_0_tilda=theta_0;
	//Step 2.1: in the variable is not in p(alpha_FSMH=0.0)=theta_0_tilda
	if(DEBUG){
	  cout << "\ttheta_0_tilda=" << theta_0_tilda << endl;
	}
	if(rand_test<theta_0_tilda){
	  doMove = true;
	  if(DEBUG){
	    cout << "\tTry to change gamma 1->0";
	  }
	  count_nb_model_10++;
	  
	  //double log_L=current_log_marg*(1.0/(*t_tun).t[current_chain]);
          double log_P=current_log_cond_post*(1.0/(*t_tun).t[current_chain]);
	  //Step 2.2: if alpha_FSMH!=0.0, alpha_FSMH=f(log_prob,theta_1...)
	  vect_gam[pos_current_chain][pos_curr_var]=1-vect_gam[pos_current_chain][pos_curr_var];
	  if(DEBUG){
	    cout << " -- modified value " << vect_gam[pos_current_chain][pos_curr_var] << endl;
	  }


	  //Calculate log_cond_post for the proposed move
	  double propLogLik,propLogPost;

      computeLogPosterior(vect_gam[pos_current_chain],current_chain,Response,
                          mat_X,propLogLik,propLogPost,mat_Y,PR,
                                gPriorFlag,indepPriorFlag,gSampleFlag,lambda,g,
                                cudaFlag,
                                Model_Data);
      prop_log_marg_condPost->data[0]=propLogLik;
	  prop_log_marg_condPost->data[1]=propLogPost;

	  n_Models_visited[sweep]++;

	  if(DEBUG){
	    cout << "\tproposed log_marg " << prop_log_marg_condPost->data[0] << endl
		 << "\tproposed log_cond_post " << prop_log_marg_condPost->data[1] << endl;
	  }
	  double temp=log(theta_1_tilda)-log(theta_0_tilda)+prop_log_marg_condPost->data[1]*tVal-log_P;
	  alpha_FSMH=min(1.0, exp(temp));
	}
	else{
	  if(DEBUG){
	    cout << "\tNothing to do" << endl;
	  }
	}

      }//end of current variable==1
      if(DEBUG){
	cout << "\talpha_FSMH=" << alpha_FSMH 
	     << " -- Test " << (alpha_FSMH>0.0) << endl;
      }
      if(doMove){
        rand_test = myrand(RandomNumberGenerator);
        if(DEBUG){
          cout << "\tu=" << rand_test << endl;
        }
        if(rand_test >= alpha_FSMH||autoReject){
	  if(DEBUG){
	    cout << "\tMove Rejected" << endl;
	  }
	  vect_gam[pos_current_chain][pos_curr_var]=1-vect_gam[pos_current_chain][pos_curr_var];
	}
	else{
	  if(DEBUG){
	    cout << "\tMove Accepted" << endl;
	  }
	  if(vect_gam[pos_current_chain][pos_curr_var]==1){//Move 0->1 accepted
	    sum_gam++;
	    count_nb_accept_01++;
	  }
	  else{//Move 1->0 accepted
	    count_nb_accept_10++;
	    sum_gam--;
	  }
	  current_log_cond_post=prop_log_marg_condPost->data[1];
	  current_log_marg=prop_log_marg_condPost->data[0];
	}
      }
      if(DEBUG){
        cout << "\tvect_gam" << endl <<"\t";
        for(unsigned int i=0;i<vect_gam[0].size();i++){
            if(vect_gam[pos_current_chain][i]==1){
            cout << i+1 << " ";
	  }
	}
	cout << endl;
	cout << "\tlogmarg " << current_log_marg << endl
	     << "\tlog_cond_post " << current_log_cond_post << endl;
      }
    }//end of for variable
    //cout << endl;//!!!!!!!!!!!!!!!!
    
    mat_log_cond_post.matrix[pos_current_chain][sweep]=current_log_cond_post;
    mat_log_marg.matrix[pos_current_chain][sweep]=current_log_marg;
    unsigned int final_n_vars_in=sum_line_std_mat(vect_gam,
						  pos_current_chain);
    if(DEBUG){
      cout << "**********************************************************************" << endl;
      cout << "End of the variables, Final Results for chain " << current_chain << endl;
      cout << "\tlog_cond_post " << mat_log_cond_post.matrix[pos_current_chain][sweep]
	   << " -- log_marg " << mat_log_marg.matrix[pos_current_chain][sweep]
	   << " -- # var in " << final_n_vars_in << endl;
      cout << "\t";
      
      for(unsigned int col=0;col<vect_gam[pos_current_chain].size();col++){
	if(vect_gam[pos_current_chain][col]==1){
	  cout << col+1 << " ";
	}
      }
      cout << endl;
      cout << "**********************************************************************" << endl << endl;
    }
  }//end of for chain
  
  count_nb_model=count_nb_model_01+count_nb_model_10; 
  count_nb_accept=count_nb_accept_01+count_nb_accept_10;
  //Storing summary in monitor
  (*My_Move_monitor).FSMH_nb_sweep++;
  (*My_Move_monitor).FSMH_nb_model_tot+=count_nb_model;
  (*My_Move_monitor).FSMH_nb_model_0_1+=count_nb_model_01;
  (*My_Move_monitor).FSMH_nb_model_1_0+=count_nb_model_10;
  (*My_Move_monitor).FSMH_nb_accept_tot+=count_nb_accept;
  (*My_Move_monitor).FSMH_nb_accept_0_1+=count_nb_accept_01;
  (*My_Move_monitor).FSMH_nb_accept_1_0+=count_nb_accept_10;
  
  (*My_Move_monitor).FSMH_list_sweep.push_back(sweep);
  (*My_Move_monitor).FSMH_nb_model_per_sweep.push_back(count_nb_model);
  (*My_Move_monitor).FSMH_nb_model_per_sweep_0_1.push_back(count_nb_model_01);
  (*My_Move_monitor).FSMH_nb_model_per_sweep_1_0.push_back(count_nb_model_10);
  (*My_Move_monitor).FSMH_nb_accept_per_sweep.push_back(count_nb_accept);
  (*My_Move_monitor).FSMH_nb_accept_per_sweep_0_1.push_back(count_nb_accept_01);
  (*My_Move_monitor).FSMH_nb_accept_per_sweep_1_0.push_back(count_nb_accept_10);

  gsl_vector_free(prop_log_marg_condPost);
}

void temp_placement(Temperatures *t_tun,
		    DR *My_DR,
		    Move_monitor *My_Move_monitor,
		    unsigned int sweep,
		    unsigned int n_vars_in_last_chain,
		    unsigned int nX,
		    bool iso_T_Flag)
{

  unsigned int n_chains=((*t_tun).t).size();
  unsigned int DR_n_accept_adj=0;
  unsigned int DR_n_proposed_adj=0;
  double DR_acceptance_rate;
  unsigned int firstChainAccept;

  // Use moves between all adjacent chains
  for(unsigned int row=0;row<(*My_DR).mat_moves_accepted.size()-1;row++){
    DR_n_accept_adj+=(*My_DR).mat_moves_accepted[row][row+1];
    DR_n_proposed_adj+=(*My_DR).mat_moves_proposed[row][row+1];
  }

  firstChainAccept=(*My_DR).mat_moves_accepted[0][1];
  if(DR_n_proposed_adj>0){
    DR_acceptance_rate=(double)(DR_n_accept_adj)/(double)(DR_n_proposed_adj);
  }else{
    DR_acceptance_rate=0.0;
  }



  if(DEBUG){
    cout << "DR_n_accept " << DR_n_accept_adj
       << " -- DR_n_try " << DR_n_proposed_adj
       << " -- firstChainAccept" << firstChainAccept
       << " -- DR acceptance rate " << DR_acceptance_rate << endl
       << "n_vars_in_last_chain " << n_vars_in_last_chain
       << " -- 10*nX " << 10*nX << endl;
  }
  //Saving the temperatures values
  (*My_Move_monitor).temperature_history[0].push_back(sweep);
  for(unsigned int chain=0;chain<n_chains;chain++){
    (*My_Move_monitor).temperature_history[chain+1].push_back((*t_tun).t[chain]);
  }
  //Adjusting temperatures


  double temp_new_b_t=(*t_tun).b_t;
  double new_b_t=(*t_tun).b_t;

  if(DR_n_accept_adj==0||firstChainAccept==0 || n_vars_in_last_chain> 10*nX ){
    temp_new_b_t=(*t_tun).b_t - ((*t_tun).b_t-1.0)/2.0;
    new_b_t=max((*t_tun).M[0],temp_new_b_t);
    if(DEBUG){
      cout << "No DR move accepted, or too many vars in last chain:" << endl
	   << "\tAcceptance rate= " << DR_acceptance_rate
	   << " -- n_vars_in_last_chain= " << n_vars_in_last_chain << endl;
      
      cout << "\tChanging b_t= " << (*t_tun).b_t
	   << " -- temp new b_t " << temp_new_b_t
	   << " -- to new_b_t = " << new_b_t
	   << endl;
    }
  }
  else if(DR_n_accept_adj==DR_n_proposed_adj){
    temp_new_b_t=(*t_tun).b_t + ((*t_tun).b_t-1.0)/2.0;
    new_b_t=min((*t_tun).M[1],temp_new_b_t);
    if(DEBUG){
      cout << "All DR move accepted:" << endl
	   << "\t-Acceptance rate= " << DR_acceptance_rate
	   << " -- n_vars_in_last_chain= " << n_vars_in_last_chain << endl;
      
      cout << "\tChanging b_t= " << (*t_tun).b_t
	   << " -- temp new b_t " << temp_new_b_t
	   << " -- to new_b_t = " << new_b_t
	   << endl;
    }

  }
  else if(DR_acceptance_rate<(*t_tun).optimal){
    double argument=log((*t_tun).b_t)/log(2.0)-(*t_tun).delta_n;
    temp_new_b_t=pow(2.0,argument);
    new_b_t=max((*t_tun).M[0],temp_new_b_t);

    if(DEBUG){
      cout << "Acceptance >0 and <1 and <optimal:" << endl
	   << "\t-Acceptance rate= " << DR_acceptance_rate
	   << " -- n_vars_in_last_chain= " << n_vars_in_last_chain 
	   << " -- optimal " << (*t_tun).optimal
	   << endl;
      
      cout << "\tChanging b_t= " << (*t_tun).b_t
	   << " -- argument " << argument
	   << " -- temp new b_t " << temp_new_b_t
	   << " -- to new_b_t = " << new_b_t
	   << endl;
    }
      

  }
  else if(DR_acceptance_rate>(*t_tun).optimal){
    double argument=log((*t_tun).b_t)/log(2.0)+(*t_tun).delta_n;
    temp_new_b_t=pow(2.0,argument);
    new_b_t=min((*t_tun).M[1],temp_new_b_t);
  
    if(DEBUG){
      cout << "Accepetance >0 and <1 and >optimal:" << endl
	   << "\t-Acceptance rate= " << DR_acceptance_rate
	   << " -- n_vars_in_last_chain= " << n_vars_in_last_chain 
	   << " -- optimal " << (*t_tun).optimal
	   << endl;
      
      cout << "\tChanging b_t= " << (*t_tun).b_t
	   << " -- argument " << argument
	   << " -- temp new b_t " << temp_new_b_t
	   << " -- to new_b_t = " << new_b_t
	   << endl;
    }
  }
  (*t_tun).b_t=new_b_t;
  for(unsigned int chain=0;chain<n_chains;chain++){
    double exponent=(*t_tun).a_t[chain];
    if(!iso_T_Flag){
      (*t_tun).t[chain]=pow(new_b_t,exponent);
    }
    else{
      (*t_tun).t[chain]=1.0;
    }
    if(DEBUG){
      cout << "Chain " << chain
	   << "-- exponent " << exponent
	   << " -- temp " << (*t_tun).t[chain]
	   << endl;
    }
  }
  (*My_DR).nb_calls_adj=0;
  (*My_DR).nb_calls=0;
  for(unsigned int row=0;row<(*My_DR).mat_moves_accepted.size();row++){
    for(unsigned int col=0;col<(*My_DR).mat_moves_accepted[row].size();col++){
    (*My_DR).mat_moves_accepted[row][col]=0;
    (*My_DR).mat_moves_proposed[row][col]=0;
    }
  }
}

void Crossover_move(unsigned int Response,
            Double_Matrices mat_log_marg,
		    Double_Matrices mat_log_cond_post,
		    gsl_matrix *mat_X,
		    gsl_matrix *mat_Y,
		    Temperatures *t_tun,
		    vector < vector <unsigned int> > &vect_gam,
		    bool gPriorFlag,
		    bool indepPriorFlag,
		    bool gSampleFlag,
		    double lambda,
		    double g,
		    Prior_param PR,
		    Move_monitor *My_Move_monitor,
		    vector <unsigned int > &chain_idx,
		    unsigned int sweep,
		    CM *My_CM,
		    vector <unsigned int > &n_Models_visited,
                    bool cudaFlag,
                    unsigned int maxPX,
                    gsl_rng *RandomNumberGenerator,
                    Model_Information &Model_Data)
{


  unsigned int n_chains=chain_idx.size();
  unsigned int c_trunc=(unsigned int)(max(ceil((double)(n_chains)*PR.Prob_sel),2.0));
  vector < double > vect_pbty_Boltz;
  vector < unsigned int > sampled_chains;
  vect_pbty_Boltz.resize(n_chains);
  sampled_chains.resize(2);
  unsigned int n_vars_tot=mat_X->size2;
  (*My_Move_monitor).CM_nb_sweep++;
  (*My_Move_monitor).CM_nb_model+=2;

  unsigned int n_pts_recomb=0;
  unsigned int n_changes=0;

  double chosenT;
  //Step 1: Sampling the two chains
  //Step 1.1: Sampling the two chains and calculating the Boltz Pbty vector.
  
  computeAndSampleCondPostBoltz(vect_pbty_Boltz,
                                sampled_chains,
                                n_chains,
                                c_trunc,
                                sweep,
                                mat_log_cond_post,
                                chain_idx,
                                t_tun,
                                chosenT,
                                RandomNumberGenerator);
  
  
  if(DEBUG){
    cout << "OUT: log_post_cond_Boltz" << endl;
    for(unsigned int current_chain=0;current_chain<c_trunc;current_chain++){
      cout << vect_pbty_Boltz[current_chain] << " ";
    }
    cout << endl;
    cout << "c1= " << sampled_chains[0] 
	 << " -- c2 = " << sampled_chains[1] << endl;
  }
  
  //Step 2: Sampling the CM type
  unsigned int type_CM_move_sampled=SampleFromDiscrete_new((*My_CM).unit_move_pbty_cum,RandomNumberGenerator);
  (*My_Move_monitor).CM_n_nb_sweep_per_CM_move[type_CM_move_sampled]++;
  vector < unsigned int > vect_gam_prop_c1;
  vector < unsigned int > vect_gam_prop_c2;
    
  vect_gam_prop_c1.resize(n_vars_tot);
  vect_gam_prop_c2.resize(n_vars_tot);

  
  //Step 3: Sampling the position of the breakpoints
  if(type_CM_move_sampled<(*My_CM).n_max_breakpoint){
    unsigned int n_bkpts=type_CM_move_sampled +1;
    n_pts_recomb=n_bkpts;
    if(DEBUG){
      cout << n_bkpts << "-point Crossover Move sampled" << endl;
    }
    double *pos_bkpts = new double[n_bkpts];
    double *full_vars = new double[n_vars_tot-1];
    //we do not sample bkpt at pos=0 or nvars
    for(unsigned current_var=0;current_var<n_vars_tot-1;current_var++){
      full_vars[current_var]=current_var;
    }
    
    My_gsl_ran_choose_double(pos_bkpts,n_bkpts,full_vars,n_vars_tot-1,RandomNumberGenerator);
    gsl_sort(pos_bkpts,1,n_bkpts);
    for(unsigned int curr_bkpt=0;curr_bkpt<n_bkpts;curr_bkpt++){
      pos_bkpts[curr_bkpt]+=1;
    }
    if(DEBUG){
      cout << "Sampled Breakpoints" << endl;
      for(unsigned int bkpt=0;bkpt<n_bkpts;bkpt++){
	cout << (unsigned int)(pos_bkpts[bkpt]) << " ";
      }
      cout << endl;
    }
    delete [] full_vars;

    //Step 4: Generating the two candidate chains
    if(DEBUG){
      cout << "c1 init" << endl;
      for(unsigned int var=0;var<n_vars_tot;var++){
	if(vect_gam[chain_idx[sampled_chains[0]]][var]==1){
	  cout << var << " ";
	}
      }
      cout << endl;
      cout << "c2 init" << endl;
      for(unsigned int var=0;var<n_vars_tot;var++){
	if(vect_gam[chain_idx[sampled_chains[1]]][var]==1){
	  cout << var << " ";
	}
      }
      cout << endl;
    }

    n_changes=recombine_chains(vect_gam_prop_c1,
			       vect_gam_prop_c2,
			       sampled_chains,
			       chain_idx,
			       vect_gam,
			       pos_bkpts,
			       n_bkpts);
    if(DEBUG){
      cout << "c1 recomb" << endl;
      for(unsigned int col=0;col<n_vars_tot;col++){
	if(vect_gam_prop_c1[col]==1){
	  cout << col << " "; 
	}
      }
      cout << endl;
      cout << "c_2 recomb" << endl;
      for(unsigned int col=0;col<n_vars_tot;col++){
	if(vect_gam_prop_c2[col]==1){
	  cout << col << " ";
	}
      }
      cout << endl;    
    }
    delete [] pos_bkpts;
    
  }
  else{
    if(DEBUG){
      cout << "Haplotype Crossover Move sampled" << endl;
    }
    //Sampling the position of the crossover
    unsigned int pos_crsv=(unsigned int)(myrand(RandomNumberGenerator)*(double)(n_vars_tot));
    if(DEBUG){    
      cout << "Pos crsv " << pos_crsv << endl;
    }
    vector < unsigned int > r_idx;
    define_haplotype_bkpts(r_idx,
			   pos_crsv,
			   mat_X,
			   PR);
    n_pts_recomb=r_idx.size();
 
    if(DEBUG){
      cout << "c1 init" << endl;
      for(unsigned int var=0;var<n_vars_tot;var++){
	if(vect_gam[chain_idx[sampled_chains[0]]][var]==1){
	  cout << var << " ";
	}
      }
      cout << endl;
      cout << "c2 init" << endl;
      for(unsigned int var=0;var<n_vars_tot;var++){
	if(vect_gam[chain_idx[sampled_chains[1]]][var]==1){
	  cout << var << " ";
	}
      }
      cout << endl;
    }
    n_changes=recombine_haplotype(vect_gam_prop_c1,
				  vect_gam_prop_c2,
				  sampled_chains,
				  chain_idx,
				  vect_gam,
				  r_idx);
    
    if(DEBUG){
      cout << "c1 recomb" << endl;
      for(unsigned int col=0;col<n_vars_tot;col++){
	if(vect_gam_prop_c1[col]==1){
	  cout << col << " "; 
	}
      }
      
      cout << endl;
      cout << "c_2 recomb" << endl;
      for(unsigned int col=0;col<n_vars_tot;col++){
	if(vect_gam_prop_c2[col]){
	  cout << col << " ";
	}
      }

      cout << endl;    
    }
    r_idx.clear();
  }
  
  bool autoReject = false;
  unsigned int nVarsInPropC1=0, nVarsInPropC2=0;
  for(unsigned int j=0;j<n_vars_tot;j++){
    if(vect_gam_prop_c1[j]==1){
      nVarsInPropC1++;
    }
    if(vect_gam_prop_c2[j]==1){
      nVarsInPropC2++;
    }
  }

  // Make sure that the move doesn't make one or the chains have more than the
  // maximum allowed
  double alpha_CM=0.0;
  gsl_vector *prop_log_marg_condPost_c1=gsl_vector_calloc(2);
  gsl_vector *prop_log_marg_condPost_c2=gsl_vector_calloc(2);
  unsigned int c1=sampled_chains[0];
  unsigned int c2=sampled_chains[1];
  unsigned int pos_c1=chain_idx[c1];
  unsigned int pos_c2=chain_idx[c2];
  gsl_vector *vect_log_cond_post_init=gsl_vector_calloc(n_chains);
  gsl_vector *vect_log_cond_post_prop=gsl_vector_calloc(n_chains);


  if(nVarsInPropC1>maxPX||nVarsInPropC2>maxPX){
    autoReject = true;
  }else{

    //Step 5: Calculating logCondPost/logMarg for both candidate chains
  
  
    get_log_cond_post_log_marg_prop(c1,Response,
                                    prop_log_marg_condPost_c1,
                                       vect_gam_prop_c1,
                                       mat_X,
                                       mat_Y,
                                       PR,
                                       gPriorFlag,
                                       indepPriorFlag,
                                       gSampleFlag,
                                       lambda,
                                       g,
                                       cudaFlag,
                                       Model_Data);

    get_log_cond_post_log_marg_prop(c2,Response,
                                    prop_log_marg_condPost_c2,
                                       vect_gam_prop_c2,
                                       mat_X,
                                       mat_Y,
                                       PR,
                                       gPriorFlag,
                                       indepPriorFlag,
                                       gSampleFlag,
                                       lambda,
                                       g,
                                       cudaFlag,
                                       Model_Data);

    if(n_changes>0){
      n_Models_visited[sweep]+=2;
    }

    for(unsigned int chain=0;chain<n_chains;chain++){
      unsigned int pos_chain=chain_idx[chain];
      vect_log_cond_post_init->data[chain]=mat_log_cond_post.matrix[pos_chain][sweep];
      if(chain==sampled_chains[0]){
        vect_log_cond_post_prop->data[chain]=prop_log_marg_condPost_c1->data[1];
      }
      else if(chain==sampled_chains[1]){
        vect_log_cond_post_prop->data[chain]=prop_log_marg_condPost_c2->data[1];
      }
      else{
        vect_log_cond_post_prop->data[chain]=vect_log_cond_post_init->data[chain];
      }
    }

    //Step 6: Calculating the acceptance pbty

    double log_T_x_y=log(vect_pbty_Boltz[c1])+log(vect_pbty_Boltz[c2])-log(1-vect_pbty_Boltz[c1]);
    double log_T_y_x=computeLogCondPostBoltz(vect_log_cond_post_prop,
                                                 chosenT,
                                                 c1,
                                                 c2,
                                                 c_trunc);


    if(DEBUG){
      cout << "log_T_x_y " << log_T_x_y
          << " -- log_T_y_x " << log_T_y_x
          << endl;
    }
    double R_CM=((vect_log_cond_post_prop->data[c1]- vect_log_cond_post_init->data[c1]) * (1 / (*t_tun).t[c1])+
                (vect_log_cond_post_prop->data[c2]- vect_log_cond_post_init->data[c2]) * (1 / (*t_tun).t[c2])+
                 log_T_y_x-log_T_x_y);
    if(DEBUG){
      cout << "R_CM " << R_CM << endl;
    }
    alpha_CM=min(1.0,exp(R_CM));
    if(DEBUG){
      cout << "alpha_CM " << alpha_CM << endl;
    }
  }
  double accept_test=myrand(RandomNumberGenerator);
  //Step 7: Updating X_gam, and log_condPost/log_marg
  if(accept_test< alpha_CM&&!autoReject){

    (*My_Move_monitor).CM_nb_accept++;
    (*My_Move_monitor).CM_n_nb_accept_per_CM_move[type_CM_move_sampled]++;
    (*My_Move_monitor).CM_history[0].push_back(sweep);
    (*My_Move_monitor).CM_history[1].push_back(type_CM_move_sampled);
    (*My_Move_monitor).CM_history[2].push_back(n_pts_recomb);
    (*My_Move_monitor).CM_history[3].push_back(c1);
    (*My_Move_monitor).CM_history[4].push_back(c2);
    if(DEBUG){
    cout << "Move accepted" << endl;
    }
    //Updating log_cond_post and logmarg
    mat_log_marg.matrix[pos_c1][sweep]=prop_log_marg_condPost_c1->data[0];
    mat_log_cond_post.matrix[pos_c1][sweep]=prop_log_marg_condPost_c1->data[1];

    mat_log_marg.matrix[pos_c2][sweep]=prop_log_marg_condPost_c2->data[0];
    mat_log_cond_post.matrix[pos_c2][sweep]=prop_log_marg_condPost_c2->data[1];
    //Updating the vect_gam
    for(unsigned int var=0;var<n_vars_tot;var++){
      vect_gam[pos_c1][var]=vect_gam_prop_c1[var];
      vect_gam[pos_c2][var]=vect_gam_prop_c2[var];
    }
  }
  else{
    if(DEBUG){
      cout << "Move Rejected" << endl;
    }
  }
  gsl_vector_free(prop_log_marg_condPost_c1);
  gsl_vector_free(prop_log_marg_condPost_c2);

  gsl_vector_free(vect_log_cond_post_init);
  gsl_vector_free(vect_log_cond_post_prop);

  vect_gam_prop_c1.clear();
  vect_gam_prop_c2.clear();
  sampled_chains.clear();
  vect_pbty_Boltz.clear();
}

void computeAndSampleCondPostBoltz(vector<double>& condPostBoltzProb,
                                      vector<unsigned int>& sampledChains,
                                      const unsigned int& nChains,
                                      const unsigned int& cTrunc,
                                      const unsigned int& sweep,
                                      const Double_Matrices& matLogCondPost,
                                      const vector<unsigned int>& chainIdx,
                                      Temperatures* const tTun,
                                      double& chosenT,
                                      gsl_rng *RandomNumberGenerator)
{


  // Step 1. Sort the log conditional posterior of the chains
  gsl_vector* logCondPostInit=gsl_vector_calloc(nChains);
  double* chainTemps = new double[nChains];
  for(unsigned int currChain=0;currChain<nChains;currChain++){
    unsigned int posChain = chainIdx[currChain];
    logCondPostInit->data[currChain]=matLogCondPost.matrix[posChain][sweep];
    chainTemps[currChain]=tTun->t[posChain];
  }
  gsl_permutation* logCondPostSortIdx=gsl_permutation_calloc(nChains);
  gsl_sort_vector_index(logCondPostSortIdx,logCondPostInit);
  gsl_vector* IdxVec=gsl_vector_calloc(nChains);
  gsl_permutation_reverse(logCondPostSortIdx);
  for(unsigned int i=0;i<nChains;i++){
    IdxVec->data[i]=logCondPostSortIdx->data[i];
  }
  gsl_permutation* logCondPostSortRank=gsl_permutation_calloc(nChains);
  gsl_sort_vector_index(logCondPostSortRank,IdxVec);

  // Step 2. Get the temperature for the move
  gsl_sort(chainTemps,1,nChains);
  chosenT=gsl_stats_median_from_sorted_data(chainTemps,1,nChains);

  // Step 3. Calculate F_t (the Boltzmann normalising const) and p_t the Bolztmann probs
  double maxLogCondPost = logCondPostInit->data[logCondPostSortIdx->data[0]];
  double Z=0.0;
  for(unsigned int i=0;i<nChains;i++){
    Z+=exp((logCondPostInit->data[i]-maxLogCondPost)/chosenT);
  }

  double normConst=0.0;
  for(unsigned int currChain=0;currChain<nChains;currChain++){
    unsigned int currRank = logCondPostSortRank->data[currChain];
    if(currRank<cTrunc){
      condPostBoltzProb[currChain]=1.0e-10+exp((logCondPostInit->data[currChain]-maxLogCondPost)/chosenT)/Z;
    }else{
      condPostBoltzProb[currChain]=1.0e-10;
    }
    normConst+=condPostBoltzProb[currChain];
  }
  for(unsigned int i=0;i<nChains;i++){
    condPostBoltzProb[i]/=normConst;
  }


  //Sampling first chain
  unsigned int sampledChain1=SampleFromDiscrete_non_cum(condPostBoltzProb,RandomNumberGenerator);
  
  // Now remove the sampled chain (set it's prob to zero and renormalise)
  vector<double> condPostBoltzProbNew(nChains);
  for(unsigned int i=0;i<nChains;i++){
    if(i==sampledChain1){
      condPostBoltzProbNew[i]=0.0;
    }else{
      condPostBoltzProbNew[i]=condPostBoltzProb[i]/(1.0-condPostBoltzProb[sampledChain1]);
    }
  }
  //Sampling second chain
  unsigned int sampledChain2=SampleFromDiscrete_non_cum(condPostBoltzProbNew,RandomNumberGenerator);

  if(DEBUG){  
    cout << " -- sampledChain#1 " << sampledChain1 << endl;
  }
  if(DEBUG){
    cout << " -- sampledChain#2 " << sampledChain2 << endl;
  }
  sampledChains[0]=sampledChain1;
  sampledChains[1]=sampledChain2;

  delete [] chainTemps;
  gsl_permutation_free(logCondPostSortRank);
  gsl_permutation_free(logCondPostSortIdx);
  gsl_vector_free(IdxVec);
  gsl_vector_free(logCondPostInit);
 
}

unsigned int recombine_chains(vector < unsigned int > &vect_gam_prop_c1,
			      vector < unsigned int > &vect_gam_prop_c2,
			      vector < unsigned int > &sampled_chains,
			      vector <unsigned int > &chain_idx,
			      vector < vector <unsigned int> > &vect_gam,
			      double *pos_bkpts,
			      unsigned int n_bkpts)
{
  unsigned int n_vars_tot=vect_gam[0].size();
  unsigned int c1=sampled_chains[0];
  unsigned int pos_c1=chain_idx[c1];
  unsigned int c2=sampled_chains[1];
  unsigned int pos_c2=chain_idx[c2];
  unsigned int lower_bound=0;
  unsigned int upper_bound=0;
  unsigned int nb_interval=n_bkpts+1;
  unsigned int how_many_changes=0;
 
  for(unsigned int interval=0;interval<nb_interval;interval++){
    if(interval==0){
      lower_bound=0;
      upper_bound=(unsigned int)(pos_bkpts[interval])-1;
    }
    else if(interval==nb_interval-1){
      lower_bound=(unsigned int)(pos_bkpts[interval-1]);
      upper_bound=n_vars_tot-1;
    }
    else{
      lower_bound=(unsigned int)(pos_bkpts[interval-1]);
      upper_bound=(unsigned int)(pos_bkpts[interval])-1;
    }
    if(DEBUG){
      cout << "Interval " << interval+1
	   << "/" << nb_interval
	   << " -- lower bound " << lower_bound
	   << " -- upper bound " << upper_bound
	   << endl;
    }
    if(interval%2!=0){
      for(unsigned int pos=lower_bound;pos<=upper_bound;pos++){
	vect_gam_prop_c1[pos]=vect_gam[pos_c2][pos];
	vect_gam_prop_c2[pos]=vect_gam[pos_c1][pos];
	if(vect_gam[pos_c2][pos]!=vect_gam[pos_c1][pos]){
	  how_many_changes++;
	}
      }
    }
    else{
      for(unsigned int pos=lower_bound;pos<=upper_bound;pos++){
	vect_gam_prop_c1[pos]=vect_gam[pos_c1][pos];
	vect_gam_prop_c2[pos]=vect_gam[pos_c2][pos];
      }
    }
  }
  return how_many_changes;
}

void define_haplotype_bkpts(vector < unsigned int > &r_idx,
			    unsigned int pos_crsv,
			    gsl_matrix *mat_X,
			    Prior_param PR)
{
  unsigned int n_vars_tot=mat_X->size2;
  unsigned int n_obs=mat_X->size1;
  //Getting the correlation matrix between X[:,crsv] and all columns of X
  vector < double > vect_r;
  vect_r.resize(n_vars_tot);
  unsigned int n_neg_r=0;
  
  for(unsigned int var=0;var<n_vars_tot;var++){
    double current_r=gsl_stats_correlation(&mat_X->data[pos_crsv],
					   n_vars_tot,
					   &mat_X->data[var],
					   n_vars_tot,
					   n_obs);
    if(current_r<=0.0){
      n_neg_r++;
    }
    vect_r[var]=current_r;
    
  }
  //Defining the sign of the correlation
  unsigned int selected_neg_r=0;
  if(DEBUG){
    cout << "selected_neg_r " << selected_neg_r << endl;
  }
  for(unsigned int var=0;var<n_vars_tot;var++){
    if(selected_neg_r==0){//postiive r
      if(vect_r[var]>=PR.Prob_crsv_r){
	r_idx.push_back(var);
      }
    }
    else{//negative r
      if(vect_r[var]<= -PR.Prob_crsv_r){
	r_idx.push_back(var);
      }
    }
  }
  if(DEBUG){
    cout << "r_idx, size " << r_idx.size() << endl;
    for(unsigned int col=0;col<r_idx.size();col++){
      cout << r_idx[col] << " ";
    }
    cout << endl;
  }
}

unsigned int recombine_haplotype(vector < unsigned int > &vect_gam_prop_c1,
				 vector < unsigned int > &vect_gam_prop_c2,
				 vector < unsigned int > &sampled_chains,
				 vector <unsigned int > &chain_idx,
				 vector < vector <unsigned int> > &vect_gam,
				 vector <unsigned int > &r_idx)
{
  unsigned int n_vars_tot=vect_gam[0].size();
  unsigned int c1=sampled_chains[0];
  unsigned int pos_c1=chain_idx[c1];
  unsigned int c2=sampled_chains[1];
  unsigned int pos_c2=chain_idx[c2];
  unsigned int rank_in_r_idx=0;
  unsigned int how_many_changes=0;
  for(unsigned int pos=0;pos<n_vars_tot;pos++){
    if(rank_in_r_idx<r_idx.size()){
      if(r_idx[rank_in_r_idx]==pos){
	vect_gam_prop_c1[pos]=vect_gam[pos_c2][pos];
	vect_gam_prop_c2[pos]=vect_gam[pos_c1][pos];
	rank_in_r_idx++;
	if(vect_gam[pos_c2][pos]!=vect_gam[pos_c1][pos]){
	  how_many_changes++;
	}
      }
      else{
	vect_gam_prop_c1[pos]=vect_gam[pos_c1][pos];
	vect_gam_prop_c2[pos]=vect_gam[pos_c2][pos];
      }
    }
    else{
	vect_gam_prop_c1[pos]=vect_gam[pos_c1][pos];
	vect_gam_prop_c2[pos]=vect_gam[pos_c2][pos];
    }
  }
  return(how_many_changes);
}

void get_log_cond_post_log_marg_prop(unsigned int Chain,unsigned int Response,
                     gsl_vector *prop_log_marg_condPost,
				     vector < unsigned int > vect_gam_prop,
				     gsl_matrix *mat_X,
				     gsl_matrix *mat_Y,
				     Prior_param PR,
				     bool gPriorFlag,
				     bool indepPriorFlag,
				     bool gSampleFlag,
				     double lambda,
				     double g,
                     bool cudaFlag,
                     Model_Information &Model_Data)
{


  //Step 2: Calculate log_cond_post if move
  double propLogLik,propLogPost;
//  computeLogPosterior(vect_gam_prop,propLogLik,propLogPost,mat_X_gam,mat_Y,PR,
    computeLogPosterior(vect_gam_prop,Chain,Response,
                        mat_X,propLogLik,propLogPost,mat_Y,PR,
                        gPriorFlag,indepPriorFlag,gSampleFlag,lambda,g,
                        //pX,nX,pY,
                        cudaFlag,
                        Model_Data);

  prop_log_marg_condPost->data[0]=propLogLik;
  prop_log_marg_condPost->data[1]=propLogPost;

}

double computeLogCondPostBoltz(gsl_vector* const logCondPostVec,
		     const double& chosenT,
		     const unsigned int& c1,
		     const unsigned int& c2,
		     const unsigned int& cTrunc)
{
  // Step 1. Rank the chains according to their logCondPost
  unsigned int nChains=logCondPostVec->size;

  gsl_permutation* logCondPostSortIdx=gsl_permutation_calloc(nChains);
  gsl_sort_vector_index(logCondPostSortIdx,logCondPostVec);
  gsl_vector* IdxVec=gsl_vector_calloc(nChains);
  gsl_permutation_reverse(logCondPostSortIdx);
  for(unsigned int i=0;i<nChains;i++){
    IdxVec->data[i]=logCondPostSortIdx->data[i];
  }
  gsl_permutation* logCondPostSortRank=gsl_permutation_calloc(nChains);
  gsl_sort_vector_index(logCondPostSortRank,IdxVec);

  double maxLogCondPost=logCondPostVec->data[logCondPostSortIdx->data[0]];

  double Z=0.0;
  for(unsigned int i=0;i<nChains;i++){
    Z+=exp((logCondPostVec->data[i]-maxLogCondPost)/chosenT);
  }
  
  double normConst=0.0;
  for(unsigned int currChain=0;currChain<nChains;currChain++){
    unsigned int currRank = logCondPostSortRank->data[currChain];
    if(currRank<cTrunc){
      normConst+=1.0e-10+exp((logCondPostVec->data[currChain]-maxLogCondPost)/chosenT)/Z;
    }else{
      normConst+=1.0e-10;
    }
  }

  double result=0.0;
  // Contribution from chain c1
  if(logCondPostSortRank->data[c1]<cTrunc){
    result+=log(1.0e-10+exp((logCondPostVec->data[c1]-maxLogCondPost)/chosenT)/Z)-log(normConst);
  }else{
    result+=log(1.0e-10)-log(normConst);
  }


  // Contribution from chain c2
  double chain1Prob = exp(result);
  if(logCondPostSortRank->data[c2]<cTrunc){
    result+=log(1.0e-10+exp((logCondPostVec->data[c2]-maxLogCondPost)/chosenT)/Z)-log(normConst);
  }else{
    result+=log(1.0e-10)-log(normConst);
  }
  // Renormalise c2 prob as without replacement
  result-=log(1-chain1Prob);

  gsl_permutation_free(logCondPostSortIdx);
  gsl_permutation_free(logCondPostSortRank);
  gsl_vector_free(IdxVec);
  return(result);
 
}
