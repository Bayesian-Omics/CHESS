/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
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

#include "Kernel_Single_Gamma.h"
#include "Model_Information.h"

#include <vector>

#include "../Kernel/Routines/matrix_handling.h"
#include "../Kernel/Routines/rand.h"
#include "../Kernel/Routines/moves.h"
#include "../Kernel/Routines/regression.h"
#include "../Kernel/Routines/cond_post.h"


#include "../Kernel/Classes/Double_Matrices.h"

#include "../Kernel/Classes/Temperatures.h"
#include "../Kernel/Classes/Move_monitor.h"

#include "../Kernel/Classes/DR.h"
#include "../Kernel/Classes/CM.h"

using namespace std;

#define DEBUG 0

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF



Kernel_Single_Gamma::Kernel_Single_Gamma()
{
}

void Kernel_Single_Gamma::free()
{
#if 1

    delete My_DR;
    delete My_CM;
    delete My_Move_monitor;
    delete t_tun;

    mat_log_marg.Free_double_matrix();
    mat_log_cond_post.Free_double_matrix();


    //gsl_matrix_free(matXGam); matXGam is freed outside the class kernel_Single_Gamma
    gsl_matrix_free(description_exchange_move);
    gsl_permutation_free(MyPerm);
#endif
}


Kernel_Single_Gamma::Kernel_Single_Gamma(bool postProcessOnly,
                                         unsigned int nb_chains,
                                         bool resumeRun,
                                         unsigned int resumeDRNCalls,
                                         unsigned int resumeDRNCallsAdj,
                                         vector<vector<unsigned int> > resumeDRAccepted,
                                         vector<vector<unsigned int> > resumeDRProposed,
                                         unsigned int k_max_from_read,
                                         unsigned int pX,
                                         unsigned int n_sweeps,
                                         Model_Information &Model_Data)
{



    //cout << "Enter creation of Kernel_Single_Gamma object" << endl;

#if 1

    //Regression_Model=_Regression_Model;
    Regression_Model=Model_Data.Regression_Model;

    //////////////////////////////////
    //  Setting up DR parameters
    //////////////////////////////////

    My_DR=new DR();
    if(!postProcessOnly){
      (*My_DR).set_DR(nb_chains);
      if(resumeRun){
        (*My_DR).nb_calls=resumeDRNCalls;
        (*My_DR).nb_calls_adj=resumeDRNCallsAdj;
        for(unsigned int j=0;j<nb_chains;j++){
          for(unsigned int k=0;k<nb_chains;k++){
            (*My_DR).mat_moves_accepted[j][k]=resumeDRAccepted[j][k];
            (*My_DR).mat_moves_proposed[j][k]=resumeDRProposed[j][k];
          }
        }
      }
      if (Model_Data.Model_Tag=="ESS")
      {
        (*My_DR).display_DR();
      }

    }

    //////////////////////////////////
    //  Setting up CM parameters
    //////////////////////////////////
    vector <unsigned int > list_CM_moves_enabled_from_read;
    My_CM=new CM();
    if(!postProcessOnly){
      list_CM_moves_enabled_from_read.push_back(1);
      list_CM_moves_enabled_from_read.push_back(2);


      (*My_CM).set_CM(k_max_from_read,
            list_CM_moves_enabled_from_read);
      if (Model_Data.Model_Tag=="ESS")
      {
        (*My_CM).display_CM();
      }

      list_CM_moves_enabled_from_read.clear();
    }

    //////////////////////////////////
    //  Defining Move Monitor Object
    //////////////////////////////////

    My_Move_monitor=new Move_monitor(nb_chains,
                           (*My_CM).n_possible_CM_moves);

    MyPerm=gsl_permutation_calloc(pX);
    t_tun=new Temperatures();
    description_exchange_move=description_exch_moves(nb_chains);

    mat_log_marg.Alloc_double_matrix(nb_chains,
                   n_sweeps+1);

    mat_log_cond_post.Alloc_double_matrix(nb_chains,
                    n_sweeps+1);

#endif



}

void Kernel_Single_Gamma::Initialise(bool postProcessOnly,
                                     unsigned int nb_chains,
                                     unsigned int pX,
                                     bool resumeRun,
                                     vector<vector<unsigned int> > resumeGam,
                                     unsigned int nConfounders,
                                     bool fixInit,
                                     vector < vector <unsigned int> > Gam_step_regr,
                                     bool iso_T_Flag,
                                     unsigned int maxPGamma,
                                     gsl_rng *RandomNumberGenerator,
                                     vector<unsigned int> initGam,
                                     double b_t_input,
                                     double a_t_den_inf_5k,
                                     double a_t_den_5_10k,
                                     double a_t_den_sup_10k,
                                     unsigned int temp_n_batch,
                                     vector <double> M_input,
                                     unsigned int burn_in,
                                     double temp_optimal_input,
                                     vector <double> resumeT,
                                     double resumeBT,
                                     unsigned int resumeSweep,
                                     unsigned int n_sweeps,
                                     Model_Information &Model_Data)
{

    //cout << "inside Initialise" << endl;

    if(!postProcessOnly){

        if (Model_Data.Model_Tag=="ESS")
        {
            cout << "Initialising chains" << endl;
        }

      vect_gam.resize(nb_chains);
      for(unsigned int chain=0;chain<nb_chains;chain++){
        vect_gam[chain].resize(pX);
      }

      if(resumeRun){
        for(unsigned int chain=0;chain<nb_chains;chain++){
          unsigned int k=0;
          for(unsigned int j=0;j<pX;j++)
          {
            if(resumeGam[chain].size()>0)
            {
              if(j<nConfounders)
              {
                vect_gam[chain][j]=1;
              }
              else if(j==resumeGam[chain][k])
              {
                k++;
                vect_gam[chain][j]=1;
              }
              else
              {
                vect_gam[chain][j]=0;
              }
            }
            else
            {
              if(j<nConfounders)
              {
                vect_gam[chain][j]=1;
              }
              else
              {
                vect_gam[chain][j]=0;
              }
            }
          }
        }

      }
      else
      {
        if(!fixInit)
        {
          get_vect_gam_init(vect_gam,
              Gam_step_regr,
              iso_T_Flag,
              maxPGamma,
              RandomNumberGenerator);
        }
        else
        {
          for(unsigned int chain=0;chain<nb_chains;chain++){
            unsigned int k=0;
            for(unsigned int j=0;j<pX;j++)
            {
              if(j<nConfounders)
              {
                vect_gam[chain][j]=1;
              }
              else if (initGam.size()>0) // In case there is no variable
              {
                  if(j==initGam[k])
                  {
                    //cout << "chain = " << chain << ", j = " << j << ", k = " << k << endl;
                    k++;
                    vect_gam[chain][j]=1;
                  }
              }
              else
              {
                vect_gam[chain][j]=0;
              }
            }
          }
        }
      }


      Gam_step_regr.clear();
    }


    if(!postProcessOnly){
      //////////////////////////////////
      // Intializing chain temparature
      //////////////////////////////////

      (*t_tun).set_Temp_param(nb_chains,
                pX,
                b_t_input,
                a_t_den_inf_5k,
                a_t_den_5_10k,
                a_t_den_sup_10k,
                temp_n_batch,
                M_input,
                burn_in,
                temp_optimal_input,
                iso_T_Flag);
      M_input.clear();
      if(resumeRun){
        for(unsigned int j=0;j<nb_chains;j++){
          (*t_tun).t[j]=resumeT[j];
        }
        (*t_tun).b_t=resumeBT;
      }
      if (Model_Data.Model_Tag=="ESS")
      {
            (*t_tun).display_Temp_param();
      }

      intialize_chain_idx(chain_idx,
                nb_chains);
      if(resumeRun){
        for(unsigned int c=0;c<nb_chains;c++){
          chain_idx[c]=resumeChainIndex[c];
        }
      }

      //description_exchange_move: each move is presented in columns
      // -line 1: first chain
      // -line 2: second chain
      // -line 3: Pbty of the move

    }

    if(postProcessOnly)
    {
      n_Models_visited.assign(resumeSweep+1,0);
    }else
    {
      n_Models_visited.assign(n_sweeps+1,0);
    }
    if(resumeRun||postProcessOnly)
    {
      if(resumeSweep>burn_in)
      {
        List_Models.resize(resumeSweep-burn_in);

      }

      // There seems to be a disagreement between ESS and HESS about when to start saving models...
      // ESS: save after burnin
      // HESS: save from first sweep
      for(unsigned int i=0;i<resumeSweep+1;i++){
        n_Models_visited[i]=pastNModelsVisited[i];
        if(i>burn_in){
          unsigned int tmp = pastModels[i][0];
          List_Models[i-1-burn_in].resize(tmp+1);
          List_Models[i-1-burn_in][0]=tmp;
          for(unsigned int j=0;j<tmp;j++){
            List_Models[i-1-burn_in][j+1]=pastModels[i][j+1];
          }
        }
      }

    }

}


void Kernel_Single_Gamma::First_Sweep(unsigned int Response,
                                      unsigned int nb_chains,
                                      gsl_matrix *mat_X_work2,
                                      gsl_matrix *mat_Y_work2,
                                      Prior_param PR,
                                      bool gPriorFlag,
                                      bool indepPriorFlag,
                                      bool gSampleFlag,
                                      double lambda,
                                      double g,
                                      //unsigned int pX,
                                      //unsigned int nX,
                                      //unsigned int pY,
                                      bool cudaFlag,
                                      unsigned int sweep,
                                      Model_Information &Model_Data)
{


    for(unsigned int chain=0;chain<nb_chains;chain++){

      //Initial Calculation of the logMarg and log_cond_post;

      double logMargLik,logPost;

      computeLogPosterior(vect_gam[chain_idx[chain]],
                          chain,
                          Response,
                          mat_X_work2,
                          logMargLik,
                          logPost,
                          //matXGam,
                          mat_Y_work2,
                          PR,
                          gPriorFlag,
                          indepPriorFlag,
                          gSampleFlag,
                          lambda,
                          g,
                          //pX,
                          //n_vars_in,
                          //nX,
                          //pY,
                          cudaFlag,
                          Model_Data);

      mat_log_marg.matrix[chain_idx[chain]][sweep]=logMargLik;
      mat_log_cond_post.matrix[chain_idx[chain]][sweep]=logPost;


      if (Model_Data.Model_Tag=="ESS")
      {

          vector < unsigned int > list_columns_X_gam;
          vector < unsigned int > list_columns_X_gam_bar;

          get_list_var_in_and_out(list_columns_X_gam,
                    list_columns_X_gam_bar,
                    vect_gam[chain_idx[chain]]);

          cout << "List_columns_X" << endl;
          for(unsigned int col=0;col<list_columns_X_gam.size();col++){
            cout << list_columns_X_gam[col] << " ";
          }
          cout << endl;
          unsigned int n_vars_in=list_columns_X_gam.size();
          cout << "n_vars_in " << n_vars_in << endl;


          cout << endl << "**************Results, chain " << chain << " -- sweep " << sweep << "***************" << endl;
          cout << "\tlog_marg " << mat_log_marg.matrix[chain_idx[chain]][sweep] << endl
              << "\tlog_cond_post " << mat_log_cond_post.matrix[chain_idx[chain]][sweep] << endl;
          cout << "********************************************************" << endl;

          list_columns_X_gam.clear();
          list_columns_X_gam_bar.clear();

      }


    }//end of for chain.

}


void Kernel_Single_Gamma::MCMC_Update(unsigned int Response,
                                      unsigned int sweep,
                                      unsigned int nb_chains,
                                      unsigned int Gibbs_n_batch,
                                      bool Log_Flag,
                                      gsl_matrix *mat_X_work2,
                                      gsl_matrix *mat_Y_work2,
                                      bool gPriorFlag,
                                      bool indepPriorFlag,
                                      bool gSampleFlag,
                                      double lambda,
                                      double g,
                                      Prior_param PR,
                                      bool cudaFlag,
                                      unsigned int nConfounders,
                                      unsigned int maxPGamma,
                                      gsl_rng *RandomNumberGenerator,
                                      unsigned int burn_in,
                                      //unsigned int nX,
                                      bool iso_T_Flag,
                                      Model_Information &Model_Data)
{

    unsigned int nX=mat_X_work2->size1;

    // Reset the permutation object
    gsl_permutation_init(MyPerm);

    n_Models_visited[sweep]=n_Models_visited[sweep-1];


    ////////////////////////////////////////////////////
    /////////////////// Local Moves  ///////////////////
    ////////////////////////////////////////////////////
    for(unsigned int chain=0;chain<nb_chains;chain++){
      mat_log_marg.matrix[chain][sweep]=mat_log_marg.matrix[chain][sweep-1];
      mat_log_cond_post.matrix[chain][sweep]=mat_log_cond_post.matrix[chain][sweep-1];
    }

    //Gibbs Moves
    if(sweep%Gibbs_n_batch==0){
      if(Log_Flag){
        cout << "Gibbs" << endl;
      }

      Gibbs_move(Response,
                mat_log_marg,
                mat_log_cond_post,
                mat_X_work2,
                mat_Y_work2,
                t_tun,
                MyPerm,
                vect_gam,
                sweep,
                gPriorFlag,
                indepPriorFlag,
                gSampleFlag,
                lambda,
                g,
                PR,
                My_Move_monitor,
                chain_idx,
                n_Models_visited,
                cudaFlag,
                nConfounders,
                maxPGamma,
                RandomNumberGenerator,
                Model_Data);

    }

    ////Fast Scan Metropolis Hastings (FSMH)

    double local_move_rand=myrand(RandomNumberGenerator);

    if(nb_chains==1){
      local_move_rand=0.0;
    }

    if(local_move_rand<PR.Prob_mut){
      if(Log_Flag){
        cout << "FSMH" << endl;
      }
      FSMH_move(Response,
                  mat_log_marg,
                  mat_log_cond_post,
                  mat_X_work2,
                  mat_Y_work2,
                  t_tun,
                  MyPerm,
                  vect_gam,
                  sweep,
                  gPriorFlag,
                  indepPriorFlag,
                  gSampleFlag,
                  lambda,
                  g,
                  PR,
                  My_Move_monitor,
                  chain_idx,
                  n_Models_visited,
                  cudaFlag,
                  nConfounders,
                  maxPGamma,
                  RandomNumberGenerator,
                  Model_Data);

    }
    else{
      ///////////////////////////////////////////////////////
      /////////////////////// CM Moves  /////////////////////
      ///////////////////////////////////////////////////////

      //Defining the number of CO move to simulate
      if(Log_Flag){
        cout << "CM" << endl;
      }
      Crossover_move(Response,
           mat_log_marg,
           mat_log_cond_post,
           mat_X_work2,
           mat_Y_work2,
           t_tun,
           vect_gam,
           gPriorFlag,
           indepPriorFlag,
           gSampleFlag,
           lambda,
           g,
           PR,
           My_Move_monitor,
           chain_idx,
           sweep,
           My_CM,
           n_Models_visited,
           cudaFlag,
           maxPGamma,
           RandomNumberGenerator,
           Model_Data);
    }

    ///////////////////////////////////////////////////////
    ///////////////////// Global Moves  ///////////////////
    ///////////////////////////////////////////////////////
    if(nb_chains>1){
      double global_move_rand=myrand(RandomNumberGenerator);
      if(sweep<=burn_in){
        global_move_rand=0.0;//during burn-in, DR is the only global move
      }

      if(global_move_rand<PR.Prob_DR){

        if(Log_Flag){
          cout << "DR" << endl;
        }


        DR_move(My_DR,
      chain_idx,
      mat_log_cond_post,
      t_tun,
      sweep,
      My_Move_monitor,
      RandomNumberGenerator,
      Model_Data);

      }else{
        if(Log_Flag){
          cout << "AE" << endl;
        }
        All_exchange_move(description_exchange_move,
            chain_idx,
            mat_log_cond_post,
            t_tun,
            sweep,
            My_Move_monitor,
            RandomNumberGenerator,
            Model_Data);
      }

    }



    ///////////////////////////////////////////////////
    ///////////////// Temp placement //////////////////
    ///////////////////////////////////////////////////
    if(nb_chains>1){
      if(sweep<=burn_in){
        if(((*My_DR).nb_calls_adj==(*t_tun).nbatch) || ((*My_DR).nb_calls==5*(*t_tun).nbatch)){
          unsigned int n_vars_in_last_chain=sum_line_std_mat(vect_gam,
                               chain_idx[nb_chains-1]);


          //cout << "t_placement" << endl;

          temp_placement(t_tun,
           My_DR,
           My_Move_monitor,
           sweep,
           n_vars_in_last_chain,
           nX,
           iso_T_Flag);
          //cout << "end -- t_placement" << endl;

        }
      }
    }

}


void Kernel_Single_Gamma::Stepwise_Regression(vector < vector <unsigned int> > & Gam_step_regr,
                                              unsigned int nConfounders,
                                              gsl_matrix * mat_X,
                                              gsl_matrix * mat_Y,
                                              gsl_vector * vect_RMSE,
                                              double Pvalue_enter,
                                              double Pvalue_remove,
                                              Model_Information &Model_Data)
{
    if (DEBUG)
    { cout << "enter stepwise regression function" << endl << endl; }

    unsigned int nX=mat_X->size1;
    unsigned int pX=mat_X->size2;

    unsigned int Nb_Outcomes=mat_Y->size2;

    Double_Matrices Gam_step_regr_pvals;
    Gam_step_regr_pvals.Alloc_double_matrix(Nb_Outcomes,
                    pX);
    Double_Matrices Gam_step_regr_SE;
    Gam_step_regr_SE.Alloc_double_matrix(Nb_Outcomes,
                     pX);

    Double_Matrices Gam_step_regr_beta;
    Gam_step_regr_beta.Alloc_double_matrix(Nb_Outcomes,
                   pX);

    double tolerance= 6.6835e-14;

    gsl_vector *current_outcome =gsl_vector_calloc(nX);
    gsl_vector *vect_residuals=gsl_vector_calloc(nX);
    gsl_vector *vect_p_value=gsl_vector_calloc(pX);
    gsl_vector *vect_beta_full=gsl_vector_calloc(pX);
    gsl_vector *vect_SE_full=gsl_vector_calloc(pX);

    if (Model_Data.Model_Tag=="ESS")
    {
        cout << "***************************************************" << endl
            << "*************  Stepwise regression   *************" << endl
            << "***************************************************" << endl << endl;
    }

    for(unsigned int outcome=0;outcome<Nb_Outcomes;outcome++)
    {
      vector < unsigned int > list_columns_X_gam;
      vector < unsigned int > list_columns_X_gam_bar;
      vector < unsigned int > is_var_in;

      is_var_in.resize(pX);
      for(unsigned int jj=0;jj<nConfounders;jj++){
        // Set up so confounders are always in
        is_var_in[jj]=1;
      }


      gsl_matrix_get_col(current_outcome,
             mat_Y,
             outcome);
      int stop=0;
      int loop=0;
      int n_loop_max=100;

      while(stop==0){

        get_list_var_in_and_out(list_columns_X_gam,
                list_columns_X_gam_bar,
                is_var_in);

        gsl_matrix *mat_X_gam=get_X_reduced_and_constant(list_columns_X_gam,
                          mat_X);

        gsl_matrix *mat_X_gam_bar=get_X_reduced(list_columns_X_gam_bar,
                         mat_X);

        //Calculating p-values for all variables
        get_full_pvalues_and_beta(vect_p_value,
              vect_beta_full,
              vect_SE_full,
              mat_X_gam,
              mat_X_gam_bar,
              current_outcome,
              vect_residuals,
              vect_RMSE,
              tolerance,
              list_columns_X_gam,
              list_columns_X_gam_bar,
              outcome);

        //Updating the list of var in
        stop=update_is_var_in(is_var_in,
              list_columns_X_gam,
              list_columns_X_gam_bar,
              vect_p_value,
              Pvalue_enter,
              Pvalue_remove,
              loop,
              n_loop_max,
              nConfounders);

        gsl_matrix_free(mat_X_gam);
        gsl_matrix_free(mat_X_gam_bar);

        loop++;

        list_columns_X_gam.clear();
        list_columns_X_gam_bar.clear();

        if(stop==1){
          get_list_var_in_and_out(list_columns_X_gam,
              list_columns_X_gam_bar,
              is_var_in);
        }

      }//end of while
      //Filling the output file

      //Inserted line:
      Gam_step_regr.resize(Nb_Outcomes);
      if (DEBUG)
      {cout << "Matrix resized correctly" << endl; }
      store_model_per_outcome(Gam_step_regr,
              list_columns_X_gam,
              vect_p_value,
              vect_beta_full,
              vect_SE_full,
              Gam_step_regr_pvals,
              Gam_step_regr_SE,
              Gam_step_regr_beta,
              outcome);
      if (DEBUG)
      {cout << "After store_model_per_outcome" << endl; }


      list_columns_X_gam.clear();
      is_var_in.clear();
    }

    if (Model_Data.Model_Tag=="ESS")
    {
        cout << "Result From Step-Wise Regression" << endl;
        display_matrix_var_dim(Gam_step_regr);

        cout << endl;

        cout << "***************************************************" << endl
            << "**********  End of Stepwise regression  **********" << endl
            << "***************************************************" << endl << endl;
    }

    gsl_vector_free(vect_residuals);
    gsl_vector_free(current_outcome);
    gsl_vector_free(vect_p_value);
    gsl_vector_free(vect_beta_full);
    gsl_vector_free(vect_SE_full);

    Gam_step_regr_pvals.Free_double_matrix();
    Gam_step_regr_SE.Free_double_matrix();
    Gam_step_regr_beta.Free_double_matrix();

}




void Kernel_Single_Gamma::save_model_per_sweep()
{
    // FOR THE MOMENT, THIS FUNCTION IS ONLY USED IN HESS.
    // NEED TO BE CAREFUL BEFORE USING IT IN ESS

    // We save both the model, and its posterior probability conditional on Omega and g,
    // evaluated at the current sweep.

    /*
    ########################################################################
                               0)   PRELIMINARIES
    ########################################################################
    */

    unsigned int pos_chain=chain_idx[0];

    List_Models.push_back(vector<unsigned int>());

    unsigned int last_pos_in_List=List_Models.size()-1;

    /*
    ########################################################################
                               1)   SAVE MODEL
    ########################################################################
    */

    List_Models[last_pos_in_List].push_back(0);//room for the #of variables in

    unsigned int count_n_vars_in=0;

    for(unsigned int var=0;var<vect_gam[0].size();var++)
    {
        if(vect_gam[pos_chain][var]==1)
        {
            //unsigned int last_pos_in_List=List_Models.size()-1;
            List_Models[last_pos_in_List].push_back(var);
            count_n_vars_in++;
        }
    }
    //unsigned int last_pos_in_List=List_Models.size()-1;
    List_Models[last_pos_in_List][0]=count_n_vars_in;

    /*
    ########################################################################
                               2)   SAVE...
    ########################################################################
    */


}
