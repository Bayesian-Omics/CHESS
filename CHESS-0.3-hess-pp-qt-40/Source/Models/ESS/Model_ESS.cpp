/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from ESS.cc in the ESS++ program
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

#include "Start_Up_ESS.h"
#include "../../General_Classes/Model_Information.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>

#include <ctime>

#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <sstream>

#include "../../General_Classes/Kernel_Single_Gamma.h"
#include "../../General_Classes/Command_Line.h"
#include "../../General_Classes/Model_Generic.h"
#include "Model_ESS.h"

#include "../../Kernel/Routines/struc.h"
#include "../../Kernel/Routines/dyn_name.h"
#include "../../Kernel/Routines/matrix_handling.h"
#include "../../Kernel/Routines/rand.h"
#include "../../Kernel/Routines/moves.h"
#include "../../Kernel/Routines/regression.h"
#include "../../Kernel/Routines/cond_post.h"
#include "../../Kernel/Routines/post_processing.h"
#include "../../Kernel/Routines/xml_file_read.h"
#include "../../Kernel/Classes/String_Matrices.h"
#include "../../Kernel/Classes/Int_Matrices.h"
#include "../../Kernel/Classes/Int_Matrices_var_dim.h"
#include "../../Kernel/Classes/Double_Matrices.h"
#include "../../Kernel/Classes/Double_Matrices_cont.h"
#include "../../Kernel/Classes/Prior_param.h"
#include "../../Kernel/Classes/Temperatures.h"
#include "../../Kernel/Classes/Move_monitor.h"
#include "../../Kernel/Classes/AdMH.h"
#include "../../Kernel/Classes/DR.h"
#include "../../Kernel/Classes/CM.h"
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula.h>
#endif


#define DEBUG 0
#define Missing_data -1

Model_ESS::Model_ESS()
{

    Model_Tag="ESS";
    Regression_Model="Undefined";

}

void Model_ESS::Run()
{

    Model_Information ESS_Information(Model_Tag, Regression_Model);

    time_t beginTime, currTime;
    beginTime = time(NULL);

    std::clock_t startTime = std::clock();
    std::clock_t tmpTime,endTime;
    double setupTime=0.0,mainLoopTime=0.0,postProcessTime=0.0;

    gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_mt19937 );


    /*   File Streams
       -------------*/
        ofstream f_out;
        ofstream f_out_n_vars_in;
        ofstream f_out_n_models_visited;
        ofstream f_out_log_cond_post_per_chain;
        ofstream f_out_FSMH;
        ofstream f_out_CM;
        ofstream f_out_AE;
        ofstream f_out_DR;
        ofstream f_out_g;
        ofstream f_out_g_adapt;
        ofstream f_out_Gibbs;
        ofstream f_out_t_tun;
        ostringstream strStrOut;
        ostringstream strStrOutNVarsIn;
        ostringstream strStrOutNModelsVisited;
        ostringstream strStrOutLogCondPostPerChain;
        ostringstream strStrOutFSMH;
        ostringstream strStrOutCM;
        ostringstream strStrOutAE;
        ostringstream strStrOutDR;
        ostringstream strStrOutG;
        ostringstream strStrOutGAdapt;
        ostringstream strStrOutGibbs;
        ostringstream strStrOutTTun;
        fstream f_resume;
        FILE *f_rng;
        ofstream f_out_time;
        ostringstream strStrOutTime;
        ofstream f_out_best_models;
        ofstream f_out_marg_gam;

    /* Other modeling parameters
       -----------------------------------*/

        double temp_optimal_input;
        double b_t_input;
        double a_t_den_inf_5k;
        double a_t_den_5_10k;
        double a_t_den_sup_10k;
        unsigned int temp_n_batch;
        unsigned int maxPGamma;
        unsigned int nb_chains;
        unsigned int nX;
        unsigned int pX;
        unsigned int nY;
        unsigned int pY;

        Prior_param PR;

        unsigned int Gibbs_n_batch;

        string resumeName;
        string rngFileName;
        string OutputName_time;


    /* Initial state
       -------------*/
        vector<unsigned int> initGam;
        vector<unsigned int> resumeChainIndex;
        vector <double> resumeT;
        double resumeCumG;
        double resumeCountG;
        unsigned int resumeSweep;
        vector<vector<unsigned int> > resumeGam;
        double resumeBT;
        unsigned int resumeDRNCalls;
        unsigned int resumeDRNCallsAdj;
        vector<vector<unsigned int> > resumeDRAccepted;
        vector<vector<unsigned int> > resumeDRProposed;
        unsigned int k_max_from_read;

    /* Stepwise Regression object
       --------------------------*/
        vector < vector <unsigned int> > Gam_step_regr;


    /* Internal PARAMETERS
       -------------*/
        vector< unsigned int > shuffleYIndex;
        vector <double> M_input;

    /* MCMC state
       -------------*/
        double g=0;
        vector<vector<unsigned int> > pastModels;
        vector<unsigned int> pastNModelsVisited;


    /**********************************************************
    ***********************************************************

                        START PROCESSING

    ***********************************************************
    ***********************************************************/



      Start_Up_ESS* Initial_Parameters=new Start_Up_ESS();

      Initial_Parameters->Initialise(
                    RandomNumberGenerator,
                    f_out,
                    f_out_n_vars_in,
                    f_out_n_models_visited,
                    f_out_log_cond_post_per_chain,
                    f_out_FSMH,
                    f_out_CM,
                    f_out_AE,
                    f_out_DR,
                    f_out_g,
                    f_out_g_adapt,
                    f_out_Gibbs,
                    f_out_t_tun,
                    f_resume,
                    f_rng,
                    f_out_time,
                    f_out_best_models,
                    f_out_marg_gam,
                    postProcessOnly,
                    resumeRun,
                    fixInit,
                    Gam_step_regr,
                    iso_T_Flag,
                    maxPGamma,
                    initGam,
                    nb_chains,
                    nConfounders,
                    pX,
                    nY,
                    shuffleYIndex,
                    b_t_input,
                    a_t_den_inf_5k,
                    a_t_den_5_10k,
                    a_t_den_sup_10k,
                    temp_n_batch,
                    M_input,
                    burn_in,
                    temp_optimal_input,
                    resumeT,
                    resumeChainIndex,
                    PR,
                    gPriorFlag,
                    indepPriorFlag,
                    gSampleFlag,
                    lambda,
                    g,
                    nX,
                    pY,
                    cudaFlag,
                    Log_Flag,
                    timeLimit,
                    Time_monitorFlag,
                    Out_full_Flag,
                    HistoryFlag,
                    resumeCumG,
                    resumeCountG,
                    resumeSweep,
                    pastModels,
                    pastNModelsVisited,
                    Gibbs_n_batch,
                    resumeGam,
                    resumeBT,
                    n_sweeps,
                    resumeName,
                    rngFileName,
                    OutputName_time,
                    n_top_models,
                    resumeDRNCalls,
                    resumeDRNCallsAdj,
                    resumeDRAccepted,
                    resumeDRProposed,
                    k_max_from_read,
                    path_name_out,
                    MY_SEED,
                    standardizeFlag,
                    g_init,
                    filename_init,
                    doYShuffle,
                    filename_in_mat_X,
                    filename_in_mat_Y,
                    filename_par);

      /**********************************************************
      ***********************************************************

                DECLARE NOW THE Kernel_Single_Gamma object

      ***********************************************************
      ***********************************************************/


    Kernel_Single_Gamma Single_Gamma(postProcessOnly,
                                 nb_chains,
                                 resumeRun,
                                 resumeDRNCalls,
                                 resumeDRNCallsAdj,
                                 resumeDRAccepted,
                                 resumeDRProposed,
                                 k_max_from_read,
                                 pX,
                                 n_sweeps,
                                 ESS_Information);

    gsl_matrix * mat_X_work2=Initial_Parameters->mat_X_work2;
    gsl_matrix * mat_Y_work2=Initial_Parameters->mat_Y_work2;

    gsl_vector * vect_RMSE=Initial_Parameters->vect_RMSE;

    AdMH * My_g_AdMH=Initial_Parameters->My_g_AdMH;

    delete Initial_Parameters;

    unsigned int sweep=0;

    Single_Gamma.Initialise(postProcessOnly,
                      nb_chains,
                      pX,
                      resumeRun,
                      resumeGam,
                      nConfounders,
                      fixInit,
                      Gam_step_regr,
                      iso_T_Flag,
                      maxPGamma,
                      RandomNumberGenerator,
                      initGam,
                      b_t_input,
                      a_t_den_inf_5k,
                      a_t_den_5_10k,
                      a_t_den_sup_10k,
                      temp_n_batch,
                      M_input,
                      burn_in,
                      temp_optimal_input,
                      resumeT,
                      resumeBT,
                      resumeSweep,
                      n_sweeps,
                      ESS_Information);


    if(!postProcessOnly)
    {
      if(resumeRun){
        sweep=resumeSweep;
      }
    }

    /*
    ####################################################################################

                  End of initialisation

    ####################################################################################
    */

    endTime = std::clock();
    tmpTime = startTime;
    setupTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);


    if(!postProcessOnly){

      cout << "***************************************" << endl
          << "             FIRST SWEEP" << endl
          <<  "***************************************" << endl << endl;

      unsigned int Response_Parameter=0;
      Single_Gamma.First_Sweep(Response_Parameter,
                               nb_chains,
                               mat_X_work2,
                               mat_Y_work2,
                               PR,
                               gPriorFlag,
                               indepPriorFlag,
                               gSampleFlag,
                               lambda,
                               g,
                               //pX,
                               //nX,
                               //pY,
                               cudaFlag,
                               sweep,
                               ESS_Information);


    cout << "***************************************" << endl
        << "          END OF FIRST SWEEP" << endl
        <<  "***************************************" << endl << endl;


        if(!resumeRun){
          print_main_results_per_sweep(f_out,
                      Single_Gamma.vect_gam,
                      Single_Gamma.chain_idx,
                      Single_Gamma.mat_log_marg,
                      Single_Gamma.mat_log_cond_post,
                      0);
        }

      }

      double cum_g=0.0;
      unsigned int count_g=0.0;
      if(resumeRun||postProcessOnly){
        cum_g=resumeCumG;
        count_g=resumeCountG;
      }


      if(!postProcessOnly){
        for(sweep=resumeSweep+1;sweep<n_sweeps+1;sweep++)
        {

          if(Log_Flag){
            cout << "***************************************" << endl
                << "             SWEEP #" << sweep << "/" << n_sweeps << endl
                << "***************************************" << endl << endl;
          }

          clock_t Time_start,Time_end;
          Time_start=clock();

          unsigned int Response_Parameter=0;

          Single_Gamma.MCMC_Update(Response_Parameter,
                      sweep,
                      nb_chains,
                      Gibbs_n_batch,
                      Log_Flag,
                      mat_X_work2,
                      mat_Y_work2,
                      gPriorFlag,
                      indepPriorFlag,
                      gSampleFlag,
                      lambda,
                      g,
                      PR,
                      cudaFlag,
                      nConfounders,
                      maxPGamma,
                      RandomNumberGenerator,
                      burn_in,
                      //nX,
                      iso_T_Flag,
                      ESS_Information);


          /////////////////////////////////////////////////////
          ///////////////////// Sampling g  ///////////////////
          /////////////////////////////////////////////////////

          if(gSampleFlag){
            sample_g_ESS(Response_Parameter,
               Single_Gamma.mat_log_marg,
               Single_Gamma.mat_log_cond_post,
               My_g_AdMH,
               Single_Gamma.t_tun,
               Single_Gamma.vect_gam,
               mat_X_work2,
               mat_Y_work2,
               gPriorFlag,
               indepPriorFlag,
               gSampleFlag,
               lambda,
               g,
               PR,
               sweep,
               Single_Gamma.My_Move_monitor,
               Single_Gamma.chain_idx,
               Single_Gamma.n_Models_visited,
               cudaFlag,
               RandomNumberGenerator,
               //Model_Tag,
               //Regression_Model,
               ESS_Information);

            if(Log_Flag){
              cout << "g = " << g << endl;
            }

          }

          if(Log_Flag){
              display_summary_result_per_sweep(Single_Gamma.vect_gam,
                            Single_Gamma.chain_idx,
                            Single_Gamma.mat_log_marg,
                            Single_Gamma.mat_log_cond_post,
                            sweep,
                            Single_Gamma.t_tun,
                            nConfounders);
          }

          if(timeLimit>0){
              print_and_save_main_results_per_sweep(f_out,
                           f_out_n_vars_in,
                           f_out_n_models_visited,
                           f_out_log_cond_post_per_chain,
                           Single_Gamma.vect_gam,
                           Single_Gamma.List_Models,
                           Single_Gamma.chain_idx,
                           Single_Gamma.mat_log_marg,
                           Single_Gamma.mat_log_cond_post,
                           sweep,
                           burn_in,
                           Single_Gamma.n_Models_visited[sweep],
                           HistoryFlag);

              if(HistoryFlag){
                (*Single_Gamma.My_Move_monitor).print_move_monitor_per_sweep(f_out_FSMH,
                                                         f_out_CM,
                                                         f_out_AE,
                                                         f_out_DR,
                                                         f_out_g,
                                                         f_out_g_adapt,
                                                         f_out_Gibbs,
                                                         f_out_t_tun,
                                                         gSampleFlag,
                                                         iso_T_Flag,
                                                         sweep,
                                                         nb_chains);

              }
          }else{
            print_and_save_main_results_per_sweep(strStrOut,
                                                     strStrOutNVarsIn,
                                                     strStrOutNModelsVisited,
                                                     strStrOutLogCondPostPerChain,
                                                     Single_Gamma.vect_gam,
                                                     Single_Gamma.List_Models,
                                                     Single_Gamma.chain_idx,
                                                     Single_Gamma.mat_log_marg,
                                                     Single_Gamma.mat_log_cond_post,
                                                     sweep,
                                                     burn_in,
                                                     Single_Gamma.n_Models_visited[sweep],
                                                     HistoryFlag);

            if(HistoryFlag){
              (*Single_Gamma.My_Move_monitor).print_move_monitor_per_sweep(strStrOutFSMH,                                                           strStrOutCM,
                                                               strStrOutAE,
                                                               strStrOutDR,
                                                               strStrOutG,
                                                               strStrOutGAdapt,
                                                               strStrOutGibbs,
                                                               strStrOutTTun,
                                                               gSampleFlag,
                                                               iso_T_Flag,
                                                               sweep,
                                                               nb_chains);

            }
          }

          if(sweep>burn_in){
            cum_g+=g;
            count_g++;
          }

          if(timeLimit>0){
            f_resume.open(resumeName.c_str(),ios::out);
            f_rng=fopen(rngFileName.c_str(),"wb");

            saveResumeFile(f_resume,f_rng,sweep,g,shuffleYIndex,nY,cum_g,count_g,vect_RMSE,Single_Gamma.t_tun,My_g_AdMH->ls,Single_Gamma.My_DR,
                              Single_Gamma.vect_gam,Single_Gamma.chain_idx,nConfounders,pY,RandomNumberGenerator);

            fclose(f_rng);
            f_resume.close();
          }

          Time_end=clock();
          double time_taken=((double)Time_end-Time_start)/CLOCKS_PER_SEC;
          if(Log_Flag){
            cout << "Sweep\tTime\t#models\tTime/model" << endl;
            cout << sweep << "\t"
                << time_taken << "\t"
                << Single_Gamma.n_Models_visited[sweep]-Single_Gamma.n_Models_visited[sweep-1] << "\t"
                << time_taken/(double)(Single_Gamma.n_Models_visited[sweep]-Single_Gamma.n_Models_visited[sweep-1])
                <<  endl;

            cout << "***************************************" << endl
                << "           END OF SWEEP #" << sweep << "/" << n_sweeps << endl
                <<  "***************************************" << endl << endl;
          }

          if(Time_monitorFlag){
            if(timeLimit>0){
              f_out_time << sweep << "\t" << time_taken << "\t" << time_taken/(double)(Single_Gamma.n_Models_visited[sweep]-Single_Gamma.n_Models_visited[sweep-1]) << endl;
            }else{
              strStrOutTime << sweep << "\t" << time_taken << "\t" << time_taken/(double)(Single_Gamma.n_Models_visited[sweep]-Single_Gamma.n_Models_visited[sweep-1]) << endl;
            }
          }

          currTime = time(NULL);
          double timeElapsed=((double)currTime-beginTime)/3600;
          //NEW
          if(Time_monitorFlag)
          {
            cout << "Time elapsed: " << timeElapsed << endl;
          }
          if(timeLimit>0&&timeElapsed>timeLimit){
            // Time has passed the limit so we need to stop
            if(HistoryFlag){
              f_out.close();
              f_out_n_vars_in.close();
              f_out_log_cond_post_per_chain.close();
              f_out_n_models_visited.close();
              f_out_FSMH.close();
              f_out_CM.close();
              f_out_AE.close();
              f_out_DR.close();
              if(gSampleFlag){
                f_out_g.close();
                f_out_g_adapt.close();
              }
              f_out_Gibbs.close();
              if(iso_T_Flag==false){
                f_out_t_tun.close();
              }
            }
            if(Time_monitorFlag){
              f_out_time.close();
            }

            //n_Models_visited.clear();
            Single_Gamma.n_Models_visited.clear();

            for(unsigned int chain=0;chain<nb_chains;chain++){
              //vect_gam[chain].clear();
              Single_Gamma.vect_gam[chain].clear();
            }
            //vect_gam.clear();
            Single_Gamma.vect_gam.clear();
            gsl_matrix_free(mat_X_work2);
            gsl_matrix_free(mat_Y_work2);

            gsl_vector_free(vect_RMSE);
            //mat_log_marg.Free_double_matrix();
            //mat_log_cond_post.Free_double_matrix();
            Single_Gamma.mat_log_marg.Free_double_matrix();
            Single_Gamma.mat_log_cond_post.Free_double_matrix();
            //gsl_permutation_free(MyPerm);
            //chain_idx.clear();
            //gsl_matrix_free(description_exchange_move);
            //delete (My_g_AdMH);
            //delete (My_Move_monitor);
            //delete (My_CM);
            //delete (My_DR);
            //delete (t_tun);
            gsl_permutation_free(Single_Gamma.MyPerm);
            Single_Gamma.chain_idx.clear();
            gsl_matrix_free(Single_Gamma.description_exchange_move);
            delete (My_g_AdMH);
            delete (Single_Gamma.My_Move_monitor);
            delete (Single_Gamma.My_CM);
            delete (Single_Gamma.My_DR);
            delete (Single_Gamma.t_tun);

    #if _CUDA_
            if(cudaFlag){
              cublasShutdown();
              culaShutdown();
            }
    #endif

            cout << "Run curtailed as time limit exceeded. Please run with -resume flag." << endl;
            //return(0);

          }


        }//end of for sweep
      }


      // Now if we weren't writing as we went along, write the output
      if(!postProcessOnly){
        if(timeLimit<0){
          if(HistoryFlag){
            f_out << strStrOut.str();
            f_out_n_vars_in << strStrOutNVarsIn.str();
            f_out_log_cond_post_per_chain << strStrOutLogCondPostPerChain.str();
            f_out_n_models_visited << strStrOutNModelsVisited.str();
            f_out_Gibbs << strStrOutGibbs.str();
            f_out_FSMH << strStrOutFSMH.str();
            f_out_CM << strStrOutCM.str();
            f_out_AE << strStrOutAE.str();
            f_out_DR << strStrOutDR.str();
            if(gSampleFlag){
              f_out_g << strStrOutG.str();
              f_out_g_adapt << strStrOutGAdapt.str();
            }
            if(!iso_T_Flag){
              f_out_t_tun << strStrOutTTun.str();
            }
            if(Time_monitorFlag){
              f_out_time << strStrOutTime.str();
            }
          }
        }

        if(HistoryFlag){
          f_out.close();
          f_out_n_vars_in.close();
          f_out_log_cond_post_per_chain.close();
          f_out_n_models_visited.close();
          f_out_FSMH.close();
          f_out_CM.close();
          f_out_AE.close();
          f_out_DR.close();
          if(gSampleFlag){
            f_out_g.close();
            f_out_g_adapt.close();
          }
          f_out_Gibbs.close();
          if(!iso_T_Flag){
            f_out_t_tun.close();
          }
        }
        if(Time_monitorFlag){
          f_out_time.close();
        }

        cout << "***************************************" << endl
         << "***************************************" << endl
         << "           END OF SWEEPS" << endl
         //<< "       # models evaluated: " << n_Models_visited[n_sweeps] << endl
         << "       # models evaluated: " << Single_Gamma.n_Models_visited[n_sweeps] << endl
         << "***************************************" << endl
         <<  "***************************************" << endl << endl;

        tmpTime = endTime;
        endTime = std::clock();
        mainLoopTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);
      }


      //////////////////
      //Post-Processing
      //////////////////

      //Step1: Calculating E(g);
      double mean_g=cum_g/(double)(count_g);
      //Step2: Get the Unique list of models and integrate log posterior

      vector< vector<double> >vect_marg_log_post;
      int pos_null_model;
      bool Include_Single_Models=1;
              pos_null_model=getUniqueList(Single_Gamma.Unique_List_Models,
                                                       Single_Gamma.List_Models,
                                                       Single_Gamma.n_Models_visited,
                                                       burn_in,
                                                       pX,
                                                       nConfounders,
                                                       Include_Single_Models);

              getLogPost(vect_marg_log_post,
                        Single_Gamma.Unique_List_Models,
                        mat_X_work2,
                        mat_Y_work2,
                        PR,
                        gPriorFlag,
                        indepPriorFlag,
                        gSampleFlag,
                        lambda,
                        mean_g,
                        cudaFlag,
                        ESS_Information);

      //Step4: Get the posterior of the model
      //unsigned int n_unique = Unique_List_Models.size();
        unsigned int n_unique = Single_Gamma.Unique_List_Models.size();
        gsl_permutation *idx_post_gam_sort=gsl_permutation_calloc(vect_marg_log_post.size());
        gsl_vector *vect_post_gam = gsl_vector_calloc(vect_marg_log_post.size());

      unsigned int n_retained_models=min(n_unique,n_top_models);


      getAndSortPostGam(vect_post_gam,
                       idx_post_gam_sort,
                       vect_marg_log_post);

      combineAndPrintBestModel(f_out_best_models,
                       idx_post_gam_sort,
                       vect_post_gam,
                       vect_marg_log_post,
                       //Unique_List_Models,
                       Single_Gamma.Unique_List_Models,
                       n_retained_models,
                       pos_null_model,
                       Out_full_Flag,
                       nConfounders);


      getAndPrintMargGam(f_out_marg_gam,
                 //Unique_List_Models,
                 Single_Gamma.Unique_List_Models,
                 vect_post_gam,
                 pX,
                 nConfounders);

      tmpTime = endTime;
      endTime = std::clock();
      postProcessTime = (endTime-tmpTime)/(double)(CLOCKS_PER_SEC);
      cout << "Setup Time: " << setupTime  << endl;
      if(!postProcessOnly){
        cout << "Main Loop Time: " << mainLoopTime  << endl;
      }
      cout << "Post Processing Time: " << postProcessTime  << endl;

      f_out_marg_gam.close();
      f_out_best_models.close();

      gsl_permutation_free(idx_post_gam_sort);
      gsl_vector_free(vect_post_gam);

    //  for(unsigned int line=0;line<List_Models.size();line++){
    //    List_Models[line].clear();
      for(unsigned int line=0;line<Single_Gamma.List_Models.size();line++){
        Single_Gamma.List_Models[line].clear();
      }
      //List_Models.clear();
      Single_Gamma.List_Models.clear();

      //for(unsigned int line=0;line<Unique_List_Models.size();line++){
      //  Unique_List_Models[line].clear();
      for(unsigned int line=0;line<Single_Gamma.Unique_List_Models.size();line++){
        Single_Gamma.Unique_List_Models[line].clear();
        vect_marg_log_post[line].clear();
      }

      //Unique_List_Models.clear();
      Single_Gamma.Unique_List_Models.clear();


      //n_Models_visited.clear();
      Single_Gamma.n_Models_visited.clear();

      if(!postProcessOnly){
        for(unsigned int chain=0;chain<nb_chains;chain++){
          //vect_gam[chain].clear();
          Single_Gamma.vect_gam[chain].clear();
        }
      }
      //vect_gam.clear();
      Single_Gamma.vect_gam.clear();

      gsl_matrix_free(mat_X_work2);
      gsl_matrix_free(mat_Y_work2);

      gsl_vector_free(vect_RMSE);
      if(!postProcessOnly){
    //    mat_log_marg.Free_double_matrix();
    //    mat_log_cond_post.Free_double_matrix();
          Single_Gamma.mat_log_marg.Free_double_matrix();
          Single_Gamma.mat_log_cond_post.Free_double_matrix();
      }
    //  gsl_permutation_free(MyPerm);
    //  chain_idx.clear();
    //  gsl_matrix_free(description_exchange_move);
    //  delete (My_g_AdMH);
    //  delete (My_Move_monitor);
    //  delete (My_CM);
    //  delete (My_DR);
    //  delete (t_tun);
      gsl_permutation_free(Single_Gamma.MyPerm);
      Single_Gamma.chain_idx.clear();
      gsl_matrix_free(Single_Gamma.description_exchange_move);
      delete (My_g_AdMH);
      delete (Single_Gamma.My_Move_monitor);
      delete (Single_Gamma.My_CM);
      delete (Single_Gamma.My_DR);
      delete (Single_Gamma.t_tun);

      gsl_rng_free(RandomNumberGenerator);


#if _CUDA_
  if(cudaFlag){
    cublasShutdown();
    culaShutdown();
  }
#endif
}



