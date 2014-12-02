/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modofied from ESS.cc in the ESS++ program
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
/*
#include <../Routines/struc.h>
#include <../Routines/dyn_name.h>
#include <../Routines/matrix_handling.h>
#include <../Routines/rand.h>
#include <../Routines/moves.h>
#include <../Routines/regression.h>
#include <../Routines/cond_post.h>
#include <../Routines/post_processing.h>
#include <../Routines/xml_file_read.h>
#include <../Classes/String_Matrices.h>
#include <../Classes/Int_Matrices.h>
#include <../Classes/Int_Matrices_var_dim.h>
#include <../Classes/Double_Matrices.h>
#include <../Classes/Double_Matrices_cont.h>
#include <../Classes/Prior_param.h>
#include <../Classes/Temperatures.h>
#include <../Classes/Move_monitor.h>
#include <../Classes/AdMH.h>
#include <../Classes/DR.h>
#include <../Classes/CM.h>
*/
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

using namespace std;

#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula.h>
#endif

#define DEBUG 0

Start_Up_ESS::Start_Up_ESS()
{
}



/*
************************************************************************************************
************************************************************************************************
                |================================================================|
                |                                                                |
                |      STEPWISE REGRESSION NEEDS TO BE ADAPTED FOR LOGISTIC      |
                |                                                                |
                |================================================================|
************************************************************************************************
************************************************************************************************
*/



void Start_Up_ESS::Initialise(
                            gsl_rng *RandomNumberGenerator,
                            ofstream &f_out,
                            ofstream &f_out_n_vars_in,
                            ofstream &f_out_n_models_visited,
                            ofstream &f_out_log_cond_post_per_chain,
                            ofstream &f_out_FSMH,
                            ofstream &f_out_CM,
                            ofstream &f_out_AE,
                            ofstream &f_out_DR,
                            ofstream &f_out_g,
                            ofstream &f_out_g_adapt,
                            ofstream &f_out_Gibbs,
                            ofstream &f_out_t_tun,
                            fstream &f_resume,
                            FILE *f_rng,
                            ofstream &f_out_time,
                            ofstream &f_out_best_models,
                            ofstream &f_out_marg_gam,
                            bool &postProcessOnly,
                            bool &resumeRun,
                            bool &fixInit,
                            vector < vector <unsigned int> > &Gam_step_regr,
                            bool &iso_T_Flag,
                            unsigned int &maxPGamma,
                            vector<unsigned int> &initGam,
                            unsigned int &nb_chains,
                            unsigned int &nConfounders,
                            unsigned int &pX,
                            unsigned int &nY,
                            vector< unsigned int > &shuffleYIndex,
                            double &b_t_input,
                            double &a_t_den_inf_5k,
                            double &a_t_den_5_10k,
                            double &a_t_den_sup_10k,
                            unsigned int &temp_n_batch,
                            vector <double> &M_input,
                            unsigned int &burn_in,
                            double &temp_optimal_input,
                            vector <double> &resumeT,
                            vector<unsigned int> &resumeChainIndex,
                            Prior_param &PR,
                            bool &gPriorFlag,
                            bool &indepPriorFlag,
                            bool &gSampleFlag,
                            double &lambda,
                            double &g,
                            unsigned int &nX,
                            unsigned int &pY,
                            bool &cudaFlag,
                            bool &Log_Flag,
                            double &timeLimit,
                            bool &Time_monitorFlag,
                            bool &Out_full_Flag,
                            bool &HistoryFlag,
                            double &resumeCumG,
                            double &resumeCountG,
                            unsigned int &resumeSweep,
                            vector<vector<unsigned int> > &pastModels,
                            vector<unsigned int> &pastNModelsVisited,
                            unsigned int &Gibbs_n_batch,
                            vector<vector<unsigned int> > &resumeGam,
                            double &resumeBT,
                            unsigned int &n_sweeps,
                            string &resumeName,
                            string &rngFileName,
                            string &OutputName_time,
                            unsigned int &n_top_models,
                            unsigned int &resumeDRNCalls,
                            unsigned int &resumeDRNCallsAdj,
                            vector<vector<unsigned int> > &resumeDRAccepted,
                            vector<vector<unsigned int> > &resumeDRProposed,
                            unsigned int &k_max_from_read,
                            string path_name_out,
                            long &MY_SEED,
                            bool &standardizeFlag,
                            double &g_init,
                            string filename_init,
                            bool &doYShuffle,
                            string filename_in_mat_X,
                            string filename_in_mat_Y,
                            string filename_par)
{

    resumeSweep=0;
    double resumeG;
    string Extension_out=".txt";
    string Name_number="sweeps";

    resumeName=Get_stddzed_name(path_name_out,
                                       n_sweeps,
                                       Name_number,
                                       "resume",
                                       Extension_out);

    rngFileName=Get_stddzed_name(path_name_out,
                                         n_sweeps,
                                         Name_number,
                                         "random",
                                         ".rng");

    if(MY_SEED<0){
      MY_SEED=(long)time(0);
    }
    smyrand((long)(MY_SEED),RandomNumberGenerator);

    if(resumeRun||postProcessOnly){
      // Get the current sweep
      f_resume.open(resumeName.c_str(),ios::in);
      f_resume >> resumeSweep;
      f_resume >> resumeG;
      // Need to read in the state of random number generator from file
      f_rng=fopen(rngFileName.c_str(),"rb");

      readRNG(f_rng,RandomNumberGenerator);
      fclose(f_rng);
    }

    if(postProcessOnly){
      if(resumeSweep < burn_in){
        cout << "Post processing cannot be applied because current number of sweeps less than burn in -- stopping run" << endl;
  #if _CUDA_
        if(cudaFlag){
          cublasShutdown();
          culaShutdown();
        }
  #endif
        exit(-1);
      }
    }

    //Reading X matrix.

    fstream f_X;
    //f_X.open(filename_in_mat_X,ios::in);
    f_X.open(filename_in_mat_X.c_str(),ios::in);
    f_X >> nX;
    f_X >> pX;
    /*
    cout << "nX=" << nX << endl;
    cout << "pX=" << pX << endl;
    cout << "filename_in_mat_X=" << filename_in_mat_X << endl;
    */

    gsl_matrix *mat_X=gsl_matrix_alloc(nX,pX);
    //cout << "Created gsl_matrix *mat_X" << endl;

    for(unsigned int i=0;i<nX;i++){
      for(unsigned int j=0;j<pX;j++){
        double tmp;
        f_X >> tmp;
        gsl_matrix_set(mat_X,i,j,tmp);
      }
    }
    f_X.close();

    //Reading Y matrix.
    fstream f_Y;
    //f_Y.open(filename_in_mat_Y);
    f_Y.open(filename_in_mat_Y.c_str());

    //unsigned int nY;
    //unsigned int pY;
    f_Y >> nY;
    f_Y >> pY;

    shuffleYIndex.resize(nY);

    gsl_matrix *mat_Y=gsl_matrix_alloc(nY,pY);
    if(resumeRun||postProcessOnly){
      for(unsigned int i=0;i<nY;i++){
        f_resume >> shuffleYIndex[i];
      }
    }else{
      for(unsigned int i=0;i<nY;i++){
        shuffleYIndex[i]=i;
      }
      if(doYShuffle){
          gsl_ran_shuffle(RandomNumberGenerator,&shuffleYIndex,nY,sizeof(unsigned int));
      }
    }
    for(unsigned int i=0;i<nY;i++){
      for(unsigned int j=0;j<pY;j++){
        double tmp;
        f_Y >> tmp;
        gsl_matrix_set(mat_Y,shuffleYIndex[i],j,tmp);
      }
    }
    f_Y.close();


    /////Testing the number of variables in Y
    if(nX<pY){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "There are too many outcomes ("<< pY
       << ") compared to the number of observations (" << nX
       << "), run stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }


    //////////////////////////////////
    //  Running options
    //////////////////////////////////

    if(n_sweeps==0 || burn_in==0){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "The Number of sweeps and/or the burn-in has not been specified" << endl
       << "Use -iter and/or -burn-in option(s) in the command line" << endl
       << "Run stopped" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }
    if(n_sweeps <= burn_in){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "The Number of sweeps: " << n_sweeps << " is lower than " << endl
       << "(or equal to) the burn-in: " << burn_in << " -- Run stopped" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }
     //Setting up modelling parameters.

    double E_p_gam_from_read=2;
    double Sd_p_gam_from_read=1;
    nb_chains=3;
    g=0;
    if(resumeRun){
      g=resumeG;
    }else{
      if(!gSampleFlag){
        g=g_init;
      }else{
        g=pow((double)(pX),2);
      }
    }

    double maxPGammaFactor = 10;
    double P_mutation_from_read=0.5;
    double P_sel_from_read=0.5;
    double P_csvr_r_from_read=0.375;
    double P_DR_from_read=0.5;

    //Reading the parameter File
    FILE *fparameter=NULL;
    char str[256]; // string used to read parameters

    //fparameter = fopen(filename_par,"r");
    fparameter = fopen(filename_par.c_str(),"r");
    MaXmlTagRead XML_PAR(fparameter); // assign class
    //Reading elements of the par file
    if(XML_PAR.ReadTag("E_P_GAM", 0, 0, str,256)){
      sscanf(str,"%lf",&E_p_gam_from_read);
    }
    if(XML_PAR.ReadTag("SD_P_GAM", 0, 0, str,256)){
      sscanf(str,"%lf",&Sd_p_gam_from_read);
    }
    if(XML_PAR.ReadTag("MAX_P_GAM_FACTOR", 0, 0, str,256)){
        sscanf(str,"%lf",&maxPGammaFactor);
    }
    if(XML_PAR.ReadTag("NB_CHAINS", 0, 0, str,256)){
      sscanf(str,"%u",&nb_chains);
    }
    if(XML_PAR.ReadTag("P_MUTATION", 0, 0, str,256)){
      sscanf(str,"%lf",&P_mutation_from_read);
    }
    if(XML_PAR.ReadTag("P_SEL", 0, 0, str,256)){
      sscanf(str,"%lf",&P_sel_from_read);
    }
    if(XML_PAR.ReadTag("P_CSRV_R", 0, 0, str,256)){
      sscanf(str,"%lf",&P_csvr_r_from_read);
    }
    if(XML_PAR.ReadTag("P_DR", 0, 0, str,256)){
      sscanf(str,"%lf",&P_DR_from_read);
    }


    // Add in the number of confounders
    E_p_gam_from_read = E_p_gam_from_read + nConfounders;
    // Set the limit on the maximum number of factors
    maxPGamma = (unsigned int)(E_p_gam_from_read+Sd_p_gam_from_read*maxPGammaFactor);
    if(maxPGamma>nX){
      maxPGamma=nX;
    }

    //Setting up regression parameters.
    double n_Pvalue_enter=0.01;
    double n_Pvalue_remove=0.01;

    if(XML_PAR.ReadTag("N_P_VALUE_ENTER", 0, 0, str,256)){
      sscanf(str,"%lf",&n_Pvalue_enter);
    }
    if(XML_PAR.ReadTag("N_P_VALUE_REMOVE", 0, 0, str,256)){
      sscanf(str,"%lf",&n_Pvalue_remove);
    }


    //Regression Setting up parameters.
    double Pvalue_enter = 1.0 - pow((1.0 - n_Pvalue_enter),(1.0/(double)(pX)));
    double Pvalue_remove = 1.0 - pow((1.0 - n_Pvalue_remove),(1.0/(double)(pX)));


    //Moves Parameters
    //g Adaptative M-H
    g_n_batch_from_read=100;
    g_AdMH_optimal_from_read=0.44;
    g_AdMH_ls_from_read=0.0;
    g_M_min_input=0.0;
    g_M_max_input=0.0;
    if(XML_PAR.ReadTag("G_N_BATCH", 0, 0, str,256)){
      sscanf(str,"%u",&g_n_batch_from_read);
    }
    if(XML_PAR.ReadTag("G_ADMH_OPTIMAL", 0, 0, str,256)){
      sscanf(str,"%lf",&g_AdMH_optimal_from_read);
    }
    if(XML_PAR.ReadTag("G_ADMH_LS", 0, 0, str,256)){
      sscanf(str,"%lf",&g_AdMH_ls_from_read);
    }
    if(XML_PAR.ReadTag("G_M_MIN", 0, 0, str,256)){
      sscanf(str,"%lf",&g_M_min_input);
    }
    if(XML_PAR.ReadTag("G_M_MAX", 0, 0, str,256)){
      sscanf(str,"%lf",&g_M_max_input);
    }
    //Crossover Move
    k_max_from_read=2;
    if(XML_PAR.ReadTag("K_MAX", 0, 0, str,256)){
      sscanf(str,"%u",&k_max_from_read);
    }
    //Gibbs Move
    Gibbs_n_batch=500;
    if(XML_PAR.ReadTag("GIBBS_N_BATCH", 0, 0, str,256)){
      sscanf(str,"%u",&Gibbs_n_batch);
    }

    //Temperature Placement
    b_t_input=2.0;
    a_t_den_inf_5k=2.0;
    a_t_den_5_10k=4.0;
    a_t_den_sup_10k=2.0;
    temp_n_batch=50;
    temp_optimal_input=0.5;
    M_input.resize(2);
    double M_min_input=1.0;
    double M_max_input=4.0;

    if(XML_PAR.ReadTag("B_T", 0, 0, str,256)){
      sscanf(str,"%lf",&b_t_input);
    }
    if(XML_PAR.ReadTag("A_T_DEN_INF_5K", 0, 0, str,256)){
      sscanf(str,"%lf",&a_t_den_inf_5k);
    }
    if(XML_PAR.ReadTag("A_T_DEN_5_10K", 0, 0, str,256)){
      sscanf(str,"%lf",&a_t_den_5_10k);
    }
    if(XML_PAR.ReadTag("A_T_DEN_SUP_10K", 0, 0, str,256)){
      sscanf(str,"%lf",&a_t_den_sup_10k);
    }
    if(XML_PAR.ReadTag("TEMP_N_BATCH", 0, 0, str,256)){
      sscanf(str,"%u",&temp_n_batch);
    }
    if(XML_PAR.ReadTag("TEMP_OPTIMAL", 0, 0, str,256)){
      sscanf(str,"%lf",&temp_optimal_input);
    }
    if(XML_PAR.ReadTag("M_MIN", 0, 0, str,256)){
      sscanf(str,"%lf",&M_min_input);
    }
    if(XML_PAR.ReadTag("M_MAX", 0, 0, str,256)){
      sscanf(str,"%lf",&M_max_input);
    }

    M_input[0]=M_min_input;
    M_input[1]=M_max_input;

    fclose(fparameter);



    cout << "**********************************************************" << endl
         << "***************** Setup options **********************" << endl;
    if(postProcessOnly){
      cout << "Post processing previous run only" << endl;
    }else{
      if(resumeRun){
        cout << "Resuming previous run" << endl;
      }else{
        cout << "Random seed " << MY_SEED << endl;
      }
      cout << "nb_chains " << nb_chains << endl
         << "n_sweeps " << n_sweeps << endl
         << "n_top_models " << n_top_models << endl
         << "burn_in " << burn_in << endl
         << "E_p_gam_input " << E_p_gam_from_read-nConfounders << endl
         << "Sd_p_gam_input " << Sd_p_gam_from_read << endl
         << "Max_p_gam_factor " << maxPGammaFactor << endl
         << "Max p_gam " << maxPGamma-nConfounders << endl;
      if(gSampleFlag){
        cout << "Sampling g" << endl;
      }else{
        cout << "Not sampling g" << endl;
        cout << "g " << g << endl;
      }
      if(gPriorFlag){
        cout << "Using g prior" << endl;
      }else{
        if(indepPriorFlag){
          cout << "Using independence prior" << endl;
        }else{
          cout << "Using powered prior with lambda = " << lambda << endl;
        }
      }
      cout << "CUDA " << cudaFlag << endl;
      cout << "No. confounders " << nConfounders << endl;
    }
    cout << "**********************************************************" << endl
         << "**********************************************************" << endl << endl;

    if(!postProcessOnly&&!resumeRun&&!fixInit){
      cout << "**********************************************************" << endl
         << "*************** Regression parameters ********************" << endl
         << "n_Pvalue_enter " << n_Pvalue_enter << endl
         << "n_Pvalue_enter " << n_Pvalue_enter << endl
         << "Pvalue_enter stepwise " << Pvalue_enter << endl
         << "Pvalue_remove stepwise " << Pvalue_remove << endl
         << "**********************************************************" << endl
         << "**********************************************************" << endl << endl;
    }

    //Testing PATH-file names for output

    Extension_out=".txt";
    Name_number="sweeps";
    string Main_Output_name="output_models_history";
    string OutputName=Get_stddzed_name(path_name_out,
                       n_sweeps,
                       Name_number,
                       Main_Output_name,
                       Extension_out);
    fstream f_in;

    string OutputName_n_vars_in=Get_stddzed_name(path_name_out,
                             n_sweeps,
                             Name_number,
                             "output_model_size_history",
                             Extension_out);

    string OutputName_n_models_visited=Get_stddzed_name(path_name_out,
                                                 n_sweeps,
                                                 Name_number,
                                                 "output_n_models_visited_history",
                                                 Extension_out);
    fstream f_in_n_models_visited;

    string OutputName_log_cond_post_per_chain=Get_stddzed_name(path_name_out,
                                   n_sweeps,
                                   Name_number,
                                   "output_log_cond_post_prob_history",
                                   Extension_out);

    string OutputName_FSMH=Get_stddzed_name(path_name_out,
                                            n_sweeps,
                                            Name_number,
                                            "output_fast_scan_history",
                                            Extension_out);

    string OutputName_CM=Get_stddzed_name(path_name_out,
                                          n_sweeps,
                                          Name_number,
                                          "output_cross_over_history",
                                          Extension_out);

    string OutputName_AE=Get_stddzed_name(path_name_out,
                                            n_sweeps,
                                            Name_number,
                                            "output_all_exchange_history",
                                            Extension_out);

    string OutputName_DR=Get_stddzed_name(path_name_out,
                                           n_sweeps,
                                           Name_number,
                                           "output_delayed_rejection_history",
                                           Extension_out);

    string OutputName_g=Get_stddzed_name(path_name_out,
                                          n_sweeps,
                                          Name_number,
                                          "output_g_history",
                                          Extension_out);

    string OutputName_g_adapt=Get_stddzed_name(path_name_out,
                                                n_sweeps,
                                                Name_number,
                                                "output_g_adaptation_history",
                                                Extension_out);

    string OutputName_Gibbs=Get_stddzed_name(path_name_out,
                                              n_sweeps,
                                              Name_number,
                                              "output_gibbs_history",
                                              Extension_out);

    string OutputName_t_tun=Get_stddzed_name(path_name_out,
                                              n_sweeps,
                                              Name_number,
                                              "output_temperature_history",
                                              Extension_out);

    ios_base::openmode fileMode;
    if(resumeRun){
      fileMode=ios::app;
    }else{
      fileMode=ios::out;
    }


    if(HistoryFlag){

    /*
    #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

       if true: Read files from previous run (resume)

    #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
    */
      if(resumeRun||postProcessOnly){
        pastModels.resize(resumeSweep+1);
        f_in.open(OutputName.c_str(),ios::in);
        if(f_in.fail()){
          cout << "Trying to resume a run where no history was written -- stopping run" << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        // Remove the first line
        string strtmp;
        for(unsigned int i=0;i<5;i++){
          f_in >> strtmp;
        }
        for(unsigned int i=0;i<resumeSweep+1;i++){
          unsigned int tmp,modelSize;
          double tmp1;
          f_in >> tmp;
          f_in >> modelSize;
          f_in >> tmp1;
          f_in >> tmp1;
          pastModels[i].resize(modelSize+1);
          pastModels[i][0]=modelSize;
          for(unsigned int j=0;j<modelSize;j++){
            f_in >> tmp;
            pastModels[i][j+1]=tmp-1;
          }
        }
        f_in.close();
        if(resumeRun){
          f_out.open(OutputName.c_str(),fileMode);
        }
      }else{
        f_out.open(OutputName.c_str(),fileMode);
      }
  /*
  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

     Various initialisations for the Output

  #=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
  */
      if(!postProcessOnly){
        if(f_out.fail()){
          cout << "Invalid Path and/or permission rights for " << OutputName << " -- run stopped." << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        else{
          if(!resumeRun){
            f_out << "Sweep\tModel_size\tlog_marg\tlog_cond_post\tModel"<<endl;
          }
        }
      }

      if(!postProcessOnly){
        f_out_n_vars_in.open(OutputName_n_vars_in.c_str(),fileMode);
        if(f_out_n_vars_in.fail()){
          cout << "Invalid Path and/or permission rights for " << OutputName_n_vars_in << " -- run stopped." << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        else{
          if(!resumeRun){
            f_out_n_vars_in << "Sweep\t";
            for(unsigned int tmp_chain=0;tmp_chain<nb_chains;tmp_chain++){
              f_out_n_vars_in << "Chain_"<< tmp_chain+1 << "\t";
            }
            f_out_n_vars_in << endl;
          }
        }
      }

      if(resumeRun||postProcessOnly){
        pastNModelsVisited.assign(resumeSweep+1,0);
        f_in_n_models_visited.open(OutputName_n_models_visited.c_str(),ios::in);
        if(f_in_n_models_visited.fail()){
          cout << "Trying to resume a run where no history was written -- stopping run" << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        string tmpstr;
        f_in_n_models_visited >> tmpstr;
        f_in_n_models_visited >> tmpstr;
        for(unsigned int i=1;i<resumeSweep+1;i++){
          unsigned int tmp;
          f_in_n_models_visited >> tmp;
          f_in_n_models_visited >> pastNModelsVisited[i];
        }
        f_in_n_models_visited.close();
        if(resumeRun){
          f_out_n_models_visited.open(OutputName_n_models_visited.c_str(),fileMode);
        }
      }else{
        f_out_n_models_visited.open(OutputName_n_models_visited.c_str(),fileMode);
      }
      if(!postProcessOnly){
        if(f_out_n_models_visited.fail()){
          cout << "Invalid Path and/or permission rights for " << OutputName_n_models_visited << " -- run stopped." << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        else{
          if(!resumeRun){
            f_out_n_models_visited << "Sweep\tn_models_visited" << endl;
          }
        }
      }

      if(!postProcessOnly){
        f_out_log_cond_post_per_chain.open(OutputName_log_cond_post_per_chain.c_str(),fileMode);
        if(f_out_log_cond_post_per_chain.fail()){
          cout << "Invalid Path and/or permission rights for " << OutputName_log_cond_post_per_chain << " -- run stopped." << endl;
  #if _CUDA_
          if(cudaFlag){
            cublasShutdown();
            culaShutdown();
          }
  #endif
          exit(-1);
        }
        else{
          if(!resumeRun){
            f_out_log_cond_post_per_chain << "Sweep"<< "\t";
            for(unsigned int tmp_chain=0;tmp_chain<nb_chains;tmp_chain++){
              f_out_log_cond_post_per_chain << "Chain_"<< tmp_chain+1 << "\t";
            }
            f_out_log_cond_post_per_chain << endl;
          }
        }
      }

      if(!postProcessOnly){
        f_out_Gibbs.open(OutputName_Gibbs.c_str(),fileMode);
        f_out_FSMH.open(OutputName_FSMH.c_str(),fileMode);
        f_out_CM.open(OutputName_CM.c_str(),fileMode);
        f_out_AE.open(OutputName_AE.c_str(),fileMode);
        f_out_DR.open(OutputName_DR.c_str(),fileMode);
        if(gSampleFlag){
          f_out_g.open(OutputName_g.c_str(),fileMode);
          f_out_g_adapt.open(OutputName_g_adapt.c_str(),fileMode);
        }
        if(iso_T_Flag==false){
          f_out_t_tun.open(OutputName_t_tun.c_str(),fileMode);
        }
      }
    }else{
      if(resumeRun||postProcessOnly){
        cout << "Trying to resume a run where no history was written -- stopping run" << endl;
  #if _CUDA_
        if(cudaFlag){
          cublasShutdown();
          culaShutdown();
        }
  #endif
        exit(-1);

      }
    }

    if(!postProcessOnly){
      if(Time_monitorFlag){
        OutputName_time=Get_stddzed_name(path_name_out,
                          n_sweeps,
                          Name_number,
                          "output_time_monitor",
                          Extension_out);
        f_out_time.open(OutputName_time.c_str(),fileMode);
        if(!resumeRun){
          f_out_time << "Sweep\tTime\tTime_per_eval_model"<<endl;
        }
      }
    }

    string OutputName_best_models=Get_stddzed_name(path_name_out,
                           n_sweeps,
                           Name_number,
                           "output_best_visited_models",
                           Extension_out);

    f_out_best_models.open(OutputName_best_models.c_str(),ios::out);
    if(f_out_best_models.fail()){
      cout << "Invalid Path and/or permission rights for " << OutputName_best_models << " -- run stopped." << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }
    else{
      if(Out_full_Flag){
        f_out_best_models << "Rank\t#Visits\tSweep_1st_visit\t#models_eval_before_1st_visit\tModel_size\tlog_Post_Prob\tModel_Post_Prob\tJeffreys_scale\tModel"<<endl;
      }
      else{
        f_out_best_models << "Rank\t#Visits\tModel_size\tlog_Post_Prob\tModel_Post_Prob\tJeffreys_scale\tModel"<<endl;
      }
    }
    string OutputName_marg_gam=Get_stddzed_name(path_name_out,
                           n_sweeps,
                           Name_number,
                           "output_marg_prob_incl",
                           Extension_out);

    f_out_marg_gam.open(OutputName_marg_gam.c_str(),ios::out);
    if(f_out_marg_gam.fail()){
      cout << "Invalid Path and/or permission rights for " << OutputName_marg_gam << " -- run stopped." << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }
    else{
      f_out_marg_gam << "Predictor\tMarg_Prob_Incl"<< endl;
    }


/*
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

        Again: Read files from previous run (resume)

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
*/

    vect_RMSE = gsl_vector_calloc(pY);

    double resumeLS;

    if(resumeRun||postProcessOnly){
      // Read the last state in from file
      f_resume >> resumeCumG;
      f_resume >> resumeCountG;
      for(unsigned int j=0;j<pY;j++){
        f_resume >> vect_RMSE->data[j];
      }
      resumeT.resize(nb_chains);
      for(unsigned int j=0;j<nb_chains;j++){
        f_resume >> resumeT[j];
      }
      f_resume >> resumeBT;
      f_resume >> resumeLS;
      resumeDRAccepted.resize(nb_chains);
      resumeDRProposed.resize(nb_chains);
      f_resume >> resumeDRNCalls;
      f_resume >> resumeDRNCallsAdj;
      for(unsigned int j=0;j<nb_chains;j++){
        resumeDRProposed[j].resize(nb_chains);
        for(unsigned int k=0;k<nb_chains;k++){
          f_resume >> resumeDRProposed[j][k];
        }
      }
      for(unsigned int j=0;j<nb_chains;j++){
        resumeDRAccepted[j].resize(nb_chains);
        for(unsigned int k=0;k<nb_chains;k++){
          f_resume >> resumeDRAccepted[j][k];
        }
      }
      resumeGam.resize(nb_chains);
      resumeChainIndex.resize(nb_chains);
      for(unsigned int j=0;j<nb_chains;j++){
        unsigned int gamSize;
        f_resume >> gamSize;
        f_resume >> resumeChainIndex[j];
        resumeGam[resumeChainIndex[j]].resize(gamSize);
        for(unsigned int k=0;k<gamSize;k++){
          f_resume >> resumeGam[resumeChainIndex[j]][k];
        }
      }
      f_resume.close();
    }

/*
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

        New Run, but initialise Gamma by reading from a file

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
*/

    if((!resumeRun)&&(!postProcessOnly))
    {
        standardize_matrix_gsl(mat_X);
    }

    if((!resumeRun)&&(!postProcessOnly))
    {
        if(fixInit){
        // Read the initial gamma from file
        ifstream inputFile;
        //inputFile.open(filename_init);
        inputFile.open(filename_init.c_str());
        if(!inputFile.is_open()){
          cout << "Input file not found" << endl;
          exit(-1);
        }
        unsigned int initSize;
        inputFile >> initSize;
        initGam.resize(initSize);
        for(unsigned int i=0;i<initSize;i++){
          inputFile >> initGam[i];
        }
        inputFile.close();

        cout << "********************************************************" << endl
            << "***************** Initial indices entered **************" << endl;
        for(unsigned int j=0;j<initGam.size();j++){
          cout << initGam[j] << " ";
        }
        cout << endl;
        cout << "********************************************************" << endl
            << "********************************************************" << endl;


        // Need to get the RMSE for prior estimate of k

        gsl_vector *current_outcome =gsl_vector_calloc(nX);

        vector < unsigned int > list_columns_X_gam;
        vector < unsigned int > is_var_in(pX);

        unsigned int k=0;
        for(unsigned int j=0;j<pX;j++){
          // Set up so confounders are always in
          if(j<nConfounders){
            is_var_in[j]=1;
          }else if(j==initGam[k]){
            k++;
            is_var_in[j]=1;
          }else{
            is_var_in[j]=0;
          }
        }
        get_list_var_in(list_columns_X_gam,is_var_in);

        gsl_matrix* mat_X_gam=get_X_reduced_and_constant(list_columns_X_gam,
                                                              mat_X);

        double tolerance= 6.6835e-14;

        for(unsigned int outcome=0;outcome<pY;outcome++){
          gsl_matrix_get_col(current_outcome,mat_Y,outcome);
          getEstimateRMSE(mat_X_gam,
                        current_outcome,
                        outcome,
                        vect_RMSE,
                        tolerance);
        }

        gsl_matrix_free(mat_X_gam);
        gsl_vector_free(current_outcome);
      }
    }

/*
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

        New Run, Stepwise Regression

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
*/

/*
************************************************************************************************
************************************************************************************************
                |================================================================|
                |                                                                |
                |      STEPWISE REGRESSION NEEDS TO BE ADAPTED FOR LOGISTIC      |
                |                                                                |
                |================================================================|
************************************************************************************************
************************************************************************************************
*/



    if((!resumeRun)&&(!postProcessOnly))
    {

      if(!fixInit){

    /*
    ************************************************************************************************
    ************************************************************************************************
        Stepwise regression: specific to each response
    ************************************************************************************************
    ************************************************************************************************
    */
        Double_Matrices Gam_step_regr_pvals;
        Gam_step_regr_pvals.Alloc_double_matrix(pY,
                        pX);
        Double_Matrices Gam_step_regr_SE;
        Gam_step_regr_SE.Alloc_double_matrix(pY,
                         pX);

        Double_Matrices Gam_step_regr_beta;
        Gam_step_regr_beta.Alloc_double_matrix(pY,
                       pX);

        double tolerance= 6.6835e-14;

        gsl_vector *current_outcome =gsl_vector_calloc(nX);
        gsl_vector *vect_residuals=gsl_vector_calloc(nX);
        gsl_vector *vect_p_value=gsl_vector_calloc(pX);
        gsl_vector *vect_beta_full=gsl_vector_calloc(pX);
        gsl_vector *vect_SE_full=gsl_vector_calloc(pX);

        cout << "***************************************************" << endl
            << "*************  Stepwise regression   *************" << endl
            << "***************************************************" << endl << endl;

        for(unsigned int outcome=0;outcome<pY;outcome++)
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
          Gam_step_regr.resize(pY);
          store_model_per_outcome(Gam_step_regr,
                  list_columns_X_gam,
                  vect_p_value,
                  vect_beta_full,
                  vect_SE_full,
                  Gam_step_regr_pvals,
                  Gam_step_regr_SE,
                  Gam_step_regr_beta,
                  outcome);


          list_columns_X_gam.clear();
          is_var_in.clear();
        }

        cout << "Result From Step-Wise Regression" << endl;
        display_matrix_var_dim(Gam_step_regr);

        cout << endl;


        cout << "***************************************************" << endl
            << "**********  End of Stepwise regression  **********" << endl
            << "***************************************************" << endl << endl;

        gsl_vector_free(vect_residuals);
        gsl_vector_free(current_outcome);
        gsl_vector_free(vect_p_value);
        gsl_vector_free(vect_beta_full);
        gsl_vector_free(vect_SE_full);

        Gam_step_regr_pvals.Free_double_matrix();
        Gam_step_regr_SE.Free_double_matrix();
        Gam_step_regr_beta.Free_double_matrix();

      }
    }


    if((!resumeRun)&&(!postProcessOnly))
    {
        gsl_matrix_free(mat_X);
    }



    // Re-read in the data (we need to do this to remove the standardisation)
    //f_X.open(filename_in_mat_X,ios::in);
    f_X.open(filename_in_mat_X.c_str(),ios::in);
    f_X >> nX;
    f_X >> pX;

    mat_X_work2=gsl_matrix_alloc(nX,pX);

    for(unsigned int i=0;i<nX;i++){
      for(unsigned int j=0;j<pX;j++){
        double tmp;
        f_X >> tmp;
        gsl_matrix_set(mat_X_work2,i,j,tmp);
      }
    }
    f_X.close();

    mat_Y_work2=mat_Y;

    //Centering X and Y
    center_matrix_gsl(mat_X_work2);
    center_matrix_gsl(mat_Y_work2);
    if(standardizeFlag){
      // Standardize X matrix if selected
      standardize_matrix_gsl(mat_X_work2);
    }


    if(!postProcessOnly){
      cout << "**********************************************************" << endl
         << "****************** MOVES parameters **********************" << endl
         << "g-adaptative M-H" << endl
         << "\t-g_n_batch: " << g_n_batch_from_read << endl
         << "\t-g_AdMH_optimal: " << g_AdMH_optimal_from_read << endl
         << "\t-g_AdMH_ls: " << g_AdMH_ls_from_read << endl
         << "Crossover Move" << endl
         << "\tk_max: " <<k_max_from_read << endl
         << "Gibbs Move" << endl
         << "\tGibbs_n_batch: " <<Gibbs_n_batch << endl
         << "**********************************************************" << endl
         << "**********************************************************" << endl << endl;

      cout << "**********************************************************" << endl
          << "****************** TEMP parameters **********************" << endl
          << "b_t " << b_t_input << endl
          << "a_t_den_inf_5k " << a_t_den_inf_5k << endl
          << "a_t_den_5_10k " << a_t_den_5_10k << endl
          << "a_t_den_sup_10k " << a_t_den_sup_10k << endl
          << "temp_n_batch " << temp_n_batch << endl
          << "temp_optimal " << temp_optimal_input << endl
          << " M= [" << M_input[0] << " - " << M_input[1] << "]" << endl
          << "**********************************************************" << endl
          << "**********************************************************" << endl << endl;
    }

    //////////////////////////////////
    //  Setting up prior parameters
    //////////////////////////////////

    PR.set_PR_param(E_p_gam_from_read,
            Sd_p_gam_from_read,
            pX,
            pY,
            nX,
            lambda,
            vect_RMSE,
            P_mutation_from_read,
            P_sel_from_read,
            P_csvr_r_from_read,
            P_DR_from_read);

    if(!postProcessOnly){
      PR.display_prior_param();
    }

    //////////////////////////////////
    //  Setting up AdMH parameters
    //////////////////////////////////

    My_g_AdMH=new AdMH;
    if(!postProcessOnly){
      (*My_g_AdMH).set_AdMH(gSampleFlag,
                g_n_batch_from_read,
                g_AdMH_optimal_from_read,
                g_AdMH_ls_from_read,
                pX,
                burn_in,
                g_M_min_input,
                g_M_max_input);

      if(resumeRun){
        (*My_g_AdMH).ls=resumeLS;
      }
      (*My_g_AdMH).display_AdMH();
    }

}
