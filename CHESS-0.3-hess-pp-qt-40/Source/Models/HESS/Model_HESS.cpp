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

#include "../../General_Classes/Kernel_Single_Gamma.h"
#include "../../General_Classes/Command_Line.h"
#include "../../General_Classes/Model_Generic.h"
#include "Model_HESS.h"
#include "../../General_Classes/Model_Information.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <iomanip>

#include <ctime>

#include <vector>
#include <algorithm>
#include <memory>
#include <fstream>
#include <sstream>

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
//#define Missing_data -1
#define Print 1

#define IGNORE 1

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF


using namespace std;
#define Fixed_Random 0

#if Fixed_Random

    //string Path_Ext_Random_Nbrs=Path_Input+"/Ext_Random_Nbrs";
    string Path_Ext_Random_Nbrs=" ";
    fstream Random_Nbrs_File;
    vector < vector < double > > Ext_Random_Nbrs;

    unsigned int Line_Ext_Rnd_Nmbr=0;

#endif

Model_HESS::Model_HESS()
{

    Model_Tag="HESS";
    Regression_Model="Undefined";

}


void Model_HESS::Run()
{



    //cout << "This is version 35 of HESS." << endl;
    //cout << "\n\nWe test the implementation of one g per response:" << endl;

#if Fixed_Random

  Path_Ext_Random_Nbrs=(string)path_name_out+(string)"_Ext_Random_Nbrs";
  cout << "Random numbers from external file " << Path_Ext_Random_Nbrs << endl;

  // Input: Random Number
  Random_Nbrs_File.open(Path_Ext_Random_Nbrs.c_str(),ios::in);
  if(!Random_Nbrs_File.is_open())
  {
    cout << "Input file for Random_Nbrs_File not found" << endl;
    exit(-1);
  }

  Read_File_Double(Random_Nbrs_File, Ext_Random_Nbrs);
  Random_Nbrs_File.close();


#endif

    cout << "Single g = " << Single_g << endl;

    bool DEBUG_sample_g=false;
    if (DEBUG_sample_g==true)
    {cout << "in Model_HESS::Run:  DEBUG_sample_g = " << DEBUG_sample_g << endl;}


    /*
      ##############################################################################

                                     Verify Output and GPU

     ##############################################################################
   */

    //Define output

    ios_base::openmode fileMode;
     //       if(resumeRun){
     //         fileMode=ios::app;
     //       }else{
     //         fileMode=ios::out;
     //       }

    fileMode=ios::out;

    cout << endl;




    bool WRITE_BY_1_LINE = true;

    string Field_Separator="\t";


    string Output_Name_Av_Rho_j=(string)path_name_out+(string)"_Av_Rho_j.txt";
    ofstream Output_File_Av_Rho_j;
    Test_Output_File(Output_Name_Av_Rho_j,Output_File_Av_Rho_j,fileMode);

    string Output_Name_Av_Omega_k=(string)path_name_out+(string)"_Av_Omega_k.txt";
    ofstream Output_File_Av_Omega_k;
    Test_Output_File(Output_Name_Av_Omega_k,Output_File_Av_Omega_k,fileMode);

    string Output_Name_History_Rho=(string)path_name_out+(string)"_History_Rho_j.txt";
    ofstream Output_File_History_Rho;
    Test_Output_File(Output_Name_History_Rho,Output_File_History_Rho,fileMode);

    string Output_Name_History_Omega=(string)path_name_out+(string)"_History_Omega_k.txt";
    ofstream Output_File_History_Omega;
    Test_Output_File(Output_Name_History_Omega,Output_File_History_Omega,fileMode);

    string Output_Name_Av_Omega=(string)path_name_out+(string)"_Av_Matrix_Omega_kj.txt";
    ofstream Output_File_Av_Omega;
    Test_Output_File(Output_Name_Av_Omega,Output_File_Av_Omega,fileMode);

    string Output_Name_History_g=(string)path_name_out+(string)"_History_g.txt";
    ofstream Output_File_History_g;
    Test_Output_File(Output_Name_History_g,Output_File_History_g,fileMode);

    string Output_Name_Counts_Gamma_kj=(string)path_name_out+(string)"_Counts_Gamma_kj.txt";
    ofstream Output_File_Counts_Gamma_kj;
    Test_Output_File(Output_Name_Counts_Gamma_kj,Output_File_Counts_Gamma_kj,fileMode);

    // FILES FOR THE POST-PROCESSING

    string Output_Name_TP_Rho_j=(string)path_name_out+(string)"_Tail_Prob_Rho_j.txt";
    ofstream Output_File_TP_Rho_j;
    Test_Output_File(Output_Name_TP_Rho_j,Output_File_TP_Rho_j,fileMode);

    string Output_Name_Best_Models_Marg=(string)path_name_out+(string)"_Best_Models_List.txt";
    ofstream Output_File_Best_Models_Marg;
    Test_Output_File(Output_Name_Best_Models_Marg,Output_File_Best_Models_Marg,fileMode);

    string Output_Name_Prob_Best_Models_Marg=(string)path_name_out+(string)"_Best_Models_Post_Prob.txt";
    ofstream Output_File_Prob_Best_Models_Marg;
    Test_Output_File(Output_Name_Prob_Best_Models_Marg,Output_File_Prob_Best_Models_Marg,fileMode);

    string Output_Name_Marg_Prob_Incl=(string)path_name_out+(string)"_Matrix_Marg_Prob_Incl.txt";
    ofstream Output_File_Marg_Prob_Incl;
    Test_Output_File(Output_Name_Marg_Prob_Incl,Output_File_Marg_Prob_Incl,fileMode);

    string Output_Name_Models_Prob=(string)path_name_out+(string)"_List_Models_Post_Prob.txt";
    ofstream Output_File_Models_Prob;
    Test_Output_File(Output_Name_Models_Prob,Output_File_Models_Prob,fileMode);



    if (Print==1)
    {
     cout << endl;
    }

    /*
      ##############################################################################

                               Define general objects

     ##############################################################################
   */


    Model_Information HESS_Information(Model_Tag, Regression_Model);

    std::clock_t startTime = std::clock();
    std::clock_t endTime;
    std::clock_t Loop_Start_Time;
    std::clock_t PostProc_Start_Time;


    gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_mt19937 );

    if(MY_SEED<0){
      MY_SEED=(long)time(0);
    }
    smyrand((long)(MY_SEED),RandomNumberGenerator);



    /*
      ##############################################################################

                                    DATA

      ##############################################################################
    */
    unsigned int nX=0;
    unsigned int pX=0;
    unsigned int nY=0;
    unsigned int Nb_Resp=0;
    unsigned int q_total_Nb_Yvar=0;


    gsl_matrix * Data_X;
    // Vectorise Data_Y
    vector <gsl_matrix *> Vector_Data_Y;

    fstream f_X;
    fstream f_Y;

    // Read Data X:

    cout << "Reading Data_X..." << endl;
    f_X.open(filename_in_mat_X.c_str(),ios::in);
    Data_X=Read_Data_X(f_X,nX,pX);
    f_X.close();

    // Read Data Y
    cout << "Reading Data_Y..." << endl;
    f_Y.open(filename_in_mat_Y.c_str(),ios::in);
    Vector_Data_Y=Read_Data_Y(f_Y,nY,q_total_Nb_Yvar,Nb_Resp);
    f_Y.close();

    //Centering X and Y

    //cout << "The centering seems to result in a loss of precision at the 5th digit" << endl;
    center_matrix_gsl(Data_X);
    for (unsigned int Response=0;Response<Nb_Resp;Response++)
    {
        center_matrix_gsl(Vector_Data_Y[Response]);
    }

    if(standardizeFlag)
    {
      // Standardize X matrix if selected
      //cout << "standardizeFlag = " << standardizeFlag << endl;
      standardize_matrix_gsl(Data_X);
    }

    /*
      ##############################################################################

                                     Test validity of data size

     ##############################################################################
   */



    Preliminary_Tests(nX,
                      nY,
                      Nb_Resp,
                      n_sweeps,
                      burn_in,
                      gamSampleFlag,
                      fixInit,
                      OmegaKSampleFlag,
                      fixInit_omega_k,
                      RhoJSampleFlag,
                      fixInit_rho_j);


    /*
      ##############################################################################

                            Kernel Parameters

      ##############################################################################
    */

    //fixInit=false;
    unsigned int maxPGamma=0;
    unsigned int nb_chains=0;

    // Now we need one PR object per response
    //Prior_param PR
    vector<Prior_param> PR_per_Resp;
        PR_per_Resp.resize(Nb_Resp);
    Prior_param Omega_k_PR;

    unsigned int Gibbs_n_batch=0;

    double b_t_input=0;
    double a_t_den_inf_5k=0;
    double a_t_den_5_10k=0;
    double a_t_den_sup_10k=0;

    unsigned int temp_n_batch=0;
    vector <double> M_input;
    double temp_optimal_input=0;

    unsigned int k_max_from_read=0;

    //Stepwise Regression
    double Pvalue_enter=0;
    double Pvalue_remove=0;

    //Adaptation
    unsigned int g_n_batch_from_read=0;
    double g_AdMH_optimal_from_read=0;
    double g_AdMH_ls_from_read=0;
    double g_M_min_input=0;
    double g_M_max_input=0;

    //Temporary variables
    double E_p_gam_from_read=2;
    double Sd_p_gam_from_read=1;
    double P_mutation_from_read=0.5;
    double P_sel_from_read=0.5;
    double P_csvr_r_from_read=0.375;
    double P_DR_from_read=0.5;


    // Default value for g_init
    // Used only if we sample g and no initial value is provided
    if(gSampleFlag && !gInitFlag)
    {
      //g_init=pow((double)(pX),2);
      g_init=nX;
    }


    n_top_models=n_sweeps; // This parameter is not used, only for display.
    // Read the xml file and setup the corresponding parameters
    Setup_Model_Parameters_HESS(postProcessOnly,
                                resumeRun,
                                fixInit,
                                maxPGamma,
                                nb_chains,
                                nConfounders,
                                pX,
                                b_t_input,
                                a_t_den_inf_5k,
                                a_t_den_5_10k,
                                a_t_den_sup_10k,
                                temp_n_batch,
                                M_input,
                                burn_in,
                                temp_optimal_input,
                                gPriorFlag,
                                indepPriorFlag,
                                gSampleFlag,
                                lambda,
                                nX,
                                cudaFlag,
                                Gibbs_n_batch,
                                n_sweeps,
                                n_top_models,
                                k_max_from_read,
                                MY_SEED,
                                g_init,
                                filename_par,
                                Pvalue_enter,
                                Pvalue_remove,
                                g_n_batch_from_read,
                                g_AdMH_optimal_from_read,
                                g_AdMH_ls_from_read,
                                g_M_min_input,
                                g_M_max_input,
                                E_p_gam_from_read,
                                Sd_p_gam_from_read,
                                P_mutation_from_read,
                                P_sel_from_read,
                                P_csvr_r_from_read,
                                P_DR_from_read);

    //cout << "Model parameters fixed." << endl << endl;



    vector< vector< unsigned int > > initGamMat;
    vector < vector<vector<unsigned int> > > Gam_step_regr_Vect;


    /*
        Parametres NOT used for the moment (concern resume option)
        ----------------------------------------------------------
    */
    unsigned int resumeDRNCalls=0;
    unsigned int resumeDRNCallsAdj=0;
    vector<vector<unsigned int> > resumeDRAccepted;
    vector<vector<unsigned int> > resumeDRProposed;
    vector < vector<vector<unsigned int> > > resumeGamMat;
    vector <double> resumeT;

    double resumeBT=0;
    unsigned int resumeSweep;

   /*
     ##############################################################################

                                    Define MCMC state

     ##############################################################################
   */

    // Old code:
    //double g=g_init;

    vector <double> g;
    if (Single_g==true)
        g.resize(1);
    else
        g.resize(Nb_Resp);


    for(unsigned int k=0;k<g.size();k++)
    {
        //cout << "g_init = " << g_init << endl;
        g[k]=g_init;

    }


    // Adaptation for g
    // Old code:
    /*
    AdMH * My_g_AdMH=new AdMH();
    (*My_g_AdMH).set_AdMH(gSampleFlag,
              g_n_batch_from_read,
              g_AdMH_optimal_from_read,
              g_AdMH_ls_from_read,
              pX,
              burn_in,
              g_M_min_input,
              g_M_max_input);
    */
    vector < AdMH * > My_g_AdMH;
    if (Single_g==true)
        My_g_AdMH.resize(1);
    else
        My_g_AdMH.resize(Nb_Resp);

    for(unsigned int k=0;k<My_g_AdMH.size();k++)
    {
        My_g_AdMH[k]=new AdMH();

        (*My_g_AdMH[k]).set_AdMH(gSampleFlag,
                  g_n_batch_from_read,
                  g_AdMH_optimal_from_read,
                  g_AdMH_ls_from_read,
                  pX,
                  burn_in,
                  g_M_min_input,
                  g_M_max_input);
    }



    // Adaptation for Omega_k, need a matrix of AdMH objects

    vector< vector < AdMH *> > Omega_k_AdMH;
    Omega_k_AdMH.resize(Nb_Resp);
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        Omega_k_AdMH[k].resize(nb_chains);
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {
            Omega_k_AdMH[k][Chain]=new AdMH();
        }
    }

    // We prefer not to use the method set_AdMH,
    // it doesn't do the same thing as in the Matlab code
    Set_Omega_k_AdMH(Omega_k_AdMH,
                     Nb_Resp,
                     nb_chains,
                     burn_in,
                     g_AdMH_ls_from_read,
                     g_n_batch_from_read,
                     g_AdMH_optimal_from_read);

    // Adaptation for Rho_j, need a matrix of AdMH objects

    vector< vector < AdMH *> > Rho_j_AdMH;
    Rho_j_AdMH.resize(pX);
    for (unsigned int j=0;j<pX;j++)
    {
        Rho_j_AdMH[j].resize(nb_chains);
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {
            Rho_j_AdMH[j][Chain]=new AdMH();
        }
    }

    // We prefer not to use the method set_AdMH,
    // it doesn't do the same thing as in the Matlab code
    Set_Rho_j_AdMH(Rho_j_AdMH,
                   pX,
                   nb_chains,
                   burn_in,
                   g_AdMH_ls_from_read,
                   g_n_batch_from_read,
                   g_AdMH_optimal_from_read);


    // Kernel_Single_Gamma:
    //#####################

    //cout << "Right before vector of Kernel_Single_Gamma objects creation" << endl;

    vector<Kernel_Single_Gamma > Gammas;
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        //cout << "Response = " << k+1 << endl;
        Gammas.push_back(Kernel_Single_Gamma(postProcessOnly,
                                             nb_chains,
                                             resumeRun,
                                             resumeDRNCallsAdj,
                                             resumeDRNCalls,
                                             resumeDRAccepted,
                                             resumeDRProposed,
                                             k_max_from_read,
                                             pX,
                                             n_sweeps,
                                             HESS_Information));

    }


    /*
      VERY IMPORTANT:
      ==============

        The position of the chain c for response k is Gammas[k]->chain_idx[c]
        --> This is important if we want to access gamma
        --> in the class Kernel_Single_Gamma, chain_idx and vect_gam should be made private, and
             a method should give access to the correct line of vect_gam.
    */


    vector < vector < double > > Omega_PerLine;
    vector < vector < double > > Rho_PerCol;
    vector < unsigned int > Active_Yk; //Subset of responses to be visited

    Omega_PerLine.resize(nb_chains);

    Rho_PerCol.resize(nb_chains);

    HESS_Information.Omega_PerLine=&Omega_PerLine;
    HESS_Information.Rho_PerCol=&Rho_PerCol;
    HESS_Information.Active_Yk=&Active_Yk;


/*
     ##############################################################################

                       Creation of the History and Output

     ##############################################################################
*/

    //We keep only the output of the first chain

    vector < vector < double > > History_Rho;
    History_Rho.resize(n_sweeps+1);
    for (unsigned int sweep=0;sweep<n_sweeps+1;sweep++)
    {
        History_Rho[sweep].resize(pX);
    }

    vector < vector < double > > History_Omega;
    History_Omega.resize(n_sweeps+1);
    for (unsigned int sweep=0;sweep<n_sweeps+1;sweep++)
    {
        History_Omega[sweep].resize(Nb_Resp);
    }

    vector < vector < double > > Average_Omega_k_j;
    Average_Omega_k_j.resize(Nb_Resp);
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        Average_Omega_k_j[k].resize(pX);
        for (unsigned int j=0;j<pX;j++)
        {
            Average_Omega_k_j[k][j]=0;
        }
    }

    vector < vector < double > > History_g;
    History_g.resize(n_sweeps+1);
    for (unsigned int sweep=0;sweep<n_sweeps+1;sweep++)
    {
        //History_g[sweep].resize(1);
        History_g[sweep].resize(g.size());
    }

    vector < vector < unsigned int > > Count_Visits_Gamma_kj;
    Count_Visits_Gamma_kj.resize(Nb_Resp);
    for(unsigned int k=0;k<Nb_Resp;k++)
    {
      Count_Visits_Gamma_kj[k].resize(pX);
      for(unsigned int j=0;j<pX;j++)
      {
        Count_Visits_Gamma_kj[k][j]=0;
      }
    }

/*
     ##############################################################################

                            Initialisations

     ##############################################################################
*/


    // Auxiliary object to set priors later:
    vector <gsl_vector *>  vect_vect_RMSE;
    vect_vect_RMSE.resize(Nb_Resp);
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        unsigned int n_outcomes=Vector_Data_Y[k]->size2;
        vect_vect_RMSE[k]=gsl_vector_calloc(n_outcomes);
        if (DEBUG)
        {   cout << "Response " << k << endl;
            cout << "n_outcomes=" << n_outcomes << endl;
        }
    }

    // Initialise Gammas

    if (Print==1)
    {   cout << "Initial value for Gamma" << endl;}
    //cout << "fixInit=" << fixInit << endl;

    // For the moment there is no Resume
    resumeGamMat.resize(Nb_Resp);

    Gam_step_regr_Vect.resize(Nb_Resp);
    if (fixInit==false)
    {   initGamMat.resize(Nb_Resp);

        //It is not going to be used, just to be able to pass the argument
        // initGamMat[k] to Kernel_Single_Gamma.initialise(...);
        //Same thing for resumeGamMat
        for (unsigned int k=0;k<Nb_Resp;k++)
        {
            cout << "Stepwise Regression number " << k+1 << " out of " << Nb_Resp << endl;
            //cout << "****************************" << endl;

            Gammas[k].Stepwise_Regression(Gam_step_regr_Vect[k],
                                          nConfounders,
                                          Data_X,
                                          Vector_Data_Y[k],
                                          //vect_RMSE,
                                          vect_vect_RMSE[k],
                                          Pvalue_enter,
                                          Pvalue_remove,
                                          HESS_Information);
            //cout << endl;
            //unsigned int Size_Gam_step_regr=Gam_step_regr_Vect[k].size();
            //cout <<  "Gam_step_regr_Vect[" << k << "].size()=" << Size_Gam_step_regr << endl;
            //cout << endl;
            //Display_Matrices(Gam_step_regr_Vect[k]);
            //cout << endl;
        }
        cout << endl;

        if (Print==1)
        {
            cout << "End of Stepwise Regression" << endl;
            cout << endl;
        }

    }
    else
    {
        // Read the initial gamma from file
        fstream inputFile;
        inputFile.open(filename_init.c_str(),ios::in);
        if(!inputFile.is_open()){
          cout << "Input file not found" << endl;
          exit(-1);
        }

        Read_File_Int(inputFile, initGamMat);
        inputFile.close();

        if (Print==1)
        {
            cout << endl;
            cout << "Initial Gamma read from the file " << filename_init << " :" << endl;
            cout << endl;

            if (DEBUG)
            {
              Write_Matrix(cout,
                           initGamMat,
                           "\t");
              cout << endl;
            }
        }


    }

    //if (Print==1){cout << "Initialise Gammas" << endl;}

    for (unsigned int k=0;k<Nb_Resp;k++)
    {   Gammas[k].Initialise(postProcessOnly,
                             nb_chains,
                             pX,
                             resumeRun,
                             resumeGamMat[k],
                             nConfounders,
                             fixInit,
                             Gam_step_regr_Vect[k],
                             iso_T_Flag,
                             maxPGamma,
                             RandomNumberGenerator,
                             initGamMat[k],
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
                             HESS_Information);

        if (DEBUG)
        {   cout << "Content of Gammas for response nb " << k << " :" << endl;
            cout << endl;

            Write_Matrix(cout,
                         Gammas[k].vect_gam,
                         "\t");
            cout << endl;
        }

    }


    //if (Print==1){cout << "End of initialisation for the Gammas" << endl << endl;}


    // Initial value for g, already done;
    //g=g_init;
    if (Print==1)
    {
        cout << "Initial value for g: " << g[0] << endl;
        cout << endl;
    }

    //Initial value for Omega_k
    if (Print==1) {cout << "Initial value for Omega_k:" << endl;}
    //if (OmegaKSampleFlag==true)
    if (fixInit_omega_k==false)
    {
        for (unsigned int C=0;C<nb_chains;C++)
        {
            (Omega_PerLine[C]).resize(Nb_Resp);
            for (unsigned int k=0;k<Nb_Resp;k++)
            {
                double new_Omega=myrand(RandomNumberGenerator);
                Omega_PerLine[C][k]=new_Omega;
            }
        }
        if (Print==1)
        {
            cout << "   Random initial values for Omega_k" << endl << endl;
        }
    }
    else
    {
        vector < double > init_omega_k;

        fstream inputFile;
        inputFile.open(filename_init_omega_k.c_str(),ios::in);
        if(!inputFile.is_open()){
          cout << " Input file for Omega_k not found" << endl;
          exit(-1);
        }

        Read_File_Double(inputFile, init_omega_k);
        inputFile.close();

        for (unsigned int C=0;C<nb_chains;C++)
        {
            Omega_PerLine[C]=init_omega_k;
        }
        if (Print==1)
        {
            cout << "   Omega_k initialised from the file " << filename_init_omega_k << endl << endl;
        }

    }
    if (DEBUG)
    {   cout << endl;
        cout << "Content of Omega_PerLine :" << endl;
        cout << endl;
        Write_Matrix(cout,
                     Omega_PerLine,
                     "\t");
        cout << endl;
    }
    if (Print==1) {cout << "Initial value for Rho_PerCol: " << endl;}

    //if (RhoJSampleFlag==true)
    if (fixInit_rho_j==false)
    {   for (unsigned int C=0;C<nb_chains;C++)
        {   Rho_PerCol[C].resize(pX);
            for (unsigned int j=0;j<pX;j++)
            {   double new_Rho=myrand(RandomNumberGenerator);
                Rho_PerCol[C][j]=new_Rho;
            }
        }
        if (Print==1)
        {   cout << "   Random initial values for Rho_j" << endl << endl;}
    }
    else
    {   vector < double > init_rho_j;

        fstream inputFile;
        inputFile.open(filename_init_rho_j.c_str(),ios::in);
        if(!inputFile.is_open()){
          cout << " Input file for Rho_j not found" << endl;
          exit(-1);
        }

        Read_File_Double(inputFile, init_rho_j);
        inputFile.close();

        for (unsigned int C=0;C<nb_chains;C++)
        {
            Rho_PerCol[C]=init_rho_j;
        }
        if (Print==1)
        {   cout << "   Rho_j initialised from the file " << filename_init_rho_j << endl << endl;}

    }

    if (DEBUG)
    {   cout << endl;
        cout << "Content of Rho_PerCol :" << endl;
        cout << endl;
        Write_Matrix(cout,
                     Rho_PerCol,
                     "\t");
        cout << endl;
    }


    Save_Rho_History(History_Rho,Rho_PerCol,0);
    Save_Omega_History(History_Omega,Omega_PerLine,0);
    Save_Average_Omega_k_j(Average_Omega_k_j,Omega_PerLine,Rho_PerCol,
                           0,burn_in);

    //History_g[0][0]=g;
    History_g[0]=g;

    /*
      ##############################################################################

                            PRIOR PARAMETERS SETUP

                         WHEN SEVERAL RESPONSES (later versions), NEEDS TO BE DONE
                                AFTER STEPWISE REGRESSION

      ##############################################################################
    */

    // We need one PR object per response, so that the prior parameters can vary and adapt
    // to each local data
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        PR_per_Resp[k].set_PR_param(E_p_gam_from_read,
                                    Sd_p_gam_from_read,
                                    pX,
                                    Vector_Data_Y[k]->size2, // Nb of outcomes per response
                                    nX,
                                    lambda,
                                    vect_vect_RMSE[k],
                                    P_mutation_from_read,
                                    P_sel_from_read,
                                    P_csvr_r_from_read,
                                    P_DR_from_read);
        /*==============================================================
                        DIFFERENCE WITH ESS
        ============================================================*/

       // : Formula used in the Matlab version of HESS
       // It concerns the prior for g, it is the same for all k

       if (Pr_g_a_flag==true)
       {
        PR_per_Resp[k].alpha=Pr_g_a;
       }
       else
       {
         if (Single_g==true)
         {
           PR_per_Resp[k].alpha=0.5*(double)Nb_Resp+((double)Nb_Resp-1.0);
         }
         else
         {
           PR_per_Resp[k].alpha=0.5;
         }
       }

       if (Pr_g_b_flag==true)
       {
        PR_per_Resp[k].beta=Pr_g_b;
       }
       else
       {
         if (Single_g==true)
         {
           PR_per_Resp[k].beta=0.5*(double)nY * (double)Nb_Resp;
         }
         else
         {
           PR_per_Resp[k].beta=0.5*(double)nY;
         }
       }


       /*==============================================================
       ============================================================*/
       if(!postProcessOnly){
           //if (Print==1)
           // Debug
           if (1==2)
           {
               cout << "Prior Parameters for Response " << k+1 << ":" << endl;
               cout << endl;
               PR_per_Resp[k].display_prior_param();
               cout << endl;
           }
       }

    }


   //Priors for OMEGAS
       // We just reuse the code, we are going to use only 2 parameters,
       // namely a_pi and b_pi, which depend on
       // E_p_gam_from_read and Sd_p_gam_from_read
   Omega_k_PR.set_PR_param(E_p_gam_from_read,
                           Sd_p_gam_from_read,
                           pX,
                           1, //This parameter is not used here
                           nX,
                           lambda,
                           vect_vect_RMSE[0], // not used, just to pass a parameter
                           P_mutation_from_read,
                           P_sel_from_read,
                           P_csvr_r_from_read,
                           P_DR_from_read);

   for (unsigned int k=0;k<Nb_Resp;k++)
   {
           gsl_vector_free(vect_vect_RMSE[k]);
   }

   if (!a_om_k_Flag){
       a_om_k=Omega_k_PR.a_pi;
   }
   if (!b_om_k_Flag){
       b_om_k=Omega_k_PR.b_pi;
   }
   if (!c_rh_j_Flag){
       c_rh_j=1.2;
   }
   if (!d_rh_j_Flag){
       d_rh_j=1.2;
   }

   cout << endl << "****************************************************************" << endl
        << "******************** Omega parameters ********************" << endl
        << "\ta_om_k = " << a_om_k << endl
        << "\tb_om_k = " << b_om_k << endl
        << "\tc_rh_j = " << c_rh_j << endl
        << "\td_rh_j = " << d_rh_j << endl;
   cout << endl;
   cout << endl << "****************************************************************" << endl
        << "******************** g parameters ********************" << endl
        << "\tPr.g.a = " << PR_per_Resp[0].alpha << endl
        << "\tPr.g.b = " << PR_per_Resp[0].beta << endl;
   cout << endl;
   cout << "\t################" << endl
        << "****************************************************************" << endl
        << "****************************************************************" << endl << endl;

    /*
         ##############################################################################

                                Preliminary sweep

         ##############################################################################
    */

   //DEBUG the log-posterior computations !!
  /*
           unsigned int k=0;
           double tmplogMargLik,tmplogPost=0.0;

           cout << endl;
           cout << "Preliminary Computation with Response " << k+1 << endl;
           unsigned int chain_pos=Gammas[k].chain_idx[0];
           cout << "chain_pos = " << chain_pos << endl;
           computeLogPosterior(Gammas[k].vect_gam[chain_pos],
                               0,//Chain
                               k,//Response
                               Data_X,
                               tmplogMargLik,
                               tmplogPost,
                               Vector_Data_Y[k],
                               PR_per_Resp[k],
                               gPriorFlag,
                               indepPriorFlag,
                               gSampleFlag,
                               lambda,
                               g[0],
                               cudaFlag,
                               HESS_Information);

           cout << "Log_Marg = " << tmplogMargLik << endl;
           cout << "Log_Post = " << tmplogPost << endl;
           cout << "g = " << g[0] << endl;
           //cout << "Pr.delta = " << PR.delta << endl;
           cout << "Pr.k = " << (PR_per_Resp[k]).k << endl;


           cout << endl;

           return;
    */


    if (Print==1)
    {
        cout << "Preliminary sweep..." << endl;
        cout << endl;
    }

    unsigned int sweep=0;
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        //We need only the k^th response from the data

        gsl_matrix * Data_Y_k=Vector_Data_Y[k];

        double g_value;
        if (Single_g==true)
            g_value=g[0];
        else
            g_value=g[k];
        Gammas[k].First_Sweep(k,
                              nb_chains,
                              Data_X,
                              Data_Y_k,
                              //PR,
                              PR_per_Resp[k],
                              gPriorFlag,
                              indepPriorFlag,
                              gSampleFlag,
                              lambda,
                              //g,
                              g_value,
                              cudaFlag,
                              sweep,
                              HESS_Information);

    }

    endTime = std::clock();
    double Setup_Time = (endTime-startTime)/(double)(CLOCKS_PER_SEC);
    cout << "Setup Time: ";
    Show_Time(Setup_Time);
    cout << endl << endl;

    /*
    // This saves the models for later post-processing
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        // CHANGE: PUT THE FIRST LINE INSIDE THE BODY OF save_model_per_sweep,
        // AND MAKE THE LATTER TAKE sweep AS A PARAMETER
        Gammas[k].Current_Sweep=sweep;
        Gammas[k].save_model_per_sweep();
    }
    */

    /*
         ##############################################################################

                                DEBUG_sample_g

         ##############################################################################
    */
    if (DEBUG_sample_g==true)
    {
      cout << setprecision(10);

      sample_g_HESS(nb_chains,
              My_g_AdMH,
              Gammas,
              Data_X,
              Vector_Data_Y,
              gPriorFlag,
              indepPriorFlag,
              gSampleFlag,
              lambda,
              g,
              PR_per_Resp,
              sweep,
              cudaFlag,
              RandomNumberGenerator,
              HESS_Information,
              Single_g);

      return;
    }

    /*
         ##############################################################################

                                MCMC Iterations

         ##############################################################################
    */

    Loop_Start_Time = std::clock();


    cout << "Begin MCMC loop..." << endl << endl;
    for (unsigned int sweep=1;sweep<n_sweeps+1;sweep++)
    {

        //update   Active_Yk
        if(gamSampleFlag)
        {
            Update_Active_Yk(Active_Yk,
                             Prop_Active,
                             Nb_Resp,
                             RandomNumberGenerator);
        }

        unsigned int Nb_Active_Resp=Active_Yk.size();

        // PRELIMINARY SWEEP: Computes the likelihoods

        for (unsigned int k=0;k<Nb_Resp;k++)
        {
            gsl_matrix * Data_Y_k=Vector_Data_Y[k];

            double g_value;
            if (Single_g==true)
                g_value=g[0];
            else
                g_value=g[k];

            Gammas[k].First_Sweep(k,
                                  nb_chains,
                                  Data_X,
                                  Data_Y_k,
                                  //PR,
                                  PR_per_Resp[k],
                                  gPriorFlag,
                                  indepPriorFlag,
                                  gSampleFlag,
                                  lambda,
                                  //g,
                                  g_value,
                                  cudaFlag,
                                  sweep-1,
                                  HESS_Information);

            //Copy the previous value of the marginal and posterior, as well as the previous models
            Gammas[k].n_Models_visited[sweep]=Gammas[k].n_Models_visited[sweep-1];

            for(unsigned int chain=0;chain<nb_chains;chain++)
            {
              Gammas[k].mat_log_marg.matrix[chain][sweep]=Gammas[k].mat_log_marg.matrix[chain][sweep-1];
              Gammas[k].mat_log_cond_post.matrix[chain][sweep]=Gammas[k].mat_log_cond_post.matrix[chain][sweep-1];
            }
        }

        //update Gammas
        if(gamSampleFlag)
        {

            for (unsigned int count=0;count<Nb_Active_Resp;count++)
            {

                unsigned int k=Active_Yk[count];
                if (DEBUG)
                {cout << "Update Response " << k << endl;}

                gsl_matrix * Data_Y_k=Vector_Data_Y[k];

                double g_value;
                if (Single_g==true)
                    g_value=g[0];
                else
                    g_value=g[k];

                Gammas[k].MCMC_Update( k,     //Response number k
                                       sweep,
                                       nb_chains,
                                       Gibbs_n_batch,
                                       Log_Flag,
                                       Data_X,
                                       Data_Y_k,
                                       gPriorFlag,
                                       indepPriorFlag,
                                       gSampleFlag,
                                       lambda,
                                       //g,
                                       g_value,
                                       PR_per_Resp[k],
                                       cudaFlag,
                                       nConfounders,
                                       maxPGamma,
                                       RandomNumberGenerator,
                                       burn_in,
                                       iso_T_Flag,
                                       HESS_Information);
            }
        }


        // update g...

        if(gSampleFlag)
        {

            //cout << endl;

            sample_g_HESS(nb_chains,
                      My_g_AdMH,
                      Gammas,
                      Data_X,
                      Vector_Data_Y,
                      gPriorFlag,
                      indepPriorFlag,
                      gSampleFlag,
                      lambda,
                      g,
                      PR_per_Resp,
                      sweep,
                      cudaFlag,
                      RandomNumberGenerator,
                      HESS_Information,
                      Single_g);
        }

        //update  Omega_PerLine
        if(OmegaKSampleFlag)
        {
            //cout << endl;
            //cout << "Going to update Omega_k" << endl;

            Sample_Omega_k( sweep,
                            RandomNumberGenerator,
                            Omega_k_AdMH,
                            a_om_k,
                            b_om_k,
                            nb_chains,
                            pX,
                            Nb_Resp,
                            Omega_PerLine,
                            Rho_PerCol,
                            Gammas,
                            HESS_Information);

        }

        //update   Rho_PerCol
        if(RhoJSampleFlag)
        {

            Sample_Rho_j(sweep,
                         RandomNumberGenerator,
                         Rho_j_AdMH,
                         c_rh_j,
                         d_rh_j,
                         nb_chains,
                         pX,
                         Nb_Resp,
                         Omega_PerLine,
                         Rho_PerCol,
                         Gammas,
                         HESS_Information);
        }


        // Save Results of the sweep

        if (sweep>burn_in)
        {
            // This saves the models for later post-processing
            // CHANGE: PUT THE FIRST LINE INSIDE THE BODY OF save_model_per_sweep,
            // AND MAKE THE LATTER TAKE sweep AS A PARAMETER
            for (unsigned int k=0;k<Nb_Resp;k++)
            {
                Gammas[k].Current_Sweep=sweep;
                Gammas[k].save_model_per_sweep();
            }
        }


        // Rho
        Save_Rho_History(History_Rho,Rho_PerCol,sweep);

        // Omega_k
        Save_Omega_History(History_Omega,Omega_PerLine,sweep);

        //Average Omega
        if (sweep>burn_in)
        {
            Save_Average_Omega_k_j(Average_Omega_k_j,Omega_PerLine,Rho_PerCol,
                                   sweep,burn_in);

            // Save Average Gamma_kj
            for (unsigned int k=0;k<Nb_Resp;k++)
            {
              //Index corresponding to the current sweep:
              unsigned int last_pos_in_List=Gammas[k].List_Models.size()-1;
              for(unsigned int Ind_Var=1;Ind_Var<Gammas[k].List_Models[last_pos_in_List].size();Ind_Var++)
              {
                unsigned int j=Gammas[k].List_Models[last_pos_in_List][Ind_Var];
                Count_Visits_Gamma_kj[k][j]=Count_Visits_Gamma_kj[k][j]+1;
              }
            }
        }

        //g
        //History_g[sweep][0]=g;
        History_g[sweep]=g;

        std::clock_t PreviousEnd=endTime;
        endTime = std::clock();
        double Elapsed = (endTime-Loop_Start_Time)/(double)(CLOCKS_PER_SEC);
        double remainingTime=((double) n_sweeps-(double)sweep)/((double) sweep)*Elapsed;
        double durationSweep=(endTime-PreviousEnd)/(double)(CLOCKS_PER_SEC);

        cout << "End of sweep " << sweep << "/" << n_sweeps << endl;
        cout << "Lasted ";
        Show_Time(durationSweep);
        cout << endl;
        //cout << "Elapsed Time: ";
        //Show_Time(Elapsed);
        //cout << endl;
        cout << "Estimated Remaining Time: ";
        Show_Time(remainingTime);
        cout << endl << endl;


        // WRITE HISTORY: only chain 0
        if (WRITE_BY_1_LINE==true)
        {
            // All the Rho's
            Write_1Line(Output_File_History_Rho,
                         Rho_PerCol[0],
                         Field_Separator);

            // All the Omega_k's
            Write_1Line(Output_File_History_Omega,
                         Omega_PerLine[0],
                         Field_Separator);


            // g
            Write_1Line(Output_File_History_g,
                         g,
                         Field_Separator);

        }



    }
    cout << endl;

    endTime = std::clock();
    double MCMCLoopTime = (endTime-Loop_Start_Time)/(double)(CLOCKS_PER_SEC);



    /*
         ##############################################################################

                                    Output MCMC

         ##############################################################################
    */


    // Writing History
    if (WRITE_BY_1_LINE==false)
    {
        // All the Rho's
        Write_Matrix(Output_File_History_Rho,
                     History_Rho,
                     Field_Separator);

        // All the Omega_k's
        Write_Matrix(Output_File_History_Omega,
                     History_Omega,
                     Field_Separator);

        // All the Omega_k's
        Write_Matrix(Output_File_History_g,
                     History_g,
                     Field_Separator);


    }

    // The average of the Omega_k_j
    Write_Matrix(Output_File_Av_Omega,
                 Average_Omega_k_j,
                 Field_Separator);

    Write_Matrix(Output_File_Counts_Gamma_kj,
                 Count_Visits_Gamma_kj,
                 Field_Separator);


    Output_File_History_Rho.close();
    Output_File_History_Omega.close();
    Output_File_Av_Omega.close();
    Output_File_History_g.close();
    Output_File_Counts_Gamma_kj.close();


    cout << endl << "MCMC loop time: ";
    Show_Time(MCMCLoopTime);
    cout << endl;



/*
     ##############################################################################

                                Post-processing

     ##############################################################################
*/


    PostProc_Start_Time = std::clock();



    // 1: Averages of Omega_k and Rho_j:

    vector<double> Average_Omega_k;
    Average_Omega_k.resize(Nb_Resp);
    Compute_Average_Histo(Average_Omega_k,
                               History_Omega,
                               burn_in);

    vector<double> Average_Rho_j;
    Average_Rho_j.resize(pX);
    Compute_Average_Histo(Average_Rho_j,
                               History_Rho,
                               burn_in);

    // 2: Tail probbilities for Rho_j:

    // for each j, count how many Rho_j after burnin are larger than 1
    vector<double> TailProb_Rho_j;
    TailProb_Rho_j.resize(pX);
    double Thresh_TP_Rho_j=1.0;
    for (unsigned int j=0;j<pX;j++)
    {
        double Number=0.0;
        for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
        {
            if (History_Rho[sweep][j]>Thresh_TP_Rho_j)
            {
                Number=Number+1.0;
            }
        }
        TailProb_Rho_j[j]=Number/((double)(n_sweeps+1-burn_in));

    }


    /*
         ##############################################################################

                                    Output basic post-processing

         ##############################################################################
    */


    // Average of Omega_k
    Write_Vector(Output_File_Av_Omega_k,
                 Average_Omega_k);

    // Average of Rho_j
    Write_Vector(Output_File_Av_Rho_j,
                 Average_Rho_j);
    // Tail Probabilities for Rho_j
    Write_Vector(Output_File_TP_Rho_j,
                 TailProb_Rho_j);


    Output_File_Av_Rho_j.close();
    Output_File_Av_Omega_k.close();
    Output_File_TP_Rho_j.close();




/*
     ##############################################################################

                                Compute MPPI

     ##############################################################################
*/



    if (HESS_do_Postproc)
    {

      cout << "Starting Post-Processing..." << endl;
      cout << endl;



        // Steps:
        // 1) Find unique models
        // 2) For each model, loop through the iterations and compute the
        //    marginal likelihood.
        // 2')Compute normalising constant for each iteration
        // 3) Find best marginal model
        // 4) Compute frequencies...
        // 5) Probability of inclusion: Sum over all the visited models



        //***************************************************************************
        //
        //                           Find unique models
        //
        //***************************************************************************
        cout << "       Listing unique models..." << endl;

        for (unsigned int k=0;k<Nb_Resp;k++)
        {

          int pos_null_model;
          pos_null_model=getUniqueList(Gammas[k].Unique_List_Models,
                                       Gammas[k].List_Models,
                                       Gammas[k].n_Models_visited,
                                       burn_in,
                                       pX,
                                       nConfounders,
                                       Include_Single_Models);
        }

        // We can stop write here the list of unique models and stop the main program.
        // That file could then be read by the post-processing program.

        //***************************************************************************
        //
        //                      Compute marginal posterior probabilities,
        //                         including the normalising constant
        //
        //***************************************************************************


        // Try and parallelise over all responses the expensive computations in one go.


        vector < vector <double> > Matrix_Mod_Marg_Post_Prob;
        Matrix_Mod_Marg_Post_Prob.resize(Nb_Resp);

        vector < vector < unsigned int> > Best_Models_Marginal;
        Best_Models_Marginal.resize(Nb_Resp);
        vector < double> Prob_Best_Models_Marginal;
        Prob_Best_Models_Marginal.resize(Nb_Resp);

        cout << endl <<"Post-processing computations of the log-posterior:" << endl;


            //***************************************************************************
            //
            //                           Decide here subsampling strategy for
            //                                  MPPI computations
            //
            //***************************************************************************

            // This will impact all the code below



            //***************************************************************************
            //
            //                           Precompute logs (1)
            //
            //***************************************************************************

        vector <double> Total_for_logs;
        Total_for_logs.resize(n_sweeps-burn_in+1);

        if (!HESS_slow_Postproc) // Approximate computation of the logs
        {
            // With the approximation, the sum can be computed once and for all
            unsigned Record_index=0;  // not great implementation

            for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
            {
                Total_for_logs[Record_index]=0.0;
                for (unsigned int j=0;j<pX;j++)
                {

                    Total_for_logs[Record_index]=Total_for_logs[Record_index]+
                                                    History_Rho[sweep][j];

                }
                Record_index=Record_index+1;
            }
        }


        vector < vector <double> > log_rho_j_M;

        cout << "Memory_Limited = " << Memory_Limited << endl;

        if (!Memory_Limited)
        {
            log_rho_j_M.resize(n_sweeps-burn_in+1);
            unsigned Record_index=0;  // not great implementation

            for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
            {

                log_rho_j_M[Record_index].resize(pX);

                for (unsigned int j=0;j<pX;j++)
                {
                    //cout << "j = " << j << endl;
                    //cout << "Dimensions of log_rho_j_M: " << log_rho_j_M.size() << " , " <<
                    //        log_rho_j_M[0].size() << endl;
                    //cout << "Dimensions of History_Rho: " << History_Rho.size() << " , " <<
                    //        History_Rho[0].size() << endl;

                    log_rho_j_M[Record_index][j]=log(History_Rho[sweep][j]);

                }
                Record_index=Record_index+1;
            }
        }



        for (unsigned int k=0;k<Nb_Resp;k++)
        {
            cout << "Response " << k+1 << " / " << Nb_Resp << " ..." << endl;

            //***************************************************************************
            //
            //                           Precompute logs  (2)
            //
            //***************************************************************************
            vector <double> log_omega_k_V;
            log_omega_k_V.resize(n_sweeps-burn_in+1);

            vector < vector <double> > log_1_omega_k_rho_j_M;
            log_1_omega_k_rho_j_M.resize(n_sweeps-burn_in+1);


            if (!Memory_Limited)
            {

                unsigned int Record_index=0;  // not great implementation
                for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
                {
                    log_omega_k_V[Record_index]=log(History_Omega[sweep][k]);
                    Record_index=Record_index+1;
                }

                if (HESS_slow_Postproc)
                {
                    Record_index=0;  // not great implementation
                    for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
                    {
                        log_1_omega_k_rho_j_M[Record_index].resize(pX);
                        Total_for_logs[Record_index]=0.0;
                        for (unsigned int j=0;j<pX;j++)
                        {
                            log_1_omega_k_rho_j_M[Record_index][j]=log(1-History_Omega[sweep][k] * History_Rho[sweep][j]);
                            Total_for_logs[Record_index]=Total_for_logs[Record_index]+
                                                            log_1_omega_k_rho_j_M[Record_index][j];
                        }
                        Record_index=Record_index+1;
                    }
                }
            }
            else // when Limited Memory, only compute the sum of logs
            {
                if (HESS_slow_Postproc)
                {
                    unsigned int Record_index=0;  // not great implementation
                    for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
                    {
                        Total_for_logs[Record_index]=0.0;
                        for (unsigned int j=0;j<pX;j++)
                        {
                            //log_1_omega_k_rho_j_M[Record_index][j]=log(1-History_Omega[sweep][k] * History_Rho[sweep][j]);
                            Total_for_logs[Record_index]=Total_for_logs[Record_index]+
                                                            log(1-History_Omega[sweep][k] * History_Rho[sweep][j]);
                        }
                        Record_index=Record_index+1;
                    }
                }

            }




            //***************************************************************************
            //
            //          Introduce normalising constants for each MCMC iteration
            //
            //***************************************************************************

            vector < double > Normalising_Const_per_Iter;
            // Add Shift_Constants to help compute exponentials of small numbers
            vector < double > Shift_per_Iter_Max;

            //Initialise Shift_per_Iter_Max
            unsigned int Nb_Records=n_sweeps-burn_in+1;

            Shift_per_Iter_Max.resize(Nb_Records);

            for(unsigned int Index_Record=0;Index_Record<Nb_Records;Index_Record++)
            {
                Shift_per_Iter_Max[Index_Record]=-1e10;
            }


            //gsl_matrix * Data_Y_k=Vector_Data_Y[k];

            unsigned int nUniqueModels=Gammas[k].Unique_List_Models.size();

            // For parallelisation:
            //vector < double> log_Post;

            //In the meantime:
            Matrix_Mod_Marg_Post_Prob[k].resize(nUniqueModels);


            // Temporary matrix to store the logs before computing exponentials
            // and normalising constants
            vector < vector < double > > logPost_tmp_Mat;
            logPost_tmp_Mat.resize(nUniqueModels);


            //cout << "Before Loop 1" << endl;
            //cout << "       Compute log posteriors..." << endl;

            //***************************************************************************
            //
            //                             LOOP 1:
            //
            //              For each model, compute its MC LOG-posterior
            //
            //***************************************************************************

                                //**********************************\\
                                //  Preliminary:                    \\
                                //  Compute Y'Y and store it        \\
                                //**********************************\\


            unsigned int pY=Vector_Data_Y[k]->size2;
            gsl_matrix * matYTY=gsl_matrix_calloc(pY,pY);
            gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,Vector_Data_Y[k],Vector_Data_Y[k],0.0,matYTY);

            // Memory LEAK IN THE FOLLOWING LOOP
            for(unsigned int model=0;model<nUniqueModels;model++)
            {
                //***************************************************************************
                //
                //                              Build Gamma
                //
                //***************************************************************************

                unsigned int nVarsIn,offset;

                nVarsIn=Gammas[k].Unique_List_Models[model][3];
                offset=4;

                vector < unsigned int > gamma;
                vector < unsigned int > gamma_list_vars;
                gamma.resize(pX);
                for (unsigned int j=0;j<gamma.size();j++)
                {
                    gamma[j]=0;
                }

                if(nVarsIn>0)
                {
                  for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++)
                  {
                      unsigned int Variable_Explicit=Gammas[k].Unique_List_Models[model][currVar];
                      gamma[Variable_Explicit]=1;
                      gamma_list_vars.push_back(Variable_Explicit);
                  }
                }


                //***************************************************************************
                //
                //          Store the R2 used in the marginal likelihood computations,
                //                                  as well as Y'Y
                //
                //***************************************************************************

                gsl_matrix * R2_Mat_GPriors=gsl_matrix_calloc(pY,pY);


                //cout << "Compute Residual Matrix..." << endl;
                Get_R2_Mat_GPriors(R2_Mat_GPriors,
                                   gamma,
                                   Data_X,
                                   Vector_Data_Y[k],
                                   cudaFlag);



                //cout << "Start Postprocessing MCMC Loop..." << endl;
                unsigned int Index_Record=0;
                // We run only after burnin
                for (unsigned int sweep=burn_in;sweep<n_sweeps+1;sweep++)
                {
                    Index_Record=sweep-burn_in;

                    //***************************************************************************
                    //
                    //                  Compute p(Y|...).p(gamma|Omega)
                    //              For the moment, only coded for the g-priors
                    //
                    //***************************************************************************


                        double log_posterior_term=0.0;

                        //log_posterior_term=3.14*(double)Index_Record;

                        // We precompute the sum of logs once and for all


                        double Total;

                        /*
                        cout << "sweep = " << sweep << endl;
                        cout << "Index_Record = " << Index_Record << endl;
                        cout << "log_rho_j_M.size() = " << log_rho_j_M.size() << endl;
                        cout << "log_rho_j_M[" << Index_Record << "].size() = " << log_rho_j_M[sweep].size() << endl;
                        */

                        double g_value;
                        if (Single_g==true)
                            g_value=History_g[sweep][0];
                        else
                            g_value=History_g[sweep][k];

                        log_posterior_term=PostProc_Log_Marg_g_prior(//gamma, // Old algorithm
                                                                     gamma_list_vars, // new algorithm
                                                                     R2_Mat_GPriors,
                                                                     matYTY,
                                                                     Data_X,
                                                                     Vector_Data_Y[k],
                                                                     PR_per_Resp[k],
                                                                     gPriorFlag,
                                                                     indepPriorFlag,
                                                                     //History_g[sweep][0],
                                                                     g_value,
                                                                     lambda,
                                                                     cudaFlag,
                                                                     History_Omega[sweep][k],
                                                                     History_Rho[sweep],
                                                                     HESS_slow_Postproc,
                                                                     Memory_Limited,
                                                                     log_omega_k_V[Index_Record], // CAREFUL !
                                                                     log_rho_j_M[Index_Record], // CAREFUL !
                                                                     log_1_omega_k_rho_j_M[Index_Record], //CAREFUL !
                                                                     Total_for_logs[Index_Record]); //CAREFUL !


                        // DEBUG
                        //cout << "log_posterior_term = " << log_posterior_term << endl;

                        logPost_tmp_Mat[model].push_back(log_posterior_term);

                        //Update the overall maximum
                        Shift_per_Iter_Max[Index_Record]=max(Shift_per_Iter_Max[Index_Record],
                                                             logPost_tmp_Mat[model][Index_Record]);

                        //cout << "Shift_per_Iter_Max[" << 0 << "] = " << Shift_per_Iter_Max[0] << endl;


                }// end for MCMC iterations
                //cout << "End Postprocessing MCMC Loop." << endl;




                gsl_matrix_free(R2_Mat_GPriors);



            }// end LOOP 1: we have computed all the marginal likelihood terms
            gsl_matrix_free(matYTY);

            //***************************************************************************
            //
            //                              LOOP 2:
            //
            //            For each MCMC iteration, compute the normalising constant,
            //               using the shift factor at the corresponding iteration.
            //
            //***************************************************************************

            Normalising_Const_per_Iter.resize(Nb_Records);
            for (unsigned int Record=0;Record<Nb_Records;Record++)
            {
                Normalising_Const_per_Iter[Record]=0.0;
                for(unsigned int model=0;model<nUniqueModels;model++)
                {
                    Normalising_Const_per_Iter[Record]=
                            Normalising_Const_per_Iter[Record]+
                            exp(logPost_tmp_Mat[model][Record]-Shift_per_Iter_Max[Record]);

                    //DEBUG
                    //cout << "Normalising_Const_per_Iter[" << Record << "] = " << Normalising_Const_per_Iter[Record] << endl;
                }

            }


            //***************************************************************************
            //
            //                                  LOOP 3:
            //
            //              For each model, compute its posterior probability,
            //             using the normalising constant AND the shift constant
            //                      also, find the Best marginal model
            //
            //***************************************************************************

            for(unsigned int model=0;model<nUniqueModels;model++)
            {
                double post_prob=0.0;
                for (unsigned int Record=0;Record<Nb_Records;Record++)
                {

                    post_prob=post_prob+
                            (exp(logPost_tmp_Mat[model][Record]-Shift_per_Iter_Max[Record])
                             /
                            Normalising_Const_per_Iter[Record]
                             );
                    //DEBUG
                    //cout << "post_prob = " << post_prob << endl;
                }
                // MCMC average:
                post_prob=post_prob/((double)Nb_Records);

                Matrix_Mod_Marg_Post_Prob[k][model]=post_prob;

            }


            // Find the best marginal model
            //*****************************
            unsigned int Best_Mod_Marg_ind=0;
            double Best_Prob=0.0;
            for (unsigned int model=0;model<nUniqueModels;model++)
            {
                if (Matrix_Mod_Marg_Post_Prob[k][model]>Best_Prob)
                {
                        Best_Prob=Matrix_Mod_Marg_Post_Prob[k][model];
                        Best_Mod_Marg_ind=model;
                }
            }

            //Save the best model:
            //*******************
            unsigned int nVarsIn=Gammas[k].Unique_List_Models[Best_Mod_Marg_ind][3];
            unsigned int offset=4;
            if(nVarsIn>0)
            {
                Best_Models_Marginal[k].resize(nVarsIn);

                for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++)
                {
                  /*
                    Best_Models_Marginal[k][currVar-offset]=
                            Gammas[k].Unique_List_Models[Best_Mod_Marg_ind][currVar];
                  */
                  //CORRECTION:
                  Best_Models_Marginal[k][currVar-offset]=
                          Gammas[k].Unique_List_Models[Best_Mod_Marg_ind][currVar]+1;

                }
            }

            Prob_Best_Models_Marginal[k]=Best_Prob;

        }// end for k


        //***************************************************************************
        //
        //              Compute frequencies. reweighted ?!?
        //              ...
        //
        //***************************************************************************



        //***************************************************************************
        //
        //                  Compute Marginal probabilities of inclusion
        //
        //***************************************************************************

        cout << endl << "Compute Marginal probabilities of inclusion..." << endl;

        vector < vector < double > > Marg_Prob_Inclusion;
        Marg_Prob_Inclusion.resize(Nb_Resp);
        for (unsigned int k=0;k<Nb_Resp;k++)
        {
            Marg_Prob_Inclusion[k].resize(pX);
            for (unsigned int j=0;j<pX;j++)
            {
                Marg_Prob_Inclusion[k][j]=0;
            }

            unsigned int nUniqueModels=Gammas[k].Unique_List_Models.size();

            for (unsigned int model=0;model<nUniqueModels;model++)
            {
                unsigned int nVarsIn=Gammas[k].Unique_List_Models[model][3];
                unsigned int offset=4;
                if(nVarsIn>0)
                {
                    for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++)
                    {
                        unsigned int Name_Var=Gammas[k].Unique_List_Models[model][currVar];
                        Marg_Prob_Inclusion[k][Name_Var]=Marg_Prob_Inclusion[k][Name_Var]+
                                Matrix_Mod_Marg_Post_Prob[k][model];
                    }
                }

            }

        }


        //***************************************************************************
        //
        //                       Other desirable output ?!?
        //              * for each sweep, for each j, number of gammas=1 in column j
        //              * Temperatures
        //              * Nb of variables in the models
        //
        //***************************************************************************





/*
     ##############################################################################

                                Output MPPI

     ##############################################################################
*/


        cout << endl << "Writing post-processing output" << endl << endl;

        // Best model for each response (best in the marginal sense)
        Write_Matrix(Output_File_Best_Models_Marg,
                     Best_Models_Marginal,
                     Field_Separator);

        // Posterior Probability for each the best model for each k
        Write_Vector(Output_File_Prob_Best_Models_Marg,
                     Prob_Best_Models_Marginal);

        // Marginal probability of inclusion for each gamma_kj
        Write_Matrix(Output_File_Marg_Prob_Incl,
                     Marg_Prob_Inclusion,
                     Field_Separator);

        Output_File_Best_Models_Marg.close();
        Output_File_Prob_Best_Models_Marg.close();
        Output_File_Marg_Prob_Incl.close();

/*
     ##############################################################################

                                Output All Models Posterior Probabilities

     ##############################################################################
*/

        // Content of the file:
        // 1st column: index of the response
        // 2nd column: posterior probability of the model in question
        // other columns: list of variables in the model

        for (unsigned int k=0;k<Nb_Resp;k++)
        {
          for(unsigned int model=0;model<Gammas[k].Unique_List_Models.size();model++)
          {
            unsigned int nVarsIn,offset;

            nVarsIn=Gammas[k].Unique_List_Models[model][3];
            offset=4;

            Output_File_Models_Prob << k+1 << Field_Separator
                                       << Matrix_Mod_Marg_Post_Prob[k][model];

            if(nVarsIn>0)
            {
              for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++)
              {
                Output_File_Models_Prob << Field_Separator << Gammas[k].Unique_List_Models[model][currVar]+1;
              }
            }
            Output_File_Models_Prob << endl;
          }
        }

        Output_File_Models_Prob.close();


    }// End if(HESS_do_Postproc)



    endTime = std::clock();
    double PostProcessTime = (endTime-PostProc_Start_Time)/(double)(CLOCKS_PER_SEC);

    cout << endl << endl;

    double Total_Time = (endTime-startTime)/(double)(CLOCKS_PER_SEC);

    cout << "Setup Time: ";
    Show_Time(Setup_Time);
    cout << endl;

    cout << "MCMC loop time: ";
    Show_Time(MCMCLoopTime);
    cout << endl;

    if (HESS_do_Postproc)
    {
        cout << "Post-process time: ";
    } else {
        cout << "Frequencies computation time: ";
    }
    Show_Time(PostProcessTime);
    cout << endl;

    cout << "Total time: ";
    Show_Time(Total_Time);
    cout << endl;



    cout << endl << endl;


    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        Gammas[k].free();
    }


    //delete My_g_AdMH;
    for (unsigned int k=0;k<My_g_AdMH.size();k++)
    {
        delete My_g_AdMH[k];
    }
    //delete Omega_k_AdMH;
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        for (unsigned int c=0;c<nb_chains;c++)
        {
            delete Omega_k_AdMH[k][c];
        }
    }
    //delete Rho_j_AdMH;
    for (unsigned int j=0;j<pX;j++)
    {
        for (unsigned int c=0;c<nb_chains;c++)
        {
            delete Rho_j_AdMH[j][c];
        }
    }



    gsl_matrix_free(Data_X);

    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        gsl_matrix_free(Vector_Data_Y[k]);
    }



    gsl_rng_free(RandomNumberGenerator);



}




void
sample_g_HESS(unsigned int n_chains,
          vector <AdMH *> My_g_AdMH,
          vector <Kernel_Single_Gamma > Gammas,
          gsl_matrix *mat_X,
          vector <gsl_matrix *> &Vector_Data_Y,
          bool gPriorFlag,
          bool indepPriorFlag,
          bool gSampleFlag,
          double lambda,
          vector <double> &g,
          vector <Prior_param> PR_per_Resp,
          unsigned int sweep,
              bool cudaFlag,
              gsl_rng *RandomNumberGenerator,
              Model_Information &Model_Data,
          bool Single_g)
{

    //cout << "Just entered sample_g_HESS" << endl;

#if Fixed_Random



#endif



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
  */

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

      sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);
      g_prop=exp(sample_tmp);

      if (DEBUG_sample_g==true)
      {
        cout << "mean proposal: " << norm_mean << endl;
        cout << "sd proposal: " << norm_sd << endl;

#if Fixed_Random
        cout << "Ext_Random_Nbrs = " << Ext_Random_Nbrs[Line_Ext_Rnd_Nmbr][0] << endl;
        sample_tmp=norm_mean+norm_sd*Ext_Random_Nbrs[Line_Ext_Rnd_Nmbr][0];
        //Line_Ext_Rnd_Nmbr++;

        //g_prop=g[0]+norm_sd;
        g_prop=exp(sample_tmp);
#endif

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

        double current_log_cond_post=(Gammas[k].mat_log_cond_post.matrix)[pos_current_chain][sweep];
        double current_log_marg=(Gammas[k].mat_log_marg.matrix)[pos_current_chain][sweep];

        gsl_matrix *Data_Y_k=Vector_Data_Y[k];

        double propLogMarg,propLogPost;

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
        */


        store_log_marg[k].push_back(propLogMarg);
        store_log_cond_post[k].push_back(propLogPost);

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
              //if(DEBUG_HESS){cout << "alpha_g is NaN" << endl;}
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
    //cout << "With single g, nb_chains = " << n_chains << endl;
    //cout << "cum_diff = " << cum_diff << endl;

    //cout << endl << endl;
  }


  if (Single_g==true)
  {
    // Begin Change from ESS
    logPG=getPriorG(PR_per_Resp[0],gSampleFlag,g[0]);
    logPG_Prop=getPriorG(PR_per_Resp[0],gSampleFlag,g_prop);

    if (DEBUG_sample_g==true)
    {
      cout << "g = " << g[0]
         << ", g_prop = " << g_prop
         << ", logP(g) = " << logPG
         << ", logP(g_prop) = " << logPG_Prop
         << ", cum_diff = " <<  cum_diff;
    }




    cum_diff=cum_diff+logPG_Prop-logPG;
    // End Change from ESS

    cum_diff=cum_diff+log(g_prop)-log(g[0]);

    alpha_g=min(1.0,exp(cum_diff));

    if (DEBUG_sample_g==true)
    {cout <<", alpha_g = " << alpha_g << endl;}

    if(!isnan(alpha_g))
    {
        double rand_test = myrand(RandomNumberGenerator);
#if Fixed_Random

        rand_test=Ext_Random_Nbrs[Line_Ext_Rnd_Nmbr][1];
        Line_Ext_Rnd_Nmbr++;

#endif

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
        //if(DEBUG_HESS){cout << "alpha_g is NaN" << endl;}
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
      /*
      if(DEBUG_HESS){
        cout << "G-adaptation" << endl;
        cout << "sweep=" << sweep
         << " -- My_g_AdMH.n_batch " << (*My_g_AdMH[k]).n_batch
         << " -- Test " << sweep%(*My_g_AdMH[k]).n_batch << endl;
      }
      */
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
      /*
      if(DEBUG_HESS){
        cout << "Sweep " << sweep << " -- Update the AdMH parameters" << endl;
        cout << "Updated g_ls " << (*My_g_AdMH[k]).ls << endl;
      }
      */

    }
  }


}















void Set_Omega_k_AdMH(vector< vector < AdMH *> > Omega_k_AdMH,
                      unsigned int Nb_Resp,
                      unsigned int nb_chains,
                      unsigned int burn_in,
                      double g_AdMH_ls_from_read,
                      unsigned int g_n_batch_from_read,
                      double g_AdMH_optimal_from_read)
{
    double Omega_k_M_min=-log(Nb_Resp)/(2*log(10));
    double Omega_k_M_max=log(Nb_Resp)/(2*log(10));
    double Omega_k_ls=g_AdMH_ls_from_read;
    double Omega_k_delta_n=max(fabs(Omega_k_M_min-Omega_k_ls),fabs(Omega_k_M_max-Omega_k_ls));
    Omega_k_delta_n=Omega_k_delta_n/((double)(burn_in)/(double)(g_n_batch_from_read));
    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        //Omega_k_AdMH[k].resize(nb_chains);
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {
            //Omega_k_AdMH[k][Chain]=new AdMH();

            // We prefer not to use the method set_AdMH,
            // it doesn't do what we want in general
            //Omega_k_AdMH[k][Chain]->set_AdMH(true,
            //          g_n_batch_from_read,
            //          g_AdMH_optimal_from_read,
            //          g_AdMH_ls_from_read,
            //          pX,
            //          burn_in,
            //          g_M_min_input,
            //          g_M_max_input);
            Omega_k_AdMH[k][Chain]->n_batch= g_n_batch_from_read;
            Omega_k_AdMH[k][Chain]->optimal=g_AdMH_optimal_from_read;
            Omega_k_AdMH[k][Chain]->ls=g_AdMH_ls_from_read;
            Omega_k_AdMH[k][Chain]->M[0]=Omega_k_M_min;
            Omega_k_AdMH[k][Chain]->M[1]=Omega_k_M_max;
            Omega_k_AdMH[k][Chain]->delta_n=Omega_k_delta_n;

        }
    }
}


void Set_Rho_j_AdMH(vector< vector < AdMH *> > Rho_j_AdMH,
                    unsigned int pX,
                    unsigned int nb_chains,
                    unsigned int burn_in,
                    double g_AdMH_ls_from_read,
                    unsigned int g_n_batch_from_read,
                    double g_AdMH_optimal_from_read
                    )
{
    double Rho_j_M_min=-log(pX)/(2*log(10));
    double Rho_j_M_max=log(pX)/(2*log(10));
    double Rho_j_ls=g_AdMH_ls_from_read;
    double Rho_j_delta_n=max(fabs(Rho_j_M_min-Rho_j_ls),fabs(Rho_j_M_max-Rho_j_ls));
    Rho_j_delta_n=Rho_j_delta_n/((double)(burn_in)/(double)(g_n_batch_from_read));
    for (unsigned int j=0;j<pX;j++)
    {
        //Rho_j_AdMH[j].resize(nb_chains);
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {
            //Rho_j_AdMH[j][Chain]=new AdMH();

            // We prefer not to use the method set_AdMH,
            // it doesn't do what we want in general
            //Rho_j_AdMH[j][Chain]->set_AdMH(true,
            //          g_n_batch_from_read,
            //          g_AdMH_optimal_from_read,
            //          g_AdMH_ls_from_read,
            //          pX,
            //          burn_in,
            //          g_M_min_input,
            //          g_M_max_input);
            Rho_j_AdMH[j][Chain]->n_batch= g_n_batch_from_read;
            Rho_j_AdMH[j][Chain]->optimal=g_AdMH_optimal_from_read;
            Rho_j_AdMH[j][Chain]->ls=g_AdMH_ls_from_read;
            Rho_j_AdMH[j][Chain]->M[0]=Rho_j_M_min;
            Rho_j_AdMH[j][Chain]->M[1]=Rho_j_M_max;
            Rho_j_AdMH[j][Chain]->delta_n=Rho_j_delta_n;
        }
    }
}


void Sample_Omega_k(unsigned int sweep,
                    gsl_rng *RandomNumberGenerator,
                    //AdMH * Omega_k_AdMH,
                    vector< vector <AdMH *> > Omega_k_AdMH,
                    double a_om_k,
                    double b_om_k,
                    unsigned int nb_chains,
                    unsigned int pX,
                    unsigned int Nb_Resp,
                    vector < vector < double > > &Omega_PerLine,
                    vector < vector < double > > &Rho_PerCol,
                    vector<Kernel_Single_Gamma > &Gammas,
                    Model_Information &Model_Data)
{
    double Curr_Omega_k;
    double Prop_Omega_k;
    double norm_mean;
    double norm_sd;

    double sample_tmp;

    double Alpha_Omega_k;

    vector<double> * Rho_V;

    for (unsigned int Response=0;Response<Nb_Resp;Response++)
    {
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {
        // 1) For each chain and reponse, propose an omega. If violates conditions, reject it

            if (DEBUG)
            {cout << "Chain " << Chain << " and Response " << Response << endl;}

            Curr_Omega_k=Omega_PerLine[Chain][Response];
            Rho_V=&(Rho_PerCol[Chain]);

            norm_mean=Logistic_Trans(Curr_Omega_k);

            norm_sd=exp(Omega_k_AdMH[Response][Chain]->ls);
            sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);

            Prop_Omega_k=Logistic_Trans_Inv(sample_tmp);

            double Max_Omega=0;
            for (unsigned int j=0;j<pX;j++)
            {
                Max_Omega=max(Max_Omega,Prop_Omega_k * (*Rho_V)[j]);

            }

            bool Failed=false;
            if (Max_Omega>1)
            {
                //Reject
                Failed=true;
            }

            if (DEBUG)
            {
                cout << "Failed=" << Failed << endl;
                cout << "Max_Omega=" << Max_Omega << endl;
                cout << "Prop_Omega_k=" << Prop_Omega_k << endl;
            }

            if (Failed==false)
            {

                // 2) Compute the acceptance probability

                double Log_P_Prop_Ok=Log_Pr_Omega_k(Prop_Omega_k,a_om_k,b_om_k,Gammas,Chain,Response,Model_Data);
                double Log_P_Curr_Ok=Log_Pr_Omega_k(Curr_Omega_k,a_om_k,b_om_k,Gammas,Chain,Response,Model_Data);

                Alpha_Omega_k=exp(Log_P_Prop_Ok-Log_P_Curr_Ok);
                Alpha_Omega_k=min((double) 1,Alpha_Omega_k);

                // 3) Accept / Reject
                double U = 1;
                U = myrand(RandomNumberGenerator);

                if (DEBUG)
                {cout << "U=" << U << endl;}

                if (U>Alpha_Omega_k)
                {
                    Failed=true;
                }
            }

        // 4) Update

            if (DEBUG)
            {cout << "Alpha_Omega_k=" << Alpha_Omega_k << "   Prop_Omega_k="
                 << Prop_Omega_k << endl;}

            if (Failed==false)
            {
                Omega_PerLine[Chain][Response]=Prop_Omega_k;

                if (DEBUG)
                {cout << "Move accepted, New Omega=" << Omega_PerLine[Chain][Response] << endl;}

            }
            else
            {
                if (DEBUG)
                {cout << "Move refused, New Omega=" << Omega_PerLine[Chain][Response] << endl;}
            }
            if (DEBUG)
            {cout << endl;}

        // 5) Adapt

            // ************************************************************************************************************

                   // CHANGES: we copy the code from Matlab to perform the adaptation
                   // Corresponding Matlab code:
                   // omega_k_AdMH.ls(k, c) = min(omega_k_AdMH.ls(k, c), (sweep / omega_k_AdMH.nbatch) ^(-1 /2));
           double batchNumber = (double)sweep/(double)(Omega_k_AdMH[Response][Chain]->n_batch);
           Omega_k_AdMH[Response][Chain]->ls=min(Omega_k_AdMH[Response][Chain]->ls,pow(batchNumber,-0.5));

           // ************************************************************************************************************


            if (sweep%Omega_k_AdMH[Response][Chain]->n_batch==0)
            {
                Omega_k_AdMH[Response][Chain]->Adapt_AdMH(sweep);
                //(*Omega_k_AdMH).Adapt_AdMH(sweep);
            }

        }
        if (DEBUG)
        {cout << endl;}

    }


}

void Sample_Rho_j(unsigned int sweep,
                  gsl_rng *RandomNumberGenerator,
                  //AdMH * Rho_j_AdMH,
                  vector< vector <AdMH *> > Rho_j_AdMH,
                  double c_rh_j,
                  double d_rh_j,
                  unsigned int nb_chains,
                  unsigned int pX,
                  unsigned int Nb_Resp,
                  vector < vector < double > > &Omega_PerLine,
                  vector < vector < double > > &Rho_PerCol,
                  vector<Kernel_Single_Gamma > &Gammas,
                  Model_Information &Model_Data)
{

    double Curr_Rho_j;
    double Prop_Rho_j;
    double norm_mean;
    double norm_sd;

    double sample_tmp;

    double Alpha_Rho_j;

    vector<double> * Omega_V;

    for (unsigned int Variable=0;Variable<pX;Variable++)
    {
        for (unsigned int Chain=0;Chain<nb_chains;Chain++)
        {

        // 1) For each chain and Covariate, propose a Rho. If violates conditions, reject it

            if (DEBUG)
            {cout << "Chain " << Chain << " and Variable " << Variable << endl;}

            Curr_Rho_j=Rho_PerCol[Chain][Variable];
            Omega_V=&(Omega_PerLine[Chain]);

            norm_mean=log(Curr_Rho_j);
            norm_sd=exp(Rho_j_AdMH[Variable][Chain]->ls);

            sample_tmp=gennor(norm_mean,norm_sd,RandomNumberGenerator);

            Prop_Rho_j=exp(sample_tmp);

            double Max_Omega=0;
            for (unsigned int k=0;k<Nb_Resp;k++)
            {
                Max_Omega=max(Max_Omega,Prop_Rho_j * (*Omega_V)[k]);
            }

            bool Failed=false;
            if (Max_Omega>1)
            {
                //Reject
                Failed=true;
            }

            if (DEBUG)
            {
                cout << "Failed=" << Failed << endl;
                cout << "Max_Omega=" << Max_Omega << endl;
                cout << "Prop_Rho_j=" << Prop_Rho_j << endl;
            }

            if (Failed==false)
            {
                // 2) Compute the acceptance probability

                double Log_P_Curr_Rj=Log_Pr_Rho_j(Curr_Rho_j,c_rh_j,d_rh_j,Omega_V,Gammas,Chain,Variable,Nb_Resp,Model_Data);
                double Log_P_Prop_Rj=Log_Pr_Rho_j(Prop_Rho_j,c_rh_j,d_rh_j,Omega_V,Gammas,Chain,Variable,Nb_Resp,Model_Data);


                Alpha_Rho_j=exp(Log_P_Prop_Rj-Log_P_Curr_Rj);
                Alpha_Rho_j=min(1.0,Alpha_Rho_j);

                // 3) Accept / Reject

                double U = 1;
                U = myrand(RandomNumberGenerator);
                if (DEBUG)
                {cout << "U=" << U << endl;}

                if (U>Alpha_Rho_j)
                {
                    Failed=true;
                }

            }

        // 4) Update

            if (DEBUG)
            {cout << "Alpha_Rho_j=" << Alpha_Rho_j << "   Prop_Rho_j="
                 << Prop_Rho_j << endl;}

            if (Failed==false)
            {
                Rho_PerCol[Chain][Variable]=Prop_Rho_j;

                if (DEBUG)
                {cout << "Move accepted, New Rho=" << Rho_PerCol[Chain][Variable] << endl;}
            }
            else
            {
                if (DEBUG)
                {cout << "Move refused, New Rho=" << Rho_PerCol[Chain][Variable] << endl;}
            }


        // 5) Adapt


            // ************************************************************************************************************

                   // CHANGES: we copy the code from Matlab to perform the adaptation
                   // Corresponding Matlab code:
                   // rho_j_AdMH.ls(c, j) = min(rho_j_AdMH.ls(c, j), (sweep / rho_j_AdMH.nbatch) ^(-1 /2));   % Condition (1)
           double batchNumber = (double)sweep/(double)(Rho_j_AdMH[Variable][Chain]->n_batch);
           Rho_j_AdMH[Variable][Chain]->ls=min(Rho_j_AdMH[Variable][Chain]->ls,pow(batchNumber,-0.5));

           // ************************************************************************************************************

            if (sweep%Rho_j_AdMH[Variable][Chain]->n_batch==0)
            {
                Rho_j_AdMH[Variable][Chain]->Adapt_AdMH(sweep);
            }

            if (DEBUG)
            {cout << endl;}
        }

    }

}

double Logistic_Trans(double Omega)
{
    double Epsilon=pow(10.0,-10.0);

    double Cut_Omega=Omega;
    Cut_Omega=max(Cut_Omega,Epsilon);
    Cut_Omega=min(Cut_Omega,1-Epsilon);


    double norm_mean=log(Cut_Omega/(1-Cut_Omega));

    return norm_mean;
}

double Logistic_Trans_Inv(double sample_tmp)
{
    double Epsilon=pow(10.0,-10.0);

    double Omega=exp(sample_tmp)/(1+exp(sample_tmp));

    double Cut_Omega=Omega;
    Cut_Omega=max(Cut_Omega,Epsilon);
    Cut_Omega=min(Cut_Omega,1-Epsilon);

    return Cut_Omega;

}


double Log_Pr_Omega_k(double x,
                      double a,
                      double b,
                      vector<Kernel_Single_Gamma > &Gammas,
                      unsigned int Chain,
                      unsigned int Response,
                      Model_Information Model_Data)
{

    if (x>1 || x<0)
    {
        cout << "Error in Log_Pr_Omega_k, x outside range [0,1]" << endl;
        exit(-1);
    }



    unsigned int Pos_Chain=Gammas[Response].chain_idx[Chain];
    vector < unsigned int > & gamma=(Gammas[Response]).vect_gam[Pos_Chain];
    double Temp=Gammas[Response].t_tun->t[Chain];
    vector<double> * Rho_V=&(*(Model_Data.Rho_PerCol))[Chain];
    unsigned int pX=Rho_V->size();

    double Log_P_Gam=0;
    for (unsigned int j=0;j<pX;j++)
    {
        if (gamma[j]==1)
        {
            Log_P_Gam=Log_P_Gam + log(x * (*Rho_V)[j]);
        }
        else
        {
            Log_P_Gam=Log_P_Gam + log(1-x * (*Rho_V)[j]);
        }

    }
    Log_P_Gam=Log_P_Gam/Temp;

    double Log_J=log(x*(1-x));

    double Log_P_Omega_k=(a-1)*log(x)+(b-1)*log(1-x)-gsl_sf_lnbeta(a,b);

    double Res=Log_P_Gam+(Log_P_Omega_k+Log_J)/Temp;


    return Res;


}

double Log_Pr_Rho_j(double x,
                    double c,
                    double d,
                    vector<double> * Omega_V,
                    vector<Kernel_Single_Gamma > &Gammas,
                    unsigned int Chain,
                    unsigned int Variable,
                    unsigned int Nb_Resp,
                    Model_Information Model_Data)
{

    if (x<0)
    {
        cout << "Error in Log_Pr_Omega_k, x outside range [0,infinity)" << endl;
        exit(-1);
    }

    unsigned int Pos_Chain;
    double Omega_k;
    double Temp;
    unsigned int gamma_kj;
    double Log_P_Gam_j=0;
    for (unsigned int Response=0;Response<Nb_Resp;Response++)
    {
        Pos_Chain=Gammas[Response].chain_idx[Chain];

        Omega_k=(*Omega_V)[Response];
        gamma_kj=(Gammas[Response]).vect_gam[Pos_Chain][Variable];

        Temp=Gammas[Response].t_tun->t[Chain];

        if (gamma_kj==1)
        {
            Log_P_Gam_j=Log_P_Gam_j+log(Omega_k*x)/Temp;
        }
        else
        {
            Log_P_Gam_j=Log_P_Gam_j+log(1-Omega_k*x)/Temp;
        }
    }


    double Log_J=log(x);

    double Log_P_Rho_j=(c-1)*log(x)-x/d;
    //double Log_P_Rho_j=(c-1)*log(x)-x/d-gsl_sf_lngamma(c)-c*log(d);

    double Res=Log_P_Gam_j+Log_P_Rho_j+Log_J;

    return Res;

}


void Update_Active_Yk(vector < unsigned int > &Active_Yk,
                      double Proportion,
                      unsigned int Nb_Resp,
                      gsl_rng * RandomNumberGenerator)
{
    Active_Yk.resize(0);

    if (DEBUG)
    {cout << "Active responses:" << endl;}

    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        double U = myrand(RandomNumberGenerator);
        if (U<Proportion)
        {
            Active_Yk.push_back(k);

            if (DEBUG)
            {cout << k << endl;}
        }
    }
    if (DEBUG)
    {cout << endl;}

}

void Preliminary_Tests(unsigned int nX,
                       unsigned int nY,
                       unsigned int Nb_Resp,
                       unsigned int n_sweeps,
                       unsigned int burn_in,
                       bool gamSampleFlag,
                       bool fixInit,
                       bool OmegaKSampleFlag,
                       bool fixInit_omega_k,
                       bool RhoJSampleFlag,
                       bool fixInit_rho_j)
{

    if(nX!=nY){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "The number of observations in data X ("<< nX
       << ") differs from the number of observations in data Y (" << nY
       << "), run stopped" << endl;
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

      exit(-1);
    }


//    if(nX<Nb_Resp){
//      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
//       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
//       << "There are too many outcomes ("<< Nb_Resp
//       << ") compared to the number of observations (" << nX
//       << "), run stopped" << endl;
//      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
//       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

//      exit(-1);
//    }

    if(n_sweeps==0 || burn_in==0){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "The Number of sweeps and/or the burn-in has not been specified" << endl
       << "Use -iter and/or -burn-in option(s) in the command line" << endl
       << "Run stopped" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(-1);
    }
    if(n_sweeps <= burn_in){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "The Number of sweeps: " << n_sweeps << " is lower than " << endl
       << "(or equal to) the burn-in: " << burn_in << " -- Run stopped" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
      exit(-1);
    }



    if (gamSampleFlag==false && fixInit==false) // This test is not enough, we have to test whether a file has been provided...
        // But not necessary to worry about gamSampleFlag, it is not documented for the moment
    {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "Incompatible Options:" << endl
         << "-gam_set is on, and there is no input file with the option -init" << endl
         << " to define the initial Gammas" << " -- Run stopped" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        exit(-1);
    }

    if (OmegaKSampleFlag==false && fixInit_omega_k==false)
    {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "Incompatible Options:" << endl
         << "-omega_k_set is on, and there is no input file with the option -init_omega_k" << endl
         << " to define the initial Omega_k" << " -- Run stopped" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        exit(-1);
    }

    if (RhoJSampleFlag==false && fixInit_rho_j==false)
    {
        cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "Incompatible Options:" << endl
         << "-rho_j_set is on, and there is no input file with the option -init_rho_j" << endl
         << " to define the initial Rho_j" << " -- Run stopped" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
         << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
        exit(-1);
    }

    /*

            FOR GPU:

            Transfer here the initialisation and verifications associated with the GPU.
            For the moment, this code is in:
            - Start_Up_ESS.cpp and
            - Command_Line.cpp

    */


}


void Setup_Model_Parameters_HESS(bool &postProcessOnly,
                                 bool &resumeRun,
                                 bool &fixInit,
                                 unsigned int &maxPGamma,
                                 unsigned int &nb_chains,
                                 unsigned int &nConfounders,
                                 unsigned int &pX,
                                 double &b_t_input,
                                 double &a_t_den_inf_5k,
                                 double &a_t_den_5_10k,
                                 double &a_t_den_sup_10k,
                                 unsigned int &temp_n_batch,
                                 vector <double> &M_input,
                                 unsigned int &burn_in,
                                 double &temp_optimal_input,
                                 bool &gPriorFlag,
                                 bool &indepPriorFlag,
                                 bool &gSampleFlag,
                                 double &lambda,
                                 unsigned int &nX,
                                 bool &cudaFlag,
                                 unsigned int &Gibbs_n_batch,
                                 unsigned int &n_sweeps,
                                 unsigned int &n_top_models,
                                 unsigned int &k_max_from_read,
                                 long &MY_SEED,
                                 double &g_init,
                                 string filename_par,
                                 double &Pvalue_enter,
                                 double &Pvalue_remove,
                                 unsigned int &g_n_batch_from_read,
                                 double &g_AdMH_optimal_from_read,
                                 double &g_AdMH_ls_from_read,
                                 double &g_M_min_input,
                                 double &g_M_max_input,
                                 double &E_p_gam_from_read,
                                 double &Sd_p_gam_from_read,
                                 double &P_mutation_from_read,
                                 double &P_sel_from_read,
                                 double &P_csvr_r_from_read,
                                 double &P_DR_from_read
                                 )
{

    nb_chains=3;

    double maxPGammaFactor = 10;

    //Reading the parameter File
    FILE *fparameter=NULL;
    char str[256]; // string used to read parameters

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
    //double
    Pvalue_enter = 1.0 - pow((1.0 - n_Pvalue_enter),(1.0/(double)(pX)));
    //double
    Pvalue_remove = 1.0 - pow((1.0 - n_Pvalue_remove),(1.0/(double)(pX)));


    //Moves Parameters

    //g Adaptative M-H
    if(XML_PAR.ReadTag("G_N_BATCH", 0, 0, str,256)){
      sscanf(str,"%u",&g_n_batch_from_read);
    }
    if(XML_PAR.ReadTag("G_ADMH_OPTIMAL", 0, 0, str,256)){
      sscanf(str,"%lf",&g_AdMH_optimal_from_read);
    }
    if(XML_PAR.ReadTag("G_ADMH_LS", 0, 0, str,256)){
      sscanf(str,"%lf",&g_AdMH_ls_from_read);
    }
    if(XML_PAR.ReadTag("M_MIN", 0, 0, str,256)){
      sscanf(str,"%lf",&g_M_min_input);
    }
    if(XML_PAR.ReadTag("M_MAX", 0, 0, str,256)){
      sscanf(str,"%lf",&g_M_max_input);
    }
    //Crossover Move
    //k_max_from_read=2;
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


    if (Print==1)
    {

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
            //cout << "g " << g_init << endl;
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
             << "n_Pvalue_remove " << n_Pvalue_remove << endl
             << "Pvalue_enter stepwise " << Pvalue_enter << endl
             << "Pvalue_remove stepwise " << Pvalue_remove << endl
             << "**********************************************************" << endl
             << "**********************************************************" << endl << endl;
        }
    }

}


//  Post-processing


void Save_Rho_History(vector<vector<double> > &History_Rho,
                      vector<vector<double> > Rho_PerCol,
                      unsigned int sweep)
{

    // We save only the first chain
    unsigned int Chain=0;

    unsigned int pX=Rho_PerCol[Chain].size();
    for (unsigned int j=0;j<pX;j++)
    {
        History_Rho[sweep][j]=Rho_PerCol[Chain][j];
    }

}

void Save_Omega_History(vector<vector<double> > &History_Rho,
                      vector<vector<double> > Rho_PerCol,
                      unsigned int sweep)
{

    Save_Rho_History(History_Rho,
                     Rho_PerCol,
                     sweep);
}


void Save_Average_Omega_k_j(vector<vector<double> > &Average_Omega_k_j,
                            vector<vector<double> > Omega_PerLine,
                            vector<vector<double> > Rho_PerCol,
                            unsigned int sweep,
                            unsigned int burn_in)
{
    // We save only the first chain
    unsigned int Chain=0;

    unsigned int pX=Rho_PerCol[Chain].size();
    unsigned int Nb_Resp=Omega_PerLine[Chain].size();

    double N=(double) (sweep-burn_in);
    double Coeff_1=N/(N+1);
    double Coeff_2=1.0/(N+1);

    for (unsigned int k=0;k<Nb_Resp;k++)
    {
        for (unsigned int j=0;j<pX;j++)
        {
            double New_Omega_kj=Omega_PerLine[Chain][k]*Rho_PerCol[Chain][j];

            Average_Omega_k_j[k][j]=Coeff_1*Average_Omega_k_j[k][j]
                                        +Coeff_2*New_Omega_kj;
        }
    }

}

void Compute_Average_Histo(vector<double> &Average,
                            vector< vector< double> > Histo,
                            unsigned int burn_in)
{
    // We save only the first chain

    unsigned int Nb_Values=Histo.size();
    unsigned int NCol=Histo[0].size();

    Average.resize(NCol); // should be done beforehand to avoid memory errors

    for (unsigned int j=0;j<NCol;j++)
    {
        Average[j]=0.0;

        for (unsigned int sweep=burn_in;sweep<Nb_Values;sweep++)
        {
            Average[j]=Average[j]+
                    Histo[sweep][j]/((double)(Nb_Values-burn_in));
        }

    }
}

void Save_g_History(vector<vector<double> > &History_Rho,
                      vector<vector<double> > Rho_PerCol,
                      unsigned int sweep)
{

    Save_Rho_History(History_Rho,
                     Rho_PerCol,
                     sweep);
}




double PostProc_Log_Marg_g_prior(//vector <unsigned int> &Gamma, // Old algorithm
                                 vector <unsigned int> &gamma_list_vars, // New algorithm
                                 gsl_matrix * R2_Mat_GPriors,
                                 gsl_matrix *matYTY,
                                 gsl_matrix *mat_X,
                                 gsl_matrix *matY,
                                 Prior_param PR,
                                 bool gPriorFlag,
                                 bool indepPriorFlag,
                                 double g,
                                 double lambda,
                                 bool cudaFlag,
                                 double omega_k,
                                 vector<double> &rho_j_V,
                                 bool HESS_slow_Postproc,
                                 bool Memory_Limited,
                                 double log_omega_k,
                                 vector<double> &log_rho_j_V,
                                 vector<double> &log_1_omega_k_rho_j_V,
                                 double Total)

{

                        // We compute p(Y|...)p(gamma|Omega)

    double log_posterior_term=0;
    double logMargLik=0;


    unsigned int nX=mat_X->size1;
    unsigned int pX=mat_X->size2;
    unsigned int pY=matY->size2;

//*********************************
//  Compute p(gamma|Omega)
//*********************************

    double logPGam=0;
    unsigned int pXGam=0;
    pXGam=gamma_list_vars.size(); // number of variables in the model


    // New algorithm
    if(HESS_slow_Postproc)
    {
        //Total is the sum of the logarithms
        logPGam=Total;
    }
    else
    {
        //Total is the sum of the rho_j's
        logPGam=-omega_k*Total;
    }


    //for (unsigned int j=0;j<pX;j++) //Old algorithm: we loop through all the vector
    // New loop: only the variables in are visited
    if (pXGam>0)
    {
        for (unsigned int count=0;count<pXGam;count++)
        {
            unsigned int j=gamma_list_vars[count];

            // 1) Additive term
            // ****************
            if (Memory_Limited)
            {
                logPGam=logPGam + log(omega_k * rho_j_V[j]);
            }
            else
            {
                // When we use pre-computed logs, we may run onto memory problems if p is very large...
                logPGam=logPGam + log_omega_k + log_rho_j_V[j];
            }

            // 2) Substracted term
            // *******************
            if (HESS_slow_Postproc) // Exact computations
            {
                if (Memory_Limited)
                {
                    // Old algorithm
                    //logPGam=logPGam + log(1-omega_k * rho_j_V[j]);

                    // New algorithm: substract from the whole total
                    logPGam=logPGam - log(1-omega_k * rho_j_V[j]);
                }
                else
                {
                    // Old algorithm
                    //logPGam=logPGam + log_1_omega_k_rho_j_V[j];
                    // New algorithm: substract from the whole total
                    logPGam=logPGam - log_1_omega_k_rho_j_V[j];
                }
            }
            else // Approximation log(1+u)=u+o(u) to speed up computations
            {
                // Old algorithm
                //logPGam=logPGam - omega_k * rho_j_V[j];

                // New algorithm: substract from the whole total
                logPGam=logPGam + omega_k * rho_j_V[j];
            }

        }
    }


    //Old algorithm:
//    for (unsigned int j=0;j<pX;j++)
//    {
//        if (Gamma[j]==1)
//        {
//            if (Memory_Limited)
//            {
//                logPGam=logPGam + log(omega_k * rho_j_V[j]);
//            }
//            else
//            {
//                // When we use pre-computed logs, we may run onto memory problems if p is very large...
//                logPGam=logPGam + log_omega_k + log_rho_j_V[j];
//            }
//            pXGam=pXGam+1;
//        }
//        else
//        {
//            if (HESS_slow_Postproc)
//            {
//                //debug:
//                //cout << "HESS_slow_Postproc = " << HESS_slow_Postproc << " when it should be true" << endl;
//                //
//                if (Memory_Limited)
//                {
//                    logPGam=logPGam + log(1-omega_k * rho_j_V[j]);
//                }
//                else
//                {
//                // When we use pre-computed logs, we may run onto memory problems if p is very large...
//                    logPGam=logPGam + log_1_omega_k_rho_j_V[j];
//                }
//            }
//            else
//            {
//                //cout << "HESS_slow_Postproc = " << HESS_slow_Postproc << " when it should be wrong" << endl;

//                // Approximation log(1+u)=u+o(u) to speed up computations
//                // Significant gain of performance
//                logPGam=logPGam - omega_k * rho_j_V[j];
//            }
//        }
//    }





//*************************************************
//  Compute the marginal likelihood: p(Y|...)
//*************************************************

    if (pY==1)
    {
        double SGamma=0.0;
        SGamma=1.0/(g+1.0)*gsl_matrix_get(matYTY,0,0)+g/(g+1.0)*gsl_matrix_get(R2_Mat_GPriors,0,0);

        double detQk=PR.k+SGamma;

        double logDetInnerMatrix=(double)pXGam*(log(1+g)-log(g));

        // THE FIRST TERM IN log(g) IS AN EXTRA TERM TO ACCOMODATE THE CODE BETWEEN
        // INDEPENDENT PRIORS AND THE OTHER PRIORS...
        // NOT EASY TO READ THOUGH, QUITE CONFUSING


        logMargLik=-((double)pY/2.0)*(double)(pXGam)*log(g)-
                    ((double)pY/2.0)*logDetInnerMatrix-
                    ((PR.delta+(double)nX+(double)pY-1.0)/2.0)*log(detQk);

    }
    else
    {

        //****************************************************************
        //  Construct matXGam: we don't need it in the case of g-priors
        //****************************************************************

        gsl_matrix *matXGam;//=gsl_matrix_calloc(nX,pXGam);

        //********************************
        //  Build matSGamma
        //********************************

        gsl_matrix *matSGamma=gsl_matrix_alloc(pY,pY);


        Quick_getSGamma_GPriors(matSGamma,
                                R2_Mat_GPriors,
                                matYTY,
                                nX,
                                pXGam,
                                pY,
                                g);

        //*************************************************
        //  Compute the marginal likelihood: p(Y|...)
        //*************************************************

            float *matrixEigenVecs,*vectorEigenVals; // Those are not used
            //matrixEigenVecs = (float*) malloc(pXGam*pXGam*sizeof(float));
            //vectorEigenVals = (float*) malloc(pXGam*sizeof(float));

            logMargLik=getLogMarg(PR,matXGam,matSGamma,matrixEigenVecs,vectorEigenVals,
                                    lambda,g,gPriorFlag,indepPriorFlag,pXGam,nX,pY);

            gsl_matrix_free(matSGamma);

            //gsl_matrix_free(matXGam);


    }




//*********************************
//  Gather results
//*********************************
    log_posterior_term=logMargLik+logPGam;

    //DEBUG
    //cout << "logPGam = " << logPGam << endl;
    //cout << "logMargLik = " << logMargLik << endl;
    //cout << "pXGam = " << pXGam << endl;
    //cout << "matYTY = " << gsl_matrix_get(matYTY,0,0) << endl;
    //cout << "R2 = " << gsl_matrix_get(R2_Mat_GPriors,0,0) << endl;
    //cout << "matSGamma = " << gsl_matrix_get(matSGamma,0,0) << endl;
    //cout << endl;


    return log_posterior_term;

}



void Test_Output_File(string Name,ofstream & File,ios_base::openmode fileMode)
{

    File.open(Name.c_str(),fileMode);
    if(File.fail())
    {
      cout << "Invalid Path and/or permission rights for " << Name << " -- run stopped." << endl;
      exit(-1);
    }
    else
    {
        if (Print==1)
        {
            cout << "File " << Name << " successefuly opened" << endl;
        }
    }

}


gsl_matrix * Read_Data_X(fstream & f,unsigned int &n,unsigned int &p)
{
    f >> n;
    f >> p;

    cout << "n = " << n << " observations" << endl;
    cout << "p = " << p << " predictors" << endl;

    gsl_matrix *mat=gsl_matrix_alloc(n,p);

    for(unsigned int i=0;i<n;i++)
    {
      for(unsigned int j=0;j<p;j++)
      {
        double tmp;
        f >> tmp;
        gsl_matrix_set(mat,i,j,tmp);
      }
    }
    cout << endl;

    return mat;
}

vector <gsl_matrix *> Read_Data_Y(fstream & f,unsigned int &n, unsigned int & q,unsigned int &Nb_Resp)
{
    //Read first two lines
    f >> n; //Nb of lines
    // DEBUG
    cout << "n = " << n << " observations" << endl;
    f >> q; // Nb of Columns
    // DEBUG
    cout << "q = " << q << " responses" << endl;

    // 3rd line:
    // list of nb of variables for each response, they must sum up to q
    unsigned int Count_Vars=0;
    vector <unsigned int > Size_Responses;
    cout << "Number of variables per response: " << endl;
    cout << "The test of data dimension can miss some errors...In twhich case results are meaningless !!" << endl;
    //Pb: The while condition can stop before the end of the line
    while (Count_Vars<q)
    {
        unsigned int q_this_resp;
        f >> q_this_resp;
        // DEBUG
        cout << q_this_resp << " ";

        Size_Responses.push_back(q_this_resp);
        Count_Vars=Count_Vars+q_this_resp;
    }
    // DEBUG
    cout << endl << endl;

    if (Count_Vars!=q)
    {
        // Problem: we haven't recovered the correct number of responses
        cout << "********************************" << endl;
        cout << "********************************" << endl;
        cout << "   Problem reading data Y" << endl;
        cout << "   Sum of variables per response does not add up" << endl;
        cout << "   to total number of variables " << endl;
        cout << "   Check first 3 lines of data Y file" << endl;
        cout << "********************************" << endl;
        cout << "********************************" << endl;
        exit(-1);
    }

    Nb_Resp=Size_Responses.size();

    // Now build the vector of gsl matrices to vectorize the data
    vector <gsl_matrix *> Vector_Data_Y;
    Vector_Data_Y.resize(Nb_Resp);
    for (unsigned int Response=0;Response<Nb_Resp;Response++)
    {
        Vector_Data_Y[Response]=gsl_matrix_alloc(n,Size_Responses[Response]);
    }

    // Read the rest of the file, and store it in the corresponding
    // place in the corresponding matrix!

    // read each line
    for(unsigned int i=0;i<n;i++)
    {
      unsigned int Curr_Resp=0;
      unsigned int Floor_count_vars=0;
      for(unsigned int j=0;j<q;j++)
      {
        double tmp;
        f >> tmp;

        // Find the correct Response and the correct variable inside this response,
        // and index it accordingly!
        if (j==Floor_count_vars+Size_Responses[Curr_Resp])
        {
            Floor_count_vars=Floor_count_vars+Size_Responses[Curr_Resp];
            Curr_Resp++;
        }
        unsigned int j_in_this_Resp=j-Floor_count_vars;

        // Record the data point
        gsl_matrix_set(Vector_Data_Y[Curr_Resp],i,j_in_this_Resp,tmp);
      }
    }

    return Vector_Data_Y;
}

// 2D version
void Read_File_Double(fstream & File, vector < vector <double > > & Mat)
{

    string Line;
    vector <double> Temp_Line;

    while ( getline(File, Line) ) //Reads the file line by line, and stores each line in Line
    {
        Temp_Line.resize(0);

        //Now we read each line Ligne into separate entries
        string Word;
        istringstream Read_Line(Line);
        while (Read_Line >> Word)
        {
            Temp_Line.push_back(atof(Word.c_str()));
        }
        Mat.push_back(Temp_Line);

    }
}

// 1D version
void Read_File_Double(fstream & File, vector <double > & Vec)
{

    string Line;

    while ( getline(File, Line) ) //Reads the file line by line, and stores each line in Line
    {
        //Now we read each line Ligne into separate entries
        string Word;
        istringstream Read_Line(Line);
        while (Read_Line >> Word)
        {
            Vec.push_back(atof(Word.c_str()));
        }

    }
}

void Read_File_Int(fstream & File, vector < vector <unsigned int > > & Mat)
{

    string Line;
    vector <unsigned int > Temp_Line;

    while ( getline(File, Line) ) //Reads the file line by line, and stores each line in Line
    {
        Temp_Line.resize(0);

        //Now we read each line Ligne into separate entries
        string Word;
        istringstream Read_Line(Line);
        while (Read_Line >> Word)
        {
            Temp_Line.push_back(atoi(Word.c_str()));
        }
        Mat.push_back(Temp_Line);

    }
}

// 1D version
void Read_File_Int(fstream & File, vector <unsigned int > & Vec)
{

    string Line;

    while ( getline(File, Line) ) //Reads the file line by line, and stores each line in Line
    {
        //Now we read each line Ligne into separate entries
        string Word;
        istringstream Read_Line(Line);
        while (Read_Line >> Word)
        {
            Vec.push_back(atoi(Word.c_str()));
        }

    }
}


void Write_Vector(ostream & File,
                  vector < unsigned int > & Vector)
{
    unsigned int Nb_Lines=Vector.size();
    for (unsigned int line=0;line<Nb_Lines;line++)
    {
        File << Vector[line]<< endl;
    }

}

void Write_Vector(ostream & File,
                  vector < double > & Vector)
{
    unsigned int Nb_Lines=Vector.size();
    for (unsigned int line=0;line<Nb_Lines;line++)
    {
        File << Vector[line]<< endl;
    }

}


void Write_Matrix(ostream & File,
                  vector <vector < double > > & Matrix,
                  string Field_Separator)
{
    unsigned int Nb_Lines=Matrix.size();
    for (unsigned int line=0;line<Nb_Lines;line++)
    {
        unsigned int Nb_Col=Matrix[line].size();
        for (unsigned int col=0;col<Nb_Col;col++)
        {
            File << Matrix[line][col] << Field_Separator;
        }
        File << endl;
    }

}



//void Mean_Matrix(vector <vector < double > > & Matrix,
//                 vector <double> &Mean,
//                 unsigned int burnin,
//                 unsigned int nsweep)
//{
//    for (unsigned int col=0;col<Mean.size();col++)
//    {   Mean[col]=0;}

//    for (unsigned int sweep=burnin;sweep<nsweep;sweep++)
//    {
//        unsigned int Nb_Col=Matrix[sweep].size();
//        for (unsigned int col=0;col<Nb_Col;col++)
//        {
//            Mean[col]=Mean[col]+Matrix[sweep][col]/(double)(nsweep-burnin);
//        }
//    }


//}

void Write_Matrix(ostream & File,
                  vector <vector < unsigned int > > & Matrix,
                  string Field_Separator)
{
    unsigned int Nb_Lines=Matrix.size();
    for (unsigned int line=0;line<Nb_Lines;line++)
    {
        unsigned int Nb_Col=Matrix[line].size();
        for (unsigned int col=0;col<Nb_Col;col++)
        {
            File << Matrix[line][col] << Field_Separator;
        }
        File << endl;
    }


}


void Write_1Line(ostream & File,
                  vector < double > & Line,
                  string Field_Separator)
{

    unsigned int Nb_Col=Line.size();
    for (unsigned int col=0;col<Nb_Col;col++)
    {
      File << Line[col] << Field_Separator;
    }
    File << endl;


}


void Write_1Line(ostream & File,
                  vector < unsigned int > & Line,
                  string Field_Separator)
{

    unsigned int Nb_Col=Line.size();
    for (unsigned int col=0;col<Nb_Col;col++)
    {
      File << Line[col] << Field_Separator;
    }
    File << endl;


}


void Show_Time(double remainingTime)
{
        double Time_Hr=0;
        double Time_Mn=0;
        double Time_s=0;

        double Residue=remainingTime;
        Time_Hr=floor(Residue/3600.0);
        Residue=Residue-3600.0*Time_Hr;

        Time_Mn=floor(Residue/60.0);
        Residue=Residue-60.0*Time_Mn;

        Time_s=Residue;

        if (Time_Hr==0)
        {
            if (Time_Mn==0)
            {
                if (Time_s==0)
                {
                    cout << " 0 s";
                }
                else
                {
                    cout << setprecision(2) << Time_s << " s";
                }
            }
            else
            {
                cout << Time_Mn << " mn " << Time_s << " s";
            }
        }
        else
        {
            cout << Time_Hr << " hr " << Time_Mn << " mn " << Time_s << " s";

        }


}



void Display_Matrices(vector < vector < unsigned int > > M)
{
    unsigned int Nb_l=M.size();
    for (unsigned int l=0;l<Nb_l;l++)
    {
        cout << "Line " << l << endl;
        unsigned int Nb_c=M[l].size();

        for (unsigned int c=0;c<Nb_c;c++)
        {
            cout << "\t" << M[l][c];
            //cout << "\t["<< l << "," << c << "]";
        }
        cout << endl;
    }
    cout << endl;

}


void Display_Matrices(vector < vector < double > > M)
{
    vector < vector < double > > *p_M=& M;
    Display_Matrices(p_M);
}

void Display_Matrices(vector < vector < double > > * M)
{
    unsigned int Nb_l=M->size();
    for (unsigned int l=0;l<Nb_l;l++)
    {
        cout << "Line " << l << endl;
        unsigned int Nb_c=(*M)[l].size();

        for (unsigned int c=0;c<Nb_c;c++)
        {
            cout << "\t" << (*M)[l][c];
            //cout << "\t["<< l << "," << c << "]";
        }
        cout << endl;
    }
    cout << endl;

}

void Display_Vector(vector < double > M)
{
    vector < double > *p_M=& M;
    Display_Vector(p_M);
}

void Display_Vector(vector < double > * M)
{
    unsigned int Nb_l=M->size();
    for (unsigned int l=0;l<Nb_l;l++)
    {
        cout << "\t" << (*M)[l];
        cout << endl;
    }
    //cout << endl;
}

void Display_Vector(vector < unsigned int > M)
{
    vector < unsigned int > *p_M=& M;
    Display_Vector(p_M);
}


void Display_Vector(vector < unsigned int > * M)
{
    unsigned int Nb_l=M->size();
    for (unsigned int l=0;l<Nb_l;l++)
    {
        cout << "\t" << (*M)[l];
        cout << endl;
    }
    //cout << endl;
}
