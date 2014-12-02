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

#include "Command_Line.h"

#include <iostream>
//#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula.h>
#endif


using namespace std;

Command_Line::Command_Line()
{
    Model_Tag="";
    Regression_Model="Linear"; // Default setting
    Model_Flag=false;

    n_sweeps=0;
    burn_in=0;

    postProcessOnly=false;
    resumeRun=false;
    timeLimit=false;
    Time_monitorFlag=false;
    Out_full_Flag=false;
    HistoryFlag=false;
    cudaFlag=false;
    Log_Flag=false;

    fixInit=false;
    fixInit_omega_k=false;
    fixInit_rho_j=false;

    gPriorFlag=false;
    indepPriorFlag=false;
    gSampleFlag=true;
    gamSampleFlag=false;
    lambda=0;

    iso_T_Flag=false;

    nConfounders=0;
    n_top_models=0;

    MY_SEED=0;
    standardizeFlag=false;
    g_init=0;
    Pr_g_a=0.0;
    Pr_g_b=0.0;
    Pr_g_a_flag=false;
    Pr_g_b_flag=false;

    doYShuffle=false;

    a_om_k_Flag=false;
    b_om_k_Flag=false;
    c_rh_j_Flag=false;
    d_rh_j_Flag=false;

    Prop_Active=.25; // Default setting

    //HESS_do_Postproc=true; // OLD Default setting
    HESS_do_Postproc=false; // NEW Default setting
    HESS_slow_Postproc=false; // Default setting
    Memory_Limited=false; // Default setting
    Single_g=true; //Default setting

    //For post-processing in HESS
    Include_Single_Models=true;

}


void Command_Line::Process_Comm_Line(int argc, char* argv[])
{

    int na=0;
    MY_SEED=-1;
    n_sweeps=0;

    burn_in=0;
    unsigned int n_top_models_from_read=100;
    n_top_models=0;
    nConfounders=0;

    // If gSampleFlag we are also sampling g
    gSampleFlag=true;
    gInitFlag=false;
    gamSampleFlag=true;

    OmegaKSampleFlag=true;
    RhoJSampleFlag=true;


    // If gPriorFlag we are using g-priors otherwise powered priors
    // If indepPriorFlag we are using independent priors
    // Even though g prior and independent prior are special cases of powered
    // priors, computations can be done more cheaply in those special cases

    gPriorFlag=true;
    indepPriorFlag=false;
    lambda=-1.0;
    g_init=1.0;
    Pr_g_a=0.0;
    Pr_g_b=0.0;
    Pr_g_a_flag=false;
    Pr_g_b_flag=false;

    standardizeFlag=true;
    resumeRun=false;
    doYShuffle=false;
    HistoryFlag=false;
    Time_monitorFlag=false;

    bool X_Flag=false;
    bool Y_Flag=false;
    bool Par_Flag=false;
    Out_full_Flag=false;
    Log_Flag=false;
    iso_T_Flag=false;
    cudaFlag=false;
    fixInit=false;
    postProcessOnly=false;

    timeLimit = -1;

    na++;
    while(na < argc)
      {
        if ( 0 == strcmp(argv[na],"-model") )
        {
            Model_Flag=true;
            Model_Tag=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-regression") )
        {
            Regression_Model=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-X") )
        {
            X_Flag=true;
            filename_in_mat_X=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-Y") )
        {
            Y_Flag=true;
            filename_in_mat_Y=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-par") )
        {
            Par_Flag=true;
            filename_par=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-nsweep") )
        {
            n_sweeps=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-burn_in") )
        {
            burn_in=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-seed") )
        {
            MY_SEED=(long)((atoi(argv[++na])));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-out") )
        {
            path_name_out=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-out_full") )
        {
            Out_full_Flag=true;
            path_name_out=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-history") )
        {
            HistoryFlag=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-time") )
        {
            Time_monitorFlag=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-top") )
        {
            n_top_models_from_read=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-log") )
        {
            Log_Flag=true;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-isoT") )
        {
            iso_T_Flag=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-noRescaleX") )
        {
            standardizeFlag=false;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-lambda") )
        {
            gPriorFlag=false;
            lambda=(double)(atof(argv[++na]));
            if(fabs(lambda)<0.000000000001)
            {
              indepPriorFlag=true;
            }
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-g_set") )
        {
            gSampleFlag=false;
            g_init=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-g_init") )
        {
            gInitFlag=true;
            g_init=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-pr_g_a") )
        {
            Pr_g_a_flag=true;
            Pr_g_a=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-pr_g_b") )
        {
            Pr_g_b_flag=true;
            Pr_g_b=(double)(atof(argv[++na]));
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-gam_set") )
        {
            gamSampleFlag=false;
            fixInit=true;
            filename_init=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-init") )
        {
            fixInit = true;
            filename_init=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-init_omega_k") )
        {
            fixInit_omega_k = true;
            filename_init_omega_k=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-omega_k_set") )
        {
            OmegaKSampleFlag=false;
            fixInit_omega_k=true;
            filename_init_omega_k=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-init_rho_j") )
        {
            fixInit_rho_j = true;
            filename_init_rho_j=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-rho_j_set") )
        {
            RhoJSampleFlag=false;
            fixInit_rho_j=true;
            filename_init_rho_j=argv[++na];
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-a_om_k") )
        {
            a_om_k_Flag=true;
            a_om_k=atof(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-b_om_k") )
        {
            b_om_k_Flag=true;
            b_om_k=atof(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-c_rh_j") )
        {
            c_rh_j_Flag=true;
            c_rh_j=atof(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-d_rh_j") )
        {
            d_rh_j_Flag=true;
            d_rh_j=atof(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-Prop_hess") )
        {
            Prop_Active=atof(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-DO_HESS_Postproc") )
        {
            //NEW DEFAULT
            HESS_do_Postproc=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-HESS_slow_Postproc") )
        {
            HESS_slow_Postproc=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-Memory_Limited") )
        {
            Memory_Limited=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-Multiple_g") )
        {
            Single_g=false;
            cout << "Using multiple independent g" << endl;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-No_Single_Models") )
        {
            Include_Single_Models=false;
            cout << "Excluding single models from post-processing" << endl;
            if (na+1==argc) break;
            na++;
        }

        else if ( 0 == strcmp(argv[na],"-cuda") )
        {
  #if _CUDA_
            cudaFlag=true;

            int nDevices;
            culaStatus culaStat;
            culaStat = culaGetDeviceCount(&nDevices);
            if(culaStat!=culaNoError){
              cout << "Error detecting how many devices" << endl;
              exit(-1);
            }else{
              cout << "Detected " << nDevices << " gpu devices" << endl;
            }

            cudaThreadExit();
            for(int i=0;i<nDevices;i++){
              culaStat=culaSelectDevice(i);
              if(culaStat!=culaNoError){
                cout << culaGetStatusString(culaStat) << endl;
                cout << culaGetErrorInfo() << endl;
                cudaThreadExit();
                continue;
              }

              culaStat=culaInitialize();
              if(culaStat!=culaNoError){
                cout << culaGetStatusString(culaStat) << endl;
                cout << culaGetErrorInfo() << endl;
                cublasShutdown();
                culaShutdown();
                cudaThreadExit();
                continue;
              }

              break;
            }

            if(culaStat!=culaNoError){
              cout << "Unable to initialise GPU" << endl;
              cublasShutdown();
              culaShutdown();
              exit(1);
            }else{
              int device;
              culaGetExecutingDevice(&device);
              cout << "Successfully initialised GPU" << endl;
              cout << "Using device " << device << endl;
            }

  #else
            cudaFlag=false;
  #endif
            if (na+1==argc) break;
            na++;

        }
        else if ( 0 == strcmp(argv[na],"-nconf") )
        {// Confounders must be in the first columns of X
            nConfounders=atoi(argv[++na]);
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-resume") )
        {// Confounders must be in the first columns of X
            resumeRun=true;
            if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-shuffle") )
        {// Confounders must be in the first columns of X
          doYShuffle=true;
          if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-postProcess") )
        {// Confounders must be in the first columns of X
          postProcessOnly=true;
          if (na+1==argc) break;
            na++;
        }
        else if ( 0 == strcmp(argv[na],"-timeLimit") )
        {// Confounders must be in the first columns of X
          timeLimit=atof(argv[++na]);
          if (na+1==argc) break;
            na++;
        }
        else
        {
            cout << "Unknown option: " << argv[na] << endl;
  #if _CUDA_
            if(cudaFlag){
                cublasShutdown();
                culaShutdown();
            }
  #endif
            exit(1);
        }
    }

//    cout << "Command_Line::Process... HESS_slow_Postproc = " << HESS_slow_Postproc << endl;
//    cout << "Pr_g_a = " << Pr_g_a << endl;
//    cout << "Pr_g_b = " << Pr_g_b << endl;
//    cout << "Pr_g_a_flag = " << Pr_g_a_flag << endl;
//    cout << "Pr_g_b_flag = " << Pr_g_b_flag << endl;

    cout << endl;

    //cout << "Single_g = " << Single_g << endl << endl;

    if(!X_Flag){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "   The predictor matrix X has not been specified, RUN STOPPED" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
          cublasShutdown();
          culaShutdown();
      }
  #endif

      exit(1);
    }

    if(!Model_Flag){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "   Model has not been specified, RUN STOPPED" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }


    if(!Y_Flag){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "   The outcome matrix Y has not been specified, RUN STOPPED" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }

    if(!Par_Flag){
      cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "   The parameters matrix has not been specified, RUN STOPPED" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
       << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  #if _CUDA_
      if(cudaFlag){
        cublasShutdown();
        culaShutdown();
      }
  #endif
      exit(-1);
    }

    if(postProcessOnly){
      resumeRun=false;
      fixInit=false;
    }

    if(n_top_models_from_read!=0){
      n_top_models=n_top_models_from_read;
    }
    else
    {
        n_top_models=n_sweeps;
    }

}
