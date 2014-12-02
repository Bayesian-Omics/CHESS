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

#ifndef START_UP_ESS_H
#define START_UP_ESS_H


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

using namespace std;

#if _CUDA_
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <cula.h>
#endif

#define DEBUG 0

class Start_Up_ESS
{
public:
    Start_Up_ESS();


    gsl_matrix * mat_X_work2;
    gsl_matrix * mat_Y_work2;

    gsl_vector * vect_RMSE;
    AdMH * My_g_AdMH;

    unsigned int g_n_batch_from_read;
    double g_AdMH_optimal_from_read;
    double g_AdMH_ls_from_read;
    double g_M_min_input;
    double g_M_max_input;


    void Initialise(   
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
                    string filename_par);


};

#endif // START_UP_ESS_H
