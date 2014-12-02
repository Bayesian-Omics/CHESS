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

#ifndef KERNEL_SINGLE_GAMMA_H
#define KERNEL_SINGLE_GAMMA_H


#include <vector>

#include "../Kernel/Routines/matrix_handling.h"
#include "../Kernel/Routines/rand.h"


#include "../Kernel/Routines/cond_post.h"

#include "../Kernel/Classes/Double_Matrices.h"

#include "../Kernel/Classes/Temperatures.h"
#include "../Kernel/Classes/Move_monitor.h"

#include "../Kernel/Classes/DR.h"
#include "../Kernel/Classes/CM.h"


using namespace std;

//class

class Kernel_Single_Gamma
{
public:
    Kernel_Single_Gamma();
    Kernel_Single_Gamma(bool postProcessOnly,
                        unsigned int nb_chains,
                        bool resumeRun,
                        unsigned int resumeDRNCalls,
                        unsigned int resumeDRNCallsAdj,
                        vector <vector <unsigned int> > resumeDRAccepted,
                        vector <vector <unsigned int> > resumeDRProposed,
                        unsigned int k_max_from_read,
                        unsigned int pX,
                        unsigned int n_sweeps,
                        Model_Information &Model_Data);

    void free();


    void Initialise(bool postProcessOnly,
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
                    Model_Information &Model_Data);

    void First_Sweep(unsigned int Response,
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
                     Model_Information &Model_Data);

    void MCMC_Update(unsigned int Response,
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
                     Model_Information &Model_Data);

    void Stepwise_Regression(vector < vector <unsigned int> > & Gam_step_regr,
                              //unsigned int nX,
                              //unsigned int pX,
                              //unsigned int Nb_Outcomes,
                              unsigned int nConfounders,
                              gsl_matrix * mat_X,
                              gsl_matrix * mat_Y,
                              gsl_vector *vect_RMSE,
                              double Pvalue_enter,
                              double Pvalue_remove,
                              Model_Information &Model_Data);

    void save_model_per_sweep();


    string Regression_Model;

    unsigned int Current_Sweep;
    // for the moment, used only in HESS to store data for the post-processing
    // see Kernel_Single_Gamma::save_model_per_sweep()



    vector < vector <unsigned int> > vect_gam;
    Temperatures *t_tun;
    Double_Matrices mat_log_marg;
    Double_Matrices mat_log_cond_post;

    gsl_permutation *MyPerm;
    vector < unsigned int > chain_idx;

    /*
      VERY IMPORTANT:
      ==============

        The position of the chain c for response k is Gammas[k]->chain_idx[c]
        --> This is important if we want to access gamma
    */


    vector<vector<unsigned int> > pastModels;
    vector<unsigned int> pastNModelsVisited;

    Move_monitor * My_Move_monitor;
    CM * My_CM;
    DR *My_DR;
    gsl_matrix *description_exchange_move;

    vector <unsigned int> n_Models_visited;
    vector < vector <unsigned int> > List_Models;

    vector < vector <unsigned int> > Unique_List_Models;

    vector<unsigned int> resumeChainIndex;
    gsl_matrix *matXGam;


};

#endif // KERNEL_SINGLE_GAMMA_H
