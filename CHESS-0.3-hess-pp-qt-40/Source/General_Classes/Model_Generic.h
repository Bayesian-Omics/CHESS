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

#ifndef MODEL_GENERIC_H
#define MODEL_GENERIC_H

#include <algorithm>
#include <string>

#include "Command_Line.h"

using namespace std;

class Model_Generic
{
public:
    Model_Generic();
    Model_Generic(Command_Line My_Command_Line);

    string Model_Tag;
    string Regression_Model;

    unsigned int n_sweeps;
    unsigned int burn_in;
    long MY_SEED;
    bool postProcessOnly;
    bool resumeRun;
    double timeLimit;
    bool Time_monitorFlag;
    bool Out_full_Flag;
    bool HistoryFlag;
    bool cudaFlag;
    bool Log_Flag;

    bool fixInit;
    bool fixInit_rho_j;
    bool fixInit_omega_k;

    bool iso_T_Flag;
    unsigned int nConfounders;
    unsigned int n_top_models;
    bool standardizeFlag;
    bool gPriorFlag;
    bool indepPriorFlag;
    bool gSampleFlag;
    bool gInitFlag;
    bool gamSampleFlag;
    bool OmegaKSampleFlag;
    bool RhoJSampleFlag;

    double lambda;
    bool doYShuffle;

    string path_name_out;
    string filename_init;
    string filename_in_mat_X;
    string filename_in_mat_Y;
    string filename_par;

    double g_init;
    double Pr_g_a;
    double Pr_g_b;
    bool Pr_g_a_flag;
    bool Pr_g_b_flag;


    string filename_init_omega_k;
    string filename_init_rho_j;

    //Post-processing option:
    bool Include_Single_Models;

    // Parameters for the Omega Prior in HESS
    double a_om_k;
    bool a_om_k_Flag;

    double b_om_k;
    bool b_om_k_Flag;

    double c_rh_j;
    bool c_rh_j_Flag;

    double d_rh_j;
    bool d_rh_j_Flag;

    // Run parameter for HESS
    double Prop_Active;

    bool HESS_do_Postproc;
    bool HESS_slow_Postproc;
    bool Memory_Limited;
    bool Single_g;



    void Run();

};

#endif // MODEL_GENERIC_H
