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

#ifndef COMMAND_LINE_H
#define COMMAND_LINE_H

#include <algorithm>
#include <string>

using namespace std;


class Command_Line
{
public:
    Command_Line();
    void Process_Comm_Line(int argc, char* argv[]);

    string Model_Tag;
    string Regression_Model;
    bool Model_Flag;

    unsigned int n_sweeps;
    unsigned int burn_in;

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

    bool gPriorFlag;
    bool indepPriorFlag;
    double lambda;

    bool gSampleFlag;
    bool gInitFlag;
    bool gamSampleFlag;
    bool OmegaKSampleFlag;
    bool RhoJSampleFlag;


    bool iso_T_Flag;

    unsigned int nConfounders;
    unsigned int n_top_models;

    long MY_SEED;
    bool standardizeFlag;
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

    //Parameters for the priors of the Matrix Omega in HESS:
    double a_om_k;
    bool a_om_k_Flag;

    double b_om_k;
    bool b_om_k_Flag;

    double c_rh_j;
    bool c_rh_j_Flag;

    double d_rh_j;
    bool d_rh_j_Flag;

    double Prop_Active;

    bool HESS_do_Postproc;
    bool HESS_slow_Postproc;
    bool Memory_Limited;
    bool Single_g;

    // For post-processing
    bool Include_Single_Models;





};

#endif // COMMAND_LINE_H
