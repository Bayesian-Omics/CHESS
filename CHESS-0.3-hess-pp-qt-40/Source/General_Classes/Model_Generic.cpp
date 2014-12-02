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

#include "Model_Generic.h"
#include "Command_Line.h"

#include <iostream>

Model_Generic::Model_Generic()
{
    Model_Tag="Generic";
    Regression_Model="Undefined";

}

Model_Generic::Model_Generic(Command_Line My_Command_Line)
{

    Model_Tag=My_Command_Line.Model_Tag;
    Regression_Model=My_Command_Line.Regression_Model;

    n_sweeps=My_Command_Line.n_sweeps;
    burn_in=My_Command_Line.burn_in;

    MY_SEED=My_Command_Line.MY_SEED;
    postProcessOnly=My_Command_Line.postProcessOnly;
    resumeRun=My_Command_Line.resumeRun;
    timeLimit=My_Command_Line.timeLimit;
    Time_monitorFlag=My_Command_Line.Time_monitorFlag;
    Out_full_Flag=My_Command_Line.Out_full_Flag;
    HistoryFlag=My_Command_Line.HistoryFlag;
    cudaFlag=My_Command_Line.cudaFlag;
    Log_Flag=My_Command_Line.Log_Flag;

    fixInit=My_Command_Line.fixInit;
    fixInit_rho_j=My_Command_Line.fixInit_rho_j;
    fixInit_omega_k=My_Command_Line.fixInit_omega_k;


    iso_T_Flag=My_Command_Line.iso_T_Flag;
    nConfounders=My_Command_Line.nConfounders;
    n_top_models=My_Command_Line.n_top_models;

    standardizeFlag=My_Command_Line.standardizeFlag;
    gPriorFlag=My_Command_Line.gPriorFlag;
    indepPriorFlag=My_Command_Line.indepPriorFlag;

    gSampleFlag=My_Command_Line.gSampleFlag;
    gInitFlag=My_Command_Line.gInitFlag;;
    gamSampleFlag=My_Command_Line.gamSampleFlag;
    OmegaKSampleFlag=My_Command_Line.OmegaKSampleFlag;
    RhoJSampleFlag=My_Command_Line.RhoJSampleFlag;


    lambda=My_Command_Line.lambda;    
    doYShuffle=My_Command_Line.doYShuffle;

    path_name_out=My_Command_Line.path_name_out;
    filename_init=My_Command_Line.filename_init;
    filename_in_mat_X=My_Command_Line.filename_in_mat_X;
    filename_in_mat_Y=My_Command_Line.filename_in_mat_Y;
    filename_par=My_Command_Line.filename_par;

    g_init=My_Command_Line.g_init;
    Pr_g_a=My_Command_Line.Pr_g_a;
    Pr_g_b=My_Command_Line.Pr_g_b;
    Pr_g_a_flag=My_Command_Line.Pr_g_a_flag;
    Pr_g_b_flag=My_Command_Line.Pr_g_b_flag;


    filename_init_omega_k=My_Command_Line.filename_init_omega_k;
    filename_init_rho_j=My_Command_Line.filename_init_rho_j;


    a_om_k=My_Command_Line.a_om_k;
    b_om_k=My_Command_Line.b_om_k;
    c_rh_j=My_Command_Line.c_rh_j;
    d_rh_j=My_Command_Line.d_rh_j;

    a_om_k_Flag=My_Command_Line.a_om_k_Flag;
    b_om_k_Flag=My_Command_Line.b_om_k_Flag;
    c_rh_j_Flag=My_Command_Line.c_rh_j_Flag;
    d_rh_j_Flag=My_Command_Line.d_rh_j_Flag;

    Prop_Active=My_Command_Line.Prop_Active;

    HESS_do_Postproc=My_Command_Line.HESS_do_Postproc;
    HESS_slow_Postproc=My_Command_Line.HESS_slow_Postproc;
    Memory_Limited=My_Command_Line.Memory_Limited;
    Single_g=My_Command_Line.Single_g;

    Include_Single_Models=My_Command_Line.Include_Single_Models;


}
