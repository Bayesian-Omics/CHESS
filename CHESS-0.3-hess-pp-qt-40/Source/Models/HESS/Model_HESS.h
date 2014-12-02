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

#ifndef MODEL_HESS_H
#define MODEL_HESS_H

#include "../../General_Classes/Model_Generic.h"

#include <algorithm>
#include <string>
#include <vector>

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


using namespace std;

class Model_HESS : public Model_Generic
{
public:
    Model_HESS();
    //Inherited constructor:
    Model_HESS(Command_Line My_Command_Line) : Model_Generic(My_Command_Line){}


    /*_______________
        Overloaded:
      ---------------*/
    void Run();

};

void sample_g_HESS(unsigned int n_chains,
          //AdMH *My_g_AdMH,
          vector<AdMH *> My_g_AdMH,
          vector <Kernel_Single_Gamma> Gammas,
          gsl_matrix *mat_X,
          vector <gsl_matrix *> &Vector_Data_Y,
          bool gPriorFlag,
          bool indepPriorFlag,
          bool gSampleFlag,
          double lambda,
          //double &g,
          vector<double> &g,
          vector<Prior_param> PR_per_Resp,
          unsigned int sweep,
              bool cudaFlag,
              gsl_rng *RandomNumberGenerator,
              Model_Information &Model_Data,
          bool Single_g);

void Set_Omega_k_AdMH(vector< vector < AdMH *> > Omega_k_AdMH,
                      unsigned int Nb_Resp,
                      unsigned int nb_chains,
                      unsigned int burn_in,
                      double g_AdMH_ls_from_read,
                      unsigned int g_n_batch_from_read,
                      double g_AdMH_optimal_from_read);

void Set_Rho_j_AdMH(vector< vector < AdMH *> > Rho_j_AdMH,
                    unsigned int pX,
                    unsigned int nb_chains,
                    unsigned int burn_in,
                    double g_AdMH_ls_from_read,
                    unsigned int g_n_batch_from_read,
                    double g_AdMH_optimal_from_read);


void Sample_Omega_k(unsigned int sweep,
                    gsl_rng *RandomNumberGenerator,
                    vector< vector <AdMH *> > Omega_k_AdMH,
                    double a_om_k,
                    double b_om_k,
                    unsigned int nb_chains,
                    unsigned int pX,
                    unsigned int Nb_Resp,
                    vector < vector < double > > &Omega_PerLine,
                    vector<vector<double> > &Rho_PerCol,
                    vector<Kernel_Single_Gamma> &Gammas,
                    Model_Information &Model_Data);

void Sample_Rho_j(unsigned int sweep,
                  gsl_rng *RandomNumberGenerator,
                  vector< vector <AdMH *> > Rho_j_AdMH,
                  double c_rh_j,
                  double d_rh_j,
                  unsigned int nb_chains,
                  unsigned int pX,
                  unsigned int Nb_Resp,
                  vector<vector<double> > &Omega_PerLine,
                  vector<vector<double> > &Rho_PerCol,
                  vector<Kernel_Single_Gamma> &Gammas,
                  Model_Information &Model_Data);

double Logistic_Trans(double Omega);
double Logistic_Trans_Inv(double sample_tmp);


double Log_Pr_Omega_k(double x,
                      double a,
                      double b,
                      vector<Kernel_Single_Gamma > &Gammas,
                      unsigned int Chain,
                      unsigned int Response,
                      Model_Information Model_Data);

double Log_Pr_Rho_j(double x,
             double c,
             double d,
             vector<double> * Omega_V,
             vector<Kernel_Single_Gamma > &Gammas,
             unsigned int Chain,
             unsigned int Variable,
             unsigned int Nb_Resp,
             Model_Information Model_Data);

void Update_Active_Yk(vector < unsigned int > &Active_Yk,
                      double Proportion,
                      unsigned int Nb_Resp,
                      gsl_rng * RandomNumberGenerator);

void Preliminary_Tests(unsigned int nX,
                       unsigned int nY,
                       unsigned int Nb_Resp,
                       unsigned int n_sweeps,
                       unsigned int burn_in,
                       bool gamSampleFlag,
                       bool fixInit,
                       bool OmegaKSampleFlag,
                       bool fixInit_omega_k,
                       bool RhoJSampleFlag, bool fixInit_rho_j);

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
                                 );

void Mean_Matrix(vector <vector < double > > & Matrix,
                 vector <double> & Mean,
                 unsigned int burnin,
                 unsigned int nsweep);


void Save_Rho_History(vector <vector < double > > & History_Rho,
                      vector<vector<double> > Rho_PerCol,
                      unsigned int sweep);

void Save_Omega_History(vector<vector<double> > &History_Rho,
                      vector<vector<double> > Rho_PerCol,
                      unsigned int sweep);


void Save_Average_Omega_k_j(vector<vector<double> > &Average_Omega_k_j,
                            vector<vector<double> > Omega_PerLine,
                            vector<vector<double> > Rho_PerCol,
                            unsigned int sweep,
                            unsigned int burn_in);

void Compute_Average_Histo(vector<double> &Average,
                            vector< vector< double> > Histo,
                            unsigned int burn_in);

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
                                 double Total);

void Test_Output_File(string Name,ofstream & File,ios_base::openmode fileMode);


gsl_matrix * Read_Data_X(fstream & f,
                         unsigned int &n,
                         unsigned int &p);

vector<gsl_matrix *> Read_Data_Y(fstream & f,
                         unsigned int &n,
                         unsigned int &q,
                         unsigned int &Nb_Resp);

void Read_File_Double(fstream & File,
                      vector<vector<double> > &Mat);
void Read_File_Double(fstream & File,
                      vector <double > & Vec);

void Read_File_Int(fstream & File,
                   vector < vector <unsigned int > > & Mat);

void Read_File_Int(fstream & File, vector <unsigned int > & Vec);

void Write_Vector(ostream & File,
                  vector < unsigned int > & Vector);

void Write_Vector(ostream & File,
                  vector < double > & Vector);

void Write_Matrix(ostream & File,
                  vector<vector<double> > &Matrix,
                  string Field_Separator);

void Write_Matrix(ostream & File,
                  vector<vector<unsigned int> > &Matrix,
                  string Field_Separator);

void Write_1Line(ostream & File,
                  vector < double > & Line,
                  string Field_Separator);

void Write_1Line(ostream & File,
                  vector < unsigned int > & Line,
                  string Field_Separator);

void Show_Time(double remainingTime);


void Display_Matrices(vector < vector < unsigned int > > M);

void Display_Matrices(vector < vector < double > > * M);
void Display_Matrices(vector < vector < double > > M);

void Display_Vector(vector < double > * M);
void Display_Vector(vector < double > M);

void Display_Vector(vector < unsigned int > M);
void Display_Vector(vector < unsigned int > * M);


#endif // MODEL_HESS_H
