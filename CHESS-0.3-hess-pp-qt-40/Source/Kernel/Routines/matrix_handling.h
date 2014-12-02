/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from matrix_handling.h in the ESS++ program
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

#ifndef MATRIX_HANDLING_H
#define MATRIX_HANDLING_H

#include <iostream>
#include <iomanip>
#include <sstream> 
#include <fstream> 
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include "../Classes/String_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Double_Matrices.h"
#include "../Classes/Int_Matrices.h"
#include "../Classes/Temperatures.h"
#include "../Classes/DR.h"
#include "../Routines/rand.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>
#include <vector>
#include <algorithm>

using namespace std; 

void standardize_matrix_gsl(gsl_matrix* M);

gsl_matrix* Double_matrices_cont_2_gsl_matrix(Double_Matrices_cont source);

gsl_matrix* Double_matrices_2_gsl_matrix(Double_Matrices source);

unsigned int sum_vector_int(vector <unsigned int> &Myvector);

void get_list_var_in_and_out(vector <unsigned int> &list_columns_var_in,
			     vector <unsigned int> &list_columns_var_out,
			     vector <unsigned int> &is_var_in);

gsl_matrix* get_X_gam(vector <unsigned int> &list_columns_var_in,
		      gsl_matrix *mat_X);

gsl_matrix* get_X_reduced(vector <unsigned int> &list_columns_var,
			    gsl_matrix *mat_X);

gsl_matrix* get_X_reduced_and_constant(vector <unsigned int> &list_columns_var,
				       gsl_matrix *mat_X);

gsl_matrix* get_sub_matrix_col(gsl_matrix *mat_X,
			       size_t first_col,
			       size_t last_col);

gsl_matrix* get_sub_matrix_row(gsl_matrix *mat_X,
			       size_t first_row,
			       size_t last_row);

void display_gsl_matrix(gsl_matrix *M);

void display_gsl_vector(gsl_vector *V);

void display_gsl_perm(gsl_permutation *P);

void fill_sub_matrix_row(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_row,
			 size_t last_row);

void fill_sub_matrix_col(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_col,
			 size_t last_col);


void display_matrix_var_dim(vector < vector <unsigned int> > &M);

void center_matrix_gsl(gsl_matrix *M);

void display_vector_int(vector < unsigned int> &vector);

unsigned int sum_line_std_mat(vector < vector <unsigned int> > &M,
			      unsigned int line);

void get_list_var_in(vector <unsigned int> &list_columns_var_in,
		     vector <unsigned int> &is_var_in);

void display_result_per_sweep(vector < vector <unsigned int> > vect_gam,
			     vector < unsigned int > chain_idx,
			     Double_Matrices mat_log_marg,
			     Double_Matrices mat_log_cond_post,
			     unsigned int sweep,
			     Temperatures *t_tun);

void display_summary_result_per_sweep(vector < vector <unsigned int> > vect_gam,
				     vector < unsigned int > chain_idx,
				     Double_Matrices mat_log_marg,
				     Double_Matrices mat_log_cond_post,
				     unsigned int sweep,
				     Temperatures *t_tun,
				     const unsigned int& nConfounders);

void print_main_results_per_sweep(ofstream &f_out,
				 vector < vector <unsigned int> > vect_gam,
				 vector < unsigned int > chain_idx,
				 Double_Matrices mat_log_marg,
				 Double_Matrices mat_log_cond_post,
				 unsigned int sweep);

void print_and_save_main_results_per_sweep(ofstream &f_out,
					  ofstream &f_out_n_vars_in,
                                          ofstream &f_out_n_models_visited,
					  ofstream &f_out_log_cond_post,
					  vector < vector <unsigned int> > vect_gam,
					  vector < vector <unsigned int> > &List_models,
					  vector < unsigned int > chain_idx,
					  Double_Matrices mat_log_marg,
					  Double_Matrices mat_log_cond_post,
					  unsigned int sweep,
					  unsigned int burn_in,
					  unsigned int nModelsVisited,
					  bool HistoryFlag);

void print_and_save_main_results_per_sweep(ostringstream &ss_out,
                                          ostringstream &ss_out_n_vars_in,
                                          ostringstream &ss_out_n_models_visited,
                                          ostringstream &ss_out_log_cond_post,
                                          vector < vector <unsigned int> > vect_gam,
                                          vector < vector <unsigned int> > &List_models,
                                          vector < unsigned int > chain_idx,
                                          Double_Matrices mat_log_marg,
                                          Double_Matrices mat_log_cond_post,
                                          unsigned int sweep,
                                          unsigned int burn_in,
                                          unsigned int nModelsVisited,
                                          bool HistoryFlag);


void saveResumeFile(fstream &fResume, FILE *fRNG, unsigned int sweep,
                    double g,
                    //unsigned int *shuffleYIndex,
                    vector<unsigned int> shuffleYIndex,
                    unsigned int nY,double cumG, unsigned int countG,
                    gsl_vector *vectRMSE,
                    Temperatures *tTun,double gLs,DR *currDR,
                    vector<vector<unsigned int> > vectGam,
                    vector < unsigned int > chainIndex,
                    unsigned int nConfounders,
                    unsigned int pY,
                    gsl_rng *RandomNumberGenerator);

#endif
