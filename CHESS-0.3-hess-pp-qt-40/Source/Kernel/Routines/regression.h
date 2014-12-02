/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from regression.h in the ESS++ program
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

#ifndef REGRESSION_H
#define REGRESSION_H

#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include "../Classes/Double_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Int_Matrices.h"
#include "../Routines/matrix_handling.h"
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <cstring>
#include <iostream>
#include <sstream> 
#include <cstdarg>
#include <math.h>
#include <cmath>
#include <gsl/gsl_statistics_double.h>

void perform_stepwise_regr(vector < vector <unsigned int> > &Gam_stepregr,
			   Double_Matrices mat_Y,
			   Double_Matrices mat_X,
			   bool scaling);

unsigned int get_rank_from_R(gsl_matrix *M,
			     double tolerance);

gsl_vector* get_vect_elem_square(gsl_matrix *M);

void get_vect_beta_var_out(gsl_vector *vect_beta_var_out,
			   gsl_vector *xx,
			   gsl_matrix *xr,
			   gsl_vector *residuals);

gsl_vector *get_se(gsl_matrix *mat_R_red,
		   double RMSE);

void get_mat_Q_red_R_red_beta_in_resid(gsl_matrix *mat_X_gam,
				       gsl_vector *current_outcome,
				       gsl_vector *vect_residuals,
				       gsl_vector *vect_beta_var_in,
				       gsl_matrix *mat_Q_red,
				       gsl_matrix *mat_R_red);

void get_X_residuals(gsl_matrix *mat_X_gam_bar,
		     gsl_matrix *mat_Q_red);

void get_se_var_out(gsl_vector *S2,
		    gsl_matrix *xr,
		    gsl_vector *vect_residuals,
		    gsl_vector *vect_beta,
		    gsl_vector *xx,
		    double df2);

void compute_beta_SE_vars_out(gsl_vector *vect_beta_var_out,
			      gsl_vector *vect_SE_var_out,
			      gsl_matrix *mat_X_gam_bar,
			      gsl_matrix *mat_Q_red,
			      gsl_vector *vect_residuals,
			      double df2);

void get_full_pvalues_and_beta(gsl_vector *vect_p_value,
			       gsl_vector *vect_beta_full,
			       gsl_vector *vect_SE_full,
			       gsl_matrix *mat_X_gam,
			       gsl_matrix *mat_X_gam_bar,
			       gsl_vector *current_outcome,
			       gsl_vector *vect_residuals,
			       gsl_vector *vect_RMSE,
			       double tolerance,
			       vector < unsigned int > &list_columns_X_gam,
			       vector < unsigned int > &list_columns_X_gam_bar,
			       unsigned int outcome);

void getEstimateRMSE(gsl_matrix *mat_X_gam,
                     gsl_vector *current_outcome,
                     unsigned int whichOutcome,
                     gsl_vector *vect_RMSE,
                     double tolerance);

int update_is_var_in(vector < unsigned int > &is_var_in,
		     vector < unsigned int > &list_columns_X_gam,
		     vector < unsigned int > &list_columns_X_gam_bar,
		     gsl_vector *vect_p_value,
		     double Pvalue_enter,
		     double Pvalue_remove,
		     unsigned int loop,
		     unsigned int n_loop_max,
		     unsigned int nConfounders);

void store_model_per_outcome(vector < vector <unsigned int> > &Gam_step_regr,
			     vector < unsigned int > &list_columns_X_gam,
			     gsl_vector *vect_p_value,
			     gsl_vector *vect_beta_full,
			     gsl_vector *vect_SE_full,
			     Double_Matrices Gam_step_regr_pvals,
			     Double_Matrices Gam_step_regr_SE,
			     Double_Matrices Gam_step_regr_beta,
			     unsigned int outcome);

#endif /*  */
