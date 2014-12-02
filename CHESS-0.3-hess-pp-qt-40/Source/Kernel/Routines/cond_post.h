/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from cond_post.h in the ESS++ program
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


#ifndef COND_POST_H
#define COND_POST_H

#include "../../General_Classes/Model_Information.h"

// DEBUG
//#define _CUDA_ 0


#include <gsl/gsl_fit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include "../Classes/Double_Matrices.h"
#include "../Classes/Double_Matrices_cont.h"
#include "../Classes/Int_Matrices.h"
#include "../Routines/matrix_handling.h"
#include "../Routines/rand.h"
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
#include <math.h>
#include <ctype.h>


#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_gamma.h>
#include "../Classes/Prior_param.h"

void get_vect_gam_init(vector < vector <unsigned int> > &vect_gam,
		       vector < vector <unsigned int> > &Gam_step_regr,
		       bool iso_T_Flag,
		       unsigned int maxPGamma,
		       gsl_rng *RandomNumberGenerator);

void getEigenDecomposition(gsl_matrix *matXGam,
                           gsl_matrix* gslEigenVecs,
                           gsl_vector* gslEigenVals,
                           unsigned int pXGam);

void getTildeResidualsMatrix(gsl_matrix * TildeResidualsMatrix,
                            gsl_matrix *QSubTYTilde);

void getQSubTYTilde(gsl_matrix * QSubTYTilde,
                            gsl_matrix *matXGamTilde,
                            gsl_matrix *matYTilde,
                            gsl_matrix *matY);

void getQSubTYTildeCula(gsl_matrix * QSubTYTilde,
                        float *matrixXGamTilde,
                        float *matrixYTilde,
                        gsl_matrix *matY,
                        unsigned int nQ,
                        unsigned int pXGam);


gsl_matrix *getSGamma(gsl_matrix *mat_X,
                   gsl_matrix *mat_Y,
                   double lambda,
                   double g,
                   gsl_matrix *gslEigenVecs,
                   gsl_vector *gslEigenVals,
                   bool gPriorFlag,
                   bool indepPriorFlag,
                   unsigned int pXGam,
                   unsigned int nX,
                   unsigned int pY
                   //,gsl_matrix *R2_Mat_GPriors
                   );

void Quick_getSGamma_GPriors(gsl_matrix *matSGamma,
                             gsl_matrix * R2_Mat_GPriors,
                             gsl_matrix *matYTY,
                             unsigned int nX,
                             unsigned int pXGam,
                             unsigned int pY,
                             double g);

void Get_R2_Mat_GPriors(gsl_matrix * R2_Mat_GPriors,
                        vector<unsigned int> &Gamma,
                        gsl_matrix *mat_X,
                        gsl_matrix * matY,
                        bool cudaFlag);


#if _CUDA_
void getEigenDecompositionCula(gsl_matrix *matXGam,
                               float* matrixEigenVecs,
                               float* vectorEigenVals,
                               unsigned int pXGam);

gsl_matrix *getSGammaCula(gsl_matrix *matXGam,
                             gsl_matrix *matY,
                             double lambda,
                             double g,
                             float* eigenVecs,
                             float* eigenVals,
                             bool gPriorFlag,
                             bool indepPriorFlag,
                             unsigned int pXGam,
                             unsigned int nX,
                             unsigned int pY
                             //,gsl_matrix * R2_Mat_GPriors
                             );


#endif

double getPriorGam_ESS(Prior_param PR,
                    unsigned int pX,
                    unsigned int pXGam);

double getPriorGam_HESS(vector <unsigned int> &Gamma,
                        unsigned int Chain,
                        unsigned int Response,
                        Model_Information &Model_Data);

double getPriorG(Prior_param PR,
			  bool gSampleFlag,
			  double g);

double invGammaPdf(double x,
			double alpha,
			double beta);

double getLogMarg(Prior_param PR,
                    gsl_matrix *matXGam,
                    gsl_matrix *matSGamma,
                    float* matrixEigenVecs,
                    float* vectorEigenVals,
                    double lambda,
                    double g,
                    bool gPriorFlag,
                    bool indepPriorFlag,
                    unsigned int pXGam,
                    unsigned int nX,
                    unsigned int pY);

void computeLogPosterior(vector <unsigned int> &Gamma,
                         unsigned int Chain,
                         unsigned int Response,
                         gsl_matrix *mat_X,
                         double& logMargLik,
                         double& logPosterior,
                         gsl_matrix *matY,
                         Prior_param PR,
                         bool gPriorFlag,
                         bool indepPriorFlag,
                         bool gSampleFlag,
                         double lambda,
                         double g,
                         bool cudaFlag,
                         Model_Information &Model_Data);

void computeLogPosterior_Linear(vector <unsigned int> &Gamma,
                                unsigned int Chain,
                                unsigned int Response,
                                gsl_matrix *mat_X,
                                double& logMargLik,
                                double& logPosterior,
                                gsl_matrix *matY,
                                Prior_param PR,
                                bool gPriorFlag,
                                bool indepPriorFlag,
                                bool gSampleFlag,
                                double lambda,
                                double g,
                                bool cudaFlag,
                                Model_Information &Model_Data);

double LinearLogMarg(vector <unsigned int> &Gamma,
                 gsl_matrix *mat_X,
                 gsl_matrix *matY,
                 Prior_param PR,
                 bool gPriorFlag,
                 bool indepPriorFlag,
                 double lambda,
                 double g,
                 bool cudaFlag
                 //,gsl_matrix *R2_Mat_GPriors
                 );

void computeLogPosterior_Logistic(vector <unsigned int> &Gamma,
                                  unsigned int Chain,
                                  unsigned int Response,
                                  gsl_matrix *mat_X,
                                  double& logMargLik,
                                  double& logPosterior,
                                  gsl_matrix *matY,
                                  Prior_param PR,
                                  bool gPriorFlag,
                                  bool indepPriorFlag,
                                  bool gSampleFlag,
                                  double lambda,
                                  double g,
                                  bool cudaFlag,
                                  Model_Information &Model_Data);

void computeLogPosterior_PostProc();


#endif /*  */
