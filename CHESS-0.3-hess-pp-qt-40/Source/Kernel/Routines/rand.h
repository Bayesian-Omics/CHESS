/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from rand.h in the ESS++ program
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

#ifndef RAND_H_
#define RAND_H_


#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdint.h> 

using namespace std;

double myrand(gsl_rng *RandomNumberGenerator);

void smyrand(long seed, gsl_rng *RandomNumberGenerator);

void readRNG(FILE *fRNG,gsl_rng *RandomNumberGenerator);

void writeRNG(FILE *fRNG,gsl_rng *RandomNumberGenerator);

double gennor(double av,double sd, gsl_rng *RandomNumberGenerator);

int genBernoulli(double p, gsl_rng *RandomNumberGenerator);

void My_Permut_unsigned_int(gsl_permutation *MyPerm,gsl_rng *RandomNumberGenerator);

int SampleFromDiscrete_All_exchange(gsl_matrix *description_all_exchange,gsl_rng *RandomNumberGenerator);

unsigned int SampleFromDiscrete_non_cum(vector<double> &pbty,gsl_rng *RandomNumberGenerator);

unsigned int SampleFromDiscrete_new( vector<double> &cdf,gsl_rng *RandomNumberGenerator);

void My_gsl_ran_choose_double(void * dest, size_t k, void * source, size_t n,gsl_rng *RandomNumberGenerator);

void My_gsl_ran_choose_u_int(void * dest, size_t k, void * source, size_t n,gsl_rng *RandomNumberGenerator);

#endif /*RAND_H_*/
