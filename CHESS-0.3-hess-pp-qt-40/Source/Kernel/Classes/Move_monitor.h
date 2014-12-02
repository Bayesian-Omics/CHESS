/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Move_monitor.h in the ESS++ program
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

#ifndef MOVE_MONITOR_H
#define MOVE_MONITOR_H 1

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <math.h>
#include <vector>
using namespace std;

class Move_monitor
{
 public:
  Move_monitor(unsigned int n_chains,
	       unsigned int n_possible_CM_moves);
  ~Move_monitor(){};

  //Monitoring FSMH Moves

  unsigned int FSMH_nb_model_tot;
  unsigned int FSMH_nb_sweep;
  unsigned int FSMH_nb_model_0_1;
  unsigned int FSMH_nb_model_1_0;

  unsigned int FSMH_nb_accept_tot;
  unsigned int FSMH_nb_accept_0_1;
  unsigned int FSMH_nb_accept_1_0;

  vector < unsigned int > FSMH_list_sweep;
  vector < unsigned int > FSMH_nb_model_per_sweep;
  vector < unsigned int > FSMH_nb_model_per_sweep_0_1;
  vector < unsigned int > FSMH_nb_model_per_sweep_1_0;
  vector < unsigned int > FSMH_nb_accept_per_sweep;
  vector < unsigned int > FSMH_nb_accept_per_sweep_0_1;
  vector < unsigned int > FSMH_nb_accept_per_sweep_1_0;

  //Monitoring All Exchange Moves

  unsigned int All_exchange_nb_sweep;
  vector < unsigned int > All_exchange_freq;
  unsigned int All_exchange_n_accept;
  vector < vector < unsigned int > > All_exchange_move_history;

  //Monitoring DR Moves
  unsigned int DR_nb_sweep;
  vector < unsigned int > DR_freq;
  unsigned int DR_n_accept;
  vector < vector < unsigned int > > DR_move_history;

  //Monitoring the sampling of g
  unsigned int g_sample_nb_sweep;
  unsigned int g_sample_nb_accept;
  vector < double > g_sample_history;
  vector < vector < double > > g_adapt_history;

  //Monitoring Gibbs Moves

  unsigned int Gibbs_nb_sweep;
  unsigned int Gibbs_nb_model;
  unsigned int Gibbs_nb_0_1;
  unsigned int Gibbs_nb_1_0;

  vector < vector < unsigned int > >Gibbs_move_history;

  //Monitoring Temperature
  vector < vector < double > > temperature_history;

  //Monitoring Crossover Moves
  unsigned int CM_nb_sweep;
  unsigned int CM_nb_model;
  unsigned int CM_nb_accept;
  vector < unsigned int > CM_n_nb_sweep_per_CM_move;
  vector < unsigned int > CM_n_nb_accept_per_CM_move;

  vector < vector < unsigned int > > CM_history;


  void display_move_monitor_full();
  void display_move_monitor();
  void print_move_monitor_full(ofstream &f_out_FSMH,
			       ofstream &f_out_CM,
			       ofstream &f_out_AE,
			       ofstream &f_out_DR,
			       ofstream &f_out_g,
			       ofstream &f_out_g_adapt,
			       ofstream &f_out_Gibbs,
			       ofstream &f_out_t_tun,
			       unsigned int g_sample,
			       bool iso_T_Flag);

  void print_move_monitor_per_sweep(ofstream &f_out_FSMH,
                                 ofstream &f_out_CM,
                                 ofstream &f_out_AE,
                                 ofstream &f_out_DR,
                                 ofstream &f_out_g,
                                 ofstream &f_out_g_adapt,
                                 ofstream &f_out_Gibbs,
                                 ofstream &f_out_t_tun,
                                 unsigned int g_sample,
                                 bool iso_T_Flag,
                                 unsigned int sweep,
                                 unsigned int nChains);

  void print_move_monitor_per_sweep(ostringstream &ss_out_FSMH,
                                 ostringstream &ss_out_CM,
                                 ostringstream &ss_out_AE,
                                 ostringstream &ss_out_DR,
                                 ostringstream &ss_out_g,
                                 ostringstream &ss_out_g_adapt,
                                 ostringstream &ss_out_Gibbs,
                                 ostringstream &ss_out_t_tun,
                                 unsigned int g_sample,
                                 bool iso_T_Flag,
                                 unsigned int sweep,
                                 unsigned int nChains);
  
};

#endif /* !defined MOVE_MONITOR_H */
