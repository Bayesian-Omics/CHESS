/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from Move_monitor.cc in the ESS++ program
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

#include "Move_monitor.h"
#define DEBUG 0

using namespace std;

Move_monitor::Move_monitor(unsigned int n_chains,
			   unsigned int n_possible_CM_moves)
{
  FSMH_nb_model_tot=0;
  FSMH_nb_sweep=0;
  FSMH_nb_model_0_1=0;
  FSMH_nb_model_1_0=0;

  FSMH_nb_accept_tot=0;
  FSMH_nb_accept_0_1=0;
  FSMH_nb_accept_1_0=0;

  All_exchange_move_history.resize(3);
  All_exchange_freq.resize(n_chains);
  All_exchange_nb_sweep=0;
  All_exchange_n_accept=0;

  DR_move_history.resize(3);
  DR_freq.resize(n_chains);
  DR_nb_sweep=0;
  DR_n_accept=0;

  g_sample_nb_sweep=0;
  g_sample_nb_accept=0;
  g_adapt_history.resize(4);

  Gibbs_nb_sweep=0;
  Gibbs_nb_model=0;
  Gibbs_nb_0_1=0;
  Gibbs_nb_1_0=0;

  Gibbs_move_history.resize(4);

  temperature_history.resize(n_chains+1);

  CM_nb_sweep=0;
  CM_nb_model=0;
  CM_nb_accept=0;
  CM_n_nb_sweep_per_CM_move.resize(n_possible_CM_moves);
  CM_n_nb_accept_per_CM_move.resize(n_possible_CM_moves);
  for(unsigned int col=0;col<CM_n_nb_sweep_per_CM_move.size();col++){
    CM_n_nb_sweep_per_CM_move[col]=0;
    CM_n_nb_accept_per_CM_move[col]=0;
  }

  CM_history.resize(5);

}

void Move_monitor::display_move_monitor_full()
{
  cout << endl << "****************************************************************************************************" << endl
       << "****************************************** Moves monitor *******************************************" << endl
       << endl << "####################" << endl
       << "  -Fast Scan M-H" << endl
       << "####################" << endl
       << "\tnb_sweep = " << FSMH_nb_sweep << endl
       << "\tnb_model_tot = " << FSMH_nb_model_tot << endl
       << "\tnb_accept = " << FSMH_nb_accept_tot << endl
       << "\tnb_model_0_1 = " << FSMH_nb_model_0_1 << endl
       << "\tnb_accept_0_1 = " << FSMH_nb_accept_0_1 << endl
       << "\tnb_model_1_0 = " << FSMH_nb_model_1_0 << endl
       << "\tnb_accept_1_0 = " << FSMH_nb_accept_1_0 << endl;
  if(FSMH_list_sweep.size()>0){
    cout << "\tFSMH Per-sweep, monitoring:" << endl
	 << "\tSweep\tn_mod\tn_accept\tn_mod_0_1\tn_accept_0_1\tn_mod_1_0\tn_accept_1_0" << endl;
    for(unsigned int col=0;col<FSMH_list_sweep.size();col++){
      cout << "\t" << FSMH_list_sweep[col] << "\t"
	   << FSMH_nb_model_per_sweep[col] <<"\t"
	   << FSMH_nb_accept_per_sweep[col] <<"\t\t"
	   << FSMH_nb_model_per_sweep_0_1[col] <<"\t\t"
	   << FSMH_nb_accept_per_sweep_0_1[col] <<"\t\t"
	   << FSMH_nb_model_per_sweep_1_0[col] <<"\t\t"
	   << FSMH_nb_accept_per_sweep_1_0[col] << endl;
    }
  }

  cout << endl << "####################" << endl
       << "  -CM Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << CM_nb_sweep << endl
       << "\tnb_model = " << CM_nb_model << endl
       << "\tnb_accept = " << CM_nb_accept << endl
       << "\tCM_nb_sweep_per_CM_move" << endl
       << "\t";
  for(unsigned int col=0;col<CM_n_nb_sweep_per_CM_move.size();col++){
    cout << CM_n_nb_sweep_per_CM_move[col] << " ";
  }
  cout << endl;
  cout << "\tCM_nb_accept_per_CM_move" << endl
       << "\t";
  for(unsigned int col=0;col<CM_n_nb_accept_per_CM_move.size();col++){
    cout << CM_n_nb_accept_per_CM_move[col] << " ";
  }
  cout << endl;
  if(CM_history[0].size()>0){
    cout << "\tCM per-sweep monitoring:" << endl;
    cout << "\tSweep\ttype\t#Bkpoints\tchain1\tchain2" << endl;
    for(unsigned int col=0;col<CM_history[0].size();col++){
      cout << "\t" << CM_history[0][col]
	   << "\t" << CM_history[1][col]
	   << "\t" << CM_history[2][col]
	   << "\t" << CM_history[3][col] 
	   << "\t" << CM_history[4][col] 
	   << endl;
    }
  }



  cout << endl << "####################" << endl
       << "  -All Exchange" << endl
       << "####################" << endl
       << "\tnb_sweep = " << All_exchange_nb_sweep << endl
       << "\tnb_accept = " << All_exchange_n_accept << endl
       << "\tVect Freq:" << endl;
  cout << "\t";
  for(unsigned int col=0;col<All_exchange_freq.size();col++){
    cout << All_exchange_freq[col] << "\t";
  }
  cout << endl;
  if(All_exchange_move_history[0].size()>0){
    cout << "\tAll exchange per-sweep monitoring:" << endl;
    cout << "\tSweep\tChain_i\tChain_f" << endl;
    for(unsigned int col=0;col<All_exchange_move_history[0].size();col++){
      cout << "\t" << All_exchange_move_history[0][col]
	   << "\t" << All_exchange_move_history[1][col]
	   << "\t" << All_exchange_move_history[2][col] << endl;
    }
  }
  cout << endl << "####################" << endl
       << "  -DR Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << DR_nb_sweep << endl
       << "\tnb_accept = " << DR_n_accept << endl
       << "\tVect Freq:" << endl;
  cout << "\t";
  for(unsigned int col=0;col<DR_freq.size();col++){
    cout << DR_freq[col] << "\t";
  }
  cout << endl;
  if(DR_move_history[0].size()>0){
    cout << "\tDR per-sweep monitoring:" << endl;
    cout << "\tSweep\tChain_i\tChain_f" << endl;
    for(unsigned int col=0;col<DR_move_history[0].size();col++){
      cout << "\t" << DR_move_history[0][col]
	   << "\t" << DR_move_history[1][col]
	   << "\t" << DR_move_history[2][col] << endl;
    }
  }
  cout << endl << "####################" << endl
       << "  -g Sample" << endl
       << "####################" << endl
       << "\tnb_sweep = " << g_sample_nb_sweep << endl
       << "\tnb_accept = " << g_sample_nb_accept << endl;
  if(g_sample_history.size()>0){
    cout << "\tHistory:" << endl;
    cout << "\t";
    for(unsigned int col=0;col<g_sample_history.size();col++){
      cout << g_sample_history[col] << "\t";
    }
    cout << endl;
  }
  if(g_adapt_history[0].size()>0){
    cout << endl
	 << "\tg_adapt History:" << endl;
    cout << "\tSweep\t%accept\tls\tls_new" << endl;
    for(unsigned int row=0;row<g_adapt_history[0].size();row++){
      cout << "\t" << g_adapt_history[0][row]
	   << "\t" << g_adapt_history[1][row]
	   << "\t" << g_adapt_history[2][row]
	   << "\t" << g_adapt_history[3][row] << endl;
      
    }
  }

  cout << endl << "####################" << endl
       << "  -Gibbs Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << Gibbs_nb_sweep << endl
       << "\tnb_model = " << Gibbs_nb_model << endl
       << "\tnb_0_1 moves = " << Gibbs_nb_0_1 << endl
       << "\tnb_1_0 moves = " << Gibbs_nb_1_0 << endl;
  if(Gibbs_move_history[0].size()>0){
    cout << "\tGibbs per-sweep monitoring:" << endl;
    cout << "\tSweep\tn_0_1\tn_1_0\tn_unchanged" << endl;
    for(unsigned int col=0;col<Gibbs_move_history[0].size();col++){
      cout << "\t" << Gibbs_move_history[0][col]
	   << "\t" << Gibbs_move_history[1][col]
	   << "\t" << Gibbs_move_history[2][col]
	   << "\t" << Gibbs_move_history[3][col] << endl;
    }
  }
  



  cout << endl << "####################" << endl
       << "  -Temperatures" << endl
       << "####################" << endl;
  
  cout << "\tSweep ";
  for(unsigned int chain=1;chain<temperature_history.size();chain++){
    cout << "\tchain # " << chain;
  }
  cout << endl;
  for(unsigned int sweep=0;sweep<temperature_history[0].size();sweep++){
    for(unsigned int row=0;row<temperature_history.size();row++){
      cout <<"\t " << temperature_history[row][sweep];
    }
    cout << endl;
  }
  
  cout << endl << "####################" << endl << endl;
 
  cout << "****************************************************************************************************" << endl
       << "****************************************************************************************************" << endl << endl;
  
}

void Move_monitor::display_move_monitor()
{
  cout << endl << "****************************************************************************************************" << endl
       << "****************************************** Moves monitor *******************************************" << endl
       << endl << "####################" << endl
       << "  -Fast Scan M-H" << endl
       << "####################" << endl
       << "\tnb_sweep = " << FSMH_nb_sweep << endl
       << "\tnb_model_tot = " << FSMH_nb_model_tot << endl
       << "\tnb_accept = " << FSMH_nb_accept_tot << endl
       << "\tnb_model_0_1 = " << FSMH_nb_model_0_1 << endl
       << "\tnb_accept_0_1 = " << FSMH_nb_accept_0_1 << endl
       << "\tnb_model_1_0 = " << FSMH_nb_model_1_0 << endl
       << "\tnb_accept_1_0 = " << FSMH_nb_accept_1_0 << endl;

  cout << endl << "####################" << endl
       << "  -CM Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << CM_nb_sweep << endl
       << "\tnb_model = " << CM_nb_model << endl
       << "\tnb_accept = " << CM_nb_accept << endl
       << "\tCM_nb_sweep_per_CM_move" << endl
       << "\t";
  for(unsigned int col=0;col<CM_n_nb_sweep_per_CM_move.size();col++){
    cout << CM_n_nb_sweep_per_CM_move[col] << " ";
  }
  cout << endl;
  cout << "\tCM_nb_accept_per_CM_move" << endl
       << "\t";
  for(unsigned int col=0;col<CM_n_nb_accept_per_CM_move.size();col++){
    cout << CM_n_nb_accept_per_CM_move[col] << " ";
  }
  cout << endl;

  cout << endl << "####################" << endl
       << "  -All Exchange" << endl
       << "####################" << endl
       << "\tnb_sweep = " << All_exchange_nb_sweep << endl
       << "\tnb_accept = " << All_exchange_n_accept << endl
       << "\tVect Freq:" << endl;
  cout << "\t";
  for(unsigned int col=0;col<All_exchange_freq.size();col++){
    cout << All_exchange_freq[col] << "\t";
  }
  cout << endl;
  cout << endl << "####################" << endl
       << "  -DR Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << DR_nb_sweep << endl
       << "\tnb_accept = " << DR_n_accept << endl
       << "\tVect Freq:" << endl;
  cout << "\t";
  for(unsigned int col=0;col<DR_freq.size();col++){
    cout << DR_freq[col] << "\t";
  }
  cout << endl;
  cout << endl << "####################" << endl
       << "  -g Sample" << endl
       << "####################" << endl
       << "\tnb_sweep = " << g_sample_nb_sweep << endl
       << "\tnb_accept = " << g_sample_nb_accept << endl;

  cout << endl << "####################" << endl
       << "  -Gibbs Move" << endl
       << "####################" << endl
       << "\tnb_sweep = " << Gibbs_nb_sweep << endl
       << "\tnb_model = " << Gibbs_nb_model << endl
       << "\tnb_0_1 moves = " << Gibbs_nb_0_1 << endl
       << "\tnb_1_0 moves = " << Gibbs_nb_1_0 << endl;
 
  cout << "****************************************************************************************************" << endl
       << "****************************************************************************************************" << endl << endl;
  
}

void Move_monitor::print_move_monitor_full(ofstream &f_out_FSMH,
					   ofstream &f_out_CM,
					   ofstream &f_out_AE,
					   ofstream &f_out_DR,
					   ofstream &f_out_g,
					   ofstream &f_out_g_adapt,
					   ofstream &f_out_Gibbs,
					   ofstream &f_out_t_tun,
					   unsigned int g_sample,
					   bool iso_T_Flag)
{
  f_out_FSMH << "Sweep\tn_mod\tn_accept\tn_mod_0_1\tn_accept_0_1\tn_mod_1_0\tn_accept_1_0" << endl;
  if(FSMH_list_sweep.size()>0){
    for(unsigned int col=0;col<FSMH_list_sweep.size();col++){
      f_out_FSMH << FSMH_list_sweep[col] << "\t"
		 << FSMH_nb_model_per_sweep[col] <<"\t"
		 << FSMH_nb_accept_per_sweep[col] <<"\t"
		 << FSMH_nb_model_per_sweep_0_1[col] <<"\t"
		 << FSMH_nb_accept_per_sweep_0_1[col] <<"\t"
		 << FSMH_nb_model_per_sweep_1_0[col] <<"\t"
		 << FSMH_nb_accept_per_sweep_1_0[col] << endl;
    }
  }
  
  f_out_CM << "Sweep\tMove_type\t#Breakpoints\tChain_l\tChain_r" << endl;
  if(CM_history[0].size()>0){
    for(unsigned int col=0;col<CM_history[0].size();col++){
      f_out_CM << CM_history[0][col]
	       << "\t" << CM_history[1][col]+1
	       << "\t" << CM_history[2][col]
	       << "\t" << CM_history[3][col]+1 
	       << "\t" << CM_history[4][col]+1 
	       << endl;
    }
  }
  
  f_out_AE << "Sweep\tChain_l\tChain_r" << endl;
  if(All_exchange_move_history[0].size()>0){
    for(unsigned int col=0;col<All_exchange_move_history[0].size();col++){
      f_out_AE << All_exchange_move_history[0][col]
	       << "\t" << All_exchange_move_history[1][col]+1
	       << "\t" << All_exchange_move_history[2][col]+1 << endl;
    }
  }

  f_out_DR << "Sweep\tChain_l\tChain_r" << endl;
  
  if(DR_move_history[0].size()>0){
    for(unsigned int col=0;col<DR_move_history[0].size();col++){
      f_out_DR << DR_move_history[0][col]
	       << "\t" << DR_move_history[1][col]+1
	       << "\t" << DR_move_history[2][col]+1 << endl;
    }
  }
  if(g_sample!=0){
    f_out_g << "Sweep\tg" << endl;
    if(g_sample_history.size()>0){
      for(unsigned int col=0;col<g_sample_history.size();col++){
	f_out_g << col+1 << "\t" 
		<< g_sample_history[col] << endl;
      }
      f_out_g << endl;
    }
    f_out_g_adapt << "Sweep\tAcceptance_rate\tlog_proposal_std" << endl;
    if(g_adapt_history[0].size()>0){
      for(unsigned int row=0;row<g_adapt_history[0].size();row++){
	f_out_g_adapt << g_adapt_history[0][row]
		      << "\t" << g_adapt_history[1][row]
		      << "\t" << g_adapt_history[3][row] << endl;
      }
    }
  }
  f_out_Gibbs << "Sweep\tn_0->1\tn_1->0" << endl;
  if(Gibbs_move_history[0].size()>0){
    for(unsigned int col=0;col<Gibbs_move_history[0].size();col++){
      f_out_Gibbs << Gibbs_move_history[0][col]
		  << "\t" << Gibbs_move_history[1][col]
		  << "\t" << Gibbs_move_history[2][col] << endl;
    }
  }
  if(!iso_T_Flag){
    f_out_t_tun << "Sweep\t";
    for(unsigned int chain=1;chain<temperature_history.size();chain++){
      f_out_t_tun << "Chain_" << chain << "\t";
    }
    f_out_t_tun << endl;
    if(temperature_history[0].size()>0){
      for(unsigned int sweep=0;sweep<temperature_history[0].size();sweep++){
	f_out_t_tun << temperature_history[0][sweep] << "\t";
	for(unsigned int row=1;row<temperature_history.size();row++){
	  f_out_t_tun << temperature_history[row][sweep] << "\t";
	}
	f_out_t_tun << endl;
      }
    }
  }
}

void Move_monitor::print_move_monitor_per_sweep(ofstream &f_out_FSMH,
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
                                           unsigned int nChains)
{

  unsigned int lastSweep=0;
  unsigned int col=0;
  unsigned int row=0;

  if(sweep==1){
    f_out_FSMH << "Sweep\tn_mod\tn_accept\tn_mod_0_1\tn_accept_0_1\tn_mod_1_0\tn_accept_1_0" << endl;
  }
  lastSweep=FSMH_list_sweep.size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(FSMH_list_sweep[col]==sweep){
      f_out_FSMH << FSMH_list_sweep[col] << "\t"
                 << FSMH_nb_model_per_sweep[col] <<"\t"
                 << FSMH_nb_accept_per_sweep[col] <<"\t"
                 << FSMH_nb_model_per_sweep_0_1[col] <<"\t"
                 << FSMH_nb_accept_per_sweep_0_1[col] <<"\t"
                 << FSMH_nb_model_per_sweep_1_0[col] <<"\t"
                 << FSMH_nb_accept_per_sweep_1_0[col] << endl;

    }
  }

  if(sweep==1){
    f_out_CM << "Sweep\tMove_type\tnBreakpoints\tChain_l\tChain_r" << endl;
  }
  lastSweep=CM_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(CM_history[0][col]==sweep){
      f_out_CM << CM_history[0][col]
               << "\t" << CM_history[1][col]+1
               << "\t" << CM_history[2][col]
               << "\t" << CM_history[3][col]+1
               << "\t" << CM_history[4][col]+1
               << endl;

    }
  }

  if(sweep==1){
    f_out_AE << "Sweep\tChain_l\tChain_r" << endl;
  }
  lastSweep=All_exchange_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(All_exchange_move_history[0][col]==sweep){
      f_out_AE << All_exchange_move_history[0][col]
               << "\t" << All_exchange_move_history[1][col]+1
               << "\t" << All_exchange_move_history[2][col]+1 << endl;

    }
  }

  if(sweep==1){
    f_out_DR << "Sweep\tChain_l\tChain_r" << endl;
  }
  lastSweep=DR_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(DR_move_history[0][col]==sweep){
      f_out_DR << DR_move_history[0][col]
               << "\t" << DR_move_history[1][col]+1
               << "\t" << DR_move_history[2][col]+1 << endl;

    }
  }

  if(g_sample!=0){
    if(sweep==1){
      f_out_g << "Sweep\tg" << endl;
    }
    lastSweep=g_sample_history.size();
    col=lastSweep-1;
    f_out_g << sweep << "\t"
            << g_sample_history[col] << endl;


    if(sweep==1){
      f_out_g_adapt << "Sweep\tAcceptance_rate\tlog_proposal_std" << endl;
    }
    lastSweep=g_adapt_history[0].size();
    row=lastSweep-1;
    if(lastSweep>0){
      if(g_adapt_history[0][row]==sweep){
        f_out_g_adapt << g_adapt_history[0][row]
                      << "\t" << g_adapt_history[1][row]
                      << "\t" << g_adapt_history[3][row] << endl;

      }
    }
  }

  if(sweep==1){
    f_out_Gibbs << "Sweep\tn_0->1\tn_1->0" << endl;
  }
  lastSweep=Gibbs_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(Gibbs_move_history[0][col]==sweep){
      f_out_Gibbs << Gibbs_move_history[0][col]
                  << "\t" << Gibbs_move_history[1][col]
                  << "\t" << Gibbs_move_history[2][col] << endl;

    }
  }
  if(!iso_T_Flag){
    if(sweep==1){
      f_out_t_tun << "Sweep\t";
      for(unsigned int chain=1;chain<nChains;chain++){
        f_out_t_tun << "Chain_" << chain << "\t";
      }
      f_out_t_tun << endl;
    }
    lastSweep=temperature_history[0].size();
    col=lastSweep-1;
    if(lastSweep>0){
      if(temperature_history[0][col]==sweep){
        f_out_t_tun << temperature_history[0][col] << "\t";
        for(unsigned int row=1;row<=nChains;row++){
          f_out_t_tun << temperature_history[row][col] << "\t";
        }
        f_out_t_tun << endl;
      }
    }
  }
}

void Move_monitor::print_move_monitor_per_sweep(ostringstream &ss_out_FSMH,
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
                                           unsigned int nChains)
{

  unsigned int lastSweep=0;
  unsigned int col=0;
  unsigned int row=0;

  if(sweep==1){
    ss_out_FSMH << "Sweep\tn_mod\tn_accept\tn_mod_0_1\tn_accept_0_1\tn_mod_1_0\tn_accept_1_0" << endl;
  }
  lastSweep=FSMH_list_sweep.size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(FSMH_list_sweep[col]==sweep){
      ss_out_FSMH << FSMH_list_sweep[col] << "\t"
                 << FSMH_nb_model_per_sweep[col] <<"\t"
                 << FSMH_nb_accept_per_sweep[col] <<"\t"
                 << FSMH_nb_model_per_sweep_0_1[col] <<"\t"
                 << FSMH_nb_accept_per_sweep_0_1[col] <<"\t"
                 << FSMH_nb_model_per_sweep_1_0[col] <<"\t"
                 << FSMH_nb_accept_per_sweep_1_0[col] << endl;

    }
  }

  if(sweep==1){
    ss_out_CM << "Sweep\tMove_type\tnBreakpoints\tChain_l\tChain_r" << endl;
  }
  lastSweep=CM_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(CM_history[0][col]==sweep){
      ss_out_CM << CM_history[0][col]
               << "\t" << CM_history[1][col]+1
               << "\t" << CM_history[2][col]
               << "\t" << CM_history[3][col]+1
               << "\t" << CM_history[4][col]+1
               << endl;

    }
  }

  if(sweep==1){
    ss_out_AE << "Sweep\tChain_l\tChain_r" << endl;
  }
  lastSweep=All_exchange_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(All_exchange_move_history[0][col]==sweep){
      ss_out_AE << All_exchange_move_history[0][col]
               << "\t" << All_exchange_move_history[1][col]+1
               << "\t" << All_exchange_move_history[2][col]+1 << endl;

    }
  }

  if(sweep==1){
    ss_out_DR << "Sweep\tChain_l\tChain_r" << endl;
  }
  lastSweep=DR_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(DR_move_history[0][col]==sweep){
      ss_out_DR << DR_move_history[0][col]
               << "\t" << DR_move_history[1][col]+1
               << "\t" << DR_move_history[2][col]+1 << endl;

    }
  }

  if(g_sample!=0){
    if(sweep==1){
      ss_out_g << "Sweep\tg" << endl;
    }
    lastSweep=g_sample_history.size();
    col=lastSweep-1;
    ss_out_g << sweep << "\t"
            << g_sample_history[col] << endl;


    if(sweep==1){
      ss_out_g_adapt << "Sweep\tAcceptance_rate\tlog_proposal_std" << endl;
    }
    lastSweep=g_adapt_history[0].size();
    row=lastSweep-1;
    if(lastSweep>0){
      if(g_adapt_history[0][row]==sweep){
        ss_out_g_adapt << g_adapt_history[0][row]
                      << "\t" << g_adapt_history[1][row]
                      << "\t" << g_adapt_history[3][row] << endl;

      }
    }
  }

  if(sweep==1){
    ss_out_Gibbs << "Sweep\tn_0->1\tn_1->0" << endl;
  }
  lastSweep=Gibbs_move_history[0].size();
  col=lastSweep-1;
  if(lastSweep>0){
    if(Gibbs_move_history[0][col]==sweep){
      ss_out_Gibbs << Gibbs_move_history[0][col]
                  << "\t" << Gibbs_move_history[1][col]
                  << "\t" << Gibbs_move_history[2][col] << endl;

    }
  }
  if(!iso_T_Flag){
    if(sweep==1){
      ss_out_t_tun << "Sweep\t";
      for(unsigned int chain=1;chain<nChains;chain++){
        ss_out_t_tun << "Chain_" << chain << "\t";
      }
      ss_out_t_tun << endl;
    }
    lastSweep=temperature_history[0].size();
    col=lastSweep-1;
    if(lastSweep>0){
      if(temperature_history[0][col]==sweep){
        ss_out_t_tun << temperature_history[0][col] << "\t";
        for(unsigned int row=1;row<=nChains;row++){
          ss_out_t_tun << temperature_history[row][col] << "\t";
        }
        ss_out_t_tun << endl;
      }
    }
  }
}

