/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is copied from matrix_handling.cc in the ESS++ program
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

#include "matrix_handling.h"

#define DEBUG 0

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF


using namespace std;

void standardize_matrix_gsl(gsl_matrix* M)
{

  int nb_rows=M->size1;
  int nb_columns=M->size2;
  for(int col=0;col<nb_columns;col++){
    double mean=0.0;
    for(int row=0;row<nb_rows;row++){
      mean+=gsl_matrix_get(M,row,col);
    }
    mean/=(double)nb_rows;
    double temp_var=0.0;
    for(int row=0;row<nb_rows;row++){
      temp_var+=pow((gsl_matrix_get(M,row,col)-mean),2.0);
    }
    if(nb_rows>1){
      temp_var/=nb_rows-1;
    }
    else{
      cout << "USAGE::standardize_matrix :: matrix size <2 to compute variance!!!!! run_stopped " << endl;
      exit(1);
    }
    double Mystd=sqrt(temp_var);
    if(Mystd>0){
      for(int row=0;row<nb_rows;row++){
        gsl_matrix_set(M,row,col,gsl_matrix_get(M,row,col)/Mystd);
      }
    }
    else{
      cout << "USAGE::standardize_matrix::  estimated std=0; matrix unchanged. " << endl;
    }
  }
}


void center_matrix_gsl(gsl_matrix *M)
{
  
  for(unsigned int col=0;col<M->size2;col++){
    double mean=0.0;
    for(unsigned int row=0;row<M->size1;row++){
      mean+=gsl_matrix_get(M,row,col);
    }
    mean/=(double)M->size1;
    for(unsigned int row=0;row<M->size1;row++){
      gsl_matrix_set(M,row,col,gsl_matrix_get(M,row,col)-mean);
    }
  }
}

gsl_matrix* Double_matrices_cont_2_gsl_matrix(Double_Matrices_cont source)
{
  gsl_matrix* M=gsl_matrix_calloc(source.nb_rows,source.nb_columns);
  M->data=&source.matrix[0];
  return M;
}

gsl_matrix* Double_matrices_2_gsl_matrix(Double_Matrices source)
{
  gsl_matrix* M=gsl_matrix_calloc(source.nb_rows,source.nb_columns);
  for(int j=0;j<source.nb_rows;j++){
    for(int k=0;k<source.nb_columns;k++){
      M->data[j*source.nb_columns+k]=source.matrix[j][k];
    }
  }
  return M;
}

unsigned int sum_vector_int(vector <unsigned int> &Myvector)
{
  unsigned int Mysum=0.0;
  for(unsigned int i=0;i<Myvector.size();i++){
    Mysum+=Myvector[i];
  }
  return Mysum;
}

void get_list_var_in_and_out(vector <unsigned int> &list_columns_var_in,
			     vector <unsigned int> &list_columns_var_out,
			     vector <unsigned int> &is_var_in)
{

  for(unsigned int i=0;i<is_var_in.size();i++){
    if(is_var_in[i]==1){
      list_columns_var_in.push_back(i);
    }
    else{
      list_columns_var_out.push_back(i);
    }
  }
}

void get_list_var_in(vector <unsigned int> &list_columns_var_in,
		     vector <unsigned int> &is_var_in)
{
  
  unsigned int numVars = 0;
  list_columns_var_in.resize(is_var_in.size());
  vector<unsigned int>::iterator it,it2,it3;
  it2=list_columns_var_in.begin();
  it3=is_var_in.end();
  unsigned int i=0;
  for(it=is_var_in.begin();it<it3;it++){
    if(*it){
      *it2=i;
      it2++;
      numVars++;
    }
    i++;
  }
  list_columns_var_in.resize(numVars);
}

gsl_matrix* get_X_gam(vector <unsigned int> &list_columns_var_in,
		      gsl_matrix *mat_X)
{
  unsigned int n_vars_in=list_columns_var_in.size();
  if (DEBUG){
    cout << "\tIn get_X_gam" << endl
	 << "\t#variables in " << list_columns_var_in.size() << endl;
    for(unsigned int i=0;i<list_columns_var_in.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var_in.size() << ") : " << list_columns_var_in[i] << endl;
    }
  }
  unsigned int nX=mat_X->size1;
  gsl_matrix *X_gam=gsl_matrix_alloc(nX,n_vars_in);
  for(unsigned int col=0;col<n_vars_in;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var_in[col]);
    //Setting the col^th column of X_gam
    gsl_matrix_set_col (X_gam,col,&current_col.vector);
  }
  return X_gam;
}

gsl_matrix* get_X_reduced_and_constant(vector <unsigned int> &list_columns_var,
				       gsl_matrix *mat_X)
{
  unsigned int n_vars=list_columns_var.size();
  if (DEBUG){
    cout << "\tIn get_X_reduced_and_constant" << endl
	 << "\t#variables in " << list_columns_var.size() << endl;
    for(unsigned int i=0;i<list_columns_var.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var.size() << ") : " << list_columns_var[i] << endl;
    }
  }

  unsigned int nX=mat_X->size1;
  gsl_matrix *X_red=gsl_matrix_alloc(nX,(n_vars)+1);
  
  //The first column is set to 1: the constant term of the regression
  gsl_vector *temp_vector=gsl_vector_alloc(nX);
  gsl_vector_set_all(temp_vector,1.0);
  gsl_matrix_set_col (X_red,0,temp_vector);

  gsl_vector_free(temp_vector);

  //Setting the other columns to the list of vars in
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var[col]);
    //Setting the col^th column of X_red
    gsl_matrix_set_col (X_red,col+1,&current_col.vector);
  }
  return X_red;
}

gsl_matrix* get_X_reduced(vector <unsigned int> &list_columns_var,
			  gsl_matrix *mat_X)
{
  unsigned int n_vars=list_columns_var.size();
  if (DEBUG){
    cout << "\tIn get_X_reduced" << endl
	 << "\t#variables in " << list_columns_var.size() << endl;
    for(unsigned int i=0;i<list_columns_var.size();i++){
      cout << "\tVariable position (" << i+1 << "/ " << list_columns_var.size() << ") : " << list_columns_var[i] << endl;
    }
  }

  if (DEBUG)
  {
      cout << "Before accessing mat_X->size1" << endl;
  }
  unsigned int nX=mat_X->size1;
  if (DEBUG)
  {
      cout << mat_X->size1 << endl;
  }
  gsl_matrix *X_red=gsl_matrix_alloc(nX,(n_vars));
  if (DEBUG)
  {
      cout << "X_red allocated" << endl;
  }
  //Setting the other columns to the list of vars in
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,list_columns_var[col]);
    if (DEBUG)
    {
        cout << "Column copied in current_col: " << col << endl;
    }
    //Setting the col^th column of X_red
    gsl_matrix_set_col (X_red,col,&current_col.vector);
    if (DEBUG)
    {
        cout << "Column set to X_red: " << col << endl;
    }

  }
  return X_red;
}

gsl_matrix* get_sub_matrix_col(gsl_matrix *mat_X,
			       size_t first_col,
			       size_t last_col)
{
  unsigned int n_vars=last_col-first_col;
  unsigned int nX=mat_X->size1;
  gsl_matrix *X_red=gsl_matrix_alloc(nX,n_vars);
  
  for(unsigned int col=0;col<n_vars;col++){
    //Copying the the col^th variable (at position list_columns_var_in[col] in mat_X)
    gsl_vector_view current_col = gsl_matrix_column (mat_X,first_col+col);
    gsl_matrix_set_col (X_red,col,&current_col.vector);
  }
  return X_red;
}

void fill_sub_matrix_col(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_col,
			 size_t last_col)
{
  unsigned int n_vars=last_col-first_col;
  
  for(unsigned int col=0;col<n_vars;col++){
    gsl_vector_view current_col = gsl_matrix_column (mat_X,first_col+col);
    gsl_matrix_set_col (X_red,col,&current_col.vector);
  }
}

gsl_matrix* get_sub_matrix_row(gsl_matrix *mat_X,
			       size_t first_row,
			       size_t last_row)
{
  unsigned int n_row=last_row-first_row;
  unsigned int nY=mat_X->size2;
  gsl_matrix *X_red=gsl_matrix_alloc(n_row,nY);
  
  for(unsigned int row=0;row<n_row;row++){
    gsl_vector_view current_row = gsl_matrix_row(mat_X,first_row+row);
    gsl_matrix_set_row(X_red,row,&current_row.vector);
  }
  return X_red;
}

void fill_sub_matrix_row(gsl_matrix *X_red,
			 gsl_matrix *mat_X,
			 size_t first_row,
			 size_t last_row)
{
  unsigned int n_row=last_row-first_row;
  
  for(unsigned int row=0;row<n_row;row++){
    gsl_vector_view current_row = gsl_matrix_row(mat_X,first_row+row);
    gsl_matrix_set_row(X_red,row,&current_row.vector);
  }
}

void display_gsl_matrix(gsl_matrix *M)
{
  cout << "nb_rows " << M->size1 << " -- nb_columns " << M->size2 << endl;
    for(unsigned int i=0;i<M->size1;i++){
      for(unsigned int j=0;j<M->size2;j++){
	cout << M->data[i*(M->size2)+j] << " ";
      }
      cout << endl;
    }
}

void display_gsl_vector(gsl_vector *V)
{
  cout << "vector size " << V->size << endl;
  for(unsigned int i=0;i<V->size;i++){
    cout << V->data[i] << " ";
  }
  cout << endl;
}

void display_gsl_perm(gsl_permutation *P)
{
  cout << "perm size " << P->size << endl;
  for(unsigned int i=0;i<P->size;i++){
    cout << P->data[i] << " ";
  }
  cout << endl;
}

void display_vector_int(vector < unsigned int> &vector)
{
  cout << "Vector size " << vector.size() << endl;
  for(unsigned int i=0;i<vector.size();i++){
    cout << vector[i] << " ";
  }
  cout << endl;
}

void display_matrix_var_dim(vector < vector <unsigned int> > &M)
{
  for(unsigned int row=0;row<M.size();row++){
    cout << "Row #" << row+1 << " -- nb columns " << M[row].size() << endl;
    for(unsigned int col=0;col<M[row].size();col++){
      cout << M[row][col] << " ";
    }
    if(M[row].size()>0){
      cout << endl;
    }
  }


}

unsigned int sum_line_std_mat(vector < vector <unsigned int> > &M,
			      unsigned int line)
{
  unsigned int tmp_sum=0;
  vector<unsigned int>::iterator it,itEnd;
  itEnd=M[line].end();
  for(it=M[line].begin();it<itEnd;++it){
    tmp_sum+=(*it);
  }
  return tmp_sum;
}

void display_result_per_sweep(vector < vector <unsigned int> > vect_gam,
			     vector < unsigned int > chain_idx,
			     Double_Matrices mat_log_marg,
			     Double_Matrices mat_log_cond_post,
			     unsigned int sweep,
			     Temperatures *t_tun)
{
  for(unsigned int chain=0;chain<vect_gam.size();chain++){
    unsigned int pos_chain=chain_idx[chain];
    unsigned int n_vars_in=sum_line_std_mat(vect_gam,
					    pos_chain);
    cout << "Chain #" << chain+1
	 << " -- position: " << pos_chain+1
	 << " -- log_marg " << mat_log_marg.matrix[pos_chain][sweep]
	 << " -- log_cond_post " << mat_log_cond_post.matrix[pos_chain][sweep]
	 << " -- n_vars_in " << n_vars_in
	 << " -- t_tun.t " << (*t_tun).t[chain]
	 << endl;
    for(unsigned int var=0;var<vect_gam[chain].size();var++){
      if(vect_gam[pos_chain][var]==1){
	cout << var << " ";
      }
    }
    cout << endl;
  }
}

void display_summary_result_per_sweep(vector < vector <unsigned int> > vect_gam,
				     vector < unsigned int > chain_idx,
				     Double_Matrices mat_log_marg,
				     Double_Matrices mat_log_cond_post,
				     unsigned int sweep,
				     Temperatures *t_tun,
				     const unsigned int& nConfounders)
{
  for(unsigned int chain=0;chain<vect_gam.size();chain++){
    unsigned int pos_chain=chain_idx[chain];
    unsigned int n_vars_in=sum_line_std_mat(vect_gam,
					    pos_chain)-nConfounders;
    cout << "Chain #" << chain+1
	 << " -- position: " << pos_chain+1
	 << " -- log_marg " << mat_log_marg.matrix[pos_chain][sweep]
	 << " -- log_cond_post " << mat_log_cond_post.matrix[pos_chain][sweep]
	 << " -- n_vars_in " << n_vars_in
	 << " -- t_tun.t " << (*t_tun).t[chain]
	 << endl;  
  }
}

void print_main_results_per_sweep(ofstream &f_out,
				 vector < vector <unsigned int> > vect_gam,
				 vector < unsigned int > chain_idx,
				 Double_Matrices mat_log_marg,
				 Double_Matrices mat_log_cond_post,
				 unsigned int sweep)
{
  unsigned int pos_chain=chain_idx[0];
  unsigned int n_vars_in=sum_line_std_mat(vect_gam,
					  pos_chain);
  f_out << sweep << "\t"
	<< n_vars_in << "\t";
  f_out << setprecision(13) << mat_log_marg.matrix[pos_chain][sweep] << "\t";
  f_out << setprecision(13) << mat_log_cond_post.matrix[pos_chain][sweep] << "\t\t";
  
  
  for(unsigned int var=0;var<vect_gam[0].size();var++){
    if(vect_gam[pos_chain][var]==1){
      f_out << var+1 << " ";
    }
  }
  f_out << endl;
}

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
					  bool HistoryFlag)
{
  unsigned int pos_chain=chain_idx[0];
  unsigned int n_vars_in=sum_line_std_mat(vect_gam,
					  pos_chain);

  if(HistoryFlag){
    f_out << sweep << "\t"
	  << n_vars_in << "\t";
    f_out << setprecision(13) << mat_log_marg.matrix[pos_chain][sweep] << "\t";
    f_out << setprecision(13) << mat_log_cond_post.matrix[pos_chain][sweep] << "\t";
  }
  if(sweep>burn_in){
    List_models.push_back(vector<unsigned int>());
    unsigned int last_pos_in_List=List_models.size()-1;
    List_models[last_pos_in_List].push_back(0);//room for the #of variables in
  }
  unsigned int count_n_vars_in=0;

  for(unsigned int var=0;var<vect_gam[0].size();var++){
    if(vect_gam[pos_chain][var]==1){
      if(HistoryFlag){
	f_out << var+1 << " ";
      }
      if(sweep>burn_in){
	unsigned int last_pos_in_List=List_models.size()-1;
	List_models[last_pos_in_List].push_back(var);
	count_n_vars_in++;
      }
    }
  }
  if(sweep>burn_in){
    unsigned int last_pos_in_List=List_models.size()-1;
    List_models[last_pos_in_List][0]=count_n_vars_in;
  }
  if(HistoryFlag){
    f_out << endl;
    f_out_log_cond_post << sweep << "\t";
    f_out_n_vars_in << sweep << "\t";
    f_out_n_models_visited << sweep << "\t" << nModelsVisited << endl;
    for(unsigned int chain=0;chain<vect_gam.size();chain++){
      unsigned int current_chain=chain_idx[chain];
      unsigned int n_vars_in_per_chain=sum_line_std_mat(vect_gam,
							current_chain);
      f_out_n_vars_in << n_vars_in_per_chain << "\t";
      f_out_log_cond_post << setprecision(13) << mat_log_cond_post.matrix[current_chain][sweep] << "\t";
    }
    f_out_n_vars_in << endl;
    f_out_log_cond_post << endl;
  }

}

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
                                          bool HistoryFlag)
{
  unsigned int pos_chain=chain_idx[0];
  unsigned int n_vars_in=sum_line_std_mat(vect_gam,
                                          pos_chain);

  if(HistoryFlag){
    ss_out << sweep << "\t"
          << n_vars_in << "\t";
    ss_out << setprecision(13) << mat_log_marg.matrix[pos_chain][sweep] << "\t";
    ss_out << setprecision(13) << mat_log_cond_post.matrix[pos_chain][sweep] << "\t";
  }
  if(sweep>burn_in){
    List_models.push_back(vector<unsigned int>());
    unsigned int last_pos_in_List=List_models.size()-1;
    List_models[last_pos_in_List].push_back(0);//room for the #of variables in
  }
  unsigned int count_n_vars_in=0;

  for(unsigned int var=0;var<vect_gam[0].size();var++){
    if(vect_gam[pos_chain][var]==1){
      if(HistoryFlag){
        ss_out << var+1 << " ";
      }
      if(sweep>burn_in){
        unsigned int last_pos_in_List=List_models.size()-1;
        List_models[last_pos_in_List].push_back(var);
        count_n_vars_in++;
      }
    }
  }
  if(sweep>burn_in){
    unsigned int last_pos_in_List=List_models.size()-1;
    List_models[last_pos_in_List][0]=count_n_vars_in;
  }
  if(HistoryFlag){
    ss_out << endl;
    ss_out_log_cond_post << sweep << "\t";
    ss_out_n_vars_in << sweep << "\t";
    ss_out_n_models_visited << sweep << "\t" << nModelsVisited << endl;
    for(unsigned int chain=0;chain<vect_gam.size();chain++){
      unsigned int current_chain=chain_idx[chain];
      unsigned int n_vars_in_per_chain=sum_line_std_mat(vect_gam,
                                                        current_chain);
      ss_out_n_vars_in << n_vars_in_per_chain << "\t";
      ss_out_log_cond_post << setprecision(13) << mat_log_cond_post.matrix[current_chain][sweep] << "\t";
    }
    ss_out_n_vars_in << endl;
    ss_out_log_cond_post << endl;
  }

}


void saveResumeFile(fstream &fResume, FILE *fRNG, unsigned int sweep,
                    double g,
                    //unsigned int* shuffleYIndex,
                    vector<unsigned int> shuffleYIndex,unsigned int nY,
                    double cumG, unsigned int countG,
                    gsl_vector *vectRMSE,
                    Temperatures *tTun,double gLs,DR *currDR,
                    vector<vector<unsigned int> > vectGam,
                    vector < unsigned int > chainIndex,
                    unsigned int nConfounders,
                    unsigned int pY,
                    gsl_rng *RandomNumberGenerator){

  // Store the random number state
  writeRNG(fRNG,RandomNumberGenerator);
  fResume << sweep << endl;
  fResume << setiosflags(ios::fixed) << setprecision(10) << g << endl;
  for(unsigned int i=0;i<nY;i++){
    fResume << shuffleYIndex[i] << endl;
  }
  fResume << cumG << endl;
  fResume << countG << endl;
  for(unsigned int j=0;j<pY;j++){
    fResume << vectRMSE->data[j] << endl;
  }
  unsigned int nChains = vectGam.size();
  for(unsigned int j=0;j<nChains;j++){
    fResume << setiosflags(ios::fixed) << setprecision(10) << tTun->t[j] << endl;
  }
  fResume << setiosflags(ios::fixed) << setprecision(10) << tTun->b_t << endl;
  fResume << setiosflags(ios::fixed) << setprecision(10) << gLs << endl;

  fResume << currDR->nb_calls << endl;
  fResume << currDR->nb_calls_adj << endl;
  for(unsigned int j=0;j<nChains;j++){
    for(unsigned int k=0;k<nChains;k++){
      fResume << currDR->mat_moves_proposed[j][k] << endl;
    }
  }

  for(unsigned int j=0;j<nChains;j++){
    for(unsigned int k=0;k<nChains;k++){
      fResume << currDR->mat_moves_accepted[j][k] << endl;
    }
  }

  for(unsigned int c=0;c<nChains;c++){
    vector<unsigned int> tmpVec;
    for(unsigned int j=0;j<vectGam[0].size();j++){
      if(j<nConfounders){
        continue;
      }
      if(vectGam[chainIndex[c]][j]==1){
        tmpVec.push_back(j);
      }
    }
    fResume << tmpVec.size() << endl;
    fResume << chainIndex[c] << endl;
    for(unsigned int j=0;j<tmpVec.size();j++){
      fResume << tmpVec[j] << endl;
    }
  }
}
