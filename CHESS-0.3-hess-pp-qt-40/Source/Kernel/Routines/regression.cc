/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from regression.cc in the ESS++ program
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

#include "regression.h"
#include "struc.h"

#define DEBUG 0

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF


using namespace std;

void perform_Univariate_regr(vector < vector <unsigned int> > &Gam_regr,
			     Double_Matrices mat_Y,
			     Double_Matrices mat_X,
			     double alpha)
{
  unsigned int nX=mat_X.nb_rows;// # observations
  unsigned int pX=mat_X.nb_columns;// # SNPs
  unsigned int nY=mat_Y.nb_rows;// # observations 
  unsigned int pY=mat_Y.nb_columns;// # outcomes
  
  double *vect_X= new double[nX];
  double *vect_Y= new double[nY];
  double *vect_Y_sim= new double[nY];
  double nu=max(0.0,(double)(nX)-2.0);
  double t_val=gsl_cdf_tdist_Pinv((1.0-alpha),nu);
  
  vector <unsigned int> retained_idx;
  vector <double> retained_pvalue;

  cout << "nu " << nu 
       << " -- alpha " << alpha
       << " -- tval " << t_val << endl;
  
  for(unsigned int Current_outcome=0;Current_outcome<pY;Current_outcome++){
    if(DEBUG){
      cout << "Current outcome " << Current_outcome << endl << endl;
    }
    //Setting the Y vector & computing mean
    
    double tmp_sum_Y=0.0;
    for(unsigned int Myrow=0;Myrow<nY;Myrow++){
      vect_Y[Myrow]=mat_Y.matrix[Myrow][Current_outcome];
      tmp_sum_Y+=vect_Y[Myrow];
    }
    double Y_bar=tmp_sum_Y/nY;
    if(DEBUG){
      cout << "\tY_bar= " << Y_bar << endl;
    }    
    for(unsigned int Current_SNP=0;Current_SNP<pX;Current_SNP++){
      if(DEBUG){
	cout << "\tCurrent SNP " << Current_SNP << endl << endl;
      }
      
      //Setting the X vector
      for(unsigned int Myrow=0;Myrow<nY;Myrow++){
	vect_X[Myrow]=mat_X.matrix[Myrow][Current_SNP];
      }
      //Defining output variables for the univariate regression

      double constant,coeff,cov00,cov01,cov11,sumsq;
      //Fitting the regression model on each XxY combonation

      gsl_fit_linear(vect_X,1,
		     vect_Y,1,
		     nY,
		     &constant,&coeff,
		     &cov00,&cov01,&cov11,
		     &sumsq);

      double RSS=0.0;
      double TSS=0.0;
      
      for(unsigned int Myrow=0;Myrow<nY;Myrow++){
	vect_Y_sim[Myrow]=constant+coeff*(double)(vect_X[Myrow]);
	RSS+=pow((vect_Y_sim[Myrow]-Y_bar),2);
	TSS+=pow((vect_Y[Myrow]-Y_bar),2);
      } 
      
      double S=sumsq/((double)(nX)-2.0);//S: variance of the error. 
      double se=sqrt(cov11);
      double F_stat=RSS/S;
 
      if(DEBUG){
	cout << "\tConstant " << constant << endl
	     << "\tcoeff " << coeff << endl
	     << "\tcov00 " << cov00 << " -- cov01 " << cov01 << " -- cov11 " << cov11 << endl
	     << "\tsum of square " << sumsq << endl
	     << "\tvariance of the error " << S << endl
	     << "\tY_bar " << Y_bar << endl
	     << "\tRSS= " << RSS << endl
	     << "\tTSS= " << TSS << endl
	     << "\tR2= " << 1-(sumsq/TSS) << endl
	     << "\tse= " << se << endl
	     << "\tF_stat=" << F_stat << endl
	     << "\tp_val= " << gsl_cdf_fdist_Q(F_stat,1,nu) << endl;
      }

 
      //Focus is made on SNP whose CI does not include 0
      
      double coeff_inf=coeff-t_val*se;
      double coeff_sup=coeff+t_val*se;
      double current_p_value=gsl_cdf_fdist_Q(F_stat,1,nu);
      
      if(coeff_inf*coeff_sup>0.0){
	retained_pvalue.push_back(current_p_value);
 	retained_idx.push_back(Current_SNP);
      }
      
    }//end of for SNP
    unsigned int nb_retained=retained_idx.size();
    if(DEBUG){
      cout << "Current_outcome " << Current_outcome 
	   << " -- Permutation size " << nb_retained << endl;
    }
    if( nb_retained>0){
      double * temp_pointer =&retained_pvalue[0];
      size_t *idx_sorted= new size_t[nb_retained];
      gsl_sort_index (idx_sorted,
		      temp_pointer,
		      1,
		      nb_retained);
      
      for(unsigned int i=0;i<nb_retained;i++){
	Gam_regr[Current_outcome].push_back(retained_idx[idx_sorted[i]]);
	if(DEBUG){
	  cout << "pos_in_sorted= " << i
	       << " -- real_idx= " << retained_idx[idx_sorted[i]] 
	       << " -- p_val= " << retained_pvalue[idx_sorted[i]]  <<  endl;
	}
      }
      
      delete [] idx_sorted;
    }

    retained_pvalue.clear();
    retained_idx.clear();

  }//end of for outcome
  
  
  delete []  vect_X;
  delete []  vect_Y;
  delete []  vect_Y_sim;

}

unsigned int get_rank_from_R(gsl_matrix *M,
			     double tolerance)
{
  unsigned int rank=0;
  unsigned int min_size=GSL_MIN(M->size1,M->size2); 
  double ref_val=tolerance*fabs(M->data[0]);
  for(unsigned int i=0;i<min_size;i++){
    if(fabs(M->data[i*M->size2+i])>ref_val){
	rank++;
    }
  }
  return(rank);
}

gsl_vector* get_vect_elem_square(gsl_matrix *M)
{
  gsl_vector *M2=gsl_vector_calloc(M->size2);

  for(unsigned int row=0;row<M->size1;row++){
    for(unsigned int col=0;col<M->size2;col++){
      M2->data[col]+=pow(M->data[row*M->size2+col],2.0);
    }
  }
  return M2;
}

void get_vect_beta_var_out(gsl_vector *vect_beta_var_out,
			   gsl_vector *xx,
			   gsl_matrix *xr,
			   gsl_vector *residuals)
{
  gsl_blas_dgemv(CblasTrans,
		 1.0,
		 xr,
		 residuals,
		 0.0,
		 vect_beta_var_out);

  for(unsigned int col=0;col<vect_beta_var_out->size;col++){
    vect_beta_var_out->data[col]/=xx->data[col];
  }  
}

void get_se_var_out(gsl_vector *S2,
		    gsl_matrix *xr,
		    gsl_vector *vect_residuals,
		    gsl_vector *vect_beta,
		    gsl_vector *xx,
		    double df2)
{


  if(DEBUG){
    cout << "Size Check: xr " << xr->size1 << "x" << xr->size2 << endl
	 << "S2 " << S2->size << endl
	 << "residuals " << vect_residuals->size << endl
	 << "vect_beta " << vect_beta->size << endl;
  }
  for(unsigned int col=0;col<xr->size2;col++){
    if(DEBUG){
      cout << "col " << col << endl;
    }
    for(unsigned int row=0;row<xr->size1;row++){
      if(DEBUG){
	cout << "\trow " << row << " -- S2init= " << S2->data[col];
      }
      S2->data[col]+=pow((vect_residuals->data[row]-
			  (xr->data[row*xr->size2+col]*vect_beta->data[col])),2)/df2;
      if(DEBUG){
	cout << " -- xr= " <<  xr->data[row*xr->size2+col]
	     << " -- b= " << vect_beta->data[col]
	     << " -- Test1= " << (xr->data[row*xr->size2+col]*vect_beta->data[col])
	     << " -- Test2= " << (vect_residuals->data[row]- (xr->data[row*xr->size2+col]*vect_beta->data[col]))
	     << " -- Test3= " << pow((vect_residuals->data[row]-(xr->data[row*xr->size2+col]*vect_beta->data[col])),2)
	     << " -- Test4= " << S2->data[col] << endl;
      }
    }
  }
  for(unsigned int col=0;col<xr->size2;col++){
    S2->data[col]=sqrt( S2->data[col]/xx->data[col]); 
  }
}

gsl_vector *get_se(gsl_matrix *mat_R_red,
		   double RMSE)
{
  
  gsl_matrix *mat_R_inv=gsl_matrix_calloc(mat_R_red->size2,mat_R_red->size2);
  gsl_matrix_set_identity(mat_R_inv);
  
  //Calculating R_inv
  gsl_blas_dtrsm(CblasLeft,
		 CblasUpper,
		 CblasNoTrans,
		 CblasNonUnit,
		 1.0,
		 mat_R_red,
		 mat_R_inv);
 
  if(DEBUG){ 
    cout << "Mat R_inv" << endl;
    display_gsl_matrix(mat_R_inv);
  }
  gsl_vector *SE=gsl_vector_calloc(mat_R_inv->size1);
  for(unsigned int row=0;row<mat_R_inv->size1;row++){
    double temp_sum=0.0;
    for(unsigned int col=0;col<mat_R_inv->size2;col++){
      temp_sum+=pow(mat_R_inv->data[row*mat_R_inv->size2+col],2);
    }
    SE->data[row]=sqrt(temp_sum)*RMSE;
  }
  gsl_matrix_free(mat_R_inv);
  return SE;
}

void get_mat_Q_red_R_red_beta_in_resid(gsl_matrix *mat_X_gam,
				       gsl_vector *current_outcome,
				       gsl_vector *vect_residuals,
				       gsl_vector *vect_beta_var_in,
				       gsl_matrix *mat_Q_red,
				       gsl_matrix *mat_R_red)
{
  

  //Declaring the vectors for QR decomposition

  gsl_vector *tau=gsl_vector_alloc(GSL_MIN(mat_X_gam->size1,mat_X_gam->size2));
  gsl_matrix *mat_Q=gsl_matrix_alloc(mat_X_gam->size1,mat_X_gam->size1);
  gsl_matrix *mat_R=gsl_matrix_alloc(mat_X_gam->size1,mat_X_gam->size2);
  
  //Step 1. QR decomposition
  //Mat_X_gam is the synthetic QR
  gsl_linalg_QR_decomp(mat_X_gam,tau);
  
  //Step 2. Getting the LS regression coefficients
  gsl_linalg_QR_lssolve(mat_X_gam,
			tau,
			current_outcome,
			vect_beta_var_in,
			vect_residuals);
  //Step 3. Extracting the full Q and R matrices
  gsl_linalg_QR_unpack(mat_X_gam,
		       tau,
		       mat_Q,
		       mat_R);
  //step 4. get economy size Q and R.
  fill_sub_matrix_col(mat_Q_red,
		      mat_Q,
		      0,
		      mat_X_gam->size2);

  fill_sub_matrix_row(mat_R_red,
		      mat_R,
		      0,
		      mat_X_gam->size2);


  gsl_matrix_free(mat_Q);
  gsl_matrix_free(mat_R);
  gsl_vector_free(tau);
  

}

void get_X_residuals(gsl_matrix *mat_X_gam_bar,
		     gsl_matrix *mat_Q_red)
{
  
  //Step 1: computing Q.QT
  gsl_matrix *Q_QT=gsl_matrix_calloc(mat_Q_red->size1,mat_Q_red->size1);
  
  gsl_blas_dgemm(CblasNoTrans,
		 CblasTrans,
		 1.0,
		 mat_Q_red,
		 mat_Q_red,
		 0.0,
		 Q_QT);
  if(DEBUG){ 
    cout << "Q_QT" << endl;
    display_gsl_matrix(Q_QT);
  }
  //Step 2: computing Q.QT.X_gam_bar
  gsl_matrix *Q_QT_X_gam_bar=gsl_matrix_calloc(Q_QT->size1,mat_X_gam_bar->size2);
  

  // VERY SLOW:
  //cout << "Before gsl_blas_dgemm 2" << endl;
  //cout << "Q_QT: (" << Q_QT->size1 << "," << Q_QT->size2 << ")" << endl;
  //cout << "mat_X_gam_bar: (" << mat_X_gam_bar->size1 << "," << mat_X_gam_bar->size2 << ")" << endl;
  //cout << "Cible: Q_QT_X_gam_bar: (" << Q_QT_X_gam_bar->size1 << "," << Q_QT_X_gam_bar->size2 << ")" << endl;
  gsl_blas_dgemm(CblasNoTrans,
		 CblasNoTrans,
		 1.0,
		 Q_QT,
		 mat_X_gam_bar,
		 0.0,
		 Q_QT_X_gam_bar);

  gsl_matrix_free(Q_QT);
   if(DEBUG){ 
     cout << "QxQTxX_gam_bar" << endl; 
     display_gsl_matrix(Q_QT_X_gam_bar);
   }
  //step 3: get X_gam_bar-Q.QT.X_gam_bar: result store in X_gam_bar

  gsl_matrix_sub(mat_X_gam_bar,
		 Q_QT_X_gam_bar);
  
  gsl_matrix_free(Q_QT_X_gam_bar);
  
  

}

void compute_beta_SE_vars_out(gsl_vector *vect_beta_var_out,
			      gsl_vector *vect_SE_var_out,
			      gsl_matrix *mat_X_gam_bar,
			      gsl_matrix *mat_Q_red,
			      gsl_vector *vect_residuals,
			      double df2)
{

  //Step1: get X_residuals= X_gam_bar-(Q.Q_T.X_gam_bar)
  //!!!!!!!!!!!!!!!!!!! WARNING !!!!!!!!!!!!!!!!!!!!!!!
  // RESULT STORED IN mat_X_gam_bar
  
  // VERY SLOW:
  get_X_residuals(mat_X_gam_bar,
		  mat_Q_red);

  if(DEBUG){ 
    cout << "mat X_residuals" << endl; 
    display_gsl_matrix(mat_X_gam_bar);
  }

  gsl_vector *vect_X_res_square=get_vect_elem_square(mat_X_gam_bar);
  if(DEBUG){ 
      cout << "vect xx" << endl; 
      display_gsl_vector(vect_X_res_square);
  }


  get_vect_beta_var_out(vect_beta_var_out,
			vect_X_res_square,
			mat_X_gam_bar,
			vect_residuals);
  if(DEBUG){ 
    cout << "vect beta var out" << endl; 
    display_gsl_vector(vect_beta_var_out);
  }
  
  //Step 2: get SE for variables out
  
  get_se_var_out(vect_SE_var_out,
		 mat_X_gam_bar,
		 vect_residuals,
		 vect_beta_var_out,
		 vect_X_res_square,
		 df2);

  if(DEBUG){ 
    cout << "vect SE var out" << endl; 
    display_gsl_vector(vect_SE_var_out);
  }

  gsl_vector_free(vect_X_res_square);

}

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
			       unsigned int outcome)
{

  

   ///////////////////////////////////////////
    //***************************************//
    // Step 1: calculating coefficients and  //
    //   Standard errors for variables in    //
    //         working with mat_X_gam        //
    //***************************************//
    ///////////////////////////////////////////

    //cout << "Before gsl_matrix_calloc" << endl;
    //Step 1.1: Getting regression coefficients/residuals and QR decompostion
    gsl_matrix *mat_Q_red=gsl_matrix_calloc(mat_X_gam->size1,mat_X_gam->size2);
    gsl_matrix *mat_R_red=gsl_matrix_calloc(mat_X_gam->size2,mat_X_gam->size2);
    gsl_vector *vect_beta_var_in=gsl_vector_calloc(mat_X_gam->size2);
  

    // A bit slow:
    get_mat_Q_red_R_red_beta_in_resid(mat_X_gam,
				      current_outcome,
				      vect_residuals,
				      vect_beta_var_in,
				      mat_Q_red,
				      mat_R_red);

    size_t actual_rank=get_rank_from_R(mat_X_gam,
				       tolerance);
    
    double dfe=mat_X_gam->size1-(double)(actual_rank);

    double RSS=gsl_stats_tss(vect_residuals->data,
			     1,
			     vect_residuals->size);
    double RMSE=sqrt(RSS/dfe);

    vect_RMSE->data[outcome]=RMSE;


    //Step 1.2: calculating SE vars_in.
    gsl_vector *vect_SE_var_in=get_se(mat_R_red,
				      RMSE);
    ///////////////////////////////////////////
    //***************************************//
    // Step 2: calculating coefficients and  //
    //   Standard errors for variables out   //
    //      working with mat_X_gam_bar       //
    //***************************************//
    ///////////////////////////////////////////
    double df2=GSL_MAX(0.0,dfe-1.0);
   
    gsl_vector *vect_SE_var_out=gsl_vector_calloc(mat_X_gam_bar->size2);
    gsl_vector *vect_beta_var_out=gsl_vector_calloc(mat_X_gam_bar->size2);


    // REALLY SLOW:
    compute_beta_SE_vars_out(vect_beta_var_out,
			     vect_SE_var_out,
			     mat_X_gam_bar,
			     mat_Q_red,
			     vect_residuals,
			     df2);
    gsl_matrix_free(mat_Q_red);
    gsl_matrix_free(mat_R_red);
    
    ///////////////////////////////////////////
    //***************************************//
    //     Step 3: calculating p-values      //
    //            for all variables          //
    //                                       //
    //***************************************//
    ///////////////////////////////////////////

    //For variables in
    for(unsigned int col=0;col<list_columns_X_gam.size();col++){
      unsigned int pos_current_var_in=list_columns_X_gam[col];
      double temp_t_stat=-fabs((vect_beta_var_in->data[col+1])/(vect_SE_var_in->data[col+1]));
      vect_p_value->data[pos_current_var_in]=2.0*gsl_cdf_tdist_P(temp_t_stat,dfe);
      vect_beta_full->data[pos_current_var_in]=vect_beta_var_in->data[col+1];
      vect_SE_full->data[pos_current_var_in]=vect_SE_var_in->data[col+1];
    }

    //For Variables out
    for(unsigned int col=0;col<list_columns_X_gam_bar.size();col++){
      unsigned int pos_current_var_out=list_columns_X_gam_bar[col];
      double temp_t_stat=-abs((vect_beta_var_out->data[col])/(vect_SE_var_out->data[col]));
      vect_p_value->data[pos_current_var_out]=2.0*gsl_cdf_tdist_P(temp_t_stat,dfe-1.0);
      vect_beta_full->data[pos_current_var_out]=vect_beta_var_out->data[col];
      vect_SE_full->data[pos_current_var_out]=vect_SE_var_out->data[col];
     }
    
        



    if(DEBUG){
      cout << endl << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl
	   << "Multivariate regression, outcome # " << outcome << endl;
      cout << "Best fit coefficient estimate var in" << endl;
      display_gsl_vector(vect_beta_var_in);
      cout << endl;
      cout << "SE var in" << endl;
      display_gsl_vector(vect_SE_var_in);
      cout << endl;
      cout << "Best fit coefficient estimate var out" << endl;
      display_gsl_vector(vect_beta_var_out);
      cout << endl;
      cout << "SE var out" << endl;
      display_gsl_vector(vect_SE_var_out);
      cout << endl;
      cout << "p-value full" << endl;
      display_gsl_vector(vect_p_value);
      cout << endl;
      cout << "beta full" << endl;
      display_gsl_vector(vect_beta_full);
      cout << endl;
      cout << endl << endl; 
    }
    
    
    gsl_vector_free(vect_beta_var_in);
    gsl_vector_free(vect_SE_var_in);
    gsl_vector_free(vect_beta_var_out);
    gsl_vector_free(vect_SE_var_out);
    
}

int update_is_var_in(vector < unsigned int > &is_var_in,
		     vector < unsigned int > &list_columns_X_gam,
		     vector < unsigned int > &list_columns_X_gam_bar,
		     gsl_vector *vect_p_value,
		     double Pvalue_enter,
		     double Pvalue_remove,
		     unsigned int loop,
		     unsigned int n_loop_max,
		     unsigned int nConfounders)
{
  size_t idx_min=0;
  size_t idx_max=0;
  size_t pos_in_X_min=0;
  size_t pos_in_X_max=0;
  double MyMin=0.0;
  double MyMax=0.0;
  unsigned int n_vars_in=list_columns_X_gam.size();
  unsigned int n_vars_out=list_columns_X_gam_bar.size();
  int included=0;
  int removed=0;
  int stop=0;

  if(n_vars_out>0 && n_vars_in>0){//some vars are out
    //Getting p_vals_in
    gsl_vector *p_value_in=gsl_vector_calloc(n_vars_in);
      for(unsigned int current_var_in=0;current_var_in<n_vars_in;current_var_in++){
	unsigned int rank_var_in=list_columns_X_gam[current_var_in];
	p_value_in->data[current_var_in]=vect_p_value->data[rank_var_in];
      }
   //Getting p_vals_out
    gsl_vector *p_value_out=gsl_vector_calloc(n_vars_out);
    for(unsigned int current_var_out=0;current_var_out<n_vars_out;current_var_out++){
      unsigned int rank_var_out=list_columns_X_gam_bar[current_var_out];
      p_value_out->data[current_var_out]=vect_p_value->data[rank_var_out];
    }
    idx_min=gsl_vector_min_index(p_value_out);
    pos_in_X_min=list_columns_X_gam_bar[idx_min];
    idx_max=gsl_vector_max_index(p_value_in);
    pos_in_X_max=list_columns_X_gam[idx_max];
    MyMin=p_value_out->data[idx_min];
    MyMax=p_value_in->data[idx_max];
    if(DEBUG){
      cout << "idx_min " << idx_min 
	   << " -- min " << MyMin 
	   << " -- pos_in_X_min " << pos_in_X_min 
	   << " -- idx_max " << idx_max 
	   << " -- pos_in_X_max " << pos_in_X_max
	   << " -- max " << MyMax
	   << endl;
    }
    if(MyMin<Pvalue_enter){//the corresponding variable is included
      included=1;
      is_var_in[pos_in_X_min]=1;
    }
    else if(MyMax>Pvalue_remove&&pos_in_X_max>=nConfounders){ // do not remove confounders
      removed=1;
      is_var_in[pos_in_X_max]=0;
    }
    else{
    }
    gsl_vector_free(p_value_in);
    gsl_vector_free(p_value_out);
    
  }
  else{//all are in OR all are out
    if(n_vars_in==0){//all are out
      gsl_vector *p_value_out=gsl_vector_calloc(n_vars_out);
      for(unsigned int current_var_out=0;current_var_out<n_vars_out;current_var_out++){
	unsigned int rank_var_out=list_columns_X_gam_bar[current_var_out];
	p_value_out->data[current_var_out]=vect_p_value->data[rank_var_out];
      }
      idx_min=gsl_vector_min_index(p_value_out);
      pos_in_X_min=list_columns_X_gam_bar[idx_min];
      MyMin=p_value_out->data[idx_min];
      
      if(DEBUG){
	cout << "idx_min " << idx_min 
	     << " -- min " << MyMin 
	     << " -- pos_in_X_min " << pos_in_X_min 
	     << endl;
      }
      if(MyMin<Pvalue_enter){//the corresponding variable is included
	included=1;
	is_var_in[pos_in_X_min]=1;
      }
      gsl_vector_free(p_value_out);
    }
    else{//all are in, n_vars_out==0
      
      gsl_vector *p_value_in=gsl_vector_calloc(n_vars_in);
      for(unsigned int current_var_in=0;current_var_in<n_vars_in;current_var_in++){
	unsigned int rank_var_in=list_columns_X_gam[current_var_in];
	p_value_in->data[current_var_in]=vect_p_value->data[rank_var_in];
      }
      idx_max=gsl_vector_max_index(p_value_in);
      pos_in_X_max=list_columns_X_gam[idx_max];
      MyMax=p_value_in->data[idx_max];
      if(DEBUG){
	cout << "idx_max " << idx_max 
	     << " -- pos_in_X_max " << pos_in_X_max
	     << " -- max " << MyMax
	     << endl;
      }
      if(MyMax>Pvalue_remove&&pos_in_X_max>=nConfounders){// do not remove confounders
	removed=1;
	is_var_in[pos_in_X_max]=0;
      }     
      gsl_vector_free(p_value_in);
      
    }
  }
  if((included==0 && removed==0) || (loop>=n_loop_max-1)){
    stop=1;
  }
  return stop;
}

void store_model_per_outcome(vector < vector <unsigned int> > &Gam_step_regr,
			     vector < unsigned int > &list_columns_X_gam,
			     gsl_vector *vect_p_value,
			     gsl_vector *vect_beta_full,
			     gsl_vector *vect_SE_full,
			     Double_Matrices Gam_step_regr_pvals,
			     Double_Matrices Gam_step_regr_SE,
			     Double_Matrices Gam_step_regr_beta,
			     unsigned int outcome)
{
  for(unsigned int col=0;col<list_columns_X_gam.size();col++){
    Gam_step_regr[outcome].push_back(list_columns_X_gam[col]);
  }
  for(unsigned int col=0;col<vect_p_value->size;col++){
    Gam_step_regr_pvals.matrix[outcome][col]=vect_p_value->data[col];
    Gam_step_regr_SE.matrix[outcome][col]=vect_SE_full->data[col];
    Gam_step_regr_beta.matrix[outcome][col]=vect_beta_full->data[col];
  }
}

void getEstimateRMSE(gsl_matrix *mat_X_gam,
                     gsl_vector *current_outcome,
                     unsigned int whichOutcome,
                     gsl_vector *vect_RMSE,
                     double tolerance)
{



   ///////////////////////////////////////////
    //***************************************//
    // Step 1: calculating coefficients and  //
    //   Standard errors for variables in    //
    //         working with mat_X_gam        //
    //***************************************//
    ///////////////////////////////////////////

    //Step 1.1: Getting regression coefficients/residuals and QR decompostion
    gsl_matrix *mat_Q_red=gsl_matrix_calloc(mat_X_gam->size1,mat_X_gam->size2);
    gsl_matrix *mat_R_red=gsl_matrix_calloc(mat_X_gam->size2,mat_X_gam->size2);
    gsl_vector *vect_beta_var_in=gsl_vector_calloc(mat_X_gam->size2);
    gsl_vector *vect_residuals=gsl_vector_calloc(mat_X_gam->size1);

    get_mat_Q_red_R_red_beta_in_resid(mat_X_gam,
                                      current_outcome,
                                      vect_residuals,
                                      vect_beta_var_in,
                                      mat_Q_red,
                                      mat_R_red);

    size_t actual_rank=get_rank_from_R(mat_X_gam,
                                       tolerance);

    double dfe=mat_X_gam->size1-(double)(actual_rank);
    double RSS=gsl_stats_tss(vect_residuals->data,
                             1,
                             vect_residuals->size);
    double RMSE=sqrt(RSS/dfe);
    vect_RMSE->data[whichOutcome]=RMSE;

    gsl_vector_free(vect_beta_var_in);
    gsl_vector_free(vect_residuals);
    gsl_matrix_free(mat_Q_red);
    gsl_matrix_free(mat_R_red);


}
