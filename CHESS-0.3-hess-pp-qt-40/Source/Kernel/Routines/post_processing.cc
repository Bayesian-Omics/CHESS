/* This file is part of CHESS.
 *      Copyright (c) Habib Saadi (h.saadi@imperial.ac.uk)
 *      2013
 *
 * The file is modified from post_processing.cc in the ESS++ program
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

#include "post_processing.h"

#define DEBUG 0

//To remove range check in gsl
//#define GSL_RANGE_CHECK_OFF


using namespace std;

int getUniqueList(vector<vector<unsigned int> > &uniqueListModels,
                        const vector<vector<unsigned int> >& listModels,
                        const vector<unsigned int>& nModelsVisited,
                        const unsigned int& burnIn,
                        const unsigned int& pX,
                        const unsigned int& nConfounders,
                        bool Include_Single_Models)
{
    int posNullModel=-1;

    unsigned int nModels=listModels.size();
    uniqueListModels.resize(nModels+pX-nConfounders+1);// Add null models and all size 1 models

    vector <unsigned int> alreadyVisited;
    vector <unsigned int> singleEffectsPresent;
    singleEffectsPresent.assign(pX-nConfounders,0);

    alreadyVisited.assign(nModels,0);
    vector<vector<unsigned int> >::const_iterator itCur,itCmp;
    unsigned int j=0,k=0,jj=0;

    for(itCur=listModels.begin();itCur<listModels.end();++itCur)
    {
        // Set the model we wish to compare
        if(alreadyVisited[j]==0)
        {
            // Mark the model as having been visited
            alreadyVisited[j]=1;
            //Only the models that have not been matched are considered
            unsigned int refModelSize=itCur->size()-1;
            if(refModelSize==nConfounders)
            {
                posNullModel=k;
            }
            k++;
            if(refModelSize==nConfounders+1)
            {
                unsigned int tmpPosn = (*itCur)[nConfounders+1]-nConfounders;
                singleEffectsPresent[tmpPosn] = 1;
            }
            unsigned int nRepeats = 1;

            // Define a working vector
            vector<unsigned int> workVector;

            // Work out how many times this model was visited in the output
            unsigned int r = j;
            for(itCmp=itCur;itCmp<listModels.end();++itCmp)
            {
                if(alreadyVisited[r]==0)
                {
                  if((*itCmp)==(*itCur))
                  {
                    alreadyVisited[r]=1;
                    nRepeats++;
                  }
                }
                r++;
            }

            // Populate the temporary vector
            workVector.push_back(nRepeats);
            workVector.push_back(j+1);
            workVector.push_back(nModelsVisited[j+burnIn+1]);
            workVector.push_back(refModelSize);
            workVector.insert(workVector.begin()+4,itCur->begin()+1,itCur->end());

            // Add the working vector to the output
            uniqueListModels[jj].insert(uniqueListModels[jj].begin(),workVector.begin(),workVector.end());
            jj++;
        }
        j++;
    }

    if(posNullModel<0)
    {
        // Add the null model to the list of unique models

        // Define a working vector
        vector<unsigned int> workVector;
        workVector.assign(4,0);
        workVector[3]=nConfounders;
        for(unsigned int i=0;i<nConfounders;i++)
        {
            workVector.push_back(i);
        }
        // Push this working vector back to the output
        //uniqueListModels.push_back(workVector);
        uniqueListModels[jj].insert(uniqueListModels[jj].begin(),workVector.begin(),workVector.end());
        jj++;
        posNullModel=k;
    }
    if (Include_Single_Models==1)
    {
      for(unsigned int i=0;i<pX-nConfounders;i++)
      {
          // Add the single effects models that weren't visited
          if(singleEffectsPresent[i]==0){
              vector<unsigned int> workVector;
              workVector.assign(4,0);
              workVector[3]=1+nConfounders;
              for(unsigned int kk=0;kk<nConfounders;kk++)
              {
                  workVector.push_back(kk);
              }
              workVector.push_back(i+nConfounders);
              //uniqueListModels.push_back(workVector);
              uniqueListModels[jj].insert(uniqueListModels[jj].begin(),workVector.begin(),workVector.end());
              jj++;
          }
      }
    }
    uniqueListModels.resize(jj);
    return(posNullModel);

}

void getLogPost(vector<vector<double> >& margLogPostVec,
                  vector<vector<unsigned int> >& uniqueListModels,
                  gsl_matrix *matX,
                  gsl_matrix *matY,
                  Prior_param PR,
                  bool gPriorFlag,
                  bool indepPriorFlag,
                  bool gSampleFlag,
                  double lambda,
                  double gMean,
                  bool cudaFlag,
                  Model_Information &Model_Data)
{
  unsigned int nUniqueModels=0;
  unsigned int margLogVecSize;

  nUniqueModels=uniqueListModels.size();
  margLogVecSize=nUniqueModels;
  margLogPostVec.resize(nUniqueModels);


  gsl_vector *propLogMargCondPost=gsl_vector_calloc(2);
  unsigned int pX=matX->size2;

  //Step3: Calculating the models' posterior:
  //We plug-in E(g)
  vector<vector<double> > logPost(margLogVecSize);

  for(unsigned int model=0;model<margLogVecSize;model++){
    logPost[model].resize(2);
    unsigned int nVarsIn,offset;

    nVarsIn=uniqueListModels[model][3];
    offset=4;

//    gsl_matrix *matXGam;
//    if(nVarsIn>0)
//    {
//      vector<unsigned int> listColumnsXGam;
//      for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++){
//        listColumnsXGam.push_back(uniqueListModels[model][currVar]);
//      }
//      matXGam=get_X_reduced(listColumnsXGam,matX);
//    }

    //Build Gamma
    vector <unsigned int> Gamma;
    Gamma.resize(pX);
    if(nVarsIn>0)
    {
        for(unsigned int currVar=offset;currVar<nVarsIn+offset;currVar++)
        {
            Gamma[uniqueListModels[model][currVar]]=1;
        }
    }

    double propLogLik,propLogPost;
//  computeLogPosterior_PostProc(propLogLik,propLogPost,matXGam,matY,PR,gPriorFlag,
    unsigned int Chain=0;
    unsigned int Response=0;
    computeLogPosterior( Gamma,Chain,Response,
                         matX,propLogLik,propLogPost,matY,PR,gPriorFlag,
                         indepPriorFlag,gSampleFlag,lambda,gMean,
                         //pX,nX,pY,
                         cudaFlag,
                         Model_Data);
    propLogMargCondPost->data[0]=propLogLik;
    propLogMargCondPost->data[1]=propLogPost;

//    if(nVarsIn>0){
//      gsl_matrix_free(matXGam);
//    }

    logPost[model][0]=propLogMargCondPost->data[0];
    logPost[model][1]=propLogMargCondPost->data[1];
  }

  for(unsigned int i=0;i<nUniqueModels;i++){
    margLogPostVec[i].assign(logPost[i].begin(),logPost[i].end());
  }

  gsl_vector_free(propLogMargCondPost);

}

void getAndSortPostGam(gsl_vector* postGamVec,
                              gsl_permutation* idxPostGamSort,
                              const vector<vector<double> >& margLogPostVec)
{

  //Initialize post_gam
  unsigned int nUnique=margLogPostVec.size();

  double meanVal=0.0;
  for(unsigned int model=0;model<nUnique;model++){
    postGamVec->data[model]=margLogPostVec[model][1];
    meanVal+=postGamVec->data[model];
  }
  meanVal=meanVal/(double)nUnique;

  gsl_sort_vector_index(idxPostGamSort,postGamVec);
  gsl_permutation_reverse(idxPostGamSort);

  double cumSum=0.0;
  for(unsigned int model=0;model<nUnique;model++){
    cumSum+=exp(postGamVec->data[model]-meanVal);
  }
 
  for(unsigned int model=0;model<nUnique;model++){
    postGamVec->data[model]=exp(postGamVec->data[model]-meanVal)/cumSum;
  }
  

  
}

void combineAndPrintBestModel(ofstream& fOut,
				  gsl_permutation* const idxPostGamSort,
				  gsl_vector* const postGamVec,
				  const vector<vector<double> >& margLogPostVec,
				  const vector<vector<unsigned int> >& uniqueListModels,
				  const unsigned int& nRetainedModels,
				  const unsigned int& posNullModel,
				  const bool& fullOutput,
				  const unsigned int& nConfounders)
{
  double nullMargLik=margLogPostVec[posNullModel][0];

  for(unsigned int model=0;model<nRetainedModels;model++){
    unsigned int origPosn=(unsigned int)(idxPostGamSort->data[model]);
    fOut << model+1 << "\t"
	  << uniqueListModels[origPosn][0] << "\t";//#repeats
    if(fullOutput){
	  fOut << uniqueListModels[origPosn][1] << "\t"//sweep first visit
	      << uniqueListModels[origPosn][2] << "\t";//#models visited before 1st visit
    }
    fOut << uniqueListModels[origPosn][3]-nConfounders << "\t"//#vars in
	  << setprecision(16) 
	  << margLogPostVec[origPosn][1] << "\t";//logPost
    fOut << setprecision(16)
	  << postGamVec->data[origPosn] << "\t";//postGam;

    fOut << setprecision(16)
	  << log10(exp(margLogPostVec[origPosn][0]-nullMargLik)) << "\t";//JS;
    for(unsigned int varIn=nConfounders;varIn<uniqueListModels[origPosn][3];varIn++){
      fOut << uniqueListModels[origPosn][varIn+4]-nConfounders+1 << "\t";
    }
    fOut << endl;
  }
}

void getAndPrintMargGam(ofstream& fOut,
                const vector<vector<unsigned int> >& uniqueListModels,
                gsl_vector* const postGamVec,
                const unsigned int& pX,
                const unsigned int& nConfounders)
{
    // New version

    vector <double> margGamVec;
    vector <double> normConstVec;
    margGamVec.assign(pX-nConfounders,0.0);
    normConstVec.assign(pX-nConfounders,0.0);
    unsigned int nUnique=uniqueListModels.size();

    for(unsigned int model=0;model<nUnique;model++){

        unsigned int nRepeats=uniqueListModels[model][0];
        double currWeight=nRepeats*postGamVec->data[model];

        const vector<unsigned int> *Current_Model=&uniqueListModels[model];
        unsigned int Nb_Vars_in=(unsigned int) Current_Model->size()-4;

        if(uniqueListModels[model][3]>nConfounders)
        {
            for(unsigned int j=0;j<Nb_Vars_in;j++)
            {
                unsigned int currVar=uniqueListModels[model][j+4];
                if(currVar>=nConfounders)
                {
                    if(DEBUG)
                    {
                        cout << setw(10) << j << setw(10) << currVar << setw(30) << endl;
                    }
                    margGamVec[currVar-nConfounders]+=currWeight;
                }
            }
        }
        for(unsigned int currVar=nConfounders;currVar<pX;currVar++)
        {
            normConstVec[currVar-nConfounders]+=currWeight;
        }
    }
    for(unsigned int currVar=0;currVar<pX-nConfounders;currVar++){
      fOut << currVar+1 << "\t"
            << setprecision(15) << margGamVec[currVar]/normConstVec[currVar] << endl;
    }
}
