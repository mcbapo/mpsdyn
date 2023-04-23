/** 
 *\file splitLocalOp.h
 * One often encounters the situation that one has a small k-local Operator \f$ O= A_1\otimes A_2 \otimes \dots \otimes A_k \f$ and one wants to split it into MPO form via SVD. This function takes the mwArray containing the operator \f$ O \f$, a vector containing the  physical dimensions \f$ d_1,\dots,d_k \f$ and a vector which after the function finished contains the MPO matrices.
 * 
 *
 \author Stefan Kühn
 \date 08.07.2015*/

#ifndef SPLIT_LOCAL_OP_H
#define SPLIT_LOCAL_OP_H

#include <vector>
#include <numeric>

#include "wrapper.h"
#include "mwArray.h"
#include "Operator.h"

using namespace shrt;

inline void splitLocalOp(mwArray Op,vector<int> dims,vector<mwArray> &matrices)
{
  //Delete potential content of the vector
  matrices.clear();
  
  //Matrices for the SVDs
  mwArray U,S,V;  
  
  ////////////////////////////////////////////////
  //Some error checking
  ////////////////////////////////////////////////
  
  //Is the given Operator a matrix?
  if(!Op.isMatrix())
  {
    cout << "Error in splitLocalOp(), the given operator is not a matrix" << endl;
    exit(666);
  }
  
  //Check if the product of the given physical dimensions is compliant with the given mwArray
  int proddim=accumulate (dims.begin(),dims.end(),1,multiplies<int>());
  Indices matrix_dim=Op.getDimensions();
  if(matrix_dim!=Indices(proddim,proddim))
  {
    cout << "Error in splitLocalOp(), the given physical dimensions are not compliant with the given operator" << endl;
    cout << "Dimension of matrix " << matrix_dim << ", product of dimensions " << proddim << endl;
    exit(666);
  }
  
  ////////////////////////////////////////////////
  //Check for special cases
  ////////////////////////////////////////////////
  if(dims.size()==1)
  {
    //Single site operator given, all I have to do is to reshape it accordingly to get an object with bond dimension one
    cout << "Single site operator given, will only reshape it" << endl;
    matrices.push_back(reshape(Op,Indices(dims[0],1,dims[0],1)));
    return;
  }
  if(dims.size()==2)
  {
    //I only need a single SVD
    cout << "Two site operator given, only computing a single SVD" << endl;
    Op.reshape(Indices(dims[0],dims[1],dims[0],dims[1]));
    Op.permute(Indices(1,3,2,4));
    Op.reshape(Indices(dims[0]*dims[0],dims[1]*dims[1]));
    wrapper::svd(Op,U,S,V);  
    //Redistribute the singular values
    S=sqrt(S);
    U.multiplyRight(S);
    V.multiplyLeft(S);
    //Reshape the matrix and save results
    U.reshape(Indices(dims[0],1,dims[0],-1));
    matrices.push_back(U);
    V.reshape(Indices(-1,dims[1],dims[1],1));
    V.permute(Indices(2,1,3,4));
    matrices.push_back(V);  
    return;
  }  
  ////////////////////////////////////////////////
  //First permute the indices in the right order  
  ////////////////////////////////////////////////  
  //Dimensions of the input matrix
  vector <int> indices;
  //As an operator has two physical legs and the kronecker product in the code orders the indices as described above, I simply double the physical dimensions
  indices = dims;	indices.insert(indices.end(), dims.begin(), dims.end());  
  //Reshape
  Op.reshape(Indices(indices));  
  //Now compute the order of the indices
  indices.clear();
  vector<int> opdims;
  for(int i=1; i<=dims.size(); i++)
  {
    indices.push_back(i);
    indices.push_back(i+dims.size());
    opdims.push_back(dims[i-1]);
    opdims.push_back(dims[i-1]);
  }  
  //Permute the indices accordingly, such that two physical legs of a single site operator are successive
  Op.permute(Indices(indices));
  
  ////////////////////////////////////////////////
  //Compute the SVDs  
  ////////////////////////////////////////////////   
  //Now take the first two indices together and do a SVD  
  Op.reshape(Indices(accumulate(opdims.begin(),opdims.begin()+2,1,multiplies<int>()),accumulate(opdims.begin()+2,opdims.end(),1,multiplies<int>())));
  wrapper::svd(Op,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(opdims[0],1,opdims[1],-1));
  matrices.push_back(U);
  Op=V;  
  //Proceed with the SVDs up to the last one
  int count = 2;
  Indices dim;
  for (int i=2; i<(dims.size()-1); i++)
  {
    dim=Op.getDimensions();
    Op.reshape(Indices(dim[0]*accumulate(opdims.begin()+count,opdims.begin()+count+2,1,multiplies<int>()),accumulate(opdims.begin()+count+2,opdims.end(),1,multiplies<int>())));
    wrapper::svd(Op,U,S,V);  
    //Redistribute the singular values
    S=sqrt(S);
    U.multiplyRight(S);
    V.multiplyLeft(S);
    //Reshape the matrix and save results
    U.reshape(Indices(dim[0],opdims[count],opdims[count+1],-1));
    U.permute(Indices(2,1,3,4));
    matrices.push_back(U);
    Op=V;
    count+=2;
  }  
  //Last SVD
  dim=Op.getDimensions();
  Op.reshape(Indices(dim[0]*accumulate(opdims.begin()+count,opdims.begin()+count+2,1,multiplies<int>()),accumulate(opdims.begin()+count+2,opdims.end(),1,multiplies<int>())));
  wrapper::svd(Op,U,S,V);  
  //Redistribute the singular values
  S=sqrt(S);
  U.multiplyRight(S);
  V.multiplyLeft(S);
  //Reshape the matrix and save results
  U.reshape(Indices(dim[0],opdims[count],opdims[count+1],-1));
  U.permute(Indices(2,1,3,4));
  matrices.push_back(U);
  V.reshape(Indices(-1,opdims[count+2],opdims[count+3],1));
  V.permute(Indices(2,1,3,4));
  matrices.push_back(V);  
  cout << "Successfully computed SVDs" << endl;
}

#endif //SPLIT_LOCAL_OP_H