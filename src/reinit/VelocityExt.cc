/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include "VelocityExt.h"

namespace reinit {

void 
VelocityExt::calcVelocityBoundary(DegreeOfFreedom *lokInd, const int indexV)
{
  double tempV;

  for (int nV=0; nV<nVelDOFs; ++nV) {

    tempV = 0.0;
    for(int i=0; i<=dim; i++)
      {
	tempV += (*(origVelDOF[nV]))[lokInd[i]]*lamVec[indexV][i];
      }
    (*(velDOF[nV]))[lokInd[indexV]]=tempV;
  }
}

void
VelocityExt::calcVelocity(DegreeOfFreedom *lokInd, const int indexV)
{
  double tempV;
  int tempIndex;
  if(indexFace>=0) //if there were several updates over element faces,
                   //take the coordinates that belong to the element face
                   //with the shortest distance to the interface
    {
      tempIndex=indexFace;
    }
  else //right coordinates are in lamVec[0]
    {
      tempIndex=0;
    }
  indexFace=-1;

// //   // test
//   int zero = permutation[0];
//   int one = permutation[1];
//   int two = permutation[2];
//   // int three = permutation[3];
// //   // test

  for (int nV=0; nV<nVelDOFs; ++nV) {

    tempV = 0.0;
    for(int i=0; i<=dim; i++)
      {
	tempV += (*(velDOF[nV]))[lokInd[permutation[i]]]*lamVec[tempIndex][i];
      }
    (*(velDOF[nV]))[lokInd[indexV]]=tempV;
  }
}

void
VelocityExt::setBarycentricCoords_2D_boundary(const double &l_0, const double &l_1, const double &l_2, const int index)
{
  lamVec[index][0]=l_0;
  lamVec[index][1]=l_1;
  lamVec[index][2]=l_2;
}

void
VelocityExt::setBarycentricCoords_2D(const double &l_0, const double &l_1, const double &l_2)
{
  setBarycentricCoords_2D_boundary(l_0, l_1, l_2, 0);
}

void
VelocityExt::calcBarycentricCoords_2D(const double &c_delta, const double &c_alpha, const double &norm_zhminusyh, const double &norm_xhminusyh)
{
  double s_delta=sqrt(1-c_delta*c_delta);
  double c_delta_alpha=c_delta*c_alpha+sqrt((1-c_alpha*c_alpha)*(1-c_delta*c_delta));
// double c_delta_alpha=sqrt(c_delta*c_alpha+sqrt((1-c_alpha*c_alpha)*(1-c_delta*c_delta)));
  double s_delta_alpha=sqrt(1-c_delta_alpha*c_delta_alpha);
  double c=(s_delta_alpha*norm_xhminusyh)/s_delta;
  double lambda=c/norm_zhminusyh;
  setBarycentricCoords_2D(1-lambda, lambda, 0);
}

void
VelocityExt::setBarycentricCoords_3D_boundary(const double &l_0, const double &l_1, const double &l_2, const double &l_3, const int index)
{
  lamVec[index][0]=l_0;
  lamVec[index][1]=l_1;
  lamVec[index][2]=l_2;
  lamVec[index][3]=l_3;
}

void
VelocityExt::calcBarycentricCoords_3D_boundary(const DimVec<double> sp1, const DimVec<double> sp2, const double lambda, int i)
{
  DimVec<double> tempVec(dim, NO_INIT);
  tempVec[0] = sp1[0]+lambda*(sp2[0]-sp1[0]);
  tempVec[1] = sp1[1]+lambda*(sp2[1]-sp1[1]);
  tempVec[2] = sp1[2]+lambda*(sp2[2]-sp1[2]);
  tempVec[3] = sp1[3]+lambda*(sp2[3]-sp1[3]);
  setBarycentricCoords_3D_boundary(tempVec[0], tempVec[1], tempVec[2], tempVec[3], i);
}

void
VelocityExt::copyAndExpandFaceCoords_3D(int vertNum,int index)
{
  double t0, t1, t2;
  t0=lamVec[0][0];
  t1=lamVec[0][1];
  t2=lamVec[0][2];
  switch(vertNum)
    {
    case 0: //coordinates belong to element face 0
      setBarycentricCoords_3D_boundary(0,t0,t1,t2,index);
      break;

    case 1: //coordinates belong to element face 1
      setBarycentricCoords_3D_boundary(t0,0,t1,t2,index);
      break;
      
    case 2: //coordinates belong to element face 2
      setBarycentricCoords_3D_boundary(t0,t1,0,t2,index);
      break;

    default:
      break;
    }

}

void
VelocityExt::setBarycentricCoords_3D(const double &l_0, const double &l_1, const double &l_2, const double &l_3)
{
    indexFace = -1;
  setBarycentricCoords_3D_boundary(l_0, l_1, l_2, l_3, 0);
}

void
VelocityExt::setIndexFaceWithShortestDist(int indexFace_)
{
  indexFace=indexFace_;
}

void
VelocityExt::setPermutation(int vertNum, int mTrav)
{
  switch(dim)
    {
    case 2:
      switch(mTrav)
	{
	case 0: //level
	  setPermutation(vertNum);
	  if(vertNum==1)
	    {
	      swapVertices(0,1);
	    }
	  break;

	case 1: //traverse
	  setPermutation(vertNum);
	  break;

	default:
	  break;
	}
      break;

    case 3:
      switch(mTrav)
	{
	case 1: //traverse
	  setPermutation(vertNum);
	  break;

	default:
	  break;
	}
      break;

    default:
      break;
    }
}

void
VelocityExt::setPermutation(int vertNum)
{
  switch(dim)
    {
    case 2:
      switch(vertNum)
	{
	case 0: //update is for vertex 0
	  setPermutation_2D(1,2,0);
	  break;
	  
	case 1: //update is for vertex 1
	  setPermutation_2D(2,0,1);
	  break;
	  
	case 2: //update is for vertex 2
	  setPermutation_2D(0,1,2);
	  break;
	  
	default:
	  break;
	}
      break;

    case 3:
      switch(vertNum)
	{
	case 0: //update is for vertex 0
	  setPermutation_3D(1,2,3,0);
	  break;
	  
	case 1: //update is for vertex 1
	  setPermutation_3D(2,3,0,1);
	  break;
	  
	case 2: //update is for vertex 2
	  setPermutation_3D(3,0,1,2);
	  break;

	case 3: //update is for vertex 3
	  setPermutation_3D(0,1,2,3);
	  break;
	  
	default:
	  break;
	}
      break;

    default:
      break;
    }
}

void
VelocityExt::setPermutation_2D(int i_0, int i_1, int i_2)
{
  permutation[0]=i_0;
  permutation[1]=i_1;
  permutation[2]=i_2;
}


void
VelocityExt::setPermutation_3D(int i0, int i1, int i2, int i3)
{
  permutation[0]=i0;
  permutation[1]=i1;
  permutation[2]=i2;
  permutation[3]=i3;
}

void
VelocityExt::swapVertices(int i1, int i2)
{
  int temp=permutation[i1];
  permutation[i1]=permutation[i2];
  permutation[i2]=temp;
}

}
