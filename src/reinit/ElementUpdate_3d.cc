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


#include "ElementUpdate_3d.h"

namespace reinit {

using namespace AMDiS;

double 
ElementUpdate_3d::calcElementUpdate(
		  const FixVec<WorldVector<double> *, VERTEX> &vert_,
		  FixVec<double, VERTEX> &uhVal)
{
  FUNCNAME("ElementUpdate_3d::calcElementUpdate");
  
  int dim = 3;
  double update = 0.0;
  WorldVector<double> B1;
  WorldVector<double> C1;
  WorldVector<double> D1;
  double uhVal0_orig;       // original value of vertex 0 after sort and
                            // translation of vertices but before
                            // translation of function values to zero
  FixVec<WorldVector<double> *, VERTEX> vert(dim, NO_INIT);
  for (int i=0; i<=dim; ++i) {
    vert[i] = vert_[i];
  } 
  
  // ===== Sort vertices and translate element to point of origin. =====
  sortAndTranslateElement(vert, uhVal, uhVal0_orig);
  
  // ===== Rotate element. 
  //       The rotated element is stored in B1, C1, D1 (A = zero vector
  //       after sort and translation and is fix under rotation).
  //       Return value: true - rotation exists
  //                     false - rotation does not exist
  bool rotationFlag = rotateElement(vert, uhVal, B1, C1, D1);
  
  // ===== rotation exists ==> calculate projection and update value
  //                           depending on the position of the
  //                           projection point
  //       rotation does not exist ==> calculate update values on the
  //                                   faces of element
 
   if (rotationFlag) {
   
     // Calculate projection of B1, C1 and D1 on x-y-plane -> B2, C2, D2.
    WorldVector<double> B2;
    WorldVector<double> C2;
    WorldVector<double> D2;
    
    B2[0] = B1[0]; B2[1] = B1[1]; B2[2] = 0.0;
    C2[0] = C1[0]; C2[1] = C1[1]; C2[2] = 0.0;
    D2[0] = D1[0]; D2[1] = D1[1]; D2[2] = 0.0;
    
    // Calculate position of D2 with respect to A, B2 and C2. 
    // If D2 lies in the triangle AB2C2, this is the minimizing point 
    // of the update formula. Else, the minimizing point lies on edges 
    // of this triangle and the 2-dimensional update formula is applied
    // on faces of element to get the minimizing point.
    int posFlag = calcPosition(B2, C2, D2);
    
    WorldVector<double> tmpVec;
    switch(posFlag) {
      
    case IN_ABC: // update value is the distance of D1 to the x-y-plane
      update = D1[2];
      break;
      
    case VERT_A: // minimizing point is vertex A
      update = uhVal[0] + sqrt(*(vert[3]) * *(vert[3]));
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->setBarycentricCoords_3D(1,0,0,0);
	}
      break;
      
    case VERT_B: // minimizing point is vertex B
      tmpVec = (*(vert[3])) - (*(vert[1]));
      update = uhVal[1] + sqrt(tmpVec * tmpVec);
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->setBarycentricCoords_3D(0,1,0,0);
	}
      break;
      
    case VERT_C: // minimizing point is vertex C
      tmpVec = (*(vert[3])) - (*(vert[2]));
      update = uhVal[2] + sqrt(tmpVec * tmpVec);
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->setBarycentricCoords_3D(0,0,1,0);
	}
      break;
      
    case EDGE_AB: // minimizing point lies on edge AB
      update = calcFaceUpdate(vert[0], vert[1], vert[3], 
			      uhVal[0], uhVal[1]);
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->copyAndExpandFaceCoords_3D(2,0);
	}
      break;
      
    case EDGE_AC: // minimizing point lies on edge AC
      update = calcFaceUpdate(vert[0], vert[2], vert[3], 
			      uhVal[0], uhVal[2]);
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->copyAndExpandFaceCoords_3D(1,0);
	}
      break;
      
    case EDGE_BC: // minimizing point lies on edge BC
      update = calcFaceUpdate(vert[1], vert[2], vert[3], 
			      uhVal[1], uhVal[2]);
      //set barycentric coordinats for calculation of the velocity
      if(velExt != NULL)
	{
	  velExt->copyAndExpandFaceCoords_3D(0,0);
	}
      break;
      
    default: ERROR_EXIT("illegal position flag !\n");
      break; 
    }
    
  }
  else {  // rotation does not exist
    
    // ===== Calculate the 2-dimensional update for D on all faces 
    //       adjacent with D. The update value is the minimum of these
    //       updates. =====
    double tmp_update;
    
    update = calcFaceUpdate(vert[0], vert[1], vert[3],
			    uhVal[0], uhVal[1]);
    //set barycentric coordinats for calculation of the velocity
    //save index of element face
    if(velExt != NULL)
      {
	velExt->copyAndExpandFaceCoords_3D(2,1);
	velExt->setIndexFaceWithShortestDist(1);
      }
    tmp_update = calcFaceUpdate(vert[0], vert[2], vert[3],
				uhVal[0], uhVal[2]);
    //set barycentric coordinats for calculation of the velocity
    if(velExt != NULL)
      {
	velExt->copyAndExpandFaceCoords_3D(1,2);
      }
    if (tmp_update < update) { 
      update = tmp_update;
      //save index of element face, if diatance is shorter
      if(velExt != NULL)
	{
	  velExt->setIndexFaceWithShortestDist(2);
	}
    }
    tmp_update = calcFaceUpdate(vert[1], vert[2], vert[3],
				uhVal[1], uhVal[2]);
    //set barycentric coordinats for calculation of the velocity
    if(velExt != NULL)
      {
	velExt->copyAndExpandFaceCoords_3D(0,3);
      }
    if (tmp_update < update) { 
      update = tmp_update;
      //save index of element face, if distance is shorter
      if(velExt != NULL)
	{
	  velExt->setIndexFaceWithShortestDist(3);
	}
    }
  }
  
  // ===== Correct initial translation of function values to zero. =====
  update += uhVal0_orig;
  
  return update;
}

void 
ElementUpdate_3d::sortAndTranslateElement(
                  FixVec<WorldVector<double> *, VERTEX> &vert,
		  FixVec<double, VERTEX> &uhVal,
		  double &uhVal0_orig)
{
  WorldVector<double> *tmpVec;
  double tmpVal;
  
  // ===== Sort vertices. =====
  if (uhVal[0] > uhVal[1]) {
    
    if (uhVal[1] > uhVal[2]) {
      
      tmpVec = vert[0];
      vert[0] = vert[2];
      vert[2] = tmpVec;
      
      tmpVal = uhVal[0];
      uhVal[0] = uhVal[2];
      uhVal[2] = tmpVal;
      //save permutation of vertices
      if(velExt != NULL)
	{
	  velExt->swapVertices(0,2);
	}
    }
    else if (uhVal[0] > uhVal[2]) {
      
      tmpVec = vert[0];
      vert[0] = vert[1];
      vert[1] = vert[2];
      vert[2] = tmpVec;
      
      tmpVal = uhVal[0];
      uhVal[0] = uhVal[1];
      uhVal[1] = uhVal[2];
      uhVal[2] = tmpVal;
      //save permutation of vertices
      if(velExt != NULL)
	{
	  velExt->swapVertices(0,1);
	  velExt->swapVertices(2,1);
	}
    }
    else {
      
      tmpVec = vert[0];
      vert[0] = vert[1];
      vert[1] = tmpVec;
      
      tmpVal = uhVal[0];
      uhVal[0] = uhVal[1];
      uhVal[1] = tmpVal;
      //save permutation of vertices
      if(velExt != NULL)
	{
	  velExt->swapVertices(0,1);
	}
    }
  }
  else if (uhVal[1] > uhVal[2]) {
    
    if (uhVal[0] > uhVal[2]) {
      
      tmpVec = vert[0];
      vert[0] = vert[2];
      vert[2] = vert[1];
      vert[1] = tmpVec;
      
      tmpVal = uhVal[0];
      uhVal[0] = uhVal[2];
      uhVal[2] = uhVal[1];
      uhVal[1] = tmpVal;
      //save permutation of vertices
      if(velExt != NULL)
	{
	  velExt->swapVertices(0,2);
	  velExt->swapVertices(1,2);
	}
    }
    else {
      
      tmpVec = vert[1];
      vert[1] = vert[2];
      vert[2] = tmpVec;
      
      tmpVal = uhVal[1];
      uhVal[1] = uhVal[2];
      uhVal[2] = tmpVal;
      //save permutation of vertices
      if(velExt != NULL)
	{
	  velExt->swapVertices(1,2);
	}
    }
  }    
  
  // ===== Translation: A is moved to point of origin. =====
  *(vert[1]) -= *(vert[0]);
  *(vert[2]) -= *(vert[0]);
  *(vert[3]) -= *(vert[0]);
  *(vert[0]) = 0.0;
  
  // ===== Adapt \ref uhVal. =====
  uhVal[1] -= uhVal[0];
  uhVal[2] -= uhVal[0];
  uhVal0_orig = uhVal[0];
  uhVal[0] = 0.0;
}

bool 
ElementUpdate_3d::rotateElement(
		  const FixVec<WorldVector<double> *, VERTEX> &vert,
		  const FixVec<double, VERTEX> &uhVal,
		  WorldVector<double> &B1,
		  WorldVector<double> &C1,
		  WorldVector<double> &D1) 
{
  int orient = 1;              // orientation of AB, AC and AD
  WorldVector<double> B = *(vert[1]);
  WorldVector<double> C = *(vert[2]);
  WorldVector<double> D = *(vert[3]);
  WorldVector<double> N;       // unit normal for ABC 
                               // (orientation of AB, AC and N is 1)
  WorldVector<double> P;       // projection of D on ABC
  WorldVector<double> unit_B;  // normalized AB
  double norm_B;               // norm of vector AB
  double norm_C;               // norm of vector AC
  double norm_P;               // norm of vector AP
  double dD;                   // (signed) distance of D to ABC
  double lambda, nue;          // scalars
  WorldVector<double> N1;      // unit normal for A1B1C1
  WorldVector<double> unit_B1; // normalized B1
  double dB1 = uhVal[1];       // distance of B1 to x-y-plane
  double dC1 = uhVal[2];       // distance of C1 to x-y-plane
  WorldVector<double> E_21;    // vector of rotated coordinate system
  int dim = 3;
  WorldVector<double> tmpVec;
  double tmp;


  // ===== Calculate orientation of ABCD. =====
  
  // Get unit normal for ABC.
  tmpVec[0] = B[1]*C[2] - B[2]*C[1];
  tmpVec[1] = B[2]*C[0] - B[0]*C[2];
  tmpVec[2] = B[0]*C[1] - B[1]*C[0];
  tmp = sqrt(tmpVec * tmpVec);
  for (int i=0; i<dim; ++i) { N[i] = tmpVec[i] / tmp; }
  
  // Calculate distance of D to ABC.
  dD = D * N;
  
  // Orientation.
  orient = (dD < 0) ? -1 : 1;
  
  
  // ===== Calculate B1. 
  //       We construct B1 in x-z-plane with distance dB1 to 
  //       x-y-plane. =====
  
  norm_B = sqrt(B * B);
  
  tmp = norm_B*norm_B - dB1*dB1;
  if (tmp < 1.e-15) {
    return false;   // rotation is not possible
  }
  B1[0] = sqrt(tmp);
  B1[1] = 0.0;
  B1[2] = dB1;
  
  
  // ===== Calculate C1.
  //       If (E_1, E_2, E_3) is a coordinate system with 
  //       AC = lambda*E_1 + nue*E_2, then for the rotated coordinate
  //       system (E_11, E_21, E_31) we have 
  //       A1C1 = lambda*E_11 + nue*E_21.
  //
  //       The coordinate system (E_1, E_2, E_3) does not need to 
  //       be calculated, but we assume
  //          --> E_1 = unit_B,
  //          --> nue > 0. =====
  
  // Calculate lambda and nue.
  for (int i=0; i<dim; ++i) { unit_B[i] = B[i] / norm_B;  }
  norm_C = sqrt(C * C);
  
  lambda = C * unit_B;      // corresponds to norm(AC)*cos(AC,AB)
  nue =  sqrt(norm_C*norm_C - lambda*lambda);
                            // corresponds to norm(AC)*sin(AC,AB)
  
  // Calculate E_11 (= unit_B1).
  for (int i=0; i<dim; ++i) { unit_B1[i] = B1[i] / norm_B; }
  
  // Calculate E_21.
  E_21[2] = (dC1 - lambda*unit_B1[2])/nue; 
                            // A1C1 = lambda*E_11 + nue*E_21 and C1[2] = dC1 

  E_21[0] = -unit_B1[2]*E_21[2]/unit_B1[0];  
                            // E_11*E_21 = 0 and unit_B[1]=0
  
  tmp = 1 - E_21[0]*E_21[0] - E_21[2]*E_21[2];
  if (tmp < 1.e-15) {
    return false;    // rotation is not possible
  }
  E_21[1] = orient * sqrt(tmp);  // norm(E_21) = 1 and
                                      // in further construction: D1 
                                      // above x-y-plane

  // Calculate C1.
  C1[0] = lambda*unit_B1[0] + nue*E_21[0];
  C1[1] = lambda*unit_B1[1] + nue*E_21[1];
  C1[2] = lambda*unit_B1[2] + nue*E_21[2];
  
  
  // ===== Calculate D1. 
  //       If we have AP = lambda*E_1 + nue*E_2 for the projection P of 
  //       D on ABC , then D1 = lambda*E_11 + nue*E_21 + dD*N1. =====
  
  // Calculate P.
  P[0] = D[0] - dD*N[0];
  P[1] = D[1] - dD*N[1];
  P[2] = D[2] - dD*N[2];
  
  // Calculate lambda and nue.
  norm_P = sqrt(P * P);
  lambda = P * unit_B;      // corresponds to norm(AP)*cos(AP,AB)
  nue =  sqrt(norm_P*norm_P - lambda*lambda);
             // corresponds to norm(AP)*sin(AP,AB)
             // Fix problem for angle(E_2,P) > Pi/2:
             //   project P onto E_2-axes (-> tmpVec), 
             //   if angle(tmpVec,AC) > Pi/2 the scalar nue has negative sign.
  tmpVec = P - unit_B*lambda;
  if (tmpVec * C < 0) { nue *= -1; }
  
  // Calculate N1.
  tmpVec[0] = B1[1]*C1[2] - B1[2]*C1[1];
  tmpVec[1] = B1[2]*C1[0] - B1[0]*C1[2];
  tmpVec[2] = B1[0]*C1[1] - B1[1]*C1[0];
  tmp = sqrt(tmpVec * tmpVec);
  for (int i=0; i<dim; ++i) { N1[i] = tmpVec[i] / tmp; }
  
  // Calculate D1.
  D1[0] = lambda*unit_B1[0] + nue*E_21[0] + dD*N1[0];
  D1[1] = lambda*unit_B1[1] + nue*E_21[1] + dD*N1[1];
  D1[2] = lambda*unit_B1[2] + nue*E_21[2] + dD*N1[2];

//   // Check projection of triangle A1B1C1 (-> triangle A2B2C2) to x-y-plane.
//   // Rotation is said to be not possible if
//   //   -> angles are too small: | cos(A2B2,A2C2) | < 1-1.e-3 
//   //   -> area of A2B2C2 is too small: 
//   //         | det(A2B2 A2C2) | > max(1.e-10, 0.01*norm_B*norm_B)
//   double abs_det_B2C2 = fabs(B1[0]*C1[1]-B1[1]*C1[0]);
//   double limit = (0.01*norm_B*norm_B < 1.e-10) ? 
//     1.e-10 : 0.01*norm_B*norm_B;
//   if (abs_det_B2C2 < limit)
//     return false;

//   double norm_AB2 = sqrt(B1[0]*B1[0]+B1[1]*B1[1]);
//   double norm_AC2 = sqrt(C1[0]*C1[0]+C1[1]*C1[1]);
//   double cos_B2AC2 = (B1[0]*C1[0]+B1[1]*C1[1]) / (norm_AB2*norm_AC2);
//   if (fabs(cos_B2AC2) > 0.999)
//     return false;
  
  return true;
}

int 
ElementUpdate_3d::calcPosition(const WorldVector<double> &B2,
			       const WorldVector<double> &C2,
			       const WorldVector<double> &D2)
{
  FUNCNAME("ElementUpdate_3d::calcPosition");
  
  double l_A, l_B2, l_C2;
  double tmp_B2, tmp_C2;
  double norm_AB2 = sqrt(B2 * B2);
  double norm_AC2 = sqrt(C2 * C2);
  double sP_AB2AC2 = B2 * C2;
  double cos_B2AC2 = sP_AB2AC2 / (norm_AB2*norm_AC2);
  int posFlag = -1;
  
  // ===== Calculate barycentric coordinates of D2 with respect to 
  //       triangle AB2C2. =====
  if (fabs(cos_B2AC2) < 1.e-15) {
    double sP_AB2AD2 = B2 * D2;
    double sP_AC2AD2 = C2 * D2;
    
    tmp_B2 = (sP_AB2AD2 < 0) ? -D2[0] : D2[0]; 
    tmp_C2 = (sP_AC2AD2 < 0) ? -D2[1] : D2[1]; 

    if (C2[1] < 0)
      tmp_C2 *= -1;
  }
  else {
    double sin_B2AC2 = sqrt(1-cos_B2AC2*cos_B2AC2);

    if (C2[1] < 0)
      sin_B2AC2 *= -1;
    
    tmp_C2 = D2[1]/sin_B2AC2;
    tmp_B2 = D2[0] - D2[1]*cos_B2AC2/sin_B2AC2;
  }
  
  l_B2 = tmp_B2 / norm_AB2;
  l_C2 = tmp_C2 / norm_AC2;

//   // for test purposes: barycentric coordinates with Cramer's rule
//   double det_B2C2 = B2[0]*C2[1]-B2[1]*C2[0];

//   TEST_EXIT(fabs(det_B2C2) > 1.e-15)("illegal projection triangle !\n");

//   double det_D2C2 = D2[0]*C2[1]-D2[1]*C2[0];
//   double det_B2D2 = B2[0]*D2[1]-B2[1]*D2[0];

//   double l_B2_Cramer = 1.0*det_D2C2/det_B2C2;
//   double l_C2_Cramer = 1.0*det_B2D2/det_B2C2;

// //   TEST_EXIT(fabs(l_B2-l_B2_Cramer) < 1.e-8 && 
// // 	    fabs(l_C2-l_C2_Cramer) < 1.e-8)
//   TEST_EXIT(fabs(l_B2*l_B2_Cramer) >= 0 && 
// 	    fabs(l_C2*l_C2_Cramer) >= 0)
//     ("problems in barycentric coordinates of projection point !\n");
//   // end: for test purposes
  
  /*     // Corrections. */
  /*     if (fabs(l_B2) < 1.e-15) l_B2 = 0.0; */
  /*     if (fabs(1-l_B2) < 1.e-15) l_B2 = 1.0; */
  /*     if (fabs(l_C2) < 1.e-15) l_C2 = 0.0; */
  /*     if (fabs(1-l_C2) < 1.e-15) l_C2 = 1.0; */
  
  l_A = 1-l_B2-l_C2;
 
  //save barycentric coordinates for calculation of the velocity
  if(velExt != NULL)
    {
      velExt->setBarycentricCoords_3D(l_A,l_B2,l_C2,0);
    }
  
  // ===== Determine position of D2. =====
  if (l_A < 0) {
    if (l_B2 < 0) {
      if (l_C2 < 0) {
	ERROR_EXIT("calculation error: sum of barycentric coordinates is equal to 1 !\n");
      }
      else {
	posFlag = VERT_C;
      }
    }      
    else if (l_C2 < 0) {
      posFlag = VERT_B;
    }
    else {  // EDGE_BC
      posFlag = EDGE_BC;
      /* 	if (l_B2 > 0) { */
      /* 	  posFlag = (l_C2 > 0) ? EDGE_BC : VERT_B; */
      /* 	} */
      /* 	else { */
      /* 	  if (l_C2 > 0) { */
      /* 	    posFlag = VERT_C; */
      /* 	  } */
      /* 	  else { */
      /* 	    ERROR_EXIT("no barycentric coordinates !\n"); */
      /* 	  } */
      /* 	} */
    }
  }
  else if (l_B2 < 0) {
    if (l_C2 < 0) {
      posFlag = VERT_A;
    }
    else {  // EDGE_AC 
      posFlag = EDGE_AC;
      /* 	if (l_A > 0) { */
      /* 	  posFlag = (l_C2 > 0) ? EDGE_AC : VERT_A; */
      /* 	} */
      /* 	else { */
      /* 	  if (l_C2 > 0) { */
      /* 	    posFlag = VERT_C; */
      /* 	  } */
      /* 	  else { */
      /* 	    ERROR_EXIT("no barycentric coordinates !\n"); */
      /* 	  } */
      /* 	} */
    }
  }      
  else if (l_C2 < 0) {  // EDGE_AB
    posFlag = EDGE_AB;
    /*       if (l_A > 0) { */
    /* 	posFlag = (l_B2 > 0) ? EDGE_AB : VERT_A; */
    /*       } */
    /*       else { */
    /* 	if (l_B2 > 0) { */
    /* 	  posFlag = VERT_B; */
    /* 	} */
    /* 	else { */
    /* 	  ERROR_EXIT("no barycentric coordinates !\n"); */
    /* 	} */
    /*       } */
  }
  else {
    posFlag = IN_ABC;
  }
  
  return posFlag;
}

double 
ElementUpdate_3d::calcFaceUpdate(WorldVector<double> *A2d,
				 WorldVector<double> *B2d,
				 WorldVector<double> *C2d,
				 double &uhValA2d,
				 double &uhValB2d)
{
  int dim = A2d->getSize();
  FixVec<WorldVector<double> *, VERTEX> vert2d(dim, NO_INIT);
  FixVec<double, VERTEX> uhVal2d(dim, NO_INIT);

  vert2d[0] = A2d;
  vert2d[1] = B2d;
  vert2d[2] = C2d;
  
  uhVal2d[0] = uhValA2d;
  uhVal2d[1] = uhValB2d;
  
  return elUpdate2d->calcElementUpdate(vert2d, uhVal2d);
}

}
