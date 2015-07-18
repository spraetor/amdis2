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


#include <vector>
#include <boost/numeric/mtl/mtl.hpp>
#include "Assembler.h"
#include "SecondOrderAssembler.h"
#include "Operator.h"
#include "QPsiPhi.h"
#include "FiniteElemSpace.h"
#include "Quadrature.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS {

  ThreadPrivate<vector<SubAssembler*> > 
  SecondOrderAssembler::optimizedSubAssemblers;
  
  ThreadPrivate<vector<SubAssembler*> > 
  SecondOrderAssembler::standardSubAssemblers;


  SecondOrderAssembler::SecondOrderAssembler(Operator *op,
					     Assembler *assembler,
					     Quadrature *quad,
					     bool optimized)
    : SubAssembler(op, assembler, quad, 2, optimized)
  {}


  SecondOrderAssembler* SecondOrderAssembler::getSubAssembler(Operator* op,
							      Assembler *assembler,
							      Quadrature *quad,
							      bool optimized)
  {
    // check if a assembler is needed at all
    if (op->secondOrder.size() == 0)
      return NULL;

    SecondOrderAssembler *newAssembler;

    vector<SubAssembler*> &subAssemblers =
      optimized ? optimizedSubAssemblers.get() : standardSubAssemblers.get();

    vector<OperatorTerm*> opTerms = op->zeroOrder;

    sort(opTerms.begin(), opTerms.end());

    // check if a new assembler is needed
    for (unsigned int i = 0; i < subAssemblers.size(); i++) {
      vector<OperatorTerm*> assTerms = *(subAssemblers[i]->getTerms());
    
      sort(assTerms.begin(), assTerms.end());

      if (opTerms == assTerms && subAssemblers[i]->getQuadrature() == quad)	
	return dynamic_cast<SecondOrderAssembler*>(subAssemblers[i]);
    }

    // check if all terms are pw_const
    bool pwConst = true;
    for (unsigned int i = 0; i < op->secondOrder.size(); i++) {
      if (!op->secondOrder[i]->isPWConst()) {
	pwConst = false;
	break;
      }
    }  
   
    // create new assembler
    if (!optimized) {
      newAssembler = new Stand2(op, assembler, quad);
    } else {
      if (pwConst) {
    	newAssembler = new Pre2(op, assembler, quad);
      } else {
    	newAssembler = new Quad2(op, assembler, quad);
      }
    }

    subAssemblers.push_back(newAssembler);

    return newAssembler;
  }


  Pre2::Pre2(Operator *op, Assembler *assembler, Quadrature *quad) 
    : SecondOrderAssembler(op, assembler, quad, true)
  {
    name = "precalculated second order assembler";

    q11 = Q11PsiPhi::provideQ11PsiPhi(rowFeSpace->getBasisFcts(), 
				      colFeSpace->getBasisFcts(), 
				      quadrature);    
    LALt.resize(1);
    LALt[0].change_dim(dim + 1, dim + 1);
  }


  void Pre2::calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat)
  {
    const int **nEntries;
    const int *k, *l;
    const double *values;

    mtl::dense2D<double> &tmpMat = LALt[0];
    tmpMat = 0.0;

    for (unsigned int i = 0; i < terms.size(); i++)
      (static_cast<SecondOrderTerm*>(terms[i]))->getLALt(elInfo, LALt);    

    // Compute: 	LALt[0] *= elInfo->getDet();
    {
      for (int i = 0; i <= dim; i++)
	for (int j = 0; j <= dim; j++)
	  tmpMat[i][j] *= elInfo->getDet();
    }

    nEntries = q11->getNumberEntries();

    if (symmetric) {
      TEST_EXIT_DBG(nCol == nRow)("nCol != nRow, but symmetric assembling!\n");

      for (int i = 0; i < nRow; i++) {
	k = q11->getKVec(i, i);
	l = q11->getLVec(i, i);
	values = q11->getValVec(i, i);
	double val = 0.0;
	for (int m = 0; m < nEntries[i][i]; m++)
	  val += values[m] * tmpMat[k[m]][l[m]];	
	mat[i][i] += val;

	for (int j = i + 1; j < nCol; j++) {
	  k = q11->getKVec(i, j);
	  l = q11->getLVec(i, j);
	  values = q11->getValVec(i, j);
	  val = 0.0;
	  for (int m = 0; m < nEntries[i][j]; m++)
	    val += values[m] * tmpMat[k[m]][l[m]];	  
	  mat[i][j] += val;
	  mat[j][i] += val;
	}
      }
    } else {  /*  A not symmetric or psi != phi        */
      for (int i = 0; i < nRow; i++) {
	for (int j = 0; j < nCol; j++) {
	  k = q11->getKVec(i, j);
	  l = q11->getLVec(i, j);
	  values = q11->getValVec(i, j);
	  double val = 0.0;
	  for (int m = 0; m < nEntries[i][j]; m++)
	    val += values[m] * tmpMat[k[m]][l[m]];
	  mat[i][j] += val;
	}
      }
    }
  }


  Quad2::Quad2(Operator *op, Assembler *assembler, Quadrature *quad) 
    : SecondOrderAssembler(op, assembler, quad, true)
  {
    name = "fast quadrature second order assembler";
  }


  void Quad2::calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat)
  {
    const int nPoints = quadrature->getNumPoints();

    if (firstCall) {
      dimVec.change_dim(dim + 1);
      LALt.resize(nPoints);
      for (int j = 0; j < nPoints; j++)
	LALt[j].change_dim(dim + 1, dim + 1);

      psiFast = 
	updateFastQuadrature(psiFast, rowFeSpace->getBasisFcts(), INIT_GRD_PHI);
      phiFast = 
	updateFastQuadrature(phiFast, rowFeSpace->getBasisFcts(), INIT_GRD_PHI);
      
      firstCall = false;
    }

    for (int i = 0; i < nPoints; i++)
      LALt[i] = 0.0;

    for (unsigned int i = 0; i < terms.size(); i++)
      (static_cast<SecondOrderTerm*>(terms[i]))->getLALt(elInfo, LALt);

    if (symmetric) {
      // === Symmetric assembling. ===
      TEST_EXIT_DBG(nCol == nRow)("nCol != nRow, but symmetric assembling!\n");

      for (int iq = 0; iq < nPoints; iq++) {
	// Compute: 	LALt[iq] *= elInfo->getDet();	
	for (int i = 0; i <= dim; i++)
	  for (int j = 0; j <= dim; j++)
	    LALt[iq][i][j] *= elInfo->getDet();	

  	const vector<mtl::dense_vector<double> >& grdPsi = psiFast->getGradient(iq);
  	const vector<mtl::dense_vector<double> >& grdPhi = phiFast->getGradient(iq);
	double weight = quadrature->getWeight(iq);

	for (int i = 0; i < nCol; i++) {
	  // Compute: dimVec = quadrate->getWeight(iq) * (LALt[iq] * grdPhi[i])
	  // 	  dimVec = LALt[iq] * grdPhi[i];
	  // 	  dimVec *= quadrature->getWeight(iq);	  
	  double v = 0.0;
	  for (int j = 0; j <= dim; j++) {
	    v = 0.0;
	    for (int k = 0; k <= dim; k++)
	      v += LALt[iq][j][k] * grdPhi[i][k];
	    dimVec[j] = weight * v;
	  }	  

 	  mat[i][i] += dot(dimVec, grdPsi[i]);
 	  for (int j = i + 1; j < nRow; j++) {
 	    double tmp = dot(dimVec, grdPhi[j]);
 	    mat[i][j] += tmp;
 	    mat[j][i] += tmp;
 	  }
	}
      }
    } else {      
      // === Non symmetric assembling. ===

      for (int iq = 0; iq < nPoints; iq++) {
	// Compute: 	LALt[iq] *= elInfo->getDet();
	{
	  for (int i = 0; i <= dim; i++)
	    for (int j = 0; j <= dim; j++)
	      LALt[iq][i][j] *= elInfo->getDet();
	}

	const vector<mtl::dense_vector<double> >& grdPsi = psiFast->getGradient(iq);
	const vector<mtl::dense_vector<double> >& grdPhi = phiFast->getGradient(iq);

	for (int i = 0; i < nRow; i++) {
	  const mtl::dense_vector<double>& grdPsi_i = grdPsi[i];

	  for (int j = 0; j < nCol; j++) {
	    const mtl::dense_vector<double>& grdPhi_j = grdPhi[j];

	    // Compute: mat[i][j] += quadrature->getWeight(iq) * (grdPsi[i] * (LALt[iq] * grdPhi[j]))
	    //	    dimVec = LALt[iq] * grdPhi[j];
	    //	    mat[i][j] += quadrature->getWeight(iq) * dot(grdPsi[i], dimVec);
	    {
	      double v = 0.0;
	      for (int k = 0; k <= dim; k++) {
		double w = 0.0;
		for (int l = 0; l <= dim; l++)
		  w += LALt[iq][k][l] * grdPhi_j[l];
		v += grdPsi_i[k] * w;
	      }
	      mat[i][j] += quadrature->getWeight(iq) * v;
	    }
	  }
	}
      }
    }
  }


  Stand2::Stand2(Operator *op, Assembler *assembler, Quadrature *quad) 
    : SecondOrderAssembler(op, assembler, quad, false)
  {
    name = "standard second order assembler";
  }


  void Stand2::calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat)
  {
    mtl::dense_vector<double> grdPsi(dim + 1);
    vector<mtl::dense_vector<double> > grdPhi(nCol);
    for (int i = 0; i < nCol; i++)
      grdPhi[i].change_dim(dim + 1);

    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    int nPoints = quadrature->getNumPoints();

    vector<mtl::dense2D<double> > LALt(nPoints);
    for (int iq = 0; iq < nPoints; iq++) {
      LALt[iq].change_dim(dim + 1, dim + 1);
      LALt[iq] = 0.0;
    }

    mtl::dense_vector<double> tmpVec;

    for (unsigned int i = 0; i < terms.size(); i++)
      (static_cast<SecondOrderTerm*>(terms[i]))->getLALt(elInfo, LALt);

    if (symmetric) {
      TEST_EXIT_DBG(nCol == nRow)("nCol != nRow, but symmetric assembling!\n");

      for (int iq = 0; iq < nPoints; iq++) {
	LALt[iq] *= elInfo->getDet();

	for (int i = 0; i < nCol; i++)
	  (*(phi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPhi[i]);

	for (int i = 0; i < nRow; i++) {
	  (*(psi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPsi);

	  tmpVec = LALt[iq] * grdPhi[i];
	  mat[i][i] += quadrature->getWeight(iq) * dot(grdPsi, tmpVec);

	  for (int j = i + 1; j < nCol; j++) {
	    tmpVec = (LALt[iq] * grdPhi[j]);
	    double val = quadrature->getWeight(iq) * dot(grdPsi, tmpVec);
	    mat[i][j] += val;
	    mat[j][i] += val;
	  }
	}
      }
    } else {      /*  non symmetric assembling   */
      for (int iq = 0; iq < nPoints; iq++) {
	LALt[iq] *= elInfo->getDet();

	for (int i = 0; i < nCol; i++)
	  (*(phi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPhi[i]);

	for (int i = 0; i < nRow; i++) {
	  (*(psi->getGrdPhi(i)))(quadrature->getLambda(iq), grdPsi);
	  for (int j = 0; j < nCol; j++) {
	    tmpVec = (LALt[iq] * grdPhi[j]);	 
	    mat[i][j] += quadrature->getWeight(iq) * dot(grdPsi, tmpVec);
	  }
	}
      }
    }
  }

}
