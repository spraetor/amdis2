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


namespace AMDiS {

  /****************************************************************************/
  /*  standard routines for getting values of a finite element function and   */
  /*  its first (and second) derivatives at the quadrature points             */
  /****************************************************************************/

  template<typename T> 
  const T *uhAtQp(const FastQuadrature *quadFast,
		  const T *uhLoc, T *vec)
  {
    FUNCNAME("uhAtQp()");
    static T *quadVec = NULL;
    static int size = 0;
    T *val;
    const double *phi;
    int i, j, k;

    int nPoints = quadFast->getQuadrature()->getNumPoints();
    int nBasFcts = quadFast->getBasisFunctions()->getNumber();
    int dim = quadFast->getDim();

    int vecSize = uhLoc[0].size();

    if (vec) {
      val = vec;
    } else {
      if (size < nPoints)  {
	int newSize = std::max(Quadrature::maxNQuadPoints[dim], nPoints);
	if(quadVec) delete [] quadVec;
	quadVec = new T[size](vecSize);
	size = newSize;
      }
      val = quadVec;
    }

    for (i = 0; i < nPoints; i++) {
      phi = quadFast->getPhi(i);
      for(k = 0; k < vecSize; k++) {
	for (val[i][k] = j = 0; j < nBasFcts; j++) {
	  val[i][k] += uhLoc[j][k] * phi[j];
	}
      }
    }
    return(const_cast<const T*>(val));
  }

  template<typename T, typename GrdT> 
  const GrdT *grdUhAtQp(const FastQuadrature *quadFast,
			const DimVec<WorldVector<double> >& grdLambda,
			const T *uhLoc, GrdT *vec)
  {
    FUNCNAME("grdUhAtQp()");

    int i, j, k, l, m;

    int dim = quadFast->getDim();
    int dow = Global::getGeo(WORLD);
    int parts = Global::getGeo(PARTS, dim);

    int nPoints = quadFast->getQuadrature()->getNumPoints();
    int nBasFcts = quadFast->getBasisFunctions()->getNumber();

    static GrdT *quadVec = NULL;
    static int size = 0;
    GrdT *val;
    VectorOfFixVecs<DimVec<double> > *gradPhi;
    DimVec<double> grd1(dim, NO_INIT);

    int vecSize = uhLoc[0].size();

    if (vec) {
      val = vec;
    } else {
      if(size < nPoints) {
	int newSize = std::max(Quadrature::maxNQuadPoints[dim], nPoints);
	if(quadVec) delete [] quadVec; 
	quadVec = new GrdT[newSize](vecSize);
	size = newSize;
      }
      val = quadVec;
    }

    for (i = 0; i < nPoints; i++) {
      gradPhi = quadFast->getGradient(i);
      for(m = 0; m < vecSize; m++) {
	grd1 = 0.0;
	for (j = 0; j < nBasFcts; j++) {
	  for (k = 0; k < parts; k++) {
	    grd1[k] += (*gradPhi)[j][k] * uhLoc[j][m];
	  }
	}

	for(l=0; l < dow; l++) {
	  for (val[i][m][l] = k = 0; k < parts; k++) {
	    val[i][m][l] += grdLambda[k][l] * grd1[k];
	  }
	}
      }
    }

    return(const_cast<const GrdT*>(val));
  }

  template<typename T, typename D2T> 
  const D2T* D2UhAtQp(const FastQuadrature *quadFast,
		      const DimVec<WorldVector<double> >& grdLambda,
		      const T *uhLoc, D2T *vec)
  {
    FUNCNAME("D2UhAtQp()");

    int i, j, k, l, m, iq;

    int nPoints = quadFast->getQuadrature()->getNumPoints();
    int nBasFcts = quadFast->getBasisFunctions()->getNumber();
    int dim = quadFast->getDim();

    int parts = Global::getGeo(PARTS, dim);
    int dow   = Global::getGeo(WORLD);

    static D2T *quadVec = NULL;
    static int size = 0;
    D2T *val;
    const VectorOfFixVecs<DimMat<double> > *D2Phil;
    DimMat<double>  D2Tmp(dim, DEFAULT_VALUE, 0.);

    int vecSize = uhLoc[0].size();

    if (vec) {
      val = vec;
    } else {
      if(size < nPoints)  {
	int newSize = std::max(Quadrature::maxNQuadPoints[dim], nPoints);
	if(quadVec) delete [] quadVec;
	quadVec = new D2T[newSize](vecSize);
	size = newSize;
      }
      val = quadVec;
    }

    for (iq = 0; iq < nPoints; iq++) {
      for(m = 0; m < vecSize; m++) {
	D2Tmp = 0.0;
	//for (k = 0; k < parts; k++)
	//  for (l = 0; l < parts; l++)
	//    D2Tmp[k][l] = 0.0;

	D2Phil = quadFast->getSecDer(iq);
	for (i = 0; i < nBasFcts; i++) {
	  for (k = 0; k < parts; k++)
	    for (l = 0; l < parts; l++)
	      D2Tmp[k][l] += uhLoc[i][m]*(*D2Phil)[i][k][l];
	}

	for (i = 0; i < dow; i++)
	  for (j = 0; j < dow; j++) {
	    val[iq][m][i][j] = 0.0;
	    for (k = 0; k < parts; k++)
	      for (l = 0; l < parts; l++)
		val[iq][m][i][j] += grdLambda[k][i]*grdLambda[l][j]*D2Tmp[k][l];
	  }
      }
    }

    return(const_cast<const D2T*>(val));
  }

}
