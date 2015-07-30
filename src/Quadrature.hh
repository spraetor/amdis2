namespace AMDiS 
{
  /****************************************************************************/
  /*  standard routines for getting values of a finite element function and   */
  /*  its first (and second) derivatives at the quadrature points             */
  /****************************************************************************/

  template <class T> 
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

  
  template <class T, class GrdT> 
  const GrdT *grdUhAtQp(const FastQuadrature *quadFast,
			const DimVec<WorldVector<double> >& grdLambda,
			const T *uhLoc, GrdT *vec)
  {
    FUNCNAME("grdUhAtQp()");

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

    if (vec)
      val = vec;
    else {
      if (size < nPoints) {
      	int newSize = std::max(Quadrature::maxNQuadPoints[dim], nPoints);
      	if(quadVec) delete [] quadVec; 
      	quadVec = new GrdT[newSize](vecSize);
      	size = newSize;
      }
      val = quadVec;
    }

    for (int i = 0; i < nPoints; i++) {
      gradPhi = quadFast->getGradient(i);
      for (int m = 0; m < vecSize; m++) {
      	grd1 = 0.0;
      	for (int j = 0; j < nBasFcts; j++) {
      	  for (int k = 0; k < parts; k++) {
      	    grd1[k] += (*gradPhi)[j][k] * uhLoc[j][m];
      	  }
      	}
      
      	for (int l = 0; l < dow; l++) {
          int k;
      	  for (val[i][m][l] = k = 0; k < parts; k++) {
      	    val[i][m][l] += grdLambda[k][l] * grd1[k];
      	  }
      	}
      }
    }

    return(const_cast<const GrdT*>(val));
  }

  
  template <class T, class D2T> 
  const D2T* D2UhAtQp(const FastQuadrature *quadFast,
            		      const DimVec<WorldVector<double> >& grdLambda,
            		      const T *uhLoc, D2T *vec)
  {
    FUNCNAME("D2UhAtQp()");

    int nPoints = quadFast->getQuadrature()->getNumPoints();
    int nBasFcts = quadFast->getBasisFunctions()->getNumber();
    int dim = quadFast->getDim();

    int parts = Global::getGeo(PARTS, dim);
    int dow   = Global::getGeo(WORLD);

    static D2T *quadVec = NULL;
    static int size = 0;
    D2T *val;
    const VectorOfFixVecs<DimMat<double> > *D2Phil;
    DimMat<double>  D2Tmp(dim, dim, 0.0);

    int vecSize = uhLoc[0].size();

    if (vec)
      val = vec;
    else {
      if (size < nPoints)  {
      	int newSize = std::max(Quadrature::maxNQuadPoints[dim], nPoints);
      	if (quadVec) delete [] quadVec;
      	quadVec = new D2T[newSize](vecSize);
      	size = newSize;
      }
      val = quadVec;
    }

    for (int iq = 0; iq < nPoints; iq++) {
      for (int m = 0; m < vecSize; m++) {
      	D2Tmp = 0.0;
      	//for (k = 0; k < parts; k++)
      	//  for (l = 0; l < parts; l++)
      	//    D2Tmp[k][l] = 0.0;
      
      	D2Phil = quadFast->getSecDer(iq);
      	for (int i = 0; i < nBasFcts; i++) {
      	  for (int k = 0; k < parts; k++)
      	    for (int l = 0; l < parts; l++)
      	      D2Tmp[k][l] += uhLoc[i][m] * (*D2Phil)[i][k][l];
      	}
      
      	for (int i = 0; i < dow; i++)
      	  for (int j = 0; j < dow; j++) {
      	    val[iq][m][i][j] = 0.0;
      	    for (int k = 0; k < parts; k++)
      	      for (int l = 0; l < parts; l++)
      		      val[iq][m][i][j] += grdLambda[k][i] * grdLambda[l][j] * D2Tmp[k][l];
      	  }
      }
    }

    return(const_cast<const D2T*>(val));
  }

} // end namespace AMDiS
