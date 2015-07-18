#include "QPInfo.h"
#include "Quadrature.h"
#include "DOFVector.h"

namespace AMDiS 
{
  std::map<const Quadrature*, QPInfo*> QPInfo::qpInfos_;

  QPInfo::QPInfo(const Quadrature *quad)
    : quadrature_(quad),
      currentElInfo_(NULL),
      coordsAtQPs_(NULL),
      coordsNumPointsValid_(0),
      elementNormalAtQPs_(NULL),
      elementNormalConst_(NULL),
      elementNormalNumPointsValid_(0),
      grdLambdaAtQPs_(NULL),
      grdLambdaConst_(NULL),
      grdLambdaNumPointsValid_(0)
  {
    numPoints_ = quadrature_->getNumPoints();
  }


  QPInfo::~QPInfo()
  {
    if (coordsAtQPs_) 
      delete [] coordsAtQPs_;

    if (elementNormalAtQPs_) {
      for (int i = 0; i < numPoints_; i++)
	delete elementNormalAtQPs_[i];

      delete [] elementNormalAtQPs_;
    }

    if (elementNormalConst_)
      delete [] elementNormalConst_;

    if (grdLambdaAtQPs_) {
      for (int i = 0; i < numPoints_; i++) 
	delete grdLambdaAtQPs_[i];

      delete [] grdLambdaAtQPs_;
    }

    if (grdLambdaConst_)
      delete [] grdLambdaConst_;

    std::map<const DOFVector<double>*, VecQPInfo*>::iterator it, itEnd = vecQPInfos_.end();

    for (it = vecQPInfos_.begin(); it != itEnd; ++it) {
      delete it->second;
    }
  }


  void QPInfo::initElement(const ElInfo *elInfo)
  {
    currentElInfo_ = elInfo;
    coordsNumPointsValid_ = 0;
    elementNormalNumPointsValid_ = 0;
    grdLambdaNumPointsValid_ = 0;

    std::map<const DOFVector<double>*, VecQPInfo*>::iterator it, itEnd = vecQPInfos_.end();

    for (it = vecQPInfos_.begin(); it != itEnd; ++it) {
      it->second->valNumPointsValid_ = 0;
      it->second->grdNumPointsValid_ = 0;
      it->second->D2NumPointsValid_ = 0;
    }  
  }


  WorldVector<double> *QPInfo::getCoordsAtQPs(int numPoints) 
  {
    if (coordsNumPointsValid_ < numPoints) {
      if (!coordsAtQPs_)
	coordsAtQPs_ = new WorldVector<double>[numPoints_];
      
      for (int i = 0; i < numPoints; i++) {
	const DimVec<double>& lambda = quadrature_->getLambda(i);
	TEST_EXIT_DBG(currentElInfo_)("currentElInfo_ not set\n");
	currentElInfo_->coordToWorld(lambda, coordsAtQPs_[i]);
      }
      coordsNumPointsValid_ = numPoints;
    }
    return coordsAtQPs_;
  }


  mtl::dense_vector<double>& QPInfo::getVecAtQPs(const DOFVector<double>* vec, 
						 int numPoints,
						 const FastQuadrature *quadFast)
  {
    // check fast quadrature
    if (quadFast) {
      TEST_EXIT_DBG(quadrature_ == quadFast->getQuadrature())
	("quadrature_ != quadFast->quadrature\n");
    }
  
    // create new info if necessary
    if (vecQPInfos_[vec] == NULL)
      vecQPInfos_[vec] = new VecQPInfo;    

    VecQPInfo *localVecQPInfo = vecQPInfos_[vec];

    // update values if necessary
    if (localVecQPInfo->valNumPointsValid_ < numPoints) {

      // allocate memory if necessary
      localVecQPInfo->valAtQPs.change_dim(numPoints_);

      // fill memory
      vec->getVecAtQPs(currentElInfo_,
		       quadrature_,
		       quadFast,
		       localVecQPInfo->valAtQPs);
    
      localVecQPInfo->valNumPointsValid_ = numPoints;
    }

    // return values
    return localVecQPInfo->valAtQPs;
  }


  mtl::dense_vector<WorldVector<double> >& QPInfo::getGrdAtQPs(const DOFVector<double>* vec,
					   int numPoints,
					   const FastQuadrature *quadFast)
  {
    // check fast quadrature
    if (quadFast) {
      TEST_EXIT_DBG(quadrature_ == quadFast->getQuadrature())
	("quadrature_ != quadFast->quadrature\n");
    }
  
    // create new info if necessary
    if (vecQPInfos_[vec] == NULL)
      vecQPInfos_[vec] = new VecQPInfo;    

    VecQPInfo *localVecQPInfo = vecQPInfos_[vec];

    // update values if necessary
    if (localVecQPInfo->grdNumPointsValid_ < numPoints) {

      // allocate memory if necessary
      localVecQPInfo->grdAtQPs_.change_dim(numPoints_);

      // fill memory
      vec->getGrdAtQPs(currentElInfo_,
		       quadrature_,
		       quadFast,
		       localVecQPInfo->grdAtQPs_);

      localVecQPInfo->grdNumPointsValid_ = numPoints;
    }

    // return values
    return localVecQPInfo->grdAtQPs_;
  }


  mtl::dense_vector<WorldMatrix<double> >& QPInfo::getD2AtQPs(const DOFVector<double>* vec,
					  int numPoints,
					  const FastQuadrature *quadFast)
  {
    // check fast quadrature
    if (quadFast) {
      TEST_EXIT_DBG(quadrature_ == quadFast->getQuadrature())
	("quadrature_ != quadFast->quadrature\n");
    }
  
    // create new info if necessary
    if (vecQPInfos_[vec] == NULL)
      vecQPInfos_[vec] = new VecQPInfo;    

    VecQPInfo *localVecQPInfo = vecQPInfos_[vec];

    // update values if necessary
    if (localVecQPInfo->D2NumPointsValid_ < numPoints) {

      // allocate memory if necessary
      localVecQPInfo->D2AtQPs_.change_dim(numPoints_);

      // fill memory
      vec->getD2AtQPs(currentElInfo_,
		      quadrature_,
		      quadFast,
		      localVecQPInfo->D2AtQPs_);

      localVecQPInfo->D2NumPointsValid_ = numPoints;
    }

    // return values
    return localVecQPInfo->D2AtQPs_;
  }


  WorldVector<double> **QPInfo::getElementNormalAtQPs(int numPoints) 
  {
    TEST_EXIT_DBG(currentElInfo_)("currentElInfo_ not set\n");

    if (currentElInfo_->getParametric()) {
      if (!elementNormalAtQPs_) {
	elementNormalAtQPs_ = new WorldVector<double>*[numPoints_];
	for (int i = 0; i < numPoints_; i++)
	  elementNormalAtQPs_[i] = new WorldVector<double>;
      }

      if (elementNormalNumPointsValid_ < numPoints) {
	for (int i = 0; i < numPoints; i++) {
	  ERROR_EXIT("This does not work in the current AMDiS version!\n");
	  // Note on the following two lines of code:
	  //   There is no AMDiS version with getElementNormal with these parameters.
	  //   Maybe there were some AMDiS version before that implemented this
	  //   functionality. Must be reimplemented, if we want to make it working.
	  //      const DimVec<double>& lambda = quadrature_->getLambda(i);
	  //	  currentElInfo_->getElementNormal(*(elementNormalAtQPs_[i]), &lambda);
	}
	elementNormalNumPointsValid_ = numPoints;
      }
      return elementNormalAtQPs_;
    } else {
      if (!elementNormalConst_) {
	elementNormalConst_ = new WorldVector<double>*[numPoints_];
	elementNormalConst_[0] = new WorldVector<double>;

	for (int i = 1; i < numPoints_; i++)
	  elementNormalConst_[i] = elementNormalConst_[0];
      }
      currentElInfo_->getElementNormal(*(elementNormalConst_[0]));
      return elementNormalConst_;
    }
  }


  DimVec<WorldVector<double> > **QPInfo::getGrdLambdaAtQPs(int numPoints) 
  {
    TEST_EXIT_DBG(currentElInfo_)("currentElInfo_ not set\n");

    if (currentElInfo_->getParametric()) {
      if (!grdLambdaAtQPs_) {
	grdLambdaAtQPs_ = new DimVec<WorldVector<double> >*[numPoints_];
	for (int i = 0; i < numPoints_; i++)
	  grdLambdaAtQPs_[i] = new DimVec<WorldVector<double> >(quadrature_->getDim());
      }

      if (grdLambdaNumPointsValid_ < numPoints) {
	for (int i = 0; i < numPoints; i++) {
	  ERROR_EXIT("This does not work in the current AMDiS version!\n");
	  // Note on the following two lines of code:
	  //   There is no AMDiS version with getElementNormal with these parameters.
	  //   Maybe there were some AMDiS version before that implemented this
	  //   functionality. Must be reimplemented, if we want to make it working.
	  // 	  const DimVec<double>& lambda = quadrature_->getLambda(i);
	  // 	  currentElInfo_->calcGrdLambda(*(grdLambdaAtQPs_[i]), &lambda);
	}
	grdLambdaNumPointsValid_ = numPoints;
      }
      return grdLambdaAtQPs_;
    } else {
      if (!grdLambdaConst_) {
	grdLambdaConst_ = new DimVec<WorldVector<double> >*[numPoints_];
	grdLambdaConst_[0] = new DimVec<WorldVector<double> >(quadrature_->getDim());

	for (int i = 1; i < numPoints_; i++)
	  grdLambdaConst_[i] = grdLambdaConst_[0];
      }
      const_cast<ElInfo*>(currentElInfo_)->calcGrdLambda(*(grdLambdaConst_[0]));
      return grdLambdaConst_;
    }
  }


  QPInfo *QPInfo::provideQPInfo(const Quadrature *quad, 
				const FastQuadrature *quadFast)
  {
    if (quadFast) {
      if (quad && (quad != quadFast->getQuadrature())) {
	ERROR_EXIT("quad != quadFast->quadrature\n");
      } else {
	quad = quadFast->getQuadrature();
      }    
    }
    if (quad) {
      if (qpInfos_[quad]) 
	return qpInfos_[quad];
      QPInfo *newQPInfo = new QPInfo(quad);
      qpInfos_[quad] = newQPInfo;
      return newQPInfo;
    } else {
      return NULL;
    }
  }


  void QPInfo::clearQPInfo(const Quadrature *quad,
			   const FastQuadrature *quadFast)
  {
    if (quadFast) {
      if (quad && (quad != quadFast->getQuadrature())) {
	ERROR_EXIT("quad != quadFast->quadrature\n");
      } else {
	quad = quadFast->getQuadrature();
      }    
    }
    TEST_EXIT_DBG(quad)("no quadrature\n");
    if (qpInfos_[quad]) {
      delete qpInfos_[quad];
      qpInfos_.erase(quad);
    }
  }


  void QPInfo::clearAllQPInfos()
  {
    std::map<const Quadrature*, QPInfo*>::iterator 
      it, itEnd = qpInfos_.end();

    for (it = qpInfos_.begin(); it != itEnd; ++it) {
      delete it->second;
      qpInfos_.erase(it);
    }
  }

} // end namespace AMDiS
