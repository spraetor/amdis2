#include <list>
#include <algorithm>
#include <math.h>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/utility/category.hpp>
#include <boost/numeric/linear_algebra/identity.hpp>

#include "FixVec.h"
#include "Boundary.h"
#include "DOFAdmin.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Global.h"
#include "Mesh.h"
#include "Quadrature.h"
#include "BoundaryManager.h"
#include "Assembler.h"
#include "Operator.h"
#include "Initfile.h"
#include "Traverse.h"
#include "DualTraverse.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MpiHelper.h"
#include "parallel/MeshDistributor.h"
#endif

// Defining the interface for MTL4
namespace mtl 
{
  // Let MTL4 know that DOFVector it is a column vector
  namespace traits 
  {
    template <class T>
    struct category< AMDiS::DOFVector<T> > 
    {
      typedef tag::dense_col_vector type;
    };
  } // end namepsace traits

  
  namespace ashape 
  {
    template <class T>
    struct ashape< AMDiS::DOFVector<T> > 
    {
      typedef cvec<typename ashape<T>::type> type;
    };
  } // end namespace ashape

  // Modelling Collection and MutableCollection
  template <class T>
  struct Collection< AMDiS::DOFVector<T> > 
  {
    typedef T               value_type;
    typedef const T&        const_reference;
    typedef std::size_t     size_type;
  };

  template <typename T>
  struct MutableCollection< AMDiS::DOFVector<T> > 
    : public Collection< AMDiS::DOFVector<T> > 
  {
    typedef T&              reference;
  };

} // namespace mtl



namespace AMDiS 
{
  template <class T>
  DOFVectorBase<T>::DOFVectorBase(const FiniteElemSpace *f, std::string n)
    : feSpace(f),
      name(n),
      elementVector(f->getBasisFcts()->getNumber()),
      boundaryManager(NULL)
  {    
    nBasFcts = feSpace->getBasisFcts()->getNumber();
    dim = feSpace->getMesh()->getDim();
  }
  

  template <class T>
  DOFVectorBase<T>::~DOFVectorBase()
  {}

 
  template <class T>
  DOFVector<T>::DOFVector(const FiniteElemSpace* f, std::string n, bool addToSynch)
    : DOFVectorBase<T>(f, n)
  {
    vec.resize(0);
    init(f, n, addToSynch);
  } 


  template <class T>
  void DOFVector<T>::init(const FiniteElemSpace* f, std::string n, bool addToSynch)
  {
    this->name = n;
    this->feSpace = f;
    if (this->feSpace && this->feSpace->getAdmin())
      (this->feSpace->getAdmin())->addDOFIndexed(this);
    this->boundaryManager = new BoundaryManager(f);
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    if (addToSynch && Parallel::MeshDistributor::globalMeshDistributor != NULL) 
      Parallel::MeshDistributor::globalMeshDistributor->addInterchangeVector(this);
#endif
  }

  template <class T>
  DOFVector<T>::~DOFVector()
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    if ( Parallel::MeshDistributor::globalMeshDistributor != NULL) 
      Parallel::MeshDistributor::globalMeshDistributor->removeInterchangeVector(this);
#endif
    if (this->feSpace && this->feSpace->getAdmin())
      (this->feSpace->getAdmin())->removeDOFIndexed(this);
    
    if (this->boundaryManager)
      delete this->boundaryManager;
    
    vec.clear();
  }

  template <class T>
  void DOFVectorBase<T>::addElementVector(T factor, 
					  const ElementVector &elVec, 
					  const BoundaryType *bound,
					  ElInfo *elInfo,
					  bool add)
  {
    std::vector<DegreeOfFreedom> indices(nBasFcts);
    feSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(), 
					     feSpace->getAdmin(),
					     indices);

    for (int i = 0; i < nBasFcts; i++) {
      BoundaryCondition *condition = 
	bound ? this->getBoundaryManager()->getBoundaryCondition(bound[i]) : NULL;

      if (!(condition && condition->isDirichlet())) {
	DegreeOfFreedom irow = indices[i];

	if (add)
	  (*this)[irow] += factor * elVec[i];
	else  
	  (*this)[irow] = factor * elVec[i];
      }
    }
  }


  template <class T>
  double DOFVector<T>::squareNrm2() const
  {
    checkFeSpace(this->feSpace, vec);
    
    double nrm = 0.0;
    Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(this)), 
			 USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator)
      nrm += (*vecIterator) * (*vecIterator);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(nrm);
#endif

    return nrm;
  }

  template <class T>
  double DOFVector<T>::nrm2() const
  {
    return std::sqrt(squareNrm2());
  }


  template <class T>
  T DOFVector<T>::asum() const
  {
    checkFeSpace(this->feSpace, vec);

    double nrm = 0.0;
    Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(this)), 
			 USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator)
      nrm += std::abs(*vecIterator);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(nrm);
#endif

    return nrm;
  }


  template <class T>
  T DOFVector<T>::sum() const
  {
    checkFeSpace(this->feSpace, vec);

    double nrm = 0.0;    
    Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(this)), 
			 USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator)
      nrm += *vecIterator;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(nrm);
#endif

    return nrm;
  }


  template <class T>
  void DOFVector<T>::set(T alpha)
  {
    checkFeSpace(this->feSpace, vec);

    Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(this), USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator)
      *vecIterator = alpha ; 
  }


  template <class T>
  void DOFVector<T>::copy(const DOFVector<T>& x)
  {
    FUNCNAME_DBG("DOFVector<T>::copy()");

    checkFeSpace(this->feSpace, vec);

    TEST_EXIT_DBG(static_cast<int>(x.vec.size()) >= 
		  this->feSpace->getAdmin()->getUsedSize())
      ("x.size = %d too small: admin->sizeUsed = %d\n", x.vec.size(), 
       this->feSpace->getAdmin()->getUsedSize());
    
    Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(this), USED_DOFS);
    Iterator xIterator(dynamic_cast<DOFVector<T>*>(const_cast<DOFVector<T>*>(&x)), 
		       USED_DOFS);
    for (vecIterator.reset(), xIterator.reset(); !vecIterator.end();
	 ++vecIterator, ++xIterator)      
      *vecIterator = *xIterator; 
  }


  template <class T>
  T DOFVector<T>::min() const
  {    
    checkFeSpace(this->feSpace, vec);

    T m;    
    Iterator vecIterator(const_cast<DOFIndexed<T>*>(dynamic_cast<const DOFIndexed<T>*>(this)), USED_DOFS);
    for (vecIterator.reset(), m = *vecIterator; !vecIterator.end(); ++vecIterator)
      m = std::min(m, *vecIterator);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMin(m);
#endif

    return m;
  }


  template <class T>
  T DOFVector<T>::max() const 
  {    
    checkFeSpace(this->feSpace, vec);

    T m;    
    Iterator vecIterator(const_cast<DOFIndexed<T>*>(dynamic_cast<const DOFIndexed<T>*>(this)), USED_DOFS);
    for (vecIterator.reset(), m = *vecIterator; !vecIterator.end(); ++vecIterator)
      m = std::max(m, *vecIterator);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMax(m);
#endif

    return m;
  }


  template <class T>
  T DOFVector<T>::absMax() const
  {
    return std::max(std::abs(max()), std::abs(min()));
  }


  template <class T>
  T DOFVector<T>::average() const
  {    
    checkFeSpace(this->feSpace, vec);

    int count = 0;
    T m = 0;
    Iterator vecIterator(const_cast<DOFIndexed<T>*>(dynamic_cast<const DOFIndexed<T>*>(this)), USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator) {
      m += *vecIterator;
      count++;
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(m);
    Parallel::mpi::globalAdd(count);
#endif

    return m / count;
  }


  template <class T>
  size_t DOFVector<T>::calcMemoryUsage() const
  {
    size_t result = 0;
    result += sizeof(DOFVector<T>);
    result += vec.size() * sizeof(T);

    return result;
  }


  template <class T>
  T DOFVectorBase<T>::evalUh(const DimVec<double>& lambda, 
			     DegreeOfFreedom* dof_indices)
  {
    BasisFunction* phi = const_cast<BasisFunction*>(this->getFeSpace()->getBasisFcts());
    int nBasisFcts = phi->getNumber();
    T val = 0.0;

    for (int i = 0; i < nBasisFcts; i++)
      val += (*this)[dof_indices[i]]*(*phi->getPhi(i))(lambda);

    // TODO: ist das im Parallelen so richtig???
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double localVal = val;
    MPI::COMM_WORLD.Allreduce(&localVal, &val, 1, MPI_DOUBLE, MPI_SUM);
#endif

    return val;
  }


  template <class T>
  double DOFVector<T>::Int(int meshLevel, Quadrature* q) const
  {
    Mesh* mesh = this->feSpace->getMesh();

    if (!q) {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(this->dim, deg);
    }

    FastQuadrature *quadFast = 
      FastQuadrature::provideFastQuadrature(this->feSpace->getBasisFcts(), *q, INIT_PHI);

    Flag flag = Mesh::FILL_COORDS | Mesh::FILL_DET;
    if (meshLevel == -1)
      flag |= Mesh::CALL_LEAF_EL;
    else
      flag |= Mesh::CALL_EL_LEVEL;

    T result; nullify(result);
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<T> uh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo =  stack.traverseFirst(mesh, meshLevel, flag);
    while (elInfo) {
      double det = elInfo->getDet();
      T normT; nullify(normT);
      this->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      for (int iq = 0; iq < nPoints; iq++)
	normT += quadFast->getWeight(iq) * (uh_vec[iq]);
      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(result);
#endif
    
    return result;
  }
  
  
  template <class T>
  double DOFVector<T>::L1Norm(Quadrature* q) const
  {  
    Mesh* mesh = this->feSpace->getMesh();

    if (!q) {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(this->dim, deg);
    }

    FastQuadrature *quadFast = 
      FastQuadrature::provideFastQuadrature(this->feSpace->getBasisFcts(), *q, INIT_PHI);

    double result = 0.0;
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<T> uh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET);
    while (elInfo) {
      double det = elInfo->getDet();
      double normT = 0.0;
      this->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      for (int iq = 0; iq < nPoints; iq++)
	normT += quadFast->getWeight(iq) * abs(uh_vec[iq]);
      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double localResult = result;
    MPI::COMM_WORLD.Allreduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM);
#endif
    
    return result;
  }


  template <class T>
  double DOFVector<T>::L2NormSquare(Quadrature* q) const
  {
    using math::sqr;
    Mesh* mesh = this->feSpace->getMesh();

    if (!q) {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(this->dim, deg);
    }

    FastQuadrature *quadFast = 
      FastQuadrature::provideFastQuadrature(this->feSpace->getBasisFcts(), *q, INIT_PHI);

    double result = 0.0;
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<T> uh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET);
    while (elInfo) {
      double det = elInfo->getDet();
      double normT = 0.0;
      this->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      for (int iq = 0; iq < nPoints; iq++)
	normT += quadFast->getWeight(iq) * sqr(uh_vec[iq]);
      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double localResult = result;
    MPI::COMM_WORLD.Allreduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM);
#endif
    
    return result;
  }


  template <class T>  
  double DOFVector<T>::H1NormSquare(Quadrature *q) const
  {
    using math::sqr;
    Mesh* mesh = this->feSpace->getMesh();

    if (!q) {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree() - 2;
      q = Quadrature::provideQuadrature(this->dim, deg);
    }

    FastQuadrature *quadFast = 
      FastQuadrature::provideFastQuadrature(this->feSpace->getBasisFcts(), *q, INIT_GRD_PHI);

    double result = 0.0;
    int nPoints = quadFast->getNumPoints();
    int dimOfWorld = Global::getGeo(WORLD);
    mtl::dense_vector<WorldVector<T> > grduh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | 
			  Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA);
    while (elInfo) {
      double det = elInfo->getDet();
      double normT = 0.0;
      this->getGrdAtQPs(elInfo, NULL, quadFast, grduh_vec);

      for (int iq = 0; iq < nPoints; iq++) {
	double norm2 = 0.0;
	for (int j = 0; j < dimOfWorld; j++)
	  norm2 += sqr(grduh_vec[iq][j]);
	normT += quadFast->getWeight(iq) * norm2;
      }

      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double localResult = result;
    MPI::COMM_WORLD.Allreduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM);
#endif
    
    return result;
  }



  template <class T>
  bool DOFVector<T>::getDofIdxAtPoint(WorldVector<double> &p, 
				      DegreeOfFreedom &idx, 
				      ElInfo *oldElInfo, 
				      bool useOldElInfo) const
  { 
    Mesh *mesh = this->feSpace->getMesh();
    const BasisFunction *basFcts = this->feSpace->getBasisFcts();

    int dim = mesh->getDim();
    int numBasFcts = basFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DimVec<double> lambda(dim, NO_INIT);

    ElInfo *elInfo = mesh->createNewElInfo();
    idx = 0;
    bool inside = false;

    if (oldElInfo && useOldElInfo && oldElInfo->getMacroElement()) {
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, oldElInfo->getMacroElement(), NULL, NULL);
      delete oldElInfo;
    } else {
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, NULL, NULL, NULL);
    }
    
    if (oldElInfo)
      oldElInfo = elInfo;
    
    if (!inside){
      delete elInfo; 
      return false;
    }
    
    basFcts->getLocalIndices(elInfo->getElement(), this->feSpace->getAdmin(), localIndices);
    
    WorldVector<double> coord;
    int minIdx = -1;
    double minDist = 1.e15;
    
    for (int i = 0; i < numBasFcts; i++) {
      elInfo->coordToWorld(*(basFcts->getCoords(i)), coord);
      WorldVector<double> dist = coord - p;
      double newDist = norm(dist);
      if (newDist < minDist) {
	minIdx = i;
	minDist = newDist;
	if (minDist < 1.e-10)
	  break;
      }
    }
    
    if (minIdx >= 0)
      idx = localIndices[minIdx];
    
    if(!oldElInfo) delete elInfo;
    return inside;
  }


  template <class T>
  void DOFVector<T>::compressDOFIndexed(int first, int last, 
					std::vector<DegreeOfFreedom> &newDOF)
  {
    for (int i = first; i <= last; i++)
      if (newDOF[i] >= 0)
	vec[newDOF[i]] = vec[i];
  }


  template <class T>
  Flag DOFVectorBase<T>::getAssembleFlag()
  {
    Flag fillFlag(0);
    
    for (std::vector<Operator*>::iterator op = this->operators.begin(); 
	 op != this->operators.end(); ++op)
      fillFlag |= (*op)->getFillFlag();

    return fillFlag;
  }


  template <class T>
  void DOFVectorBase<T>::finishAssembling()
  {
    // call the operatos cleanup procedures
    for (std::vector<Operator*>::iterator it = this->operators.begin();
	 it != this->operators.end(); ++it)
      (*it)->finishAssembling();
  }


  template <class T>
  DOFVector<T>& DOFVector<T>::operator=(const DOFVector<T>& rhs)
  {
    this->feSpace = rhs.feSpace;
    this->dim = rhs.dim;
    this->nBasFcts = rhs.nBasFcts;
    vec = rhs.vec;
    this->elementVector.change_dim(this->nBasFcts);
    this->operators = rhs.operators;
    this->operatorFactor = rhs.operatorFactor;
    
    if (rhs.boundaryManager) {
      if (this->boundaryManager) 
	delete this->boundaryManager; 

      this->boundaryManager = new BoundaryManager(*rhs.boundaryManager);
    } else {
      this->boundaryManager = NULL;
    }

    return *this;
  }


  template <class T>
  void DOFVectorBase<T>::getLocalVector(const Element *el, 
					mtl::dense_vector<T>& d) const
  {
    FUNCNAME_DBG("DOFVectorBase<T>::getLocalVector()");

    TEST_EXIT_DBG(feSpace->getMesh() == el->getMesh())
      ("Element is defined on a different mesh than the DOF vector!\n");

    std::vector<DegreeOfFreedom> localIndices(nBasFcts);
    const DOFAdmin* admin = feSpace->getAdmin();
    feSpace->getBasisFcts()->getLocalIndices(el, admin, localIndices);

    for (int i = 0; i < nBasFcts; i++)
      d[i] = (*this)[localIndices[i]];
  }


  template <class T>
  void DOFVectorBase<T>::getVecAtQPs(const ElInfo *elInfo, 
				     const Quadrature *quad,
				     const FastQuadrature *quadFast,
				     mtl::dense_vector<T>& vecAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getVecAtQPs()");
    
    TEST_EXIT_DBG(quad || quadFast)
      ("Neither quad nor quadFast defined!\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
      ("quad != quadFast->quadrature\n");    
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
      ("Invalid basis functions!");

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    mtl::dense_vector<T> localVec(nBasFcts);
    getLocalVector(elInfo->getElement(), localVec);
    
    if (quadFast) {
      // using precalculated values at QPs
      const mtl::dense2D<double>& phi = quadFast->getPhi();
      vecAtQPs.change_dim(num_rows(phi)); 	// = quadrature->getNumPoints()
      vecAtQPs = phi * localVec; 		// Matrix<double> * Vector<T>
    } else {
      // evaluate basisFunctions at QPs
      int nPoints = quad->getNumPoints();
      vecAtQPs.change_dim(nPoints);

      for (int iq = 0; iq < nPoints; iq++) {
       	nullify(vecAtQPs[iq]);
       	for (int j = 0; j < nBasFcts; j++)
       	  vecAtQPs[iq] += 
       	    localVec[j] * (*(basFcts->getPhi(j)))(quad->getLambda(iq));
      }
    }
  }


  template <class T>
  void DOFVectorBase<T>::getGrdAtQPs(const ElInfo *elInfo,
				     const Quadrature *quad,
				     const FastQuadrature *quadFast,
				     mtl::dense_vector<typename GradientType<T>::type> &grdAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGrdAtQPs()");

    TEST_EXIT_DBG(quad || quadFast)("neither quad nor quadFast defined\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
      ("quad != quadFast->quadrature\n");
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
      ("invalid basis functions");

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    int dow = Global::getGeo(WORLD);
    int nPoints = quadFast ? quadFast->getQuadrature()->getNumPoints() : quad->getNumPoints();

    mtl::dense_vector<T> localVec(nBasFcts);
    this->getLocalVector(elInfo->getElement(), localVec);

    mtl::dense_vector<T> grd1(dim + 1);
    int parts = Global::getGeo(PARTS, dim);
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();

    grdAtQPs.change_dim(nPoints);
    if (quadFast) {
      for (int iq = 0; iq < nPoints; iq++) {
	nullify(grd1);

	for (int j = 0; j < nBasFcts; j++) { // #BasisFunctions 
	  for (int k = 0; k < parts; k++)  // #edges (2d) or #faces (3d)
	    grd1[k] += quadFast->getGradient(iq, j, k) * localVec[j];
	}

	for (int l = 0; l < dow; l++) {
	  nullify(grdAtQPs[iq][l]);
	  for (int k = 0; k < parts; k++)
	    grdAtQPs[iq][l] += grdLambda[k][l] * grd1[k];
	}
      }

    } else {
      mtl::dense_vector<double> grdPhi(dim + 1);

      for (int iq = 0; iq < nPoints; iq++) {
	nullify(grd1);

	for (int j = 0; j < nBasFcts; j++) {
	  (*(basFcts->getGrdPhi(j)))(quad->getLambda(iq), grdPhi);
	  for (int k = 0; k < parts; k++)
	    grd1[k] += grdPhi[k] * localVec[j];
	}

	for (int l = 0; l < dow; l++) {
	  nullify(grdAtQPs[iq][l]);
	  for (int k = 0; k < parts; k++)
	    grdAtQPs[iq][l] += grdLambda[k][l] * grd1[k];
	}
      }
    }
  }


  template <class T>
  void DOFVectorBase<T>::getDerivativeAtQPs(const ElInfo *elInfo,
					    const Quadrature *quad,
					    const FastQuadrature *quadFast,
					    int comp,
					    mtl::dense_vector<T> &derivativeAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGrdAtQPs()");

    TEST_EXIT_DBG(quad || quadFast)("neither quad nor quadFast defined\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
      ("quad != quadFast->quadrature\n");
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
      ("invalid basis functions");

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    int nPoints = quadFast ? quadFast->getQuadrature()->getNumPoints() : quad->getNumPoints();

    mtl::dense_vector<T> localVec(nBasFcts);
    this->getLocalVector(elInfo->getElement(), localVec);

    mtl::dense_vector<T> grd1(dim + 1);
    int parts = Global::getGeo(PARTS, dim);
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();

    derivativeAtQPs.change_dim(nPoints);
    if (quadFast) {
      for (int iq = 0; iq < nPoints; iq++) {
	nullify(grd1);

	for (int j = 0; j < nBasFcts; j++) // #BasisFunctions
	  for (int k = 0; k < parts; k++)  // #edges (2d) or #faces (3d)
	    grd1[k] += quadFast->getGradient(iq, j, k) * localVec[j];

	  nullify(derivativeAtQPs[iq]);
	  for (int k = 0; k < parts; k++)
	    derivativeAtQPs[iq] += grdLambda[k][comp] * grd1[k];
      }

    } else {
      mtl::dense_vector<double> grdPhi(dim + 1);

      for (int iq = 0; iq < nPoints; iq++) {
	nullify(grd1);

	for (int j = 0; j < nBasFcts; j++) {
	  (*(basFcts->getGrdPhi(j)))(quad->getLambda(iq), grdPhi);
	  for (int k = 0; k < parts; k++)
	    grd1[k] += grdPhi[k] * localVec[j];
	}

	  nullify(derivativeAtQPs[iq]);
	  for (int k = 0; k < parts; k++)
	    derivativeAtQPs[iq] += grdLambda[k][comp] * grd1[k];
      }
    }
  }


  template <class T>
  inline void set_to_zero(AMDiS::DOFVector<T>& v)
  {
    T my_zero; nullify(my_zero);
    std::fill(v.begin(), v.end(), my_zero);
  }


  template <class T>
  double DOFVector<T>::DoubleWell(Quadrature* q) const
  {
    using math::sqr;
    Mesh* mesh = this->feSpace->getMesh();

    if (!q) {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(this->dim, deg);
    }

    FastQuadrature *quadFast = 
      FastQuadrature::provideFastQuadrature(this->feSpace->getBasisFcts(), *q, INIT_PHI);

    double result = 0.0;
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<T> uh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, -1,
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET);
    while (elInfo) {
      double det = elInfo->getDet();
      double normT = 0.0;
      this->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      for (int iq = 0; iq < nPoints; iq++)
	normT += quadFast->getWeight(iq) * sqr(uh_vec[iq]) * sqr(1.0 - uh_vec[iq]);
      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double localResult = result;
    MPI::COMM_WORLD.Allreduce(&localResult, &result, 1, MPI_DOUBLE, MPI_SUM);
#endif
    
    return result;
  }


  template <class T>
  DOFVector<typename GradientType<T>::type>*
  DOFVector<T>::getGradient(DOFVector<typename GradientType<T>::type> *grad) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGradient()");
    const FiniteElemSpace *feSpace = DOFVector<T>::feSpace;

    // define result vector
    static DOFVector<typename GradientType<T>::type> *result = NULL; // TODO: REMOVE STATIC

    if (grad) {
      result = grad;
    } else {
      if(result && result->getFeSpace() != feSpace) {
	delete result;
	result = new DOFVector<typename GradientType<T>::type>(feSpace, "gradient");
      }
    }

    Mesh *mesh = feSpace->getMesh();
    int dim = mesh->getDim();
    const BasisFunction *basFcts = feSpace->getBasisFcts();
    DOFAdmin *admin = feSpace->getAdmin();

    // count number of nodes and dofs per node
    std::vector<int> nNodeDOFs;
    std::vector<int> nNodePreDofs;
    std::vector<DimVec<double>*> bary;

    int nNodes = 0;
    int nDofs = 0;

    for (int i = 0; i < dim + 1; i++) {
      GeoIndex geoIndex = INDEX_OF_DIM(i, dim);
      int nPositions = mesh->getGeo(geoIndex);
      int numPreDofs = admin->getNumberOfPreDofs(i);
      for (int j = 0; j < nPositions; j++) {
	int dofs = basFcts->getNumberOfDofs(geoIndex);
	nNodeDOFs.push_back(dofs);
	nDofs += dofs;
	nNodePreDofs.push_back(numPreDofs);
      }
      nNodes += nPositions;
    }

    TEST_EXIT_DBG(nDofs == basFcts->getNumber())
      ("number of dofs != number of basis functions\n");

    for (int i = 0; i < nDofs; i++)
      bary.push_back(basFcts->getCoords(i));

    mtl::dense_vector<T> localUh(basFcts->getNumber());

    // traverse mesh
    std::vector<bool> visited(getUsedSize(), false);
    TraverseStack stack;
    Flag fillFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while (elInfo) {
      const DegreeOfFreedom **dof = elInfo->getElement()->getDof();
      const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
      this->getLocalVector(elInfo->getElement(), localUh);

      int localDOFNr = 0;
      for (int i = 0; i < nNodes; i++) { // for all nodes
	for (int j = 0; j < nNodeDOFs[i]; j++) { // for all dofs at this node
	  DegreeOfFreedom dofIndex = dof[i][nNodePreDofs[i] + j];
	  if (!visited[dofIndex]) {
	    basFcts->evalGrdUh(*(bary[localDOFNr]), grdLambda,
			       localUh, ((*result)[dofIndex]));

	    visited[dofIndex] = true;
	  }
	  localDOFNr++;
	}
      }

      elInfo = stack.traverseNext(elInfo);
    }

    return result;
  }


  template <class T>
  DOFVector<typename GradientType<T>::type>*
  DOFVector<T>::getRecoveryGradient(DOFVector<typename GradientType<T>::type> *grad) const
  {
    const FiniteElemSpace *feSpace = DOFVector<T>::feSpace;
    int dim = DOFVector<T>::dim;

    // define result vector
    static DOFVector<typename GradientType<T>::type> *vec = NULL; // TODO: REMOVE STATIC

    DOFVector<typename GradientType<T>::type> *result = grad;

    if (!result) {
      if (vec && vec->getFeSpace() != feSpace) {
	delete vec;
	vec = NULL;
      }
      if (!vec)
	vec = new DOFVector<typename GradientType<T>::type>(feSpace, "gradient");

      result = vec;
    }

    typename GradientType<T>::type grd;
    nullify(grd);

    result->set(grd);

    DOFVector<double> volume(feSpace, "volume");
    volume.set(0.0);

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int nBasisFcts = basFcts->getNumber();
    DimVec<double> bary(dim, DEFAULT_VALUE, (1.0 / (dim + 1.0)));

    // traverse mesh
    Mesh *mesh = feSpace->getMesh();
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1,
					 Mesh::CALL_LEAF_EL | Mesh::FILL_DET |
					 Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS);

    mtl::dense_vector<T> localUh(nBasisFcts);
    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);

    while (elInfo) {
      double det = elInfo->getDet();
      const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
      this->getLocalVector(elInfo->getElement(), localUh);
      basFcts->evalGrdUh(bary, grdLambda, localUh, grd);
      basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), localIndices);

      for (int i = 0; i < nBasisFcts; i++) {
	(*result)[localIndices[i]] += grd * det;
	volume[localIndices[i]] += det;
      }

      elInfo = stack.traverseNext(elInfo);
    
    }
   
    // NOTE: We have to synchronize the vectors in PARALLEL_DOMAIN_AMDIS mode
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->checkMeshChange(false);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(*result);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(volume);
#endif
    
    DOFVector<double>::Iterator volIt(&volume, USED_DOFS);
    typename DOFVector<typename GradientType<T>::type>::Iterator grdIt(result, USED_DOFS);

    for (volIt.reset(), grdIt.reset(); !volIt.end(); ++volIt, ++grdIt)
      if (*volIt != 0.0)
	*grdIt *= 1.0 / (*volIt);
      
    return result;
  }


  template <class T>
  std::vector<DOFVector<double>*> *transform(DOFVector<typename GradientType<T>::type> *vec,
					     std::vector<DOFVector<double>*> *res)
  {
    FUNCNAME_DBG("DOFVector<T>::transform()");

    TEST_EXIT_DBG(vec)("no vector\n");

    static std::vector<DOFVector<double>*> *result = NULL; // TODO: REMOVE STATIC

    int len = num_rows(GradientType<T>::getValues((*vec)[0]));
    
    if (!res && !result) {
      result = new std::vector<DOFVector<double>*>(len);
      for (int i = 0; i < len; i++)
	(*result)[i] = new DOFVector<double>(vec->getFeSpace(), "transform");
    } else if (res->size() == 0 || (*res)[0] == NULL) {
      res->resize(len, NULL);
      for (int i = 0; i < len; i++)
	(*res)[i] = new DOFVector<double>(vec->getFeSpace(), "transform");
    }

    std::vector<DOFVector<double>*> *r = res ? res : result;

    int vecSize = vec->getSize();
    for (int i = 0; i < vecSize; i++)
      for (int j = 0; j < len; j++)
	(*((*r)[j]))[i] = (GradientType<T>::getValues((*vec)[i]))[j];

    return r;
  }
  
  
  
  

  template<typename T>
  const DOFVector<T>& operator*=(DOFVector<T>& x, T scal)
  {
    FUNCNAME_DBG("DOFVector<T>::operator*=(DOFVector<T>& x, T scal)");

    TEST_EXIT_DBG(x.getFeSpace() && x.getFeSpace()->getAdmin())
      ("pointer is NULL: %8X, %8X\n", x.getFeSpace(), x.getFeSpace()->getAdmin());

    typename DOFVector<T>::Iterator vecIterator(dynamic_cast<DOFIndexed<T>*>(&x), 
						USED_DOFS);
    for (vecIterator.reset(); !vecIterator.end(); ++vecIterator)
      (*vecIterator) *= scal; 

    return x;
  }


  template<typename T>
  const DOFVector<T>& operator+=(DOFVector<T>& x, const DOFVector<T>& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator+=(DOFVector<T>& x, const DOFVector<T>& y)");
    
    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    typename DOFVector<T>::Iterator xIterator(dynamic_cast<DOFIndexed<T>*>(&x), USED_DOFS);
    typename DOFVector<T>::Iterator yIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(&y)), USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator)
      *xIterator += *yIterator; 

    return x;
  }


  template<typename T>
  const DOFVector<T>& operator-=(DOFVector<T>& x, const DOFVector<T>& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator-=(DOFVector<T>& x, const DOFVector<T>& y)");

    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    typename DOFVector<T>::Iterator xIterator(dynamic_cast<DOFIndexed<T>*>(&x), USED_DOFS);
    typename DOFVector<T>::Iterator yIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(&y)), USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator)
      *xIterator -= *yIterator; 

    return x;
  }


  template<typename T>
  const DOFVector<T>& operator*=(DOFVector<T>& x, const DOFVector<T>& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator*=(DOFVector<T>& x, const DOFVector<T>& y)");
    
    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    typename DOFVector<T>::Iterator xIterator(dynamic_cast<DOFIndexed<T>*>(&x), USED_DOFS);
    typename DOFVector<T>::Iterator yIterator(dynamic_cast<DOFIndexed<T>*>(const_cast<DOFVector<T>*>(&y)), USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator)
      *xIterator *= *yIterator; 

    return x;
  }


  template<typename T>
  T operator*(DOFVector<T>& x, DOFVector<T>& y)
  {
    FUNCNAME("DOFVector<T>::operator*(DOFVector<T>& x, DOFVector<T>& y)");
    const DOFAdmin *admin = NULL;

    TEST_EXIT(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT((admin = x.getFeSpace()->getAdmin()) && (admin == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT(x.getSize() == y.getSize())("different sizes\n");

    T dot = 0;

    typename DOFVector<T>::Iterator xIterator(dynamic_cast<DOFIndexed<T>*>(&x), USED_DOFS);
    typename DOFVector<T>::Iterator yIterator(dynamic_cast<DOFIndexed<T>*>(&y), USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end();
	 ++xIterator, ++yIterator) 
      dot += (*xIterator) * (*yIterator);

    return dot;
  }
  

  template<typename T>
  const DOFVector<T>& operator*(const DOFVector<T>& v, double d)
  {
    static DOFVector<T> result; // TODO: REMOVE STATIC
    return mult(d, v, result); 
  }


  template<typename T>
  const DOFVector<T>& operator*(double d, const DOFVector<T>& v)
  {
    static DOFVector<T> result; // TODO: REMOVE STATIC
    return mult(d, v, result);
  }


  template<typename T>
  const DOFVector<T>& operator+(const DOFVector<T> &v1 , const DOFVector<T> &v2)
  {
    static DOFVector<T> result; // TODO: REMOVE STATIC
    return add(v1, v2, result);
  }

  
} // end namespace AMDiS
