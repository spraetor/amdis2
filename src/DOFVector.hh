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
// #include "Expressions.h"
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
#include "MatrixVectorOperations.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MpiHelper.h"
#include "parallel/MeshDistributor.h"
#endif

#include <expressions2/_Expressions.h>

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
  void DOFVector<T>::copy(DOFVector<T> const& x)
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
  T DOFVector<T>::Int(int meshLevel, Quadrature* q) const
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
    DenseVector<T> uh_vec(nPoints);
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
    DenseVector<T> uh_vec(nPoints);
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET);
    while (elInfo) {
      double det = elInfo->getDet();
      double normT = 0.0;
      this->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      for (int iq = 0; iq < nPoints; iq++)
        normT += quadFast->getWeight(iq) * std::abs(uh_vec[iq]);
      result += det * normT;

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(result);
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
    DenseVector<T> uh_vec(nPoints);
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
    Parallel::mpi::globalAdd(result);
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
    DenseVector<WorldVector<T> > grduh_vec(nPoints);
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
    Parallel::mpi::globalAdd(result);
#endif
    
    return result;
  }



  template <class T>
  bool DOFVector<T>::getDofIdxAtPoint(WorldVector<double>& p, 
                                      DegreeOfFreedom& idx, 
                                      ElInfo* oldElInfo, 
                                      bool useOldElInfo) const
  { 
    Mesh *mesh = this->feSpace->getMesh();
    const BasisFunction *basFcts = this->feSpace->getBasisFcts();

    int dim = mesh->getDim();
    int numBasFcts = basFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(numBasFcts);
    DimVec<double> lambda(dim);

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
      double newDist = norm(coord - p);
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
                                        std::vector<DegreeOfFreedom>& newDOF)
  {
    for (int i = first; i <= last; i++)
      if (newDOF[i] >= 0)
        vec[newDOF[i]] = vec[i];
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
    DenseVector<T> uh_vec(nPoints);
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
    Parallel::mpi::globalAdd(result);
#endif
    
    return result;
  }


  template <class T>
  DOFVector<Gradient_t<T>>*
  DOFVector<T>::getGradient(DOFVector<Gradient_t<T>> *grad) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGradient()");
    const FiniteElemSpace *feSpace = DOFVector<T>::feSpace;

    // define result vector
    static DOFVector<Gradient_t<T>> *result = NULL; // TODO: REMOVE STATIC

    if (grad) {
      result = grad;
    } else {
      if(result && result->getFeSpace() != feSpace) {
        delete result;
        result = new DOFVector<Gradient_t<T>>(feSpace, "gradient");
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

    DenseVector<T> localUh(basFcts->getNumber());

    // traverse mesh
    std::vector<bool> visited(getUsedSize(), false);
    TraverseStack stack;
    Flag fillFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    while (elInfo) {
      const DegreeOfFreedom **dof = elInfo->getElement()->getDof();
      auto& grdLambda = elInfo->getGrdLambda();
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
  DOFVector<Gradient_t<T>>*
  DOFVector<T>::getRecoveryGradient(DOFVector<Gradient_t<T>> *grad) const
  {
    const FiniteElemSpace *feSpace = DOFVector<T>::feSpace;
    int dim = DOFVector<T>::dim;

    // define result vector
    static DOFVector<Gradient_t<T>> *vec = NULL; // TODO: REMOVE STATIC

    DOFVector<Gradient_t<T>> *result = grad;

    if (!result) {
      if (vec && vec->getFeSpace() != feSpace) {
        delete vec;
        vec = NULL;
      }
      if (!vec)
        vec = new DOFVector<Gradient_t<T>>(feSpace, "gradient");

      result = vec;
    }

    Gradient_t<T> grd;
    nullify(grd);

    result->set(grd);

    DOFVector<double> volume(feSpace, "volume");
    volume.set(0.0);

    const BasisFunction *basFcts = feSpace->getBasisFcts();
    int nBasisFcts = basFcts->getNumber();
    DimVec<double> bary(dim, (1.0 / (dim + 1.0)));

    // traverse mesh
    Mesh *mesh = feSpace->getMesh();
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1,
					 Mesh::CALL_LEAF_EL | Mesh::FILL_DET |
					 Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS);

    DenseVector<T> localUh(nBasisFcts);
    std::vector<DegreeOfFreedom> localIndices(nBasisFcts);

    while (elInfo) {
      double det = elInfo->getDet();
      auto& grdLambda = elInfo->getGrdLambda();
      this->getLocalVector(elInfo->getElement(), localUh);
      basFcts->evalGrdUh(bary, grdLambda, localUh, grd);
      basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), localIndices);

      for (int i = 0; i < nBasisFcts; i++) {
        (*result)[localIndices[i]] += grd * det;
        volume[localIndices[i]] += det;
      }

      elInfo = stack.traverseNext(elInfo);
    
    }
   
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::MeshDistributor::globalMeshDistributor->checkMeshChange(false);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(*result);
    Parallel::MeshDistributor::globalMeshDistributor->synchAddVector(volume);
#endif
    
    DOFVector<double>::Iterator volIt(&volume, USED_DOFS);
    typename DOFVector<Gradient_t<T>>::Iterator grdIt(result, USED_DOFS);

    for (volIt.reset(), grdIt.reset(); !volIt.end(); ++volIt, ++grdIt)
      if (*volIt != 0.0)
        *grdIt *= 1.0 / (*volIt);
      
    return result;
  }


  template <class T>
  std::vector<DOFVector<double>*> *transform(DOFVector<Gradient_t<T>>* vec,
                                             std::vector<DOFVector<double>*>* res)
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
  DOFVector<T>& operator*=(DOFVector<T>& x, T scal)
  {
    FUNCNAME_DBG("DOFVector<T>::operator*=(DOFVector<T>& x, T scal)");

    TEST_EXIT_DBG(x.getFeSpace() && x.getFeSpace()->getAdmin())
      ("pointer is NULL: %8X, %8X\n", x.getFeSpace(), x.getFeSpace()->getAdmin());

    DOFIterator<T> xIterator(&x, USED_DOFS);
    for (xIterator.reset(); !xIterator.end(); ++xIterator)
      (*xIterator) *= scal; 

    return x;
  }


  template<typename T>
  DOFVector<T>& operator+=(DOFVector<T>& x, DOFVector<T> const& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator+=(DOFVector<T>& x, const DOFVector<T>& y)");
    
    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    DOFIterator<T> xIterator(&x, USED_DOFS);
    DOFConstIterator<T> yIterator(&y, USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end(); ++xIterator, ++yIterator)
      *xIterator += *yIterator; 

    return x;
  }


  template<typename T>
  DOFVector<T>& operator-=(DOFVector<T>& x, DOFVector<T> const& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator-=(DOFVector<T>& x, const DOFVector<T>& y)");

    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    DOFIterator<T> xIterator(&x, USED_DOFS);
    DOFConstIterator<T> yIterator(&y, USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end(); ++xIterator, ++yIterator)
      *xIterator -= *yIterator; 

    return x;
  }


  template<typename T>
  DOFVector<T>& operator*=(DOFVector<T>& x, DOFVector<T> const& y)
  {
    FUNCNAME_DBG("DOFVector<T>::operator*=(DOFVector<T>& x, const DOFVector<T>& y)");
    
    TEST_EXIT_DBG(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT_DBG(x.getFeSpace()->getAdmin() &&
	      (x.getFeSpace()->getAdmin() == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT_DBG(x.getSize() == y.getSize())("different sizes\n");
    
    DOFIterator<T> xIterator(&x, USED_DOFS);
    DOFConstIterator<T> yIterator(&y, USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end(); ++xIterator, ++yIterator)
      *xIterator *= *yIterator; 

    return x;
  }


  template<typename T>
  T operator*(DOFVector<T> const& x, DOFVector<T> const& y)
  {
    FUNCNAME("DOFVector<T>::operator*(DOFVector<T>& x, DOFVector<T>& y)");
    const DOFAdmin *admin = NULL;

    TEST_EXIT(x.getFeSpace() && y.getFeSpace())
      ("feSpace is NULL: %8X, %8X\n", x.getFeSpace(), y.getFeSpace());
    TEST_EXIT((admin = x.getFeSpace()->getAdmin()) && (admin == y.getFeSpace()->getAdmin()))
      ("no admin or different admins: %8X, %8X\n",
       x.getFeSpace()->getAdmin(), y.getFeSpace()->getAdmin());
    TEST_EXIT(x.getSize() == y.getSize())("different sizes\n");

    T dot; nullify(dot);
    
    DOFConstIterator<T> xIterator(&x, USED_DOFS);
    DOFConstIterator<T> yIterator(&y, USED_DOFS);
    for (xIterator.reset(), yIterator.reset(); !xIterator.end(); ++xIterator, ++yIterator) 
      dot += (*xIterator) * (*yIterator);

    return dot;
  }
  

  template<typename T>
  DOFVector<T> operator*(DOFVector<T> v, T d)
  {
    v *= d;
    return v; 
  }


  template<typename T>
  DOFVector<T>& operator*(T d, DOFVector<T> v)
  {
    v *= d;
    return v; 
  }


  template<typename T>
  DOFVector<T> operator+(DOFVector<T> v1 , const DOFVector<T> &v2)
  {
    v1 += v2;
    return v1;
  }

  template<typename T>
  DOFVector<T> operator-(DOFVector<T> v1 , const DOFVector<T> &v2)
  {
    v1 -= v2;
    return v1;
  }
  
} // end namespace AMDiS
