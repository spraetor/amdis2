/** \file DOFVectorBase.hh */

#include <vector>

#include "BasisFunction.h"
#include "Boundary.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Global.h"
#include "Mesh.h"
#include "Operator.h"
#include "Quadrature.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MpiHelper.h"
#endif

namespace AMDiS
{
  template <class T>
  DOFVectorBase<T>::DOFVectorBase(const FiniteElemSpace* f, std::string n)
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
  void DOFVectorBase<T>::addElementVector(T factor,
                                          const DenseVector<double>& elVec,
                                          const BoundaryType* bound,
                                          ElInfo* elInfo,
                                          bool add)
  {
    std::vector<DegreeOfFreedom> indices(nBasFcts);
    feSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(),
        indices);

    for (int i = 0; i < nBasFcts; i++)
    {
      BoundaryCondition* condition =
        bound ? this->getBoundaryManager()->getBoundaryCondition(bound[i]) : NULL;

      if (!(condition && condition->isDirichlet()))
      {
        DegreeOfFreedom irow = indices[i];

        if (add)
          (*this)[irow] += factor * elVec[i];
        else
          (*this)[irow] = factor * elVec[i];
      }
    }
  }


  template <class T>
  T DOFVectorBase<T>::evalUh(DimVec<double> const& lambda,
                             DegreeOfFreedom* dof_indices)
  {
    BasisFunction const* phi = this->getFeSpace()->getBasisFcts();
    int nBasisFcts = phi->getNumber();
    T val = 0.0;

    for (int i = 0; i < nBasisFcts; i++)
      val += (*this)[dof_indices[i]]*(*phi->getPhi(i))(lambda);

    // TODO: ist das im Parallelen so richtig???
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(val);
#endif

    return val;
  }


  template <class T>
  Flag DOFVectorBase<T>::getAssembleFlag()
  {
    Flag fillFlag(0);

    for (Operator* op : operators)
      fillFlag |= op->getFillFlag();

    return fillFlag;
  }


  template <class T>
  void DOFVectorBase<T>::finishAssembling()
  {
    // call the operatos cleanup procedures
    for (Operator* op : operators)
      op->finishAssembling();
  }


  template <class T>
  void DOFVectorBase<T>::getLocalVector(Element const* el,
                                        DenseVector<T>& d) const
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
  void DOFVectorBase<T>::getVecAtQPs(ElInfo const* elInfo,
                                     Quadrature const* quad,
                                     FastQuadrature const* quadFast,
                                     DenseVector<T>& vecAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getVecAtQPs()");

    TEST_EXIT_DBG(quad || quadFast)
    ("Neither quad nor quadFast defined!\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
    ("quad != quadFast->quadrature\n");
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
    ("Invalid basis functions!");

    const BasisFunction* basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    DenseVector<T> localVec(nBasFcts);
    getLocalVector(elInfo->getElement(), localVec);

    if (quadFast)
    {
      // using precalculated values at QPs
      auto& phi = quadFast->getPhi();
      vecAtQPs.change_dim(num_rows(phi));   // = quadrature->getNumPoints()
      vecAtQPs = phi * localVec;    // Matrix<double> * Vector<T>
    }
    else
    {
      // evaluate basisFunctions at QPs
      int nPoints = quad->getNumPoints();
      vecAtQPs.change_dim(nPoints);

      for (int iq = 0; iq < nPoints; iq++)
      {
        nullify(vecAtQPs[iq]);
        for (int j = 0; j < nBasFcts; j++)
          vecAtQPs[iq] +=
            localVec[j] * (*(basFcts->getPhi(j)))(quad->getLambda(iq));
      }
    }
  }


  template <class T>
  void DOFVectorBase<T>::getGrdAtQPs(ElInfo const* elInfo,
                                     Quadrature const* quad,
                                     FastQuadrature const* quadFast,
                                     DenseVector<Gradient_t<T>>& grdAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGrdAtQPs()");

    TEST_EXIT_DBG(quad || quadFast)("neither quad nor quadFast defined\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
    ("quad != quadFast->quadrature\n");
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
    ("invalid basis functions");

    const BasisFunction* basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    int dow = Global::getGeo(WORLD);
    int nPoints = quadFast ? quadFast->getQuadrature()->getNumPoints() : quad->getNumPoints();

    DenseVector<T> localVec(nBasFcts);
    this->getLocalVector(elInfo->getElement(), localVec);

    DenseVector<T> grd1(dim + 1);
    int parts = Global::getGeo(PARTS, dim);
    auto& grdLambda = elInfo->getGrdLambda();

    grdAtQPs.change_dim(nPoints);
    if (quadFast)
    {
      for (int iq = 0; iq < nPoints; iq++)
      {
        nullify(grd1);

        for (int j = 0; j < nBasFcts; j++)   // #BasisFunctions
        {
          for (int k = 0; k < parts; k++)  // #edges (2d) or #faces (3d)
            grd1[k] += quadFast->getGradient(iq, j, k) * localVec[j];
        }

        for (int l = 0; l < dow; l++)
        {
          nullify(grdAtQPs[iq][l]);
          for (int k = 0; k < parts; k++)
            grdAtQPs[iq][l] += grdLambda[k][l] * grd1[k];
        }
      }

    }
    else
    {
      DenseVector<double> grdPhi(dim + 1);

      for (int iq = 0; iq < nPoints; iq++)
      {
        nullify(grd1);

        for (int j = 0; j < nBasFcts; j++)
        {
          (*(basFcts->getGrdPhi(j)))(quad->getLambda(iq), grdPhi);
          for (int k = 0; k < parts; k++)
            grd1[k] += grdPhi[k] * localVec[j];
        }

        for (int l = 0; l < dow; l++)
        {
          nullify(grdAtQPs[iq][l]);
          for (int k = 0; k < parts; k++)
            grdAtQPs[iq][l] += grdLambda[k][l] * grd1[k];
        }
      }
    }
  }


  template <class T>
  void DOFVectorBase<T>::getDerivativeAtQPs(ElInfo const* elInfo,
      Quadrature const* quad,
      FastQuadrature const* quadFast,
      int comp,
      DenseVector<T>& derivativeAtQPs) const
  {
    FUNCNAME_DBG("DOFVector<T>::getGrdAtQPs()");

    TEST_EXIT_DBG(quad || quadFast)("neither quad nor quadFast defined\n");
    TEST_EXIT_DBG(!(quad && quadFast) || quad == quadFast->getQuadrature())
    ("quad != quadFast->quadrature\n");
    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
    ("invalid basis functions");

    const BasisFunction* basFcts = feSpace->getBasisFcts();
    int nBasFcts  = basFcts->getNumber();
    int nPoints = quadFast ? quadFast->getQuadrature()->getNumPoints() : quad->getNumPoints();

    DenseVector<T> localVec(nBasFcts);
    this->getLocalVector(elInfo->getElement(), localVec);

    DenseVector<T> grd1(dim + 1);
    int parts = Global::getGeo(PARTS, dim);
    auto& grdLambda = elInfo->getGrdLambda();

    derivativeAtQPs.change_dim(nPoints);
    if (quadFast)
    {
      for (int iq = 0; iq < nPoints; iq++)
      {
        nullify(grd1);

        for (int j = 0; j < nBasFcts; j++) // #BasisFunctions
          for (int k = 0; k < parts; k++)  // #edges (2d) or #faces (3d)
            grd1[k] += quadFast->getGradient(iq, j, k) * localVec[j];

        nullify(derivativeAtQPs[iq]);
        for (int k = 0; k < parts; k++)
          derivativeAtQPs[iq] += grdLambda[k][comp] * grd1[k];
      }

    }
    else
    {
      DenseVector<double> grdPhi(dim + 1);

      for (int iq = 0; iq < nPoints; iq++)
      {
        nullify(grd1);

        for (int j = 0; j < nBasFcts; j++)
        {
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

} // end namespace AMDiS
