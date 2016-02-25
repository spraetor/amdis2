/** \file DOFVectorBase.cc */

#include "DOFVectorBase.hpp"
#include "Traverse.hpp"
#include "DualTraverse.hpp"
#include "FixVec.hpp"
#include "ElementFunction.hpp"

namespace AMDiS
{
  template<>
  void DOFVectorBase<double>::getD2AtQPs( const ElInfo* elInfo,
                                          const Quadrature* quad,
                                          const FastQuadrature* quadFast,
                                          DenseVector<D2Type<double>::type>& d2AtQPs) const
  {
    FUNCNAME("DOFVector<double>::getD2AtQPs()");

    TEST_EXIT_DBG(quad || quadFast)("neither quad nor quadFast defined\n");

    if (quad && quadFast)
    {
      TEST_EXIT_DBG(quad == quadFast->getQuadrature())
      ("quad != quadFast->quadrature\n");
    }

    TEST_EXIT_DBG(!quadFast || quadFast->getBasisFunctions() == feSpace->getBasisFcts())
    ("invalid basis functions");

    Element* el = elInfo->getElement();

    int dow = Global::getGeo(WORLD);
    int nPoints = quadFast ? quadFast->getQuadrature()->getNumPoints() : quad->getNumPoints();

    DenseVector<double> localVec(nBasFcts);
    getLocalVector(el, localVec);

    DimMat<double> D2Tmp(dim, dim, 0.0);
    int parts = Global::getGeo(PARTS, dim);
    const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();

    d2AtQPs.change_dim(nPoints);
    if (quadFast)
    {
      for (int iq = 0; iq < nPoints; iq++)
      {
        for (int k = 0; k < parts; k++)
          for (int l = 0; l < parts; l++)
            D2Tmp[k][l] = 0.0;

        for (int i = 0; i < nBasFcts; i++)
        {
          for (int k = 0; k < parts; k++)
            for (int l = 0; l < parts; l++)
              D2Tmp[k][l] += localVec[i] * quadFast->getSecDer(iq, i, k, l);
        }

        for (int i = 0; i < dow; i++)
          for (int j = 0; j < dow; j++)
          {
            d2AtQPs[iq][i][j] = 0.0;
            for (int k = 0; k < parts; k++)
              for (int l = 0; l < parts; l++)
                d2AtQPs[iq][i][j] += grdLambda[k][i]*grdLambda[l][j]*D2Tmp[k][l];
          }
      }
    }
    else
    {
      const BasisFunction* basFcts = feSpace->getBasisFcts();
      DimMat<double> D2Phi(dim, dim);

      for (int iq = 0; iq < nPoints; iq++)
      {
        for (int k = 0; k < parts; k++)
          for (int l = 0; l < parts; l++)
            D2Tmp[k][l] = 0.0;

        for (int i = 0; i < nBasFcts; i++)
        {
          WARNING("not tested after index correction\n");
          (*(basFcts->getD2Phi(i)))(quad->getLambda(iq), D2Phi);

          for (int k = 0; k < parts; k++)
            for (int l = 0; l < parts; l++)
              D2Tmp[k][l] += localVec[i] * D2Phi[k][l];
        }

        for (int i = 0; i < dow; i++)
          for (int j = 0; j < dow; j++)
          {
            d2AtQPs[iq][i][j] = 0.0;
            for (int k = 0; k < parts; k++)
              for (int l = 0; l < parts; l++)
                d2AtQPs[iq][i][j] += grdLambda[k][i] * grdLambda[l][j] * D2Tmp[k][l];
          }
      }
    }
  }


  template<>
  void DOFVectorBase<double>::assemble(double factor, ElInfo* elInfo,
                                       const BoundaryType* bound,
                                       Operator* op)
  {
    if (!(op || operators.size()))
      return;

    set_to_zero(elementVector);
    bool addVector = false;

    if (op)
    {
      op->getElementVector(elInfo, elementVector);
      addVector = true;
    }
    else
    {
      std::vector<double*>::iterator factorIt = operatorFactor.begin();

      for (Operator* op : operators)
      {
        if (!op->getNeedDualTraverse())
        {
          op->getElementVector(elInfo, elementVector, (*factorIt ? **factorIt : 1.0) );
          addVector = true;
        }
        ++factorIt;
      }
    }

    if (addVector)
      addElementVector(factor, this->elementVector, bound, elInfo);
  }


  template<>
  void DOFVectorBase<double>::assembleOperator(Operator& op)
  {
    FUNCNAME("DOFVectorBase::assembleOperator()");

    TEST_EXIT(op.getRowFeSpace() == feSpace)
    ("Row FE spaces do not fit together!\n");

    Mesh* mesh = feSpace->getMesh();
    mesh->dofCompress();
    const BasisFunction* basisFcts = feSpace->getBasisFcts();

    Flag assembleFlag = getAssembleFlag() |
                        Mesh::CALL_LEAF_EL                        |
                        Mesh::FILL_COORDS                         |
                        Mesh::FILL_DET                            |
                        Mesh::FILL_GRD_LAMBDA |
                        Mesh::FILL_NEIGH |
                        Mesh::FILL_BOUND;

    BoundaryType* bound = new BoundaryType[basisFcts->getNumber()];

    if (getBoundaryManager())
      getBoundaryManager()->initVector(this);


    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, assembleFlag);
    while (elInfo)
    {
      basisFcts->getBound(elInfo, bound);

      assemble(1.0, elInfo, bound, &op);

      if (getBoundaryManager())
        getBoundaryManager()->fillBoundaryConditions(elInfo, this);

      elInfo = stack.traverseNext(elInfo);
    }

    finishAssembling();
    getBoundaryManager()->exitVector(this);

    delete [] bound;
  }

} // end namespace AMDiS
