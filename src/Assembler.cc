#include <vector>
#include <algorithm>

#include <boost/numeric/mtl/mtl.hpp>

#include "Assembler.h"
#include "Operator.h"
#include "Element.h"
#include "QPsiPhi.h"
#include "DOFVector.h"
#include "ElInfo.h"

namespace AMDiS
{
  Assembler::Assembler(not_null<Operator*> op,
                       not_null<const FiniteElemSpace*> row,
				const FiniteElemSpace*  col)
    : operat(op),
      rowFeSpace(row),
      colFeSpace(col ? col : row.get()),
      nRow(rowFeSpace->getBasisFcts()->getNumber()),
      nCol(colFeSpace->getBasisFcts()->getNumber()),
      remember(true),
      rememberElMat(false),
      rememberElVec(false),
      elementMatrix(nRow, nCol),
      elementVector(nRow),
      tmpMat(nRow, nCol),
      lastMatEl(NULL),
      lastVecEl(NULL),
      lastTraverseId(-1)
  {}


  void Assembler::calculateElementMatrix(const ElInfo* elInfo,
                                         ElementMatrix& userMat,
                                         double factor)
  {
    if (remember && (factor != 1.0 || operat->getUhOld()))
      rememberElMat = true;

    Element* el = elInfo->getElement();

    if (el != lastMatEl || !operat->isOptimized())
    {
      initElement(elInfo);

      if (rememberElMat)
        set_to_zero(elementMatrix);

      lastMatEl = el;
    }
    else
    {
      // Only possible in single mesh case when one operator
      // is used more than twice?
      if (rememberElMat)
      {
        if (&userMat != &elementMatrix)
          userMat += factor * elementMatrix;
        return;
      }
    }

    ElementMatrix& mat = rememberElMat ? elementMatrix : userMat;

    if (secondOrderAssembler)
      secondOrderAssembler->calculateElementMatrix(elInfo, mat);
    if (firstOrderAssemblerGrdPsi)
      firstOrderAssemblerGrdPsi->calculateElementMatrix(elInfo, mat);
    if (firstOrderAssemblerGrdPhi)
      firstOrderAssemblerGrdPhi->calculateElementMatrix(elInfo, mat);
    if (zeroOrderAssembler)
      zeroOrderAssembler->calculateElementMatrix(elInfo, mat);

    if (rememberElMat && &userMat != &elementMatrix)
      userMat += factor * elementMatrix;
  }


  void Assembler::calculateElementVector(const ElInfo* elInfo,
                                         DenseVector<double>& userVec,
                                         double factor)
  {
    if (remember && factor != 1.0)
      rememberElVec = true;

    Element* el = elInfo->getElement();

    if ((el != lastMatEl && el != lastVecEl) || !operat->isOptimized())
      initElement(elInfo);

    if (el != lastVecEl || !operat->isOptimized())
    {
      if (rememberElVec)
        set_to_zero(elementVector);

      lastVecEl = el;
    }
    else
    {
      // Only possible in single mesh case when one operator
      // is used more than twice at dof vector?
      if (rememberElVec)
      {
        userVec += factor * elementVector;
        return;
      }
    }

    DenseVector<double>& vec = rememberElVec ? elementVector : userVec;

    if (operat->getUhOld() && remember)
    {
      matVecAssemble(elInfo, vec);
      if (rememberElVec)
        userVec += factor * elementVector;

      return;
    }

    if (firstOrderAssemblerGrdPsi)
      firstOrderAssemblerGrdPsi->calculateElementVector(elInfo, vec);
    if (zeroOrderAssembler)
      zeroOrderAssembler->calculateElementVector(elInfo, vec);

    if (rememberElVec)
      userVec += factor * elementVector;
  }


  void Assembler::matVecAssemble(const ElInfo* elInfo, DenseVector<double>& vec)
  {

    Element* el = elInfo->getElement();
    DenseVector<double> uhOldLoc(operat->getUhOld()->getFeSpace() == rowFeSpace ?
                                 nRow : nCol);
    operat->getUhOld()->getLocalVector(el, uhOldLoc);

    if (el != lastMatEl)
    {
      set_to_zero(elementMatrix);
      calculateElementMatrix(elInfo, elementMatrix);
    }


    vec += elementMatrix * uhOldLoc;
  }


  void Assembler::initElement(const ElInfo* smallElInfo,
                              const ElInfo* largeElInfo,
                              Quadrature* quad)
  {
    if (secondOrderAssembler)
      secondOrderAssembler->initElement(smallElInfo, largeElInfo, quad);
    if (firstOrderAssemblerGrdPsi)
      firstOrderAssemblerGrdPsi->initElement(smallElInfo, largeElInfo, quad);
    if (firstOrderAssemblerGrdPhi)
      firstOrderAssemblerGrdPhi->initElement(smallElInfo, largeElInfo, quad);
    if (zeroOrderAssembler)
      zeroOrderAssembler->initElement(smallElInfo, largeElInfo, quad);
  }


  void Assembler::checkQuadratures()
  {
    if (secondOrderAssembler)
    {
      // create quadrature
      if (!secondOrderAssembler->getQuadrature())
      {
        int dim = rowFeSpace->getMesh()->getDim();
        int degree = operat->getQuadratureDegree(2);
        Quadrature* quadrature = Quadrature::provideQuadrature(dim, degree);
        secondOrderAssembler->setQuadrature(quadrature);
      }
    }
    if (firstOrderAssemblerGrdPsi)
    {
      // create quadrature
      if (!firstOrderAssemblerGrdPsi->getQuadrature())
      {
        int dim = rowFeSpace->getMesh()->getDim();
        int degree = operat->getQuadratureDegree(1, GRD_PSI);
        Quadrature* quadrature = Quadrature::provideQuadrature(dim, degree);
        firstOrderAssemblerGrdPsi->setQuadrature(quadrature);
      }
    }
    if (firstOrderAssemblerGrdPhi)
    {
      // create quadrature
      if (!firstOrderAssemblerGrdPhi->getQuadrature())
      {
        int dim = rowFeSpace->getMesh()->getDim();
        int degree = operat->getQuadratureDegree(1, GRD_PHI);
        Quadrature* quadrature = Quadrature::provideQuadrature(dim, degree);
        firstOrderAssemblerGrdPhi->setQuadrature(quadrature);
      }
    }
    if (zeroOrderAssembler)
    {
      // create quadrature
      if (!zeroOrderAssembler->getQuadrature())
      {
        int dim = rowFeSpace->getMesh()->getDim();
        int degree = operat->getQuadratureDegree(0);
        Quadrature* quadrature = Quadrature::provideQuadrature(dim, degree);
        zeroOrderAssembler->setQuadrature(quadrature);
      }
    }
  }


  void Assembler::finishAssembling()
  {
    lastVecEl = NULL;
    lastMatEl = NULL;
  }


  OptimizedAssembler::OptimizedAssembler(Operator*  op,
                                         Quadrature* quad2,
                                         Quadrature* quad1GrdPsi,
                                         Quadrature* quad1GrdPhi,
                                         Quadrature* quad0,
                                         const FiniteElemSpace* rowFeSpace,
                                         const FiniteElemSpace* colFeSpace)
    : Assembler(op, rowFeSpace, colFeSpace)
  {
    bool opt = (rowFeSpace->getBasisFcts() == colFeSpace->getBasisFcts());

    // create sub assemblers
    secondOrderAssembler =
      SecondOrderAssembler::getSubAssembler(op, this, quad2, opt);
    firstOrderAssemblerGrdPsi =
      FirstOrderAssembler::getSubAssembler(op, this, quad1GrdPsi, GRD_PSI, opt);
    firstOrderAssemblerGrdPhi =
      FirstOrderAssembler::getSubAssembler(op, this, quad1GrdPhi, GRD_PHI, opt);
    zeroOrderAssembler =
      ZeroOrderAssembler::getSubAssembler(op, this, quad0, opt);

    checkQuadratures();
  }


  StandardAssembler::StandardAssembler(Operator* op,
                                       Quadrature* quad2,
                                       Quadrature* quad1GrdPsi,
                                       Quadrature* quad1GrdPhi,
                                       Quadrature* quad0,
                                       const FiniteElemSpace* rowFeSpace,
                                       const FiniteElemSpace* colFeSpace)
    : Assembler(op, rowFeSpace, colFeSpace)
  {
    remember = false;

    // create sub assemblers
    secondOrderAssembler =
      SecondOrderAssembler::getSubAssembler(op, this, quad2, false);
    firstOrderAssemblerGrdPsi =
      FirstOrderAssembler::getSubAssembler(op, this, quad1GrdPsi, GRD_PSI, false);
    firstOrderAssemblerGrdPhi =
      FirstOrderAssembler::getSubAssembler(op, this, quad1GrdPhi, GRD_PHI, false);
    zeroOrderAssembler =
      ZeroOrderAssembler::getSubAssembler(op, this, quad0, false);

    checkQuadratures();
  }

} // end namespace AMDiS
