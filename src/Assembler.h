/** \file Assembler.h */

/**
 * \defgroup Assembler Assembler module
 *
 * \brief
 * Contains the operator and assembler classes:
 * @{ <img src="assembler.png"> @}
 */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
#include "FixVec.h"
#include "ZeroOrderAssembler.h"
#include "FirstOrderAssembler.h"
#include "SecondOrderAssembler.h"
#include <traits/not_null.hpp>

namespace AMDiS
{
  /**
   * \ingroup Assembler
   *
   * \brief
   * Assembles element matrices and vectors for a given Operator. Uses
   * one SubAssembler for all second order terms of the Operator, one for all
   * first order terms, and one for all zero order terms.
   */
  class Assembler
  {
  public:
    /// Constructor
    Assembler(not_null<Operator*> op,
              not_null<const FiniteElemSpace*> rowFeSpace,
		       const FiniteElemSpace*  colFeSpace = NULL);

    /// Destructor
    virtual ~Assembler() {};

    /// Assembles the element matrix for the given ElInfo
    void calculateElementMatrix(const ElInfo* elInfo,
                                ElementMatrix& userMat,
                                double factor = 1.0);

    /// Assembles the element vector for the given ElInfo
    void calculateElementVector(const ElInfo* elInfo,
                                DenseVector<double>& userVec,
                                double factor = 1.0);

    /// Returns \ref rowFeSpace.
    const FiniteElemSpace* getRowFeSpace() const
    {
      return rowFeSpace;
    }

    /// Returns \ref colFeSpace.
    const FiniteElemSpace* getColFeSpace() const
    {
      return colFeSpace;
    }

    /// Returns \ref nRow.
    int getNRow() const
    {
      return nRow;
    }

    /// Returns \ref nCol.
    int getNCol() const
    {
      return nCol;
    }

    /// Sets \ref rememberElMat.
    void rememberElementMatrix(bool rem)
    {
      rememberElMat = rem;
    }

    /// Sets \ref rememberElVec.
    void rememberElementVector(bool rem)
    {
      rememberElVec = rem;
    }

    /// Returns \ref zeroOrderAssembler.
    ZeroOrderAssembler* getZeroOrderAssembler() const
    {
      return zeroOrderAssembler;
    }

    /** \brief
     * Returns \ref firstOrderAssemblerGrdPsi or \ref firstOrderAssemblerGrdPhi
     * depending on type.
     */
    FirstOrderAssembler* getFirstOrderAssembler(FirstOrderType type = GRD_PSI) const
    {
      return (type == GRD_PSI)
             ? firstOrderAssemblerGrdPsi
             : firstOrderAssemblerGrdPhi;
    }

    /// Returns \ref secondOrderAssembler.
    SecondOrderAssembler* getSecondOrderAssembler() const
    {
      return secondOrderAssembler;
    }

    /// Returns \ref operat;
    Operator* getOperator() const
    {
      return operat;
    }

    /// Initialisation for the given ElInfo. The call is deligated to the sub assemblers.
    void initElement(const ElInfo* smallElInfo,
                     const ElInfo* largeElInfo = NULL,
                     Quadrature* quad = NULL);

    /// Sets quadratures of all sub assemblers.
    void setQuadratures(Quadrature* quad2,
                        Quadrature* quad1GrdPsi,
                        Quadrature* quad1GrdPhi,
                        Quadrature* quad0)
    {
      if (secondOrderAssembler)
      {
        TEST_EXIT(!secondOrderAssembler->getQuadrature())
        ("quadrature already existing\n");
        secondOrderAssembler->setQuadrature(quad2);
      }
      if (firstOrderAssemblerGrdPsi)
      {
        TEST_EXIT(!firstOrderAssemblerGrdPsi->getQuadrature())
        ("quadrature already existing\n");
        firstOrderAssemblerGrdPsi->setQuadrature(quad1GrdPsi);
      }
      if (firstOrderAssemblerGrdPhi)
      {
        TEST_EXIT(!firstOrderAssemblerGrdPhi->getQuadrature())
        ("quadrature already existing\n");
        firstOrderAssemblerGrdPhi->setQuadrature(quad1GrdPhi);
      }
      if (zeroOrderAssembler)
      {
        TEST_EXIT(!zeroOrderAssembler->getQuadrature())
        ("quadrature already existing\n");
        zeroOrderAssembler->setQuadrature(quad0);
      }
    }

    /// That function must be called after one assembling cycle has been finished.
    void finishAssembling();

  protected:
    /** \brief
     * Vector assembling by element matrix-vector multiplication.
     * Usefull if an element matrix was already calculated.
     */
    void matVecAssemble(const ElInfo* elInfo, DenseVector<double>& vec);

    /** \brief
     * Checks whether quadratures for subassemblers are already set.
     * If not they will be created.
     */
    void checkQuadratures();

  protected:
    /// Operator this Assembler belongs to.
    Operator* operat;

    /// Row FiniteElemSpace.
    const FiniteElemSpace* rowFeSpace;

    /// Column FiniteElemSpace.
    const FiniteElemSpace* colFeSpace;

    /// Number of rows.
    int nRow;

    /// Number of columns.
    int nCol;

    /// SubAssembler for the second order terms
    SecondOrderAssembler* secondOrderAssembler;

    /// SubAssembler for the first order terms (grdPsi)
    FirstOrderAssembler* firstOrderAssemblerGrdPsi;

    /// SubAssembler for the first order terms (grdPhi)
    FirstOrderAssembler* firstOrderAssemblerGrdPhi;

    /// SubAssembler for the zero order terms
    ZeroOrderAssembler* zeroOrderAssembler;

    ///
    bool remember;

    /// Determines whether the element matrix should be stored locally.
    bool rememberElMat;

    /// Determines whether the element vector should be stored locally.
    bool rememberElVec;

    /// Locally stored element matrix
    ElementMatrix elementMatrix;

    /// Locally stored element vector
    DenseVector<double> elementVector;

    ///
    ElementMatrix tmpMat;

    /// Used to check whether \ref initElement() must be called, because
    /// a new Element is visited.
    Element* lastMatEl;

    /// Used to check whether \ref initElement() must be called, because
    /// a new Element is visited.
    Element* lastVecEl;

    /// Used to check for new traverse.
    int lastTraverseId;
  };

  /**
   * \ingroup Assembler
   *
   * \brief
   * Assembler using non optimized sub assemblers.
   */
  class StandardAssembler : public Assembler
  {
  public:
    /// Constructor
    StandardAssembler(Operator* op,
                      Quadrature* quad2,
                      Quadrature* quad1GrdPsi,
                      Quadrature* quad1GrdPhi,
                      Quadrature* quad0,
                      const FiniteElemSpace* rowFeSpace,
                      const FiniteElemSpace* colFeSpace = NULL);
  };

  /**
   * \ingroup Assembler
   *
   * \brief
   * Assembler using optimized sub assemblers.
   */
  class OptimizedAssembler : public Assembler
  {
  public:
    /// Constructor
    OptimizedAssembler(Operator* op,
                       Quadrature* quad2,
                       Quadrature* quad1GrdPsi,
                       Quadrature* quad1GrdPhi,
                       Quadrature* quad0,
                       const FiniteElemSpace* rowFeSpace,
                       const FiniteElemSpace* colFeSpace = NULL);
  };

} // end namespace AMDiS
