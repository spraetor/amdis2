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



/** \file ZeroOrderAssembler.h */

#ifndef AMDIS_ZERO_ORDER_ASSEMBLER_H
#define AMDIS_ZERO_ORDER_ASSEMBLER_H

#include <vector>
#include "AMDiS_fwd.h"
#include "SubAssembler.h"

namespace AMDiS {

  /**
   * \ingroup Assembler
   * 
   * \brief
   * SubAssembler for zero order terms.
   */
  class ZeroOrderAssembler : public SubAssembler
  {
  public:
    /** \brief
     * Creates and returns the ZeroOrderAssembler for Operator op and
     * the given assembler. If all terms are piecewise constant precalculated 
     * integrals can be used while assembling and the returned 
     * ZeroOrderAssembler is of type PrecalcZOA. Otherwise a FastQuadZOA
     * object will be returned.
     */
    static ZeroOrderAssembler* getSubAssembler(Operator *op,
					       Assembler *assembler,
					       Quadrature *quadrat,
					       bool optimized);

    /// Destructor.
    virtual ~ZeroOrderAssembler() {}

  protected:
    /// Constructor.
    ZeroOrderAssembler(Operator *op, 
		       Assembler *assembler, 
		       Quadrature *quadrat,
		       bool optimized);

  protected:
    /// List of all yet created optimized SubAssembler objects.
    static ThreadPrivate<std::vector<SubAssembler*> > optimizedSubAssemblers;

    /// List of all yet created standard SubAssembler objects.
    static ThreadPrivate<std::vector<SubAssembler*> > standardSubAssemblers;
  };


  /**
   * \ingroup Assembler
   * 
   * \brief
   * Standard zero order assembler.
   */
  class StandardZOA :  public ZeroOrderAssembler 
  {
  public:
    /// Constructor.
    StandardZOA(Operator *op, Assembler *assembler, Quadrature *quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(const ElInfo *elInfo, ElementVector& vec) override;
  };


  /** 
   * \ingroup Assembler
   *
   * \brief
   * Zero order assembler using fast quadratures.
   */
  class FastQuadZOA :  public ZeroOrderAssembler
  {
  public:
    /// Constructor.
    FastQuadZOA(Operator *op, Assembler *assembler, Quadrature *quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(const ElInfo *elInfo, ElementVector& vec) override;

  protected:
    ElementVector c;
  };


  /**
   * \ingroup Assembler
   * 
   * \brief
   * Zero order assembler using precaculated integrals.
   */
  class PrecalcZOA : public ZeroOrderAssembler
  {
  public:
    /// Constructor.
    PrecalcZOA(Operator *op, Assembler *assembler, Quadrature *quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(const ElInfo *elInfo, ElementMatrix& mat) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(const ElInfo *elInfo, ElementVector& vec) override;

  protected:
    /// Integral of the product of psi and phi.
    const Q00PsiPhi *q00;

    /// Integral of psi.
    const Q0Psi *q0;
 
#if 0
    friend class ZeroOrderAssembler;
#endif
  };

}

#endif // AMDIS_ZERO_ORDER_ASSEMBLER_H
