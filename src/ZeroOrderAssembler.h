/** \file ZeroOrderAssembler.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
#include "SubAssembler.h"

namespace AMDiS 
{
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

} // end namespace AMDiS
