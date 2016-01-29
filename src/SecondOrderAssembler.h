/** \file SecondOrderAssembler.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
#include "SubAssembler.h"

namespace AMDiS
{
  // forward declaration
  class Q11PsiPhi;

  /**
   * \ingroup Assembler
   *
   * \brief
   * SubAssembler for second order terms.
   */
  class SecondOrderAssembler : public SubAssembler
  {
  public:
    /** \brief
     * Creates and returns the SecondOrderAssembler for Operator op and
     * the given assembler. If all terms are piecewise constant precalculated
     * integrals can be used while assembling and the returned
     * ZeroOrderAssembler is of type Pre0. Otherwise a Quad0 object will
     * be returned.
     */
    static SecondOrderAssembler* getSubAssembler(Operator* op,
						 Assembler* assembler,
						 Quadrature* quadrat,
						 bool optimized);

  protected:
    /// Constructor.
    SecondOrderAssembler(Operator* op,
                         Assembler* assembler,
                         Quadrature* quadrat,
                         bool optimized);

  protected:
    /// List of all yet created optimized second order assemblers.
    static ThreadPrivate<std::vector<SubAssembler*>> optimizedSubAssemblers;

    /// List of all yet created standard second order assemblers.
    static ThreadPrivate<std::vector<SubAssembler*>> standardSubAssemblers;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Standard second order assembler
   */
  class Stand2 : public SecondOrderAssembler
  {
  public:
    /// Constructor.
    Stand2(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, 
					    ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, 
					    DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Second order assembler using fast quadratures.
   */
  class Quad2 : public SecondOrderAssembler
  {
  public:
    /// Constructor.
    Quad2(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, 
					    ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, 
					    DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }

  protected:
    DenseVector<double> dimVec;

    std::vector<DenseMatrix<double>> LALt;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Second order assembler using predefined integrals.
   */
  class Pre2 : public SecondOrderAssembler
  {
  public:
    /// Constructor.
    Pre2(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, 
					    ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, 
					    DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }

  protected:
    /// Integral of the product of the derivative of psi and the derivative of phi.
    const Q11PsiPhi* q11;

    std::vector<DenseMatrix<double>> LALt;
  };

} // end namespace AMDiS
