/** \file FirstOrderAssembler.h */

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
   * SubAssembler for first order terms.
   */
  class FirstOrderAssembler : public SubAssembler
  {
  public:
    /// Creates and returns the FirstOrderAssembler for Operator op and
    /// the given assembler. If all terms are piecewise constant precalculated
    /// integrals can be used while assembling and the returned
    /// ZeroOrderAssembler is of type Pre0. Otherwise a Quad0 object will
    /// be returned.
    static FirstOrderAssembler* getSubAssembler(Operator* op,
        Assembler* assembler,
        Quadrature* quadrat,
        FirstOrderType type,
        bool optimized);

    /// Destructor.
    virtual ~FirstOrderAssembler() {}

  protected:
    /// Constructor.
    FirstOrderAssembler(Operator* op,
                        Assembler* assembler,
                        Quadrature* quadrat,
                        bool optimized,
                        FirstOrderType type);


  protected:
    /// Vector of DimMats for calculation in function calculateElementMatrix().
    std::vector<DenseVector<double>> Lb;

    /// List of all yet created optimized zero order assemblers for grdPsi.
    static ThreadPrivate<std::vector<SubAssembler*>> optimizedSubAssemblersGrdPsi;

    /// List of all yet created standard zero order assemblers for grdPsi.
    static ThreadPrivate<std::vector<SubAssembler*>> standardSubAssemblersGrdPsi;

    /// List of all yet created optimized zero order assemblers for grdPhi.
    static ThreadPrivate<std::vector<SubAssembler*>> optimizedSubAssemblersGrdPhi;

    /// List of all yet created standard zero order assemblers for grdPhi.
    static ThreadPrivate<std::vector<SubAssembler*>> standardSubAssemblersGrdPhi;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Standard first order assembler for grdPsi.
   */
  class Stand10 : public FirstOrderAssembler
  {
  public:
    /// Constructor
    Stand10(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override;

  protected:
    BasisFunction const *psi, *phi;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Standard first order assembler for grdPhi.
   */
  class Stand01 : public FirstOrderAssembler
  {
  public:
    /// Constructor.
    Stand01(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }

  protected:
    std::vector<DenseVector<double>> grdPhi;

    BasisFunction const *psi, *phi;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * First order assembler for grdPsi using fast quadratures.
   */
  class Quad10 : public FirstOrderAssembler
  {
  public:
    /// Constructor.
    Quad10(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * First order assembler for grdPhi using fast quadratures.
   */
  class Quad01 : public FirstOrderAssembler
  {
  public:
    /// Constructor.
    Quad01(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * First order assembler for grdPsi using precalculated integrals
   */
  class Pre10 : public FirstOrderAssembler
  {
  public:
    /// Constructor.
    Pre10(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override;

  protected:
    /// Integral of the product of the derivative of psi and phi.
    Q10PsiPhi const* q10;

    /// Integral of the derivative of psi.
    Q1Psi const* q1;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   *  First order assembler for grdPhi using precalculated integrals
   */
  class Pre01 : public FirstOrderAssembler
  {
  public:
    /// Constructor.
    Pre01(Operator* op, Assembler* assembler, Quadrature* quad);

  private:
    /// Implements SubAssembler::calculateElementMatrix().
    virtual void calculateElementMatrixImpl(ElInfo const*, ElementMatrix&) override;

    /// Implements SubAssembler::calculateElementVector().
    virtual void calculateElementVectorImpl(ElInfo const*, DenseVector<double>&) override
    {
      ERROR_EXIT("should not be called\n");
    }

  protected:
    /// Integral of the product of psi and the derivative of phi.
    Q01PsiPhi const* q01;

    /// Integral of the derivative of phi.
    Q1Psi const* q1;
  };

} // end namespace AMDiS
