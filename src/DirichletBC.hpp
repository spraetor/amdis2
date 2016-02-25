/** \file DirichletBC.h */

#pragma once

#include "AMDiS_fwd.hpp"
#include "BoundaryCondition.hpp"
#include "FixVec.hpp"
#include "OperatorTerm.hpp"
#include "DOFVector.hpp"

#include "expressions/TermGenerator.hpp"

namespace AMDiS
{
  namespace detail
  {
    class DirichletBC : public BoundaryCondition
    {
    public:

      /// Constructor.
      DirichletBC(BoundaryType type,
                  const FiniteElemSpace* rowFeSpace,
                  const FiniteElemSpace* colFeSpace,
                  bool apply)
        : BoundaryCondition(type, rowFeSpace, colFeSpace),
          applyBC(apply)
      { }

      /// Implementation of BoundaryCondition::fillBoundaryCondition().
      virtual void fillBoundaryCondition(DOFMatrix* matrix,
                                         ElInfo* elInfo,
                                         const DegreeOfFreedom* dofIndices,
                                         const BoundaryType* localBound,
                                         int nBasFcts) override;

      /// Implementation of BoundaryCondition::initVector()
      virtual void initVector(DOFVectorBase<double>* vec) override;

      /// Implementation of BoundaryCondition::boundResidual().
      virtual double boundResidual(ElInfo*, DOFMatrix*,
                                   const DOFVectorBase<double>*) override
      {
        return 0.0;
      }

      /// Because this is a Dirichlet boundary condition, always return true.
      virtual bool isDirichlet() const override
      {
        return true;
      }

      /// Returns \ref applyBC.
      virtual bool applyBoundaryCondition() const override
      {
        return applyBC;
      }

    protected:

      /// Defines, if the boundary condition must be applied to the matrix. See
      /// comment of \ref BoundaryCondition::applyBoundaryCondition.
      bool applyBC;
    };

  } // end namespace detail


  /**
  * \ingroup Assembler
  *
  * \brief
  * Sub class of BoundaryCondition. Implements Dirichlet boundary conditions.
  * A DOFVectors is set to a given value at a Dirichlet dof and in a DOFMatrix
  * the row corresponding to a Dirichlet dof is replaced by a row containing
  * only a 1.0 in the diagonal.
  */
  template <class Expr>
  class DirichletBC : public detail::DirichletBC
  {
    using Super = detail::DirichletBC;
    using TermType = ToTerm_t<Expr>;

  public:
    /// Constructor.
    template <class Expr_>
    DirichletBC(BoundaryType type,
                Expr_&& fct_,
                const FiniteElemSpace* rowFeSpace,
                const FiniteElemSpace* colFeSpace = NULL,
                bool apply = true)
      : Super(type, rowFeSpace, colFeSpace, apply),
        fct(toTerm(std::forward<Expr_>(fct_)))
    {}


    /// Implementation of BoundaryCondition::fillBoundaryCondition().
    virtual void fillBoundaryCondition(DOFVectorBase<double>* vector,
                                       ElInfo* elInfo,
                                       const DegreeOfFreedom* dofIndices,
                                       const BoundaryType* localBound,
                                       int nBasFcts) override
    {
      const BasisFunction* basFcts = rowFeSpace->getBasisFcts();
      // initialize expression on ElInfo
      fct.initElement(elInfo, NULL, NULL, basFcts);
      for (int i = 0; i < nBasFcts; i++)
        if (localBound[i] == boundaryType)
        {
          Value_t<TermType> value = fct[i];
          vector->setDirichletDofValue(dofIndices[i], value);
          (*vector)[dofIndices[i]] = value;
        }
    }

  protected:
    TermType fct;
  };

} // end namespace AMDiS
