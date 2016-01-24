/** \file TriangularPreconditioner.h */

#pragma once

#include "solver/CombinedPreconditioner.h"

namespace AMDiS
{
  /// preconditioner structure that combines various block-preconditioners
  template <class MatrixType>
  struct TriangularPreconditioner : public CombinedPreconditioner<MatrixType>
  {
    using Super = CombinedPreconditioner<MatrixType>;
    using Self  = TriangularPreconditioner;
    
    using precon_base = ITL_PreconditionerBase<MatrixType, MTLTypes::MTLVector>;

    struct Creator : public CreatorInterfaceName<precon_base>
    {
      /// Constructor for a list of integer arguments, representing the sizes
      /// of the sub-blocks
      template <class... Ints,
	class = Requires_t< and_< std::is_convertible<Ints, int>... > >>
      Creator(Ints... ls) : l{int(ls)...} {}

      /// Constructor for a vector of integers, representing the sizes
      /// of the sub-blocks
      Creator(std::vector<int> const& l) : l(l) {}

      virtual precon_base* create() override
      {
        return new Self(l);
      }

    private:
      std::vector<int> l;
    };

    /// Constructor
    TriangularPreconditioner(std::vector<int> const& parts_)
      : Super(parts_)
    {}

    /// Apply the preconditioners block-wise
    /**
     * solve Px = b, with
     * P = diag(P1, P2, ...)
     **/
    virtual void solve(MTLTypes::MTLVector const& b,
                       MTLTypes::MTLVector& x) const override
    {
      FUNCNAME("CombinedPreconditioner::solve()");

      x.change_dim(num_rows(b));
      y.change_dim(num_rows(b));
      c.change_dim(num_rows(b));
      set_to_zero(c);

      for (int i = int{super::precon.size()}-1; i >= 0; --i)
      {
        MTLTypes::MTLVector y_i(b[super::rows[i]]);
        if (i < super::precon.size()-1)
          y_i -= y[super::rows[i]];

        MTLTypes::MTLVector x_i(x[super::rows[i]]);
        super::precon[i]->solve(y_i, x_i);

        if (i > 0)
        {
          c[super::rows[i]] = x_i;
          y = (*super::fullMatrix) * c;
        }
      }
    }

  private:
    mutable MTLTypes::MTLVector y;
    mutable MTLTypes::MTLVector c;
  };

} // end namespace AMDiS
