/** \file CombinedPreconditioner.h */

#pragma once

#include <solver/BlockPreconditioner.h>

namespace AMDiS
{
  /// preconditioner structure that combines various block-preconditioners
  template <class MatrixType>
  struct CombinedPreconditioner : public BlockPreconditioner<MatrixType>
  {
    using Self        = CombinedPreconditioner; 
    using Super       = BlockPreconditioner<MatrixType>;
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
    CombinedPreconditioner(std::vector<int> const& parts_)
      : Super(), 
	parts(parts_)
    {
      precon.resize(parts.size(), &identity);
      subA.resize(parts.size());
    }

    /// Implementation of \ref ITL_PreconditionerBase::init
    /// Extract iranges from SolverMatrix to be used to extract sub-vectors and sub-matrices.
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A_,
                      MatrixType const& fullMatrix_) override
    {
      Super::A = &A_;
      Super::fullMatrix = &fullMatrix_;

      BlockMapper mapper(A_);
      rows.resize(mapper.getNumComponents());
      int start = 0;
      int sum = 0;
      for (size_t p = 0; p < parts.size(); p++)
      {
        sum += parts[p];
        TEST_EXIT(sum <= mapper.getNumComponents())("range out of bound!");
        mapper.setRow(sum);
        int finish = mapper.row(0);
        TEST_EXIT(finish > start)("ranges overlap!");
        rows[p].set(start, finish);
        start = finish;
      }

      Matrix<DOFMatrix*> const& mat = *A_.getOriginalMat();
      sum = 0;
      for (size_t p = 0; p < parts.size(); p++)
      {
        Matrix<DOFMatrix*>* subMat = new Matrix<DOFMatrix*>(parts[p],parts[p]);
        for (int i = 0; i < parts[p]; i++)
          for (int j = 0; j < parts[p]; j++)
            (*subMat)[i][j] = mat[i+sum][j+sum];
        subA[p].setMatrix(*subMat);
        precon[p]->init(subA[p], fullMatrix_);
      }
    }


    /// Implementation of \ref PreconditionerInterface::init
    virtual void exit() override
    {
      for (size_t p = 0; p < parts.size(); p++)
      {
        precon[p]->exit();
        delete subA[p].getOriginalMat();
      }
    }


    /// Implementation of \ref ITL_PreconditionerBase::solve
    /// Apply the preconditioners block-wise
    /**
     * solve Px = b, with
     * P = diag(P1, P2, ...)
     **/
    virtual void solve(MTLTypes::MTLVector const& b, MTLTypes::MTLVector& x) const override
    {
      FUNCNAME("CombinedPreconditioner::solve()");

      x.change_dim(num_rows(b));

      for (size_t i = 0; i < precon.size(); i++)
      {
        const MTLTypes::MTLVector b_i(b[rows[i]]);
        MTLTypes::MTLVector x_i(x[rows[i]]);
        precon[i]->solve(b_i, x_i);
      }
    }

    /// Sets a preconditioner for part i. If non is set, an identity preconditioner is used.
    void setPreconditioner(size_t i, precon_base& p)
    {
      precon[i] = &p;
    }

  protected:
    std::vector<mtl::irange> rows;

    /// lengths of blocks assigned to seperate preconditioners
    std::vector<int> parts;

    std::vector<precon_base*> precon;
    std::vector<SolverMatrix<Matrix<DOFMatrix*>>> subA;

    ITL_Preconditioner<itl::pc::identity<MatrixType>, MatrixType, MTLTypes::MTLVector> identity;
  };

} // end namespace AMDiS
