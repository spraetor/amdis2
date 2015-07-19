/** \file CombinedPreconditioner.h */

#pragma once

#include "solver/BlockPreconditioner.h"

namespace AMDiS 
{
  /// preconditioner structure that combines various block-preconditioners
  template <class MatrixType>
  struct CombinedPreconditioner : public BlockPreconditioner<MatrixType>
  {
    typedef BlockPreconditioner<MatrixType>                         super;
    typedef ITL_PreconditionerBase<MatrixType, MTLTypes::MTLVector> precon_base;
    typedef CombinedPreconditioner<MatrixType>                      self;
 
    class Creator : public CreatorInterfaceName<precon_base>
    {
    public:
      // TODO: replace with initializer lists
      Creator(int l0) 
      { l.push_back(l0); };
      
      Creator(int l0, int l1) 
      { l.push_back(l0); l.push_back(l1); };
      
      Creator(int l0, int l1, int l2) 
      { l.push_back(l0); l.push_back(l1); l.push_back(l2); };
      
      Creator(const std::vector<int>& l) : l(l) {}
      
      virtual precon_base* create() override { 
	       return new self(l);
      }
      
    private:
      std::vector<int> l;
    };  
    
    /// Constructor
    CombinedPreconditioner(const std::vector<int>& parts_) 
      : super(), parts(parts_)
    { 
      precon.resize(parts.size(), &identity);
      subA.resize(parts.size());
    }
    
    /// Extract iranges from SolverMatrix to be used to extract sub-vectors and sub-matrices.
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*> >& A_, 
                      const MatrixType& fullMatrix_) override
    {
      super::A = &A_;
      super::fullMatrix = &fullMatrix_;
            
      BlockMapper mapper(A_);
      rows.resize(mapper.getNumComponents());
      int start = 0;
      int sum = 0;
      for (int p = 0; p < parts.size(); p++) {
      	sum += parts[p];
      	TEST_EXIT(sum <= mapper.getNumComponents())("range out of bound!");
      	mapper.setRow(sum);
      	int finish = mapper.row(0);
      	TEST_EXIT(finish > start)("ranges overlap!");
      	rows[p].set(start, finish);
      	start = finish;
      }
      
      const Matrix<DOFMatrix*>& mat = *A_.getOriginalMat();
      sum = 0;
      for (int p = 0; p < parts.size(); p++) {
      	Matrix<DOFMatrix*>* subMat = new Matrix<DOFMatrix*>(parts[p],parts[p]);
      	for (int i = 0; i < parts[p]; i++)
      	  for (int j = 0; j < parts[p]; j++)
      	    (*subMat)[i][j] = mat[i+sum][j+sum];
      	subA[p].setMatrix(*subMat);
      	precon[p]->init(subA[p], fullMatrix_);
      }
    }
    
    
    virtual void exit() override
    {
      for (int p = 0; p < parts.size(); p++) {
      	precon[p]->exit();
      	delete subA[p].getOriginalMat();
      }
    }
    
    
    /// Apply the preconditioners block-wise
    /**
     * solve Px = b, with
     * P = diag(P1, P2, ...)
     **/
    void solve(const MTLTypes::MTLVector& b, MTLTypes::MTLVector& x) const override
    { FUNCNAME("CombinedPreconditioner::solve()");

      x.change_dim(num_rows(b));
      
      for (size_t i = 0; i < precon.size(); i++) {
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
    
    // length of blocks assigned to a seperate preconditioners
    std::vector<int> parts;
    
    std::vector<precon_base*> precon;
    std::vector<SolverMatrix<Matrix<DOFMatrix*> > > subA;
    
    DOFMatrix::base_matrix_type dummy;
    ITL_Preconditioner<itl::pc::identity<MatrixType>, MatrixType, MTLTypes::MTLVector > identity;
  };
  
} // end namespace AMDiS
