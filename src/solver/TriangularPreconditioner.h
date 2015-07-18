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
 

/** \file TriangularPreconditioner.h */

#ifndef AMDIS_TRIANGULAR_PRECONDITIONER_H
#define AMDIS_TRIANGULAR_PRECONDITIONER_H

#include "solver/CombinedPreconditioner.h"

namespace AMDiS {

  /// preconditioner structure that combines various block-preconditioners
  template<typename MatrixType>
  struct TriangularPreconditioner : CombinedPreconditioner<MatrixType>
  {
    typedef CombinedPreconditioner<MatrixType>                         super;
    typedef ITL_PreconditionerBase<MatrixType, MTLTypes::MTLVector>    precon_base;
    typedef TriangularPreconditioner<MatrixType>                       self;
 
    class Creator : public CreatorInterfaceName<precon_base>
    {
    public:
      Creator(int l0) 
      { l.push_back(l0); };
      
      Creator(int l0, int l1) 
      { l.push_back(l0); l.push_back(l1); };
      
      Creator(int l0, int l1, int l2) 
      { l.push_back(l0); l.push_back(l1); l.push_back(l2); };
      
      Creator(const std::vector<int>& l) : l(l) {}
      
      virtual ~Creator() {}
      
      precon_base* create() { 
	return new self(l);
      }
    private:
      std::vector<int> l;
    };  
    
    /// Constructor
    TriangularPreconditioner(const std::vector<int>& parts_) 
      : super(parts_)
    { }
    
    virtual ~TriangularPreconditioner() {}
    
    /// Apply the preconditioners block-wise
    /**
     * solve Px = b, with
     * P = diag(P1, P2, ...)
     **/
    void solve(const MTLTypes::MTLVector& b, MTLTypes::MTLVector& x) const override
    { FUNCNAME("CombinedPreconditioner::solve()");
      
      x.change_dim(num_rows(b));
      y.change_dim(num_rows(b));
      c.change_dim(num_rows(b));
      set_to_zero(c);
      
      for (size_t i = super::precon.size()-1; i >= 0; i--) {
	MTLTypes::MTLVector y_i(b[super::rows[i]]);
	if (i < super::precon.size()-1)
	  y_i -= y[super::rows[i]];
	
	MTLTypes::MTLVector x_i(x[super::rows[i]]);
	super::precon[i]->solve(y_i, x_i);
	
	if (i > 0) {
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

#endif // AMDIS_COMBINED_PRECONDITIONER_H
