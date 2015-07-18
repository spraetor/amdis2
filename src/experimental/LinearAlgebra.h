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

 /// This file is work in progress and may replace the assembling now part of
 /// DOFMatrix directly. The idea is to assemble matrices to a data-structure
 /// provided by a linear algebra backend. You have to specify how an element-
 /// matrix can be added to the global matrix, how a Dirichlet-row is inserted
 /// and how the matrix is initialized/finalized.


/** \file LinearAlgebra.h */

#ifndef AMDIS_LINEAR_ALGEBRA_H
#define AMDIS_LINEAR_ALGEBRA_H

#include "utility/calculate_nnz.hpp"

namespace AMDiS {
    
  namespace tag 
  {
    struct distributed {};
    struct non_distributed {};
    
    // backends
    struct mtl {};
    struct pmtl {};
    struct petsc {};
    struct p_petsc {};
  }
  
  namespace traits
  {
    template< class T >
    struct distributed_tag : boost::mpl::if_< 
      mtl::traits::is_distributed< T >,
      tag::distributed,
      tag::non_distributed
    > {};
  }
  
  namespace detail
  {
    template< class Dummy = void >
    struct DefaultParameters {};
  }
  
  /** fill for each backend a linear algebra class that
   *  contains the following typedefs and functions
   **/
  template< class Backend = tag::mtl, class Value = double, class Parameters = detail::DefaultParameters<> >
  struct LinearAlgebra 
  {
  public: // typedefs
    
    typedef Value         ValueType;
    typedef unsigned int  SizeType;
    
    typedef mtl::matrix::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, SizeType> MatrixPara;
    typedef mtl::matrix::compressed2D<ValueType, MatrixPara>   MatrixBase;
    typedef mtl::vector::dense_vector<ValueType>               VectorBase;
    
    typedef mtl::matrix::inserter<MatrixBase, mtl::operations::update_plus<ValueType> >  Inserter;
    typedef MPI::Intracomm   Communicator;
    
    static const int minSlotSize = 20;
    
  public: // methods
    
    /// constructor
    LinearAlgebra(Communicator& comm)
      : comm(comm), ins(NULL)
    { }
    
    
    /// prepare matrix for insertion
    template< class M >
    void init_matrix(MatrixBase& matrix, MapperBase<M>& mapper, bool add_to_matrix = false)
    {
      int slot_size = estimate_slot_size(matrix);
      prepare_matrix(matrix, mapper, traits::distributed_tag<MatrixBase>::type() );
      
      if (!add_to_matrix)
	set_to_zero(matrix);
      
      // prepare the inserter
      if (ins) {
	delete ins;
	ins = NULL; 
      }
      ins = new Inserter(matrix, slot_size);
    }
    
    
    /// finish insertion
    template< class M >
    void exit_matrix(MatrixBase& matrix, MapperBase<M>& mapper)
    {
      // destroy inserter, to start insertion
      assert( ins );
      delete ins;
      inserter = NULL;
    }
    
    
    /// insert values of element_matrix into matrix
    template< class ElementMatrix, class VectorType >
    void add_element_matrix(MatrixBase& matrix, const ElementMatrix& elMatrix, 
			    const VectorType& rowIndices, const VectorType& colIndices)
    {
      assert( ins );
      *ins << mtl::matrix::element_matrix(elMatrix, rowIndices, colIndices);
    }
    
  protected: // methods
    
    /// resize matrix to appropriate size
    template< class M >
    void prepare_matrix(MatrixBase& matrix, MapperBase<M>& mapper, tag::non_distributed)
    {
      matrix.change_dim(mapper.getNumRows(), mapper.getNumCols());
    }
    
    
    /// estimate the number of nonzeros per row
    int estimate_slot_size(const MatrixBase& matrix)
    {
      int nnzPerRow = 0;
      if (num_rows(matrix) != 0)
	nnzPerRow = static_cast<int>((matrix.nnz() * 1.2) / num_rows(matrix)); 
      return std::max( nnzPerRow, minSlotSize );
    }
    
  private: // member variables
    Communicator& comm;
    Inserter* ins;
  };

#ifdef HAVE_PARALLEL_MTL4
  template< class Value, class Parameters >
  struct LinearAlgebra< tag::pmtl, Value, Parameters > : public LinearAlgebra< tag::mtl, Value, Parameters >
  {
    typedef LinearAlgebra< tag::mtl, Value, Parameters > super;
    
    /// constructor
    LinearAlgebra(Communicator& comm)
      : super(comm)
    { }
    
  protected:
    // eventuell als freie funktion schreiben
    template< class M >
    void prepare_matrix(MatrixBase& matrix, MapperBase<M>& mapper, tag::distributed)
    {
      mtl::par::block_distribution dist(mapper.getNumRows());
      dist.setup_from_local_size(mapper.getMap().getLocalDofs());
      matrix.change_dim(0, 0);
      matrix.init_distribution(dist, dist, mapper.getNumRows(), mapper.getNumRows());
    }
  };
#endif
  
  
#ifdef HAVE_SEQ_PETSC
  template< class Value, class Parameters >
  struct LinearAlgebra< tag::petsc, Value, Parameters >
  {
  public: // typedefs
    
    typedef PetscValue  ValueType;
    typedef PetscInt    SizeType;
    
    typedef PetscMatrix   MatrixBase;
    typedef PetscVector   VectorBase;
    
    typedef MPI::Intracomm   Communicator;
    
  public: // methods
    
    /// constructor
    LinearAlgebra(Communicator& comm)
      : comm(comm)
    { }
    
    
    /// prepare matrix for insertion
    template< class M >
    void init_matrix(MatrixBase& mat, MapperBase<M>& mapper, bool add_to_matrix = false)
    {
      std::vector<size_t> nnz_per_row(mapper.getNumRows());
      calculate_nnz(dof_matrix /* ?? */, nnz_per_row); // siehe auch MatrixNnzStructure fuer parallele NNZ

      // eventuell nur max_nnz angeben
      MatCreateSeqAIJ(comm, mapper.getNumRows(), mapper.getNumCols(), 0, &(nnz_per_row[0]), &mat.matrix);
    }
    
    
    /// finish insertion
    template< class M >
    void exit_matrix(MatrixBase& mat, MapperBase<M>& mapper)
    {
      MatAssemblyBegin(mat.matrix, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat.matrix, MAT_FINAL_ASSEMBLY);
    }
    
    
    /// insert values of element_matrix into matrix
    template< class ElementMatrix, class VectorType >
    void add_element_matrix(MatrixBase& mat, const ElementMatrix& elMatrix, 
			    const VectorType& rowIndices, const VectorType& colIndices)
    {
      // delete dirichlet-dofs by changing the sign of these indices
      MatSetValues(mat.matrix, rowIndices.size(), &rowIndices[0], colIndices.size(), &colIndices[0], elMatrix.address_data(), ADD_VALUES);
    }
    
  private: // member variables
    Communicator& comm;
  };
#endif
}

#endif  // AMDIS_DOFMATRIX_H
