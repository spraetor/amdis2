/** \file DOFMatrix.h */

#pragma once

#include <vector>
#include <set>
#include <memory>
#include <list>

#include "AMDiS_fwd.h"
#include "Flag.h"
#include "RCNeighbourList.h"
#include "DOFAdmin.h"
#include "DOFIndexed.h"
#include "Boundary.h"
#include "MTL4Types.h"

namespace AMDiS 
{
  /** \ingroup DOFAdministration
   * \brief
   * A DOFMatrix is a sparse matrix representation for matrices that work
   * on DOFVectors. The underlying matrix type 
   * is a CRS matrix of double.
   * TODO: make generic with value_type
   */
  class DOFMatrix
  {
  public:
    /// Type of scalars in the underlying matrix
    typedef MTLTypes::value_type value_type; 	// double
    typedef MTLTypes::size_type size_type;	// unsigned (long)

    /// Type of underlying matrix
    typedef mtl::matrix::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
    typedef mtl::compressed2D<value_type, para> base_matrix_type;

    /// Type of inserter for the base matrix;
    typedef mtl::matrix::inserter<base_matrix_type, mtl::operations::update_plus<value_type> >  inserter_type;

  private:
    /// Symbolic constant for an unused matrix entry
    static const int UNUSED_ENTRY = -1;
    
    /// Symbolic constant for an unused entry which is not followed by any
    /// used entry in this row
    static const int NO_MORE_ENTRIES = -2;
    
    /// Minimal slot-size used for the inserter. See \ref startInsertion for details.
    static const int MIN_NNZ_PER_ROW = 20;

  public:    
    DOFMatrix();

    /// Constructs a DOFMatrix with name n and the given row and column FeSpaces.
    DOFMatrix(const FiniteElemSpace* rowFeSpace, 
      	      const FiniteElemSpace* colFeSpace,
      	      std::string n = "");

    /// Copy-Constructor
    DOFMatrix(const DOFMatrix& rhs);

    /// Destructor
    virtual ~DOFMatrix();
  
    /// Assignment operator.
    DOFMatrix& operator=(const DOFMatrix& rhs);

    void copy(const DOFMatrix& rhs);

    /// Access underlying matrix directly
    base_matrix_type& getBaseMatrix()
    {
      return matrix;
    }

    /// Access underlying matrix directly (const)
    const base_matrix_type& getBaseMatrix() const
    {
      return matrix;
    }

    /// Returns \ref coupleMatrix.
    bool isCoupleMatrix() const
    { 
      return coupleMatrix; 
    }

    /// Returns \ref coupleMatrix.
    void setCoupleMatrix(bool c) 
    { 
      coupleMatrix = c; 
    }

    /// Adds an operator to the DOFMatrix. A factor, that is multipled to the 
    /// operator, and a multilier factor for the estimator may be also given.
    void addOperator(Operator *op, 
            		     double* factor = NULL, 
            		     double* estFactor = NULL);

    ///
    void clearOperators();

    std::vector<double*>::iterator getOperatorFactorBegin()
    {
      return operatorFactor.begin();
    }

    std::vector<double*>::iterator getOperatorFactorEnd() 
    {
      return operatorFactor.end();
    }

    std::vector<double*>::iterator getOperatorEstFactorBegin()
    {
      return operatorEstFactor.begin();
    }

    std::vector<double*>::iterator getOperatorEstFactorEnd()
    {
      return operatorEstFactor.end();
    }

    std::vector<Operator*>::iterator getOperatorsBegin()
    {
      return operators.begin();
    }

    std::vector<Operator*>::iterator getOperatorsEnd()
    {
      return operators.end();
    }

    Flag getAssembleFlag();

    /** \brief
     * Updates the matrix matrix by traversing the underlying mesh and assembling
     * the element contributions into the matrix. Information about the 
     * computation of element matrices and connection of local and global DOFs is
     * stored in minfo; the flags for the mesh traversal are stored at 
     * minfo->fill_flags which specifies the elements to be visited and 
     * information that should be present on the elements for the calculation of 
     * the element matrices and boundary information (if minfo->boundBas is not
     * NULL). On the elements, information about the row DOFs is accessed by 
     * minfo->rowBas using info->row_admin; this vector is also used for the 
     * column DOFs if minfo->nCol is less or equal zero, or minfo->col_admin or 
     * minfo->colBas is a NULL pointer; if row and column DOFs are the same, the 
     * boundary type of the DOFs is accessed by minfo->boundBas if 
     * minfo->boundBas is not a NULL pointer; then the element matrix is 
     * computed by minfo->fillElementMatrix(el info, minfo->myFill); these 
     * contributions, multiplied by minfo->factor, are eventually added to matrix
     * by a call of addElementMatrix() with all information about row and column 
     * DOFs, the element matrix, and boundary types, if available;
     * updateMatrix() only adds element contributions; this makes several calls 
     * for the assemblage of one matrix possible; before the first call, the 
     * matrix should be cleared by calling clear dof matrix().
     */
  
    void assemble(double factor, ElInfo *elInfo, const BoundaryType *bound);

    void assemble(double factor, ElInfo *elInfo, const BoundaryType *bound,
	                Operator *op);

    /// Adds an element matrix to \ref matrix
    void addElementMatrix(const ElementMatrix& elMat, 
                  			  const BoundaryType *bound,
                  			  ElInfo* rowElInfo,
                  			  ElInfo* colElInfo);

    ///
    void assembleOperator(Operator &op);

    /// That function must be called after the matrix assembling has been 
    /// finished. This makes it possible to start some cleanup or matrix 
    /// data compressing procedures.
    void finishAssembling();

    /// Enable insertion for assembly. You can optionally give an upper limit for
    /// entries per row (per column for CCS matrices).  Choosing this parameter
    /// too small can induce perceivable overhead for compressed matrices. Thus,
    /// it's better to choose a bit too large than too small.
    void startInsertion(int nnz_per_row = MIN_NNZ_PER_ROW);

    /// Finishes insertion. For compressed matrix types, this is where the
    /// compression happens.
    void finishInsertion()
    {
      FUNCNAME("DOFMatrix::finishInsertion()");

      TEST_EXIT(inserter)("Inserter wasn't used or is already finished.\n");
      
      delete inserter;
      inserter= 0;
    }

    /// Returns const \ref rowFeSpace
    const FiniteElemSpace* getRowFeSpace() const 
    { 
      return rowFeSpace; 
    }

    /// Returns const \ref colFeSpace
    const FiniteElemSpace* getColFeSpace() const 
    { 
      return colFeSpace; 
    }

    /// Returns const \ref rowFeSpace
    const FiniteElemSpace* getFeSpace() const 
    { 
      return rowFeSpace; 
    }

    /// Returns number of rows (\ref matrix.size())
    int getSize() const 
    { 
      return num_rows(matrix);
    }

    /// Returns the number of used rows (equal to number of used DOFs in
    /// the row FE space).
    int getUsedSize() const;

    std::set<DegreeOfFreedom>& getDirichletRows()
    {
      return dirichletDofs;
    }

    /// Returns \ref name
    std::string getName() const 
    { 
      return name; 
    }

    void clearDirichletRows();

    /// Prints \ref matrix to stdout
    void print() const;

    /// Removes all matrix entries
    void clear()
    {
	set_to_zero(matrix);
    }

    std::vector<Operator*>& getOperators() 
    { 
      return operators; 
    }
    
    std::vector<double*>& getOperatorFactor() 
    { 
      return operatorFactor; 
    }

    std::vector<double*>& getOperatorEstFactor() 
    { 
      return operatorEstFactor; 
    }

    BoundaryManager* getBoundaryManager() const 
    { 
      return boundaryManager; 
    }

    void setBoundaryManager(BoundaryManager *bm) 
    {
      boundaryManager = bm;
    }

    /// Calculates the average of non zero entries per row in matrix.
    void calculateNnz()
    {
      nnzPerRow = 0;

      if (num_rows(matrix) != 0)
        nnzPerRow = int(double(matrix.nnz()) / num_rows(matrix) * 1.2); 
      if (nnzPerRow < MIN_NNZ_PER_ROW) 
        nnzPerRow = MIN_NNZ_PER_ROW;
    }

    /// Returns \ref nnzPerRow.
    int getNnz() const
    {
      return nnzPerRow;
    }
    
    std::vector<DegreeOfFreedom>& getRowIndices()
    {
      return rowIndices;
    }
    
    std::vector<DegreeOfFreedom>& getColIndices()
    {
      return colIndices;
    }
    
    ///
    int memsize() const;

  protected:
    /// Pointer to a FiniteElemSpace with information about corresponding row DOFs
    /// and basis functions
    const FiniteElemSpace *rowFeSpace;

    /// Pointer to a FiniteElemSpace with information about corresponding 
    /// column DOFs and basis functions
    const FiniteElemSpace *colFeSpace;

    /// Name of the DOFMatrix
    std::string name;

    /// Sparse matrix, type is a template parameter by 
    /// default compressed2D<double>
    base_matrix_type matrix;

    /// Pointers to all operators of the equation systems. Are used in the
    /// assembling process.
    std::vector<Operator*> operators;
    
    /// Defines for each operator a factor which is used to scal the element
    /// matrix after the assembling process of the operator.
    std::vector<double*> operatorFactor;

    ///
    std::vector<double*> operatorEstFactor;

    ///
    BoundaryManager *boundaryManager;

    /// If false, the matrix is a diagonal matrix within a matrix of DOF matrices.
    /// Otherwise the value is true, and the matrix is an off-diagonal matrix.
    bool coupleMatrix;

    /// Temporary variable used in assemble()
    ElementMatrix elementMatrix;

    /// Number of basis functions in the row fe space.
    int nRow;

    /// Number of basis functions in the col fe space.
    int nCol;

    /// Maps local row indices of an element to global matrix indices.
    std::vector<DegreeOfFreedom> rowIndices;

    /// Maps local col indices of an element to global matrix indices.
    std::vector<DegreeOfFreedom> colIndices;

    /// A set of row indices. When assembling the DOF matrix, all rows, that
    /// correspond to a dof at a dirichlet boundary, are ignored and the row is
    /// left blank. After assembling, the diagonal element of the matrix must be
    /// set to 1. The indices of all rows, where this should be done, are stored
    /// in this set.
    std::set<DegreeOfFreedom> dirichletDofs;

    /// Number of non zero entries per row (average). For instationary problems this
    /// information may be used in the next timestep to accelerate insertion of
    /// elemnts during assembling.
    int nnzPerRow;

    /// Inserter object: implemented as pointer, allocated and deallocated as needed
    inserter_type *inserter;
      
    friend class DOFAdmin;
    friend class DOFVector<double>;
    friend class DOFVector<unsigned char>;
    friend class DOFVector<int>;
    friend class DOFVector<WorldVector<double> >;
  };


  namespace traits
  {
    template <>
    struct category< AMDiS::DOFMatrix > 
    {
      typedef tag::matrix  tag;
      typedef MTLTypes::value_type  value_type;
      typedef MTLTypes::size_type   size_type;
    };
  }
  
  inline size_t size(DOFMatrix const& A) { return size(A.getBaseMatrix()); }
  
  inline size_t num_rows(DOFMatrix const& A) { return num_rows(A.getBaseMatrix()); }
  
  inline size_t num_cols(DOFMatrix const& A) { return num_cols(A.getBaseMatrix()); } 
  
  inline void set_to_zero(DOFMatrix& A) { A.clear(); }
  
} // end namespace AMDiS
