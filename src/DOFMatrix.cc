#include <algorithm>

#include <boost/numeric/mtl/mtl.hpp>

#include "DOFMatrix.h"
#include "QPsiPhi.h"
#include "BasisFunction.h"
#include "Boundary.h"
#include "DOFAdmin.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "DOFVector.h"
#include "Operator.h"
#include "BoundaryCondition.h"
#include "BoundaryManager.h"
#include "Assembler.h"

namespace AMDiS 
{
  using namespace mtl;

  DOFMatrix::DOFMatrix()
    : rowFeSpace(NULL),
      colFeSpace(NULL),
      elementMatrix(3, 3),
      nRow(0),
      nCol(0),
      nnzPerRow(0),
      inserter(NULL)
  {}


  DOFMatrix::DOFMatrix(const FiniteElemSpace* rowSpace,
		       const FiniteElemSpace* colSpace,
		       std::string n)
    : rowFeSpace(rowSpace),
      colFeSpace(colSpace),
      name(n), 
      coupleMatrix(false),
      nnzPerRow(0),
      inserter(NULL)
  {
    FUNCNAME("DOFMatrix::DOFMatrix()");

    TEST_EXIT(rowFeSpace)("No fe space for row!\n");

    if (!colFeSpace)
      colFeSpace = rowFeSpace;
    boundaryManager = new BoundaryManager(rowFeSpace);

    nRow = rowFeSpace->getBasisFcts()->getNumber();
    nCol = colFeSpace->getBasisFcts()->getNumber();
    elementMatrix.change_dim(nRow, nCol);
    rowIndices.resize(nRow);
    colIndices.resize(nCol);
  }


  DOFMatrix::DOFMatrix(const DOFMatrix& rhs)
    : name(rhs.name + "copy")
  {
    FUNCNAME("DOFMatrix::DOFMatrix()");

    *this = rhs;
    TEST_EXIT(rhs.inserter == 0)("Cannot copy during insertion!\n");
    inserter = 0;
  }


  DOFMatrix::~DOFMatrix()
  {
    if (boundaryManager) 
      delete boundaryManager;
    if (inserter) 
      delete inserter;
  }


  void DOFMatrix::print() const
  {
    if (inserter) 
      inserter->print();
  }


  DOFMatrix& DOFMatrix::operator=(const DOFMatrix& rhs)
  {
    rowFeSpace = rhs.rowFeSpace;
    colFeSpace = rhs.colFeSpace;
    operators = rhs.operators;
    operatorFactor = rhs.operatorFactor;
    coupleMatrix = rhs.coupleMatrix;

    /// The matrix values may only be copyed, if there is no active insertion.
    if (rhs.inserter == 0 && inserter == 0)
      matrix = rhs.matrix;

    if (rhs.boundaryManager)
      boundaryManager = new BoundaryManager(*rhs.boundaryManager);
    else
      boundaryManager = NULL;
    
    nRow = rhs.nRow;
    nCol = rhs.nCol;
    elementMatrix.change_dim(nRow, nCol);

    return *this;
  }


  void DOFMatrix::addElementMatrix(const ElementMatrix& elMat, 
				   const BoundaryType *bound,
				   ElInfo* rowElInfo,
				   ElInfo* colElInfo)
  {
    FUNCNAME_DBG("DOFMatrix::addElementMatrix()");

    TEST_EXIT_DBG(inserter)("DOFMatrix is not in insertion mode\n");
    TEST_EXIT_DBG(rowFeSpace)("Have now rowFeSpace!\n");

    inserter_type &ins= *inserter;
 
    // === Get indices mapping from local to global matrix indices. ===

    rowFeSpace->getBasisFcts()->getLocalIndices(rowElInfo->getElement(),
						rowFeSpace->getAdmin(),
						rowIndices);
    if (rowFeSpace == colFeSpace) {
      colIndices = rowIndices;
    } else {
      if (colElInfo) {
	colFeSpace->getBasisFcts()->getLocalIndices(colElInfo->getElement(),
						    colFeSpace->getAdmin(),
						    colIndices);
      } else {
	// If there is no colElInfo pointer, use rowElInfo the get the indices.
	colFeSpace->getBasisFcts()->getLocalIndices(rowElInfo->getElement(),
						    colFeSpace->getAdmin(),
						    colIndices);
      }
    }

    for (int i = 0; i < nRow; i++)  {
      DegreeOfFreedom row = rowIndices[i];

      BoundaryCondition *condition = 
	bound ? boundaryManager->getBoundaryCondition(bound[i]) : NULL;

      if (condition && condition->isDirichlet()) {	
	if (condition->applyBoundaryCondition())
	  dirichletDofs.insert(row);
      } else {
	for (int j = 0; j < nCol; j++) {
	  DegreeOfFreedom col = colIndices[j];
	  ins[row][col] += elMat[i][j];
	}
      }
    }
  }


  void DOFMatrix::assembleOperator(Operator &op)
  {
    FUNCNAME("DOFMatrix::assembleOperator()");

    TEST_EXIT(rowFeSpace->getMesh() == colFeSpace->getMesh())
      ("This function does not support for multi mesh procedure!\n");
    TEST_EXIT(op.getRowFeSpace() == rowFeSpace)
      ("Row FE spaces do not fit together!\n");
    TEST_EXIT(op.getColFeSpace() == colFeSpace)
      ("Column FE spaces do not fit together!\n");

    clearOperators();
    addOperator(&op);

    matrix.change_dim(rowFeSpace->getAdmin()->getUsedSize(),
		      colFeSpace->getAdmin()->getUsedSize());

    Mesh *mesh = rowFeSpace->getMesh();
    mesh->dofCompress();
    const BasisFunction *basisFcts = rowFeSpace->getBasisFcts();

    Flag assembleFlag = getAssembleFlag() | 
      Mesh::CALL_LEAF_EL                        | 
      Mesh::FILL_COORDS                         |
      Mesh::FILL_DET                            |
      Mesh::FILL_GRD_LAMBDA |
      Mesh::FILL_NEIGH |
      Mesh::FILL_BOUND;

    BoundaryType *bound = new BoundaryType[basisFcts->getNumber()];

    calculateNnz();
    if (getBoundaryManager())
      getBoundaryManager()->initMatrix(this);

    startInsertion(getNnz());

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, assembleFlag);
    while (elInfo) {
      basisFcts->getBound(elInfo, bound);

      assemble(1.0, elInfo, bound);

      if (getBoundaryManager())
	getBoundaryManager()->fillBoundaryConditions(elInfo, this);

      elInfo = stack.traverseNext(elInfo);
    }

    clearDirichletRows();
    finishAssembling();
    finishInsertion();
    getBoundaryManager()->exitMatrix(this);

    delete [] bound;
  }
  

  void DOFMatrix::assemble(double factor, 
			   ElInfo *elInfo, 
			   const BoundaryType *bound)
  {
    set_to_zero(elementMatrix);

    std::vector<Operator*>::iterator it = operators.begin();
    std::vector<double*>::iterator factorIt = operatorFactor.begin();
    for (; it != operators.end(); ++it, ++factorIt)
      if ((*it)->getNeedDualTraverse() == false && 
	  (*factorIt == NULL || **factorIt != 0.0))
	(*it)->getElementMatrix(elInfo,	elementMatrix, *factorIt ? **factorIt : 1.0);

    if (factor != 1.0)
      elementMatrix *= factor;

    if (operators.size())
      addElementMatrix(elementMatrix, bound, elInfo, NULL); 
  }


  void DOFMatrix::assemble(double factor, 
			   ElInfo *elInfo, 
			   const BoundaryType *bound,
			   Operator *op)
  {
    FUNCNAME_DBG("DOFMatrix::assemble()");
    
    TEST_EXIT_DBG(op)("No operator!\n");
    
    set_to_zero(elementMatrix);
    op->getElementMatrix(elInfo, elementMatrix, factor);
    
    if (factor != 1.0)
	elementMatrix *= factor;
    
    addElementMatrix(elementMatrix, bound, elInfo, NULL);
  }

  
  void DOFMatrix::finishAssembling()
  {
    // call the operators cleanup procedures
    for (std::vector<Operator*>::iterator it = operators.begin();
	 it != operators.end(); ++it)
      (*it)->finishAssembling();
  }


  // Should work as before
  Flag DOFMatrix::getAssembleFlag()
  {
    Flag fillFlag(0);
    for (std::vector<Operator*>::iterator op = operators.begin(); 
	 op != operators.end(); ++op)
      fillFlag |= (*op)->getFillFlag();

    return fillFlag;
  }

  
  void DOFMatrix::addOperator(Operator *op, double* factor, double* estFactor) 
  { 
    operators.push_back(op);
    operatorFactor.push_back(factor);
    operatorEstFactor.push_back(estFactor);
  }


  void DOFMatrix::clearOperators()
  {
    operators.clear();
    operatorFactor.clear();
    operatorEstFactor.clear();
  }


  void DOFMatrix::copy(const DOFMatrix& rhs) 
  {
    matrix = rhs.matrix;
  }


  void DOFMatrix::clearDirichletRows()
  {      
    // Do the following only in sequential code. In parallel mode, the specific
    // solver method must care about dirichlet boundary conditions.
    inserter_type &ins = *inserter;  
    for (std::set<DegreeOfFreedom>::iterator it = dirichletDofs.begin(); 
	 it != dirichletDofs.end(); ++it)
      ins[(*it)][(*it)] = 1.0;
  }


  int DOFMatrix::memsize() const
  {   
    return (num_rows(matrix) + matrix.nnz()) * sizeof(base_matrix_type::size_type)
      + matrix.nnz() * sizeof(base_matrix_type::value_type);
  }


  void DOFMatrix::startInsertion(int nnz_per_row)
  {
    if (inserter) {
      delete inserter;
      inserter = NULL; 
    }

    inserter = new inserter_type(matrix, nnz_per_row);

    dirichletDofs.clear();
  }

} // end namespace AMDiS
