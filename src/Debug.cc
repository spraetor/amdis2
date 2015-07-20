#include <iostream>
#include <fstream>
#include <string>

#include "Debug.h"
#include "DOFVector.h"
#include "MacroElement.h"
#include "ElementDofIterator.h"
#include "io/VtkWriter.h"
#include "io/ElementFileWriter.h"

namespace AMDiS {

  namespace debug {

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    void writeLocalElementDofs(int rank, int elIdx, 
			       const FiniteElemSpace *feSpace)
    {      
      if (MPI::COMM_WORLD.Get_rank() == rank) {
	DOFVector<double> tmp(feSpace, "tmp");
	colorDofVectorByLocalElementDofs(tmp, feSpace->getMesh(), elIdx);
	io::VtkWriter::writeFile(tmp, "tmp" + std::to_string(elIdx) + ".vtu");
      }
    }

    
    void writeDofMesh(int rank, DegreeOfFreedom dof, 
		      const FiniteElemSpace *feSpace)
    {
      if (MPI::COMM_WORLD.Get_rank() == rank) {
	DOFVector<double> tmp(feSpace, "tmp");
	tmp.set(0.0);
	tmp[dof] = 1.0;
	io::VtkWriter::writeFile(tmp, "dofmesh" + std::to_string(rank) + ".vtu");
      }    
    }

    
    void writeMesh(const FiniteElemSpace *feSpace, int rank, std::string filename)
    {
      int myRank = MPI::COMM_WORLD.Get_rank();
      if (rank == -1 || myRank == rank) {
	DOFVector<double> tmp(feSpace, "tmp");
	io::VtkWriter::writeFile(tmp, filename + std::to_string(myRank) + ".vtu", 
				 io::VtkWriter::ASCII,
				 false,
				 false);
      }
    }
#endif

    
    void writeDofIndexMesh(const FiniteElemSpace *feSpace,
			   std::string filename)
    {
      DOFVector<double> tmp(feSpace, "tmp");
      DOFIterator<double> it(&tmp, USED_DOFS);
      for (it.reset(); !it.end(); ++it)
	*it = it.getDOFIndex();
      io::VtkWriter::writeFile(tmp, filename);
    }


    void colorEdgeInMesh(const FiniteElemSpace *feSpace,
			 Element *el, 
			 int localEdgeNo, 
			 std::string filename)
    {
      DOFVector<double> tmp(feSpace, "tmp");
      tmp.set(0.0);
      tmp[el->getEdge(localEdgeNo).first] = 1.0;
      tmp[el->getEdge(localEdgeNo).second] = 1.0;
      io::VtkWriter::writeFile(tmp, filename);
    }


    void writeElementIndexMesh(Mesh *mesh, std::string filename, int level)
    {
      std::map<int, double> vec;    
      TraverseStack stack;
      ElInfo *elInfo = 
	level == -1 ? 
	stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL) :
	stack.traverseFirst(mesh, level, Mesh::CALL_EL_LEVEL);
      
      while (elInfo) {		  
	int index = elInfo->getElement()->getIndex();
	vec[index] = index;
	elInfo = stack.traverseNext(elInfo);
      }

      io::ElementFileWriter::writeFile(vec, mesh, filename, ".vtu", level);
    }


    void writeMacroElementIndexMesh(Mesh *mesh,
				    std::string filename)
    {
      std::map<int, double> vec;    
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);

      while (elInfo) {		  
	int index = elInfo->getElement()->getIndex();
	int macroIndex = elInfo->getMacroElement()->getIndex();
	vec[index] = macroIndex;
	elInfo = stack.traverseNext(elInfo);
      }

      io::ElementFileWriter::writeFile(vec, mesh, filename, ".vtu");
    }


    void highlightElementIndexMesh(Mesh *mesh, int idx, std::string filename)
    {
      std::map<int, double> vec;    
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      bool markChildren = false;
      int markLevel = -1;

      while (elInfo) {		  
	if (markChildren && elInfo->getLevel() <= markLevel)
	  markChildren = false;

	int index = elInfo->getElement()->getIndex();
	if (index == idx) {
	  markChildren = true;
	  markLevel = elInfo->getLevel();
	}
  
	if (elInfo->getElement()->isLeaf())
	  vec[index] = (markChildren ? 1.0 : 0.0);

	elInfo = stack.traverseNext(elInfo);
      }

      io::ElementFileWriter::writeFile(vec, mesh, filename);
    }


    void colorMeshByMacroIndex(Mesh *mesh, std::string filename)
    {
      std::map<int, double> vec;    
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      
      while (elInfo) {		  
	int index = elInfo->getElement()->getIndex();
	vec[index] = elInfo->getMacroElement()->getIndex();
	elInfo = stack.traverseNext(elInfo);
      }

      io::ElementFileWriter::writeFile(vec, mesh, filename);
    }


    void colorDofVectorByLocalElementDofs(DOFVector<double>& vec, Element *el)
    {
      // === Get local indices of the given element. ===
      
      const BasisFunction *basisFcts = vec.getFeSpace()->getBasisFcts();
      int nBasisFcts = basisFcts->getNumber();
      std::vector<DegreeOfFreedom> localDofs(nBasisFcts);
      basisFcts->getLocalIndices(el, vec.getFeSpace()->getAdmin(), localDofs);
      
      // === Set the values of the dof vector. ===
      
      vec.set(0.0);
      for (int i = 0; i < nBasisFcts; i++)
	vec[localDofs[i]] = static_cast<double>(i);
    }

    
    bool colorDofVectorByLocalElementDofs(DOFVector<double>& vec, Mesh *mesh, 
					  int elIndex)
    {      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo) {
	if (elInfo->getElement()->getIndex() == elIndex) {
	  colorDofVectorByLocalElementDofs(vec, elInfo->getElement());
	  return true;
	}
	elInfo = stack.traverseNext(elInfo);
      }
      
      return false;
    }

    
    Element* getDofIndexElement(const FiniteElemSpace *feSpace, 
				DegreeOfFreedom dof)
    {
      const BasisFunction* basFcts = feSpace->getBasisFcts();
      int nBasFcts = basFcts->getNumber();
      std::vector<DegreeOfFreedom> dofVec(nBasFcts);
      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL);
      while (elInfo) {
	basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), dofVec);
	for (int i = 0; i < nBasFcts; i++) 
	  if (dofVec[i] == dof)
	    return elInfo->getElement();	  
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;
    }

    
    Element* getLevel0ParentElement(Mesh *mesh, Element *el)
    {    
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement() == el)
	  return elInfo->getMacroElement()->getElement();      
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;
    }


    Element* getLevel0ParentElement(Mesh *mesh, int elIndex)
    {    
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement()->getIndex() == elIndex)
	  return elInfo->getMacroElement()->getElement();      
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;
    }


    Element* getParentElement(Mesh *mesh, Element *el)
    {
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement()->getChild(0) == el || 
	    elInfo->getElement()->getChild(1) == el)
	  return elInfo->getElement();      
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;
    }


    Element* getElement(Mesh *mesh, int elIndex)
    {
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement()->getIndex() == elIndex)
	  return elInfo->getElement();      
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;      
    }


    Element* getParentElement(Mesh *mesh, int elIndex)
    {
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement()->isLeaf() == false) 
	  if (elInfo->getElement()->getChild(0)->getIndex() == elIndex || 
	      elInfo->getElement()->getChild(1)->getIndex() == elIndex)
	    return elInfo->getElement();
	
	elInfo = stack.traverseNext(elInfo);
      }
      
      return NULL;

    }


    void printElementInfo(Element *el)
    {
      FUNCNAME("printElementInfo()");
      
      if (el->isLeaf()) {
	MSG("el %d is leaf!\n", el->getIndex());
      } else {
	MSG("el child0 %d    child1 %d\n", 
	    el->getChild(0)->getIndex(),
	    el->getChild(1)->getIndex());
	printElementInfo(el->getChild(0));
	printElementInfo(el->getChild(1));
      }
    }

    
    void printElementCoords(const FiniteElemSpace *feSpace, Element *el)
    {
      FUNCNAME("debug:printElementCoords()");
      
      Mesh *mesh = feSpace->getMesh();

      for (int i = 0; i <= mesh->getDim(); i++) {
	DegreeOfFreedom dof = el->getDof(i, 0);
	WorldVector<double> coords;
	mesh->getDofIndexCoords(dof, feSpace, coords);
	MSG("%d-th DOF of element %d: %f %f %f\n", 
	    i, el->getIndex(), coords[0], coords[1], coords[2]);
      }
    }


    void printInfoByDof(const FiniteElemSpace *feSpace, DegreeOfFreedom dof)
    {
      FUNCNAME("debug::printInfoByDof()");

      WorldVector<double> coords;
      feSpace->getMesh()->getDofIndexCoords(dof, feSpace, coords);

      Element *el = getDofIndexElement(feSpace, dof);
      Element *parEl = getLevel0ParentElement(feSpace->getMesh(), el);
      
      if (coords.getSize() == 2)
	MSG("[DBG] DOF-INFO:  dof = %d  coords: %e %e\n", dof, coords[0], coords[1]);
      else 
	MSG("[DBG] DOF-INFO:  dof = %d  coords: %e %e %e\n", dof, coords[0], coords[1], coords[2]);

      MSG("[DBG] DOF-INFO:  dof = %d  elidx = %d  pelidx = %d\n", 
	  dof, el->getIndex(), parEl->getIndex());
      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(), -1, 
					   Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (elInfo->getElement()->getIndex() == parEl->getIndex()) {
	  MSG("[DBG] EL INFO TO %d: type = %d\n", parEl->getIndex(), elInfo->getType());
	  break;
	}

	elInfo = stack.traverseNext(elInfo);
      }	  
    }


    void printMatValuesStatistics(Matrix<DOFMatrix*> *mat)
    {
      std::map<int, int> counter;
      int counter0 = 0;

      for (int i = 0; i < mat->getNumRows(); i++) {
	for (int j = 0; j < mat->getNumCols(); j++) {
	  DOFMatrix *dofMat = (*mat)[i][j];
	  if (dofMat) {

	    using mtl::tag::major; using mtl::tag::nz; using mtl::begin; using mtl::end;
	    namespace traits = mtl::traits;
	    typedef DOFMatrix::base_matrix_type base_matrix_type;
	    
	    traits::const_value<base_matrix_type>::type value(dofMat->getBaseMatrix());
	    typedef traits::range_generator<major, base_matrix_type>::type cursor_type;
	    typedef traits::range_generator<nz, cursor_type>::type icursor_type;
    
	    for (cursor_type cursor = begin<major>(dofMat->getBaseMatrix()), 
		   cend = end<major>(dofMat->getBaseMatrix()); cursor != cend; ++cursor)
	      for (icursor_type icursor = begin<nz>(cursor), 
		     icend = end<nz>(cursor); icursor != icend; ++icursor) {
		if (value(*icursor) == 0.0) {
		  counter0++;
		} else {
		  int a = static_cast<int>(std::floor(std::log10(std::abs(value(*icursor)))));
		  counter[a]++;
		}
	      }
		
	  }
	}
      }

      std::cout << "value = 0: " << counter0 << std::endl;
      for (std::map<int, int>::iterator it = counter.begin(); it != counter.end(); ++it)
	std::cout << std::pow(10.0, it->first) << " <= values <= " 
		  << std::pow(10.0, it->first + 1.0) << ": " 
		  << it->second << std::endl;
    }


    void printAllDofCoords(const FiniteElemSpace *feSpace)
    {
      FUNCNAME("printAllDofCoords()");

      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(), -1, 
					   Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
      while (elInfo) {
	Element *el = elInfo->getElement();
	for (int i = 0; i <= feSpace->getMesh()->getDim(); i++) {
	  MSG("Coords for DOF %d = %f %f\n", 
	      *(el->getDof(i)), (elInfo->getCoord(i))[0], (elInfo->getCoord(i))[1]);
	}
	elInfo = stack.traverseNext(elInfo);
      }      
    }


    void getAllDofs(const FiniteElemSpace *feSpace, 
		                std::set<const DegreeOfFreedom*>& dofs)
    {
      ElementDofIterator elDofIter(feSpace);
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(feSpace->getMesh(), -1, Mesh::CALL_LEAF_EL);
      while (elInfo) {	
      	elDofIter.reset(elInfo->getElement());
      	do {
      	  dofs.insert(elDofIter.getDofPtr());
      	} while (elDofIter.next());
      	elInfo = stack.traverseNext(elInfo);
      }      
    }


    void writeMatlabMatrix(DOFMatrix &mat, std::string filename)
    {
      using mtl::tag::major; using mtl::tag::nz; using mtl::begin; using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      traits::row<Matrix>::type row(mat.getBaseMatrix());
      traits::col<Matrix>::type col(mat.getBaseMatrix());
      traits::const_value<Matrix>::type value(mat.getBaseMatrix());

      typedef traits::range_generator<major, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      std::ofstream out;
      out.open(filename.c_str());
      out.precision(20);
    
      for (cursor_type cursor = begin<major>(mat.getBaseMatrix()), 
	     cend = end<major>(mat.getBaseMatrix()); cursor != cend; ++cursor)
	for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); icursor != icend; ++icursor)
	  out << row(*icursor) + 1 << " " 
	      << col(*icursor) + 1 << " " 
	      << value(*icursor) << "\n";

      out.close();
    }


    void writeMatlabMatrix(Matrix<DOFMatrix*> &mat, std::string filename)
    {
      using mtl::tag::major; using mtl::tag::nz; using mtl::begin; using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      std::ofstream out;
      out.open(filename.c_str());

      for (int i = 0; i < mat.getNumRows(); i++) {
	for (int j = 0; j < mat.getNumCols(); j++) {
	  if (!mat[i][j])
	    continue;

	  traits::row<Matrix>::type row(mat[i][j]->getBaseMatrix());
	  traits::col<Matrix>::type col(mat[i][j]->getBaseMatrix());
	  traits::const_value<Matrix>::type value(mat[i][j]->getBaseMatrix());
	  
	  typedef traits::range_generator<major, Matrix>::type cursor_type;
	  typedef traits::range_generator<nz, cursor_type>::type icursor_type;	  

	  for (cursor_type cursor = begin<major>(mat[i][j]->getBaseMatrix()), 
		 cend = end<major>(mat[i][j]->getBaseMatrix()); cursor != cend; ++cursor)
	    for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor); icursor != icend; ++icursor)
	      out << i * num_rows(mat[i][j]->getBaseMatrix()) + row(*icursor) + 1 << " " 
		  << j * num_cols(mat[i][j]->getBaseMatrix()) + col(*icursor) + 1 << " " 
		  << value(*icursor) << "\n";
	}
      }
    
      out.close();
    }


    void writeMatlabVector(DOFVector<double> &vec, std::string filename)
    {
      std::ofstream out;    
      out.open(filename.c_str());

      DOFIterator<double> it(&vec, USED_DOFS);
      for (it.reset(); !it.end(); ++it)
	out << *it << "\n";

      out.close();      
    }

    
    void writeMatlabVector(SystemVector &vec, std::string filename)
    {
      std::ofstream out;    
      out.open(filename.c_str());

      for (int i = 0; i < vec.getSize(); i++) {
	DOFIterator<double> it(vec.getDOFVector(i), USED_DOFS);
	for (it.reset(); !it.end(); ++it)
	  out << *it << "\n";
      }

      out.close();            
    }


    void writeCoordsFile(const FiniteElemSpace* feSpace, std::string filename)
    {
      DOFVector<WorldVector<double> > coords(feSpace, "tmp");
      feSpace->getMesh()->getDofIndexCoords(coords);

      std::ofstream file;
      file.open(filename.c_str());
      file << coords.getUsedSize() << "\n";
      DOFIterator<WorldVector<double> > it(&coords, USED_DOFS);
      for (it.reset(); !it.end(); ++it) {
	file << it.getDOFIndex();
	for (int i = 0; i < feSpace->getMesh()->getDim(); i++)
	  file << " " << (*it)[i];
	file << "\n";
      }
      file.close();
    }


    void printElementHierarchie(Mesh *mesh, int elIndex)
    {
      FUNCNAME("debug::printElementHierarchie()");

      bool printInfo = false;
      int elLevel = -1;

      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo) {
	if (printInfo && elInfo->getLevel() <= elLevel)
	  return;

	if (elInfo->getElement()->getIndex() == elIndex) {
	  printInfo = true;
	  elLevel = elInfo->getLevel();
	  MSG(" Hierarchie for elIdx = %d\n", elIndex);
	} else {
	  if (printInfo) {
	    std::stringstream oss;
	    for (int i = 0; i < (elInfo->getLevel() - elLevel); i++)
	      oss << "   ";
	    oss << "|--" << elInfo->getElement()->getIndex();
	    if (oss.str().length() >= 255) {
	      WARNING("String is to long! Element index is %d.\n", elInfo->getElement()->getIndex());
	    } else {
	      MSG("%s\n", oss.str().c_str());
	    }
	  }
	}

	elInfo = stack.traverseNext(elInfo);
      }

      if (!printInfo) {
	MSG("Could not find element with index %d\n", elIndex);
      }
    }


    void printElementRefinementSequence(Mesh *mesh, Element *el)
    {
      FUNCNAME("printElementRefinementSequence()");

      int elIndex = el->getIndex();
      std::vector<int> refSeq;
      Element *parent = getParentElement(mesh, el);

      while (parent) {
	if (parent->getChild(0) == el)
	  refSeq.push_back(0);
	else
	  refSeq.push_back(1);

	el = parent;
	parent = getParentElement(mesh, el);
      }

      std::stringstream oss;
      for (int i = refSeq.size() - 1; i >= 0; i--)
	oss << refSeq[i];
      
      MSG("Ref-Seq for element %d: %s\n", elIndex, oss.str().c_str());
    }


    int getLocalNeighbourIndex(Mesh *mesh, int elIndex, int neighIndex)
    {
      FUNCNAME("debug::getLocalNeighbourIndex()");

      TraverseStack stack;
      ElInfo *elInfo = 
	stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER | Mesh::FILL_NEIGH);
      while (elInfo) {
	if (elInfo->getElement()->getIndex() == elIndex) {
	  int returnValue = -1;
	  for (int i = 0; i <= mesh->getDim(); i++) {
	    if (elInfo->getNeighbour(i)) {
	      MSG("NEIGH %d of El-Idx %d: %d\n", i, elIndex, elInfo->getNeighbour(i)->getIndex());
	    } else {
	      MSG("NEIGH %d of El-Idx %d: %d\n", i, elIndex, -1);
	    }
	    
	    if (elInfo->getNeighbour(i) && 
		elInfo->getNeighbour(i)->getIndex() == neighIndex)
	      returnValue = i;
	  }

	  return returnValue;	  
	}
	elInfo = stack.traverseNext(elInfo);
      }      

      return -2;
    }


    void importDofVectorByCoords(DOFVector<double>* vec, std::string filename)
    { 
      DOFVector<WorldVector<double> > coords(vec->getFeSpace(), "dofCoords");
      vec->getFeSpace()->getMesh()->getDofIndexCoords(coords);
      int dim = vec->getFeSpace()->getMesh()->getDim();

      std::ifstream file;
      file.open(filename.c_str());
      
      int nEntries;
      file >> nEntries;

      for (int i = 0; i < nEntries; i++) {
	WorldVector<double> dofCoords;
	for (int j = 0; j < dim; j++)
	  file >> dofCoords[j];

	double value;
	file >> value;

      
	DOFIterator<WorldVector<double> > it(&coords, USED_DOFS);
	for (it.reset(); !it.end(); ++it) {
	  bool found = true;
	  for (int j = 0; j < dim; j++) {
	    if (fabs((*it)[j] - dofCoords[j]) > 1e-8) {
	      found = false;
	      break;
	    }
	  }

	  if (found) {
	    (*vec)[it.getDOFIndex()] = value;
	    break;
	  }
	}	
      }

      file.close();
    }


    void exportDofVectorByCoords(const DOFVector<double>* vec, 
				 std::string filename)
    {
      DOFVector<WorldVector<double> > coords(vec->getFeSpace(), "dofCoords");
      vec->getFeSpace()->getMesh()->getDofIndexCoords(coords);
      int dim = vec->getFeSpace()->getMesh()->getDim();

      std::ofstream file;
      file.open(filename.c_str());
      file << vec->getUsedSize() << "\n";

      DOFIterator<WorldVector<double> > it(&coords, USED_DOFS);
      for (it.reset(); !it.end(); ++it) {
        for (int i = 0; i < dim; i++)
	  file << (*it)[i] << " ";
        file << (*vec)[it.getDOFIndex()] << "\n";
      }
    
      file.close();
    }


    void createSortedDofs(Mesh *mesh, ElementIdxToDofs &elMap)
    {
      FUNCNAME("debug::dbgCreateElementMap()");
      
      int dim = mesh->getDim();
      elMap.clear();
      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo) {
	Element *el = elInfo->getElement();
	switch (dim) {
	case 2:
	  sortDofs(el->getDof(0), el->getDof(1), el->getDof(2), elMap[el->getIndex()]);
	  break;
	case 3:
	  sortDofs(el->getDof(0), el->getDof(1), el->getDof(2), el->getDof(3), elMap[el->getIndex()]);
	  break;
	default:
	  ERROR_EXIT("What is this?\n");
	}
	elInfo = stack.traverseNext(elInfo);
      }
    }
    
    void testSortedDofs(Mesh *mesh, ElementIdxToDofs &elMap)
    {
      FUNCNAME("debug::dbgTestElementMap()");
      
      int dim = mesh->getDim();
      int nVertex = Global::getGeo(VERTEX, dim);
      DofContainer vec(nVertex);
      
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo) {
	Element *el = elInfo->getElement();
	
	switch (dim) {
	case 2:
	  sortDofs(el->getDof(0), el->getDof(1), el->getDof(2), vec);
	  break;
	case 3:
	  sortDofs(el->getDof(0), el->getDof(1), el->getDof(2), el->getDof(3), vec);
	  break;
	default:
	  ERROR_EXIT("What is this?\n");
	}
	
	for (int i = 0; i < nVertex; i++) {
	  if (elMap[el->getIndex()][i] != vec[i]) {
	    MSG("[DBG] Wrong new dof numeration in element = %d\n!", el->getIndex());
	    std::cout << "[DBG]: Old numeration was: ";
	    for (int j = 0; j < nVertex; j++)
	      std::cout << elMap[el->getIndex()][j] << " = " 
			<< *(elMap[el->getIndex()][j]) << "  ";
	    std::cout << std::endl;
	    std::cout << "[DBG]: New numeration is:  ";
	    for (int j = 0; j < nVertex; j++)
	      std::cout << vec[j] << " = "  << *(vec[j]) << "  ";
	    std::cout << std::endl;
	    ERROR_EXIT("WRONG NEW DOF NUMERATION!\n");
	  }
	}
	elInfo = stack.traverseNext(elInfo);
      }
    }


    void createNodeCoords(Mesh *mesh, ElementIdxToCoords& coords)
    {
      TraverseStack stack;
      ElInfo *elInfo = 
        stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);

      while (elInfo) {
	coords.insert(std::make_pair(elInfo->getElement()->getIndex(), elInfo->getCoords()));
	  
	elInfo = stack.traverseNext(elInfo);
      } 
    }
    
    void testNodeCoords(Mesh* mesh, ElementIdxToCoords& coords)
    {
      FUNCNAME("debug::testNodeCoords()");
      
      int dim = mesh->getDim();
      
      TraverseStack stack;
      ElInfo *elInfo = 
        stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);

      while (elInfo) {
	Element *el = elInfo->getElement();
	
	FixVec<WorldVector<double>, VERTEX> elCoords = elInfo->getCoords();
	  
	for (size_t i = 0; i < (size_t)(size(elCoords)); i++) 
	  for (size_t j = 0; j < (size_t)(mesh->getDim()); j++) 
	    if(elCoords[i][j] != coords[el->getIndex()][i][j]) {
	      
	      MSG("[DBG] Wrong coordnate on element = %d, vertex index = %d\n!", el->getIndex(), i);
	      
	      switch (dim) {
	      case 2:
		std::cout << "[DBG]: one coord is: ";
		std::cout << "(" << elCoords[i][0] << ", " << elCoords[i][1]<< ") ";
		std::cout << "another is: ";
		std::cout << "(" << coords[el->getIndex()][i][0] << ", " << coords[el->getIndex()][i][1] << ")\n";
		break;
	      case 3:
		std::cout << "[DBG]: one coord is: ";
		std::cout << "(" << elCoords[i][0] << ", " << elCoords[i][1]<< ", " << elCoords[i][2] << ") ";
		std::cout << "another is: ";
		std::cout << "(" << coords[el->getIndex()][i][0] << ", " << coords[el->getIndex()][i][1]
		  << ", " << coords[el->getIndex()][i][2] << ")\n";
		break;
	      default:
		ERROR_EXIT("What is this?\n");
	      }
	      ERROR_EXIT("Mesh vertex coords don't match!\n");
	    }
	  
	elInfo = stack.traverseNext(elInfo);
      } 
    }
    
    void sortDofs(const DegreeOfFreedom* dof0,
		  const DegreeOfFreedom* dof1,
		  const DegreeOfFreedom* dof2,
		  DofContainer &vec)
    {
      DofPtrSortFct dofPtrSort;
      vec.resize(3);
      vec[0] = dof0; 
      vec[1] = dof1; 
      vec[2] = dof2;
      sort(vec.begin(), vec.end(), dofPtrSort);
    }

    void sortDofs(const DegreeOfFreedom* dof0,
		  const DegreeOfFreedom* dof1,
		  const DegreeOfFreedom* dof2,
		  const DegreeOfFreedom* dof3,
		  DofContainer &vec)
    {
      DofPtrSortFct dofPtrSort;
      vec.resize(4);
      vec[0] = dof0; 
      vec[1] = dof1; 
      vec[2] = dof2;
      vec[3] = dof3;
      sort(vec.begin(), vec.end(), dofPtrSort);
    }


    void testDofsByCoords(const FiniteElemSpace *feSpace,
			  DofContainer &dofs0, DofContainer &dofs1)
    {
      FUNCNAME("debug::testDofsByCoords()");

      TEST_EXIT(dofs0.size() == dofs1.size())
	("The dof containers have different sizes %d %d!\n",
	 dofs0.size(), dofs1.size());

      DOFVector<WorldVector<double> > coords(feSpace, "dofCorrds");
      feSpace->getMesh()->getDofIndexCoords(coords);

      for (unsigned int i = 0; i < dofs0.size(); i++) {
	WorldVector<double> tmp = coords[*(dofs0[i])];
	tmp -= coords[*(dofs1[i])];

	TEST_EXIT(norm(tmp) < 1e-13)
	  ("DOFs %d and %d (i = %d)  have different coords!\n", 
 	   *(dofs0[i]), *(dofs1[i]), i);
      }
      
    }


    void testDofsByCoords(DOFVector<WorldVector<double> > &coords,
			  DofContainer &dofs0, 
			  DofContainer &dofs1)
    {
      FUNCNAME("debug::testDofsByCoords()");

      TEST_EXIT(dofs0.size() == dofs1.size())
	("The dof containers have different sizes %d %d!\n",
	 dofs0.size(), dofs1.size());

      for (unsigned int i = 0; i < dofs0.size(); i++) {
	WorldVector<double> tmp = coords[*(dofs0[i])];
	tmp -= coords[*(dofs1[i])];
	

	TEST(norm(tmp) < 1e-13)("DOFs 0: %e : %e : %e \n DOFs 1: %e : %e : %e \n norm: %e \n", 
	   coords[*(dofs0[i])][0], coords[*(dofs0[i])][1], 
	   coords[*(dofs0[i])][2],
	   coords[*(dofs1[i])][0], coords[*(dofs1[i])][1],
				coords[*(dofs1[i])][2], norm(tmp));
	
	TEST_EXIT(norm(tmp) < 1e-13)
	  ("DOFs %d and %d (i = %d)  have different coords!\n", 
	   *(dofs0[i]), *(dofs1[i]), i);
      }
      
    }

  } // namespace debug
  
} // namespace AMDiS
