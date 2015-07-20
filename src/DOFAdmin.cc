#include <algorithm>

#include "QPsiPhi.h"
#include "BasisFunction.h"
#include "Boundary.h"
#include "DOFAdmin.h"
#include "ElInfo.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "DOFVector.h"
#include "DOFIterator.h"

namespace AMDiS 
{
  const int DOFAdmin::sizeIncrement = 10;

  DOFAdmin::DOFAdmin(Mesh* m) 
    : mesh(m), 
      nDof(mesh->getDim(), NO_INIT),
      nPreDof(mesh->getDim(), NO_INIT)
  { 
    init(); 
  }


  DOFAdmin::DOFAdmin(Mesh* m, std::string aName) 
    : name(aName), 
      mesh(m), 
      nDof(mesh->getDim(), NO_INIT),
      nPreDof(mesh->getDim(), NO_INIT)
  { 
    init(); 
  }


  DOFAdmin::~DOFAdmin() 
  {}


  void DOFAdmin::init()
  {
    firstHole = 0;
    size = 0;
    usedCount = 0;
    holeCount = 0;
    sizeUsed = 0;
    dofFree.clear();
  }


  DOFAdmin& DOFAdmin::operator=(const DOFAdmin& src) 
  {
    if (this != &src) { 
      mesh = src.mesh;
      name = src.name;
      dofFree = src.dofFree;
      firstHole = src.firstHole;
      size = src.size;
      usedCount = src.usedCount;
      holeCount = src.holeCount;
      sizeUsed = src.sizeUsed;
      for (int i = 0; i <= mesh->getDim(); i++) {
	nDof[i] = src.nDof[i];
	nPreDof[i] = src.nPreDof[i];
      }
      dofIndexedList = src.dofIndexedList;
      dofContainerList = src.dofContainerList;
    }

    return *this;
  }


  bool DOFAdmin::operator==(const DOFAdmin& ad) const
  {
    if (name != ad.name) 
      return false;
    if (mesh != ad.mesh) 
      return false;

    return true;
  }


  DOFAdmin::DOFAdmin(const DOFAdmin&)
  {
    FUNCNAME("DOFAdmin::DOFAdmin()");

    ERROR_EXIT("TODO\n");
  }


  void DOFAdmin::freeDofIndex(DegreeOfFreedom dof) 
  {    
    FUNCNAME_DBG("DOFAdmin::freeDofIndex()");

    TEST_EXIT_DBG(usedCount > 0)("No DOFs in use!\n");
    TEST_EXIT_DBG(dof >= 0 && dof < size)("Invalid DOF index %d!\n", dof);

    std::list<DOFIndexedBase*>::iterator di;
    std::list<DOFIndexedBase*>::iterator end = dofIndexedList.end();

    for (di = dofIndexedList.begin(); di != end; ++di)
      (*di)->freeDOFContent(dof);

    std::list<DOFContainer*>::iterator dc;
    std::list<DOFContainer*>::iterator dcend = dofContainerList.end();

    for (dc = dofContainerList.begin(); dc != dcend; ++dc)
      (*dc)->freeDofIndex(dof);

    dofFree[dof] = true;
    
    if (firstHole > dof) 
      firstHole = dof;

    usedCount--;
    holeCount++;
  }


  DegreeOfFreedom DOFAdmin::getDOFIndex()
  {
    FUNCNAME_DBG("DOFAdmin::getDOFIndex()");
    DegreeOfFreedom dof = 0;

    // if there is a hole
    if (firstHole < static_cast<DofIndex::size_type>(dofFree.size())) {      
      TEST_EXIT_DBG(dofFree[firstHole])("no hole at firstHole!\n");
      // its no longer a hole
      dofFree[firstHole] = false;
      dof = firstHole;
      // search new hole
      DofIndex::size_type dfsize = static_cast<DofIndex::size_type>(dofFree.size());
      DofIndex::size_type i = firstHole + 1;
      for (; i < dfsize; i++)
	if (dofFree[i])
	  break;

      firstHole = i;
    } else {                    // if there is no hole
      // enlarge dof-list
      enlargeDofLists();

      TEST_EXIT_DBG(firstHole < static_cast<DofIndex::size_type>(dofFree.size()))
	("no free entry after enlargeDofLists\n");
      TEST_EXIT_DBG(dofFree[firstHole])("no free bit at firstHole\n");
      dofFree[firstHole] = false;
      dof = firstHole;
      firstHole++;
    }

    usedCount++;
    if (holeCount > 0) 
      holeCount--;
    sizeUsed = std::max(sizeUsed, dof + 1);

    return dof;
  }


  void DOFAdmin::enlargeDofLists(int minsize)
  {  
    DofIndex::size_type old = size;
    if (minsize > 0)
      if (old > minsize) 
	return;
  
    DofIndex::size_type newval = std::max(static_cast<DofIndex::size_type>(minsize), 
					  static_cast<DofIndex::size_type>((dofFree.size() + sizeIncrement)));

    size = newval;
  
    // stl resizes dofFree to at least newval and sets all new values true
    dofFree.resize(newval, true);

    firstHole = old;
  
    // enlarge all vectors and matrices
    // but DOFVectors<int> don't have to be changed
  
    std::list<DOFIndexedBase*>::iterator di;
    for (di = dofIndexedList.begin(); di != dofIndexedList.end(); ++di)
      if ((*di)->getSize() < newval)
 	(*di)->resize(newval);
  }


  void DOFAdmin::addDOFIndexed(DOFIndexedBase* dofIndexed) 
  {
    FUNCNAME("DOFAdmin::addDOFIndexed()");

    TEST_EXIT(dofIndexed)("no dofIndexed\n");

    if (dofIndexed->getSize() < size)
      dofIndexed->resize(size);
    
    dofIndexedList.push_back(dofIndexed);    
  }


  void DOFAdmin::removeDOFIndexed(DOFIndexedBase* dofIndexed)
  {
    FUNCNAME("DOFAdmin::removeDOFIndexed()");

    bool removed = false;
    std::list<DOFIndexedBase*>::iterator it;
    std::list<DOFIndexedBase*>::iterator end = dofIndexedList.end();
    for (it = dofIndexedList.begin(); it != end; ++it) {
      if (*it == dofIndexed) {
	dofIndexedList.erase(it);
	removed = true;
	break;	
      }
    }

    TEST_EXIT(removed)("DOFIndexed not in list!\n");
  }


  void DOFAdmin::addDOFContainer(DOFContainer* cont)
  {
    FUNCNAME_DBG("DOFAdmin::addDOFContainer()");

    TEST_EXIT_DBG(cont)("no container\n");

    dofContainerList.push_back(cont);  
  }


  void DOFAdmin::removeDOFContainer(DOFContainer* cont)
  {
    FUNCNAME("DOFAdmin::removeDOFContainer()");

    std::list<DOFContainer*>::iterator it;
    std::list<DOFContainer*>::iterator end = dofContainerList.end();
    for (it = dofContainerList.begin(); it != end; ++it) {
      if (*it == cont) {
	dofContainerList.erase(it);
	return;
      }
    }

    ERROR("Container not in list!\n");
  }


  void DOFAdmin::compress(std::vector<DegreeOfFreedom> &newDofIndex)
  {
    FUNCNAME_DBG("DOFAdmin::compress()");

    // nothing to do ?
    if (size < 1 || usedCount < 1 || holeCount < 1) 
      return;

    // vector to mark used dofs
    for (int i = 0; i < size; i++)
      newDofIndex[i] = -1;

    // mark used dofs
    DOFIteratorBase it(this, USED_DOFS);
    for (it.reset(); !it.end(); ++it)
      newDofIndex[it.getDOFIndex()] = 1;   
    
    // create a MONOTONE compress
    int n = 0, last = 0;
    for (int i = 0; i < size; i++) {
      if (newDofIndex[i] == 1) {
	newDofIndex[i] = n++;
	last = i;
      }
    }
  
    TEST_EXIT_DBG(n == usedCount)("count %d != usedCount %d\n", n, usedCount);
  
    // mark used dofs in compressed dofFree
    for (int i = 0; i < n; i++)
      dofFree[i] = false;

    // mark unused dofs in compressed dofFree
    for (int i = n; i < size; i++)
      dofFree[i] = true;

    firstHole = n;  
    holeCount = 0;
    sizeUsed  = n;
  
    // get index of first changed dof
    int first = last;
    for (int i = 0; i < size; i++) {
      if (newDofIndex[i] < i && newDofIndex[i] >= 0) {
	first = i;
	break;
      }
    }
    
    for (std::list<DOFIndexedBase*>::iterator di = dofIndexedList.begin(); 
	 di != dofIndexedList.end(); ++di)
      (*di)->compressDOFIndexed(first, last, newDofIndex);

    for (std::list<DOFContainer*>::iterator dc = dofContainerList.begin(); 
	 dc != dofContainerList.end(); ++dc)
      (*dc)->compressDofContainer(n, newDofIndex);
  }


  void DOFAdmin::setNumberOfDofs(int i, int v) 
  { 
    FUNCNAME_DBG("DOFAdmin::setNumberOfDOFs()");

    TEST_EXIT_DBG(0 <= i && 4 > i)("Should not happen!\n");

    nDof[i] = v; 
  }


  void DOFAdmin::setNumberOfPreDofs(int i, int v) 
  { 
    FUNCNAME_DBG("DOFAdmin::setNumberOfPreDOFs()");

    TEST_EXIT_DBG(0 <= i && 4 > i)("Should not happen!\n"); 

    nPreDof[i] = v; 
  }


  int DOFAdmin::calcMemoryUsage() const
  {
    return sizeof(DOFAdmin);
  }

} // end namespace AMDiS
