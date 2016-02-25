#include "BoundaryManager.hpp"

// AMDiS includes
#include "BasisFunction.hpp"
#include "DOFIndexed.hpp"
#include "DOFVector.hpp"
#include "ElInfo.hpp"
#include "FiniteElemSpace.hpp"
#include "Traverse.hpp"

namespace AMDiS
{

  std::map<BoundaryType, std::vector<BoundaryCondition*>>
      BoundaryManager::globalBoundaryMap;


  BoundaryManager::BoundaryManager(not_null<const FiniteElemSpace*> feSpace)
  {
    allocatedMemoryLocalBounds = feSpace->getBasisFcts()->getNumber();
    localBound = new BoundaryType[allocatedMemoryLocalBounds];
  }


  BoundaryManager::BoundaryManager(BoundaryManager const& bm)
  {
    localBCs = bm.localBCs;
    allocatedMemoryLocalBounds = bm.allocatedMemoryLocalBounds;
    localBound = new BoundaryType[allocatedMemoryLocalBounds];
  }


  BoundaryManager::~BoundaryManager()
  {
    delete [] localBound;
  }


  void BoundaryManager::addBoundaryCondition(BoundaryCondition* localBC)
  {
    FUNCNAME("BoundaryManager::addBoundaryCondition()");

    BoundaryType type = localBC->getBoundaryType();
    TEST_EXIT(localBCs[type] == NULL)
    ("There is already a condition for this type %d.\n",type);
    localBCs[type] = localBC;

    std::vector<BoundaryCondition*>& boundMap = globalBoundaryMap[type];
    boundMap.push_back(localBC);
  }


  double BoundaryManager::boundResidual(ElInfo* elInfo,
                                        DOFMatrix* matrix,
                                        const DOFVectorBase<double>* dv)
  {
    double result = 0.0;
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second)
        result += (*it).second->boundResidual(elInfo, matrix, dv);

    return result;
  }


  void BoundaryManager::fillBoundaryConditions(ElInfo* elInfo,
      DOFVectorBase<double>* vec)
  {
    if (localBCs.size() > 0)
    {
      const FiniteElemSpace* feSpace = vec->getFeSpace();
      const BasisFunction* basisFcts = feSpace->getBasisFcts();
      int nBasFcts = basisFcts->getNumber();
      dofVec.resize(nBasFcts);

      // get boundaries of all DOFs
      basisFcts->getBound(elInfo, localBound);

      // get dof indices
      basisFcts->getLocalIndices(elInfo->getElement(),
                                 feSpace->getAdmin(),
                                 dofVec);

      // apply non dirichlet boundary conditions
      for (BoundaryIndexMap::iterator it = localBCs.begin();
           it != localBCs.end(); ++it)
        if ((*it).second && !(*it).second->isDirichlet())
          (*it).second->fillBoundaryCondition(vec, elInfo, &dofVec[0],
                                              localBound, nBasFcts);

      // apply dirichlet boundary conditions
      for (BoundaryIndexMap::iterator it = localBCs.begin();
           it != localBCs.end(); ++it)
        if ((*it).second && (*it).second->isDirichlet())
          (*it).second->fillBoundaryCondition(vec, elInfo, &dofVec[0],
                                              localBound, nBasFcts);
    }
  }


  void BoundaryManager::fillBoundaryConditions(ElInfo* elInfo, DOFMatrix* mat)
  {
    if (localBCs.size() <= 0)
      return;

    const FiniteElemSpace* feSpace = mat->getRowFeSpace();
    const BasisFunction* basisFcts = feSpace->getBasisFcts();
    int nBasFcts = basisFcts->getNumber();
    dofVec.resize(nBasFcts);

    // get boundaries of all DOFs
    basisFcts->getBound(elInfo, localBound);

    // get DOF indices
    basisFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), dofVec);

    // apply non dirichlet boundary conditions
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && !(*it).second->isDirichlet())
        (*it).second->fillBoundaryCondition(mat, elInfo, &dofVec[0],
                                            localBound, nBasFcts);

    // apply dirichlet boundary conditions
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && (*it).second->isDirichlet())
        (*it).second->fillBoundaryCondition(mat, elInfo, &dofVec[0],
                                            localBound, nBasFcts);
  }

  bool BoundaryManager::isBoundaryPeriodic(BoundaryType b)
  {
    for (BoundaryCondition* boundary_map : globalBoundaryMap[b])
      if (boundary_map->isPeriodic())
        return true;

    return false;
  }


  void BoundaryManager::initMatrix(DOFMatrix* matrix)
  {
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && !(*it).second->isDirichlet())
        (*it).second->initMatrix(matrix);

    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && (*it).second->isDirichlet())
        (*it).second->initMatrix(matrix);
  }


  void BoundaryManager::exitMatrix(DOFMatrix* matrix)
  {
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && !(*it).second->isDirichlet())
        (*it).second->exitMatrix(matrix);

    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && (*it).second->isDirichlet())
        (*it).second->exitMatrix(matrix);
  }


  void BoundaryManager::initVector(DOFVectorBase<double>* vector)
  {
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && !(*it).second->isDirichlet())
        (*it).second->initVector(vector);

    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && (*it).second->isDirichlet())
        (*it).second->initVector(vector);
  }


  void BoundaryManager::exitVector(DOFVectorBase<double>* vector)
  {
    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && !(*it).second->isDirichlet())
        (*it).second->exitVector(vector);

    for (BoundaryIndexMap::iterator it = localBCs.begin(); it != localBCs.end(); ++it)
      if ((*it).second && (*it).second->isDirichlet())
        (*it).second->exitVector(vector);
  }

} // nd namespace AMDiS
