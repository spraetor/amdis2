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


#include <map>
#include <fstream>
#include "Global.h"
#include "parallel/ElementObjectDatabase.h"
#include "parallel/PeriodicMap.h"
#include "parallel/ParallelTypes.h"

namespace AMDiS { namespace Parallel {

  using namespace std;

  void PeriodicMap::add(const FiniteElemSpace *feSpace,
			PeriodicDofMap &newMap)
  {
    for (PeriodicDofMap::iterator it = newMap.begin(); it != newMap.end(); ++it) {
      for (std::map<DegreeOfFreedom, DegreeOfFreedom>::iterator dofIt = it->second.begin();
	   dofIt != it->second.end(); ++dofIt)
	add(feSpace, it->first, dofIt->second, dofIt->first);
    }
  }


  void PeriodicMap::fillAssociations(const FiniteElemSpace* feSpace,
				     int globalDofIndex,
				     const ElementObjectDatabase &elObjDb,
				     std::set<int>& perAsc)
  {
    if (isPeriodic(feSpace, globalDofIndex) == false)
      return;

    std::set<int>& tmpAsc = getAssociations(feSpace, globalDofIndex);
    for (std::set<int>::iterator it = tmpAsc.begin(); it != tmpAsc.end(); ++it)
      if (elObjDb.isValidPeriodicType(*it))
	perAsc.insert(*it);
  }


  void PeriodicMap::mapDof(const FiniteElemSpace* feSpace,
			   int globalDofIndex, 
			   const std::set<int>& perAsc, 
			   vector<int>& mappedDofs)
  {
    mappedDofs.clear();
    mappedDofs.push_back(globalDofIndex);
    
    for (std::set<int>::const_iterator it = perAsc.begin();
	 it != perAsc.end(); ++it) {
      int nDofs = static_cast<int>(mappedDofs.size());
      
      for (int i = 0; i < nDofs; i++) 
	if (isPeriodic(feSpace, *it, mappedDofs[i]))
	  mappedDofs.push_back(map(feSpace, *it, mappedDofs[i]));
    }
  }


  void PeriodicMap::mapDof(const FiniteElemSpace* rowFeSpace,
			   const FiniteElemSpace* colFeSpace,
			   pair<int, int> globalDofIndex, 
			   const std::set<int>& perAsc, 
			   vector<pair<int, int> >& mappedDofs)
  {
    mappedDofs.clear();
    mappedDofs.push_back(globalDofIndex);

    for (std::set<int>::iterator it = perAsc.begin(); 
	 it != perAsc.end(); ++it) {
      int nDofs = static_cast<int>(mappedDofs.size());

      for (int i = 0; i < nDofs; i++) {
	int perRowDof = 0;	
	if (isPeriodic(rowFeSpace, *it, mappedDofs[i].first))
	  perRowDof = map(rowFeSpace, *it, mappedDofs[i].first);
	else
	  perRowDof = mappedDofs[i].first;
	
	int perColDof;
	if (isPeriodic(colFeSpace, *it, mappedDofs[i].second))
	  perColDof = map(colFeSpace, *it, mappedDofs[i].second);
	else
	  perColDof = mappedDofs[i].second;	      	      
	
	mappedDofs.push_back(make_pair(perRowDof, perColDof));
      }
    }    
  }


  void PeriodicMap::serialize(ostream &out, 
			      vector<const FiniteElemSpace*> feSpaces)
  {
    int nFeSpace = static_cast<int>(feSpaces.size());
    
    for (int i = 0; i < nFeSpace; i++)
      serialize(out, periodicDofMap[feSpaces[i]]);
    
    for (int i = 0; i < nFeSpace; i++)
      serialize(out, periodicDofAssociations[feSpaces[i]]);
  }


  void PeriodicMap::deserialize(istream &in, 
				vector<const FiniteElemSpace*> feSpaces)
  {
    int nFeSpace = static_cast<int>(feSpaces.size());
    
    for (int i = 0; i < nFeSpace; i++)
      deserialize(in, periodicDofMap[feSpaces[i]]);
    
    for (int i = 0; i < nFeSpace; i++)
      deserialize(in, periodicDofAssociations[feSpaces[i]]);
  }
  
  
  void PeriodicMap::serialize(ostream &out, PeriodicDofMap &data)
  {
    int mapSize = data.size();
    SerUtil::serialize(out, mapSize);
    
    for (PeriodicDofMap::iterator it = data.begin(); it != data.end(); ++it) {
      int type = it->first;
      std::map<DegreeOfFreedom, DegreeOfFreedom> dofMap = it->second;
      
      SerUtil::serialize(out, type);
      SerUtil::serialize(out, dofMap);
    }
  }
  

  void PeriodicMap::serialize(ostream &out, std::map<int, std::set<int> >& data)
  {
    int mapSize = data.size();
    SerUtil::serialize(out, mapSize);
    
    for (std::map<int, std::set<int> >::iterator it = data.begin(); 
	 it != data.end(); ++it) {
      int dof = it->first;
      std::set<int> typeSet = it->second;
      
      SerUtil::serialize(out, dof);
      SerUtil::serialize(out, typeSet);
    }
  }
  
  
  void PeriodicMap::deserialize(istream &in, PeriodicDofMap &data)
  {
    data.clear();
    
    int mapSize = 0;
    SerUtil::deserialize(in, mapSize);
    
    for (int i = 0; i < mapSize; i++) {
      int type;
      std::map<DegreeOfFreedom, DegreeOfFreedom> dofMap;
      
      SerUtil::deserialize(in, type);
      SerUtil::deserialize(in, dofMap);
      
      data[type] = dofMap;
    }
  }
  
  
  void PeriodicMap::deserialize(istream &in, std::map<int, std::set<int> >& data)
  {
    data.clear();
    
    int mapSize = 0;
    SerUtil::deserialize(in, mapSize);
    
    for (int i = 0; i < mapSize; i++) {
      int dof;
      std::set<int> typeSet;
      
      SerUtil::deserialize(in, dof);
      SerUtil::deserialize(in, typeSet);
      
      data[dof] = typeSet;
    }
  }
  
} }
