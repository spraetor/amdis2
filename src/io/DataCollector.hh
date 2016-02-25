#pragma once

#include "Traverse.hpp"
#include "DOFVector.hpp"
#include "SurfaceRegion_ED.hpp"
#include "ElementRegion_ED.hpp"
#include "Projection.hpp"

namespace AMDiS
{

  template<typename T>
  DataCollector<T>::DataCollector(const FiniteElemSpace* fe,
                                  DOFVector<T>* val,
                                  int l,
                                  Flag flag,
                                  bool (*writeElemFct)(ElInfo*))
    : values(val),
      level(l),
      traverseFlag(flag),
      feSpace(fe),
      nVertices(0),
      nElements(0),
      nInterpPoints(0),
      nConnection(0),
      elementDataCollected(false),
      valueDataCollected(false),
      periodicDataCollected(false),
      writeElem(writeElemFct)
  {
    FUNCNAME("DataCollector<T>::DataCollector()");

    TEST_EXIT(feSpace)("No finite elem space defined!\n");

    // === get mesh  ===
    mesh = feSpace->getMesh();

    // === get admin ===
    localAdmin = const_cast<DOFAdmin*>(feSpace->getAdmin());
    // === create vertex info vector ===
    vertexInfos = new DOFVector<std::list<VertexInfo>>(feSpace, "vertex infos");

    interpPointInd = new DOFVector<int>(feSpace, "interpolation point indices");
    interpPointCoords = new DOFVector<std::list<WorldVector<double>>>(feSpace, "interpolation point coordinates");
    dofCoord = new DOFVector<std::list<WorldVector<double>>>(feSpace, "dof coords");

    dim = mesh->getDim();
    nPreDofs = localAdmin->getNumberOfPreDofs(VERTEX);
  }


  template<typename T>
  DataCollector<T>::~DataCollector()
  {
    delete vertexInfos;
    delete interpPointInd;
    delete interpPointCoords;
    delete dofCoord;
  }


  template<typename T>
  void DataCollector<T>::fillAllData()
  {
    if (!elementDataCollected)
      startCollectingElementData();

    if (!periodicDataCollected)
      startCollectingPeriodicData();

    if (!valueDataCollected)
      startCollectingValueData();
  }


  template<typename T>
  void DataCollector<T>::startCollectingElementData()
  {
    Flag flag = traverseFlag;
    flag |=
      Mesh::FILL_NEIGH      |
      Mesh::FILL_COORDS     |
      Mesh::FILL_OPP_COORDS |
      Mesh::FILL_BOUND;

    elements.resize(0, dim);

    TraverseStack stack;

    // Traverse elements to create continuous element indices
    ElInfo* elInfo = stack.traverseFirst(mesh, level, flag);
    while (elInfo)
    {
      if (!writeElem || writeElem(elInfo))
        outputIndices[elInfo->getElement()->getIndex()] = nElements++;
      elInfo = stack.traverseNext(elInfo);
    }

    // Traverse elements to create element information
    elInfo = stack.traverseFirst(mesh, level, flag);

    while (elInfo)
    {
      if (!writeElem || writeElem(elInfo))
        addElementData(elInfo);

      elInfo = stack.traverseNext(elInfo);
    }

    elementDataCollected = true;
  }


  template<typename T>
  void DataCollector<T>::startCollectingValueData()
  {
    DOFVector<int>::Iterator intPointIt(interpPointInd, USED_DOFS);
    for (intPointIt.reset(); !intPointIt.end(); ++intPointIt)
      (*intPointIt) = -1;

    interpPoints.clear();

    basisFcts = const_cast<BasisFunction*>(feSpace->getBasisFcts());
    nBasisFcts = basisFcts->getNumber();

    // Traverse elements to add value information and to mark
    // interpolation points.
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, level,
                                         traverseFlag | Mesh::FILL_COORDS);
    while (elInfo)
    {
      if (!writeElem || writeElem(elInfo))
        addValueData(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }

    // If there are interpolation points, add them to the corresponding
    // data array.
    if (nInterpPoints > 0)
    {
      // Remove all interpolation marks and, instead, set to each
      // interpolation point its continous index starting from 0.
      int i = 0;
      for (intPointIt.reset(); !intPointIt.end(); ++intPointIt)
        if (*intPointIt == -3)
          *intPointIt = i++;

      // Traverse elements to create interpolation values.
      elInfo = stack.traverseFirst(mesh, level, traverseFlag | Mesh::FILL_COORDS);
      while (elInfo)
      {
        if (!writeElem || writeElem(elInfo))
          addInterpData(elInfo);
        elInfo = stack.traverseNext(elInfo);
      }
    }

    valueDataCollected = true;
  }


  template<typename T>
  void DataCollector<T>::startCollectingPeriodicData()
  {
    periodicConnections.clear();

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, level, traverseFlag);
    while (elInfo)
    {
      if (!writeElem || writeElem(elInfo))
      {
        LeafDataPeriodic* ldp = dynamic_cast<LeafDataPeriodic*>
                                (elInfo->getElement()->
                                 getElementData()->
                                 getElementData(PERIODIC));

        if (ldp)
          nConnection += dynamic_cast<LeafDataPeriodic*>(ldp)->getInfoList().size();

        periodicConnections.push_back(DimVec<bool>(dim, false));
      }
      elInfo = stack.traverseNext(elInfo);
    }

    nConnection /= 2;

    periodicInfos.clear();

    Flag flag = traverseFlag;
    flag |=
      Mesh::FILL_COORDS    |
      Mesh::FILL_OPP_COORDS|
      Mesh::FILL_NEIGH     |
      Mesh::FILL_BOUND;

    elInfo = stack.traverseFirst(mesh, level, flag);
    while (elInfo)
    {
      if (!writeElem || writeElem(elInfo))
        addPeriodicData(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }

    periodicDataCollected = true;
  }


  template<typename T>
  void DataCollector<T>::addElementData(ElInfo* elInfo)
  {
    const DegreeOfFreedom** dof = elInfo->getElement()->getDof();

    // create ElementInfo
    ElementInfo elementInfo(dim);

    // read element region
    ElementData* ed = elInfo->getElement()->getElementData(ELEMENT_REGION);

    if (ed)
      elementInfo.elementRegion = dynamic_cast<ElementRegion_ED*>(ed)->getRegion();
    else
      elementInfo.elementRegion = -1;

    // read surface regions to element info
    ed = elInfo->getElement()->getElementData(SURFACE_REGION);
    elementInfo.surfaceRegions.set(-1);
    while (ed)
    {
      SurfaceRegion_ED* sr = dynamic_cast<SurfaceRegion_ED*>(ed);
      elementInfo.surfaceRegions[sr->getSide()] = sr->getRegion();
      ed = ed->getDecorated(SURFACE_REGION);
    }

    // for all vertices
    for (int i = 0; i < mesh->getGeo(VERTEX); i++)
    {
      // get coords of this vertex
      WorldVector<double>& vertexCoords = elInfo->getCoord(i);

      // get dof index of this vertex
      DegreeOfFreedom vertexDof = dof[i][nPreDofs];

      // search for coords at this dof
      std::list<VertexInfo>::iterator it =
        find((*vertexInfos)[vertexDof].begin(), (*vertexInfos)[vertexDof].end(),
             vertexCoords);

      // coords not yet in list?
      if (it == (*vertexInfos)[vertexDof].end())
      {
        // create new vertex info and increment number of vertices
        VertexInfo newVertexInfo = {vertexCoords, nVertices};

        // add new vertex info to list
        (*vertexInfos)[vertexDof].push_front(newVertexInfo);

        // set iterator to new vertex info
        it = (*vertexInfos)[vertexDof].begin();

        nVertices++;
      }

      // fill element info
      elementInfo.vertexInfo[i] = it;
      elementInfo.boundary[i] = elInfo->getBoundary(i);
      elementInfo.projection[i] = elInfo->getProjection(i);
      elementInfo.neighbour[i] =
        elInfo->getNeighbour(i) ?
        outputIndices[elInfo->getNeighbour(i)->getIndex()] :
        -1;
    }

    if (dim == 3)
      elementInfo.type = elInfo->getType();

    // remember element info
    elements.push_back(elementInfo);
  }


  template<typename T>
  void DataCollector<T>::addValueData(ElInfo* elInfo)
  {
    std::vector<DegreeOfFreedom> localDOFs(basisFcts->getNumber());
    basisFcts->getLocalIndices(elInfo->getElement(), localAdmin, localDOFs);

    WorldVector<double> vertexCoords;
    // First, traverse all DOFs at the vertices of the element, determine
    // their coordinates and add them to the corresponding entry in dofCoords.
    for (int i = 0; i < mesh->getGeo(VERTEX); i++)
    {
      DegreeOfFreedom dofi = localDOFs[i];

      (*interpPointInd)[dofi] = -2; // mark as vertex

      // get coords of this vertex
      vertexCoords = elInfo->getCoord(i);

      // search for coords at this dof
      std::list<WorldVector<double>>::iterator it =
                                    find((*dofCoord)[dofi].begin(),
                                         (*dofCoord)[dofi].end(),
                                         vertexCoords);

      // coords not yet in list?
      if (it == (*dofCoord)[dofi].end())
      {
        // add new coords to list
        (*dofCoord)[dofi].push_back(vertexCoords);
      }
    }

    // Then, traverse all interpolation DOFs of the element, determine
    // their coordinates and add them to the corresponding entry in
    // interpPointCoords.
    for (int i = mesh->getGeo(VERTEX); i < nBasisFcts; i++)
    {
      DegreeOfFreedom dofi = localDOFs[i];

      elInfo->coordToWorld(*basisFcts->getCoords(i), vertexCoords);

      if ((*interpPointInd)[dofi] == -1)
      {
        // mark as interpolation point
        (*interpPointInd)[dofi] = -3;

        // search for interpolation point coordinates, and insert them to the
        // dof-entry, if not contained in the list
        std::list<WorldVector<double>>::iterator it =
                                      find((*interpPointCoords)[dofi].begin(),
                                           (*interpPointCoords)[dofi].end(),
                                           vertexCoords);

        if (it == (*interpPointCoords)[dofi].end())
        {
          (*interpPointCoords)[dofi].push_back(vertexCoords);
          nInterpPoints++;
        }
      }
    }
  }


  template<typename T>
  void DataCollector<T>::addInterpData(ElInfo* elInfo)
  {
    std::vector<DegreeOfFreedom> localDOFs(basisFcts->getNumber());
    basisFcts->getLocalIndices(elInfo->getElement(), localAdmin, localDOFs);

    std::vector<int> elemInterpPoints(0);
    for (int i = mesh->getGeo(VERTEX); i < nBasisFcts; i++)
      elemInterpPoints.push_back((*interpPointInd)[localDOFs[i]]);

    interpPoints.push_back(elemInterpPoints);
  }


  template<typename T>
  void DataCollector<T>::addPeriodicData(ElInfo* elInfo)
  {
    LeafDataPeriodic* ldp = dynamic_cast<LeafDataPeriodic*>
                            (elInfo->getElement()->
                             getElementData()->
                             getElementData(PERIODIC));

    if (ldp)
    {
      std::list<LeafDataPeriodic::PeriodicInfo>::iterator it;

      for (it = dynamic_cast<LeafDataPeriodic*>(ldp)->getInfoList().begin();
           it != dynamic_cast<LeafDataPeriodic*>(ldp)->getInfoList().end();
           ++it)
      {

        int outputIndex = outputIndices[elInfo->getElement()->getIndex()];
        int neighIndex  = outputIndices[elInfo->
                                        getNeighbour(it->elementSide)->
                                        getIndex()];

        if (!periodicConnections[outputIndex][it->elementSide])
        {
          PeriodicInfo periodicInfo;

          periodicInfo.mode = it->periodicMode;
          periodicInfo.type = it->type;
          periodicInfo.outputIndex = outputIndex;
          periodicInfo.neighIndex = neighIndex;
          periodicInfo.vertexMap.clear();

          periodicConnections[outputIndex][it->elementSide] = true;
          periodicConnections
          [neighIndex][elInfo->getOppVertex(it->elementSide)] = true;


          for (int i = 0; i < dim; i++)
          {
            int index1 = elInfo->getElement()->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim),
                         it->elementSide,
                         i);
            DegreeOfFreedom dof1 = elInfo->getElement()->getDof(index1, nPreDofs);

            for (int j = 0; j < dim; j++)
            {
              int index2 = elInfo->getElement()->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim),
                           elInfo->getOppVertex(it->elementSide),
                           j);
              DegreeOfFreedom dof2 = elInfo->getNeighbour(it->elementSide)->getDof(index2, nPreDofs);

              if ((dof1 == dof2) || (mesh->associated(dof1, dof2)))
              {
                periodicInfo.vertexMap[index1] = index2;
                break;
              }
            }
          }

          periodicInfos.push_back(periodicInfo);
        }
      }
    }
  }


  template<typename T>
  std::list<ElementInfo>* DataCollector<T>::getElementInfos()
  {
    if (!elementDataCollected)
      startCollectingElementData();

    return &elements;
  }


  template<typename T>
  DOFVector<std::list<VertexInfo>>* DataCollector<T>::getVertexInfos()
  {
    if (!elementDataCollected)
      startCollectingElementData();

    return vertexInfos;
  }


  template<typename T>
  std::list<PeriodicInfo>* DataCollector<T>::getPeriodicInfos()
  {
    if (!periodicDataCollected)
      startCollectingPeriodicData();

    return &periodicInfos;
  }


  template<typename T>
  int DataCollector<T>::getNumberVertices()
  {
    if (!elementDataCollected)
      startCollectingElementData();

    return nVertices;
  }


  template<typename T>
  int DataCollector<T>::getNumberElements()
  {
    if (!elementDataCollected)
      startCollectingElementData();

    return nElements;
  }


  template<typename T>
  int DataCollector<T>::getNumberInterpPoints()
  {
    if (!valueDataCollected)
      startCollectingValueData();

    return nInterpPoints;
  }


  template<typename T>
  int DataCollector<T>::getNumberConnections()
  {
    if (!periodicDataCollected)
      startCollectingPeriodicData();

    return nConnection;
  }


  template<typename T>
  DOFVector<std::list<WorldVector<double>>>* DataCollector<T>::getDofCoords()
  {
    if (!valueDataCollected)
      startCollectingValueData();

    return dofCoord;
  }


  template<typename T>
  DOFVector<int>* DataCollector<T>::getInterpPointInd()
  {
    if (!valueDataCollected)
      startCollectingValueData();

    return interpPointInd;
  }


  template<typename T>
  DOFVector<std::list<WorldVector<double>>>* DataCollector<T>::getInterpPointCoords()
  {
    if (!valueDataCollected)
      startCollectingValueData();

    return interpPointCoords;
  }


  template<typename T>
  std::vector<std::vector<int>>* DataCollector<T>::getInterpPoints()
  {
    if (!valueDataCollected)
      startCollectingValueData();

    return &interpPoints;
  }
}
