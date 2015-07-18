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



/** \file DataCollector.h */

#ifndef AMDIS_DATACOLLECTOR_H
#define AMDIS_DATACOLLECTOR_H

#include <list>
#include <vector>
#include "Global.h"
#include "ElementInfo.h"
#include "VertexInfo.h"
#include "PeriodicInfo.h"
#include "Mesh.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  /** 
   * \ingroup Output
   *
   * \brief
   * Provides data collecting of element, vertex and value data for file writer.
   */
  template<typename T=double>
  class DataCollector 
  {
  public:
    /// Constructor
    DataCollector(const FiniteElemSpace *feSpace,
		  DOFVector<T> *values = NULL,
		  int level = -1,
		  Flag traverseFlag = Mesh::CALL_LEAF_EL,
		  bool (*writeElem)(ElInfo*) = NULL);

    ~DataCollector();

    /// Fills the DataCollector with all possible datas.
    void fillAllData();

    /// Returns list with element information.
    std::list<ElementInfo>* getElementInfos();

    /// Returns vector with vertex information.
    DOFVector<std::list<VertexInfo> >* getVertexInfos();

    /// Returns the finite element space of the problem.
    const FiniteElemSpace* getFeSpace()
    {
      return feSpace;
    }

    /// Returns vector with value information.
    DOFVector<T>* getValues()
    {
      if (!valueDataCollected)
	startCollectingValueData();

      return values;
    }

    /// Returns vector with information about dof coordinates.
    DOFVector< std::list<WorldVector<double> > >* getDofCoords();

    /// Returns vector with information about interpolation point indexing.
    DOFVector<int>* getInterpPointInd();

    ///
    DOFVector< std::list<WorldVector<double> > >* getInterpPointCoords();

    /// Returns list of interpolation point information.
    std::vector< std::vector<int> >* getInterpPoints();

    /// Returns list of information about periodics.
    std::list<PeriodicInfo>* getPeriodicInfos();

    /// Returns the number of vertices.
    int getNumberVertices();

    /// Returns the number of elements.
    int getNumberElements();

    /// Returns the number of interpolation points.
    int getNumberInterpPoints();

    // Returns the number of connections.
    int getNumberConnections();
      
    /// Returns the mesh of the problem.
    Mesh* getMesh()
    {
      return mesh;
    }

    void setMesh(Mesh *m) 
    {
      mesh = m;
    }

  protected:  
    /// Start collecting element and vertex data of the problem.
    void startCollectingElementData();

    /// Start collecting value data of the problem.
    void startCollectingValueData();

    /// Start collecting periodic data of the problem.
    void startCollectingPeriodicData();

    /// Adds information about one element and its vertices.
    void addElementData(ElInfo* elInfo);

    /// Adds value information of one element.
    void addValueData(ElInfo *elInfo);

    /// Adds information about interpolation points of vertices.
    void addInterpData(ElInfo *elInfo);

    /// Adds value information of one element.
    void addPeriodicData(ElInfo *elInfo);

    /// Vector with vertex values
    DOFVector<T> *values;

    /// Level information for traversing the mesh.
    int level;

    /// Flags for traversing the mesh.
    Flag traverseFlag;

    ///
    const FiniteElemSpace *feSpace;

    /// Mesh that should be written
    Mesh *mesh;
      
    /// DOFAdmin of values
    DOFAdmin *localAdmin;

    /// vertex pre-dofs
    int nPreDofs;

    /// Number of vertices.
    int nVertices;

    /// Number of elements.
    int nElements;

    /// Total number of interpolation points.
    int nInterpPoints;

    /// Number of connections in periodic problems.
    int nConnection;

    /// Dimension of \ref mesh
    int dim;
           
    /// Maps internal element indices to global output indices.
    std::map<int, int> outputIndices;

    /// Global interpolation point indexing
    DOFVector<int> *interpPointInd;

    /// Stores for each simplex the interpolation points.
    std::vector<std::vector<int> > interpPoints;

    /// Stores for each DOF a list of its coordinates. If there are now periodic
    /// boundaries than there is also only one coordinate per DOF.
    DOFVector<std::list<WorldVector<double> > > *interpPointCoords;

    /// list of coords for each dof
    DOFVector<std::list<WorldVector<double> > > *dofCoord;

    /// List that stores an ElementInfo for each element.
    std::list<ElementInfo> elements;

    /// List stat stores information about all periodics.
    std::list<PeriodicInfo> periodicInfos;

    /// Stores a list of vertex infos for each dof.
    DOFVector<std::list<VertexInfo> > *vertexInfos;

    /// periodicConnections[i][j] stores whether the connection at side j of 
    /// the element with output index i has already been written.
    std::vector<DimVec<bool> > periodicConnections;

    /// Stores if element data was collected before.
    bool elementDataCollected;

    /// Stores if value data was collected before.   
    bool valueDataCollected;

    /// Stores if periodic data was collected before.
    bool periodicDataCollected;

    /// Pointer to a function which decides whether an element is considered.
    bool (*writeElem)(ElInfo*);

    /// Temporary variable used in functions addValueData() and addInterpData().
    BasisFunction *basisFcts;

    /// Temporary variable used in functions addValueData() and addInterpData().
    int nBasisFcts;

    /// Temporary variable used in function \ref addValueData.
    WorldVector<double> *vertexCoords;
  };
}

#include "DataCollector.hh"

#endif

