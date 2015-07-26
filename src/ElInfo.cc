#include "ElInfo.h"
#include "BasisFunction.h"
#include "Element.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "FiniteElemSpace.h"
#include "Flag.h"
#include "MacroElement.h"
#include "Mesh.h"
#include "Global.h"
#include "FixVec.h"
#include "DOFVector.h"

namespace AMDiS 
{

  std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double> > > ElInfo::subElemMatrices(5);

  std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double> > > ElInfo::subElemGradMatrices(5);

  ElInfo::ElInfo(Mesh *aMesh) 
    : mesh(aMesh),
      element(NULL),
      parent(NULL),
      macroElement(NULL),
      level(0),
      elType(0),
      iChild(0),
      coord(mesh->getDim(), NO_INIT),
      boundary(mesh->getDim(), DEFAULT_VALUE, INTERIOR),
      projection(mesh->getDim(), NO_INIT),
      oppCoord(mesh->getDim(), NO_INIT),
      neighbour(mesh->getDim(), NO_INIT),
      neighbourCoord(mesh->getDim(), NO_INIT),
      oppVertex(mesh->getDim(), NO_INIT),
      grdLambda(mesh->getDim(), NO_INIT),
      refinementPath(0),
      refinementPathLength(0)
  {
    projection.set((Projection*)(NULL));

    for (int i = 0; i < neighbourCoord.getSize(); i++)
      neighbourCoord[i].init(mesh->getDim());

    dimOfWorld = Global::getGeo(WORLD);
  }
 

  ElInfo::~ElInfo()
  {}


  void ElInfo::coordToWorld(const DimVec<double>& l, 
			    WorldVector<double>& w) const

  {
    testFlag(Mesh::FILL_COORDS);

    double c = l[0];

    for (int j = 0; j < dimOfWorld; j++)
      w[j] = c * coord[0][j];
  
    int vertices = Global::getGeo(VERTEX, l.getSize() - 1);

    for (int i = 1; i < vertices; i++) {
      c = l[i];
      for (int j = 0; j < dimOfWorld; j++)
	w[j] += c * coord[i][j];
    }
  }

  double ElInfo::calcDet() const
  {
    testFlag(Mesh::FILL_COORDS);
    return calcDet(coord);
  }

  double ElInfo::calcDet(const FixVec<WorldVector<double>, VERTEX> &coords) const
  {
    FUNCNAME("ElInfo::calcDet()");

    int dim = coords.getSize() - 1;
    int dow = Global::getGeo(WORLD);

    TEST_EXIT_DBG(dim <= dow)("dim > dow\n");

    double det = 0.0;

    if (dim == 0)
      return 1.0;

    switch (dow) {
    case 1:
      det = coords[1][0] - coords[0][0];
      break;
    case 2:
      if (dim == 1) {
      
      	WorldVector<double> e1;
      
      	e1[0] = coords[1][0] - coords[0][0];
      	e1[1] = coords[1][1] - coords[0][1];
      	det = norm(&e1);

      } else {
      	det = (coords[1][0] - coords[0][0]) * (coords[2][1] - coords[0][1]) - 
      	  (coords[1][1] - coords[0][1]) * (coords[2][0] - coords[0][0]);
      }
      break;
    case 3: 
      {
      	WorldVector<double> e1, e2, e3, n;
      	
      	for (int i = 0; i < dow; i++) {
      	  e1[i] = coords[1][i] - coords[0][i];
      	  if (dim > 1)
      	    e2[i] = coords[2][i] - coords[0][i];
      	  if (dim > 2)
      	    e3[i] = coords[3][i] - coords[0][i];
      	}
        
      	if (dim > 1) {
      	  n[0] = e1[1] * e2[2] - e1[2] * e2[1];
      	  n[1] = e1[2] * e2[0] - e1[0] * e2[2];
      	  n[2] = e1[0] * e2[1] - e1[1] * e2[0];
      	}
      	
      	if (dim == 1) {
      	  det = norm(&e1);
      	} else if (dim == 2) {
      	  det = norm(&n);
      	} else if (dim == 3) {
      	  det = n[0] * e3[0] + n[1] * e3[1] + n[2] * e3[2];
      	} else
      	  ERROR_EXIT("not yet for problem dimension = %d", dim);
      	break;
      } 
    default:
      ERROR_EXIT("not yet for Global::getGeo(WORLD) = %d", dow);
    }
    
    return math::abs(det);
  }


  double ElInfo::calcSurfaceDet(VectorOfFixVecs<DimVec<double> > &surfVert) const
  {
    double surfDet;
    int dim = surfVert[0].getSize() - 1;
    FixVec<WorldVector<double>, VERTEX> worldCoords(dim - 1, NO_INIT);
  
    // transform barycentric coords to world coords
    for (int i = 0; i < dim; i++) {
      coordToWorld(surfVert[i], worldCoords[i]);
    }
  
    // calculate determinant for surface
    surfDet = calcDet(worldCoords);

    return surfDet;
  }


  void ElInfo::fillDetGrdLambda() 
  { 
    if (fillFlag.isSet(Mesh::FILL_GRD_LAMBDA)) {
      det = calcGrdLambda(grdLambda);
    } else {
      if (fillFlag.isSet(Mesh::FILL_DET))
	det = calcDet();
    }
  }


  void ElInfo::fillElInfo(const MacroElement *mel, int refinementPathLength, unsigned long refinementPath)
  {
    if (refinementPathLength == 0) {
      fillMacroInfo(mel);
      return;
    }
    
    std::vector<ElInfo*> elInfo;
    elInfo.push_back(mesh->createNewElInfo());
    elInfo.push_back(this);

    elInfo[0]->setFillFlag(getFillFlag());
    elInfo[0]->fillMacroInfo(mel);

    int i = 0;
    for (; i < refinementPathLength; i++) {
      elInfo[(i+1)%2]->fillElInfo(static_cast<int>((refinementPath & (1<<i)) == static_cast<unsigned long>(1<<i)), elInfo[i%2]);
    }
    if (i%2 == 0)
      *this = *elInfo[0];

    delete elInfo[0];
  }


  BoundaryType ElInfo::getBoundary(GeoIndex pos, int i)
  {
    static int indexOffset[3][3] = {
      { 0, 0, 0},
      { 3, 0, 0},
      {10, 4, 0}
    };

    int dim = mesh->getDim();
    int posIndex = DIM_OF_INDEX(pos, dim);
    int offset = indexOffset[dim - 1][posIndex];
    
    return boundary[offset + i];
  }

} // end namespace AMDiS
