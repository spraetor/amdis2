#include "ElInfo.hpp"
#include "BasisFunction.hpp"
#include "Element.hpp"
#include "Line.hpp"
#include "Triangle.hpp"
#include "Tetrahedron.hpp"
#include "FiniteElemSpace.hpp"
#include "Flag.hpp"
#include "MacroElement.hpp"
#include "Mesh.hpp"
#include "Global.hpp"
#include "FixVec.hpp"
#include "DOFVector.hpp"

namespace AMDiS
{

  std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double>>> ElInfo::subElemMatrices(5);

  std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double>>> ElInfo::subElemGradMatrices(5);

  ElInfo::ElInfo(Mesh* aMesh)
    : mesh(aMesh),
      element(NULL),
      parent(NULL),
      macroElement(NULL),
      level(0),
      elType(0),
      iChild(0),
      coord(mesh->getDim()),
      boundary(mesh->getDim(), INTERIOR),
      projection(mesh->getDim()),
      oppCoord(mesh->getDim()),
      neighbour(mesh->getDim()),
      neighbourCoord(mesh->getDim(), FixVec<WorldVector<double>, VERTEX>(mesh->getDim())),
      oppVertex(mesh->getDim()),
      grdLambda(mesh->getDim()),
      refinementPath(0),
      refinementPathLength(0)
  {
    projection.set((Projection*)(NULL));

    dimOfWorld = Global::getGeo(WORLD);
  }



  ElInfo& ElInfo::operator=(ElInfo const& rhs)
  {
    mesh = rhs.mesh;
    element = rhs.element;
    parent = rhs.parent;
    macroElement = rhs.macroElement;
    fillFlag = rhs.fillFlag;
    level = rhs.level;
    elType = rhs.elType;
    iChild = rhs.iChild;
    coord = rhs.coord;
    boundary = rhs.boundary;
    projection = rhs.projection;
    oppCoord = rhs.oppCoord;
    neighbour = rhs.neighbour;
    neighbourCoord = rhs.neighbourCoord;
    oppVertex = rhs.oppVertex;
    det = rhs.det;
    grdLambda = rhs.grdLambda;
    parametric = rhs.parametric;
    dimOfWorld = rhs.dimOfWorld;

    return *this;
  }


  void ElInfo::getGrdLambda(mtl::dense2D<double>& grd_lam)
  {
    grd_lam.change_dim(grdLambda.getSize(), Global::getGeo(WORLD));
    for (size_t i = 0; i < static_cast<size_t>(grdLambda.getSize()); i++)
      for (size_t j = 0; j < static_cast<size_t>(Global::getGeo(WORLD)); j++)
        grd_lam(i,j) = grd_lam[i][j];
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

    for (int i = 1; i < vertices; i++)
    {
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

  double ElInfo::calcDet(const FixVec<WorldVector<double>, VERTEX>& coords) const
  {
    FUNCNAME("ElInfo::calcDet()");

    int dim = coords.getSize() - 1;
    int dow = Global::getGeo(WORLD);

    TEST_EXIT_DBG(dim <= dow)("dim > dow\n");

    double det = 0.0;

    if (dim == 0)
      return 1.0;

    switch (dow)
    {
    case 1:
      det = coords[1][0] - coords[0][0];
      break;
    case 2:
      if (dim == 1)
      {

        WorldVector<double> e1;

        e1[0] = coords[1][0] - coords[0][0];
        e1[1] = coords[1][1] - coords[0][1];
        det = norm(e1);

      }
      else
      {
        det = (coords[1][0] - coords[0][0]) * (coords[2][1] - coords[0][1]) -
              (coords[1][1] - coords[0][1]) * (coords[2][0] - coords[0][0]);
      }
      break;
    case 3:
    {
      WorldVector<double> e1, e2, e3, n;

      for (int i = 0; i < dow; i++)
      {
        e1[i] = coords[1][i] - coords[0][i];
        if (dim > 1)
          e2[i] = coords[2][i] - coords[0][i];
        if (dim > 2)
          e3[i] = coords[3][i] - coords[0][i];
      }

      if (dim > 1)
      {
        n[0] = e1[1] * e2[2] - e1[2] * e2[1];
        n[1] = e1[2] * e2[0] - e1[0] * e2[2];
        n[2] = e1[0] * e2[1] - e1[1] * e2[0];
      }

      if (dim == 1)
      {
        det = norm(e1);
      }
      else if (dim == 2)
      {
        det = norm(n);
      }
      else if (dim == 3)
      {
        det = n[0] * e3[0] + n[1] * e3[1] + n[2] * e3[2];
      }
      else
        ERROR_EXIT("not yet for problem dimension = %d", dim);
      break;
    }
    default:
      ERROR_EXIT("not yet for Global::getGeo(WORLD) = %d", dow);
    }

    return math::abs(det);
  }


  double ElInfo::calcSurfaceDet(VectorOfFixVecs<DimVec<double>>& surfVert) const
  {
    double surfDet;
    int dim = surfVert[0].getSize() - 1;
    FixVec<WorldVector<double>, VERTEX> worldCoords(dim - 1);

    // transform barycentric coords to world coords
    for (int i = 0; i < dim; i++)
    {
      coordToWorld(surfVert[i], worldCoords[i]);
    }

    // calculate determinant for surface
    surfDet = calcDet(worldCoords);

    return surfDet;
  }


  void ElInfo::fillDetGrdLambda()
  {
    if (fillFlag.isSet(Mesh::FILL_GRD_LAMBDA))
    {
      det = calcGrdLambda(grdLambda);
    }
    else
    {
      if (fillFlag.isSet(Mesh::FILL_DET))
        det = calcDet();
    }
  }


  void ElInfo::fillElInfo(const MacroElement* mel, int refinementPathLength, unsigned long refinementPath)
  {
    if (refinementPathLength == 0)
    {
      fillMacroInfo(mel);
      return;
    }

    std::vector<ElInfo*> elInfo;
    elInfo.push_back(mesh->createNewElInfo());
    elInfo.push_back(this);

    elInfo[0]->setFillFlag(getFillFlag());
    elInfo[0]->fillMacroInfo(mel);

    int i = 0;
    for (; i < refinementPathLength; i++)
    {
      elInfo[(i+1)%2]->fillElInfo(static_cast<int>((refinementPath & (1<<i)) == static_cast<unsigned long>(1<<i)), elInfo[i%2]);
    }
    if (i%2 == 0)
      *this = *elInfo[0];

    delete elInfo[0];
  }


  BoundaryType ElInfo::getBoundary(GeoIndex pos, int i)
  {
    static int indexOffset[3][3] =
    {
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
