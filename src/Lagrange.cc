#include <stdio.h>
#include <algorithm>
#include <list>

#include "Mesh.h"
#include "Element.h"
#include "Lagrange.h"
#include "DOFAdmin.h"
#include "RCNeighbourList.h"
#include "DOFVector.h"
#include "Traverse.h"
#include "Line.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "Parametric.h"
#include "Debug.h"

using namespace std;

namespace AMDiS
{

  std::vector<DimVec<double>*> Lagrange::baryDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  DimVec<int>* Lagrange::ndofDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  int Lagrange::nBasFctsDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  std::vector<BasFctType*> Lagrange::phiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  std::vector<GrdBasFctType*> Lagrange::grdPhiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  std::vector<D2BasFctType*> Lagrange::D2PhiDimDegree[MAX_DIM + 1][MAX_DEGREE + 1];
  std::list<Lagrange*> Lagrange::allBasFcts;


  Lagrange::Lagrange(int dim, int degree)
    : BasisFunction(std::string("Lagrange"), dim, degree)
  {
    // set name
    name += std::to_string(dim) + " " + std::to_string(degree);

    // set nDOF
    setNDOF();

    // set barycentric coordinates
    setBary();

    // set function pointer
    setFunctionPointer();
  }


  Lagrange::~Lagrange()
  {
    for (int i = 0; i < static_cast<int>(bary->size()); i++)
      if ((*bary)[i])
      {
        delete (*bary)[i];
        (*bary)[i] = NULL;
      }
  }


  Lagrange* Lagrange::getLagrange(int dim, int degree)
  {
    std::list<Lagrange*>::iterator it;
    for (it = allBasFcts.begin(); it != allBasFcts.end(); it++)
      if (((*it)->dim == dim) && ((*it)->degree == degree))
        return (*it);

    Lagrange* newLagrange = new Lagrange(dim, degree);
    allBasFcts.push_back(newLagrange);
    return newLagrange;
  }


  void Lagrange::clear()
  {
    for (std::list<Lagrange*>::iterator it = allBasFcts.begin();
         it != allBasFcts.end(); it++)
      if (*it)
      {
        delete *it;
        *it = NULL;
      }
  }


  void Lagrange::setFunctionPointer()
  {
    if (static_cast<int>(phiDimDegree[dim][degree].size()) == 0)
    {
      // for all positions
      for (int i = 0; i < dim + 1; i++)
      {
        // no vertex dofs for degree 0 ?
        if (degree == 0 && i != dim)
          continue;
        // for all vertices/edges/...
        for (int j = 0; j < Global::getGeo(INDEX_OF_DIM(i, dim),dim); j++)
        {
          // for all dofs at this position
          for (int k = 0; k < (*nDOF)[INDEX_OF_DIM(i, dim)]; k++)
          {
            // basis function
            phiDimDegree[dim][degree].push_back(new Phi(this,
                                                INDEX_OF_DIM(i ,dim), j, k));
            // gradients
            grdPhiDimDegree[dim][degree].push_back(new GrdPhi(this,
                                                   INDEX_OF_DIM(i, dim),
                                                   j, k));
            // D2-Matrices
            D2PhiDimDegree[dim][degree].push_back(new D2Phi(this,
                                                  INDEX_OF_DIM(i, dim),
                                                  j, k));
          }
        }
      }
    }
    phi = &phiDimDegree[dim][degree];
    grdPhi = &grdPhiDimDegree[dim][degree];
    d2Phi = &D2PhiDimDegree[dim][degree];

    switch (degree)
    {
    case 0:
      refineInter_fct = refineInter0;
      coarseRestr_fct = coarseRestr0;
      coarseInter_fct = coarseInter0;  // not yet implemented
      break;
    case 1:
      refineInter_fct = refineInter1;
      coarseRestr_fct = coarseRestr1;
      coarseInter_fct = NULL;  // not yet implemented
      break;
    case 2:
      switch (dim)
      {
      case 1:
        refineInter_fct = refineInter2_1d;
        coarseRestr_fct = NULL; // not yet implemented
        coarseInter_fct = coarseInter2_1d;
        break;
      case 2:
        refineInter_fct = refineInter2_2d;
        coarseRestr_fct = coarseRestr2_2d;
        coarseInter_fct = coarseInter2_2d;
        break;
      case 3:
        refineInter_fct = refineInter2_3d;
        coarseRestr_fct = coarseRestr2_3d;
        coarseInter_fct = coarseInter2_3d;
        break;
      default:
        ERROR_EXIT("invalid dim\n");
      }
      break;
    case 3:
      switch (dim)
      {
      case 1:
        refineInter_fct = refineInter3_1d;
        coarseRestr_fct = coarseRestr3_1d;
        coarseInter_fct = coarseInter3_1d;
        break;
      case 2:
        refineInter_fct = refineInter3_2d;
        coarseRestr_fct = coarseRestr3_2d;
        coarseInter_fct = coarseInter3_2d;
        break;
      case 3:
        refineInter_fct = refineInter3_3d;
        coarseRestr_fct = coarseRestr3_3d;
        coarseInter_fct = coarseInter3_3d;
        break;
      default:
        ERROR_EXIT("invalid dim\n");
      }
      break;
    case 4:
      switch (dim)
      {
      case 1:
        refineInter_fct = refineInter4_1d;
        coarseRestr_fct = coarseRestr4_1d;
        coarseInter_fct = coarseInter4_1d;
        break;
      case 2:
        refineInter_fct = refineInter4_2d;
        coarseRestr_fct = coarseRestr4_2d;
        coarseInter_fct = coarseInter4_2d;
        break;
      case 3:
        refineInter_fct = refineInter4_3d;
        coarseRestr_fct = coarseRestr4_3d;
        coarseInter_fct = coarseInter4_3d;
        break;
      default:
        ERROR_EXIT("invalid dim\n");
      }
      break;
    default:
      ERROR_EXIT("invalid degree\n");
    }
  }


  void Lagrange::setNDOF()
  {
    if (static_cast<int>(baryDimDegree[dim][degree].size()) == 0)
    {
      ndofDimDegree[dim][degree] = new DimVec<int>(dim, 0);

      if (degree != 0)
        (*ndofDimDegree[dim][degree])[VERTEX] = 1;
      else
        (*ndofDimDegree[dim][degree])[VERTEX] = 0;

      for (int i = 1; i < dim + 1; i++)
      {
        nBasFcts = getNumberOfDofs(i, degree);
        (*ndofDimDegree[dim][degree])[INDEX_OF_DIM(i, dim)] = nBasFcts;
        for (int j = 0; j < i; j++)
        {
          (*ndofDimDegree[dim][degree])[INDEX_OF_DIM(i, dim)] -=
            Global::getGeo(INDEX_OF_DIM(j, dim), i) *
            (*ndofDimDegree[dim][degree])[INDEX_OF_DIM(j, dim)];
        }
      }
      nBasFctsDimDegree[dim][degree] = nBasFcts;
    }
    nBasFcts = nBasFctsDimDegree[dim][degree];
    nDOF = ndofDimDegree[dim][degree];
  }

  DimVec<double>* Lagrange::getCoords(int i) const
  {
    return (*bary)[i];
  }

  void Lagrange::setVertices(int dim, int degree,
                             GeoIndex position, int positionIndex, int nodeIndex,
                             int** vertices)
  {
    FUNCNAME_DBG("Lagrange::setVertices()");

    TEST_EXIT_DBG(*vertices == NULL)("vertices != NULL\n");

    int dimOfPosition = DIM_OF_INDEX(position, dim);

    *vertices = new int[dimOfPosition + 1];

    if (degree == 4 && dimOfPosition == 1)
    {
      (*vertices)[(nodeIndex != 2) ? 0 : 1] =
        Global::getReferenceElement(dim)->getVertexOfPosition(position,
            positionIndex,
            0);
      (*vertices)[(nodeIndex != 2) ? 1 : 0] =
        Global::getReferenceElement(dim)->getVertexOfPosition(position,
            positionIndex,
            1);
    }
    else if (degree==4 && dimOfPosition == 2)
    {
      for (int i = 0; i < dimOfPosition + 1; i++)
      {
        (*vertices)[(i + dimOfPosition*nodeIndex) % (dimOfPosition + 1)] =
          Global::getReferenceElement(dim)->getVertexOfPosition(position,
              positionIndex,
              i);
      }
    }
    else
    {
      for (int i = 0; i < dimOfPosition + 1; i++)
      {
        (*vertices)[(i + nodeIndex) % (dimOfPosition + 1)] =
          Global::getReferenceElement(dim)->getVertexOfPosition(position,
              positionIndex,
              i);
      }
    }
  }

  Lagrange::Phi::Phi(Lagrange* owner,
                     GeoIndex position,
                     int positionIndex,
                     int nodeIndex)
    : vertices(NULL)
  {
    FUNCNAME("Lagrange::Phi::Phi()");

    // get relevant vertices
    Lagrange::setVertices(owner->getDim(),
                          owner->getDegree(),
                          position,
                          positionIndex,
                          nodeIndex,
                          &vertices);

    // set function pointer
    switch (owner->getDegree())
    {
    case 0:
      switch (position)
      {
      case CENTER:
        func = phi0c;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 1:
      switch (position)
      {
      case VERTEX:
        func = phi1v;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 2:
      switch (position)
      {
      case VERTEX:
        func = phi2v;
        break;
      case EDGE:
        TEST_EXIT_DBG(owner->getDim() > 1)("no edge in 1d\n");
        func = phi2e;
        break;
      case CENTER:
        TEST_EXIT_DBG(owner->getDim() == 1)("no center dofs for dim != 1\n");
        func = phi2e;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 3:
      switch (position)
      {
      case VERTEX:
        func = phi3v;
        break;
      case EDGE:
        func = phi3e;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = phi3f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          func = phi3e;
          break;
        case 2:
          func = phi3f;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 4:
      switch (position)
      {
      case VERTEX:
        func = phi4v;
        break;
      case EDGE:
        if (nodeIndex == 1)
          func = phi4e1;
        else
          func = phi4e0;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = phi4f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          if (nodeIndex == 1)
            func = phi4e1;
          else
            func = phi4e0;
          break;
        case 2:
          func = phi4f;
          break;
        case 3:
          func = phi4c;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    default:
      ERROR_EXIT("invalid degree\n");
    }
  }


  Lagrange::Phi::~Phi()
  {
    delete [] vertices;
  }


  Lagrange::GrdPhi::GrdPhi(Lagrange* owner,
                           GeoIndex position,
                           int positionIndex,
                           int nodeIndex)
    : vertices(NULL)
  {
    // get relevant vertices
    Lagrange::setVertices(owner->getDim(),
                          owner->getDegree(),
                          position,
                          positionIndex,
                          nodeIndex,
                          &vertices);

    // set function pointer
    switch (owner->getDegree())
    {
    case 0:
      switch (position)
      {
      case CENTER:
        func = grdPhi0c;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 1:
      switch (position)
      {
      case VERTEX:
        func = grdPhi1v;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 2:
      switch (position)
      {
      case VERTEX:
        func = grdPhi2v;
        break;
      case EDGE:
        func = grdPhi2e;
        break;
      case CENTER:
        TEST_EXIT_DBG(owner->getDim() == 1)("no center dofs for dim != 1\n");
        func = grdPhi2e;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 3:
      switch (position)
      {
      case VERTEX:
        func = grdPhi3v;
        break;
      case EDGE:
        func = grdPhi3e;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = grdPhi3f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          func = grdPhi3e;
          break;
        case 2:
          func = grdPhi3f;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 4:
      switch (position)
      {
      case VERTEX:
        func = grdPhi4v;
        break;
      case EDGE:
        if (nodeIndex == 1)
          func = grdPhi4e1;
        else
          func = grdPhi4e0;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = grdPhi4f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          if (nodeIndex == 1)
            func = grdPhi4e1;
          else
            func = grdPhi4e0;
          break;
        case 2:
          func = grdPhi4f;
          break;
        case 3:
          func = grdPhi4c;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    default:
      ERROR_EXIT("invalid degree\n");
    }
  }

  Lagrange::GrdPhi::~GrdPhi()
  {
    delete [] vertices;
  }

  Lagrange::D2Phi::D2Phi(Lagrange* owner,
                         GeoIndex position,
                         int positionIndex,
                         int nodeIndex)
    : vertices(NULL)
  {
    // get relevant vertices
    Lagrange::setVertices(owner->getDim(),
                          owner->getDegree(),
                          position,
                          positionIndex,
                          nodeIndex,
                          &vertices);

    // set function pointer
    switch (owner->getDegree())
    {
    case 0:
      switch (position)
      {
      case CENTER:
        func = D2Phi0c;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 1:
      switch (position)
      {
      case VERTEX:
        func = D2Phi1v;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 2:
      switch (position)
      {
      case VERTEX:
        func = D2Phi2v;
        break;
      case EDGE:
        TEST_EXIT_DBG(owner->getDim() > 1)("no edge in 1d\n");
        func = D2Phi2e;
        break;
      case CENTER:
        TEST_EXIT_DBG(owner->getDim() == 1)("no center dofs for dim != 1\n");
        func = D2Phi2e;
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 3:
      switch (position)
      {
      case VERTEX:
        func = D2Phi3v;
        break;
      case EDGE:
        func = D2Phi3e;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = D2Phi3f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          func = D2Phi3e;
          break;
        case 2:
          func = D2Phi3f;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    case 4:
      switch (position)
      {
      case VERTEX:
        func = D2Phi4v;
        break;
      case EDGE:
        if (nodeIndex == 1)
          func = D2Phi4e1;
        else
          func = D2Phi4e0;
        break;
      case FACE:
        TEST_EXIT_DBG(owner->getDim() >= 3)("no faces in dim < 3\n");
        func = D2Phi4f;
        break;
      case CENTER:
        switch (owner->getDim())
        {
        case 1:
          if (nodeIndex == 1)
            func = D2Phi4e1;
          else
            func = D2Phi4e0;
          break;
        case 2:
          func = D2Phi4f;
          break;
        case 3:
          func = D2Phi4c;
          break;
        default:
          ERROR_EXIT("invalid dim\n");
          break;
        }
        break;
      default:
        ERROR_EXIT("invalid position\n");
      }
      break;
    default:
      ERROR_EXIT("invalid degree\n");
    }
  }

  Lagrange::D2Phi::~D2Phi()
  {
    delete [] vertices;
  }

  void Lagrange::createCoords(int* coordInd, int numCoords, int dimIndex, int rest,
                              DimVec<double>* vec)
  {
    if (vec == NULL)
      vec = new DimVec<double>(dim, 0.0);

    if (dimIndex == numCoords - 1)
    {
      (*vec)[coordInd[dimIndex]] = double(rest) / degree;
      DimVec<double>* newCoords = new DimVec<double>(*vec);
      bary->push_back(newCoords);
    }
    else
    {
      for (int i = rest - 1; i >= 1; i--)
      {
        (*vec)[coordInd[dimIndex]] = double(i) / degree;
        createCoords(coordInd, numCoords, dimIndex + 1, rest - i, vec);
      }
    }
  }

  void Lagrange::setBary()
  {
    bary = &baryDimDegree[dim][degree];

    if (static_cast<int>(bary->size()) == 0)
    {
      for (int i = 0; i <= dim; i++)   // for all positions
      {
        int partsAtPos = Global::getGeo(INDEX_OF_DIM(i, dim), dim);
        for (int j = 0; j < partsAtPos; j++)   // for all vertices/edges/faces/...
        {
          int* coordInd = new int[i + 1];      // indices of relevant coords
          for (int k = 0; k < i + 1; k++)
            coordInd[k] = Global::getReferenceElement(dim)->
                          getVertexOfPosition(INDEX_OF_DIM(i, dim), j, k);

          createCoords(coordInd, i + 1, 0, degree);
          delete [] coordInd;
          if (static_cast<int>(bary->size()) == nBasFcts)
            return;
        }
      }
    }
  }

  int Lagrange::getNumberOfDofs(int dim, int degree)
  {
    int result = 0;
    for (int i = 0; i <= degree; i++)
      result += fac(dim - 1 + i) / (fac(i) * fac(dim - 1));

    return result;
  }


  int* Lagrange::orderOfPositionIndices(const Element* el,
                                        GeoIndex position,
                                        int positionIndex) const
  {
    static int sortedVertex = 0;
    static int sortedEdgeDeg2 = 0;
    static int sortedEdgeDeg3[2][2] = {{0, 1}, {1, 0}};
    static int sortedEdgeDeg4[2][3] = {{0, 1, 2}, {2, 1, 0}};
    static int sortedFaceDeg3 = 0;
    static int sortedFaceDeg4[7][3] = {{0,0,0}, {0,2,1}, {1,0,2}, {0,1,2},
      {2,1,0}, {1,2,0}, {2,0,1}
    };
    static int sortedCenterDeg4 = 0;

    int dimOfPosition = DIM_OF_INDEX(position, dim);

    // vertex
    if (dimOfPosition == 0)
      return &sortedVertex;

    // edge
    if (dimOfPosition == 1 && degree == 2)
      return &sortedEdgeDeg2;

    int vertex[3];
    DegreeOfFreedom** dof = const_cast<DegreeOfFreedom**>(el->getDof());
    int verticesOfPosition = dimOfPosition + 1;

    for (int i = 0; i < verticesOfPosition; i++)
      vertex[i] = Global::getReferenceElement(dim)->
                  getVertexOfPosition(position, positionIndex, i);

    if (dimOfPosition == 1)
    {
      if (degree == 3)
      {
        if (dof[vertex[0]][0] < dof[vertex[1]][0])
          return sortedEdgeDeg3[0];
        else
          return sortedEdgeDeg3[1];
      }
      else     // degree == 4
      {
        if (dof[vertex[0]][0] < dof[vertex[1]][0])
          return sortedEdgeDeg4[0];
        else
          return sortedEdgeDeg4[1];
      }
    }

    // face
    if (dimOfPosition == 2)
    {
      if (degree == 3)
      {
        return &sortedFaceDeg3;
      }
      else      // degree == 4!
      {
        int no = 0;
        const Element* refElem = Global::getReferenceElement(dim);

        if (dof[refElem->getVertexOfPosition(position, positionIndex, 0)][0] <
            dof[refElem->getVertexOfPosition(position, positionIndex, 1)][0])
          no++;

        if (dof[refElem->getVertexOfPosition(position, positionIndex, 1)][0] <
            dof[refElem->getVertexOfPosition(position, positionIndex, 2)][0])
          no += 2;

        if (dof[refElem->getVertexOfPosition(position, positionIndex, 2)][0] <
            dof[refElem->getVertexOfPosition(position, positionIndex, 0)][0])
          no += 4;

        return sortedFaceDeg4[no];
      }
    }

    // center
    if (dimOfPosition == 3 && degree == 4)
      return &sortedCenterDeg4;

    ERROR_EXIT("should not be reached\n");
    return NULL;
  }

  void Lagrange::getBound(const ElInfo* elInfo, BoundaryType* bound) const
  {
    elInfo->testFlag(Mesh::FILL_BOUND);

    // boundaries
    int index = 0;
    int offset = 0;
    BoundaryType boundaryType;
    for (int i = dim - 1; i > 0; i--)
      offset += Global::getGeo(INDEX_OF_DIM(i, dim), dim);

    for (int i = 0; i < dim; i++)
    {
      int jto = offset + Global::getGeo(INDEX_OF_DIM(i, dim), dim);
      for (int j = offset; j < jto; j++)
      {
        boundaryType = elInfo->getBoundary(j);
        int kto = (*nDOF)[INDEX_OF_DIM(i, dim)];
        for (int k = 0; k < kto; k++)
          bound[index++] = boundaryType;
      }
      offset -= Global::getGeo(INDEX_OF_DIM(i + 1, dim), dim);
    }

    // interior nodes in the center
    for (int i = 0; i < (*nDOF)[CENTER]; i++)
      bound[index++] = INTERIOR;

    TEST_EXIT_DBG(index == nBasFcts)("found not enough boundarys\n");
  }


  void Lagrange::interpol(const ElInfo* elInfo,
                          int no,
                          const int* b_no,
                          std::function<double(WorldVector<double>)> f,
                          DenseVector<double>& rvec) const
  {
    FUNCNAME_DBG("Lagrange::interpol()");

    WorldVector<double> x;

    elInfo->testFlag(Mesh::FILL_COORDS);

    if (b_no)
    {
      TEST_EXIT_DBG(no >= 0 && no < getNumber())("Something is wrong!\n");

      for (int i = 0; i < no; i++)
      {
        if (b_no[i] < Global::getGeo(VERTEX, dim))
        {
          rvec[i] = f(elInfo->getCoord(b_no[i]));
        }
        else
        {
          elInfo->coordToWorld(*(*bary)[b_no[i]], x);
          rvec[i] = f(x);
        }
      }
    }
    else
    {
      int vertices = Global::getGeo(VERTEX, dim);
      for (int i = 0; i < vertices; i++)
        rvec[i] = f(elInfo->getCoord(i));
      for (int i = vertices; i < nBasFcts; i++)
      {
        elInfo->coordToWorld(*(*bary)[i], x);
        rvec[i] = f(x);
      }
    }
  }


  void Lagrange::interpol(const ElInfo* elInfo,
                          int no,
                          const int* b_no,
                          std::function<WorldVector<double>(WorldVector<double>)> f,
                          DenseVector<WorldVector<double>>& rvec) const
  {
    FUNCNAME_DBG("*Lagrange::interpol()");

    WorldVector<double> x;

    elInfo->testFlag(Mesh::FILL_COORDS);

    int vertices = Global::getGeo(VERTEX, dim);

    if (b_no)
    {
      TEST_EXIT_DBG(no >= 0 && no < getNumber())("Something is wrong!\n");

      for (int i = 0; i < no; i++)
      {
        if (b_no[i] < Global::getGeo(VERTEX, dim))
        {
          rvec[i] = f(elInfo->getCoord(b_no[i]));
        }
        else
        {
          elInfo->coordToWorld(*(*bary)[b_no[i]], x);
          rvec[i] = f(x);
        }
      }
    }
    else
    {
      for (int i = 0; i < vertices; i++)
        rvec[i] = f(elInfo->getCoord(i));
      for (int i = vertices; i < nBasFcts; i++)
      {
        elInfo->coordToWorld(*(*bary)[i], x);
        rvec[i] = f(x);
      }
    }
  }


  void Lagrange::getLocalIndices(const Element* el,
                                 const DOFAdmin* admin,
                                 std::vector<DegreeOfFreedom>& dofs) const
  {
    if (static_cast<int>(dofs.size()) < nBasFcts)
      dofs.resize(nBasFcts);

    const DegreeOfFreedom** dof = el->getDof();

    for (int pos = 0, j = 0; pos <= dim; pos++)
    {
      GeoIndex posIndex = INDEX_OF_DIM(pos, dim);
      int nrDOFs = admin->getNumberOfDofs(posIndex);

      if (nrDOFs > 0)
      {
        int n0 = admin->getNumberOfPreDofs(posIndex);
        int node0 = admin->getMesh()->getNode(posIndex);
        int num = Global::getGeo(posIndex, dim);

        for (int i = 0; i < num; node0++, i++)
        {
          const int* indi = orderOfPositionIndices(el, posIndex, i);

          for (int k = 0; k < nrDOFs; k++)
            dofs[j++] = dof[node0][n0 + indi[k]];
        }
      }
    }
  }


  void Lagrange::getLocalDofPtrVec(const Element* el,
                                   const DOFAdmin* admin,
                                   std::vector<const DegreeOfFreedom*>& dofs) const
  {
    if (static_cast<int>(dofs.size()) < nBasFcts)
      dofs.resize(nBasFcts);

    const DegreeOfFreedom** dof = el->getDof();

    for (int pos = 0, j = 0; pos <= dim; pos++)
    {
      GeoIndex posIndex = INDEX_OF_DIM(pos, dim);
      int nrDOFs = admin->getNumberOfDofs(posIndex);

      if (nrDOFs)
      {
        int n0 = admin->getNumberOfPreDofs(posIndex);
        int node0 = admin->getMesh()->getNode(posIndex);
        int num = Global::getGeo(posIndex, dim);

        for (int i = 0; i < num; node0++, i++)
        {
          const int* indi = orderOfPositionIndices(el, posIndex, i);

          for (int k = 0; k < nrDOFs; k++)
            dofs[j++] = &(dof[node0][n0 + indi[k]]);
        }
      }
    }
  }


  void Lagrange::l2ScpFctBas(Quadrature* /*q*/,
                             std::function<WorldVector<double>(WorldVector<double>)> /*f*/,
                             DOFVector<WorldVector<double>>* /*fh*/)
  {
    ERROR_EXIT("not yet\n");
  }


  void Lagrange::l2ScpFctBas(Quadrature* quad,
                             std::function<double(WorldVector<double>)> f,
                             DOFVector<double>* fh)
  {
    FUNCNAME_DBG("Lagrange::l2ScpFctBas()");

    TEST_EXIT_DBG(fh)("no DOF_REAL_VEC fh\n");
    TEST_EXIT_DBG(fh->getFeSpace())
    ("no fe_space in DOF_REAL_VEC %s\n", fh->getName().c_str());
    TEST_EXIT_DBG(fh->getFeSpace()->getBasisFcts() == this)
    ("wrong basis fcts for fh\n");
    TEST_EXIT_DBG(!fh->getFeSpace()->getMesh()->getParametric())
    ("Not yet implemented!");

    if (!quad)
      quad = Quadrature::provideQuadrature(dim, 2 * degree - 2);

    const FastQuadrature* quad_fast =
      FastQuadrature::provideFastQuadrature(this, *quad, INIT_PHI);
    vector<double> wdetf_qp(quad->getNumPoints());
    int nPoints = quad->getNumPoints();
    DOFAdmin* admin = fh->getFeSpace()->getAdmin();
    WorldVector<double> x;
    vector<DegreeOfFreedom> dof;

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(fh->getFeSpace()->getMesh(), -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo)
    {
      getLocalIndices(elInfo->getElement(), admin, dof);
      double det = elInfo->getDet();

      for (int iq = 0; iq < nPoints; iq++)
      {
        elInfo->coordToWorld(quad->getLambda(iq), x);
        wdetf_qp[iq] = quad->getWeight(iq) * det * f(x);
      }

      for (int j = 0; j < nBasFcts; j++)
      {
        double val = 0.0;
        for (int iq = 0; iq < nPoints; iq++)
          val += quad_fast->getPhi(iq, j) * wdetf_qp[iq];

        (*fh)[dof[j]] += val;
      }

      elInfo = stack.traverseNext(elInfo);
    }
  }


  // ===== refineInter functions =====

  void  Lagrange::refineInter0(DOFIndexed<double>* drv,
                               RCNeighbourList* list,
                               int n,
                               BasisFunction* /*basFct*/)
  {
    if (n < 1)
      return;

    int n0 = drv->getFeSpace()->getAdmin()->getNumberOfPreDofs(CENTER);
    Element* el = list->getElement(0);
    int node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    // Parent center
    DegreeOfFreedom dof0 = el->getDof(node, n0);
    // Newest vertex is center
    DegreeOfFreedom dof_new = el->getChild(0)->getDof(node, n0);

    (*drv)[dof_new] = (*drv)[dof0];
    // Newest vertex is center
    dof_new = el->getChild(1)->getDof(node, n0);
    (*drv)[dof_new] = (*drv)[dof0];
  }

  void Lagrange::refineInter1(DOFIndexed<double>* drv,
                              RCNeighbourList* list,
                              int n,
                              BasisFunction* /*basFct*/)
  {
    if (n < 1)
      return;

    int dim = drv->getFeSpace()->getMesh()->getDim();
    int n0 = drv->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX);
    Element* el = list->getElement(0);
    // 1st endpoint of refinement edge
    DegreeOfFreedom dof0 = el->getDof(0, n0);
    // 2nd endpoint of refinement edge
    DegreeOfFreedom dof1 = el->getDof(1, n0);
    // newest vertex is DIM
    DegreeOfFreedom dof_new = el->getChild(0)->getDof(dim, n0);
    (*drv)[dof_new] = 0.5 * ((*drv)[dof0] + (*drv)[dof1]);
  }

  void Lagrange::refineInter2_1d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n,
                                 BasisFunction* basFct)
  {
    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);
    int n0 = admin->getNumberOfPreDofs(VERTEX);

    /****************************************************************************/
    /*  newest vertex of child[0] and child[1]                                  */
    /****************************************************************************/

    // newest vertex is DIM
    DegreeOfFreedom cdof = el->getChild(0)->getDof(node + 1, n0);
    (*drv)[cdof] = (*drv)[pdof[2]];

    node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    n0 = admin->getNumberOfPreDofs(CENTER);

    /****************************************************************************/
    /*  midpoint of edge on child[0] at the refinement edge                     */
    /****************************************************************************/

    cdof = el->getChild(0)->getDof(node, n0);
    (*drv)[cdof] =
      0.375 * (*drv)[pdof[0]] - 0.125 * (*drv)[pdof[1]] + 0.75 * (*drv)[pdof[2]];

    /****************************************************************************/
    /*  midpoint of edge on child[1] at the refinement edge                     */
    /****************************************************************************/

    cdof = el->getChild(1)->getDof(node, n0);
    (*drv)[cdof] =
      -0.125 * (*drv)[pdof[0]] + 0.375 * (*drv)[pdof[1]] + 0.75 * (*drv)[pdof[2]];
  }

  void Lagrange::refineInter2_2d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    Element* el = list->getElement(0);
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);
    int n0 = admin->getNumberOfPreDofs(VERTEX);

    /****************************************************************************/
    /*  newest vertex of child[0] and child[1]                                  */
    /****************************************************************************/

    // newest vertex is DIM
    DegreeOfFreedom cdof = el->getChild(0)->getDof(node + 2, n0);
    (*drv)[cdof] = (*drv)[pdof[5]];

    node = drv->getFeSpace()->getMesh()->getNode(EDGE);
    n0 = admin->getNumberOfPreDofs(EDGE);

    /****************************************************************************/
    /*  midpoint of edge on child[0] at the refinement edge                     */
    /****************************************************************************/

    cdof = el->getChild(0)->getDof(node, n0);
    (*drv)[cdof] =
      0.375 * (*drv)[pdof[0]] - 0.125 * (*drv)[pdof[1]] + 0.75 * (*drv)[pdof[5]];

    /****************************************************************************/
    /* node in the common edge of child[0] and child[1]                         */
    /****************************************************************************/

    cdof = el->getChild(0)->getDof(node + 1, n0);
    (*drv)[cdof] =
      -0.125 * ((*drv)[pdof[0]] + (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[5]]
      + 0.5 * ((*drv)[pdof[3]] + (*drv)[pdof[4]]);

    /****************************************************************************/
    /*  midpoint of edge on child[1] at the refinement edge                     */
    /****************************************************************************/

    cdof = el->getChild(1)->getDof(node + 1, n0);
    (*drv)[cdof] =
      -0.125 * (*drv)[pdof[0]] + 0.375 * (*drv)[pdof[1]]  + 0.75 * (*drv)[pdof[5]];

    if (n > 1)
    {
      /****************************************************************************/
      /*  adjust the value at the midpoint of the common edge of neigh's children */
      /****************************************************************************/
      el = list->getElement(1);
      basFct->getLocalIndices(el, admin, pdof);

      cdof = el->getChild(0)->getDof(node + 1, n0);
      (*drv)[cdof] =
        -0.125 * ((*drv)[pdof[0]] + (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[5]]
        + 0.5 * ((*drv)[pdof[3]] + (*drv)[pdof[4]]);
    }
  }

  void Lagrange::refineInter2_3d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME_DBG("Lagrange::refineInter2_3d()");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(10);
    basFct->getLocalIndices(el, admin, pdof);

    int node0 = drv->getFeSpace()->getMesh()->getNode(EDGE);
    int n0 = admin->getNumberOfPreDofs(EDGE);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    vector<DegreeOfFreedom> cdof(10);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    (*drv)[cdof[3]] = ((*drv)[pdof[4]]);
    (*drv)[cdof[6]] =
      (0.375 * (*drv)[pdof[0]] - 0.125 * (*drv)[pdof[1]]
       + 0.75 * (*drv)[pdof[4]]);
    (*drv)[cdof[8]] =
      (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
       + 0.5 * ((*drv)[pdof[5]] + (*drv)[pdof[7]]));
    (*drv)[cdof[9]] =
      (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
       + 0.5 * ((*drv)[pdof[6]] + (*drv)[pdof[8]]));

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    DegreeOfFreedom cdofi = el->getChild(1)->getDof(node0 + 2, n0);
    (*drv)[cdofi] =
      (-0.125 * (*drv)[pdof[0]] + 0.375 * (*drv)[pdof[1]]
       + 0.75 * (*drv)[pdof[4]]);

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);
      TEST_EXIT_DBG(el)("Should not happen!\n");

      basFct->getLocalIndices(el, admin, pdof);

      int lr_set = 0;
      if (list->getNeighbourElement(i, 0) && list->getNeighbourNr(i, 0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i, 1) && list->getNeighbourNr(i, 1) < i)
        lr_set += 2;

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      switch (lr_set)
      {
      case 0:
        cdofi = el->getChild(0)->getDof(node0 + 4, n0);
        (*drv)[cdofi] =
          (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
           + 0.5 * ((*drv)[pdof[5]] + (*drv)[pdof[7]]));
        cdofi = el->getChild(0)->getDof(node0 + 5, n0);
        (*drv)[cdofi] =
          (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
           + 0.5 * ((*drv)[pdof[6]] + (*drv)[pdof[8]]));
        break;
      case 1:
        cdofi = el->getChild(0)->getDof(node0 + 4, n0);
        (*drv)[cdofi] =
          (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
           + 0.5 * ((*drv)[pdof[5]] + (*drv)[pdof[7]]));
        break;
      case 2:
        cdofi = el->getChild(0)->getDof(node0 + 5, n0);
        (*drv)[cdofi] =
          (0.125 * (-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.25 * (*drv)[pdof[4]]
           + 0.5 * ((*drv)[pdof[6]] + (*drv)[pdof[8]]));
      }
    }
  }

  void Lagrange::refineInter3_1d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof, cdof;
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cdof[2]] =
      (-0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.5625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[3]] =
      (0.3125 * ((*drv)[pdof[0]] - (*drv)[pdof[8]]) + 0.0625 * (*drv)[pdof[1]]
       + 0.9375 * (*drv)[pdof[7]]);
    (*drv)[cdof[4]] = (*drv)[pdof[7]];
    (*drv)[cdof[5]] = (*drv)[pdof[9]];
    (*drv)[cdof[6]] =
      (0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       - 0.25 * ((*drv)[pdof[3]] + (*drv)[pdof[6]])
       + 0.5 * ((*drv)[pdof[4]] + (*drv)[pdof[5]] + (*drv)[pdof[9]])
       - 0.0625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[9]] =
      (0.0625 * (-(*drv)[pdof[0]] + (*drv)[pdof[1]]) - 0.125 * (*drv)[pdof[3]]
       + 0.375 * (*drv)[pdof[6]] + 0.1875 * ((*drv)[pdof[7]] - (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[5]] =  (*drv)[pdof[8]];
    (*drv)[cdof[6]] =
      (0.0625 * (*drv)[pdof[0]] + 0.9375 * (*drv)[pdof[8]]
       + 0.3125 * ((*drv)[pdof[1]] - (*drv)[pdof[7]]));
    (*drv)[cdof[9]] =
      (0.0625 * ((*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.375 * (*drv)[pdof[3]]
       - 0.125 * (*drv)[pdof[6]] + 0.1875 * (-(*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);

    if (n <= 1)
      return;

    /****************************************************************************/
    /*  adjust the value on the neihgbour                                       */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);

    /****************************************************************************/
    /*  values on neigh's child[0]                                              */
    /****************************************************************************/
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    (*drv)[cdof[5]] =  (*drv)[pdof[9]];
    (*drv)[cdof[6]] =
      (0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       - 0.25 * ((*drv)[pdof[3]] + (*drv)[pdof[6]])
       + 0.5 * ((*drv)[pdof[4]] + (*drv)[pdof[5]] + (*drv)[pdof[9]])
       - 0.0625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[9]] =
      (0.0625 * (-(*drv)[pdof[0]] + (*drv)[pdof[1]]) - 0.12500 * (*drv)[pdof[3]]
       + 0.375 * (*drv)[pdof[6]] + 0.1875 * ((*drv)[pdof[7]] - (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);
    /****************************************************************************/
    /*  (*drv)alues on neigh's child[0]                                         */
    /****************************************************************************/

    int node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    int n0 = admin->getNumberOfPreDofs(CENTER);
    DegreeOfFreedom dof9 = el->getChild(1)->getDof(node, n0);

    (*drv)[dof9] =
      (0.0625 * ((*drv)[pdof[0]] - (*drv)[pdof[1]]) +  0.375 * (*drv)[pdof[3]]
       - 0.125 * (*drv)[pdof[6]] + 0.1875 * (-(*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.75  *(*drv)[pdof[9]]);
  }

  void Lagrange::refineInter3_2d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n,
                                 BasisFunction* basFct)
  {
    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();

    vector<DegreeOfFreedom> pdof, cdof;
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cdof[2]] =
      (-0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.5625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[3]] =
      (0.3125 * ((*drv)[pdof[0]] - (*drv)[pdof[8]]) + 0.0625 * (*drv)[pdof[1]]
       + 0.9375 * (*drv)[pdof[7]]);
    (*drv)[cdof[4]] = (*drv)[pdof[7]];
    (*drv)[cdof[5]] = (*drv)[pdof[9]];
    (*drv)[cdof[6]] =
      (0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       - 0.25 * ((*drv)[pdof[3]] + (*drv)[pdof[6]])
       + 0.5 * ((*drv)[pdof[4]] + (*drv)[pdof[5]] + (*drv)[pdof[9]])
       - 0.0625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[9]] =
      (0.0625 * (-(*drv)[pdof[0]] + (*drv)[pdof[1]]) - 0.125 * (*drv)[pdof[3]]
       + 0.375 * (*drv)[pdof[6]] + 0.1875 * ((*drv)[pdof[7]] - (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[5]] = (*drv)[pdof[8]];
    (*drv)[cdof[6]] =
      (0.0625 * (*drv)[pdof[0]] + 0.9375 * (*drv)[pdof[8]]
       + 0.3125 * ((*drv)[pdof[1]] - (*drv)[pdof[7]]));
    (*drv)[cdof[9]] =
      (0.0625 * ((*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.375 * (*drv)[pdof[3]]
       - 0.125 * (*drv)[pdof[6]] + 0.1875 * (-(*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);

    if (n <= 1)
      return;

    /****************************************************************************/
    /*  adjust the value on the neihgbour                                       */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neigh's child[0]                                              */
    /****************************************************************************/

    (*drv)[cdof[5]] =  (*drv)[pdof[9]];
    (*drv)[cdof[6]] =
      (0.0625 * ((*drv)[pdof[0]] + (*drv)[pdof[1]])
       - 0.25 * ((*drv)[pdof[3]] + (*drv)[pdof[6]])
       + 0.5 * ((*drv)[pdof[4]] + (*drv)[pdof[5]] + (*drv)[pdof[9]])
       - 0.0625 * ((*drv)[pdof[7]] + (*drv)[pdof[8]]));
    (*drv)[cdof[9]] =
      (0.0625 * (-(*drv)[pdof[0]] + (*drv)[pdof[1]]) - 0.12500 * (*drv)[pdof[3]]
       + 0.375 * (*drv)[pdof[6]] + 0.1875 * ((*drv)[pdof[7]] - (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);
    /****************************************************************************/
    /*  values on neigh's child[0]                                              */
    /****************************************************************************/

    int node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    int n0 = admin->getNumberOfPreDofs(CENTER);
    DegreeOfFreedom dof9 = el->getChild(1)->getDof(node, n0);

    (*drv)[dof9] =
      (0.0625 * ((*drv)[pdof[0]] - (*drv)[pdof[1]]) +  0.375 * (*drv)[pdof[3]]
       - 0.125 * (*drv)[pdof[6]] + 0.1875 * (-(*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.75 * (*drv)[pdof[9]]);
  }

  void Lagrange::refineInter3_3d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME_DBG("Lagrange::refineInter3_3d()");

    if (n < 1)
      return;

    DegreeOfFreedom cdi;
    int i, lr_set, node0, n0;
    Element* el = list->getElement(0);
    int typ = list->getType(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pd(20), cd(20);
    basFct->getLocalIndices(el, admin, pd);
    basFct->getLocalIndices(el->getChild(0), admin, cd);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cd[3]] =
      (0.0625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.5625*((*drv)[pd[4]] + (*drv)[pd[5]]));
    (*drv)[cd[8]] =
      (0.3125*((*drv)[pd[0]] - (*drv)[pd[5]]) + 0.0625*(*drv)[pd[1]]
       + 0.937500*(*drv)[pd[4]]);
    (*drv)[cd[9]] = (*drv)[pd[4]];
    (*drv)[cd[12]] =
      (0.0625*((*drv)[pd[0]] + (*drv)[pd[1]] - (*drv)[pd[4]] - (*drv)[pd[5]])
       + 0.25*(-(*drv)[pd[6]] - (*drv)[pd[10]])
       + 0.5*((*drv)[pd[7]] + (*drv)[pd[11]] + (*drv)[pd[19]]));
    (*drv)[cd[13]] = (*drv)[pd[19]];
    (*drv)[cd[14]] =
      (0.0625*((*drv)[pd[0]] + (*drv)[pd[1]] - (*drv)[pd[4]] - (*drv)[pd[5]])
       + 0.25*(-(*drv)[pd[8]] - (*drv)[pd[12]])
       + 0.5*((*drv)[pd[9]] + (*drv)[pd[13]] + (*drv)[pd[18]]));
    (*drv)[cd[15]] = ((*drv)[pd[18]]);
    (*drv)[cd[16]] =
      (0.0625*((*drv)[pd[0]] + (*drv)[pd[1]] - (*drv)[pd[4]] - (*drv)[pd[5]])
       + 0.125*(-(*drv)[pd[6]]-(*drv)[pd[8]]-(*drv)[pd[10]]-(*drv)[pd[12]])
       + 0.5*((*drv)[pd[16]] + (*drv)[pd[17]])
       + 0.25*((*drv)[pd[18]] + (*drv)[pd[19]]));
    (*drv)[cd[17]] =
      (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
       + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[8]]
       - 0.125*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
    (*drv)[cd[18]] =
      (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
       + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[6]]
       - 0.125*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cd);

    if (typ == 0)
    {
      /****************************************************************************/
      /*  parent of el_type 0                                                     */
      /****************************************************************************/

      (*drv)[cd[8]] =
        (0.0625*(*drv)[pd[0]] + 0.3125*((*drv)[pd[1]] - (*drv)[pd[4]])
         + 0.9375*(*drv)[pd[5]]);
      (*drv)[cd[9]] = (*drv)[pd[5]];
      (*drv)[cd[17]] =
        (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]])
         - 0.125*(*drv)[pd[6]] + 0.375*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
      (*drv)[cd[18]] =
        (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[8]]
         + 0.375*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
    }
    else
    {
      /****************************************************************************/
      /*  parent of el_type 0                                                     */
      /****************************************************************************/

      (*drv)[cd[8]] =
        (0.0625*(*drv)[pd[0]] + 0.3125*((*drv)[pd[1]] - (*drv)[pd[4]])
         + 0.9375*(*drv)[pd[5]]);
      (*drv)[cd[9]] = (*drv)[pd[5]];
      (*drv)[cd[17]] =
        (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[8]]
         + 0.375*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
      (*drv)[cd[18]] =
        (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]]) +
         0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[6]]
         + 0.375*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
    }

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    node0 = drv->getFeSpace()->getMesh()->getNode(FACE);
    n0 = admin->getNumberOfPreDofs(FACE);

    for (i = 1; i < n; i++)
    {
      el = list->getElement(i);
      typ = list->getType(i);
      basFct->getLocalIndices(el, admin, pd);

      lr_set = 0;
      if (list->getNeighbourElement(i, 0)  &&  list->getNeighbourNr(i, 0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i, 1)  &&  list->getNeighbourNr(i, 1) < i)
        lr_set += 2;

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(0), admin, cd);

      switch (lr_set)
      {
      case 0:
        (*drv)[cd[12]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.25*(-(*drv)[pd[6]] - (*drv)[pd[10]])
           + 0.5*((*drv)[pd[7]] + (*drv)[pd[11]] + (*drv)[pd[19]]));
        (*drv)[cd[13]] = (*drv)[pd[19]];
        (*drv)[cd[14]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.25*(-(*drv)[pd[8]] - (*drv)[pd[12]])
           + 0.5*((*drv)[pd[9]] + (*drv)[pd[13]] + (*drv)[pd[18]]));
        (*drv)[cd[15]] = (*drv)[pd[18]];
        (*drv)[cd[16]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.125*(-(*drv)[pd[6]]-(*drv)[pd[8]]-(*drv)[pd[10]]-(*drv)[pd[12]])
           + 0.5*((*drv)[pd[16]] + (*drv)[pd[17]])
           + 0.25*((*drv)[pd[18]] + (*drv)[pd[19]]));
        (*drv)[cd[17]] =
          (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[8]]
           - 0.125*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
        (*drv)[cd[18]] =
          (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[6]]
           - 0.125*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
        break;
      case 1:
        (*drv)[cd[12]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.25*(-(*drv)[pd[6]] - (*drv)[pd[10]])
           + 0.5*((*drv)[pd[7]] + (*drv)[pd[11]] + (*drv)[pd[19]]));
        (*drv)[cd[13]] = (*drv)[pd[19]];
        (*drv)[cd[16]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.125*(-(*drv)[pd[6]]-(*drv)[pd[8]]-(*drv)[pd[10]]-(*drv)[pd[12]])
           + 0.5*((*drv)[pd[16]] + (*drv)[pd[17]])
           + 0.25*((*drv)[pd[18]] + (*drv)[pd[19]]));
        (*drv)[cd[18]] =
          (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[6]]
           - 0.125*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
        break;
      case 2:
        (*drv)[cd[14]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.25*(-(*drv)[pd[8]] - (*drv)[pd[12]])
           + 0.5*((*drv)[pd[9]] + (*drv)[pd[13]] + (*drv)[pd[18]]));
        (*drv)[cd[15]] = (*drv)[pd[18]];
        (*drv)[cd[16]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.125*(-(*drv)[pd[6]]-(*drv)[pd[8]]-(*drv)[pd[10]]-(*drv)[pd[12]])
           + 0.5*((*drv)[pd[16]] + (*drv)[pd[17]])
           + 0.25*((*drv)[pd[18]] + (*drv)[pd[19]]));
        (*drv)[cd[17]] =
          (0.0625*(-(*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.1875*((*drv)[pd[4]] - (*drv)[pd[5]]) + 0.375*(*drv)[pd[8]]
           - 0.125*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
        break;
      case 3:
        (*drv)[cd[16]] =
          (0.0625*((*drv)[pd[0]]+(*drv)[pd[1]]-(*drv)[pd[4]]-(*drv)[pd[5]])
           + 0.125*(-(*drv)[pd[6]]-(*drv)[pd[8]]-(*drv)[pd[10]]-(*drv)[pd[12]])
           + 0.5*((*drv)[pd[16]] + (*drv)[pd[17]])
           + 0.25*((*drv)[pd[18]] + (*drv)[pd[19]]));
        break;
      }

      /****************************************************************************/
      /*  values on child[1]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(1), admin, cd);

      if (typ == 0)
      {
        switch (lr_set)
        {
        case 1:
          cdi = el->getChild(1)->getDof(node0 + 1, n0);
          TEST_EXIT_DBG(cdi == cd[17])("cdi != cd[17]\n");
          (*drv)[cdi] =
            (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[6]]
             + 0.375*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
          break;
        case 2:
          cdi = el->getChild(1)->getDof(node0 + 2, n0);
          TEST_EXIT_DBG(cdi == cd[18])("cdi != cd[18]\n");
          (*drv)[cdi] =
            (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[8]]
             + 0.375*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
          break;
        }
      }
      else
      {
        switch (lr_set)
        {
        case 1:
          cdi = el->getChild(1)->getDof(node0 + 2, n0);
          TEST_EXIT_DBG(cdi == cd[18])("cdi != cd[18]\n");
          (*drv)[cdi] =
            (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[6]]
             + 0.375*(*drv)[pd[10]] + 0.75*(*drv)[pd[19]]);
          break;
        case 2:
          cdi = el->getChild(1)->getDof(node0 + 1, n0);
          TEST_EXIT_DBG(cdi == cd[17])("cdi != cd[17]\n");
          (*drv)[cdi] =
            (0.0625*((*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.1875*(-(*drv)[pd[4]] + (*drv)[pd[5]]) - 0.125*(*drv)[pd[8]]
             + 0.375*(*drv)[pd[12]] + 0.75*(*drv)[pd[18]]);
          break;
        }
      }
    }
  }

  void Lagrange::refineInter4_1d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    Element* el = list->getElement(0);

    vector<DegreeOfFreedom> pdof, cdof;
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cdof[2]] = (*drv)[pdof[10]];
    (*drv)[cdof[3]] =
      (0.2734375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 1.09375*(*drv)[pdof[9]] - 0.546875*(*drv)[pdof[10]]
       + 0.21875*(*drv)[pdof[11]]);
    (*drv)[cdof[4]] = (*drv)[pdof[9]];
    (*drv)[cdof[5]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.46875*(*drv)[pdof[9]] + 0.703125*(*drv)[pdof[10]]
       - 0.15625*(*drv)[pdof[11]]);
    (*drv)[cdof[6]] =
      (0.0234375*((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.0625*(-(*drv)[pdof[3]] - (*drv)[pdof[8]])
       + 0.09375*(-(*drv)[pdof[9]] - (*drv)[pdof[11]]) + 0.140625*(*drv)[pdof[10]]
       + 0.5625*((*drv)[pdof[12]] + (*drv)[pdof[13]]));
    (*drv)[cdof[7]] = (*drv)[pdof[14]];
    (*drv)[cdof[8]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.1875*((*drv)[pdof[3]] + (*drv)[pdof[8]]-(*drv)[pdof[12]]-(*drv)[pdof[13]])
       + 0.375*(-(*drv)[pdof[4]] - (*drv)[pdof[7]])
       + 0.5*((*drv)[pdof[5]] + (*drv)[pdof[6]])
       + 0.03125*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       + 0.015625*(*drv)[pdof[10]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[12]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 0.125*((*drv)[pdof[3]] - (*drv)[pdof[4]] - (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[7]] + (*drv)[pdof[12]] - (*drv)[pdof[13]])
       - 0.03125*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       + 0.09375*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[13]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.0625*(*drv)[pdof[3]]
       + 0.3125*((*drv)[pdof[8]] - (*drv)[pdof[13]])
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[12]]);
    (*drv)[cdof[14]] = (*drv)[pdof[12]];

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[6]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       - 0.15625*(*drv)[pdof[9]] + 0.703125*(*drv)[pdof[10]]
       + 0.46875*(*drv)[pdof[11]]);
    (*drv)[cdof[7]] = (*drv)[pdof[11]];
    (*drv)[cdof[8]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.2734375*(*drv)[pdof[1]]
       + 0.21875*(*drv)[pdof[9]] - 0.546875*(*drv)[pdof[10]]
       + 1.09375*(*drv)[pdof[11]]);
    (*drv)[cdof[12]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.3125*((*drv)[pdof[3]] - (*drv)[pdof[12]]) + 0.0625*(*drv)[pdof[8]]
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[13]]);
    (*drv)[cdof[13]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.125*(-(*drv)[pdof[3]] - (*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[4]] - (*drv)[pdof[12]] + (*drv)[pdof[13]])
       + 0.09375*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       - 0.03125*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[14]] = (*drv)[pdof[13]];

    if (n <= 1)
      return;

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neighbour's child[0]                                          */
    /****************************************************************************/

    (*drv)[cdof[6]] =
      (0.0234375*((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.0625*(-(*drv)[pdof[3]] - (*drv)[pdof[8]])
       + 0.09375*(-(*drv)[pdof[9]] - (*drv)[pdof[11]])
       + 0.140625*(*drv)[pdof[10]] + 0.5625*((*drv)[pdof[12]] + (*drv)[pdof[13]]));
    (*drv)[cdof[7]] = (*drv)[pdof[14]];
    (*drv)[cdof[8]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.1875*((*drv)[pdof[3]]+(*drv)[pdof[8]]-(*drv)[pdof[12]]-(*drv)[pdof[13]])
       + 0.375*(-(*drv)[pdof[4]] - (*drv)[pdof[7]])
       + 0.5*((*drv)[pdof[5]] + (*drv)[pdof[6]])
       + 0.03125*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       + 0.015625*(*drv)[pdof[10]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[12]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 0.125*((*drv)[pdof[3]] - (*drv)[pdof[4]] - (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[7]] + (*drv)[pdof[12]] - (*drv)[pdof[13]])
       - 0.03125*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       + 0.09375*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[13]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.0625*(*drv)[pdof[3]]
       + 0.3125*((*drv)[pdof[8]] - (*drv)[pdof[13]])
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[12]]);
    (*drv)[cdof[14]] = (*drv)[pdof[12]];

    /****************************************************************************/
    /*  values on neighbour's child[1]                                          */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[12]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.3125*((*drv)[pdof[3]] - (*drv)[pdof[12]])
       + 0.0625*(*drv)[pdof[8]] + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[13]]);
    (*drv)[cdof[13]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.125*(-(*drv)[pdof[3]] - (*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[4]] - (*drv)[pdof[12]] + (*drv)[pdof[13]])
       + 0.09375*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       - 0.03125*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[14]] = (*drv)[pdof[13]];
  }

  void Lagrange::refineInter4_2d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(15), cdof(15);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cdof[2]] = (*drv)[pdof[10]];
    (*drv)[cdof[3]] =
      (0.2734375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 1.09375*(*drv)[pdof[9]] - 0.546875*(*drv)[pdof[10]]
       + 0.21875*(*drv)[pdof[11]]);
    (*drv)[cdof[4]] = (*drv)[pdof[9]];
    (*drv)[cdof[5]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.46875*(*drv)[pdof[9]] + 0.703125*(*drv)[pdof[10]]
       - 0.15625*(*drv)[pdof[11]]);
    (*drv)[cdof[6]] =
      (0.0234375*((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.0625*(-(*drv)[pdof[3]] - (*drv)[pdof[8]])
       + 0.09375*(-(*drv)[pdof[9]] - (*drv)[pdof[11]]) + 0.140625*(*drv)[pdof[10]]
       + 0.5625*((*drv)[pdof[12]] + (*drv)[pdof[13]]));
    (*drv)[cdof[7]] = (*drv)[pdof[14]];
    (*drv)[cdof[8]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.1875*((*drv)[pdof[3]] +
                 (*drv)[pdof[8]]-(*drv)[pdof[12]]
                 -(*drv)[pdof[13]])
       + 0.375*(-(*drv)[pdof[4]] - (*drv)[pdof[7]])
       + 0.5*((*drv)[pdof[5]] + (*drv)[pdof[6]])
       + 0.03125*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       + 0.015625*(*drv)[pdof[10]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[12]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 0.125*((*drv)[pdof[3]] - (*drv)[pdof[4]] - (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[7]] + (*drv)[pdof[12]] - (*drv)[pdof[13]])
       - 0.03125*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       + 0.09375*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[13]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.0625*(*drv)[pdof[3]]
       + 0.3125*((*drv)[pdof[8]] - (*drv)[pdof[13]])
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[12]]);
    (*drv)[cdof[14]] = (*drv)[pdof[12]];

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[6]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       - 0.15625*(*drv)[pdof[9]] + 0.703125*(*drv)[pdof[10]]
       + 0.46875*(*drv)[pdof[11]]);
    (*drv)[cdof[7]] = (*drv)[pdof[11]];
    (*drv)[cdof[8]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.2734375*(*drv)[pdof[1]]
       + 0.21875*(*drv)[pdof[9]] - 0.546875*(*drv)[pdof[10]]
       + 1.09375*(*drv)[pdof[11]]);
    (*drv)[cdof[12]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.3125*((*drv)[pdof[3]] - (*drv)[pdof[12]]) + 0.0625*(*drv)[pdof[8]]
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[13]]);
    (*drv)[cdof[13]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.125*(-(*drv)[pdof[3]] - (*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[4]] - (*drv)[pdof[12]] + (*drv)[pdof[13]])
       + 0.09375*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       - 0.03125*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[14]] = (*drv)[pdof[13]];

    if (n <= 1)
      return;

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neighbour's child[0]                                          */
    /****************************************************************************/

    (*drv)[cdof[6]] =
      (0.0234375*((*drv)[pdof[0]] + (*drv)[pdof[1]])
       + 0.0625*(-(*drv)[pdof[3]] - (*drv)[pdof[8]])
       + 0.09375*(-(*drv)[pdof[9]] - (*drv)[pdof[11]])
       + 0.140625*(*drv)[pdof[10]] + 0.5625*((*drv)[pdof[12]] + (*drv)[pdof[13]]));
    (*drv)[cdof[7]] = (*drv)[pdof[14]];
    (*drv)[cdof[8]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.1875*((*drv)[pdof[3]]+(*drv)[pdof[8]]-(*drv)[pdof[12]]-(*drv)[pdof[13]])
       + 0.375*(-(*drv)[pdof[4]] - (*drv)[pdof[7]])
       + 0.5*((*drv)[pdof[5]] + (*drv)[pdof[6]])
       + 0.03125*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       + 0.015625*(*drv)[pdof[10]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[12]] =
      (0.0234375*(*drv)[pdof[0]] - 0.0390625*(*drv)[pdof[1]]
       + 0.125*((*drv)[pdof[3]] - (*drv)[pdof[4]] - (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[7]] + (*drv)[pdof[12]] - (*drv)[pdof[13]])
       - 0.03125*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       + 0.09375*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[13]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]]) + 0.0625*(*drv)[pdof[3]]
       + 0.3125*((*drv)[pdof[8]] - (*drv)[pdof[13]])
       + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[12]]);
    (*drv)[cdof[14]] = (*drv)[pdof[12]];

    /****************************************************************************/
    /*  values on neighbour's child[1]                                          */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[cdof[12]] =
      (0.0390625*(-(*drv)[pdof[0]] - (*drv)[pdof[1]])
       + 0.3125*((*drv)[pdof[3]] - (*drv)[pdof[12]])
       + 0.0625*(*drv)[pdof[8]] + 0.15625*((*drv)[pdof[9]] + (*drv)[pdof[11]])
       - 0.234375*(*drv)[pdof[10]] + 0.9375*(*drv)[pdof[13]]);
    (*drv)[cdof[13]] =
      (-0.0390625*(*drv)[pdof[0]] + 0.0234375*(*drv)[pdof[1]]
       + 0.125*(-(*drv)[pdof[3]] - (*drv)[pdof[7]] + (*drv)[pdof[8]])
       + 0.375*((*drv)[pdof[4]] - (*drv)[pdof[12]] + (*drv)[pdof[13]])
       + 0.09375*(*drv)[pdof[9]] - 0.046875*(*drv)[pdof[10]]
       - 0.03125*(*drv)[pdof[11]] + 0.75*(*drv)[pdof[14]]);
    (*drv)[cdof[14]] = (*drv)[pdof[13]];
  }

  void Lagrange::refineInter4_3d(DOFIndexed<double>* drv,
                                 RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    int typ = list->getType(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pd(35), cd(35);
    basFct->getLocalIndices(el, admin, pd);
    basFct->getLocalIndices(el->getChild(0), admin, cd);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[cd[3]] = ((*drv)[pd[5]]);
    (*drv)[cd[10]] =
      (0.2734375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
       + 1.09375*(*drv)[pd[4]] - 0.546875*(*drv)[pd[5]]
       + 0.21875*(*drv)[pd[6]]);
    (*drv)[cd[11]] = ((*drv)[pd[4]]);
    (*drv)[cd[12]] =
      (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
       + 0.46875*(*drv)[pd[4]] + 0.703125*(*drv)[pd[5]]
       - 0.15625*(*drv)[pd[6]]);
    (*drv)[cd[16]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
       + 0.015625*(*drv)[pd[5]] +
       0.1875*((*drv)[pd[7]]+(*drv)[pd[13]]-(*drv)[pd[31]]-(*drv)[pd[32]])
       + 0.375*(-(*drv)[pd[8]] - (*drv)[pd[14]])
       + 0.5*((*drv)[pd[9]] + (*drv)[pd[15]]) + 0.75*(*drv)[pd[33]]);
    (*drv)[cd[17]] = ((*drv)[pd[33]]);
    (*drv)[cd[18]] =
      (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
       + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]]) + 0.140625*(*drv)[pd[5]]
       + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[13]])
       + 0.5625*((*drv)[pd[31]] + (*drv)[pd[32]]));
    (*drv)[cd[19]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]]) + 0.015625*(*drv)[pd[5]]
       + 0.1875*((*drv)[pd[10]]+(*drv)[pd[16]]-(*drv)[pd[28]]-(*drv)[pd[29]])
       + 0.375*(-(*drv)[pd[11]] - (*drv)[pd[17]])
       + 0.5*((*drv)[pd[12]] + (*drv)[pd[18]]) + 0.75*(*drv)[pd[30]]);
    (*drv)[cd[20]] = ((*drv)[pd[30]]);
    (*drv)[cd[21]] =
      (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
       + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]]) + 0.140625*(*drv)[pd[5]]
       + 0.0625*(-(*drv)[pd[10]] - (*drv)[pd[16]])
       + 0.5625*((*drv)[pd[28]] + (*drv)[pd[29]]));
    (*drv)[cd[22]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]]) + 0.015625*(*drv)[pd[5]]
       + 0.125*((*drv)[pd[7]]
                -(*drv)[pd[8]]+(*drv)[pd[13]]-(*drv)[pd[14]]
                -(*drv)[pd[31]]-(*drv)[pd[32]])
       + 0.0625*((*drv)[pd[10]]+(*drv)[pd[16]]-(*drv)[pd[28]]-(*drv)[pd[29]])
       + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[33]])
       + 0.5*((*drv)[pd[23]] + (*drv)[pd[26]] + (*drv)[pd[34]]));
    (*drv)[cd[23]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]]) + 0.015625*(*drv)[pd[5]]
       + 0.0625*((*drv)[pd[7]]+(*drv)[pd[13]]-(*drv)[pd[31]]-(*drv)[pd[32]])
       + 0.125*((*drv)[pd[10]]-(*drv)[pd[11]]+(*drv)[pd[16]]-(*drv)[pd[17]]
                -(*drv)[pd[28]]-(*drv)[pd[29]])
       + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[30]])
       + 0.5*((*drv)[pd[24]] + (*drv)[pd[27]] + (*drv)[pd[34]]));
    (*drv)[cd[24]] = ((*drv)[pd[34]]);
    (*drv)[cd[25]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]]) - 0.234375*(*drv)[pd[5]]
       + 0.3125*((*drv)[pd[10]] - (*drv)[pd[29]]) + 0.0625*(*drv)[pd[16]]
       + 0.9375*(*drv)[pd[28]]);
    (*drv)[cd[26]] =
      (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
       - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
       + 0.09375*(*drv)[pd[6]]
       + 0.125*(-(*drv)[pd[10]] + (*drv)[pd[16]] - (*drv)[pd[17]])
       + 0.375*((*drv)[pd[11]] + (*drv)[pd[28]] - (*drv)[pd[29]])
       + 0.75*(*drv)[pd[30]]);
    (*drv)[cd[27]] = ((*drv)[pd[28]]);
    (*drv)[cd[28]] =
      (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
       + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]]) - 0.234375*(*drv)[pd[5]]
       + 0.3125*((*drv)[pd[7]] - (*drv)[pd[32]]) + 0.0625*(*drv)[pd[13]]
       + 0.9375*(*drv)[pd[31]]);
    (*drv)[cd[29]] =
      (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
       - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
       + 0.09375*(*drv)[pd[6]]
       + 0.125*(-(*drv)[pd[7]] + (*drv)[pd[13]] - (*drv)[pd[14]])
       + 0.375*((*drv)[pd[8]] + (*drv)[pd[31]] - (*drv)[pd[32]])
       + 0.75*(*drv)[pd[33]]);
    (*drv)[cd[30]] = ((*drv)[pd[31]]);
    (*drv)[cd[34]] =
      (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
       - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
       + 0.09375*(*drv)[pd[6]]
       + 0.0625*(-(*drv)[pd[7]]-(*drv)[pd[10]]+(*drv)[pd[13]]+(*drv)[pd[16]])
       - 0.125*(*drv)[pd[22]] + 0.375*(*drv)[pd[25]]
       + 0.1875*((*drv)[pd[28]]-(*drv)[pd[29]]+(*drv)[pd[31]]-(*drv)[pd[32]])
       + 0.75*(*drv)[pd[34]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cd);
    if (typ == 0)
    {
      /****************************************************************************/
      /*  parent of el_type 0                                                     */
      /****************************************************************************/

      (*drv)[cd[10]] =
        (-0.0390625*(*drv)[pd[0]] + 0.2734375*(*drv)[pd[1]]
         + 0.21875*(*drv)[pd[4]] - 0.546875*(*drv)[pd[5]]
         + 1.09375*(*drv)[pd[6]]);
      (*drv)[cd[11]] = (*drv)[pd[6]];
      (*drv)[cd[12]] =
        (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
         - 0.15625*(*drv)[pd[4]] + 0.703125*(*drv)[pd[5]]
         + 0.46875*(*drv)[pd[6]]);
      (*drv)[cd[25]] =
        (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
         - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[7]]
         + 0.3125*((*drv)[pd[13]] - (*drv)[pd[31]])
         + 0.9375*(*drv)[pd[32]]);
      (*drv)[cd[26]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] - (*drv)[pd[13]])
         + 0.375*((*drv)[pd[14]] - (*drv)[pd[31]] + (*drv)[pd[32]])
         + 0.75*(*drv)[pd[33]]);
      (*drv)[cd[27]] = (*drv)[pd[32]];
      (*drv)[cd[28]] =
        (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
         - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[10]]
         + 0.3125*((*drv)[pd[16]] - (*drv)[pd[28]])
         + 0.9375*(*drv)[pd[29]]);
      (*drv)[cd[29]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] - (*drv)[pd[16]])
         + 0.375*((*drv)[pd[17]] - (*drv)[pd[28]] + (*drv)[pd[29]])
         + 0.75*(*drv)[pd[30]]);
      (*drv)[cd[30]] = (*drv)[pd[29]];
      (*drv)[cd[34]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                   -(*drv)[pd[13]] - (*drv)[pd[16]])
         + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
         + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                   -(*drv)[pd[31]] + (*drv)[pd[32]])
         + 0.75*(*drv)[pd[34]]);
    }
    else
    {
      /****************************************************************************/
      /*  parent of el_type 1|2                                                   */
      /****************************************************************************/

      (*drv)[cd[10]] =
        (-0.0390625*(*drv)[pd[0]] + 0.2734375*(*drv)[pd[1]]
         + 0.21875*(*drv)[pd[4]] - 0.546875*(*drv)[pd[5]]
         + 1.09375*(*drv)[pd[6]]);
      (*drv)[cd[11]] = (*drv)[pd[6]];
      (*drv)[cd[12]] =
        (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
         - 0.15625*(*drv)[pd[4]] + 0.703125*(*drv)[pd[5]]
         + 0.46875*(*drv)[pd[6]]);
      (*drv)[cd[25]] =
        (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
         - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[10]]
         + 0.3125*((*drv)[pd[16]] - (*drv)[pd[28]])
         + 0.9375*(*drv)[pd[29]]);
      (*drv)[cd[26]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] - (*drv)[pd[16]])
         + 0.375*((*drv)[pd[17]] - (*drv)[pd[28]] + (*drv)[pd[29]])
         + 0.75*(*drv)[pd[30]]);
      (*drv)[cd[27]] = (*drv)[pd[29]];
      (*drv)[cd[28]] =
        (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
         + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
         - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[7]]
         + 0.3125*((*drv)[pd[13]] - (*drv)[pd[31]])
         + 0.9375*(*drv)[pd[32]]);
      (*drv)[cd[29]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] - (*drv)[pd[13]])
         + 0.375*((*drv)[pd[14]] - (*drv)[pd[31]] + (*drv)[pd[32]])
         + 0.75*(*drv)[pd[33]]);
      (*drv)[cd[30]] = (*drv)[pd[32]];
      (*drv)[cd[34]] =
        (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
         + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
         - 0.03125*(*drv)[pd[6]]
         + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                   - (*drv)[pd[13]] - (*drv)[pd[16]])
         + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
         + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                   - (*drv)[pd[31]] + (*drv)[pd[32]])
         + 0.75*(*drv)[pd[34]]);
    }

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);
      typ = list->getType(i);
      basFct->getLocalIndices(el, admin, pd);

      int lr_set = 0;
      if (list->getNeighbourElement(i, 0)  &&  list->getNeighbourNr(i, 0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i, 1)  &&  list->getNeighbourNr(i, 1) < i)
        lr_set += 2;

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(0), admin, cd);

      switch (lr_set)
      {
      case 0:
        (*drv)[cd[16]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.1875*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.375*(-(*drv)[pd[8]] - (*drv)[pd[14]])
           + 0.5*((*drv)[pd[9]] + (*drv)[pd[15]])
           + 0.75*(*drv)[pd[33]]);
        (*drv)[cd[17]] = (*drv)[pd[33]];
        (*drv)[cd[18]] =
          (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]])
           + 0.140625*(*drv)[pd[5]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[13]])
           + 0.5625*((*drv)[pd[31]] + (*drv)[pd[32]]));
        (*drv)[cd[19]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.1875*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.375*(-(*drv)[pd[11]] - (*drv)[pd[17]])
           + 0.5*((*drv)[pd[12]] + (*drv)[pd[18]])
           + 0.75*(*drv)[pd[30]]);
        (*drv)[cd[20]] = (*drv)[pd[30]];
        (*drv)[cd[21]] =
          (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]])
           + 0.140625*(*drv)[pd[5]]
           + 0.0625*(-(*drv)[pd[10]] - (*drv)[pd[16]])
           + 0.5625*((*drv)[pd[28]] + (*drv)[pd[29]]));
        (*drv)[cd[22]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] + (*drv)[pd[13]]
                    - (*drv)[pd[14]] - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.0625*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[33]])
           + 0.5*((*drv)[pd[23]] + (*drv)[pd[26]] + (*drv)[pd[34]]));
        (*drv)[cd[23]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.0625*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] + (*drv)[pd[16]]
                    - (*drv)[pd[17]] - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[30]])
           + 0.5*((*drv)[pd[24]] + (*drv)[pd[27]] + (*drv)[pd[34]]));
        (*drv)[cd[24]] = (*drv)[pd[34]];
        (*drv)[cd[25]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
           - 0.234375*(*drv)[pd[5]]
           + 0.3125*((*drv)[pd[10]] - (*drv)[pd[29]])
           + 0.0625*(*drv)[pd[16]] + 0.9375*(*drv)[pd[28]]);
        (*drv)[cd[26]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.125*(-(*drv)[pd[10]] + (*drv)[pd[16]] - (*drv)[pd[17]])
           + 0.375*((*drv)[pd[11]] + (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.75*(*drv)[pd[30]]);
        (*drv)[cd[27]] = (*drv)[pd[28]];
        (*drv)[cd[28]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
           - 0.234375*(*drv)[pd[5]]
           + 0.3125*((*drv)[pd[7]] - (*drv)[pd[32]])
           + 0.0625*(*drv)[pd[13]] + 0.9375*(*drv)[pd[31]]);
        (*drv)[cd[29]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.125*(-(*drv)[pd[7]] + (*drv)[pd[13]] - (*drv)[pd[14]])
           + 0.375*((*drv)[pd[8]] + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[33]]);
        (*drv)[cd[30]] = (*drv)[pd[31]];
        (*drv)[cd[34]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[10]]
                     + (*drv)[pd[13]] + (*drv)[pd[16]])
           - 0.125*(*drv)[pd[22]] + 0.375*(*drv)[pd[25]]
           + 0.1875*((*drv)[pd[28]] - (*drv)[pd[29]]
                     + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[34]]);
        break;
      case 1:
        (*drv)[cd[16]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.1875*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.375*(-(*drv)[pd[8]] - (*drv)[pd[14]])
           + 0.5*((*drv)[pd[9]] + (*drv)[pd[15]])
           + 0.75*(*drv)[pd[33]]);
        (*drv)[cd[17]] = (*drv)[pd[33]];
        (*drv)[cd[18]] =
          (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]])
           + 0.140625*(*drv)[pd[5]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[13]])
           + 0.5625*((*drv)[pd[31]] + (*drv)[pd[32]]));
        (*drv)[cd[22]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] + (*drv)[pd[13]]
                    - (*drv)[pd[14]] - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.0625*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[33]])
           + 0.5*((*drv)[pd[23]] + (*drv)[pd[26]] + (*drv)[pd[34]]));
        (*drv)[cd[23]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.0625*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] + (*drv)[pd[16]]
                    - (*drv)[pd[17]] - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[30]])
           + 0.5*((*drv)[pd[24]] + (*drv)[pd[27]] + (*drv)[pd[34]]));
        (*drv)[cd[24]] = (*drv)[pd[34]];
        (*drv)[cd[28]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
           - 0.234375*(*drv)[pd[5]]
           + 0.3125*((*drv)[pd[7]] - (*drv)[pd[32]])
           + 0.0625*(*drv)[pd[13]] + 0.9375*(*drv)[pd[31]]);
        (*drv)[cd[29]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.125*(-(*drv)[pd[7]] + (*drv)[pd[13]] - (*drv)[pd[14]])
           + 0.375*((*drv)[pd[8]] + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[33]]);
        (*drv)[cd[30]] = (*drv)[pd[31]];
        (*drv)[cd[34]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[10]]
                     + (*drv)[pd[13]] + (*drv)[pd[16]])
           - 0.125*(*drv)[pd[22]] + 0.375*(*drv)[pd[25]]
           + 0.1875*((*drv)[pd[28]] - (*drv)[pd[29]]
                     + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[34]]);
        break;
      case 2:
        (*drv)[cd[19]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.1875*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.375*(-(*drv)[pd[11]] - (*drv)[pd[17]])
           + 0.5*((*drv)[pd[12]] + (*drv)[pd[18]])
           + 0.75*(*drv)[pd[30]]);
        (*drv)[cd[20]] = (*drv)[pd[30]];
        (*drv)[cd[21]] =
          (0.0234375*((*drv)[pd[0]] + (*drv)[pd[1]])
           + 0.09375*(-(*drv)[pd[4]] - (*drv)[pd[6]])
           + 0.140625*(*drv)[pd[5]]
           + 0.0625*(-(*drv)[pd[10]] - (*drv)[pd[16]])
           + 0.5625*((*drv)[pd[28]] + (*drv)[pd[29]]));
        (*drv)[cd[22]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] + (*drv)[pd[13]]
                    - (*drv)[pd[14]] - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.0625*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[33]])
           + 0.5*((*drv)[pd[23]] + (*drv)[pd[26]] + (*drv)[pd[34]]));
        (*drv)[cd[23]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.0625*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] + (*drv)[pd[16]]
                    - (*drv)[pd[17]] - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[30]])
           + 0.5*((*drv)[pd[24]] + (*drv)[pd[27]] + (*drv)[pd[34]]));
        (*drv)[cd[24]] = (*drv)[pd[34]];
        (*drv)[cd[25]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
           - 0.234375*(*drv)[pd[5]]
           + 0.3125*((*drv)[pd[10]] - (*drv)[pd[29]])
           + 0.0625*(*drv)[pd[16]] + 0.9375*(*drv)[pd[28]]);
        (*drv)[cd[26]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.125*(-(*drv)[pd[10]] + (*drv)[pd[16]] - (*drv)[pd[17]])
           + 0.375*((*drv)[pd[11]] + (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.75*(*drv)[pd[30]]);
        (*drv)[cd[27]] = (*drv)[pd[28]];
        (*drv)[cd[34]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[10]]
                     + (*drv)[pd[13]] + (*drv)[pd[16]])
           - 0.125*(*drv)[pd[22]] + 0.375*(*drv)[pd[25]]
           + 0.1875*((*drv)[pd[28]] - (*drv)[pd[29]]
                     + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[34]]);
        break;
      case 3:
        (*drv)[cd[22]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] + (*drv)[pd[13]]
                    - (*drv)[pd[14]] - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.0625*((*drv)[pd[10]] + (*drv)[pd[16]]
                     - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[33]])
           + 0.5*((*drv)[pd[23]] + (*drv)[pd[26]] + (*drv)[pd[34]]));
        (*drv)[cd[23]] =
          (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
           + 0.03125*((*drv)[pd[4]] + (*drv)[pd[6]])
           + 0.015625*(*drv)[pd[5]]
           + 0.0625*((*drv)[pd[7]] + (*drv)[pd[13]]
                     - (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] + (*drv)[pd[16]]
                    - (*drv)[pd[17]] - (*drv)[pd[28]] - (*drv)[pd[29]])
           + 0.25*(-(*drv)[pd[22]] - (*drv)[pd[25]] + (*drv)[pd[30]])
           + 0.5*((*drv)[pd[24]] + (*drv)[pd[27]] + (*drv)[pd[34]]));
        (*drv)[cd[24]] = (*drv)[pd[34]];
        (*drv)[cd[34]] =
          (0.0234375*(*drv)[pd[0]] - 0.0390625*(*drv)[pd[1]]
           - 0.03125*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
           + 0.09375*(*drv)[pd[6]]
           + 0.0625*(-(*drv)[pd[7]] - (*drv)[pd[10]]
                     + (*drv)[pd[13]] + (*drv)[pd[16]])
           - 0.125*(*drv)[pd[22]] + 0.375*(*drv)[pd[25]]
           + 0.1875*((*drv)[pd[28]] - (*drv)[pd[29]]
                     + (*drv)[pd[31]] - (*drv)[pd[32]])
           + 0.75*(*drv)[pd[34]]);
        break;
      }

      /****************************************************************************/
      /*  values on child[1]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(1), admin, cd);

      if (typ == 0)
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[cd[25]] =
            (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
             - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[7]]
             + 0.3125*((*drv)[pd[13]] - (*drv)[pd[31]])
             + 0.9375*(*drv)[pd[32]]);
          (*drv)[cd[26]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] - (*drv)[pd[13]])
             + 0.375*((*drv)[pd[14]] - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[33]]);
          (*drv)[cd[27]] = (*drv)[pd[32]];
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        case 2:
          (*drv)[cd[28]] =
            (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
             - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[10]]
             + 0.3125*((*drv)[pd[16]] - (*drv)[pd[28]])
             + 0.9375*(*drv)[pd[29]]);
          (*drv)[cd[29]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] - (*drv)[pd[16]])
             + 0.375*((*drv)[pd[17]] - (*drv)[pd[28]] + (*drv)[pd[29]])
             + 0.75*(*drv)[pd[30]]);
          (*drv)[cd[30]] = (*drv)[pd[29]];
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        case 3:
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        }
      }
      else
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[cd[28]] =
            (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
             - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[7]]
             + 0.3125*((*drv)[pd[13]] - (*drv)[pd[31]])
             + 0.9375*(*drv)[pd[32]]);
          (*drv)[cd[29]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.125*((*drv)[pd[7]] - (*drv)[pd[8]] - (*drv)[pd[13]])
             + 0.375*((*drv)[pd[14]] - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[33]]);
          (*drv)[cd[30]] = (*drv)[pd[32]];
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        case 2:
          (*drv)[cd[25]] =
            (0.0390625*(-(*drv)[pd[0]] - (*drv)[pd[1]])
             + 0.15625*((*drv)[pd[4]] + (*drv)[pd[6]])
             - 0.234375*(*drv)[pd[5]] + 0.0625*(*drv)[pd[10]]
             + 0.3125*((*drv)[pd[16]] - (*drv)[pd[28]])
             + 0.9375*(*drv)[pd[29]]);
          (*drv)[cd[26]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.125*((*drv)[pd[10]] - (*drv)[pd[11]] - (*drv)[pd[16]])
             + 0.375*((*drv)[pd[17]] - (*drv)[pd[28]] + (*drv)[pd[29]])
             + 0.75*(*drv)[pd[30]]);
          (*drv)[cd[27]] = (*drv)[pd[29]];
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        case 3:
          (*drv)[cd[34]] =
            (-0.0390625*(*drv)[pd[0]] + 0.0234375*(*drv)[pd[1]]
             + 0.09375*(*drv)[pd[4]] - 0.046875*(*drv)[pd[5]]
             - 0.03125*(*drv)[pd[6]]
             + 0.0625*((*drv)[pd[7]] + (*drv)[pd[10]]
                       - (*drv)[pd[13]] - (*drv)[pd[16]])
             + 0.375*(*drv)[pd[22]] - 0.125*(*drv)[pd[25]]
             + 0.1875*(-(*drv)[pd[28]] + (*drv)[pd[29]]
                       - (*drv)[pd[31]] + (*drv)[pd[32]])
             + 0.75*(*drv)[pd[34]]);
          break;
        }
      }
    }
  }


  void Lagrange::coarseRestr0(DOFIndexed<double>* drv,
                              RCNeighbourList* list,
                              int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    int n0 = drv->getFeSpace()->getAdmin()->getNumberOfPreDofs(CENTER);
    Element* el = list->getElement(0);
    DegreeOfFreedom dof0 = el->getDof(0,n0);           /* 1st endpoint of refinement edge */
    DegreeOfFreedom dof1 = el->getDof(1,n0);           /* 2nd endpoint of refinement edge */
    DegreeOfFreedom dof_new = el->getChild(0)->getDof(basFct->getDim(), n0);
    /*  newest vertex is DIM */
    (*drv)[dof0] += 0.5*(*drv)[dof_new];
    (*drv)[dof1] += 0.5*(*drv)[dof_new];
  }

  void Lagrange::coarseRestr1(DOFIndexed<double>* drv,
                              RCNeighbourList* list,
                              int n, BasisFunction* basFct)
  {
    if (n < 1)
      return;

    int n0 = drv->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX);
    Element* el = list->getElement(0);

    // 1st endpoint of refinement edge
    DegreeOfFreedom dof0 = el->getDof(0, n0);
    // 2nd endpoint of refinement edge
    DegreeOfFreedom dof1 = el->getDof(1, n0);
    DegreeOfFreedom dof_new = el->getChild(0)->getDof(basFct->getDim(), n0);
    // newest vertex is DIM
    (*drv)[dof0] += 0.5 * (*drv)[dof_new];
    (*drv)[dof1] += 0.5 * (*drv)[dof_new];
  }

  void Lagrange::coarseRestr2_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr2_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin  = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof;
    basFct->getLocalIndices(el, admin, pdof);

    /****************************************************************************/
    /*  contributions of dofs located on child[0]                               */
    /****************************************************************************/

    int node = drv->getFeSpace()->getMesh()->getNode(VERTEX);
    int n0 = admin->getNumberOfPreDofs(VERTEX);
    DegreeOfFreedom cdof2 = el->getChild(0)->getDof(node + 2, n0);

    node = drv->getFeSpace()->getMesh()->getNode(EDGE);
    n0 = admin->getNumberOfPreDofs(EDGE);
    DegreeOfFreedom cdof3 = el->getChild(0)->getDof(node, n0);
    DegreeOfFreedom cdof4 = el->getChild(0)->getDof(node + 1, n0);

    (*drv)[pdof[0]] += 0.375*(*drv)[cdof3] - 0.125*(*drv)[cdof4];
    (*drv)[pdof[1]] += -0.125*((*drv)[cdof3] + (*drv)[cdof4]);
    (*drv)[pdof[3]] += 0.5*(*drv)[cdof4];
    (*drv)[pdof[4]] += 0.5*(*drv)[cdof4];
    (*drv)[pdof[5]] = (*drv)[cdof2] + 0.75*(*drv)[cdof3] + 0.25*(*drv)[cdof4];

    /****************************************************************************/
    /*  contributions of dofs located on child[1] and not on child[0]           */
    /****************************************************************************/

    cdof4 = el->getChild(1)->getDof(node + 1, n0);

    (*drv)[pdof[0]] += -0.125*(*drv)[cdof4];
    (*drv)[pdof[1]] += 0.375*(*drv)[cdof4];
    (*drv)[pdof[5]] += 0.75*(*drv)[cdof4];

    if (n > 1)
    {
      el = list->getElement(1);
      basFct->getLocalIndices(el, admin, pdof);

      /****************************************************************************/
      /*  first set those values not effected by previous element                 */
      /****************************************************************************/

      cdof4 = el->getChild(0)->getDof(node + 1, n0);
      (*drv)[pdof[3]] += 0.5*(*drv)[cdof4];
      (*drv)[pdof[4]] += 0.5*(*drv)[cdof4];

      /****************************************************************************/
      /*  and now, update the values in the refinement edge                       */
      /****************************************************************************/

      (*drv)[pdof[0]] += -0.125*(*drv)[cdof4];
      (*drv)[pdof[1]] += -0.125*(*drv)[cdof4];
      (*drv)[pdof[5]] += 0.25*(*drv)[cdof4];
    }
  }


  void Lagrange::coarseRestr2_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr2_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();

    vector<DegreeOfFreedom> pdof(10), cdof(10);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    int node0 = drv->getFeSpace()->getMesh()->getNode(EDGE);
    int n0 = admin->getNumberOfPreDofs(EDGE);
    int i, lr_set;

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.375*(*drv)[cdof[6]] + 0.125*(-(*drv)[cdof[8]] - (*drv)[cdof[9]]));
    (*drv)[pdof[1]] +=
      (0.125*(-(*drv)[cdof[6]] - (*drv)[cdof[8]] - (*drv)[cdof[9]]));
    (*drv)[pdof[4]] =
      ((*drv)[cdof[3]] + 0.75*(*drv)[cdof[6]]
       + 0.25*((*drv)[cdof[8]] + (*drv)[cdof[9]]));
    (*drv)[pdof[5]] +=
      (0.5*(*drv)[cdof[8]]);
    (*drv)[pdof[6]] +=
      (0.5*(*drv)[cdof[9]]);
    (*drv)[pdof[7]] +=
      (0.5*(*drv)[cdof[8]]);
    (*drv)[pdof[8]] +=
      (0.5*(*drv)[cdof[9]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);
    DegreeOfFreedom cdofi = el->getChild(1)->getDof(node0 + 2, n0);

    (*drv)[pdof[0]] += (-0.125*(*drv)[cdofi]);
    (*drv)[pdof[1]] += (0.375*(*drv)[cdofi]);
    (*drv)[pdof[4]] += (0.75*(*drv)[cdofi]);

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    for (i = 1; i < n; i++)
    {
      el = list->getElement(i);
      basFct->getLocalIndices(el, admin, pdof);

      lr_set = 0;
      if (list->getNeighbourElement(i,0)  &&  list->getNeighbourNr(i,0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i,1)  &&  list->getNeighbourNr(i,1) < i)
        lr_set += 2;

      TEST_EXIT_DBG(lr_set)("no values set on both neighbours\n");

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(0), admin, cdof);

      switch (lr_set)
      {
      case 1:
        cdofi = el->getChild(0)->getDof(node0 + 4, n0);
        (*drv)[pdof[0]] += (-0.125*(*drv)[cdofi]);
        (*drv)[pdof[1]] += (-0.125*(*drv)[cdofi]);
        (*drv)[pdof[4]] += (0.25*(*drv)[cdofi]);
        (*drv)[pdof[5]] += (0.5*(*drv)[cdofi]);
        (*drv)[pdof[7]] += (0.5*(*drv)[cdofi]);
        break;
      case 2:
        cdofi = el->getChild(0)->getDof(node0 + 5, n0);
        (*drv)[pdof[0]] += (-0.125*(*drv)[cdofi]);
        (*drv)[pdof[1]] += (-0.125*(*drv)[cdofi]);
        (*drv)[pdof[4]] += (0.25*(*drv)[cdofi]);
        (*drv)[pdof[6]] += (0.5*(*drv)[cdofi]);
        (*drv)[pdof[8]] += (0.5*(*drv)[cdofi]);
        break;
      }
    }
  }


  void Lagrange::coarseRestr3_1d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr3_1d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    int node, n0;
    DegreeOfFreedom dof9;

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(10), cdof(10);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.0625*(-(*drv)[cdof[2]] + (*drv)[cdof[6]] - (*drv)[cdof[9]])
       + 0.3125*(*drv)[cdof[3]]);
    (*drv)[pdof[1]] +=
      0.0625*(-(*drv)[cdof[2]] +  (*drv)[cdof[3]]
              + (*drv)[cdof[6]] + (*drv)[cdof[9]]);
    (*drv)[pdof[3]] +=
      -0.25*(*drv)[cdof[6]] - 0.125*(*drv)[cdof[9]];
    (*drv)[pdof[4]] +=
      0.5*(*drv)[cdof[6]];
    (*drv)[pdof[5]] +=
      0.5*(*drv)[cdof[6]];
    (*drv)[pdof[6]] +=
      -0.25*(*drv)[cdof[6]] + 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[7]] =
      (0.5625*(*drv)[cdof[2]] + 0.9375*(*drv)[cdof[3]] + (*drv)[cdof[4]]
       - 0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]]);
    (*drv)[pdof[8]] =
      (0.5625*(*drv)[cdof[2]] - 0.3125*(*drv)[cdof[3]]
       - 0.0625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[9]]);
    (*drv)[pdof[9]] =
      (*drv)[cdof[5]] + 0.5*(*drv)[cdof[6]] + 0.75*(*drv)[cdof[9]];

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] += 0.0625*(*drv)[cdof[6]] + 0.0625*(*drv)[cdof[9]];
    (*drv)[pdof[1]] += 0.3125*(*drv)[cdof[6]] - 0.0625*(*drv)[cdof[9]];
    (*drv)[pdof[3]] += 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[6]] += -0.125*(*drv)[cdof[9]];
    (*drv)[pdof[7]] += -0.3125*(*drv)[cdof[6]] - 0.18750*(*drv)[cdof[9]];
    (*drv)[pdof[8]] +=
      (*drv)[cdof[5]] + 0.9375*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[9]] += 0.75*(*drv)[cdof[9]];

    if (n <= 1)
      return;

    /****************************************************************************/
    /*  adjust the value on the neihgbour                                       */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neigh's child[0]                                              */
    /****************************************************************************/

    (*drv)[pdof[0]] += 0.0625*((*drv)[cdof[6]] - (*drv)[cdof[9]]);
    (*drv)[pdof[1]] += 0.0625*((*drv)[cdof[6]] + (*drv)[cdof[9]]);
    (*drv)[pdof[3]] += -0.25*(*drv)[cdof[6]] - 0.12500*(*drv)[cdof[9]];
    (*drv)[pdof[4]] += 0.5*(*drv)[cdof[6]];
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[6]];
    (*drv)[pdof[6]] += -0.25*(*drv)[cdof[6]] + 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[7]] += -0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[8]] += -0.0625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[9]] =
      (*drv)[cdof[5]] + 0.5*(*drv)[cdof[6]] + 0.75*(*drv)[cdof[9]];

    /****************************************************************************/
    /*  values on neigh's child[1]                                              */
    /****************************************************************************/

    node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    n0 = admin->getNumberOfPreDofs(CENTER);
    dof9 = el->getChild(1)->getDof(node, n0);

    (*drv)[pdof[0]] += 0.0625*(*drv)[dof9];
    (*drv)[pdof[1]] -= 0.0625*(*drv)[dof9];
    (*drv)[pdof[3]] += 0.375*(*drv)[dof9];
    (*drv)[pdof[6]] -= 0.125*(*drv)[dof9];
    (*drv)[pdof[7]] -= 0.1875*(*drv)[dof9];
    (*drv)[pdof[8]] += 0.1875*(*drv)[dof9];
    (*drv)[pdof[9]] += 0.75*(*drv)[dof9];
  }

  void Lagrange::coarseRestr3_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr3_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;


    int node, n0;
    DegreeOfFreedom dof9;

    const Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(10), cdof(10);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.0625*(-(*drv)[cdof[2]] + (*drv)[cdof[6]] - (*drv)[cdof[9]])
       + 0.3125*(*drv)[cdof[3]]);
    (*drv)[pdof[1]] +=
      0.0625*(-(*drv)[cdof[2]] +  (*drv)[cdof[3]]
              + (*drv)[cdof[6]] + (*drv)[cdof[9]]);
    (*drv)[pdof[3]] += -0.25*(*drv)[cdof[6]] - 0.125*(*drv)[cdof[9]];
    (*drv)[pdof[4]] +=  0.5*(*drv)[cdof[6]];
    (*drv)[pdof[5]] +=  0.5*(*drv)[cdof[6]];
    (*drv)[pdof[6]] += -0.25*(*drv)[cdof[6]] + 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[7]] =
      (0.5625*(*drv)[cdof[2]] + 0.9375*(*drv)[cdof[3]] + (*drv)[cdof[4]]
       - 0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]]);
    (*drv)[pdof[8]] =
      (0.5625*(*drv)[cdof[2]] - 0.3125*(*drv)[cdof[3]]
       - 0.0625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[9]]);
    (*drv)[pdof[9]] =
      (*drv)[cdof[5]] + 0.5*(*drv)[cdof[6]] + 0.75*(*drv)[cdof[9]];

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] += 0.0625*(*drv)[cdof[6]] + 0.0625*(*drv)[cdof[9]];
    (*drv)[pdof[1]] += 0.3125*(*drv)[cdof[6]] - 0.0625*(*drv)[cdof[9]];
    (*drv)[pdof[3]] += 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[6]] += -0.125*(*drv)[cdof[9]];
    (*drv)[pdof[7]] += -0.3125*(*drv)[cdof[6]] - 0.18750*(*drv)[cdof[9]];
    (*drv)[pdof[8]] +=
      (*drv)[cdof[5]] + 0.9375*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[9]] += 0.75*(*drv)[cdof[9]];

    if (n <= 1)  return;

    /****************************************************************************/
    /*  adjust the value on the neihgbour                                       */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neigh's child[0]                                              */
    /****************************************************************************/

    (*drv)[pdof[0]] += 0.0625*((*drv)[cdof[6]] - (*drv)[cdof[9]]);
    (*drv)[pdof[1]] += 0.0625*((*drv)[cdof[6]] + (*drv)[cdof[9]]);
    (*drv)[pdof[3]] += -0.25*(*drv)[cdof[6]] - 0.12500*(*drv)[cdof[9]];
    (*drv)[pdof[4]] += 0.5*(*drv)[cdof[6]];
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[6]];
    (*drv)[pdof[6]] += -0.25*(*drv)[cdof[6]] + 0.375*(*drv)[cdof[9]];
    (*drv)[pdof[7]] += -0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[8]] += -0.0625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[9]];
    (*drv)[pdof[9]] =
      (*drv)[cdof[5]] + 0.5*(*drv)[cdof[6]] + 0.75*(*drv)[cdof[9]];

    /****************************************************************************/
    /*  values on neigh's child[1]                                              */
    /****************************************************************************/

    node = drv->getFeSpace()->getMesh()->getNode(CENTER);
    n0 = admin->getNumberOfPreDofs(CENTER);
    dof9 = el->getChild(1)->getDof(node, n0);

    (*drv)[pdof[0]] += 0.0625*(*drv)[dof9];
    (*drv)[pdof[1]] -= 0.0625*(*drv)[dof9];
    (*drv)[pdof[3]] += 0.375*(*drv)[dof9];
    (*drv)[pdof[6]] -= 0.125*(*drv)[dof9];
    (*drv)[pdof[7]] -= 0.1875*(*drv)[dof9];
    (*drv)[pdof[8]] += 0.1875*(*drv)[dof9];
    (*drv)[pdof[9]] += 0.75*(*drv)[dof9];
  }

  void Lagrange::coarseRestr3_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr3_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    DegreeOfFreedom cdi;
    const Element* el = list->getElement(0);
    int typ = list->getType(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pd(20), cd(20);
    basFct->getLocalIndices(el, admin, pd);
    basFct->getLocalIndices(el->getChild(0), admin, cd);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pd[0]] +=
      (0.0625*(-(*drv)[cd[3]] + (*drv)[cd[12]] + (*drv)[cd[14]]
               + (*drv)[cd[16]] - (*drv)[cd[17]] - (*drv)[cd[18]])
       + 0.3125*(*drv)[cd[8]]);
    (*drv)[pd[1]] +=
      (0.0625*(-(*drv)[cd[3]] + (*drv)[cd[8]] + (*drv)[cd[12]] + (*drv)[cd[14]]
               + (*drv)[cd[16]] + (*drv)[cd[17]] + (*drv)[cd[18]]));
    (*drv)[pd[4]] =
      ((*drv)[cd[9]] + 0.5625*(*drv)[cd[3]] + 0.9375*(*drv)[cd[8]]
       + 0.0625*(-(*drv)[cd[12]] - (*drv)[cd[14]] - (*drv)[cd[16]])
       + 0.1875*((*drv)[cd[17]] + (*drv)[cd[18]]));
    (*drv)[pd[5]] =
      (0.5625*(*drv)[cd[3]] - 0.3125*(*drv)[cd[8]]
       + 0.0625*(-(*drv)[cd[12]] - (*drv)[cd[14]] - (*drv)[cd[16]])
       + 0.1875*(-(*drv)[cd[17]] - (*drv)[cd[18]]));
    (*drv)[pd[6]] +=
      (-0.25*(*drv)[cd[12]] - 0.125*(*drv)[cd[16]] + 0.375*(*drv)[cd[18]]);
    (*drv)[pd[7]] += 0.5*(*drv)[cd[12]];
    (*drv)[pd[8]] +=
      (-0.25*(*drv)[cd[14]] - 0.125*(*drv)[cd[16]] + 0.375*(*drv)[cd[17]]);
    (*drv)[pd[9]] += 0.5*(*drv)[cd[14]];
    (*drv)[pd[10]] +=
      (-0.25*(*drv)[cd[12]] + 0.125*(-(*drv)[cd[16]] - (*drv)[cd[18]]));
    (*drv)[pd[11]] += 0.5*(*drv)[cd[12]];
    (*drv)[pd[12]] +=
      (-0.25*(*drv)[cd[14]] + 0.125*(-(*drv)[cd[16]] - (*drv)[cd[17]]));
    (*drv)[pd[13]] += 0.5*(*drv)[cd[14]];
    (*drv)[pd[16]] += 0.5*(*drv)[cd[16]];
    (*drv)[pd[17]] += 0.5*(*drv)[cd[16]];
    (*drv)[pd[18]] =
      ((*drv)[cd[15]] + 0.5*(*drv)[cd[14]] + 0.25*(*drv)[cd[16]]
       + 0.75*(*drv)[cd[17]]);
    (*drv)[pd[19]] =
      ((*drv)[cd[13]] + 0.5*(*drv)[cd[12]] + 0.25*(*drv)[cd[16]]
       + 0.75*(*drv)[cd[18]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cd);

    if (typ == 0)
    {
      /****************************************************************************/
      /*  parent of el_type 0                                                     */
      /****************************************************************************/

      (*drv)[pd[0]] += 0.0625*((*drv)[cd[8]] + (*drv)[cd[17]] + (*drv)[cd[18]]);
      (*drv)[pd[1]] +=
        0.3125*(*drv)[cd[8]] + 0.0625*(-(*drv)[cd[17]] - (*drv)[cd[18]]);
      (*drv)[pd[4]] +=
        -0.3125*(*drv)[cd[8]] + 0.1875*(-(*drv)[cd[17]] - (*drv)[cd[18]]);
      (*drv)[pd[5]] +=
        ((*drv)[cd[9]] + 0.9375*(*drv)[cd[8]]
         + 0.1875*((*drv)[cd[17]] + (*drv)[cd[18]]));
      (*drv)[pd[6]] += -0.125*(*drv)[cd[17]];
      (*drv)[pd[8]] += -0.125*(*drv)[cd[18]];
      (*drv)[pd[10]] += 0.375*(*drv)[cd[17]];
      (*drv)[pd[12]] += 0.375*(*drv)[cd[18]];
      (*drv)[pd[18]] += 0.75*(*drv)[cd[18]];
      (*drv)[pd[19]] += 0.750000*(*drv)[cd[17]];
    }
    else
    {
      /****************************************************************************/
      /*  parent of el_type 1|2                                                   */
      /****************************************************************************/

      (*drv)[pd[0]] += 0.0625*((*drv)[cd[8]] + (*drv)[cd[17]] + (*drv)[cd[18]]);
      (*drv)[pd[1]] +=
        0.3125*(*drv)[cd[8]] + 0.0625*(-(*drv)[cd[17]] - (*drv)[cd[18]]);
      (*drv)[pd[4]] +=
        -0.3125*(*drv)[cd[8]] + 0.1875*(-(*drv)[cd[17]] - (*drv)[cd[18]]);
      (*drv)[pd[5]] +=
        ((*drv)[cd[9]] + 0.9375*(*drv)[cd[8]]
         + 0.1875*((*drv)[cd[17]] + (*drv)[cd[18]]));
      (*drv)[pd[6]] += -0.125*(*drv)[cd[18]];
      (*drv)[pd[8]] += -0.125*(*drv)[cd[17]];
      (*drv)[pd[10]] += 0.375*(*drv)[cd[18]];
      (*drv)[pd[12]] += 0.375*(*drv)[cd[17]];
      (*drv)[pd[18]] += 0.75*(*drv)[cd[17]];
      (*drv)[pd[19]] += 0.75*(*drv)[cd[18]];
    }

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    int node0 = drv->getFeSpace()->getMesh()->getNode(FACE);
    int n0 = admin->getNumberOfPreDofs(FACE);

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);
      typ = list->getType(i);
      basFct->getLocalIndices(el, admin, pd);

      int lr_set = 0;
      if (list->getNeighbourElement(i, 0) &&  list->getNeighbourNr(i, 0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i, 1) &&  list->getNeighbourNr(i, 1) < i)
        lr_set += 2;

      TEST_EXIT_DBG(lr_set)("no values set on both neighbours\n");

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(0), admin, cd);
      cdi = cd[16];

      switch (lr_set)
      {
      case 1:
        (*drv)[pd[0]] += 0.0625*((*drv)[cd[12]] + (*drv)[cdi] - (*drv)[cd[18]]);
        (*drv)[pd[1]] += 0.0625*((*drv)[cd[12]] + (*drv)[cdi] + (*drv)[cd[18]]);
        (*drv)[pd[4]] +=
          (0.0625*(-(*drv)[cd[12]] - (*drv)[cdi]) + 0.1875*(*drv)[cd[18]]);
        (*drv)[pd[5]] +=
          (0.0625*(-(*drv)[cd[12]] - (*drv)[cdi]) - 0.1875*(*drv)[cd[18]]);
        (*drv)[pd[6]] +=
          (-0.25*(*drv)[cd[12]] - 0.125*(*drv)[cdi] + 0.375*(*drv)[cd[18]]);
        (*drv)[pd[7]] += 0.5*(*drv)[cd[12]];
        (*drv)[pd[8]] += -0.125*(*drv)[cdi];
        (*drv)[pd[10]] +=
          -0.25*(*drv)[cd[12]] + 0.125*(-(*drv)[cdi] - (*drv)[cd[18]]);
        (*drv)[pd[11]] += 0.5*(*drv)[cd[12]];
        (*drv)[pd[12]] += -0.125*(*drv)[cdi];
        (*drv)[pd[16]] += 0.5*(*drv)[cdi];
        (*drv)[pd[17]] += 0.5*(*drv)[cdi];
        (*drv)[pd[18]] += 0.25*(*drv)[cdi];
        (*drv)[pd[19]] =
          ((*drv)[cd[13]] + 0.5*(*drv)[cd[12]]
           + 0.25*(*drv)[cdi] + 0.75*(*drv)[cd[18]]);
        break;
      case 2:
        (*drv)[pd[0]] += 0.0625*((*drv)[cd[14]] + (*drv)[cdi] - (*drv)[cd[17]]);
        (*drv)[pd[1]] += 0.0625*((*drv)[cd[14]] + (*drv)[cdi] + (*drv)[cd[17]]);
        (*drv)[pd[4]] +=
          0.0625*(-(*drv)[cd[14]] - (*drv)[cdi]) + 0.1875*(*drv)[cd[17]];
        (*drv)[pd[5]] +=
          0.0625*(-(*drv)[cd[14]] - (*drv)[cdi]) - 0.1875*(*drv)[cd[17]];
        (*drv)[pd[6]] += -0.125*(*drv)[cdi];
        (*drv)[pd[8]] +=
          -0.25*(*drv)[cd[14]] - 0.125*(*drv)[cdi] + 0.375*(*drv)[cd[17]];
        (*drv)[pd[9]] += 0.5*(*drv)[cd[14]];
        (*drv)[pd[10]] += -0.125*(*drv)[cdi];
        (*drv)[pd[12]] +=
          -0.25*(*drv)[cd[14]] + 0.125*(-(*drv)[cdi] - (*drv)[cd[17]]);
        (*drv)[pd[13]] += 0.5*(*drv)[cd[14]];
        (*drv)[pd[16]] += 0.5*(*drv)[cdi];
        (*drv)[pd[17]] += 0.5*(*drv)[cdi];
        (*drv)[pd[18]] =
          ((*drv)[cd[15]] + 0.5*(*drv)[cd[14]]
           + 0.25*(*drv)[cdi] + 0.75*(*drv)[cd[17]]);
        (*drv)[pd[19]] += 0.25*(*drv)[cdi];
        break;
      case 3:
        (*drv)[pd[0]] += 0.0625*(*drv)[cdi];
        (*drv)[pd[1]] += 0.0625*(*drv)[cdi];
        (*drv)[pd[4]] += -0.0625*(*drv)[cdi];
        (*drv)[pd[5]] += -0.0625*(*drv)[cdi];
        (*drv)[pd[6]] += -0.125*(*drv)[cdi];
        (*drv)[pd[8]] += -0.125*(*drv)[cdi];
        (*drv)[pd[10]] += -0.125*(*drv)[cdi];
        (*drv)[pd[12]] += -0.125*(*drv)[cdi];
        (*drv)[pd[16]] += 0.5*(*drv)[cdi];
        (*drv)[pd[17]] += 0.5*(*drv)[cdi];
        (*drv)[pd[18]] += 0.25*(*drv)[cdi];
        (*drv)[pd[19]] += 0.25*(*drv)[cdi];
        break;
      }

      /****************************************************************************/
      /*  values on child[1]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(1), admin, cd);

      if (typ == 0)
      {
        switch (lr_set)
        {
        case 1:
          cdi = el->getChild(1)->getDof(node0 + 1, n0);
          TEST_EXIT_DBG(cdi == cd[17])("cdi != cd[17]\n");
          (*drv)[pd[0]] += 0.0625*(*drv)[cdi];
          (*drv)[pd[1]] += -0.0625*(*drv)[cdi];
          (*drv)[pd[4]] += -0.1875*(*drv)[cdi];
          (*drv)[pd[5]] += 0.1875*(*drv)[cdi];
          (*drv)[pd[6]] += -0.125*(*drv)[cdi];
          (*drv)[pd[10]] += 0.375*(*drv)[cdi];
          (*drv)[pd[19]] += 0.75*(*drv)[cdi];
          break;
        case 2:
          cdi = el->getChild(1)->getDof(node0 + 2, n0);
          TEST_EXIT_DBG(cdi == cd[18])("cdi != cd[18]\n");
          (*drv)[pd[0]] += 0.0625*(*drv)[cdi];
          (*drv)[pd[1]] += -0.0625*(*drv)[cdi];
          (*drv)[pd[4]] += -0.1875*(*drv)[cdi];
          (*drv)[pd[5]] += 0.1875*(*drv)[cdi];
          (*drv)[pd[8]] += -0.125*(*drv)[cdi];
          (*drv)[pd[12]] += 0.375*(*drv)[cdi];
          (*drv)[pd[18]] += 0.75*(*drv)[cdi];
          break;
        }
      }
      else
      {
        switch (lr_set)
        {
        case 1:
          cdi = el->getChild(1)->getDof(node0 + 2, n0);
          TEST_EXIT_DBG(cdi == cd[18])("cdi != cd[18]\n");
          (*drv)[pd[0]] += 0.0625*(*drv)[cdi];
          (*drv)[pd[1]] += -0.0625*(*drv)[cdi];
          (*drv)[pd[4]] += -0.1875*(*drv)[cdi];
          (*drv)[pd[5]] += 0.1875*(*drv)[cdi];
          (*drv)[pd[6]] += -0.125*(*drv)[cdi];
          (*drv)[pd[10]] += 0.375*(*drv)[cdi];
          (*drv)[pd[19]] += 0.75*(*drv)[cdi];
          break;
        case 2:
          cdi = el->getChild(1)->getDof(node0 + 1, n0);
          TEST_EXIT_DBG(cdi == cd[17])("cdi != cd[17]\n");
          (*drv)[pd[0]] += 0.0625*(*drv)[cdi];
          (*drv)[pd[1]] += -0.0625*(*drv)[cdi];
          (*drv)[pd[4]] += -0.1875*(*drv)[cdi];
          (*drv)[pd[5]] += 0.1875*(*drv)[cdi];
          (*drv)[pd[8]] += -0.125*(*drv)[cdi];
          (*drv)[pd[12]] += 0.375*(*drv)[cdi];
          (*drv)[pd[18]] += 0.75*(*drv)[cdi];
          break;
        }
      }
    }
  }

  void Lagrange::coarseRestr4_1d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr4_1d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(15), cdof(15);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.2734375*(*drv)[cdof[3]]
       + 0.0390625*(-(*drv)[cdof[5]] - (*drv)[cdof[8]] - (*drv)[cdof[13]])
       + 0.0234375*((*drv)[cdof[6]] + (*drv)[cdof[12]]));
    (*drv)[pdof[1]] +=
      (0.0390625*(-(*drv)[cdof[3]] - (*drv)[cdof[8]]
                  - (*drv)[cdof[12]] - (*drv)[cdof[13]])
       + 0.0234375*((*drv)[cdof[5]] + (*drv)[cdof[6]]));
    (*drv)[pdof[3]] +=
      (0.0625*(-(*drv)[cdof[6]] + (*drv)[cdof[13]]) + 0.1875*(*drv)[cdof[8]]
       + 0.125*(*drv)[cdof[12]]);
    (*drv)[pdof[4]] += (-0.375*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]);
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[6]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[7]] += 0.375*(-(*drv)[cdof[8]] + (*drv)[cdof[12]]);
    (*drv)[pdof[8]] +=
      (-0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]
       + 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[9]] =
      ((*drv)[cdof[4]] + 1.09375*(*drv)[cdof[3]] + 0.46875*(*drv)[cdof[5]]
       - 0.09375*(*drv)[cdof[6]] + 0.15625*(*drv)[cdof[13]]
       + 0.03125*((*drv)[cdof[8]] - (*drv)[cdof[12]]));
    (*drv)[pdof[10]] =
      ((*drv)[cdof[2]] - 0.546875*(*drv)[cdof[3]] + 0.703125*(*drv)[cdof[5]]
       + 0.140625*(*drv)[cdof[6]] + 0.015625*(*drv)[cdof[8]]
       - 0.046875*(*drv)[cdof[12]] - 0.234375*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] =
      (0.21875*(*drv)[cdof[3]] + 0.15625*(-(*drv)[cdof[5]] + (*drv)[cdof[13]])
       + 0.09375*(-(*drv)[cdof[6]] + (*drv)[cdof[12]]) + 0.03125*(*drv)[cdof[8]]);
    (*drv)[pdof[12]] =
      ((*drv)[cdof[14]] + 0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]]
       + 0.375*(*drv)[cdof[12]] + 0.9375*(*drv)[cdof[13]]);
    (*drv)[pdof[13]] =
      (0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]] - 0.375*(*drv)[cdof[12]]
       - 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[14]] =
      ((*drv)[cdof[7]] + 0.75*((*drv)[cdof[8]] + (*drv)[cdof[12]]));

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] +=
      (0.0234375*(*drv)[cdof[6]]
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[12]] - (*drv)[cdof[13]]));
    (*drv)[pdof[1]] +=
      (0.0390625*(-(*drv)[cdof[6]] - (*drv)[cdof[12]])
       + 0.2734375*(*drv)[cdof[8]] + 0.0234375*(*drv)[cdof[13]]);
    (*drv)[pdof[3]] += 0.3125*(*drv)[cdof[12]] - 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[4]] += 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[7]] += -0.125*(*drv)[cdof[13]];
    (*drv)[pdof[8]] += 0.0625*(*drv)[cdof[12]] + 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[9]] +=
      (0.15625*(-(*drv)[cdof[6]] + (*drv)[cdof[12]]) + 0.21875*(*drv)[cdof[8]]
       + 0.09375*(*drv)[cdof[13]]);
    (*drv)[pdof[10]] +=
      (0.703125*(*drv)[cdof[6]] - 0.546875*(*drv)[cdof[8]]
       - 0.234375*(*drv)[cdof[12]] - 0.046875*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] +=
      ((*drv)[cdof[7]] + 0.46875*(*drv)[cdof[6]] + 1.09375*(*drv)[cdof[8]]
       + 0.15625*(*drv)[cdof[12]] - 0.03125*(*drv)[cdof[13]]);
    (*drv)[pdof[12]] += -0.3125*(*drv)[cdof[12]] - 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[13]] +=
      (*drv)[cdof[14]] + 0.9375*(*drv)[cdof[12]] + 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[14]] += 0.75*(*drv)[cdof[13]];

    if (n <= 1)  return;

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neighbour's child[0]                                          */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.0234375*((*drv)[cdof[6]] + (*drv)[cdof[12]])
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[13]]));
    (*drv)[pdof[1]] +=
      (0.0234375*(*drv)[cdof[6]]
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[12]] - (*drv)[cdof[13]]));
    (*drv)[pdof[3]] +=
      (0.0625*(-(*drv)[cdof[6]] + (*drv)[cdof[13]]) + 0.1875*(*drv)[cdof[8]]
       + 0.12500000*(*drv)[cdof[12]]);
    (*drv)[pdof[4]] += -0.375*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]];
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[6]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[7]] += 0.375*(-(*drv)[cdof[8]] + (*drv)[cdof[12]]);
    (*drv)[pdof[8]] +=
      (-0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]
       + 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[9]] +=
      (-0.09375*(*drv)[cdof[6]] + 0.03125*((*drv)[cdof[8]] - (*drv)[cdof[12]])
       + 0.15625*(*drv)[cdof[13]]);
    (*drv)[pdof[10]] +=
      (0.140625*(*drv)[cdof[6]] + 0.015625*(*drv)[cdof[8]]
       - 0.046875*(*drv)[cdof[12]] - 0.234375*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] +=
      (0.09375*(-(*drv)[cdof[6]] + (*drv)[cdof[12]])
       + 0.03125*(*drv)[cdof[8]] + 0.15625*(*drv)[cdof[13]]);
    (*drv)[pdof[12]] =
      ((*drv)[cdof[14]] + 0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]]
       + 0.375*(*drv)[cdof[12]] + 0.9375*(*drv)[cdof[13]]);
    (*drv)[pdof[13]] =
      (0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]] - 0.375*(*drv)[cdof[12]]
       - 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[14]] =
      (*drv)[cdof[7]] + 0.75*((*drv)[cdof[8]] + (*drv)[cdof[12]]);

    /****************************************************************************/
    /*  values on neighbour's child[1]                                          */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] += 0.0390625*(-(*drv)[cdof[12]] - (*drv)[cdof[13]]);
    (*drv)[pdof[1]] += -0.0390625*(*drv)[cdof[12]] + 0.0234375*(*drv)[cdof[13]];
    (*drv)[pdof[3]] += 0.3125*(*drv)[cdof[12]] - 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[4]] += 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[7]] += -0.125*(*drv)[cdof[13]];
    (*drv)[pdof[8]] += 0.0625*(*drv)[cdof[12]] + 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[9]] += 0.15625*(*drv)[cdof[12]] + 0.09375*(*drv)[cdof[13]];
    (*drv)[pdof[10]] += -0.234375*(*drv)[cdof[12]] - 0.046875*(*drv)[cdof[13]];
    (*drv)[pdof[11]] += 0.15625*(*drv)[cdof[12]] - 0.03125*(*drv)[cdof[13]];
    (*drv)[pdof[12]] += -0.3125*(*drv)[cdof[12]] - 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[13]] +=
      (*drv)[cdof[14]] + 0.9375*(*drv)[cdof[12]] + 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[14]] += 0.75*(*drv)[cdof[13]];
  }

  void Lagrange::coarseRestr4_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr4_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(15), cdof(15);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.2734375*(*drv)[cdof[3]]
       + 0.0390625*(-(*drv)[cdof[5]] - (*drv)[cdof[8]] - (*drv)[cdof[13]])
       + 0.0234375*((*drv)[cdof[6]] + (*drv)[cdof[12]]));
    (*drv)[pdof[1]] +=
      (0.0390625*(-(*drv)[cdof[3]] - (*drv)[cdof[8]]
                  - (*drv)[cdof[12]] - (*drv)[cdof[13]])
       + 0.0234375*((*drv)[cdof[5]] + (*drv)[cdof[6]]));
    (*drv)[pdof[3]] +=
      (0.0625*(-(*drv)[cdof[6]] + (*drv)[cdof[13]]) + 0.1875*(*drv)[cdof[8]]
       + 0.125*(*drv)[cdof[12]]);
    (*drv)[pdof[4]] += (-0.375*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]);
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[6]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[7]] += 0.375*(-(*drv)[cdof[8]] + (*drv)[cdof[12]]);
    (*drv)[pdof[8]] +=
      (-0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]
       + 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[9]] =
      ((*drv)[cdof[4]] + 1.09375*(*drv)[cdof[3]] + 0.46875*(*drv)[cdof[5]]
       - 0.09375*(*drv)[cdof[6]] + 0.15625*(*drv)[cdof[13]]
       + 0.03125*((*drv)[cdof[8]] - (*drv)[cdof[12]]));
    (*drv)[pdof[10]] =
      ((*drv)[cdof[2]] - 0.546875*(*drv)[cdof[3]] + 0.703125*(*drv)[cdof[5]]
       + 0.140625*(*drv)[cdof[6]] + 0.015625*(*drv)[cdof[8]]
       - 0.046875*(*drv)[cdof[12]] - 0.234375*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] =
      (0.21875*(*drv)[cdof[3]] + 0.15625*(-(*drv)[cdof[5]] + (*drv)[cdof[13]])
       + 0.09375*(-(*drv)[cdof[6]] + (*drv)[cdof[12]]) + 0.03125*(*drv)[cdof[8]]);
    (*drv)[pdof[12]] =
      ((*drv)[cdof[14]] + 0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]]
       + 0.375*(*drv)[cdof[12]] + 0.9375*(*drv)[cdof[13]]);
    (*drv)[pdof[13]] =
      (0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]] - 0.375*(*drv)[cdof[12]]
       - 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[14]] =
      ((*drv)[cdof[7]] + 0.75*((*drv)[cdof[8]] + (*drv)[cdof[12]]));

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] +=
      (0.0234375*(*drv)[cdof[6]]
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[12]] - (*drv)[cdof[13]]));
    (*drv)[pdof[1]] +=
      (0.0390625*(-(*drv)[cdof[6]] - (*drv)[cdof[12]])
       + 0.2734375*(*drv)[cdof[8]] + 0.0234375*(*drv)[cdof[13]]);
    (*drv)[pdof[3]] += 0.3125*(*drv)[cdof[12]] - 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[4]] += 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[7]] += -0.125*(*drv)[cdof[13]];
    (*drv)[pdof[8]] += 0.0625*(*drv)[cdof[12]] + 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[9]] +=
      (0.15625*(-(*drv)[cdof[6]] + (*drv)[cdof[12]]) + 0.21875*(*drv)[cdof[8]]
       + 0.09375*(*drv)[cdof[13]]);
    (*drv)[pdof[10]] +=
      (0.703125*(*drv)[cdof[6]] - 0.546875*(*drv)[cdof[8]]
       - 0.234375*(*drv)[cdof[12]] - 0.046875*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] +=
      ((*drv)[cdof[7]] + 0.46875*(*drv)[cdof[6]] + 1.09375*(*drv)[cdof[8]]
       + 0.15625*(*drv)[cdof[12]] - 0.03125*(*drv)[cdof[13]]);
    (*drv)[pdof[12]] += -0.3125*(*drv)[cdof[12]] - 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[13]] +=
      (*drv)[cdof[14]] + 0.9375*(*drv)[cdof[12]] + 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[14]] += 0.75*(*drv)[cdof[13]];

    if (n <= 1)  return;

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    el = list->getElement(1);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    /****************************************************************************/
    /*  values on neighbour's child[0]                                          */
    /****************************************************************************/

    (*drv)[pdof[0]] +=
      (0.0234375*((*drv)[cdof[6]] + (*drv)[cdof[12]])
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[13]]));
    (*drv)[pdof[1]] +=
      (0.0234375*(*drv)[cdof[6]]
       + 0.0390625*(-(*drv)[cdof[8]] - (*drv)[cdof[12]] - (*drv)[cdof[13]]));
    (*drv)[pdof[3]] +=
      (0.0625*(-(*drv)[cdof[6]] + (*drv)[cdof[13]]) + 0.1875*(*drv)[cdof[8]]
       + 0.12500000*(*drv)[cdof[12]]);
    (*drv)[pdof[4]] += -0.375*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]];
    (*drv)[pdof[5]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[6]] += 0.5*(*drv)[cdof[8]];
    (*drv)[pdof[7]] += 0.375*(-(*drv)[cdof[8]] + (*drv)[cdof[12]]);
    (*drv)[pdof[8]] +=
      (-0.0625*(*drv)[cdof[6]] + 0.1875*(*drv)[cdof[8]] - 0.125*(*drv)[cdof[12]]
       + 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[9]] +=
      (-0.09375*(*drv)[cdof[6]] + 0.03125*((*drv)[cdof[8]] - (*drv)[cdof[12]])
       + 0.15625*(*drv)[cdof[13]]);
    (*drv)[pdof[10]] +=
      (0.140625*(*drv)[cdof[6]] + 0.015625*(*drv)[cdof[8]]
       - 0.046875*(*drv)[cdof[12]] - 0.234375*(*drv)[cdof[13]]);
    (*drv)[pdof[11]] +=
      (0.09375*(-(*drv)[cdof[6]] + (*drv)[cdof[12]])
       + 0.03125*(*drv)[cdof[8]] + 0.15625*(*drv)[cdof[13]]);
    (*drv)[pdof[12]] =
      ((*drv)[cdof[14]] + 0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]]
       + 0.375*(*drv)[cdof[12]] + 0.9375*(*drv)[cdof[13]]);
    (*drv)[pdof[13]] =
      (0.5625*(*drv)[cdof[6]] - 0.1875*(*drv)[cdof[8]] - 0.375*(*drv)[cdof[12]]
       - 0.3125*(*drv)[cdof[13]]);
    (*drv)[pdof[14]] =
      (*drv)[cdof[7]] + 0.75*((*drv)[cdof[8]] + (*drv)[cdof[12]]);

    /****************************************************************************/
    /*  values on neighbour's child[1]                                          */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[0]] += 0.0390625*(-(*drv)[cdof[12]] - (*drv)[cdof[13]]);
    (*drv)[pdof[1]] += -0.0390625*(*drv)[cdof[12]] + 0.0234375*(*drv)[cdof[13]];
    (*drv)[pdof[3]] += 0.3125*(*drv)[cdof[12]] - 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[4]] += 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[7]] += -0.125*(*drv)[cdof[13]];
    (*drv)[pdof[8]] += 0.0625*(*drv)[cdof[12]] + 0.125*(*drv)[cdof[13]];
    (*drv)[pdof[9]] += 0.15625*(*drv)[cdof[12]] + 0.09375*(*drv)[cdof[13]];
    (*drv)[pdof[10]] += -0.234375*(*drv)[cdof[12]] - 0.046875*(*drv)[cdof[13]];
    (*drv)[pdof[11]] += 0.15625*(*drv)[cdof[12]] - 0.03125*(*drv)[cdof[13]];
    (*drv)[pdof[12]] += -0.3125*(*drv)[cdof[12]] - 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[13]] +=
      (*drv)[cdof[14]] + 0.9375*(*drv)[cdof[12]] + 0.375*(*drv)[cdof[13]];
    (*drv)[pdof[14]] += 0.75*(*drv)[cdof[13]];
  }

  void Lagrange::coarseRestr4_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseRestr4_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    vector<DegreeOfFreedom> pd(35), cd(35);
    const Element* el = list->getElement(0);
    int typ = list->getType(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    basFct->getLocalIndices(el, admin, pd);
    basFct->getLocalIndices(el->getChild(0), admin, cd);

    /****************************************************************************/
    /*  values on child[0]                                                      */
    /****************************************************************************/

    (*drv)[pd[0]] +=
      (0.2734375*(*drv)[cd[10]]
       + 0.0390625*(-(*drv)[cd[12]] - (*drv)[cd[16]] - (*drv)[cd[19]]
                    - (*drv)[cd[22]] - (*drv)[cd[23]] - (*drv)[cd[25]]
                    - (*drv)[cd[28]])
       + 0.0234375*((*drv)[cd[18]] + (*drv)[cd[21]] + (*drv)[cd[26]]
                    + (*drv)[cd[29]] + (*drv)[cd[34]]));
    (*drv)[pd[1]] +=
      (0.0390625*(-(*drv)[cd[10]] - (*drv)[cd[16]] - (*drv)[cd[19]]
                  - (*drv)[cd[22]] - (*drv)[cd[23]] - (*drv)[cd[25]]
                  - (*drv)[cd[26]] - (*drv)[cd[28]] - (*drv)[cd[29]]
                  - (*drv)[cd[34]])
       + 0.0234375*((*drv)[cd[12]] + (*drv)[cd[18]] + (*drv)[cd[21]]));
    (*drv)[pd[4]] =
      ((*drv)[cd[11]] + 1.09375*(*drv)[cd[10]] + 0.46875*(*drv)[cd[12]]
       + 0.03125*((*drv)[cd[16]] + (*drv)[cd[19]] + (*drv)[cd[22]]
                  + (*drv)[cd[23]] - (*drv)[cd[26]] - (*drv)[cd[29]]
                  - (*drv)[cd[34]])
       + 0.09375*(-(*drv)[cd[18]] - (*drv)[cd[21]])
       + 0.15625*((*drv)[cd[25]] + (*drv)[cd[28]]));
    (*drv)[pd[5]] =
      ((*drv)[cd[3]] - 0.546875*(*drv)[cd[10]]
       + 0.703125*(*drv)[cd[12]]
       + 0.015625*((*drv)[cd[16]] + (*drv)[cd[19]]
                   + (*drv)[cd[22]] + (*drv)[cd[23]])
       + 0.140625*((*drv)[cd[18]] + (*drv)[cd[21]])
       + 0.234375*(-(*drv)[cd[25]] - (*drv)[cd[28]])
       + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
    (*drv)[pd[6]] =
      (0.21875*(*drv)[cd[10]]
       + 0.15625*(-(*drv)[cd[12]] + (*drv)[cd[25]] + (*drv)[cd[28]])
       + 0.03125*((*drv)[cd[16]] + (*drv)[cd[19]]
                  + (*drv)[cd[22]] + (*drv)[cd[23]])
       + 0.09375*(-(*drv)[cd[18]] - (*drv)[cd[21]] + (*drv)[cd[26]]
                  + (*drv)[cd[29]] + (*drv)[cd[34]]));
    (*drv)[pd[7]] +=
      (0.1875*(*drv)[cd[16]]
       + 0.0625*(-(*drv)[cd[18]] + (*drv)[cd[23]] - (*drv)[cd[34]])
       + 0.125*((*drv)[cd[22]] - (*drv)[cd[29]]) + 0.3125*(*drv)[cd[28]]);
    (*drv)[pd[8]] +=
      (0.375*(-(*drv)[cd[16]] + (*drv)[cd[29]]) - 0.125*(*drv)[cd[22]]);
    (*drv)[pd[9]] += 0.5*(*drv)[cd[16]];
    (*drv)[pd[10]] +=
      (0.1875*(*drv)[cd[19]]
       + 0.0625*(-(*drv)[cd[21]] + (*drv)[cd[22]] - (*drv)[cd[34]])
       + 0.125*((*drv)[cd[23]] - (*drv)[cd[26]])
       + 0.3125*(*drv)[cd[25]]);
    (*drv)[pd[11]] +=
      (0.375*(-(*drv)[cd[19]] + (*drv)[cd[26]])
       - 0.125*(*drv)[cd[23]]);
    (*drv)[pd[12]] += 0.5*(*drv)[cd[19]];
    (*drv)[pd[13]] +=
      (0.1875*(*drv)[cd[16]]
       + 0.0625*(-(*drv)[cd[18]] + (*drv)[cd[23]]
                 + (*drv)[cd[28]] + (*drv)[cd[34]])
       + 0.125*((*drv)[cd[22]] + (*drv)[cd[29]]));
    (*drv)[pd[14]] +=
      (-0.375*(*drv)[cd[16]] + 0.125*(-(*drv)[cd[22]] - (*drv)[cd[29]]));
    (*drv)[pd[15]] += 0.5*(*drv)[cd[16]];
    (*drv)[pd[16]] +=
      (0.1875*(*drv)[cd[19]]
       + 0.0625*(-(*drv)[cd[21]] + (*drv)[cd[22]]
                 + (*drv)[cd[25]] + (*drv)[cd[34]])
       + 0.125*((*drv)[cd[23]] + (*drv)[cd[26]]));
    (*drv)[pd[17]] +=
      (-0.375*(*drv)[cd[19]] + 0.125*(-(*drv)[cd[23]] - (*drv)[cd[26]]));
    (*drv)[pd[18]] += 0.5*(*drv)[cd[19]];
    (*drv)[pd[22]] +=
      (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) - 0.125*(*drv)[cd[34]]);
    (*drv)[pd[23]] += 0.5*(*drv)[cd[22]];
    (*drv)[pd[24]] += 0.5*(*drv)[cd[23]];
    (*drv)[pd[25]] +=
      (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) + 0.375*(*drv)[cd[34]]);
    (*drv)[pd[26]] += 0.5*(*drv)[cd[22]];
    (*drv)[pd[27]] += 0.5*(*drv)[cd[23]];
    (*drv)[pd[28]] =
      ((*drv)[cd[27]] + 0.1875*(-(*drv)[cd[19]] + (*drv)[cd[34]])
       + 0.5625*(*drv)[cd[21]] - 0.0625*(*drv)[cd[22]]
       - 0.125*(*drv)[cd[23]] + 0.9375*(*drv)[cd[25]]
       + 0.375*(*drv)[cd[26]]);
    (*drv)[pd[29]] =
      (0.1875*(-(*drv)[cd[19]] - (*drv)[cd[34]])
       + 0.5625*(*drv)[cd[21]] - 0.0625*(*drv)[cd[22]]
       - 0.125*(*drv)[cd[23]] - 0.3125*(*drv)[cd[25]]
       - 0.375*(*drv)[cd[26]]);
    (*drv)[pd[30]] =
      ((*drv)[cd[20]] + 0.75*((*drv)[cd[19]] + (*drv)[cd[26]])
       + 0.25*(*drv)[cd[23]]);
    (*drv)[pd[31]] =
      ((*drv)[cd[30]] + 0.1875*(-(*drv)[cd[16]] + (*drv)[cd[34]])
       + 0.5625*(*drv)[cd[18]] - 0.125*(*drv)[cd[22]]
       - 0.0625*(*drv)[cd[23]] + 0.9375*(*drv)[cd[28]]
       + 0.375*(*drv)[cd[29]]);
    (*drv)[pd[32]] =
      (0.1875*(-(*drv)[cd[16]] - (*drv)[cd[34]])
       + 0.5625*(*drv)[cd[18]] - 0.125*(*drv)[cd[22]]
       - 0.0625*(*drv)[cd[23]] - 0.3125*(*drv)[cd[28]]
       - 0.375*(*drv)[cd[29]]);
    (*drv)[pd[33]] =
      ((*drv)[cd[17]] + 0.75*((*drv)[cd[16]] + (*drv)[cd[29]])
       + 0.25*(*drv)[cd[22]]);
    (*drv)[pd[34]] =
      ((*drv)[cd[24]] + 0.5*((*drv)[cd[22]] + (*drv)[cd[23]])
       + 0.75*(*drv)[cd[34]]);

    /****************************************************************************/
    /*  values on child[1]                                                      */
    /****************************************************************************/

    basFct->getLocalIndices(el->getChild(1), admin, cd);

    if (typ == 0)
    {
      /****************************************************************************/
      /*  parent of el_type 0                                                     */
      /****************************************************************************/

      (*drv)[pd[0]] +=
        (0.0390625*(-(*drv)[cd[10]] - (*drv)[cd[25]] - (*drv)[cd[26]]
                    - (*drv)[cd[28]] - (*drv)[cd[29]] - (*drv)[cd[34]])
         + 0.0234375*(*drv)[cd[12]]);
      (*drv)[pd[1]] +=
        (0.2734375*(*drv)[cd[10]]
         + 0.0390625*(-(*drv)[cd[12]] - (*drv)[cd[25]] - (*drv)[cd[28]])
         + 0.0234375*((*drv)[cd[26]] + (*drv)[cd[29]] + (*drv)[cd[34]]));
      (*drv)[pd[4]] +=
        (0.21875*(*drv)[cd[10]]
         + 0.15625*(-(*drv)[cd[12]] + (*drv)[cd[25]] + (*drv)[cd[28]])
         + 0.09375*((*drv)[cd[26]] + (*drv)[cd[29]] + (*drv)[cd[34]]));
      (*drv)[pd[5]] +=
        (-0.546875*(*drv)[cd[10]] + 0.703125*(*drv)[cd[12]]
         + 0.234375*(-(*drv)[cd[25]] - (*drv)[cd[28]])
         + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
      (*drv)[pd[6]] +=
        ((*drv)[cd[11]] + 1.09375*(*drv)[cd[10]] + 0.46875*(*drv)[cd[12]]
         + 0.15625*((*drv)[cd[25]] + (*drv)[cd[28]])
         + 0.03125*(-(*drv)[cd[26]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
      (*drv)[pd[7]] +=
        (0.0625*((*drv)[cd[25]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[26]]);
      (*drv)[pd[8]] += -0.125*(*drv)[cd[26]];
      (*drv)[pd[10]] +=
        (0.0625*((*drv)[cd[28]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[29]]);
      (*drv)[pd[11]] += -0.125*(*drv)[cd[29]];
      (*drv)[pd[13]] +=
        (0.3125*(*drv)[cd[25]] - 0.125*(*drv)[cd[26]]
         - 0.0625*(*drv)[cd[34]]);
      (*drv)[pd[14]] += 0.375*(*drv)[cd[26]];
      (*drv)[pd[16]] +=
        (0.3125*(*drv)[cd[28]] - 0.125*(*drv)[cd[29]]
         - 0.0625*(*drv)[cd[34]]);
      (*drv)[pd[17]] += 0.375*(*drv)[cd[29]];
      (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
      (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
      (*drv)[pd[28]] +=
        (-0.3125*(*drv)[cd[28]] - 0.375*(*drv)[cd[29]]
         - 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[29]] +=
        ((*drv)[cd[30]] + 0.9375*(*drv)[cd[28]] + 0.375*(*drv)[cd[29]]
         + 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[30]] += 0.75*(*drv)[cd[29]];
      (*drv)[pd[31]] +=
        (-0.3125*(*drv)[cd[25]] - 0.375*(*drv)[cd[26]]
         - 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[32]] +=
        ((*drv)[cd[27]] + 0.9375*(*drv)[cd[25]]
         + 0.375*(*drv)[cd[26]] + 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[33]] += 0.75*(*drv)[cd[26]];
      (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
    }
    else
    {
      /****************************************************************************/
      /*  parent of el_type 1|2                                                   */
      /****************************************************************************/

      (*drv)[pd[0]] +=
        (0.0390625*(-(*drv)[cd[10]] - (*drv)[cd[25]] - (*drv)[cd[26]]
                    - (*drv)[cd[28]] - (*drv)[cd[29]] - (*drv)[cd[34]])
         + 0.0234375*(*drv)[cd[12]]);
      (*drv)[pd[1]] +=
        (0.2734375*(*drv)[cd[10]]
         + 0.0390625*(-(*drv)[cd[12]] - (*drv)[cd[25]] - (*drv)[cd[28]])
         + 0.0234375*((*drv)[cd[26]] + (*drv)[cd[29]] + (*drv)[cd[34]]));
      (*drv)[pd[4]] +=
        (0.21875*(*drv)[cd[10]]
         + 0.15625*(-(*drv)[cd[12]] + (*drv)[cd[25]] + (*drv)[cd[28]])
         + 0.09375*((*drv)[cd[26]] + (*drv)[cd[29]] + (*drv)[cd[34]]));
      (*drv)[pd[5]] +=
        (-0.546875*(*drv)[cd[10]] + 0.703125*(*drv)[cd[12]]
         + 0.234375*(-(*drv)[cd[25]] - (*drv)[cd[28]])
         + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
      (*drv)[pd[6]] +=
        ((*drv)[cd[11]] + 1.09375*(*drv)[cd[10]] + 0.46875*(*drv)[cd[12]]
         + 0.15625*((*drv)[cd[25]] + (*drv)[cd[28]])
         + 0.03125*(-(*drv)[cd[26]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
      (*drv)[pd[7]] +=
        (0.0625*((*drv)[cd[28]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[29]]);
      (*drv)[pd[8]] += -0.125*(*drv)[cd[29]];
      (*drv)[pd[10]] +=
        (0.0625*((*drv)[cd[25]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[26]]);
      (*drv)[pd[11]] += -0.125*(*drv)[cd[26]];
      (*drv)[pd[13]] +=
        (0.3125*(*drv)[cd[28]] - 0.125*(*drv)[cd[29]]
         - 0.0625*(*drv)[cd[34]]);
      (*drv)[pd[14]] += 0.375*(*drv)[cd[29]];
      (*drv)[pd[16]] +=
        (0.3125*(*drv)[cd[25]] - 0.125*(*drv)[cd[26]]
         - 0.0625*(*drv)[cd[34]]);
      (*drv)[pd[17]] += 0.375*(*drv)[cd[26]];
      (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
      (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
      (*drv)[pd[28]] +=
        (-0.3125*(*drv)[cd[25]] - 0.375*(*drv)[cd[26]]
         - 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[29]] +=
        ((*drv)[cd[27]] + 0.9375*(*drv)[cd[25]] + 0.375*(*drv)[cd[26]]
         + 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[30]] += 0.75*(*drv)[cd[26]];
      (*drv)[pd[31]] +=
        (-0.3125*(*drv)[cd[28]] - 0.375*(*drv)[cd[29]]
         - 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[32]] +=
        ((*drv)[cd[30]] + 0.9375*(*drv)[cd[28]]
         + 0.375*(*drv)[cd[29]] + 0.1875*(*drv)[cd[34]]);
      (*drv)[pd[33]] += 0.75*(*drv)[cd[29]];
      (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
    }

    /****************************************************************************/
    /*   adjust neighbour values                                                */
    /****************************************************************************/

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);
      typ = list->getType(i);
      basFct->getLocalIndices(el, admin, pd);
      basFct->getLocalIndices(el->getChild(0), admin, cd);

      int lr_set = 0;
      if (list->getNeighbourElement(i,0) &&  list->getNeighbourNr(i,0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i,1) &&  list->getNeighbourNr(i,1) < i)
        lr_set += 2;

      TEST_EXIT_DBG(lr_set)("no values set on both neighbours\n");

      /****************************************************************************/
      /*  values on child[0]                                                      */
      /****************************************************************************/

      switch (lr_set)
      {
      case 1:
        (*drv)[pd[0]] +=
          (0.0390625*(-(*drv)[cd[16]] - (*drv)[cd[22]]
                      - (*drv)[cd[23]] - (*drv)[cd[28]])
           + 0.0234375*((*drv)[cd[18]] + (*drv)[cd[29]] + (*drv)[cd[34]]));
        (*drv)[pd[1]] +=
          (0.0390625*(-(*drv)[cd[16]] - (*drv)[cd[22]] - (*drv)[cd[23]]
                      - (*drv)[cd[28]] - (*drv)[cd[29]] - (*drv)[cd[34]])
           + 0.0234375*(*drv)[cd[18]]);
        (*drv)[pd[4]] +=
          (0.03125*((*drv)[cd[16]] + (*drv)[cd[22]] + (*drv)[cd[23]]
                    - (*drv)[cd[29]] - (*drv)[cd[34]])
           - 0.09375*(*drv)[cd[18]] + 0.15625*(*drv)[cd[28]]);
        (*drv)[pd[5]] +=
          (0.015625*((*drv)[cd[16]] + (*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.140625*(*drv)[cd[18]] - 0.234375*(*drv)[cd[28]]
           + 0.046875*(-(*drv)[cd[29]] - (*drv)[cd[34]]));
        (*drv)[pd[6]] +=
          (0.03125*((*drv)[cd[16]] + (*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.09375*(-(*drv)[cd[18]] + (*drv)[cd[29]] + (*drv)[cd[34]])
           + 0.15625*(*drv)[cd[28]]);
        (*drv)[pd[7]] +=
          (0.1875*(*drv)[cd[16]]
           + 0.0625*(-(*drv)[cd[18]] + (*drv)[cd[23]] - (*drv)[cd[34]])
           + 0.125*((*drv)[cd[22]] - (*drv)[cd[29]])
           + 0.3125*(*drv)[cd[28]]);
        (*drv)[pd[8]] +=
          (0.375*(-(*drv)[cd[16]] + (*drv)[cd[29]])
           - 0.125*(*drv)[cd[22]]);
        (*drv)[pd[9]] += 0.5*(*drv)[cd[16]];
        (*drv)[pd[10]] +=
          (0.0625*((*drv)[cd[22]] - (*drv)[cd[34]]) + 0.125*(*drv)[cd[23]]);
        (*drv)[pd[11]] += -0.125*(*drv)[cd[23]];
        (*drv)[pd[13]] +=
          (0.1875*(*drv)[cd[16]]
           + 0.0625*(-(*drv)[cd[18]] + (*drv)[cd[23]]
                     + (*drv)[cd[28]] + (*drv)[cd[34]])
           + 0.125*((*drv)[cd[22]] + (*drv)[cd[29]]));
        (*drv)[pd[14]] +=
          (-0.375*(*drv)[cd[16]]
           + 0.125*(-(*drv)[cd[22]] - (*drv)[cd[29]]));
        (*drv)[pd[15]] += 0.5*(*drv)[cd[16]];
        (*drv)[pd[16]] +=
          (0.0625*((*drv)[cd[22]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[23]]);
        (*drv)[pd[17]] += -0.125*(*drv)[cd[23]];
        (*drv)[pd[22]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) - 0.125*(*drv)[cd[34]]);
        (*drv)[pd[23]] += 0.5*(*drv)[cd[22]];
        (*drv)[pd[24]] += 0.5*(*drv)[cd[23]];
        (*drv)[pd[25]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) + 0.375*(*drv)[cd[34]]);
        (*drv)[pd[26]] += 0.5*(*drv)[cd[22]];
        (*drv)[pd[27]] += 0.5*(*drv)[cd[23]];
        (*drv)[pd[28]] +=
          (-0.0625*(*drv)[cd[22]] - 0.125*(*drv)[cd[23]]
           + 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[29]] +=
          (-0.0625*(*drv)[cd[22]] - 0.125*(*drv)[cd[23]]
           - 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[30]] += 0.25*(*drv)[cd[23]];
        (*drv)[pd[31]] =
          ((*drv)[cd[30]] + 0.1875*(-(*drv)[cd[16]] + (*drv)[cd[34]])
           + 0.5625*(*drv)[cd[18]] - 0.125*(*drv)[cd[22]]
           - 0.0625*(*drv)[cd[23]] + 0.9375*(*drv)[cd[28]]
           + 0.375*(*drv)[cd[29]]);
        (*drv)[pd[32]] =
          (0.1875*(-(*drv)[cd[16]] - (*drv)[cd[34]])
           + 0.5625*(*drv)[cd[18]] - 0.125*(*drv)[cd[22]]
           - 0.0625*(*drv)[cd[23]] - 0.3125*(*drv)[cd[28]]
           - 0.375*(*drv)[cd[29]]);
        (*drv)[pd[33]] =
          ((*drv)[cd[17]] + 0.75*((*drv)[cd[16]] + (*drv)[cd[29]])
           + 0.25*(*drv)[cd[22]]);
        (*drv)[pd[34]] =
          ((*drv)[cd[24]] + 0.5*((*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.75*(*drv)[cd[34]]);
        break;
      case 2:
        (*drv)[pd[0]] +=
          (0.0390625*(-(*drv)[cd[19]] - (*drv)[cd[22]]
                      - (*drv)[cd[23]] - (*drv)[cd[25]])
           + 0.0234375*((*drv)[cd[21]] + (*drv)[cd[26]] + (*drv)[cd[34]]));
        (*drv)[pd[1]] +=
          (0.0390625*(-(*drv)[cd[19]] - (*drv)[cd[22]] - (*drv)[cd[23]]
                      - (*drv)[cd[25]] - (*drv)[cd[26]] - (*drv)[cd[34]])
           + 0.0234375*(*drv)[cd[21]]);
        (*drv)[pd[4]] +=
          (0.03125*((*drv)[cd[19]] + (*drv)[cd[22]] + (*drv)[cd[23]]
                    - (*drv)[cd[26]] - (*drv)[cd[34]])
           - 0.09375*(*drv)[cd[21]] + 0.15625*(*drv)[cd[25]]);
        (*drv)[pd[5]] +=
          (0.015625*((*drv)[cd[19]] + (*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.140625*(*drv)[cd[21]] - 0.234375*(*drv)[cd[25]]
           + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[34]]));
        (*drv)[pd[6]] +=
          (0.03125*((*drv)[cd[19]] + (*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.09375*(-(*drv)[cd[21]] + (*drv)[cd[26]] + (*drv)[cd[34]])
           + 0.15625*(*drv)[cd[25]]);
        (*drv)[pd[7]] +=
          (0.125*(*drv)[cd[22]] + 0.0625*((*drv)[cd[23]] - (*drv)[cd[34]]));
        (*drv)[pd[8]] += -0.125*(*drv)[cd[22]];
        (*drv)[pd[10]] +=
          (0.1875*(*drv)[cd[19]]
           + 0.0625*(-(*drv)[cd[21]] + (*drv)[cd[22]] - (*drv)[cd[34]])
           + 0.125*((*drv)[cd[23]] - (*drv)[cd[26]])
           + 0.3125*(*drv)[cd[25]]);
        (*drv)[pd[11]] +=
          (0.375*(-(*drv)[cd[19]] + (*drv)[cd[26]]) - 0.125*(*drv)[cd[23]]);
        (*drv)[pd[12]] += 0.5*(*drv)[cd[19]];
        (*drv)[pd[13]] +=
          (0.125*(*drv)[cd[22]] + 0.0625*((*drv)[cd[23]] + (*drv)[cd[34]]));
        (*drv)[pd[14]] += -0.125*(*drv)[cd[22]];
        (*drv)[pd[16]] +=
          (0.1875*(*drv)[cd[19]]
           + 0.0625*(-(*drv)[cd[21]] + (*drv)[cd[22]]
                     + (*drv)[cd[25]] + (*drv)[cd[34]])
           + 0.125*((*drv)[cd[23]] + (*drv)[cd[26]]));
        (*drv)[pd[17]] +=
          (-0.375*(*drv)[cd[19]]
           + 0.125*(-(*drv)[cd[23]] - (*drv)[cd[26]]));
        (*drv)[pd[18]] += 0.5*(*drv)[cd[19]];
        (*drv)[pd[22]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) - 0.125*(*drv)[cd[34]]);
        (*drv)[pd[23]] += (0.5*(*drv)[cd[22]]);
        (*drv)[pd[24]] += (0.5*(*drv)[cd[23]]);
        (*drv)[pd[25]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) + 0.375*(*drv)[cd[34]]);
        (*drv)[pd[26]] += (0.5*(*drv)[cd[22]]);
        (*drv)[pd[27]] += (0.5*(*drv)[cd[23]]);
        (*drv)[pd[28]] =
          ((*drv)[cd[27]] + 0.1875*(-(*drv)[cd[19]] + (*drv)[cd[34]])
           + 0.5625*(*drv)[cd[21]] - 0.0625*(*drv)[cd[22]]
           - 0.125*(*drv)[cd[23]] + 0.9375*(*drv)[cd[25]]
           + 0.375*(*drv)[cd[26]]);
        (*drv)[pd[29]] =
          (0.1875*(-(*drv)[cd[19]] - (*drv)[cd[34]])
           + 0.5625*(*drv)[cd[21]] - 0.0625*(*drv)[cd[22]]
           - 0.125*(*drv)[cd[23]] - 0.3125*(*drv)[cd[25]]
           - 0.375*(*drv)[cd[26]]);
        (*drv)[pd[30]] =
          ((*drv)[cd[20]] + 0.75*((*drv)[cd[19]] + (*drv)[cd[26]])
           + 0.25*(*drv)[cd[23]]);
        (*drv)[pd[31]] +=
          (-0.125*(*drv)[cd[22]] - 0.0625*(*drv)[cd[23]]
           + 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[32]] +=
          (-0.125*(*drv)[cd[22]] - 0.0625*(*drv)[cd[23]]
           - 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[33]] += 0.25*(*drv)[cd[22]];
        (*drv)[pd[34]] =
          ((*drv)[cd[24]] + 0.5*((*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.75*(*drv)[cd[34]]);
        break;
      case 3:
        (*drv)[pd[0]] +=
          (0.0390625*(-(*drv)[cd[22]] - (*drv)[cd[23]])
           + 0.0234375*(*drv)[cd[34]]);
        (*drv)[pd[1]] +=
          (0.0390625*(-(*drv)[cd[22]] - (*drv)[cd[23]] - (*drv)[cd[34]]));
        (*drv)[pd[4]] +=
          (0.03125*((*drv)[cd[22]] + (*drv)[cd[23]] - (*drv)[cd[34]]));
        (*drv)[pd[5]] +=
          (0.015625*((*drv)[cd[22]] + (*drv)[cd[23]])
           - 0.046875*(*drv)[cd[34]]);
        (*drv)[pd[6]] +=
          (0.03125*((*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.09375*(*drv)[cd[34]]);
        (*drv)[pd[7]] +=
          (0.125*(*drv)[cd[22]] + 0.0625*((*drv)[cd[23]] - (*drv)[cd[34]]));
        (*drv)[pd[8]] += -0.125*(*drv)[cd[22]];
        (*drv)[pd[10]] +=
          (0.0625*((*drv)[cd[22]] - (*drv)[cd[34]]) + 0.125*(*drv)[cd[23]]);
        (*drv)[pd[11]] += -0.125*(*drv)[cd[23]];
        (*drv)[pd[13]] +=
          (0.125*(*drv)[cd[22]] + 0.0625*((*drv)[cd[23]] + (*drv)[cd[34]]));
        (*drv)[pd[14]] += -0.125*(*drv)[cd[22]];
        (*drv)[pd[16]] +=
          (0.0625*((*drv)[cd[22]] + (*drv)[cd[34]]) + 0.125*(*drv)[cd[23]]);
        (*drv)[pd[17]] += -0.125*(*drv)[cd[23]];
        (*drv)[pd[22]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) - 0.125*(*drv)[cd[34]]);
        (*drv)[pd[23]] += 0.5*(*drv)[cd[22]];
        (*drv)[pd[24]] += 0.5*(*drv)[cd[23]];
        (*drv)[pd[25]] +=
          (0.25*(-(*drv)[cd[22]] - (*drv)[cd[23]]) + 0.375*(*drv)[cd[34]]);
        (*drv)[pd[26]] += 0.5*(*drv)[cd[22]];
        (*drv)[pd[27]] += 0.5*(*drv)[cd[23]];
        (*drv)[pd[28]] +=
          (-0.0625*(*drv)[cd[22]] - 0.125*(*drv)[cd[23]]
           + 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[29]] +=
          (-0.0625*(*drv)[cd[22]] - 0.125*(*drv)[cd[23]]
           - 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[30]] += 0.25*(*drv)[cd[23]];
        (*drv)[pd[31]] +=
          (-0.125*(*drv)[cd[22]] - 0.0625*(*drv)[cd[23]]
           + 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[32]] +=
          (-0.125*(*drv)[cd[22]] - 0.0625*(*drv)[cd[23]]
           - 0.1875*(*drv)[cd[34]]);
        (*drv)[pd[33]] += 0.25*(*drv)[cd[22]];
        (*drv)[pd[34]] =
          ((*drv)[cd[24]] + 0.5*((*drv)[cd[22]] + (*drv)[cd[23]])
           + 0.75*(*drv)[cd[34]]);
        break;
      }

      /****************************************************************************/
      /*  values on child[1]                                                      */
      /****************************************************************************/

      basFct->getLocalIndices(el->getChild(1), admin, cd);

      if (typ == 0)
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[pd[0]] +=
            (0.0390625*(-(*drv)[cd[25]] - (*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[1]] +=
            (-0.0390625*(*drv)[cd[25]]
             + 0.0234375*((*drv)[cd[26]] + (*drv)[cd[34]]));
          (*drv)[pd[4]] +=
            (0.15625*(*drv)[cd[25]]
             + 0.09375*((*drv)[cd[26]] + (*drv)[cd[34]]));
          (*drv)[pd[5]] +=
            (-0.234375*(*drv)[cd[25]]
             + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[6]] +=
            (0.15625*(*drv)[cd[25]]
             + 0.03125*(-(*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[7]] +=
            (0.0625*((*drv)[cd[25]] + (*drv)[cd[34]])
             + 0.125*(*drv)[cd[26]]);
          (*drv)[pd[8]] += -0.125*(*drv)[cd[26]];
          (*drv)[pd[10]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[13]] +=
            (0.3125*(*drv)[cd[25]] - 0.125*(*drv)[cd[26]]
             - 0.0625*(*drv)[cd[34]]);
          (*drv)[pd[14]] += 0.375*(*drv)[cd[26]];
          (*drv)[pd[16]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[29]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[31]] +=
            (-0.3125*(*drv)[cd[25]] - 0.375*(*drv)[cd[26]]
             - 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[32]] +=
            ((*drv)[cd[27]] + 0.9375*(*drv)[cd[25]]
             + 0.375*(*drv)[cd[26]] + 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[33]] += 0.75*(*drv)[cd[26]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        case 2:
          (*drv)[pd[0]] +=
            (0.0390625*(-(*drv)[cd[28]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[1]] +=
            (-0.0390625*(*drv)[cd[28]]
             + 0.0234375*((*drv)[cd[29]] + (*drv)[cd[34]]));
          (*drv)[pd[4]] +=
            (0.15625*(*drv)[cd[28]]
             + 0.09375*((*drv)[cd[29]] + (*drv)[cd[34]]));
          (*drv)[pd[5]] +=
            (-0.234375*(*drv)[cd[28]]
             + 0.046875*(-(*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[6]] +=
            (0.15625*(*drv)[cd[28]]
             + 0.03125*(-(*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[7]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[10]] +=
            (0.0625*((*drv)[cd[28]] + (*drv)[cd[34]])
             + 0.125*(*drv)[cd[29]]);
          (*drv)[pd[11]] += -0.125*(*drv)[cd[29]];
          (*drv)[pd[13]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[16]] +=
            (0.3125*(*drv)[cd[28]] - 0.125*(*drv)[cd[29]]
             - 0.0625*(*drv)[cd[34]]);
          (*drv)[pd[17]] += 0.375*(*drv)[cd[29]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] +=
            (-0.3125*(*drv)[cd[28]] - 0.375*(*drv)[cd[29]]
             - 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[29]] +=
            ((*drv)[cd[30]] + 0.9375*(*drv)[cd[28]] + 0.375*(*drv)[cd[29]]
             + 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[30]] += 0.75*(*drv)[cd[29]];
          (*drv)[pd[31]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[32]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        case 3:
          (*drv)[pd[0]] += -0.0390625*(*drv)[cd[34]];
          (*drv)[pd[1]] += 0.0234375*(*drv)[cd[34]];
          (*drv)[pd[4]] += 0.09375*(*drv)[cd[34]];
          (*drv)[pd[5]] += -0.046875*(*drv)[cd[34]];
          (*drv)[pd[6]] += -0.03125*(*drv)[cd[34]];
          (*drv)[pd[7]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[10]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[13]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[16]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[29]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[31]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[32]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        }
      }
      else
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[pd[0]] +=
            (0.0390625*(-(*drv)[cd[28]] - (*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[1]] +=
            (-0.0390625*(*drv)[cd[28]]
             + 0.0234375*((*drv)[cd[29]] + (*drv)[cd[34]]));
          (*drv)[pd[4]] +=
            (0.15625*(*drv)[cd[28]]
             + 0.09375*((*drv)[cd[29]] + (*drv)[cd[34]]));
          (*drv)[pd[5]] +=
            (-0.234375*(*drv)[cd[28]]
             + 0.046875*(-(*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[6]] +=
            (0.15625*(*drv)[cd[28]]
             + 0.03125*(-(*drv)[cd[29]] - (*drv)[cd[34]]));
          (*drv)[pd[7]] +=
            (0.0625*((*drv)[cd[28]] + (*drv)[cd[34]])
             + 0.125*(*drv)[cd[29]]);
          (*drv)[pd[8]] += -0.125*(*drv)[cd[29]];
          (*drv)[pd[10]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[13]] +=
            (0.3125*(*drv)[cd[28]] - 0.125*(*drv)[cd[29]]
             - 0.0625*(*drv)[cd[34]]);
          (*drv)[pd[14]] += 0.375*(*drv)[cd[29]];
          (*drv)[pd[16]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[29]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[31]] +=
            (-0.3125*(*drv)[cd[28]] - 0.375*(*drv)[cd[29]]
             - 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[32]] +=
            ((*drv)[cd[30]] + 0.9375*(*drv)[cd[28]] + 0.375*(*drv)[cd[29]]
             + 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[33]] += 0.75*(*drv)[cd[29]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        case 2:
          (*drv)[pd[0]] +=
            (0.0390625*(-(*drv)[cd[25]] - (*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[1]] += (-0.0390625*(*drv)[cd[25]]
                            + 0.0234375*((*drv)[cd[26]] + (*drv)[cd[34]]));
          (*drv)[pd[4]] += (0.15625*(*drv)[cd[25]]
                            + 0.09375*((*drv)[cd[26]] + (*drv)[cd[34]]));
          (*drv)[pd[5]] += (-0.234375*(*drv)[cd[25]]
                            + 0.046875*(-(*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[6]] += (0.15625*(*drv)[cd[25]]
                            + 0.03125*(-(*drv)[cd[26]] - (*drv)[cd[34]]));
          (*drv)[pd[7]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[10]] += (0.0625*((*drv)[cd[25]] + (*drv)[cd[34]])
                             + 0.125*(*drv)[cd[26]]);
          (*drv)[pd[11]] += -0.125*(*drv)[cd[26]];
          (*drv)[pd[13]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[16]] += (0.3125*(*drv)[cd[25]] - 0.125*(*drv)[cd[26]]
                             - 0.0625*(*drv)[cd[34]]);
          (*drv)[pd[17]] += 0.375*(*drv)[cd[26]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] += (-0.3125*(*drv)[cd[25]] - 0.375*(*drv)[cd[26]]
                             - 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[29]] += ((*drv)[cd[27]] + 0.9375*(*drv)[cd[25]]
                             + 0.375*(*drv)[cd[26]] + 0.1875*(*drv)[cd[34]]);
          (*drv)[pd[30]] += 0.75*(*drv)[cd[26]];
          (*drv)[pd[31]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[32]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        case 3:
          (*drv)[pd[0]] += -0.0390625*(*drv)[cd[34]];
          (*drv)[pd[1]] += 0.0234375*(*drv)[cd[34]];
          (*drv)[pd[4]] += 0.09375*(*drv)[cd[34]];
          (*drv)[pd[5]] += -0.046875*(*drv)[cd[34]];
          (*drv)[pd[6]] += -0.03125*(*drv)[cd[34]];
          (*drv)[pd[7]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[10]] += 0.0625*(*drv)[cd[34]];
          (*drv)[pd[13]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[16]] += -0.0625*(*drv)[cd[34]];
          (*drv)[pd[22]] += 0.375*(*drv)[cd[34]];
          (*drv)[pd[25]] += -0.125*(*drv)[cd[34]];
          (*drv)[pd[28]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[29]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[31]] += -0.1875*(*drv)[cd[34]];
          (*drv)[pd[32]] += 0.1875*(*drv)[cd[34]];
          (*drv)[pd[34]] += 0.75*(*drv)[cd[34]];
          break;
        }
      }
    }
  }


  void Lagrange::coarseInter0(DOFIndexed<double>* drv, RCNeighbourList* list,
                              int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME_DBG("Lagrange::coarseInter0()");

    TEST_EXIT_DBG(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT_DBG(drv->getFeSpace()->getBasisFcts())
    ("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    const Mesh* mesh = drv->getFeSpace()->getMesh();

    // values on child[0]
    DegreeOfFreedom cdof = el->getChild(0)->getDof(mesh->getNode(CENTER) + 2,
                           admin->getNumberOfPreDofs(CENTER));
    DegreeOfFreedom pdof =
      el->getDof(mesh->getNode(CENTER) + 2, admin->getNumberOfPreDofs(CENTER));
    (*drv)[pdof] = (*drv)[cdof];
  }


  void Lagrange::coarseInter2_1d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME_DBG("Lagrange::coarseInter2_1d()");

    TEST_EXIT_DBG(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT_DBG(drv->getFeSpace()->getBasisFcts())
    ("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    Mesh* mesh = const_cast<Mesh*>(drv->getFeSpace()->getMesh());

    // values on child[0]
    DegreeOfFreedom cdof = el->getChild(0)->getDof(mesh->getNode(VERTEX) + 1,
                           admin->getNumberOfPreDofs(VERTEX));
    DegreeOfFreedom pdof = el->getDof(mesh->getNode(CENTER), admin->getNumberOfPreDofs(CENTER));
    (*drv)[pdof] = (*drv)[cdof];
  }


  void Lagrange::coarseInter2_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME("Lagrange::coarseInter2_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    const Mesh* mesh = drv->getFeSpace()->getMesh();

    // values on child[0]
    DegreeOfFreedom cdof = el->getChild(0)->getDof(mesh->getNode(VERTEX) + 2,
                           admin->getNumberOfPreDofs(VERTEX));
    DegreeOfFreedom pdof =
      el->getDof(mesh->getNode(EDGE) + 2, admin->getNumberOfPreDofs(EDGE));
    (*drv)[pdof] = (*drv)[cdof];
  }

  void Lagrange::coarseInter2_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME("Lagrange::coarseInter2_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    const Mesh* mesh = drv->getFeSpace()->getMesh();

    DegreeOfFreedom cdof = el->getChild(0)->getDof(mesh->getNode(VERTEX) + 3,
                           admin->getNumberOfPreDofs(VERTEX));
    DegreeOfFreedom pdof = el->getDof(mesh->getNode(EDGE), admin->getNumberOfPreDofs(EDGE));
    (*drv)[pdof] = (*drv)[cdof];
  }

  void Lagrange::coarseInter3_1d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME("Lagrange::coarseInter3_1d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    const Element* child = el->getChild(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    const Mesh* mesh = drv->getFeSpace()->getMesh();
    int node = mesh->getNode(EDGE);
    int n0 = admin->getNumberOfPreDofs(EDGE);

    // values on child[0]

    DegreeOfFreedom cdof, pdof;

    if (el->getDof(0, 0) < el->getDof(1, 0))
      pdof = el->getDof(node + 2, n0);
    else
      pdof = el->getDof(node + 2, n0 + 1);

    if (child->getDof(1, 0) < child->getDof(2, 0))
      cdof = child->getDof(node, n0 + 1);
    else
      cdof = child->getDof(node, n0);

    (*drv)[pdof] = (*drv)[cdof];

    if (child->getDof(2, 0) < child->getDof(0, 0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[el->getDof(mesh->getNode(CENTER), admin->getNumberOfPreDofs(CENTER))] =
      (*drv)[cdof];

    // values on child[1]

    child = el->getChild(1);

    if (el->getDof(0, 0) < el->getDof(1, 0))
      pdof = el->getDof(node + 2, n0 + 1);
    else
      pdof = el->getDof(node + 2, n0);

    if (child->getDof(2, 0) < child->getDof(0,0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[pdof] = (*drv)[cdof];

    if (n <= 1)
      return;

    // adjust neighbour values

    el = list->getElement(1);
    child = el->getChild(0);

    if (child->getDof(2,0) < child->getDof(0, 0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[el->getDof(mesh->getNode(CENTER), admin->getNumberOfPreDofs(CENTER))] =
      (*drv)[cdof];
  }

  void Lagrange::coarseInter3_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* /*basFct*/)
  {
    FUNCNAME("Lagrange::coarseInter3_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    const Element* child = el->getChild(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    const Mesh* mesh = drv->getFeSpace()->getMesh();
    int node = mesh->getNode(EDGE);
    int n0 = admin->getNumberOfPreDofs(EDGE);

    // values on child[0]

    DegreeOfFreedom cdof, pdof;

    if (el->getDof(0, 0) < el->getDof(1, 0))
      pdof = el->getDof(node + 2, n0);
    else
      pdof = el->getDof(node + 2, n0 + 1);

    if (child->getDof(1, 0) < child->getDof(2, 0))
      cdof = child->getDof(node, n0 + 1);
    else
      cdof = child->getDof(node, n0);

    (*drv)[pdof] = (*drv)[cdof];

    if (child->getDof(2, 0) < child->getDof(0, 0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[el->getDof(mesh->getNode(CENTER), admin->getNumberOfPreDofs(CENTER))] =
      (*drv)[cdof];

    // values on child[1]

    child = el->getChild(1);

    if (el->getDof(0, 0) < el->getDof(1, 0))
      pdof = el->getDof(node + 2, n0 + 1);
    else
      pdof = el->getDof(node + 2, n0);

    if (child->getDof(2, 0) < child->getDof(0, 0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[pdof] = (*drv)[cdof];

    if (n <= 1)
      return;

    // adjust neighbour values

    el = list->getElement(1);
    child = el->getChild(0);

    if (child->getDof(2, 0) < child->getDof(0, 0))
      cdof = child->getDof(node + 1, n0);
    else
      cdof = child->getDof(node + 1, n0 + 1);

    (*drv)[el->getDof(mesh->getNode(CENTER), admin->getNumberOfPreDofs(CENTER))] =
      (*drv)[cdof];
  }

  void Lagrange::coarseInter3_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseInter3_3d()");

    if (n < 1)
      return;

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    const Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();

    int node_e = drv->getFeSpace()->getMesh()->getNode(EDGE);
    int n0_e = admin->getNumberOfPreDofs(EDGE);
    int node_f = drv->getFeSpace()->getMesh()->getNode(FACE);
    int n0_f = admin->getNumberOfPreDofs(FACE);

    // values on child[0]

    const DegreeOfFreedom** pds = el->getDof();
    const DegreeOfFreedom** cds = el->getChild(0)->getDof();

    DegreeOfFreedom cd, pd;
    vector<DegreeOfFreedom> pd_o(20);
    basFct->getLocalIndices(el, admin, pd_o);
    pd = (pds[0][0] < pds[1][0]) ? pds[node_e][n0_e] : pds[node_e][n0_e + 1];
    cd = cds[0][0] < cds[3][0] ? cds[node_e + 2][n0_e + 1] : cds[node_e + 2][n0_e];
    (*drv)[pd] = (*drv)[cd];

    pd = el->getDof(node_f + 2, n0_f);
    cd = cds[2][0] < cds[3][0] ? cds[node_e + 5][n0_e + 1] : cds[node_e + 5][n0_e];
    (*drv)[pd] = (*drv)[cd];

    pd = el->getDof(node_f + 3, n0_f);
    cd = cds[1][0] < cds[3][0] ? cds[node_e + 4][n0_e + 1] : cds[node_e + 4][n0_e];
    (*drv)[pd] = (*drv)[cd];

    // values on child[1]

    cds = el->getChild(1)->getDof();
    pd = (pds[0][0] < pds[1][0]) ? pds[node_e][n0_e + 1] : pds[node_e][n0_e];
    cd = cds[0][0] < cds[3][0] ? cds[node_e + 2][n0_e + 1] : cds[node_e + 2][n0_e];
    (*drv)[pd] = (*drv)[cd];

    // adjust neighbour values

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);

      pds = el->getDof();
      cds = el->getChild(0)->getDof();

      int lr_set = 0;
      if (list->getNeighbourElement(i, 0)  &&  list->getNeighbourNr(i, 0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i, 1)  &&  list->getNeighbourNr(i, 1) < i)
        lr_set += 2;

      TEST_EXIT_DBG(lr_set)("no values set on both neighbours\n");

      switch (lr_set)
      {
      case 1:
        pd = el->getDof(node_f + 3, n0_f);
        cd = cds[1][0] < cds[3][0] ?
             cds[node_e + 4][n0_e + 1] : cds[node_e + 4][n0_e];
        (*drv)[pd] = (*drv)[cd];
        break;
      case 2:
        pd = el->getDof(node_f + 2, n0_f);
        cd = cds[2][0] < cds[3][0] ?
             cds[node_e + 5][n0_e + 1] : cds[node_e + 5][n0_e];
        (*drv)[pd] = (*drv)[cd];
        break;
      }
    }
  }

  void Lagrange::coarseInter4_1d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseInter4_1d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(15), cdof(15);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    // values on child[0]

    (*drv)[pdof[9]] = (*drv)[cdof[4]];
    (*drv)[pdof[10]] = (*drv)[cdof[2]];
    (*drv)[pdof[12]] = (*drv)[cdof[14]];
    (*drv)[pdof[14]] = (*drv)[cdof[7]];

    // values on child[1]

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[11]] = (*drv)[cdof[7]];
    (*drv)[pdof[13]] = (*drv)[cdof[14]];

    if (n <= 1)
      return;

    //   adjust neighbour values

    el = list->getElement(1);

    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    // values on neighbour's child[0]

    (*drv)[pdof[12]] = (*drv)[cdof[14]];
    (*drv)[pdof[14]] = (*drv)[cdof[7]];

    // values on neighbour's child[1]

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[13]] = (*drv)[cdof[14]];
  }

  void Lagrange::coarseInter4_2d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("Lagrange::coarseInter4_2d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pdof(15), cdof(15);
    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    // values on child[0]

    (*drv)[pdof[9]] = (*drv)[cdof[4]];
    (*drv)[pdof[10]] = (*drv)[cdof[2]];
    (*drv)[pdof[12]] = (*drv)[cdof[14]];
    (*drv)[pdof[14]] = (*drv)[cdof[7]];

    // values on child[1]

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[11]] = (*drv)[cdof[7]];
    (*drv)[pdof[13]] = (*drv)[cdof[14]];

    if (n <= 1)
      return;

    // adjust neighbour values

    el = list->getElement(1);

    basFct->getLocalIndices(el, admin, pdof);
    basFct->getLocalIndices(el->getChild(0), admin, cdof);

    // values on neighbour's child[0]

    (*drv)[pdof[12]] = (*drv)[cdof[14]];
    (*drv)[pdof[14]] = (*drv)[cdof[7]];

    // values on neighbour's child[1]

    basFct->getLocalIndices(el->getChild(1), admin, cdof);

    (*drv)[pdof[13]] = (*drv)[cdof[14]];
  }

  void Lagrange::coarseInter4_3d(DOFIndexed<double>* drv, RCNeighbourList* list,
                                 int n, BasisFunction* basFct)
  {
    FUNCNAME("void Lagrange::coarseInter4_3d()");

    TEST_EXIT(drv->getFeSpace())("No fe_space in dof_real_vec!\n");
    TEST_EXIT(drv->getFeSpace()->getBasisFcts())("No basis functions in fe_space!\n");

    if (n < 1)
      return;

    const Element* el = list->getElement(0);
    int typ = list->getType(0);
    const DOFAdmin* admin = drv->getFeSpace()->getAdmin();
    vector<DegreeOfFreedom> pd(35), cd(35);
    basFct->getLocalIndices(el, admin, pd);
    basFct->getLocalIndices(el->getChild(0), admin, cd);

    // values on child[0]

    (*drv)[pd[4]] = (*drv)[cd[11]];
    (*drv)[pd[5]] = (*drv)[cd[3]];
    (*drv)[pd[28]] = (*drv)[cd[27]];
    (*drv)[pd[30]] = (*drv)[cd[20]];
    (*drv)[pd[31]] = (*drv)[cd[30]];
    (*drv)[pd[33]] = (*drv)[cd[17]];
    (*drv)[pd[34]] = (*drv)[cd[24]];

    // values on child[1]

    basFct->getLocalIndices(el->getChild(1), admin, cd);

    if (typ == 0)
    {
      // parent of el_type 0

      (*drv)[pd[6]] = (*drv)[cd[11]];
      (*drv)[pd[29]] = (*drv)[cd[30]];
      (*drv)[pd[32]] = (*drv)[cd[27]];
    }
    else
    {
      // parent of el_type 1|2

      (*drv)[pd[6]] = (*drv)[cd[11]];
      (*drv)[pd[29]] = (*drv)[cd[27]];
      (*drv)[pd[32]] = (*drv)[cd[30]];
    }

    // adjust neighbour values

    for (int i = 1; i < n; i++)
    {
      el = list->getElement(i);
      typ = list->getType(i);
      basFct->getLocalIndices(el, admin, pd);

      int lr_set = 0;
      if (list->getNeighbourElement(i,0) &&  list->getNeighbourNr(i,0) < i)
        lr_set = 1;

      if (list->getNeighbourElement(i,1) &&  list->getNeighbourNr(i,1) < i)
        lr_set += 2;

      TEST_EXIT_DBG(lr_set)("no values set on both neighbours\n");

      // values on child[0]

      basFct->getLocalIndices(el->getChild(0), admin, cd);

      switch (lr_set)
      {
      case 1:
        (*drv)[pd[31]] = (*drv)[cd[30]];
        (*drv)[pd[33]] = (*drv)[cd[17]];
        (*drv)[pd[34]] = (*drv)[cd[24]];
        break;
      case 2:
        (*drv)[pd[28]] = (*drv)[cd[27]];
        (*drv)[pd[30]] = (*drv)[cd[20]];
        (*drv)[pd[34]] = (*drv)[cd[24]];
        break;
      case 3:
        (*drv)[pd[34]] = (*drv)[cd[24]];
        break;
      }

      // values on child[1]

      basFct->getLocalIndices(el->getChild(1), admin, cd);

      if (typ == 0)
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[pd[32]] = (*drv)[cd[27]];
          break;
        case 2:
          (*drv)[pd[29]] = (*drv)[cd[30]];
          break;
        }
      }
      else
      {
        switch (lr_set)
        {
        case 1:
          (*drv)[pd[32]] = (*drv)[cd[30]];
          break;
        case 2:
          (*drv)[pd[29]] = (*drv)[cd[27]];
          break;
        }
      }
    }
  }

}
