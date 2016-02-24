#include "PeriodicBC.h"
#include "BasisFunction.h"
#include "DOFVector.h"
#include "LeafData.h"
#include "DOFMatrix.h"
#include "Traverse.h"
#include "Boundary.h"
#include "VertexVector.h"

namespace AMDiS
{
  std::vector<PeriodicDOFMapping*> PeriodicDOFMapping::mappings;


  PeriodicDOFMapping*
  PeriodicDOFMapping::providePeriodicDOFMapping(BasisFunction const* fcts)
  {
    for (std::vector<PeriodicDOFMapping*>::iterator it = mappings.begin();
         it != mappings.end(); ++it)
      if ((*it)->basFcts == fcts)
        return *it;

    PeriodicDOFMapping* newMapping = new PeriodicDOFMapping(fcts);
    mappings.push_back(newMapping);
    return newMapping;
  }


  PeriodicDOFMapping::PeriodicDOFMapping(BasisFunction const* fcts)
    : basFcts(fcts)
  {
    FUNCNAME("PeriodicDOFMapping::PeriodicDOFMapping()");
    TEST_EXIT(basFcts->getDim() > 1)("dim == 1\n");
    int nBasFcts = basFcts->getNumber();
    DimVec<double>* lambda;
    for (int i = 0; i < nBasFcts; i++)
    {
      lambda = basFcts->getCoords(i);
      indexOfCoords[*lambda] = i;
    }
  }


  PeriodicDOFMapping::~PeriodicDOFMapping()
  {
    std::map<DimVec<int>, DegreeOfFreedom*, DimVecLess<int>>::iterator it;
    for (it = dofPermutation.begin(); it != dofPermutation.end(); ++it)
      if (it->second)
        delete [] it->second;
  }


  DegreeOfFreedom const* PeriodicDOFMapping::getDOFPermutation(DimVec<int> const& vertexPermutation)
  {
    FUNCNAME("PeriodicDOFMapping::getDOFPermutation()");

    if (dofPermutation[vertexPermutation] == NULL)
    {
      int dim = basFcts->getDim();
      int nBasFcts = basFcts->getNumber();
      int sum = 0;
      for (int i = 0; i < dim + 1; i++)
      {
        sum += i - vertexPermutation[i];
        TEST_EXIT(vertexPermutation[i] < dim + 1)
        ("invalid vertexPermuation\n");
      }
      TEST_EXIT(sum == 0)("invalid vertexPermutation\n");

      // create dof permutation
      DimVec<double>* lambda;
      DimVec<double> newLambda(dim);

      DegreeOfFreedom* mapping = new DegreeOfFreedom[nBasFcts];

      for (int i = 0; i < nBasFcts; i++)
      {
        lambda = basFcts->getCoords(i);
        for (int j = 0; j < dim + 1; j++)
          newLambda[vertexPermutation[j]] = (*lambda)[j];

        mapping[i] = indexOfCoords[newLambda];
      }

      dofPermutation[vertexPermutation] = mapping;
    }

    return dofPermutation[vertexPermutation];
  }


  PeriodicBC::PeriodicBC(BoundaryType type, FiniteElemSpace const* rowSpace)
    : BoundaryCondition(type, rowSpace, NULL),
      masterMatrix(NULL)
  {
    if (rowFeSpace->getMesh()->getDim() > 1)
      periodicDOFMapping =
        PeriodicDOFMapping::providePeriodicDOFMapping(rowFeSpace->getBasisFcts());
    else
      periodicDOFMapping = NULL;
  }


  PeriodicBC::~PeriodicBC()
  {}


  void PeriodicBC::initMatrix(DOFMatrix* matrix)
  {
    FUNCNAME("PeriodicBC::initMatrix()");

    if (!masterMatrix)
    {
      masterMatrix = matrix;
      Mesh* mesh = matrix->getRowFeSpace()->getMesh();
      associated = mesh->getPeriodicAssociations()[boundaryType];

      TEST_EXIT(associated)
      ("No associations for periodic boundary condition %d!\n", boundaryType);
    }
  }


  void PeriodicBC::fillBoundaryCondition(DOFMatrix* matrix,
                                         ElInfo* elInfo,
                                         DegreeOfFreedom const* dofIndices,
                                         BoundaryType const* /*localBound*/,
                                         int /*nBasFcts*/)
  {
    FUNCNAME_DBG("PeriodicBC::fillBoundaryCondition()");

    if (matrix != masterMatrix)
      return;

    int dim = rowFeSpace->getMesh()->getDim();
    if (dim == 1)
      return;

    DOFAdmin* admin = rowFeSpace->getAdmin();
    FixVec<int, WORLD> elFace(dim);
    FixVec<int, WORLD> neighFace(dim);
    DimVec<int> vertexPermutation(dim);
    const BasisFunction* basFcts = rowFeSpace->getBasisFcts();
    int num = basFcts->getNumber();
    Element* element = elInfo->getElement();
    DimVec<DegreeOfFreedom> periodicDOFs(dim - 1);
    GeoIndex sideGeoIndex = INDEX_OF_DIM(dim - 1, dim);
    std::vector<DegreeOfFreedom> neighIndices(num);

    for (int side = 0; side <= dim; side++)
    {
      if (elInfo->getBoundary(sideGeoIndex, side) == boundaryType)
      {
        for (int vertex = 0; vertex < dim; vertex++)
        {
          int index = element->getVertexOfPosition(sideGeoIndex, side, vertex);
          periodicDOFs[vertex] = (*associated)[dofIndices[index]];
        }

        Element* neigh = elInfo->getNeighbour(side);
        TEST_EXIT_DBG(neigh)("Wrong neighbour information at side %d!\n", side);

        basFcts->getLocalIndices(neigh, admin, neighIndices);

        int oppVertex = 0;
        for (int i = 0; i < dim + 1; i++)
        {
          // get vertex permutation
          if (i == side)
          {
            vertexPermutation[i] = 0;
          }
          else
          {
            DegreeOfFreedom periodicDOF =
              periodicDOFs[element->getPositionOfVertex(side, i)];

            int j = 0;
            for (; j < dim + 1; j++)
              if (neigh->getDof(j, 0) == periodicDOF)
                break;

            vertexPermutation[i] = j;
          }
          oppVertex += i - vertexPermutation[i];
        }
        vertexPermutation[side] = oppVertex;

        // get DOF permutation
        const DegreeOfFreedom* dofPermutation =
          periodicDOFMapping->getDOFPermutation(vertexPermutation);

        // set associated dofs
        for (int i = 0; i < num; i++)
          if ((*(basFcts->getCoords(i)))[side] == 0)
            (*associated)[dofIndices[i]] = neighIndices[dofPermutation[i]];
      }
    }
  }


  void PeriodicBC::exitMatrix(DOFMatrix* matrix)
  {
    FUNCNAME("PeriodicBC::exitMatrix()");

    TEST_EXIT(matrix)("No matrix\n");

    TEST_EXIT(associated)("No associated vector!\n");

    if (matrix == masterMatrix)
      masterMatrix = NULL;

    using namespace mtl;

    DOFAdmin* admin = rowFeSpace->getAdmin();
    std::vector<int> dofMap(admin->getUsedSize());
    for (int i = 0; i < admin->getUsedSize(); i++)
      dofMap[i] = (*associated)[i];

    // Compute reorder matrix (newRow and newCol yields transposed!!!)
    matrix::traits::reorder<>::type R= matrix::reorder(dofMap);
    DOFMatrix::base_matrix_type& A= matrix->getBaseMatrix(), C;

    C = R * A * trans(R) + A;
    A = 0.5 * C;
  }


  void PeriodicBC::exitVector(DOFVectorBase<double>* vector)
  {
    DOFIterator<double> vecIt(vector, USED_DOFS);
    Mesh* mesh = vector->getFeSpace()->getMesh();
    VertexVector* associated = mesh->getPeriodicAssociations()[boundaryType];

    for (vecIt.reset(); !vecIt.end(); ++vecIt)
    {
      DegreeOfFreedom dof = vecIt.getDOFIndex();
      DegreeOfFreedom newDOF = (*associated)[dof];

      if (dof < newDOF)
      {
        double entry = ((*vector)[dof] + (*vector)[newDOF]) * 0.5;
        (*vector)[dof] = entry;
        (*vector)[newDOF] = entry;
      }
    }
  }

} // end namespace AMDiS
