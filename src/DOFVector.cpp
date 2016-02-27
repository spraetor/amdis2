#include "DOFVector.hpp"
#include "Traverse.hpp"
// #include "DualTraverse.hpp"
#include "FixVec.hpp"
#include "ElementFunction.hpp"

namespace AMDiS
{
  template<>
  void DOFVector<double>::coarseRestrict(RCNeighbourList& list, int n)
  {
    FUNCNAME("DOFVector<double>::coarseRestrict()");

    switch (coarsenOperation)
    {
    case NO_OPERATION:
      return;
      break;
    case COARSE_RESTRICT:
      (const_cast<BasisFunction*>(feSpace->getBasisFcts()))->coarseRestr(this, &list, n);
      break;
    case COARSE_INTERPOL:
      TEST_EXIT_DBG(feSpace)("Should not happen!\n");
      TEST_EXIT_DBG(feSpace->getBasisFcts())("Shoud not happen!\n");

      (const_cast<BasisFunction*>(feSpace->getBasisFcts()))->coarseInter(this, &list, n);
      break;
    default:
      WARNING("Invalid coarsen operation \"%d\" in vector \"%s\"\n",
              coarsenOperation, this->name.c_str());
    }
  }


  template<>
  void DOFVector<double>::refineInterpol(RCNeighbourList& list, int n)
  {
    switch (refineOperation)
    {
    case NO_OPERATION:
      return;
      break;
    case REFINE_INTERPOL:
    default:
      (const_cast<BasisFunction*>(feSpace->getBasisFcts()))->refineInter(this, &list, n);
      break;
    }
  }


  template<>
  void DOFVector<WorldVector<double>>::refineInterpol(RCNeighbourList& list, int n)
  {
    if (refineOperation == NO_OPERATION)
      return;

    if (n < 1)
      return;

    Element* el = list.getElement(0);
    int n0 = feSpace->getAdmin()->getNumberOfPreDofs(VERTEX);
    DegreeOfFreedom dof0 = el->getDof(0, n0);
    DegreeOfFreedom dof1 = el->getDof(1, n0);
    DegreeOfFreedom dof_new = el->getChild(0)->getDof(feSpace->getMesh()->getDim(), n0);
    vec[dof_new] = vec[dof0];
    vec[dof_new] += vec[dof1];
    vec[dof_new] *= 0.5;
  }


  template<>
  double DOFVector<double>::evalAtPoint(
    WorldVector<double> const& p,
    ElInfo* oldElInfo) const
  {
    Mesh* mesh = feSpace->getMesh();
    const BasisFunction* basFcts = feSpace->getBasisFcts();

    int dim = mesh->getDim();
    int nBasFcts = basFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasFcts);
    DimVec<double> lambda(dim);

    ElInfo* elInfo = mesh->createNewElInfo();
    double value = 0.0;
    int inside = 0;

    if (oldElInfo && oldElInfo->getMacroElement())
    {
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, oldElInfo->getMacroElement(), NULL, NULL);
      delete oldElInfo;
    }
    else
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, NULL, NULL, NULL);

    if (oldElInfo)
      oldElInfo = elInfo;

    if (inside > 0)
    {
      basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), localIndices);
      DenseVector<double> uh(nBasFcts);
      for (int i = 0; i < nBasFcts; i++)
        uh[i] = operator[](localIndices[i]);
      value = basFcts->evalUh(lambda, uh);
    }
    else
    {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      value = 0.0;
#else
      ERROR_EXIT("Can not eval DOFVector at point p, because point is outside geometry.");
#endif
    }


    if (oldElInfo == NULL)
      delete elInfo;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(value);
#endif
    return value;
  }


  template<>
  WorldVector<double> DOFVector<WorldVector<double>>::evalAtPoint(
        WorldVector<double> const& p,
        ElInfo* oldElInfo) const
  {
    Mesh* mesh = feSpace->getMesh();
    const BasisFunction* basFcts = feSpace->getBasisFcts();

    int dim = mesh->getDim();
    int nBasFcts = basFcts->getNumber();

    std::vector<DegreeOfFreedom> localIndices(nBasFcts);
    DimVec<double> lambda(dim);
    ElInfo* elInfo = mesh->createNewElInfo();
    WorldVector<double> value(DEFAULT_SIZE, 0.0);
    int inside = 0;

    if (oldElInfo && oldElInfo->getMacroElement())
    {
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, oldElInfo->getMacroElement(), NULL, NULL);
      delete oldElInfo;
    }
    else
      inside = mesh->findElInfoAtPoint(p, elInfo, lambda, NULL, NULL, NULL);

    if (oldElInfo)
      oldElInfo = elInfo;

    if (inside > 0)
    {
      basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), localIndices);
      DenseVector<WorldVector<double>> uh(nBasFcts);
      for (int i = 0; i < nBasFcts; i++)
        uh[i] = operator[](localIndices[i]);
      value = basFcts->evalUh(lambda, uh);
    }
    else
    {
      ERROR_EXIT("Can not eval DOFVector at point p, because point is outside geometry.");
    }

    if (oldElInfo == NULL)
      delete elInfo;

    return value;
  }


  template<>
  double DOFVector<double>::IntOnBoundary(BoundaryType boundaryType,
                                          Quadrature* q) const
  {
    FUNCNAME("DOFVector<T>::IntOnBoundary()");

    Mesh* mesh = this->feSpace->getMesh();
    int dim = mesh->getDim();

    if (!q)
    {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(dim - 1, deg);
    }
    else
    {
      TEST_EXIT(q->getDim() == dim-1)
      ("Provided quadrature has wrong dimension for surfaceQuadrature!\n");
    }

    // create barycentric coords for each vertex of each side
    const Element* refElement = Global::getReferenceElement(dim);
    VectorOfFixVecs<DimVec<double>>** coords;
    coords = new VectorOfFixVecs<DimVec<double>>* [dim + 1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++)
    {
      coords[i] = new VectorOfFixVecs<DimVec<double>>(dim, dim, DEFAULT_VALUE,
          DimVec<double>(dim, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++)
      {
        int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim),
                    i, k);
        (*coords[i])[k][index] = 1.0;
      }
    }

    std::vector<SurfaceQuadrature*> quadSurfaces(dim + 1);
    for (int i = 0; i < dim + 1; i++)
      quadSurfaces[i] = new SurfaceQuadrature(q, *coords[i]);

    double result = 0.0;
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_BOUND | Mesh::FILL_COORDS);
    while (elInfo)
    {
      for (int face = 0; face < dim + 1; face++)
      {
        if (elInfo->getBoundary(face) == boundaryType)
        {
          int nPoints = quadSurfaces[face]->getNumPoints();
          DenseVector<double> uh_vec(nPoints);
          double det = elInfo->calcSurfaceDet(*coords[face]);
          double normT = 0.0;
          this->getVecAtQPs(elInfo, quadSurfaces[face], NULL, uh_vec);
          for (int iq = 0; iq < nPoints; iq++)
            normT += quadSurfaces[face]->getWeight(iq) * (uh_vec[iq]);
          result += det * normT;
        }
      }

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(result);
#endif

    return result;
  }


  template<>
  double DOFVector<WorldVector<double>>::IntOnBoundaryNormal(
                                       BoundaryType boundaryType, Quadrature* q) const
  {
    FUNCNAME("DOFVector<T>::IntOnBoundary()");

    Mesh* mesh = this->feSpace->getMesh();
    int dim = mesh->getDim();

    if (!q)
    {
      int deg = 2 * this->feSpace->getBasisFcts()->getDegree();
      q = Quadrature::provideQuadrature(dim - 1, deg);
    }
    else
    {
      TEST_EXIT(q->getDim() == dim-1)
      ("Provided quadrature has wrong dimension for surfaceQuadrature!\n");
    }

    // create barycentric coords for each vertex of each side
    const Element* refElement = Global::getReferenceElement(dim);
    VectorOfFixVecs<DimVec<double>>** coords;
    coords = new VectorOfFixVecs<DimVec<double>>* [dim + 1];

    // for all element sides
    for (int i = 0; i < dim + 1; i++)
    {
      coords[i] = new VectorOfFixVecs<DimVec<double>>(dim, dim, DEFAULT_VALUE,
          DimVec<double>(dim, 0.0));
      // for each vertex of the side
      for (int k = 0; k < dim; k++)
      {
        int index = refElement->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim),
                    i, k);
        (*coords[i])[k][index] = 1.0;
      }
    }

    std::vector<SurfaceQuadrature*> quadSurfaces(dim + 1);
    for (int i = 0; i < dim + 1; i++)
      quadSurfaces[i] = new SurfaceQuadrature(q, *coords[i]);

    double result = 0.0;
    WorldVector<double> normal;
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_BOUND | Mesh::FILL_COORDS);
    while (elInfo)
    {
      for (int face = 0; face < dim + 1; face++)
      {
        if (elInfo->getBoundary(face) == boundaryType)
        {
          elInfo->getNormal(face, normal);
          double det = elInfo->calcSurfaceDet(*coords[face]);

          int nPoints = quadSurfaces[face]->getNumPoints();
          DenseVector<WorldVector<double>> uh_vec(nPoints);
          WorldVector<double> normT;
          normT.set(0.0);
          this->getVecAtQPs(elInfo, quadSurfaces[face], NULL, uh_vec);
          for (int iq = 0; iq < nPoints; iq++)
            normT += quadSurfaces[face]->getWeight(iq) * (uh_vec[iq]);
          // scalar product between vector-valued solution and normal vector
          result += det * (normT * normal);
        }
      }

      elInfo = stack.traverseNext(elInfo);
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(result);
#endif

    return result;
  }


#if 0
  template<>
  void DOFVector<double>::interpol(DOFVector<double>* source, double factor)
  {

    const FiniteElemSpace* sourceFeSpace = source->getFeSpace();
    const BasisFunction* basisFcts = feSpace->getBasisFcts();
    const BasisFunction* sourceBasisFcts = sourceFeSpace->getBasisFcts();

    int nBasisFcts = basisFcts->getNumber();
    int nSourceBasisFcts = sourceBasisFcts->getNumber();

    this->set(0.0);

    std::vector<DegreeOfFreedom> myLocalIndices(nBasisFcts);
    DenseVector<double> sourceLocalCoeffs(nSourceBasisFcts);

    if (feSpace->getMesh() == sourceFeSpace->getMesh())
    {
      DimVec<double>* coords = NULL;
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(feSpace->getMesh(), -1,
                                           Mesh::CALL_LEAF_EL |
                                           Mesh::FILL_COORDS);

      if (basisFcts->isNodal())
      {
        while (elInfo)
        {
          Element* el = elInfo->getElement();

          basisFcts->getLocalIndices(el, feSpace->getAdmin(), myLocalIndices);
          source->getLocalVector(el, sourceLocalCoeffs);

          for (int i = 0; i < nBasisFcts; i++)
          {
            if (vec[myLocalIndices[i]] == 0.0)
            {
              coords = basisFcts->getCoords(i);
              vec[myLocalIndices[i]] = sourceBasisFcts->evalUh(*coords, sourceLocalCoeffs) * factor;
            }
          }
          elInfo = stack.traverseNext(elInfo);
        }
      }
      else
      {
        // nonnodal base:
        ElementFunctionWorld<double> F(source, factor);

        while (elInfo)
        {
          F.setElInfo(elInfo);

          Element* el = elInfo->getElement();

          basisFcts->getLocalIndices(el, feSpace->getAdmin(), myLocalIndices);

          DenseVector<double> rvec(nBasisFcts);
          std::function<double(WorldVector<double>)> f = F;
          basisFcts->interpol(elInfo, nBasisFcts, NULL, f, rvec);

          for (int i = 0; i < nBasisFcts; i++)
          {
            if (vec[myLocalIndices[i]] == 0.0)
            {
              vec[myLocalIndices[i]] = rvec[i];
            }
          }
          elInfo = stack.traverseNext(elInfo);
        }
      }
    }
    else
    {
      if (basisFcts->isNodal())
      {
        DimVec<double> coords2(feSpace->getMesh()->getDim());
        DualTraverse dualStack;
        ElInfo* elInfo1, *elInfo2;
        ElInfo* elInfoSmall, *elInfoLarge;
        WorldVector<double> worldVec;

        bool nextTraverse = dualStack.traverseFirst(feSpace->getMesh(),
                            sourceFeSpace->getMesh(),
                            -1, -1,
                            Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS,
                            Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS,
                            &elInfo1, &elInfo2,
                            &elInfoSmall, &elInfoLarge);
        while (nextTraverse)
        {
          basisFcts->getLocalIndices(elInfo1->getElement(), feSpace->getAdmin(),
                                     myLocalIndices);
          source->getLocalVector(elInfo2->getElement(), sourceLocalCoeffs);

          for (int i = 0; i < nBasisFcts; i++)
          {
            if (vec[myLocalIndices[i]] == 0.0)
            {
              elInfo1->coordToWorld(*(basisFcts->getCoords(i)), worldVec);
              elInfo2->worldToCoord(worldVec, coords2);

              bool isPositive = true;
              for (int j = 0; j < coords2.getSize(); j++)
              {
                if (coords2[j] < -0.00001)
                {
                  isPositive = false;
                  break;
                }
              }

              if (isPositive)
                vec[myLocalIndices[i]] =
                  sourceBasisFcts->evalUh(coords2, sourceLocalCoeffs);
            }
          }

          nextTraverse =
            dualStack.traverseNext(&elInfo1, &elInfo2, &elInfoSmall, &elInfoLarge);
        }
      }
      else
      {
        // nonnodal base
        ERROR_EXIT("not yet implemented\n");
      }
    }
  }


  template<>
  void DOFVector<WorldVector<double>>::interpol(DOFVector<WorldVector<double>>* v,
                                   double factor)
  {
    FUNCNAME("DOFVector<WorldVector<double> >::interpol()");
    WorldVector<double> nul(DEFAULT_VALUE,0.0);

    this->set(nul);

    DimVec<double>* coords = NULL;
    const FiniteElemSpace* vFeSpace = v->getFeSpace();

    if (feSpace == vFeSpace)
      WARNING("same FE-spaces\n");

    const BasisFunction* basFcts = feSpace->getBasisFcts();
    const BasisFunction* vBasFcts = vFeSpace->getBasisFcts();

    int nBasFcts = basFcts->getNumber();
    int vNumBasFcts = vBasFcts->getNumber();

    if (feSpace->getMesh() == vFeSpace->getMesh())
    {
      std::vector<DegreeOfFreedom> myLocalIndices(nBasFcts);
      DenseVector<WorldVector<double>> vLocalCoeffs(vNumBasFcts);
      Mesh* mesh = feSpace->getMesh();
      TraverseStack stack;
      ElInfo* elInfo =
        stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
      if (basFcts->isNodal())
      {
        while (elInfo)
        {
          Element* el = elInfo->getElement();
          basFcts->getLocalIndices(el, feSpace->getAdmin(), myLocalIndices);
          v->getLocalVector(el, vLocalCoeffs);

          for (int i = 0; i < nBasFcts; i++)
          {
            if (vec[myLocalIndices[i]] == nul)
            {
              coords = basFcts->getCoords(i);
              vec[myLocalIndices[i]] += vBasFcts->evalUh(*coords, vLocalCoeffs) * factor;
            }
          }
          elInfo = stack.traverseNext(elInfo);
        }
      }
      else
      {
        //nonnodal base:
        ElementFunctionWorld<WorldVector<double>> F(v, factor);

        while (elInfo)
        {
          F.setElInfo(elInfo);

          Element* el = elInfo->getElement();

          basFcts->getLocalIndices(el, feSpace->getAdmin(), myLocalIndices);
          v->getLocalVector(el, vLocalCoeffs);

          DenseVector<WorldVector<double>> rvec(nBasFcts);
          std::function<WorldVector<double>(WorldVector<double>)> f = F;
          basFcts->interpol(elInfo, nBasFcts, NULL, f, rvec);

          for (int i = 0; i < nBasFcts; i++)
          {
            if (vec[myLocalIndices[i]] == nul)
            {
              vec[myLocalIndices[i]] = rvec[i];
            }
          }
          elInfo = stack.traverseNext(elInfo);
        }
      }
    }
    else
    {
      ERROR_EXIT("not yet implemented for dual traverse\n");
    }
  }
#endif

  template<>
  WorldVector<DOFVector<double>*>* DOFVector<double>::getGradient(WorldVector<DOFVector<double>*>* grad) const
  {
    FUNCNAME_DBG("DOFVector<double>::getGradient()");

    Mesh* mesh = feSpace->getMesh();
    int dim = mesh->getDim();
    int dow = Global::getGeo(WORLD);

    const BasisFunction* basFcts = feSpace->getBasisFcts();

    DOFAdmin* admin = feSpace->getAdmin();

    // define result vector
    static WorldVector<DOFVector<double>*>* result = NULL; // TODO: REMOVE STATIC
    DOFVector<double>* null_ptr = NULL;

    if (grad)
    {
      result = grad;
    }
    else
    {
      if (!result)
      {
        result = new WorldVector<DOFVector<double>*>;

        result->set(null_ptr);
      }
      for (int i = 0; i < dow; i++)
      {
        if ((*result)[i] && (*result)[i]->getFeSpace() != feSpace)
        {
          delete (*result)[i];
          (*result)[i] = new DOFVector<double>(feSpace, "gradient");
        }
      }
    }

    // count number of nodes and dofs per node
    std::vector<int> nNodeDOFs;
    std::vector<int> nNodePreDofs;
    std::vector<DimVec<double>*> bary;

    int nNodes = 0;
    int nDofs = 0;

    for (int i = 0; i < dim + 1; i++)
    {
      GeoIndex geoIndex = INDEX_OF_DIM(i, dim);
      int numPositionNodes = mesh->getGeo(geoIndex);
      int numPreDofs = admin->getNumberOfPreDofs(i);
      for (int j = 0; j < numPositionNodes; j++)
      {
        int dofs = basFcts->getNumberOfDofs(geoIndex);
        nNodeDOFs.push_back(dofs);
        nDofs += dofs;
        nNodePreDofs.push_back(numPreDofs);
      }
      nNodes += numPositionNodes;
    }

    TEST_EXIT_DBG(nDofs == basFcts->getNumber())
    ("number of dofs != number of basis functions\n");

    for (int i = 0; i < nDofs; i++)
      bary.push_back(basFcts->getCoords(i));

    // traverse mesh
    std::vector<bool> visited(getUsedSize(), false);
    TraverseStack stack;
    Flag fillFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, fillFlag);
    WorldVector<double> grd;
    DenseVector<double> localUh(basFcts->getNumber());

    while (elInfo)
    {
      const DegreeOfFreedom** dof = elInfo->getElement()->getDof();
      getLocalVector(elInfo->getElement(), localUh);
      const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();

      int localDOFNr = 0;
      for (int i = 0; i < nNodes; i++)   // for all nodes
      {
        for (int j = 0; j < nNodeDOFs[i]; j++)   // for all dofs at this node
        {
          DegreeOfFreedom dofIndex = dof[i][nNodePreDofs[i] + j];

          if (!visited[dofIndex])
          {
            basFcts->evalGrdUh(*(bary[localDOFNr]), grdLambda, localUh, grd);

            for (int k = 0; k < dow; k++)
              (*(*result)[k])[dofIndex] = grd[k];

            visited[dofIndex] = true;
          }
          localDOFNr++;
        }
      }

      elInfo = stack.traverseNext(elInfo);
    }

    return result;
  }


  DOFVectorDOF::DOFVectorDOF() :
    DOFVector<DegreeOfFreedom>()
  {}


  void DOFVectorDOF::freeDOFContent(DegreeOfFreedom /*dof*/)
  {}


  WorldVector<DOFVector<double>*>* transform(DOFVector<WorldVector<double>>* vec,
      WorldVector<DOFVector<double>*>* res)
  {
    FUNCNAME_DBG("DOFVector<double>::transform()");

    TEST_EXIT_DBG(vec)("no vector\n");

    int dow = Global::getGeo(WORLD);
    static WorldVector<DOFVector<double>*>* result = NULL; // TODO: REMOVE STATIC

    if (!res && !result)
    {
      result = new WorldVector<DOFVector<double>*>;
      for (int i = 0; i < dow; i++)
        (*result)[i] = new DOFVector<double>(vec->getFeSpace(), "transform");
    }

    WorldVector<DOFVector<double>*>* r = res ? res : result;

    int vecSize = vec->getSize();
    for (int i = 0; i < vecSize; i++)
      for (int j = 0; j < dow; j++)
        (*((*r)[j]))[i] = (*vec)[i][j];

    return r;
  }

  // explicit template instatiation
  template class DOFVector<double>;
  // template class DOFVector<WorldVector<double> >;


} // end namespace AMDiS
