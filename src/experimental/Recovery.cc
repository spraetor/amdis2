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


#include "Recovery.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MeshDistributor.h"
#endif

namespace AMDiS
{

  RecoveryStructure& RecoveryStructure::operator=(const RecoveryStructure& rhs)
  {
    if (rhs.coords)
    {
      if (!coords)
        coords = new WorldVector<double>;
      *coords = *rhs.coords;
    }
    else
    {
      if (coords)
      {
        delete coords;
        coords = NULL;
      }
    }

    if (rhs.A)
    {
      if (!A)
        A = new Matrix<double>(rhs.A->getNumRows(), rhs.A->getNumCols());
      *A = *rhs.A;
    }
    else
    {
      if (A)
      {
        delete A;
        A = NULL;
      }
    }

    if (rhs.rec_uh)
    {
      if (!rec_uh)
        rec_uh = new Vector<double>(rhs.rec_uh->getSize());
      *rec_uh = *rhs.rec_uh;
    }
    else
    {
      if (rec_uh)
      {
        delete rec_uh;
        rec_uh = NULL;
      }
    }

    if (rhs.rec_grdUh)
    {
      if (!rec_grdUh)
        rec_grdUh = new Vector<WorldVector<double>>(rhs.rec_grdUh->getSize());
      *rec_grdUh = *rhs.rec_grdUh;
    }
    else
    {
      if (rec_grdUh)
      {
        delete rec_grdUh;
        rec_grdUh = NULL;
      }
    }

    if (rhs.neighbors)
    {
      if (!neighbors)
        neighbors = new std::set<DegreeOfFreedom>;
      *neighbors = *rhs.neighbors;
    }
    else
    {
      if (neighbors)
      {
        delete neighbors;
        neighbors = NULL;
      }
    }

    return *this;
  }


  void RecoveryStructure::print()
  {
    FUNCNAME("RecoveryStructure::print()");

    std::cout << std::endl;

    MSG("Coordinates of the node: ");
    std::cout << *coords << std::endl;

    if (A)
    {
      MSG("Interior vertex: printing system information.\n\n");

      int n = A->getNumRows();

      MSG("System matrix:\n");
      for (int i = 0; i < n; i++)
      {
        MSG("( ");
        for (int j = 0; j < i; j++)
          std::cout << "* ";
        for (int j = i; j < n; j++)
          std::cout << (*A)[i][j] << " ";
        std::cout << ")" << std::endl;
      }

      MSG("Right hand side:\n");
      for (int i = 0; i < n; i++)
      {
        MSG("( ");
        if (rec_grdUh)
          std::cout << (*rec_grdUh)[i];
        else
          std::cout << (*rec_uh)[i];
        std::cout << " )" << std::endl;
      }

      if (neighbors)
      {
        MSG("Printing neighbors vertices\n\n");

        MSG("Number of neighbors: ");
        std::cout << neighbors->size() << std::endl << std::endl;

        MSG("List of neighbors: ");
        std::set<DegreeOfFreedom>::const_iterator setIterator;

        for (setIterator = neighbors->begin(); setIterator != neighbors->end();
             ++setIterator)
          std::cout << " " << *setIterator;

        std::cout << std::endl << std::endl;
      }
    }
    else
    {
      MSG("Boundary vertex or not a vertex node: printing interior neighbors\n\n");

      MSG("Number of neighbors: ");
      std::cout << neighbors->size() << std::endl << std::endl;

      MSG("List of neighbors: ");
      std::set<DegreeOfFreedom>::const_iterator setIterator;

      for (setIterator  = neighbors->begin(); setIterator != neighbors->end();
           ++setIterator)
        std::cout << " " << *setIterator;

      std::cout << std::endl << std::endl;
    }

    WAIT;
  }


  void Recovery::set_feSpace(const FiniteElemSpace* fe_space)
  {
    if (!feSpace || feSpace != fe_space)
    {
      if (struct_vec)
      {
        delete struct_vec;
        struct_vec = NULL;
      }

      feSpace = fe_space;

      // create new structure vector
      struct_vec = new DOFVector<RecoveryStructure>(feSpace, "struct vec");
    }
  }


  int Recovery::set_exponents(int degree)
  {
    FUNCNAME("Recovery::set_exponents()");

    int dow = Global::getGeo(WORLD);
    int number_monomials = degree + 1;

    // Computing number of monomials.
    if (dow > 1)
    {
      number_monomials *= (degree + 2);
      number_monomials /= 2;
    }

    if (dow == 3)
    {
      number_monomials *= (degree + 3);
      number_monomials /= 3;
    }

    // Allocating memory.
    if (n_monomials != number_monomials)
    {
      n_monomials = number_monomials;
      exponents.resize(n_monomials);

      if (matrix_fcts)
        matrix_fcts->resize(n_monomials, n_monomials);
      else
        matrix_fcts = new Matrix<Monomial*>(n_monomials, n_monomials);
    }

    // Setting vector of exponents.
    int count = 0;

    switch (dow)
    {
    case 1:     // 1D monomials.
      for (int i = 0; i <= degree; i++)
        exponents[count++][0] = i;
      break;

    case 2:     // 2D monomials.
      for (int i = 0; i <= degree; i++)
        for (int j = 0; j <= i; j++)
        {
          exponents[count][0] = i - j;
          exponents[count++][1] = j;
        }
      break;

    case 3:     // 3D monomials.
      for (int i = 0; i <= degree; i++)
        for (int j = 0; j <= i; j++)
          for (int k = 0; k <= j; k++)
          {
            exponents[count][0] = i - j;
            exponents[count][1] = j - k;
            exponents[count++][2] = k;
          }
      break;

    default:
      ERROR_EXIT("Which dimension have your world???\n");
    }

    TEST_EXIT_DBG(count == n_monomials)("There must be an error!\n");

    // Setting matrix of monomials.
    WorldVector<int> sum;

    for (int i = 0; i < n_monomials; i++)
      for (int j = i; j < n_monomials; j++)
      {
        // Computing exponent vector of monomial.
        sum = exponents[i] + exponents[j];
        (*matrix_fcts)[i][j] = new Monomial(sum);
      }

    return n_monomials;
  }


  void Recovery::compute_integrals(DOFVector<double>* uh, ElInfo* elInfo,
                                   RecoveryStructure* rec_struct,
                                   AbstractFunction<double, WorldVector<double>>* f_vec,
                                   AbstractFunction<double, double>* f_scal,
                                   DOFVector<double>* aux_vec)
  {
    FUNCNAME_DBG("Recovery::compute_integrals()");

    TEST_EXIT_DBG(!(f_vec && f_scal))("Only one diffusion function, please!\n");

    WorldVector<double> vec_sum;
    WorldVector<double> quad_pts;  // For world coordinates of quadrature points.

    int deg_f = 0;
    if (f_vec)
      deg_f = f_vec->getDegree();
    if (f_scal)
      deg_f = f_scal->getDegree();

    if (gradient)
      deg_f += feSpace->getBasisFcts()->getDegree() - 1;
    else
      deg_f += feSpace->getBasisFcts()->getDegree();

    for (int i = 0; i < n_monomials; i++)
    {
      // Computing contributions to system matrix.
      for (int j = i; j < n_monomials; j++)
      {
        double sum  = 0.0;
        Quadrature* quad = Quadrature::provideQuadrature(Global::getGeo(WORLD),
                           (*matrix_fcts)[i][j]->getDegree());
        int n_points = quad->getNumPoints();

        for (int k = 0; k < n_points; k++)
        {
          elInfo->coordToWorld(quad->getLambda(k), quad_pts);
          sum += quad->getWeight(k) *
                 (*(*matrix_fcts)[i][j])(quad_pts, *rec_struct->coords);
        }
        (*(rec_struct->A))[i][j] += sum * elInfo->getDet();
      }

      Quadrature* quad = Quadrature::
                         provideQuadrature(Global::getGeo(WORLD),
                                           (*matrix_fcts)[0][i]->getDegree() + deg_f);
      int n_points = quad->getNumPoints();
      mtl::dense_vector<double> uhAtQP(n_points);

      // Computing contributions to right hand side.
      if (gradient)      // For gradient recovery.
      {
        double fAtQP = 1.0;
        if (f_scal)
        {
          if (aux_vec)
            aux_vec->getVecAtQPs(elInfo, quad, NULL, uhAtQP);
          else
            uh->getVecAtQPs(elInfo, quad, NULL, uhAtQP);
        }

        // Get gradient at quadrature points
        mtl::dense_vector<WorldVector<double>> grdAtQP(n_points);
        uh->getGrdAtQPs(elInfo, quad, NULL, grdAtQP);
        vec_sum = 0.0;
        for (int k = 0; k < n_points; k++)
        {
          elInfo->coordToWorld(quad->getLambda(k), quad_pts);
          if (f_vec)
            fAtQP = (*f_vec)(quad_pts);
          if (f_scal)
            fAtQP = (*f_scal)(uhAtQP[k]);

          vec_sum = vec_sum + grdAtQP[k] * fAtQP * quad->getWeight(k)
                    * (*(*matrix_fcts)[0][i])(quad_pts, *rec_struct->coords);
        }
        (*rec_struct->rec_grdUh)[i] = (*rec_struct->rec_grdUh)[i]
                                      + vec_sum * elInfo->getDet();
      }
      else             // For recovery of DOFVector.
      {
        // Get uh at quadrature points
        uh->getVecAtQPs(elInfo, quad, NULL, uhAtQP);
        double sum = 0.0;
        for (int k = 0; k < n_points; k++)
        {
          elInfo->coordToWorld(quad->getLambda(k), quad_pts);
          sum += uhAtQP[k] * quad->getWeight(k)
                 * (*(*matrix_fcts)[0][i])(quad_pts, *rec_struct->coords);
        }
        (*rec_struct->rec_uh)[i] += sum * elInfo->getDet();
      }
    }
  }


  void Recovery::compute_interior_sums(DOFVector<double>* uh, ElInfo* elInfo,
                                       RecoveryStructure* rec_struct, Quadrature* quad,
                                       AbstractFunction<double, WorldVector<double>>* f_vec,
                                       AbstractFunction<double, double>* f_scal,
                                       DOFVector<double>* aux_vec)
  {
    FUNCNAME_DBG("Recovery::compute_sums()");

    TEST_EXIT_DBG(gradient)("SPR of solution need computing node sums.\n");
    TEST_EXIT_DBG(!(f_vec && f_scal))("Only one diffusion function, please!\n");

    WorldVector<double> vec_sum;
    int n_points = quad->getNumPoints();
    WorldVector<double> quad_pts;  // For world coordinates of quadrature points.
    mtl::dense_vector<double> uhAtQP(n_points);
    mtl::dense_vector<WorldVector<double>> grdAtQP(n_points);

    for (int i = 0; i < n_monomials; i++)
    {
      // Computing contributions to system matrix.
      for (int j = i; j < n_monomials; j++)
      {
        double sum = 0.0;
        for (int k = 0; k < n_points; k++)
        {
          elInfo->coordToWorld(quad->getLambda(k), quad_pts);
          sum += (*(*matrix_fcts)[i][j])(quad_pts, *rec_struct->coords);
        }
        (*(rec_struct->A))[i][j] += sum;
      }

      // Computing contributions to right hand side.
      double fAtQP = 1.0;
      if (f_scal)
      {
        if (aux_vec)
          aux_vec->getVecAtQPs(elInfo, quad, NULL, uhAtQP);
        else
          uh->getVecAtQPs(elInfo, quad, NULL, uhAtQP);
      }

      // Get gradient at quadrature points
      uh->getGrdAtQPs(elInfo, quad, NULL, grdAtQP);
      vec_sum = 0.0;
      for (int k = 0; k < n_points; k++)
      {
        elInfo->coordToWorld(quad->getLambda(k), quad_pts);
        if (f_vec)
          fAtQP = (*f_vec)(quad_pts);
        if (f_scal)
          fAtQP = (*f_scal)(uhAtQP[k]);

        vec_sum = vec_sum + grdAtQP[k] * fAtQP
                  * (*(*matrix_fcts)[0][i])(quad_pts, *rec_struct->coords);
      }
      (*rec_struct->rec_grdUh)[i] = (*rec_struct->rec_grdUh)[i] + vec_sum;
    }
  }


  void Recovery::compute_node_sums(DOFVector<double>* uh, ElInfo* elInfo,
                                   RecoveryStructure* rec_struct, DimVec<int> preDofs,
                                   int n_vertices, int n_edges, int n_faces)
  {
    FUNCNAME_DBG("Recovery::compute_sums()");

    TEST_EXIT_DBG(!gradient)
    ("SPR of flux or gradient need computing interior sums\n");
    TEST_EXIT_DBG(feSpace->getMesh()->getDim() == 1)
    ("At the moment only for linear finite elements.\n");

    WorldVector<double> node;  // For world coordinates at nodes.
    const DegreeOfFreedom** dof = elInfo->getElement()->getDof();
    DenseVector<double> uh_loc(n_vertices);
    uh->getLocalVector(elInfo->getElement(), uh_loc);

    for (int l = 0; l < n_vertices; l++)
    {
      // Computing contributions of vertex nodes
      if (rec_struct->neighbors->insert(dof[l][preDofs[VERTEX]]).second)
      {
        node = elInfo->getCoord(l);
        for (int i = 0; i < n_monomials; i++)
        {
          // Computing contributions to system matrix.
          for (int j = i; j < n_monomials; j++)
            (*(rec_struct->A))[i][j] += (*(*matrix_fcts)[i][j])(node,
                                        *rec_struct->coords);

          // Computing contributions to right hand side.
          (*rec_struct->rec_uh)[i] += uh_loc[l]
                                      * (*(*matrix_fcts)[0][i])(node, *rec_struct->coords);
        }
      }
    }
  }


  void Recovery::compute_sums_linear(DOFVector<double>* uh, ElInfo* elInfo,
                                     RecoveryStructure* rec_struct,
                                     int vertex, DimVec<int> preDofs,
                                     int n_vertices)
  {
    FUNCNAME_DBG("Recovery::compute_sums_linear()");

    TEST_EXIT_DBG(!gradient)
    ("SPR of flux or gradient need computing interior sums\n");

    WorldVector<double> node;     // For world coordinates at nodes.
    const DegreeOfFreedom** dof = elInfo->getElement()->getDof();
    DenseVector<double> uh_loc(n_vertices);
    uh->getLocalVector(elInfo->getElement(), uh_loc);

    for (int l = 0;  l < n_vertices; l++)
    {
      // Computing contributions of vertex nodes
      DegreeOfFreedom k = dof[l][preDofs[VERTEX]];
      if (rec_struct->neighbors->insert(k).second)
      {
        node = elInfo->getCoord(l);
        for (int i = 0; i < n_monomials; i++)
        {
          // Computing contributions to system matrix.
          for (int j = i; j < n_monomials; j++)
            (*(rec_struct->A))[i][j] += (*(*matrix_fcts)[i][j])(node,
                                        *rec_struct->coords);

          // Computing contributions to right hand side.
          (*rec_struct->rec_uh)[i] += uh_loc[l]
                                      * (*(*matrix_fcts)[0][i])(node, *rec_struct->coords);
        }
      }
    }

    if (vertex > 1 && elInfo->getNeighbour(vertex))
    {
      int oppVertex = elInfo->getOppVertex(vertex);
      DegreeOfFreedom k = elInfo->getNeighbour(vertex)->getDof(oppVertex)[preDofs[VERTEX]];

      if (rec_struct->neighbors->insert(k).second)
      {
        node = elInfo->getOppCoord(vertex);
        for (int i = 0; i < n_monomials; i++)
        {
          // Computing contributions to system matrix.
          for (int j = i; j < n_monomials; j++)
            (*(rec_struct->A))[i][j] += (*(*matrix_fcts)[i][j])(node, *rec_struct->coords);

          // Computing contributions to right hand side.
          (*rec_struct->rec_uh)[i] += (*uh)[k] * (*(*matrix_fcts)[0][i])(node, *rec_struct->coords);
        }
      }
    }
  }


  void Recovery::fill_struct_vec(DOFVector<double>* uh,
                                 AbstractFunction<double, WorldVector<double>>* f_vec,
                                 AbstractFunction<double, double>* f_scal,
                                 DOFVector<double>* aux_vec)
  {
    FUNCNAME_DBG("Recovery::fill_struct_vec()");

    // Information on the mesh.
    Mesh* mesh = feSpace->getMesh();
    int dim = mesh->getDim();

    // Geometric information.
    int n_vertices = Global::getGeo(VERTEX, dim);
    int n_edges = Global::getGeo(EDGE, dim);
    int n_faces = Global::getGeo(FACE, dim);

    // Information concerning the finite element space.
    const BasisFunction* basis_fcts = feSpace->getBasisFcts();
    DimVec<int>* nDOFs = basis_fcts->getNumberOfDofs();

    // Information from DOFAdmin.
    const DOFAdmin* admin = feSpace->getAdmin();
    DimVec<int> preDofs = admin->getNumberOfPreDofs();

    // Variables for storing temporary information.
    DimVec<DegreeOfFreedom> interior_vertices(dim);
    WorldVector<double> coordinates;

    // Variables for passing information to integration routines.
    int degree = basis_fcts->getDegree();
    Quadrature* quad = NULL;
    if (gradient && !method)
      quad = Quadrature::provideQuadrature(Global::getGeo(WORLD), degree);

    DimVec<int> pre_dofs(dim, NO_INIT);
    if (!gradient)
      pre_dofs = uh->getFeSpace()->getAdmin()->getNumberOfPreDofs();

    // Variables for traversing the mesh.
    Flag fill_flag =
      Mesh::CALL_LEAF_EL |
      Mesh::FILL_COORDS |
      Mesh::FILL_BOUND |
      Mesh::FILL_DET |
      Mesh::FILL_GRD_LAMBDA;

    if (degree == 2 && dim > 1)
      fill_flag |= Mesh::FILL_NEIGH | Mesh::FILL_OPP_COORDS;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    DofContainer boundaryDofs;
    DofContainerSet boundaryDofSet;
    Parallel::MeshDistributor::globalMeshDistributor->getAllBoundaryDofs(feSpace, 0, boundaryDofs);
    boundaryDofSet.insert(boundaryDofs.begin(), boundaryDofs.end());
#endif

    TraverseStack stack;
    ElInfo* el_info = stack.traverseFirst(mesh, -1, fill_flag);

    while (el_info)      // traversing the mesh.
    {
      const DegreeOfFreedom** dof = el_info->getElement()->getDof();

      int n_neighbors = 0;     // counting interior vertices of element
      for (int i = 0; i < n_vertices; i++)
      {
        DegreeOfFreedom k = dof[i][preDofs[VERTEX]];
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        bool isInterior = (el_info->getBoundary(VERTEX, i) == INTERIOR) && (boundaryDofSet.count(dof[i]) == 0);
#else
        bool isInterior = el_info->getBoundary(VERTEX, i) == INTERIOR;
#endif

        if (isInterior)
          interior_vertices[n_neighbors++] = k;
      }

      TEST_EXIT_DBG(n_neighbors)
      ("Each element should have a least one interior vertex!\n");

      for (int i = 0; i < n_vertices; i++)       // Handling nodes on vertices.
      {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        bool isInterior = (el_info->getBoundary(VERTEX, i) == INTERIOR) && (boundaryDofSet.count(dof[i]) == 0);
#else
        bool isInterior = el_info->getBoundary(VERTEX, i) == INTERIOR;
#endif

        DegreeOfFreedom k = dof[i][preDofs[VERTEX]];

        // Setting world coordinates of node.
        if (!(*struct_vec)[k].coords)
        {
          (*struct_vec)[k].coords  = new WorldVector<double>;
          *(*struct_vec)[k].coords = el_info->getCoord(i);
        }

        if (isInterior)
        {
          // Allocating memory for matrix and right hand side.
          if (!(*struct_vec)[k].A)
          {
            (*struct_vec)[k].A = new Matrix<double>(n_monomials, n_monomials);
            *(*struct_vec)[k].A = 0.0;

            if (gradient)
            {
              (*struct_vec)[k].rec_grdUh = new Vector<WorldVector<double>>(n_monomials);
              for (int j = 0; j < n_monomials; j++)
                (*(*struct_vec)[k].rec_grdUh)[j] = 0.0;
            }
            else
            {
              (*struct_vec)[k].rec_uh = new Vector<double>(n_monomials);
              *(*struct_vec)[k].rec_uh = 0.0;
            }
          }

          // Computing the integrals.
          if (method)
            compute_integrals(uh, el_info, &(*struct_vec)[k],
                              f_vec, f_scal, aux_vec);
          else if (gradient)
            compute_interior_sums(uh, el_info, &(*struct_vec)[k], quad,
                                  f_vec, f_scal, aux_vec);
          else
          {
            if (!(*struct_vec)[k].neighbors)
              (*struct_vec)[k].neighbors = new std::set<DegreeOfFreedom>;

            if (degree == 2 && dim > 1)
              compute_sums_linear(uh, el_info, &(*struct_vec)[k],
                                  i, pre_dofs, n_vertices);
            else
              compute_node_sums(uh, el_info, &(*struct_vec)[k],
                                pre_dofs, n_vertices, n_edges, n_faces);
          }
        }
        else         // Setting list of adjacent interior vertices.
        {
          if (!(*struct_vec)[k].neighbors)
            (*struct_vec)[k].neighbors = new std::set<DegreeOfFreedom>;

          for (int j = 0; j < n_neighbors; j++)
            (*struct_vec)[k].neighbors->insert(interior_vertices[j]);
        }
      }

      int n_count = n_vertices;

      if (dim > 1)       // Handling nodes on edges.
      {
        for (int i = 0; i < n_edges; i++)
        {
          bool isEdgeInterior = el_info->getBoundary(EDGE, i) == INTERIOR;

          for (int j = 0; j < (*nDOFs)[EDGE]; j++)
          {
            DegreeOfFreedom k = dof[n_vertices + i][preDofs[EDGE] + j];

            if (!(*struct_vec)[k].coords)
            {
              // Setting world coordinates of node.
              el_info->coordToWorld(*basis_fcts->getCoords(n_count),
                                    coordinates);
              (*struct_vec)[k].coords = new WorldVector<double>;
              *(*struct_vec)[k].coords = coordinates;

              // Setting list of adjacent interior vertices.
              (*struct_vec)[k].neighbors = new std::set<DegreeOfFreedom>;

              if (isEdgeInterior)
              {
                for (int m = 0; m < 2; m++)
                {
                  int l = Global::getReferenceElement(dim)->getVertexOfEdge(i, m);
                  if (el_info->getBoundary(VERTEX, l) == INTERIOR)
                    (*struct_vec)[k].neighbors->insert(dof[l][preDofs[VERTEX]]);
                }
              }
              else
              {
                for (int m = 0; m < n_neighbors; m++)
                  (*struct_vec)[k].neighbors->insert(interior_vertices[m]);
              }
            }

            n_count++;
          }
        }
      }

      if (dim == 3)     // Handling nodes on faces.
        for (int i = 0; i < n_faces; i++)
          for (int j = 0; j < (*nDOFs)[FACE]; j++)
          {
            DegreeOfFreedom k = dof[n_vertices+n_edges + i][preDofs[FACE] + j];

            if (!(*struct_vec)[k].coords)
            {
              // Setting world coordinates of node.
              el_info->coordToWorld(*basis_fcts->getCoords(n_count),
                                    coordinates);
              (*struct_vec)[k].coords  = new WorldVector<double>;
              *(*struct_vec)[k].coords = coordinates;

              // Setting list of adjacent interior vertices.
              (*struct_vec)[k].neighbors = new std::set<DegreeOfFreedom>;

              if (el_info->getBoundary(FACE, i) == INTERIOR)
                for (int m = 0; m < 3; m++)
                {
                  int l = Global::getReferenceElement(dim)->getVertexOfPosition(FACE, i, m);
                  if (el_info->getBoundary(VERTEX, l) == INTERIOR)
                    (*struct_vec)[k].neighbors->insert(dof[l][preDofs[VERTEX]]);
                }
              else
                for (int m = 0; m < n_neighbors; m++)
                  (*struct_vec)[k].neighbors->insert(interior_vertices[m]);
            }

            n_count++;
          }

      if ((*nDOFs)[CENTER])    // Handling nodes on center of element.
        for (int j = 0; j < (*nDOFs)[CENTER]; j++)
        {
          DegreeOfFreedom k = dof[n_vertices+n_edges+n_faces][preDofs[CENTER] + j];

          // Setting world coordinates of node.
          el_info->coordToWorld(*basis_fcts->getCoords(n_count), coordinates);
          (*struct_vec)[k].coords = new WorldVector<double>;
          *(*struct_vec)[k].coords = coordinates;

          // Setting list of adjacent interior vertices.
          (*struct_vec)[k].neighbors = new std::set<DegreeOfFreedom>;

          for (int m = 0; m < n_neighbors; m++)
            (*struct_vec)[k].neighbors->insert(interior_vertices[m]);

          n_count++;
        }

      el_info = stack.traverseNext(el_info);
    }
  }


  void Recovery::recoveryUh(DOFVector<double>* uh, DOFVector<double>& rec_vec)
  {
    FUNCNAME("Recovery::recoveryUh()");

    clear();

    gradient = false;

    const FiniteElemSpace* fe_space = rec_vec.getFeSpace();
    set_feSpace(fe_space);                                  // Setting feSpace.
    set_exponents(feSpace->getBasisFcts()->getDegree());    // Setting exponents.
    fill_struct_vec(uh);               // Filling vector of recovery structures.

    DOFVector<RecoveryStructure>::Iterator SV_it(struct_vec, USED_DOFS);

    // Solving local systems.
    for (SV_it.reset(); !SV_it.end(); ++SV_it)
    {
      if ((*SV_it).A)
      {
        TEST(Cholesky::solve((*SV_it).A, (*SV_it).rec_uh, (*SV_it).rec_uh))
        ("There must be an error, matrix is not positive definite.\n");
      }
    }

    // define result vector
    DOFVector<double>* result = &rec_vec;

    result->set(0.0);

    DOFVector<double>::Iterator result_it(result, USED_DOFS);
    std::set<DegreeOfFreedom>::const_iterator setIterator;

    // #ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    //   DOFVector<int> neighbourSize(feSpace, "neighbourSize");
    //   neighbourSize.set(0);
    // #endif

    for (SV_it.reset(), result_it.reset(); !result_it.end(); ++SV_it, ++result_it)
    {
      if ((*SV_it).rec_uh)
      {
        *result_it = (*(*SV_it).rec_uh)[0];
      }
      else
      {
        if ((*SV_it).neighbors)
        {
          for (setIterator  = (*SV_it).neighbors->begin();
               setIterator != (*SV_it).neighbors->end();
               ++setIterator)
          {
            for (int i = 0; i < n_monomials; i++)
              *result_it += (*(*struct_vec)[*setIterator].rec_uh)[i] *
                            (*(*matrix_fcts)[0][i])(*(*SV_it).coords, *(*struct_vec)[*setIterator].coords);
          }
          // #ifdef HAVE_PARALLEL_DOMAIN_AMDIS
          // 	neighbourSize[result_it.getDOFIndex()] = (*SV_it).neighbors->size();
          // #else
          *result_it /= (*SV_it).neighbors->size();
          //#endif
        }
        else
        {
          *result_it = 0.0;
        }
      }
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    //   MeshDistributor::globalMeshDistributor->synchAddVector(rec_vec);
    //   MeshDistributor::globalMeshDistributor->synchAddVector(neighbourSize);

    //   for (result_it.reset(); !result_it.end(); ++result_it)
    //     if (neighbourSize[result_it.getDOFIndex()] > 0)
    //       *result_it /= neighbourSize[result_it.getDOFIndex()];

    Parallel::MeshDistributor::globalMeshDistributor->synchVector(rec_vec);
#endif
  }


  DOFVector<double>*
  Recovery::recoveryUh(DOFVector<double>* uh, const FiniteElemSpace* fe_space)
  {
    FUNCNAME("Recovery::recoveryUh()");

    clear();

    gradient = false;

    set_feSpace(fe_space);                                  // Setting feSpace.
    set_exponents(feSpace->getBasisFcts()->getDegree());    // Setting exponents.
    fill_struct_vec(uh);               // Filling vector of recovery structures.

    DOFVector<RecoveryStructure>::Iterator SV_it(struct_vec, USED_DOFS);

    // Solving local systems.
    for (SV_it.reset(); !SV_it.end(); ++SV_it)
    {
      if ((*SV_it).A)
      {
        TEST(Cholesky::solve((*SV_it).A, (*SV_it).rec_uh, (*SV_it).rec_uh))
        ("There must be an error, matrix is not positive definite.\n");
      }
    }

    // define result vector
    static DOFVector<double>* vec = NULL;// TODO: REMOVE STATIC
    DOFVector<double>* result = NULL;

    // Allocate memory for result vector
    if (vec && vec->getFeSpace() != feSpace)
    {
      delete vec;
      vec = NULL;
    }

    if (!vec)
      vec = new DOFVector<double>(feSpace, "gradient");
    result = vec;

    result->set(0.0);

    DOFVector<double>::Iterator result_it(result, USED_DOFS);
    std::set<DegreeOfFreedom>::const_iterator setIterator;

    for (SV_it.reset(), result_it.reset(); !result_it.end();
         ++SV_it, ++result_it)
    {
      if ((*SV_it).rec_uh)
      {
        *result_it = (*(*SV_it).rec_uh)[0];
      }
      else
      {
        if ((*SV_it).neighbors)
        {
          for (setIterator  = (*SV_it).neighbors->begin();
               setIterator != (*SV_it).neighbors->end();
               ++setIterator)
          {
            for (int i = 0; i < n_monomials; i++)
              *result_it = *result_it + (*(*struct_vec)[*setIterator].rec_uh)[i] *
                           (*(*matrix_fcts)[0][i])(*(*SV_it).coords,
                                                   *(*struct_vec)[*setIterator].coords);
          }
          *result_it /= (*SV_it).neighbors->size();
        }
        else
          *result_it=0.0;
      }
    }

    return result;
  }


  DOFVector<WorldVector<double>>*
                              Recovery::recovery(DOFVector<double>* uh, const FiniteElemSpace* fe_space,
                                  AbstractFunction<double, WorldVector<double>>* f_vec,
                                  AbstractFunction<double, double>* f_scal,
                                  DOFVector<double>* aux_vec)
  {
    FUNCNAME_DBG("Recovery::recovery()");

    clear();

    gradient = true;

    set_feSpace(fe_space);                                  // Setting feSpace.
    set_exponents(feSpace->getBasisFcts()->getDegree());    // Setting exponents.
    fill_struct_vec(uh, f_vec, f_scal, aux_vec);  // Filling vec. of rec. struct.

    DOFVector<RecoveryStructure>::Iterator SV_it(struct_vec, USED_DOFS);

    // Solving local systems.
    for (SV_it.reset(); !SV_it.end(); ++SV_it)
    {
      if ((*SV_it).A)
      {
        DBG_VAR(int error =)
        Cholesky::solve((*SV_it).A,
                        (*SV_it).rec_grdUh,
                        (*SV_it).rec_grdUh);
        TEST_EXIT_DBG(error)
        ("There must be some error, matrix is not positive definite.\n");
      }
    }

    // define result vector
    static DOFVector<WorldVector<double>>* vec = NULL;// TODO: REMOVE STATIC
    DOFVector<WorldVector<double>>* result = NULL;

    // Allocate memory for result vector
    if (vec && vec->getFeSpace() != feSpace)
    {
      delete vec;
      vec = NULL;
    }

    if (!vec)
      vec = new DOFVector<WorldVector<double>>(feSpace, "gradient");

    result = vec;

    result->set(WorldVector<double>(DEFAULT_SIZE, 0.0));

    DOFVector<WorldVector<double>>::Iterator grdIt(result, USED_DOFS);
    std::set<DegreeOfFreedom>::const_iterator setIterator;

    for (SV_it.reset(), grdIt.reset(); !grdIt.end(); ++SV_it, ++grdIt)
    {
      if ((*SV_it).rec_grdUh)
      {
        *grdIt = (*(*SV_it).rec_grdUh)[0];
      }
      else
      {
        for (setIterator  = (*SV_it).neighbors->begin();
             setIterator != (*SV_it).neighbors->end();
             ++setIterator)
        {
          for (int i = 0; i < n_monomials; i++)
            *grdIt = *grdIt + (*(*struct_vec)[*setIterator].rec_grdUh)[i] *
                     (*(*matrix_fcts)[0][i])(*(*SV_it).coords,
                                             *(*struct_vec)[*setIterator].coords);
        }
        *grdIt = *grdIt * (1.0 / (*SV_it).neighbors->size());
      }
    }

    return result;
  }


  DOFVector<WorldVector<double>>*
                              Recovery::recovery(DOFVector<double>* uh,
                                  AbstractFunction<double, WorldVector<double>>* f_vec,
                                  AbstractFunction<double, double>* f_scal,
                                  DOFVector<double>* aux_vec)
  {
    FUNCNAME_DBG("Recovery::simpleAveraging()");

    TEST_EXIT_DBG(!(f_vec && f_scal))("Only one diffusion function, please!\n");

    const FiniteElemSpace* fe_space = uh->getFeSpace();

    // define result vector
    static DOFVector<WorldVector<double>>* vec = NULL;// TODO: REMOVE STATIC
    DOFVector<WorldVector<double>>* result = NULL;

    // Allocate memory for result vector
    if (vec && vec->getFeSpace() != fe_space)
    {
      delete vec;
      vec = NULL;
    }

    if (!vec)
      vec = new DOFVector<WorldVector<double>>(fe_space, "gradient");

    result = vec;
    result->set(WorldVector<double>(DEFAULT_SIZE, 0.0));

    DOFVector<double> volume(fe_space, "volume");
    volume.set(0.0);

    Mesh* mesh = fe_space->getMesh();
    int dim = mesh->getDim();
    const BasisFunction* basFcts = fe_space->getBasisFcts();
    DOFAdmin* admin = fe_space->getAdmin();
    int numPreDofs = admin->getNumberOfPreDofs(0);
    DimVec<double> bary(dim, (1.0 / (dim + 1.0)));
    WorldVector<double> barycenter;     // For world coordinates at barycenter

    // traverse mesh
    TraverseStack stack;
    Flag fillFlag =
      Mesh::CALL_LEAF_EL | Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA | Mesh::FILL_COORDS;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, fillFlag);

    DenseVector<double> localUh(basFcts->getNumber());

    while (elInfo)
    {
      double det = elInfo->getDet();
      const DegreeOfFreedom** dof = elInfo->getElement()->getDof();
      uh->getLocalVector(elInfo->getElement(), localUh);


      const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();
      WorldVector<double> grd;
      basFcts->evalGrdUh(bary, grdLambda, localUh, grd);

      double fAtBary = 1.0;

      if (f_vec)
      {
        elInfo->coordToWorld(bary, barycenter);
        fAtBary = (*f_vec)(barycenter);
      }

      if (f_scal)
      {
        if (aux_vec)
          aux_vec->getLocalVector(elInfo->getElement(), localUh);

        fAtBary = basFcts->evalUh(bary, localUh);
        fAtBary = (*f_scal)(fAtBary);
      }

      for (int i = 0; i < dim + 1; i++)
      {
        DegreeOfFreedom dofIndex = dof[i][numPreDofs];
        (*result)[dofIndex] += grd * fAtBary * det;
        volume[dofIndex] += det;
      }

      elInfo = stack.traverseNext(elInfo);
    }

    DOFVector<double>::Iterator volIt(&volume, USED_DOFS);
    DOFVector<WorldVector<double>>::Iterator grdIt(result, USED_DOFS);

    for (volIt.reset(), grdIt.reset(); !volIt.end(); ++volIt, ++grdIt)
      *grdIt *= 1.0/(*volIt);

    return result;
  }


  void Recovery::test(DOFVector<double>* uh, const FiniteElemSpace* fe_space)
  {
    FUNCNAME("Recovery::test()");

    clear();

    set_feSpace(fe_space);                                  // Setting feSpace.
    set_exponents(feSpace->getBasisFcts()->getDegree());    // Setting exponents.
    fill_struct_vec(uh);                // Filling vector of recovery structures.

    DOFVector<RecoveryStructure>::Iterator LM_iterator(struct_vec, USED_DOFS);
    int position;
    WorldVector<double> coord;

    // for every DOFs
    for (LM_iterator.reset(); !LM_iterator.end(); ++LM_iterator)
    {
      position = LM_iterator.getDOFIndex();
      MSG("Node: ");
      std::cout << position << std::endl;
      (*struct_vec)[position].print();
    }
  }

}
