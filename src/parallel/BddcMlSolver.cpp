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


#ifdef HAVE_BDDC_ML

extern "C" {
#include <bddcml_interface_c.h>
}

#include "parallel/BddcMlSolver.hpp"
#include "parallel/MpiHelper.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    using namespace std;

    int BddcMlSolver::solveSystemImpl(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                                        SystemVector& vec,
                                        SystemVector& rhsVec,
                                        bool /* createMatrixData */,
                                        bool /* storeMatrixData */)
    {
      FUNCNAME("BddcMlSolver::solvePetscMatrix()");

      TEST_EXIT(meshDistributor)("No meshDistributor provided!\n");

      const Matrix<DOFMatrix*>* mat = A.getOriginalMat();

      int nComponents = vec.getSize();
      const FiniteElemSpace* feSpace = vec.getFeSpace(0);
      Mesh* mesh = feSpace->getMesh();


      // === First, create a continous leaf element index on each subdomain ===

      std::set<int> leafElIndex;
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        leafElIndex.insert(elInfo->getElement()->getIndex());
        elInfo = stack.traverseNext(elInfo);
      }

      map<int, int> mapElIndex;
      int nLeafEls = 0;
      for (std::set<int>::iterator it = leafElIndex.begin();
           it != leafElIndex.end(); ++it)
        mapElIndex[*it] = nLeafEls++;



      int nLevel = 2;
      int nSubdomains[nLevel];
      nSubdomains[0] = meshDistributor->getMpiComm(0).Get_size();
      nSubdomains[1] = 1;

      int nSubPerProc = 1;
      MPI_Fint c2f = MPI_Comm_c2f(meshDistributor->getMpiComm(0));
      int verboseLevel = 2;
      int numbase = 0;

      //     MSG("call to \"bddcml_init\" with the following arguments:\n");
      //     MSG("  %d, [%d, %d], %d, %d, %d, %d, %d\n", nLevel, nSubdomains[0], nSubdomains[1], nLevel, nSubPerProc, c2f, verboseLevel, numbase);

      bddcml_init(&nLevel, nSubdomains, &nLevel, &nSubPerProc,
                  &c2f, &verboseLevel, &numbase);

      // global number of elements
      int nelem = mesh->getNumberOfLeaves();
      mpi::globalAdd(nelem);

      MSG("nelem = %d\n", nelem);

      // global number of nodes
      int nnod = dofMap[feSpace].nOverallDofs;

      MSG("nnod = %d\n", nnod);

      // global number of dofs
      int ndof = nnod * nComponents;

      MSG("ndof = %d\n", ndof);

      // space dimenstion
      int ndim = mesh->getDim();

      // mesh dimension
      int meshdim = mesh->getDim();

      // global indes of subdomain
      int isub = meshDistributor->getMpiRank();

      // local number of elements
      int nelems = nLeafEls;

      MSG("nelems = %d\n", nelems);

      // local number of nodes
      int nnods = feSpace->getAdmin()->getUsedSize();

      // local number of dofs
      int ndofs = nnods * nComponents;

      MSG("local nnods %d     ndofs %d\n", nnods, ndofs);

      int nVertices = mesh->getGeo(VERTEX);

      // Length of array inet
      int linet = nelems * nVertices;

      // Local array with indices of nodes on each element
      int inet[linet];
      elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        TEST_EXIT_DBG(mapElIndex.count(elInfo->getElement()->getIndex()))
        ("Should not happen!\n");

        int localElIndex = mapElIndex[elInfo->getElement()->getIndex()];
        for (int i = 0; i < nVertices; i++)
          inet[localElIndex * nVertices + i] = elInfo->getElement()->getDof(i, 0);
        elInfo = stack.traverseNext(elInfo);
      }


      // local array of number of nodes per element
      int nnet[nelems];
      for (int i = 0; i < nelems; i++)
        nnet[i] = nVertices;

      // local array with number of DOFs per node.
      int nndf[nnods];
      for (int i = 0; i < nnods; i++)
        nndf[i] = nComponents;

      // array of indices of subdomain nodes in global numbering
      int isngn[nnods];
      for (int i = 0; i < nnods; i++)
        isngn[i] = dofMap[feSpace][i].global; //meshDistributor->mapDofToGlobal(feSpace, i);

      // array of indices of subdomain variables in global numbering
      int isvgvn[ndofs];
      for (int j = 0; j < nnods; j++)
        for (int i = 0; i < nComponents; i++)
          isvgvn[j * nComponents + i] =
            dofMap[feSpace][j].global * nComponents + i;

      // array of indices of subdomain elements in global numbering
      int isegn[nelems];
      int rStartEl, nOverallEl;
      mpi::getDofNumbering(meshDistributor->getMpiComm(0),
                           nelems, rStartEl, nOverallEl);
      MSG("rStartEl = %d\n", rStartEl);
      for (int i = 0; i < nelems; i++)
        isegn[i] = rStartEl + i;



      int lxyz1 = nnods;
      int lxyz2 = mesh->getDim();
      // local array with coordinates of nodes
      double xyz[lxyz1 * lxyz2];

      {
        DOFVector<WorldVector<double>> coordDofs(feSpace, "tmp");
        mesh->getDofIndexCoords(coordDofs);

        for (int i = 0; i < lxyz2; i++)
          for (int j = 0; j < nnods; j++)
            xyz[i * nnods + j] = coordDofs[j][i];
      }


      // === Fill for dirichlet boundary conditions. ===

      // local array of indices denoting dirichlet boundary data
      int ifix[ndofs];
      for (int i = 0; i < ndofs; i++)
        ifix[i] = 0;

      // local array of values for dirichlet boundary data
      double fixv[ndofs];
      for (int i = 0; i < ndofs; i++)
        fixv[i] = 0.0;

      for (int i = 0; i < rhsVec.getSize(); i++)
      {
        map<DegreeOfFreedom, double>& dcValues =
          rhsVec.getDOFVector(i)->getDirichletValues();

        for (map<DegreeOfFreedom, double>::iterator it = dcValues.begin();
             it != dcValues.end(); ++it)
        {
          int index = it->first * nComponents + i;
          TEST_EXIT_DBG(index < ndofs)("Should not happen!\n");
          ifix[index] = 1;
          fixv[index] = it->second;
        }
      }


      // local rhs data
      double rhs[ndofs];
      for (int i = 0; i < nComponents; i++)
      {
        DOFVector<double>& dofvec = *(rhsVec.getDOFVector(i));
        for (int j = 0; j < nnods; j++)
          rhs[j * nComponents + i] = dofvec[j];
      }

      // Completenes of the rhs vector on subdomains
      int is_rhs_complete = 0;

      // Local array with initial solution guess
      double sol[ndofs];
      for (int i = 0; i < ndofs; i++)
        sol[i] = 0.0;

      // matrix type (set here to unsymmetric)
      int matrixtype = 0;
      Parameters::get(name + "->bddcml->matrix type", matrixtype);

      // Non zero structure of matrix
      vector<int> i_sparse;
      vector<int> j_sparse;
      vector<double> a_sparse;

      for (int i = 0; i < nComponents; i++)
        for (int j = 0; j < nComponents; j++)
          if ((*mat)[i][j])
            addDofMatrix((*mat)[i][j],
                         i_sparse, j_sparse, a_sparse, nComponents, i, j);


      // Number of non-zero entries in matrix
      int la = i_sparse.size();

      //    MSG("LOCAL LA = %d\n", la);

      // Matrix is assembled
      int is_assembled_int = 0;


      // Users constraints
      double user_constraints = 0.0;

      int luser_constraints1 = 0;

      int luser_constraints2 = 0;

      /*
      string tmp="";

      MSG("call to \"bddcml_upload_subdomain_data\" with the following arguments (each in one line):\n");
      MSG("  nelem = %d\n", nelem);
      MSG("  nnod = %d\n", nnod);
      MSG("  ndof = %d\n", ndof);
      MSG("  ndim = %d\n", ndim);
      MSG("  meshdim = %d\n", meshdim);
      MSG("  isub = %d\n", isub);
      MSG("  nelems = %d\n", nelems);
      MSG("  nnods = %d\n", nnods);
      MSG("  ndofs = %d\n", ndofs);
      MSG("  inet = [%d, %d, %d, %d, %d, %d]\n", inet[0], inet[1], inet[2], inet[3], inet[4], inet[5]);
      MSG("  linet = %d\n", linet);
      MSG("  nnet = [%d, %d]\n", nnet[0], nnet[1]);
      MSG("  lnnet = %d\n", nelems);
      MSG("  nndf = [%d, %d, %d, %d]\n", nndf[0], nndf[1], nndf[2], nndf[3]);
      MSG("  lnndf = %d\n", nnods);
      MSG("  isngn = [%d, %d, %d, %d]\n", isngn[0], isngn[1], isngn[2], isngn[3]);
      MSG("  lisngn = %d\n", nnods);
      MSG("  isvgvn = [%d, %d, %d, %d]\n", isvgvn[0], isvgvn[1], isvgvn[2], isvgvn[3]);
      MSG("  lisvgvn = %d\n", nnods);
      MSG("  isegn = [%d, %d]\n", isegn[0], isegn[1]);
      MSG("  lisegn = %d\n", nelems);
      MSG("  xyz = [%f, %f, %f, %f, %f, %f, %f, %f]\n", xyz[0], xyz[1], xyz[2], xyz[3], xyz[4], xyz[5], xyz[6], xyz[7]);
      MSG("  lxyz1 = %d\n", lxyz1);
      MSG("  lxyz2 = %d\n", lxyz2);
      MSG("  ifix = [%d, %d, %d, %d]\n", ifix[0], ifix[1], ifix[2], ifix[3]);
      MSG("  lifix = %d\n", ndofs);
      MSG("  fixv = [%f, %f, %f, %f]\n", fixv[0], fixv[1], fixv[2], fixv[3]);
      MSG("  lfixv = %d\n", ndofs);
      MSG("  rhs = [%f, %f, %f, %f]\n", rhs[0], rhs[1], rhs[2], rhs[3]);
      MSG("  lrhs = %d\n", ndofs);
      MSG("  is_rhs_complete = %d\n", is_rhs_complete);
      MSG("  sol = [%f, %f, %f, %f]\n", sol[0], sol[1], sol[2], sol[3]);
      MSG("  lsol = %d\n", ndofs);
      MSG("  matrixtype = %d\n", matrixtype);


      MSG("  ispare = [%d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", i_sparse[0], i_sparse[1], i_sparse[2], i_sparse[3], i_sparse[4], i_sparse[5], i_sparse[6], i_sparse[7], i_sparse[8], i_sparse[9]);
      MSG("  jspare = [%d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", j_sparse[0], j_sparse[1], j_sparse[2], j_sparse[3], j_sparse[4], j_sparse[5], j_sparse[6], j_sparse[7], j_sparse[8], j_sparse[9]);
      MSG("  a_spare = [%f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", a_sparse[0], a_sparse[1], a_sparse[2], a_sparse[3], a_sparse[4], a_sparse[5], a_sparse[6], a_sparse[7], a_sparse[8], a_sparse[9]);
      MSG("  la = %d\n", la);
      MSG("  is_assembled_int = %d\n", is_assembled_int);
      */

      bddcml_upload_subdomain_data(&nelem,
                                   &nnod,
                                   &ndof,
                                   &ndim,
                                   &meshdim,
                                   &isub,
                                   &nelems,
                                   &nnods,
                                   &ndofs,
                                   inet,
                                   &linet,
                                   nnet,
                                   &nelems,
                                   nndf,
                                   &nnods,
                                   isngn,
                                   &nnods,
                                   isvgvn,
                                   &ndofs,
                                   isegn,
                                   &nelems,
                                   xyz,
                                   &lxyz1,
                                   &lxyz2,
                                   ifix,
                                   &ndofs,
                                   fixv,
                                   &ndofs,
                                   rhs,
                                   &ndofs,
                                   &is_rhs_complete,
                                   sol,
                                   &ndofs,
                                   &matrixtype,
                                   &(i_sparse[0]),
                                   &(j_sparse[0]),
                                   &(a_sparse[0]),
                                   &la,
                                   &is_assembled_int,
                                   &user_constraints,
                                   &luser_constraints1,
                                   &luser_constraints2);


      int use_defaults_int = 0;
      int parallel_division_int = 1;
      int use_arithmetic_int = 1;
      int use_adaptive_int = 0;
      int use_user_constraint_int = 0;

      Parameters::get(name + "->bddcml->arithmetic constraints", use_arithmetic_int);

      MSG("call to \"bddcml_setup_preconditioner\" with the following arguments (each in one line):\n");
      MSG("  matrixtype              = %d\n", matrixtype);
      MSG("  use_defaults_int        = %d\n", use_defaults_int);
      MSG("  parallel_division_int   = %d\n", parallel_division_int);
      MSG("  use_arithmetic_int      = %d\n", use_arithmetic_int);
      MSG("  use_adaptive_int        = %d\n", use_adaptive_int);
      MSG("  use_user_constraint_int = %d\n", use_user_constraint_int);

      bddcml_setup_preconditioner(&matrixtype,
                                  &use_defaults_int,
                                  &parallel_division_int,
                                  &use_arithmetic_int,
                                  &use_adaptive_int,
                                  &use_user_constraint_int);

      int method = 1;
      double tol = 1.e-6;
      int maxit = 1000;
      int ndecrmax = 30;
      int num_iter = 0;
      int converged_reason = 0;
      double condition_number = 0.0;

      bddcml_solve(&c2f,
                   &method,
                   &tol,
                   &maxit,
                   &ndecrmax,
                   &num_iter,
                   &converged_reason,
                   &condition_number);

      MSG("BDDCML converged reason: %d within %d iterations \n",
          converged_reason, num_iter);

      bddcml_download_local_solution(&isub, rhs, &ndofs);

      for (int i = 0; i < nComponents; i++)
      {
        DOFVector<double>& dofvec = *(vec.getDOFVector(i));
        for (int j = 0; j < nnods; j++)
          dofvec[j] = rhs[j * nComponents + i];
      }

      bddcml_finalize();

      return 0;
    }


    void BddcMlSolver::addDofMatrix(const DOFMatrix* dmat,
                                    vector<int>& i_sparse,
                                    vector<int>& j_sparse,
                                    vector<double>& a_sparse,
                                    int nComponents,
                                    int ithRowComponent,
                                    int ithColComponent)
    {
      FUNCNAME("BddcMlSolver::addDofMatrix()");

      TEST_EXIT(dmat)("Should not happen!\n");

      const FiniteElemSpace* feSpace = dmat->getFeSpace();

      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;

      traits::col<Matrix>::type col(dmat->getBaseMatrix());
      traits::const_value<Matrix>::type value(dmat->getBaseMatrix());

      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      for (cursor_type cursor = begin<row>(dmat->getBaseMatrix()),
           cend = end<row>(dmat->getBaseMatrix()); cursor != cend; ++cursor)
      {
        for (icursor_type icursor = begin<nz>(cursor), icend = end<nz>(cursor);
             icursor != icend; ++icursor)
        {
          i_sparse.push_back(cursor.value() * nComponents + ithRowComponent);
          j_sparse.push_back(col(*icursor) * nComponents + ithColComponent);
          a_sparse.push_back(value(*icursor));
        }
      }

    }

  }
} // end namespaces

#endif

