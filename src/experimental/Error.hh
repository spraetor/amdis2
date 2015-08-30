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


#include "Mesh.h"
#include "Parametric.h"
#include "Quadrature.h"
#include "Traverse.h"

namespace AMDiS
{

  template<typename T>
  T Error<T>::errUFct(const DimVec<double>& lambda)
  {
    WorldVector<double> x;
    elinfo->coordToWorld(lambda, x);
    return ((*pU)(x));
  }


  template<typename T>
  WorldVector<T> Error<T>::grdErrUFct(const DimVec<double>& lambda)
  {
    WorldVector<double> x;
    elinfo->coordToWorld(lambda, x);
    return ((*pGrdU)(x));
  }


  template<typename T>
  double Error<T>::maxErrAtQp(const AbstractFunction<T, WorldVector<double>>& u,
                              const DOFVector<T>& uh,
                              const Quadrature* q)
  {
    FUNCNAME("Error<T>::maxErrAtQp()");

    const FiniteElemSpace* fe_space;
    if (!(pU = &u))
    {
      ERROR("no function u specified; doing nothing\n");
      return(-1.0);
    }
    if (!(errUh = &uh) || !(fe_space = uh->getFeSpace()))
    {
      ERROR("no discrete function or no fe_space for it; doing nothing\n");
      return(-1.0);
    }
    if (!(basFct = fe_space->getBasisFcts()))
    {
      ERROR("no basis functions at discrete solution ; doing nothing\n");
      return(-1.0);
    }

    if (!q)
      q = Quadrature::provideQuadrature(fe_space->getMesh()->getDim(),
                                        2 * fe_space->getBasisFcts()->getDegree() -  2);

    quadFast = FastQuadrature::provideFastQuadrature(basFct, *q, INIT_PHI);
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<double> uh_vec(nPoints);
    double maxErr = 0.0;
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(fe_space->getMesh(), -1,
                                         Mesh::FILL_COORDS | Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      elinfo = elInfo;
      double err = 0.0;
      const double* u_vec;

      u_vec = quadFast->getQuadrature()->fAtQp(errU, NULL);
      errUh->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);

      for (int i = 0; i < nPoints; i++)
      {
        err = u_vec[i] > uh_vec[i] ? u_vec[i] - uh_vec[i] : uh_vec[i] - u_vec[i];
        maxErr = std::max(maxErr, err);
      }

      elInfo = stack.traverseNext(elInfo);
    }

    return maxErr;
  }


  template<typename T>
  double Error<T>::H1Err(const AbstractFunction<WorldVector<T>, WorldVector<double>>& grdU,
                         const DOFVector<T>& uh,
                         int relErr,
                         double* max,
                         bool writeLeafData,
                         int comp)
  {
    FUNCNAME("Error<T>::H1Err()");

    const FiniteElemSpace* fe_space;
    writeInLeafData = writeLeafData;
    component = comp;
    Quadrature* q = NULL;
    pGrdU = &grdU;
    errUh = &uh;

    if (!(fe_space = uh.getFeSpace()))
    {
      ERROR("no fe_space for uh; doing nothing\n");
      return(0.0);
    }
    if (!(basFct = fe_space->getBasisFcts()))
    {
      ERROR("no basis functions at discrete solution ; doing nothing\n");
      return(0.0);
    }

    int dim = fe_space->getMesh()->getDim();
    int deg = grdU.getDegree();
    int degree = deg ? deg : 2 * fe_space->getBasisFcts()->getDegree() - 2;

    q = Quadrature::provideQuadrature(dim, degree);
    quadFast = FastQuadrature::provideFastQuadrature(basFct,
               *q,
               INIT_GRD_PHI);

    double relative = relErr;
    double maxErr = 0.0, h1Err2 = 0.0, h1Norm2 = 0.0;


    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(fe_space->getMesh(), -1,
                                         Mesh::FILL_COORDS |
                                         Mesh::CALL_LEAF_EL |
                                         Mesh::FILL_DET |
                                         Mesh::FILL_GRD_LAMBDA);
    while (elInfo)
    {
      elinfo = elInfo;
      int i, j;
      double err, err_2, h1_err_el, norm_el, norm2, det, exact;
      const WorldVector<double>* grdu_vec;
      mtl::dense_vector<WorldVector<double>> grduh_vec;

      grdu_vec = quadFast->getQuadrature()->grdFAtQp(grdErrU, NULL);
      det = elInfo->getDet();
      errUh->getGrdAtQPs(elinfo, NULL, quadFast, grduh_vec);

      int nPoints = quadFast->getNumPoints();
      int dow = Global::getGeo(WORLD);

      for (h1_err_el = i = 0; i < nPoints; i++)
      {
        for (err_2 = j = 0; j < dow; j++)
        {
          err = grdu_vec[i][j] - grduh_vec[i][j];
          err_2 += sqr(err);
        }
        h1_err_el += quadFast->getWeight(i)*err_2;
      }

      exact = det*h1_err_el;
      h1Err2 += exact;
      maxErr = std::max(maxErr, exact);

      if (writeInLeafData)
        elInfo->getElement()->setEstimation(exact, component);

      if (relative)
      {
        for (norm_el = i = 0; i < nPoints; i++)
        {
          for (norm2 = j = 0; j < dow; j++)
            norm2 += sqr(grdu_vec[i][j]);
          norm_el += quadFast->getWeight(i)*norm2;
        }
        h1Norm2 += det*norm_el;
      }

      elInfo = stack.traverseNext(elInfo);
    }

    if (relative)
    {
      double relNorm2 = h1Norm2 + 1.e-15;

      elInfo = stack.traverseFirst(fe_space->getMesh(), -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        double exact = elInfo->getElement()->getEstimation(component) / relNorm2;
        if (writeInLeafData)
          elinfo->getElement()->setEstimation(exact, component);
        elInfo = stack.traverseNext(elInfo);
      }

      h1Err2 /= relNorm2;
      maxErr /= relNorm2;
    }

    if (max)
      *max = maxErr;

    return sqrt(h1Err2);
  }


  template<typename T>
  double Error<T>::L2Err(const AbstractFunction<T, WorldVector<double>>& u,
                         const DOFVector<T>& uh,
                         int relErr,
                         double* max,
                         bool writeLeafData,
                         int comp)
  {
    FUNCNAME("Error<T>::L2Err()");

    const FiniteElemSpace* fe_space;
    Quadrature* q = NULL;
    writeInLeafData = writeLeafData;
    component = comp;

    if (!(pU = &u))
    {
      ERROR("no function u specified; doing nothing\n");
      return(0.0);
    }

    if (!(errUh = &uh)  ||  !(fe_space = uh.getFeSpace()))
    {
      ERROR("no discrete function or no fe_space for it; doing nothing\n");
      return(0.0);
    }

    if (!(basFct = fe_space->getBasisFcts()))
    {
      ERROR("no basis functions at discrete solution ; doing nothing\n");
      return(0.0);
    }

    int dim = fe_space->getMesh()->getDim();
    int deg = u.getDegree();
    int degree = deg ? deg :  2 * fe_space->getBasisFcts()->getDegree() - 2;

    q = Quadrature::provideQuadrature(dim, degree);
    quadFast = FastQuadrature::provideFastQuadrature(basFct, *q, INIT_PHI);

    double relative = relErr;
    double maxErr = 0.0, l2Err2 = 0.0, l2Norm2 = 0.0;
    double* u_vec = new double[quadFast->getQuadrature()->getNumPoints()];
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<double> uh_vec(nPoints);

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(fe_space->getMesh(), -1,
                                         Mesh::FILL_COORDS |
                                         Mesh::CALL_LEAF_EL |
                                         Mesh::FILL_DET |
                                         Mesh::FILL_GRD_LAMBDA);
    while (elInfo)
    {
      elinfo = elInfo;

      quadFast->getQuadrature()->fAtQp(errU, u_vec);
      errUh->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      double det = elInfo->getDet();
      double l2_err_el = 0.0, err = 0.0;

      for (int i = 0; i < nPoints; i++)
      {
        err = u_vec[i] - uh_vec[i];
        l2_err_el += quadFast->getWeight(i) * sqr(err);
      }

      double exact = det * l2_err_el;
      l2Err2 += exact;
      maxErr = std::max(maxErr, exact);

      if (writeInLeafData)
        elInfo->getElement()->setEstimation(exact, component);

      if (relative)
      {
        double norm_el = 0.0;
        for (int i = 0; i < nPoints; i++)
          norm_el += quadFast->getWeight(i) * sqr(u_vec[i]);
        l2Norm2 += det * norm_el;
      }

      elInfo = stack.traverseNext(elInfo);
    }

    delete [] u_vec;

    if (relative)
    {
      double relNorm2 = l2Norm2 + 1.e-15;

      elInfo = stack.traverseFirst(fe_space->getMesh(), -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        double exact = elInfo->getElement()->getEstimation(component) / relNorm2;
        if (writeInLeafData)
          elInfo->getElement()->setEstimation(exact, component);
        elInfo = stack.traverseNext(elInfo);
      }

      l2Err2 /= relNorm2;
    }

    if (max)
      *max = maxErr;

    return (sqrt(l2Err2));
  }

  template<typename T>
  double Error<T>::L2Err_ElementWise(const AbstractFunction<T, WorldVector<double>>& u,
                                     const DOFVector<T>& uh,
                                     double* max,
                                     bool writeLeafData,
                                     int comp,
                                     int level,
                                     std::map<int, double>& estMap)
  {
    FUNCNAME("Error<T>::L2Err_ElementWise()");

    const FiniteElemSpace* fe_space;
    Quadrature* q = NULL;
    writeInLeafData = writeLeafData;
    component = comp;

    if (!(pU = &u))
    {
      ERROR("no function u specified; doing nothing\n");
      return(0.0);
    }

    if (!(errUh = &uh)  ||  !(fe_space = uh.getFeSpace()))
    {
      ERROR("no discrete function or no fe_space for it; doing nothing\n");
      return(0.0);
    }

    if (!(basFct = fe_space->getBasisFcts()))
    {
      ERROR("no basis functions at discrete solution ; doing nothing\n");
      return(0.0);
    }

    int dim = fe_space->getMesh()->getDim();
    int deg = u.getDegree();
    int degree = deg ? deg :  2 * fe_space->getBasisFcts()->getDegree() - 2;

    q = Quadrature::provideQuadrature(dim, degree);
    quadFast = FastQuadrature::provideFastQuadrature(basFct, *q, INIT_PHI);

    double maxErr = 0.0, l2Err2 = 0.0;
    double* u_vec = new double[quadFast->getQuadrature()->getNumPoints()];
    int nPoints = quadFast->getNumPoints();
    mtl::dense_vector<double> uh_vec(nPoints);

    Flag flag = Mesh::FILL_COORDS | Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA;
    if (level == -1)
      flag |= Mesh::CALL_LEAF_EL;
    else
      flag |= Mesh::CALL_EL_LEVEL;

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(fe_space->getMesh(), level, flag);

    double a = 0.0;
    double b = 0.0;

    while (elInfo)
    {
      elinfo = elInfo;

      quadFast->getQuadrature()->fAtQp(errU, u_vec);
      errUh->getVecAtQPs(elInfo, NULL, quadFast, uh_vec);
      double det = elInfo->getDet();
      double l2_err_el = 0.0;

      double int_u_vec = 0.0;
      double int_uh_vec = 0.0;

      for (int i = 0; i < nPoints; i++)
      {
        int_u_vec += quadFast->getWeight(i) * u_vec[i];
        int_uh_vec += quadFast->getWeight(i) * uh_vec[i];
      }

      a += det * int_u_vec - det * int_uh_vec;
      b += fabs(det * int_u_vec - det * int_uh_vec);
      l2_err_el = pow(fabs(det * int_u_vec - det * int_uh_vec), 2.0);
      l2Err2 += l2_err_el;
      maxErr = std::max(maxErr, l2_err_el);

      if (writeInLeafData)
        elInfo->getElement()->setEstimation(l2_err_el, component);

      estMap[elInfo->getElement()->getIndex()] = l2_err_el;

      elInfo = stack.traverseNext(elInfo);
    }

    MSG("L2ERR values %e %e %e %e\n", a, b, l2Err2, sqrt(l2Err2));

    delete [] u_vec;

    if (max)
      *max = maxErr;

    return (sqrt(l2Err2));
  }

}
