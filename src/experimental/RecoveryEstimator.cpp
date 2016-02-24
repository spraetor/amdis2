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


#include "RecoveryEstimator.h"
#include "Initfile.h"

namespace AMDiS
{

  RecoveryEstimator::RecoveryEstimator(std::string name, DOFVector<double>* uh_, int r)
    : Estimator(name, r),
      uh(uh_),
      relative(0),
      C(1.0),
      method(0),
      feSpace(NULL),
      f_vec(NULL),
      f_scal(NULL),
      aux_vec(NULL),
      rec_struct(NULL)
  {
    FUNCNAME("RecoveryEstimator::constructor()");

    Parameters::get(name + "->rec method", method); // 0, 1, or 2 (see Recovery.h)
    Parameters::get(name + "->rel error", relative); // 0 or 1
    Parameters::get(name + "->C", C);

    C = C > 1e-25 ? sqr(C) : 1.0;

    if (norm == H1_NORM)
    {
      feSpace = uh_->getFeSpace();
      degree = feSpace->getBasisFcts()->getDegree();

      if (degree <= 2 && C != 1.0)
      {
        WARNING("Recovery estimators in the H1_NORM usually very good for linear and quadratic finite element; normally you do not need to overwrite the default value of C\n");
        WAIT;
      }
    }
    else
    {
      degree = uh_->getFeSpace()->getBasisFcts()->getDegree() + 1;
      feSpace =
        FiniteElemSpace::provideFeSpace(NULL,
                                        Lagrange::getLagrange(uh_->getFeSpace()->getMesh()->getDim(),
                                            degree),
                                        uh_->getFeSpace()->getMesh(),
                                        name + "->feSpace");

      if (method == 2)
      {
        ERROR("Simple averaging only for the H1_NORM; using SPR instead\n");
        WAIT;
        method = 0;
      }
    }

    if (method == 2 && degree !=1)
    {
      ERROR("Simple averaging only for linear elements; using SPR instead\n");
      WAIT;
      method = 0;
    }

    rec_struct = new Recovery(norm, method);
  }

  void RecoveryEstimator::init(double ts)
  {
    basFcts = uh->getFeSpace()->getBasisFcts();
    int dim = mesh->getDim();
    h1Norm2 = 0.0;

    if (norm == H1_NORM)      // sets recovery gradient.
    {
      if (method == 2)
        rec_grd = rec_struct->recovery(uh, f_vec, f_scal, aux_vec);
      else
        rec_grd = rec_struct->recovery(uh, feSpace, f_vec, f_scal, aux_vec);

      rec_basFcts = rec_grd->getFeSpace()->getBasisFcts();
    }
    else                     // sets higher-order recovery solution.
    {
      rec_uh = rec_struct->recoveryUh(uh, feSpace);
      rec_basFcts = rec_uh->getFeSpace()->getBasisFcts();
    }

    int deg = 2 * std::max(basFcts->getDegree(), rec_basFcts->getDegree());
    quad = Quadrature::provideQuadrature(dim, deg);
    nPoints = quad->getNumPoints();

    grdAtQP = mtl::dense_vector<WorldVector<double>>(nPoints);
    recoveryGrdAtQP = mtl::dense_vector<WorldVector<double>>(nPoints);
    uhAtQP = mtl::dense_vector<double>(nPoints);
    recoveryUhAtQP = mtl::dense_vector<double>(nPoints);

    quadFast = FastQuadrature::provideFastQuadrature(basFcts, *quad, INIT_PHI | INIT_GRD_PHI);
    rec_quadFast = FastQuadrature::provideFastQuadrature(rec_basFcts, *quad, INIT_PHI | INIT_GRD_PHI);

    est_sum = 0.0;
    est_max = 0.0;
    est_t_sum = 0.0;
    est_t_max = 0.0;

    traverseFlag = Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA;

  }

  void RecoveryEstimator::exit(bool output)
  {
    FUNCNAME("RecoveryEstimator::exit()");

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double send_est_sum = est_sum;
    double send_est_max = est_max;

    MPI::COMM_WORLD.Allreduce(&send_est_sum, &est_sum, 1, MPI_DOUBLE, MPI_SUM);
    MPI::COMM_WORLD.Allreduce(&send_est_max, &est_max, 1, MPI_DOUBLE, MPI_MAX);
#endif

    // Computing relative errors
    if (relative)
    {
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
      while (elInfo)
      {
        double estEl = elInfo->getElement()->getEstimation(row);
        estEl /= h1Norm2;
        elInfo->getElement()->setEstimation(estEl, row);
        elInfo = stack.traverseNext(elInfo);
      }

      est_max /= h1Norm2;
      est_sum /= h1Norm2;
    }

    est_sum = sqrt(est_sum);

    if (output)
    {
      MSG("estimate for component %d = %.8e\n", row, est_sum);
    }
  }

  void RecoveryEstimator::estimateElement(ElInfo* elInfo, DualElInfo* dualElInfo)
  {
    Element* el = elInfo->getElement();
    double det = elInfo->getDet();
    double errEl = 0.0;
    double estEl = 0.0;
    int dow = Global::getGeo(WORLD);

    if (norm == H1_NORM)
    {
      // get gradient and recovery gradient at quadrature points
      uh->getGrdAtQPs(elInfo, NULL, quadFast, grdAtQP);
      rec_grd->getVecAtQPs(elInfo, NULL, rec_quadFast, recoveryGrdAtQP);
      if (f_scal)
      {
        if (aux_vec)
          aux_vec->getVecAtQPs(elInfo, NULL, quadFast, uhAtQP);
        else
          uh->getVecAtQPs(elInfo, NULL, quadFast, uhAtQP);
      }

      // calc h1 error
      for (int i = 0; i < nPoints; i++)
      {
        double  err2 = 0.0;
        double fAtQP = 1.0;
        if (f_scal)
          fAtQP = (*f_scal)(uhAtQP[i]);
        if (f_vec)
        {
          elInfo->coordToWorld(quad->getLambda(i), quad_pt);
          fAtQP = (*f_vec)(quad_pt);
        }

        for (int j = 0; j < dow; j++)
          err2 += sqr(recoveryGrdAtQP[i][j] - fAtQP * grdAtQP[i][j]);
        errEl += quad->getWeight(i) * err2;
      }
    }
    else
    {
      // get vector and recovery vector at quadrature points
      uh->getVecAtQPs(elInfo, NULL, quadFast, uhAtQP);
      rec_uh->getVecAtQPs(elInfo, NULL, rec_quadFast, recoveryUhAtQP);

      // calc l2 error
      for (int i = 0; i < nPoints; i++)
        errEl += quad->getWeight(i) * sqr(recoveryUhAtQP[i] - uhAtQP[i]);
    }

    estEl += C * det * errEl;
    el->setEstimation(estEl, row);
    est_sum += estEl;
    est_max = std::max(est_max, estEl);

    if (relative)
    {
      double normEl = 0.0;

      if (norm == H1_NORM)
      {
        for (int i = 0; i < nPoints; i++)
        {
          double norm2 = 0.0;
          for (int j = 0; j < dow; j++)
            norm2 += sqr(recoveryGrdAtQP[i][j]);
          normEl += quad->getWeight(i) * norm2;
        }
      }
      else
      {
        for (int i = 0; i < nPoints; i++)
          normEl += quad->getWeight(i) * sqr(recoveryUhAtQP[i]);
      }

      h1Norm2 += det * normEl;
    }
  }

}
