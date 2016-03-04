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



/** \file Cholesky.h */

#ifndef CHOLESKY_H
#define CHOLESKY_H

#include "FixVec.h"

namespace AMDiS
{

  using namespace std;
  using namespace AMDiS;

  /**
   * \ingroup Solvers
   *
   * \brief
   * Solves a symmetric positive definite system using Cholesky decomposition.
   * Returns true if the matrix in (numerically) positive definite.
   */

  class Cholesky
  {
  public:
    /** \brief
     * Computes the Cholesky factor of a positive definite matrix A.
     * On input, only the upper triangle of A need be given; it is not modified.
     * The Cholesky factor is returned in the lower triangle of A, except for its
     * diagonal elements which are returned in P.
     */
    static bool factorization(Matrix<double>* A, Vector<double>* p);

    /** \brief
     * Solves system A*X=B, where A is a positive definite matrix.
     * If P=NULL; A is assumed to be positive definite, and a Cholesky
     * decomposition is computed using the previous routine.
     * If P is given, A and P are assumed to be already given as the output of
     * the previous routine.
     */
    static bool solve(Matrix<double>* A, Vector<double>* b, Vector<double>* x,
                      Vector<double>* p = NULL);
    static bool solve(Matrix<double>* A, Vector<WorldVector<double>>* b,
                      Vector<WorldVector<double>>* x,
                      Vector<double>* p = NULL);
  };

}

#endif
