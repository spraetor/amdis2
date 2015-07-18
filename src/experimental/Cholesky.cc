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


#include "Cholesky.h"
#include "MatrixVectorOperations.h"

namespace AMDiS {

using namespace std;
using namespace AMDiS;

bool Cholesky::factorization(Matrix<double> *A, Vector<double> *p)
{
  FUNCNAME("Cholesky::factorization()");

  int n = A->getNumRows();

  // Checking memory for vector P of diagonal elements of factorization.
  static Vector<double> *pT = NULL;

  if (p)
  {
    if (p->getSize() != n)
      p->resize(n);
  }
  else
  {
    if (pT)
    {
      if (pT->getSize() != n)
	pT->resize(n);
    }
    else
      pT = new Vector<double>(n);

    p = pT;
  }

  // Constructs the Cholesky factorization A = L*L, with L lower triangular.
  // Only the upper triangle need be given; it is not modified.
  // The Cholesky factor L is stored in the lower triangle of A, except for
  // its diagonal which is stored in P.

  int i, j, k;
  double sum;

  for (i=0; i<n; i++)
  {
    for (j=i; j<n; j++)
    {
      for (sum=(*A)[i][j], k=i-1; k>=0; k--)
	sum -= (*A)[i][k] * (*A)[j][k];

      if (i==j)
      {
	if (sum<=0)
	{
	  MSG("Matrix is (numerically) not positive definite!\n");
	  MSG("Cholesky decomposition does not work; choose another method for solving the system.\n");
	  return false;
	}
	(*p)[i] = sqrt(sum);
      }
      else
	(*A)[j][i] = sum / (*p)[i];
    }
  }

  return true;
}

bool Cholesky::solve(Matrix<double> *A, Vector<double> *b, Vector<double> *x,
		     Vector<double> *p)
{
  FUNCNAME("Cholesky::solve");

  bool success = true;
  int        n = b->getSize();
  TEST_EXIT(n == A->getNumRows())
    ("Dimensions of matrix and vector do not match!\n");

  // Checking memory for solution vector X.
  if (x && (x->getSize() != n))
    x->resize(n);

  if (!x)
    x = new Vector<double>(n);

  // Checking vector P.
  static Vector<double> *pT = NULL;

  if (!p || (p->getSize() != n))
  {
    if (pT && pT->getSize() != n)
      delete pT;

    if (!pT)
      pT = new Vector<double>(n);

    p  = pT;

    success = factorization(A, p);
  }

  // Now solve the system by backsubstitution.
  int i, k;
  double sum;

  for (i=0; i<n; i++)        // Solve L*Y = B, storing Y in X.
  {
    for (sum=(*b)[i], k=i-1; k>=0; k--)
      sum -= (*A)[i][k] * (*x)[k];
    (*x)[i] = sum / (*p)[i];
  }

  for (i=n-1; i>=0; i--)     // Solve L^T*X = Y.
  {
    for (sum=(*x)[i], k=i+1; k<n; k++)
      sum -= (*A)[k][i] * (*x)[k];
    (*x)[i] = sum / (*p)[i];
  }

  return success;
}

bool Cholesky::solve(Matrix<double> *A, Vector<WorldVector<double> > *b,
		     Vector<WorldVector<double> > *x, Vector<double> *p)
{
  FUNCNAME("Cholesky::solve");

  bool success = true;
  int        n = b->getSize();
  TEST_EXIT(n == A->getNumRows())
    ("Dimensions of matrix and vector do not match!\n");

  // Checking memory for solution vector X.
  if (x && (x->getSize() != n))
    x->resize(n);

  if (!x)
    x = new Vector<WorldVector<double> >(n);

  // Checking vector P.
  static Vector<double> *pT = NULL;

  if (!p || (p->getSize() != n))
  {
    pT = new Vector<double>(n);
    p  = pT;

    success = factorization(A, p);
  }

  // Now solve the system by backsubstitution.
  int i, k;
  WorldVector<double> vec_sum;

  for (i=0; i<n; i++)        // Solve L*Y = B, storing L in X.
  {
    for (vec_sum=(*b)[i], k=i-1; k>=0; k--)
      vec_sum -= (*x)[k] * (*A)[i][k] ;
    (*x)[i] = vec_sum * (1.0/(*p)[i]);
  }

  for (i=n-1; i>=0; i--)     // Solve L^T*X = Y.
  {
    for (vec_sum=(*x)[i], k=i+1; k<n; k++)
      vec_sum -= (*x)[k] * (*A)[k][i] ;
    (*x)[i] = vec_sum * (1.0/(*p)[i]);
  }

  return success;
}

}
