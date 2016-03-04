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



/** \file Recovery.h */

#ifndef AMDIS_RECOVERY_H
#define AMDIS_RECOVERY_H

#include <set>
#include "Lagrange.h"
#include "Traverse.h"
#include "DOFVector.h"
#include "Cholesky.h"

namespace AMDiS
{

  class Monomial : public
    BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double>>
  {
  public:
    Monomial(WorldVector<int> expon)
      : BinaryAbstractFunction<double, WorldVector<double>, WorldVector<double>>(),
        exponent(expon)
    {
      degree_ = exponent[0];
      for (int i = 1; i<exponent.getSize(); i++)
        degree_ += exponent[i];
    }

    virtual ~Monomial() {}

    double operator()(const WorldVector<double>& y,
                      const WorldVector<double>& z) const
    {
      double result = std::pow(y[0] - z[0], double(exponent[0]));
      for (int i = 1; i < exponent.getSize(); i++)
        result *= std::pow(y[i] - z[i], double(exponent[i]));

      return result;
    }

  private:
    WorldVector<int> exponent;
  };


  class RecoveryStructure
  {
  public:
    RecoveryStructure()
      : coords(NULL),
        A(NULL),
        rec_uh(NULL),
        rec_grdUh(NULL),
        neighbors(NULL)
    {}

    ~RecoveryStructure()
    {
      clear();
    }

    RecoveryStructure(const RecoveryStructure& rhs)
      : coords(NULL),
        A(NULL),
        rec_uh(NULL),
        rec_grdUh(NULL),
        neighbors(NULL)
    {
      *this = rhs;
    }

    RecoveryStructure& operator=(const RecoveryStructure& rhs);

    /// Clear recovery structure
    inline void clear()
    {
      if (coords != NULL)
      {
        delete coords;
        coords = NULL;
      }

      if (A != NULL)
      {
        delete A;
        A = NULL;
      }

      if (rec_uh != NULL)
      {
        delete rec_uh;
        rec_uh = NULL;
      }

      if (rec_grdUh != NULL)
      {
        delete rec_grdUh;
        rec_grdUh = NULL;
      }

      if (neighbors != NULL)
      {
        delete neighbors;
        neighbors = NULL;
      }
    }

    /// Prints recovery structure (for test purposes)
    void print();


  private:
    /// World coordinates of the node.
    WorldVector<double>* coords;

    /// For interior nodes.
    Matrix<double>* A;
    Vector<double>* rec_uh;
    Vector<WorldVector<double>>* rec_grdUh;

    /// For boundary, edge, face, of center nodes: interior neighbors nodes.
    /// For interior vertices in uh-recovery: neighbors nodes.
    std::set<DegreeOfFreedom>* neighbors;

    friend class Recovery;
  };


  /**
   * \ingroup Recovery
   *
   * \brief
   * Recovering the gradient of a finite element function.
   */
  class Recovery
  {
  public:
    Recovery(int norm, int method_)
      : struct_vec(NULL),
        feSpace(NULL),
        matrix_fcts(NULL),
        method(method_)
    {
      n_monomials = 0;
      gradient = (norm == H1_NORM);
    }

    /// Clear vector of recovery structures
    inline void clear()
    {
      if (struct_vec)
      {
        DOFVector<RecoveryStructure>::Iterator SV_it(struct_vec, ALL_DOFS);
        for (SV_it.reset(); !SV_it.end(); ++SV_it)
          (*SV_it).clear();
      }
    }

    /// Recovers flux or gradient of given DOFVector.
    DOFVector<WorldVector<double>>*
                                recovery(DOFVector<double>* uh,
                                         AbstractFunction<double, WorldVector<double>>* f_vec = NULL,
                                         AbstractFunction<double, double>* f_scal = NULL,
                                         DOFVector<double>* aux_vec = NULL);

    DOFVector<WorldVector<double>>*
                                recovery(DOFVector<double>* uh, const FiniteElemSpace* fe_space,
                                         AbstractFunction<double, WorldVector<double>>* f_vec = NULL,
                                         AbstractFunction<double, double>* f_scal = NULL,
                                         DOFVector<double>* aux_vec = NULL);

    /// Computes higher order approximation of given DOFVector.
    void recoveryUh(DOFVector<double>* uh, DOFVector<double>& rec_vec);

    DOFVector<double>* recoveryUh(DOFVector<double>* uh,
                                  const FiniteElemSpace* fe_space);

    /// For test purposes
    void test(DOFVector<double>* uh, const FiniteElemSpace* fe_space);

  private:
    /// Defines a new finite element space if necessary.
    void set_feSpace(const FiniteElemSpace* fe_space);

    /// Fills vector of exponents of the monomials.
    int set_exponents(int degree);

    /// Fills vector of recovery structures.
    void fill_struct_vec(DOFVector<double>* uh,
                         AbstractFunction<double, WorldVector<double>>* f_vec = NULL,
                         AbstractFunction<double, double>* f = NULL,
                         DOFVector<double>* aux_vec = NULL);

    /// Compute integrals defining matrix and vector on elemen (continuous ZZ-recovery)
    void compute_integrals(DOFVector<double>* uh, ElInfo* elInfo,
                           RecoveryStructure* rec_struct,
                           AbstractFunction<double, WorldVector<double>>* f_vec = NULL,
                           AbstractFunction<double, double>* f_scal = NULL,
                           DOFVector<double>* aux_vec = NULL);

    /// Compute integrals defining matrix and vector on element (superconvergent patch recovery)
    void compute_interior_sums(DOFVector<double>* uh, ElInfo* elInfo,
                               RecoveryStructure* rec_struct, Quadrature* quad,
                               AbstractFunction<double, WorldVector<double>>* f_vec = NULL,
                               AbstractFunction<double, double>* f_scal = NULL,
                               DOFVector<double>* aux_vec = NULL);

    void compute_node_sums(DOFVector<double>* uh, ElInfo* elInfo,
                           RecoveryStructure* rec_struct, DimVec<int> preDOFs,
                           int n_vertices, int n_edges, int n_faces);

    void compute_sums_linear(DOFVector<double>* uh, ElInfo* elInfo,
                             RecoveryStructure* rec_struct,
                             int vertex, DimVec<int> preDOFs, int n_vertices);

  private:
    /// Structure storing needed information.
    DOFVector<RecoveryStructure>* struct_vec;

    /// Working finite element space.
    const FiniteElemSpace* feSpace;

    /// Number of monomials.
    int n_monomials;

    /// Exponents of the monomials.
    Vector<WorldVector<int>> exponents;

    /// Functions for system matrix.
    Matrix<Monomial*>* matrix_fcts;

    /** \brief
     * True if gradient or flux should be recovered.
     * False if seeking for higher order approximation of uh.
     */
    bool gradient;

    /** \brief
     * 0: superconvergent patch recovery (discrete ZZ)
     * 1: local L2-averaging (continuous ZZ-recovery)
     * 2: simple averaging
     */
    int method;
  };

}
#endif
