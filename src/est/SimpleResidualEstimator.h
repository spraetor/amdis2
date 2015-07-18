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



/** \file SimpleResidualEstimator.h */

/** \defgroup Estimator Estimator module
 * @{ <img src="estimator.png"> @}
 */

#ifndef AMDIS_SIMPLE_RESIDUAL_ESTIMATOR_H
#define AMDIS_SIMPLE_RESIDUAL_ESTIMATOR_H

#include "Estimator.h"
#include "FixVec.h"

namespace AMDiS {

  /**
   * \ingroup Estimator
   * 
   * \brief
   * Simple Residual estimator.
   *
   * The restrictions for using the simple residual estimator are the following:
   *   - no system of PDEs, thus a PDE with only one unknown
   *   - no first order terms
   *   - no time dependent terms
   */
  class SimpleResidualEstimator : public Estimator
  {
  public:
    class Creator : public EstimatorCreator
    {
    public:
      Creator() 
	: EstimatorCreator() 
      {}

      virtual ~Creator() 
      {}

      /// Returns a new ODirSolver object.
      Estimator* create() 
      { 
	return new SimpleResidualEstimator(name);
      }
    };
  
    /// Constructor.
    SimpleResidualEstimator(std::string name);

    /// Init the error estimator.
    void init(double);

    /// Run estimator on one element. 
    /// \param[in]   elInfo      Info object for the element to be estimated.
    /// \param[in]   dualElInfo  Not used here. In general, this may be used for
    ///                          estimating with the multi-mesh technique.
    void estimateElement(ElInfo *elInfo, DualElInfo *dualElInfo = NULL);

    /// Finalize the error estimator, i.e., delete all temporary data structures.
    void exit(bool output = true);

    /// Returns pow(det,2.0/dim).
    inline double h2_from_det(double det, int dim) 
    {
      return pow(det, 2.0 / dim);
    }

  protected:
    /// Computes the element residual for a given element.
    double computeElementResidual(ElInfo *elInfo);

    /// Computes the jump residual for a given element.
    double computeJumpResidual(ElInfo *elInfo);

    /** \brief
     * Returns residual square at quadrature point. Not Member of
     * Estimator to avoid multiple instantiation.
     *
     * \param[in]  elInfo   Current element information object.
     * \param[in]  nPoints  Number of quadrature points on this element.
     * \param[in]  uhIq     Solution vector on element quadrature points.
     * \param[in]  D2UhIq   Second derivations of the solution vector on
     *                      element quadrature points.
     * \param[in]  A        Matrix that contains the left hand side operators.
     * \param[in]  fh       Vector that contains the right hand side operators.
     * \param[in]  quad     Object for numerical quadrature.
     * \param[out] result   Vector containing the residual on the quadrature 
     *                      points.
     */
    void r(const ElInfo *elInfo,
	   int nPoints,
	   const mtl::dense_vector<double>& uhIq,
	   const mtl::dense_vector<WorldMatrix<double> > &D2UhIq,
	   DOFMatrix *A, 
	   DOFVector<double> *fh,
	   Quadrature *quad,
	   mtl::dense_vector<double>& result);
 
  protected:
    /// Constant in front of element residual.
    double C0;
 
    /// Constant in front of edge/face residual.
    double C1;

    /// Number of quadrature points.
    int nPoints;

    /// Dimension of the mesh.
    int dim;

    /// Polynomial degree.
    int degree;

    /// Object for numerical quadrature
    Quadrature *quad;

    /// Object for fast numerical quadrature
    FastQuadrature *quadFast;

    /// Pointer to the basis functions of the FE space.
    const BasisFunction *basFcts;

    /// Vector that stores all global DOFs of one element.
    ElementVector uhEl;

    /// Vector that stores all global DOFs of one element.
    ElementVector uhNeigh;

    /// Vector that stores values on all quadrature points (QP) on one element.
    ElementVector uhQP;

    mtl::dense_vector<WorldVector<double> > grdUhQp;

    /// Matrix that stores the second derivations on all quadrature points
    /// on one element.
    mtl::dense_vector<WorldMatrix<double> > D2UhQp;

    /// Stores the element residual computed at the quadrature points of 
    /// the element.
    ElementVector riq;

    /// Pointer to the information object of some neighbouring element.
    ElInfo *neighInfo;

    /// Surface quadrature object, used for computing the jump residual.
    Quadrature *surfaceQuad;

    /// Number of quadrature points of the surface quadrature object.
    int nPointsSurface;

    /// Stores on all surface quadrature points the gradient of a function.
    std::vector<WorldVector<double> > grdUhEl;

    /// Stores on all surface quadrature points the gradient of a function.
    std::vector<WorldVector<double> > grdUhNeigh;

    /// Stores on all surface quadrature points the jump of a function.
    std::vector<WorldVector<double> > jump;

    /// Stores on all surface quadrature points the jump of a function.
    std::vector<WorldVector<double> > localJump;

    /// Vector to store, for a given element, all its neighbouring element indices.
    WorldVector<int> faceIndEl;

    /// Vector to store, for a given element, all its neighbouring element indices.
    WorldVector<int> faceIndNeigh;

    DimVec<WorldVector<double> > *lambdaNeigh;
    
    DimVec<double> *lambda;

    /// Maximal number of neighbours an element may have in the used dimension.
    int nNeighbours;

    double kappa_inv;
  };
}

#endif // AMDIS_SIMPLE_RESIDUAL_ESTIMATOR_H
